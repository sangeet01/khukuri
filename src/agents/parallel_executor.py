"""
Async/parallel execution engine for Khukuri agent workflows.

Provides two execution modes:
  1. AsyncExecutor  – asyncio-based, for I/O-bound LLM calls
  2. ParallelExecutor – ThreadPoolExecutor-based, safe drop-in for sync code
"""

import asyncio
import logging
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from typing import Any, Callable, Dict, List, Optional

logger = logging.getLogger('khukuri')


@dataclass
class TaskResult:
    task_id: str
    success: bool
    output: Any = None
    error: Optional[str] = None
    duration_s: float = 0.0
    metadata: Dict[str, Any] = field(default_factory=dict)


@dataclass
class Task:
    task_id: str
    fn: Callable
    args: tuple = field(default_factory=tuple)
    kwargs: Dict[str, Any] = field(default_factory=dict)
    description: str = ""
    depends_on: List[str] = field(default_factory=list)


class ProgressReporter:
    """Simple callback-based progress tracker."""

    def __init__(self, total: int, on_update: Optional[Callable] = None):
        self.total = total
        self.completed = 0
        self.failed = 0
        self._on_update = on_update or (lambda msg: logger.info(msg))

    def update(self, result: TaskResult):
        self.completed += 1
        if not result.success:
            self.failed += 1
        pct = int(100 * self.completed / self.total)
        status = "✓" if result.success else "✗"
        self._on_update(
            f"[{pct:3d}%] {status} {result.task_id} ({result.duration_s:.1f}s)"
        )

    @property
    def summary(self) -> str:
        return (f"{self.completed}/{self.total} tasks — "
                f"{self.completed - self.failed} succeeded, "
                f"{self.failed} failed")


class ParallelExecutor:
    """
    Runs tasks concurrently using a thread pool.
    Respects task dependencies: tasks with `depends_on` run only after their
    prerequisites complete successfully.
    """

    def __init__(self, max_workers: int = 6,
                 on_progress: Optional[Callable] = None):
        self.max_workers = max_workers
        self.on_progress = on_progress

    def run_parallel(
        self,
        tasks: List[Task],
        timeout_per_task: float = 120.0,
    ) -> List[TaskResult]:
        """Execute tasks, honouring dependencies, return results in submission order."""
        results: Dict[str, TaskResult] = {}
        pending = {t.task_id: t for t in tasks}
        reporter = ProgressReporter(len(tasks), self.on_progress)

        while pending:
            ready = [
                t for t in pending.values()
                if all(dep in results and results[dep].success for dep in t.depends_on)
            ]
            blocked_by_failure = [
                t for t in pending.values()
                if any(dep in results and not results[dep].success for dep in t.depends_on)
            ]
            for t in blocked_by_failure:
                results[t.task_id] = TaskResult(
                    task_id=t.task_id, success=False,
                    error="Skipped: dependency failed",
                )
                reporter.update(results[t.task_id])
                del pending[t.task_id]

            if not ready:
                if pending:
                    logger.error("Dependency deadlock detected")
                break

            with ThreadPoolExecutor(max_workers=min(self.max_workers, len(ready))) as pool:
                future_to_task = {pool.submit(self._run_task, t): t for t in ready}
                for future in as_completed(future_to_task, timeout=timeout_per_task * len(ready)):
                    result = future.result()
                    results[result.task_id] = result
                    reporter.update(result)
                    del pending[result.task_id]

        logger.info(f"Parallel execution complete: {reporter.summary}")
        id_order = [t.task_id for t in tasks]
        return [results[tid] for tid in id_order if tid in results]

    @staticmethod
    def _run_task(task: Task) -> TaskResult:
        start = time.time()
        try:
            output = task.fn(*task.args, **task.kwargs)
            return TaskResult(
                task_id=task.task_id,
                success=True,
                output=output,
                duration_s=time.time() - start,
                metadata={"description": task.description},
            )
        except Exception as exc:
            logger.warning(f"Task '{task.task_id}' failed: {exc}")
            return TaskResult(
                task_id=task.task_id,
                success=False,
                error=str(exc),
                duration_s=time.time() - start,
            )


class AsyncExecutor:
    """Run tasks concurrently using asyncio."""

    def __init__(self, max_concurrency: int = 6,
                 on_progress: Optional[Callable] = None):
        self.max_concurrency = max_concurrency
        self.on_progress = on_progress

    async def run_parallel_async(
        self,
        tasks: List[Task],
        timeout_per_task: float = 120.0,
    ) -> List[TaskResult]:
        semaphore = asyncio.Semaphore(self.max_concurrency)
        reporter = ProgressReporter(len(tasks), self.on_progress)
        loop = asyncio.get_event_loop()

        async def run_one(task: Task) -> TaskResult:
            async with semaphore:
                start = time.time()
                try:
                    if asyncio.iscoroutinefunction(task.fn):
                        output = await asyncio.wait_for(
                            task.fn(*task.args, **task.kwargs),
                            timeout=timeout_per_task,
                        )
                    else:
                        output = await loop.run_in_executor(
                            None, lambda: task.fn(*task.args, **task.kwargs)
                        )
                    result = TaskResult(
                        task_id=task.task_id, success=True,
                        output=output, duration_s=time.time() - start,
                    )
                except Exception as exc:
                    result = TaskResult(
                        task_id=task.task_id, success=False,
                        error=str(exc), duration_s=time.time() - start,
                    )
                reporter.update(result)
                return result

        results = await asyncio.gather(*[run_one(t) for t in tasks])
        logger.info(f"Async execution complete: {reporter.summary}")
        return list(results)

    def run_parallel(self, tasks: List[Task], **kwargs) -> List[TaskResult]:
        """Sync wrapper — creates a new event loop if needed."""
        try:
            loop = asyncio.get_event_loop()
            if loop.is_running():
                import nest_asyncio
                nest_asyncio.apply()
        except RuntimeError:
            loop = asyncio.new_event_loop()
            asyncio.set_event_loop(loop)

        return loop.run_until_complete(self.run_parallel_async(tasks, **kwargs))


def make_executor(
    mode: str = "parallel",
    max_workers: int = 6,
    on_progress: Optional[Callable] = None,
):
    """Create an executor. mode: "parallel" (thread-based, default) | "async" (asyncio-based)"""
    if mode == "async":
        return AsyncExecutor(max_concurrency=max_workers, on_progress=on_progress)
    return ParallelExecutor(max_workers=max_workers, on_progress=on_progress)
