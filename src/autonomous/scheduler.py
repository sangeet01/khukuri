"""
Khukuri PINCER Scheduler.
Pairs cron jobs with LLM reasoning to automate discovery.
"""

import logging
from .controller import AutonomousLoopController, Mode

logger = logging.getLogger('khukuri.autonomous')

class PINCERScheduler:
    """
    Scheduler for PINCER autonomous runs.
    
    Acts as the entry point for cron-triggered discovery cycles.
    """
    def __init__(self, project: str = "default"):
        self.project = project
        self.controller = AutonomousLoopController(
            mode=Mode.LIGHTS_OUT,
            project=project
        )

    def run_cycle(self):
        """Execute a single autonomous cycle."""
        logger.info(f"PINCERScheduler: running cycle for {self.project}")
        return self.controller.run_cycle()
