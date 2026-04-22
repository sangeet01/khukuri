"""Base agent class — uses unified LLMProvider interface."""

import json
import logging
from typing import Any, Dict, Optional

from .llm_provider import LLMProvider, FallbackProvider, create_provider

logger = logging.getLogger('khukuri')


class BaseAgent:
    """
    Base class for all Khukuri agents.

    Backwards-compatible: still accepts `openai_client` for existing code,
    but internally wraps it in an OpenAIProvider.
    """

    def __init__(
        self,
        role: str,
        expertise: str,
        provider: Optional[LLMProvider] = None,
        openai_client=None,          # legacy param — auto-wrapped
    ):
        self.role = role
        self.expertise = expertise
        self.memory = []

        if provider is not None:
            self.provider = provider
        elif openai_client is not None:
            self.provider = create_provider("openai", openai_client=openai_client)
        else:
            self.provider = FallbackProvider()

        self.openai_client = openai_client   # kept for legacy attribute access

    def analyze(self, data: Dict[str, Any], question: str) -> Dict[str, Any]:
        """Analyse `data` and answer `question`. Returns a dict."""
        system = (
            f"You are a {self.role} with expertise in {self.expertise}. "
            "Return only valid JSON."
        )
        prompt = (
            f"Question: {question}\n\n"
            f"Data: {json.dumps(data, default=str)}\n\n"
            "Return JSON with keys: role, analysis, recommendations, confidence"
        )
        try:
            result = self.provider.complete_json(
                [{"role": "user", "content": prompt}],
                system=system,
                temperature=0.7,
                max_tokens=500,
            )
            self.memory.append({"question": question, "result": result})
            return result
        except Exception as exc:
            logger.warning(f"{self.role} analysis failed: {exc}")
            return self._fallback_analyze(data, question)

    def _fallback_analyze(self, data, question):
        return {
            "role": self.role,
            "analysis": f"Analysis of {question}",
            "recommendations": ["Proceed with validation", "Monitor results"],
            "confidence": "medium",
        }

    # Legacy shims
    def _ai_analyze(self, data, question):
        return self.analyze(data, question)
