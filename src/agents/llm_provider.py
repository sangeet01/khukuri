"""
Multi-LLM Provider abstraction for Khukuri agents.

Supports: Anthropic Claude, OpenAI, Google Gemini, Ollama (local).
All providers expose a unified .complete(messages, system, **kwargs) -> str interface.
"""

import json
import logging
from abc import ABC, abstractmethod
from typing import Any, Dict, List, Optional

logger = logging.getLogger('khukuri')


# ---------------------------------------------------------------------------
# Base
# ---------------------------------------------------------------------------

class LLMProvider(ABC):
    """Abstract base for all LLM backends."""

    @abstractmethod
    def complete(
        self,
        messages: List[Dict[str, str]],
        system: Optional[str] = None,
        temperature: float = 0.7,
        max_tokens: int = 1000,
    ) -> str:
        """Return the assistant's reply as a plain string."""

    def complete_json(
        self,
        messages: List[Dict[str, str]],
        system: Optional[str] = None,
        temperature: float = 0.3,
        max_tokens: int = 1500,
    ) -> Dict[str, Any]:
        """
        Return parsed JSON from the model.
        Strips markdown fences automatically.
        """
        raw = self.complete(messages, system=system,
                            temperature=temperature, max_tokens=max_tokens)
        return _parse_json(raw)

    @property
    def name(self) -> str:
        return self.__class__.__name__


# ---------------------------------------------------------------------------
# Anthropic Claude
# ---------------------------------------------------------------------------

class AnthropicProvider(LLMProvider):
    """Anthropic Claude via the official anthropic SDK."""

    DEFAULT_MODEL = "claude-sonnet-4-20250514"

    def __init__(self, api_key: Optional[str] = None, model: Optional[str] = None):
        try:
            import anthropic
        except ImportError:
            raise ImportError("pip install anthropic")

        self._client = anthropic.Anthropic(api_key=api_key)
        self.model = model or self.DEFAULT_MODEL

    def complete(self, messages, system=None, temperature=0.7, max_tokens=1000):
        kwargs: Dict[str, Any] = dict(
            model=self.model,
            max_tokens=max_tokens,
            temperature=temperature,
            messages=messages,
        )
        if system:
            kwargs["system"] = system

        response = self._client.messages.create(**kwargs)
        return response.content[0].text

    @property
    def name(self):
        return f"Anthropic({self.model})"


# ---------------------------------------------------------------------------
# OpenAI
# ---------------------------------------------------------------------------

class OpenAIProvider(LLMProvider):
    """OpenAI ChatCompletion via the openai SDK."""

    DEFAULT_MODEL = "gpt-4o-mini"

    def __init__(self, client=None, api_key: Optional[str] = None,
                 model: Optional[str] = None):
        """
        Pass an existing openai.OpenAI client *or* an api_key string.
        Accepts a pre-built client so existing code that passes openai_client
        keeps working with zero changes.
        """
        if client is not None:
            self._client = client
        else:
            try:
                import openai
            except ImportError:
                raise ImportError("pip install openai")
            self._client = openai.OpenAI(api_key=api_key)

        self.model = model or self.DEFAULT_MODEL

    def complete(self, messages, system=None, temperature=0.7, max_tokens=1000):
        full_messages = []
        if system:
            full_messages.append({"role": "system", "content": system})
        full_messages.extend(messages)

        response = self._client.chat.completions.create(
            model=self.model,
            messages=full_messages,
            temperature=temperature,
            max_tokens=max_tokens,
        )
        return response.choices[0].message.content

    @property
    def name(self):
        return f"OpenAI({self.model})"


# ---------------------------------------------------------------------------
# Google Gemini
# ---------------------------------------------------------------------------

class GeminiProvider(LLMProvider):
    """Google Gemini via google-generativeai SDK."""

    DEFAULT_MODEL = "gemini-1.5-flash"

    def __init__(self, api_key: Optional[str] = None, model: Optional[str] = None):
        try:
            import google.generativeai as genai
        except ImportError:
            raise ImportError("pip install google-generativeai")

        if api_key:
            genai.configure(api_key=api_key)
        self.model = model or self.DEFAULT_MODEL
        self._genai = genai

    def complete(self, messages, system=None, temperature=0.7, max_tokens=1000):
        model = self._genai.GenerativeModel(
            self.model,
            system_instruction=system,
        )
        history = []
        for m in messages[:-1]:
            role = "model" if m["role"] == "assistant" else m["role"]
            history.append({"role": role, "parts": [m["content"]]})

        chat = model.start_chat(history=history)
        response = chat.send_message(
            messages[-1]["content"],
            generation_config=self._genai.GenerationConfig(
                temperature=temperature,
                max_output_tokens=max_tokens,
            ),
        )
        return response.text

    @property
    def name(self):
        return f"Gemini({self.model})"


# ---------------------------------------------------------------------------
# Ollama (local)
# ---------------------------------------------------------------------------

class OllamaProvider(LLMProvider):
    """Local models via Ollama's REST API (no SDK needed)."""

    def __init__(self, model: str = "llama3", base_url: str = "http://localhost:11434"):
        self.model = model
        self.base_url = base_url.rstrip("/")

    def complete(self, messages, system=None, temperature=0.7, max_tokens=1000):
        import urllib.request

        full_messages = []
        if system:
            full_messages.append({"role": "system", "content": system})
        full_messages.extend(messages)

        payload = json.dumps({
            "model": self.model,
            "messages": full_messages,
            "stream": False,
            "options": {"temperature": temperature, "num_predict": max_tokens},
        }).encode()

        req = urllib.request.Request(
            f"{self.base_url}/api/chat",
            data=payload,
            headers={"Content-Type": "application/json"},
        )
        with urllib.request.urlopen(req, timeout=120) as resp:
            data = json.loads(resp.read())
        return data["message"]["content"]

    @property
    def name(self):
        return f"Ollama({self.model})"


# ---------------------------------------------------------------------------
# Fallback (no LLM)
# ---------------------------------------------------------------------------

class FallbackProvider(LLMProvider):
    """
    Returns deterministic rule-based responses when no LLM is configured.
    Ensures the system always works, even without any API key.
    """

    def complete(self, messages, system=None, temperature=0.7, max_tokens=1000):
        last = messages[-1]["content"] if messages else ""
        return json.dumps({
            "role": "Fallback Agent",
            "analysis": f"Rule-based response to: {last[:80]}",
            "recommendations": [
                "Proceed with standard workflow",
                "Validate computationally before synthesis",
                "Monitor for resistance markers",
            ],
            "confidence": "low",
            "note": "No LLM configured — using deterministic fallback.",
        })

    @property
    def name(self):
        return "FallbackProvider"


# ---------------------------------------------------------------------------
# Factory
# ---------------------------------------------------------------------------

def create_provider(
    provider: str = "auto",
    api_key: Optional[str] = None,
    model: Optional[str] = None,
    openai_client=None,
    base_url: Optional[str] = None,
) -> LLMProvider:
    """
    Create an LLMProvider by name.

    provider values:
        "anthropic"  → AnthropicProvider
        "openai"     → OpenAIProvider
        "gemini"     → GeminiProvider
        "ollama"     → OllamaProvider
        "fallback"   → FallbackProvider
        "auto"       → try Anthropic → OpenAI → Fallback
    """
    p = provider.lower()

    if p == "anthropic":
        return AnthropicProvider(api_key=api_key, model=model)

    if p == "openai":
        return OpenAIProvider(client=openai_client, api_key=api_key, model=model)

    if p == "gemini":
        return GeminiProvider(api_key=api_key, model=model)

    if p == "ollama":
        return OllamaProvider(
            model=model or "llama3",
            base_url=base_url or "http://localhost:11434",
        )

    if p == "fallback":
        return FallbackProvider()

    if p == "auto":
        if openai_client is not None:
            logger.info("LLMProvider: using supplied OpenAI client")
            return OpenAIProvider(client=openai_client, model=model)

        import os
        if os.getenv("ANTHROPIC_API_KEY"):
            logger.info("LLMProvider: auto-detected Anthropic")
            return AnthropicProvider(model=model)
        if os.getenv("OPENAI_API_KEY"):
            logger.info("LLMProvider: auto-detected OpenAI")
            return OpenAIProvider(api_key=os.getenv("OPENAI_API_KEY"), model=model)
        if os.getenv("GOOGLE_API_KEY"):
            logger.info("LLMProvider: auto-detected Gemini")
            return GeminiProvider(api_key=os.getenv("GOOGLE_API_KEY"), model=model)

        logger.warning("LLMProvider: no API key found — using FallbackProvider")
        return FallbackProvider()

    raise ValueError(f"Unknown provider '{provider}'. "
                     "Choose: anthropic, openai, gemini, ollama, fallback, auto")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _parse_json(text: str) -> Dict[str, Any]:
    """Strip markdown fences then parse JSON."""
    text = text.strip()
    for fence in ("```json", "```"):
        if fence in text:
            text = text.split(fence, 1)[1].rsplit("```", 1)[0].strip()
            break
    try:
        return json.loads(text)
    except json.JSONDecodeError as exc:
        logger.warning(f"JSON parse failed: {exc}. Returning raw text.")
        return {"raw_response": text, "parse_error": str(exc)}
