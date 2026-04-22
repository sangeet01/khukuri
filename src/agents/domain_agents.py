"""
Domain-specialist agents for Khukuri.

All agents inherit from DomainAgent which uses the unified LLMProvider
interface.
"""

import json
import logging
from typing import Any, Dict, List, Optional

from .llm_provider import LLMProvider, FallbackProvider

logger = logging.getLogger('khukuri')


class DomainAgent:
    """Base class for all domain-specialist agents."""

    role: str = "Domain Specialist"
    expertise: str = "general drug discovery"
    _system_template: str = (
        "You are a {role} with deep expertise in {expertise}. "
        "Respond only with valid JSON unless asked otherwise."
    )

    def __init__(self, provider: Optional[LLMProvider] = None):
        self.provider = provider or FallbackProvider()
        self.memory: List[Dict[str, Any]] = []
        self._system = self._system_template.format(
            role=self.role, expertise=self.expertise
        )

    def analyze(self, data: Dict[str, Any], question: str) -> Dict[str, Any]:
        """Analyse `data` and answer `question`. Returns a dict."""
        prompt = self._build_prompt(data, question)
        messages = [{"role": "user", "content": prompt}]

        try:
            result = self.provider.complete_json(
                messages, system=self._system, temperature=0.7
            )
            self.memory.append({"question": question, "result": result})
            return result
        except Exception as exc:
            logger.warning(f"{self.role} analysis failed: {exc}")
            return self._fallback(question)

    def critique(self, response: Dict[str, Any], task: str) -> Dict[str, Any]:
        """Provide critical scientific feedback on another agent's response."""
        prompt = (
            f"Critically review this response for scientific rigour.\n"
            f"Task: {task}\n"
            f"Response: {json.dumps(response, default=str)}\n\n"
            "Return JSON with keys: strengths, weaknesses, suggested_improvements, "
            "overall_quality (high/medium/low)."
        )
        try:
            return self.provider.complete_json(
                [{"role": "user", "content": prompt}],
                system=self._system,
                temperature=0.3,
            )
        except Exception as exc:
            logger.warning(f"Critique failed: {exc}")
            return {"overall_quality": "unknown", "weaknesses": [str(exc)]}

    def debate(self, other_view: str, topic: str) -> Dict[str, Any]:
        """Respond to another agent's view on `topic`."""
        prompt = (
            f"Topic: {topic}\n\n"
            f"A colleague said:\n{other_view}\n\n"
            f"As a {self.role}, provide your counterpoint or agreement. "
            "Return JSON: {agrees: bool, points: [...], evidence: [...], revised_position: str}"
        )
        try:
            return self.provider.complete_json(
                [{"role": "user", "content": prompt}],
                system=self._system,
                temperature=0.8,
            )
        except Exception as exc:
            logger.warning(f"Debate failed: {exc}")
            return {"agrees": False, "points": [str(exc)]}

    def _build_prompt(self, data: Dict[str, Any], question: str) -> str:
        return (
            f"Question: {question}\n\n"
            f"Data:\n{json.dumps(data, indent=2, default=str)}\n\n"
            "Return a JSON object with at minimum: "
            "{assessment, concerns, recommendations, confidence}"
        )

    def _fallback(self, question: str) -> Dict[str, Any]:
        return {
            "role": self.role,
            "assessment": f"Rule-based {self.role} analysis",
            "concerns": ["Requires manual review"],
            "recommendations": ["Proceed with standard protocol"],
            "confidence": "low",
            "question": question,
        }


class ChemistAgent(DomainAgent):
    role = "Medicinal Chemist"
    expertise = "small molecule drug design, synthetic accessibility, bioisosteres, SAR analysis"


class BiologistAgent(DomainAgent):
    role = "Pharmacologist"
    expertise = "target biology, mechanism of action, selectivity profiling, pharmacology"


class ToxicologistAgent(DomainAgent):
    role = "Toxicologist"
    expertise = "drug safety, structural alerts, hERG liability, hepatotoxicity, regulatory toxicology"


class MicrobiologyAgent(DomainAgent):
    role = "Microbiologist"
    expertise = "pathogen biology, antimicrobial resistance mechanisms, MIC analysis, efflux pumps"


class GenomicsAgent(DomainAgent):
    role = "Genomics Specialist"
    expertise = "resistance genomics, mutation analysis, horizontal gene transfer, evolutionary modelling"


class CheminformaticsAgent(DomainAgent):
    role = "Cheminformatics Specialist"
    expertise = "QSAR modelling, molecular fingerprints, chemical space analysis, scaffold hopping"


class ClinicalAgent(DomainAgent):
    role = "Clinical Scientist"
    expertise = "translational medicine, clinical trial design, PK/PD modelling, dosing strategy"


AGENT_REGISTRY: Dict[str, type] = {
    "chemist": ChemistAgent,
    "biologist": BiologistAgent,
    "toxicologist": ToxicologistAgent,
    "microbiologist": MicrobiologyAgent,
    "genomics": GenomicsAgent,
    "cheminformatics": CheminformaticsAgent,
    "clinical": ClinicalAgent,
}


def create_agent(role: str, provider: Optional[LLMProvider] = None) -> DomainAgent:
    """Instantiate a domain agent by role name (case-insensitive)."""
    cls = AGENT_REGISTRY.get(role.lower())
    if cls is None:
        raise ValueError(
            f"Unknown agent role '{role}'. "
            f"Available: {list(AGENT_REGISTRY.keys())}"
        )
    return cls(provider=provider)
