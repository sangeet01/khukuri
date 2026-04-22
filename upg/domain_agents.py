"""
Domain-specialist agents for Khukuri.

All agents inherit from DomainAgent which uses the unified LLMProvider
interface — no more hardcoded OpenAI or broken API signatures.

Agents:
    ChemistAgent          – medicinal chemistry & synthesis
    BiologistAgent        – pharmacology & target biology
    ToxicologistAgent     – safety & toxicity
    MicrobiologyAgent     – pathogen biology & AMR mechanisms
    GenomicsAgent         – resistance genomics & mutation analysis
    CheminformaticsAgent  – QSAR, property prediction, fingerprints
    ClinicalAgent         – translational medicine & trial design
"""

import json
import logging
from typing import Any, Dict, List, Optional

from .llm_provider import LLMProvider, FallbackProvider

logger = logging.getLogger('khukuri')


# ---------------------------------------------------------------------------
# Base
# ---------------------------------------------------------------------------

class DomainAgent:
    """
    Base class for all domain-specialist agents.

    Args:
        provider: An LLMProvider instance.  Defaults to FallbackProvider.
    """

    role: str = "Domain Specialist"
    expertise: str = "general drug discovery"
    _system_template: str = (
        "You are a {role} with deep expertise in {expertise}. "
        "Respond only with valid JSON unless asked otherwise."
    )

    def __init__(self, provider: Optional[LLMProvider] = None):
        self.provider = provider or FallbackProvider()
        self.memory: List[Dict[str, Any]] = []          # conversation history
        self._system = self._system_template.format(
            role=self.role, expertise=self.expertise
        )

    # ------------------------------------------------------------------
    # Primary interface
    # ------------------------------------------------------------------

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
        """
        Respond to another agent's view on `topic`.
        Enables direct peer-to-peer debate without routing through the PI.
        """
        prompt = (
            f"Topic: {topic}\n\n"
            f"A colleague said:\n{other_view}\n\n"
            f"As a {self.role}, provide your counterpoint or agreement. "
            "Return JSON: {agrees: bool, points: [...], evidence: [...], "
            "revised_position: str}"
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

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# Specialist agents
# ---------------------------------------------------------------------------

class ChemistAgent(DomainAgent):
    role = "Medicinal Chemist"
    expertise = (
        "small molecule drug design, synthetic accessibility, "
        "bioisosteres, SAR analysis, and Lipinski/Veber rules"
    )

    def analyze(self, data: Dict[str, Any], question: str) -> Dict[str, Any]:
        prompt = (
            f"As a Medicinal Chemist, analyse this molecule:\n"
            f"SMILES: {data.get('smiles', 'N/A')}\n"
            f"Properties: {json.dumps(data.get('properties', {}))}\n"
            f"Docking score: {data.get('docking_score', 'N/A')}\n\n"
            f"Question: {question}\n\n"
            "Return JSON: {drug_likeness, synthetic_accessibility, "
            "concerns, suggested_modifications, confidence}"
        )
        try:
            return self.provider.complete_json(
                [{"role": "user", "content": prompt}], system=self._system
            )
        except Exception:
            return self._fallback(question)


class BiologistAgent(DomainAgent):
    role = "Pharmacologist"
    expertise = (
        "target biology, mechanism of action, selectivity profiling, "
        "in vitro / in vivo pharmacology"
    )

    def analyze(self, data: Dict[str, Any], question: str) -> Dict[str, Any]:
        prompt = (
            f"As a Pharmacologist, assess this drug candidate:\n"
            f"Target: {data.get('target', 'N/A')}\n"
            f"Binding affinity: {data.get('binding_affinity', 'N/A')}\n"
            f"ADMET: {json.dumps(data.get('admet', {}))}\n\n"
            f"Question: {question}\n\n"
            "Return JSON: {pharmacological_viability, selectivity_concerns, "
            "recommended_assays, concerns, confidence}"
        )
        try:
            return self.provider.complete_json(
                [{"role": "user", "content": prompt}], system=self._system
            )
        except Exception:
            return self._fallback(question)


class ToxicologistAgent(DomainAgent):
    role = "Toxicologist"
    expertise = (
        "drug safety, structural alerts, hERG liability, hepatotoxicity, "
        "genotoxicity, and regulatory toxicology"
    )

    def analyze(self, data: Dict[str, Any], question: str) -> Dict[str, Any]:
        prompt = (
            f"As a Toxicologist, assess safety for:\n"
            f"SMILES: {data.get('smiles', 'N/A')}\n"
            f"Toxicity alerts: {json.dumps(data.get('toxicity', {}))}\n"
            f"ADMET: {json.dumps(data.get('admet', {}))}\n\n"
            f"Question: {question}\n\n"
            "Return JSON: {safety_assessment, structural_alerts, "
            "herg_risk, hepatotoxicity_risk, recommended_studies, confidence}"
        )
        try:
            return self.provider.complete_json(
                [{"role": "user", "content": prompt}], system=self._system
            )
        except Exception:
            return self._fallback(question)


class MicrobiologyAgent(DomainAgent):
    role = "Microbiologist"
    expertise = (
        "pathogen biology, antimicrobial resistance mechanisms, "
        "MIC analysis, efflux pumps, and biofilm"
    )

    def analyze(self, data: Dict[str, Any], question: str) -> Dict[str, Any]:
        prompt = (
            f"As a Microbiologist specialising in AMR, analyse:\n"
            f"Pathogen: {data.get('pathogen', 'N/A')}\n"
            f"Strain: {data.get('strain_id', 'N/A')}\n"
            f"Resistance profile: {data.get('resistance_profile', 'N/A')}\n"
            f"MIC data: {data.get('mic_data', 'N/A')}\n\n"
            f"Question: {question}\n\n"
            "Return JSON: {resistance_mechanisms, efficacy_assessment, "
            "cross_resistance_risk, concerns, recommendations, confidence}"
        )
        try:
            return self.provider.complete_json(
                [{"role": "user", "content": prompt}], system=self._system
            )
        except Exception:
            return self._fallback(question)


class GenomicsAgent(DomainAgent):
    role = "Genomics Specialist"
    expertise = (
        "resistance genomics, mutation analysis, horizontal gene transfer, "
        "CRISPR, and evolutionary pressure modelling"
    )

    def analyze(self, data: Dict[str, Any], question: str) -> Dict[str, Any]:
        prompt = (
            f"As a Genomics Specialist, analyse:\n"
            f"Resistance genes: {data.get('resistance_genes', [])}\n"
            f"Mutations: {json.dumps(data.get('mutations', {}))}\n"
            f"Drug pressure: {data.get('drug_pressure', 'N/A')}\n\n"
            f"Question: {question}\n\n"
            "Return JSON: {resistance_genetic_basis, evolution_potential, "
            "high_risk_mutations, multi_target_recommendation, confidence}"
        )
        try:
            return self.provider.complete_json(
                [{"role": "user", "content": prompt}], system=self._system
            )
        except Exception:
            return self._fallback(question)


class CheminformaticsAgent(DomainAgent):
    role = "Cheminformatics Specialist"
    expertise = (
        "QSAR modelling, molecular fingerprints, chemical space analysis, "
        "scaffold hopping, and virtual screening"
    )

    def analyze(self, data: Dict[str, Any], question: str) -> Dict[str, Any]:
        prompt = (
            f"As a Cheminformatics Specialist, analyse:\n"
            f"SMILES: {data.get('smiles', 'N/A')}\n"
            f"Fingerprint similarity: {data.get('similarity', 'N/A')}\n"
            f"Chemical space: {json.dumps(data.get('chemical_space', {}))}\n\n"
            f"Question: {question}\n\n"
            "Return JSON: {novelty_score, scaffold_analysis, "
            "similar_known_drugs, qsar_prediction, confidence}"
        )
        try:
            return self.provider.complete_json(
                [{"role": "user", "content": prompt}], system=self._system
            )
        except Exception:
            return self._fallback(question)


class ClinicalAgent(DomainAgent):
    role = "Clinical Scientist"
    expertise = (
        "translational medicine, clinical trial design, PK/PD modelling, "
        "dosing strategy, and regulatory pathways"
    )

    def analyze(self, data: Dict[str, Any], question: str) -> Dict[str, Any]:
        prompt = (
            f"As a Clinical Scientist, assess translational potential:\n"
            f"Compound: {data.get('name', 'N/A')}\n"
            f"ADMET: {json.dumps(data.get('admet', {}))}\n"
            f"Preclinical data: {json.dumps(data.get('preclinical', {}))}\n\n"
            f"Question: {question}\n\n"
            "Return JSON: {translational_viability, pk_pd_assessment, "
            "dosing_strategy, regulatory_pathway, confidence}"
        )
        try:
            return self.provider.complete_json(
                [{"role": "user", "content": prompt}], system=self._system
            )
        except Exception:
            return self._fallback(question)


# ---------------------------------------------------------------------------
# Agent registry
# ---------------------------------------------------------------------------

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
