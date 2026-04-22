"""
Peer debate system for Khukuri agents.

Enables direct agent-to-agent communication without routing through the PI.
Supports PeerDebate, PanelDiscussion, and SocraticDialog.
"""

import json
import logging
from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from .domain_agents import DomainAgent
from .llm_provider import LLMProvider, FallbackProvider

logger = logging.getLogger('khukuri')


@dataclass
class Turn:
    agent_role: str
    content: str
    agrees: Optional[bool] = None
    key_points: List[str] = field(default_factory=list)


@dataclass
class DebateRecord:
    topic: str
    turns: List[Turn] = field(default_factory=list)
    consensus: Optional[str] = None
    key_disagreements: List[str] = field(default_factory=list)
    final_recommendation: str = ""
    rounds_completed: int = 0

    def to_dict(self) -> Dict[str, Any]:
        return {
            "topic": self.topic,
            "rounds": self.rounds_completed,
            "consensus": self.consensus,
            "key_disagreements": self.key_disagreements,
            "final_recommendation": self.final_recommendation,
            "turns": [
                {
                    "agent": t.agent_role,
                    "content": t.content[:300],
                    "agrees": t.agrees,
                }
                for t in self.turns
            ],
        }


class PeerDebate:
    """Structured debate between exactly two agents."""

    def __init__(
        self,
        agent_a: DomainAgent,
        agent_b: DomainAgent,
        moderator_provider: Optional[LLMProvider] = None,
    ):
        self.agent_a = agent_a
        self.agent_b = agent_b
        self.mod_provider = moderator_provider or agent_a.provider

    def run(
        self,
        topic: str,
        context: Dict[str, Any],
        rounds: int = 2,
    ) -> DebateRecord:
        record = DebateRecord(topic=topic)

        opening_a = self.agent_a.analyze(context, f"State your position on: {topic}")
        opening_b = self.agent_b.analyze(context, f"State your position on: {topic}")

        record.turns.append(Turn(
            agent_role=self.agent_a.role,
            content=json.dumps(opening_a),
            key_points=opening_a.get("recommendations", []),
        ))
        record.turns.append(Turn(
            agent_role=self.agent_b.role,
            content=json.dumps(opening_b),
            key_points=opening_b.get("recommendations", []),
        ))

        last_a = json.dumps(opening_a)
        last_b = json.dumps(opening_b)

        for r in range(rounds):
            response_a = self.agent_a.debate(last_b, topic)
            record.turns.append(Turn(
                agent_role=self.agent_a.role,
                content=json.dumps(response_a),
                agrees=response_a.get("agrees"),
            ))
            last_a = json.dumps(response_a)

            response_b = self.agent_b.debate(last_a, topic)
            record.turns.append(Turn(
                agent_role=self.agent_b.role,
                content=json.dumps(response_b),
                agrees=response_b.get("agrees"),
            ))
            last_b = json.dumps(response_b)

        record.rounds_completed = rounds
        synthesis = self._moderate(record, topic)
        record.consensus = synthesis.get("consensus")
        record.key_disagreements = synthesis.get("key_disagreements", [])
        record.final_recommendation = synthesis.get("recommendation", "")

        return record

    def _moderate(self, record: DebateRecord, topic: str) -> Dict[str, Any]:
        turns_text = "\n".join(
            f"{t.agent_role}: {t.content[:200]}" for t in record.turns
        )
        prompt = (
            f"Moderate this scientific debate on: {topic}\n\n"
            f"TRANSCRIPT:\n{turns_text}\n\n"
            "Return JSON: {consensus: str|null, key_disagreements: [], recommendation: str}"
        )
        try:
            return self.mod_provider.complete_json(
                [{"role": "user", "content": prompt}],
                system="You are an impartial scientific moderator.",
                temperature=0.2,
            )
        except Exception as exc:
            logger.warning(f"Moderation failed: {exc}")
            return {
                "consensus": "Debate inconclusive",
                "key_disagreements": [],
                "recommendation": "Further discussion needed",
            }


class PanelDiscussion:
    """Round-robin panel where N agents each speak and respond."""

    def __init__(
        self,
        agents: List[DomainAgent],
        moderator_provider: Optional[LLMProvider] = None,
    ):
        if len(agents) < 2:
            raise ValueError("Panel needs at least 2 agents.")
        self.agents = agents
        self.mod_provider = moderator_provider or agents[0].provider

    def run(
        self,
        topic: str,
        context: Dict[str, Any],
        max_rounds: int = 3,
    ) -> DebateRecord:
        record = DebateRecord(topic=topic)
        panel_memory: List[str] = []

        for round_num in range(max_rounds):
            for agent in self.agents:
                prompt_ctx = {
                    **context,
                    "panel_so_far": panel_memory[-6:],
                    "round": round_num + 1,
                }
                response = agent.analyze(
                    prompt_ctx,
                    f"{'State your position' if round_num == 0 else 'Respond to your colleagues'} on: {topic}",
                )
                turn = Turn(
                    agent_role=agent.role,
                    content=json.dumps(response),
                    key_points=response.get("recommendations", []),
                )
                record.turns.append(turn)
                panel_memory.append(f"{agent.role}: {json.dumps(response)[:200]}")

            consensus = self._check_consensus(record.turns[-len(self.agents):], topic)
            if consensus.get("reached"):
                record.consensus = consensus.get("statement")
                break

        record.rounds_completed = round_num + 1

        synthesis = self._synthesise(record, topic)
        record.consensus = record.consensus or synthesis.get("consensus")
        record.key_disagreements = synthesis.get("key_disagreements", [])
        record.final_recommendation = synthesis.get("recommendation", "")

        return record

    def _check_consensus(self, recent_turns: List[Turn], topic: str) -> Dict[str, Any]:
        """Quick check — has consensus emerged?"""
        snippets = [t.content[:150] for t in recent_turns]
        prompt = (
            f"Have these agents reached consensus on '{topic}'?\n"
            + "\n".join(snippets)
            + "\nReturn JSON: {reached: bool, statement: str}"
        )
        try:
            return self.mod_provider.complete_json(
                [{"role": "user", "content": prompt}],
                temperature=0.1,
            )
        except Exception:
            return {"reached": False}

    def _synthesise(self, record: DebateRecord, topic: str) -> Dict[str, Any]:
        transcript = "\n".join(
            f"{t.agent_role}: {t.content[:200]}" for t in record.turns[-12:]
        )
        prompt = (
            f"Synthesise this panel discussion on: {topic}\n\n{transcript}\n\n"
            "Return JSON: {consensus, key_disagreements: [], recommendation}"
        )
        try:
            return self.mod_provider.complete_json(
                [{"role": "user", "content": prompt}],
                system="You are a scientific moderator synthesising a panel.",
                temperature=0.2,
            )
        except Exception as exc:
            return {"consensus": None, "recommendation": f"Synthesis failed: {exc}"}


class SocraticDialog:
    """One agent probes another with targeted questions."""

    def __init__(
        self,
        questioner: DomainAgent,
        responder: DomainAgent,
        n_questions: int = 3,
    ):
        self.questioner = questioner
        self.responder = responder
        self.n_questions = n_questions

    def run(self, claim: str, context: Dict[str, Any]) -> DebateRecord:
        record = DebateRecord(topic=f"Socratic probe: {claim}")

        initial = self.responder.analyze(context, f"Explain your reasoning: {claim}")
        record.turns.append(Turn(
            agent_role=self.responder.role,
            content=json.dumps(initial),
        ))

        current_answer = json.dumps(initial)

        for i in range(self.n_questions):
            question_resp = self.questioner.analyze(
                {"claim": claim, "current_answer": current_answer, "context": context},
                f"Ask a sharp Socratic question (round {i+1}/{self.n_questions}).",
            )
            question_text = question_resp.get("assessment", json.dumps(question_resp))
            record.turns.append(Turn(
                agent_role=f"{self.questioner.role} (questioner)",
                content=question_text,
            ))

            answer = self.responder.analyze(
                {**context, "question": question_text, "prior_answer": current_answer},
                "Answer this probing question rigorously.",
            )
            current_answer = json.dumps(answer)
            record.turns.append(Turn(
                agent_role=self.responder.role,
                content=current_answer,
            ))

        record.rounds_completed = self.n_questions
        record.final_recommendation = (
            f"Refined position after Socratic probing: {current_answer[:300]}"
        )
        return record
