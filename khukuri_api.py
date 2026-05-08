"""
Khukuri API Server — local Flask backend for the GUI.

Run: python -m khukuri_api
Then open: http://localhost:5000

Endpoints:
    GET  /api/status          — system health + mode
    GET  /api/world           — world model summary
    GET  /api/sentinel        — drug surveillance dashboard
    GET  /api/compounds       — top compounds from world model
    GET  /api/hypotheses      — open hypotheses
    GET  /api/runs            — recent cycle history
    POST /api/run             — trigger a discovery cycle
    POST /api/mode            — switch supervised/lights_out
    GET  /api/health          — watchdog system health
    GET  /api/pincer/threats  — current threat breakdown
"""

import json
import os
import sys
import logging
from datetime import datetime
from pathlib import Path

from flask import Flask, jsonify, request
from flask_cors import CORS

# Add khukuri to path
sys.path.insert(0, str(Path(__file__).parent))
sys.path.insert(0, str(Path(__file__).parent.parent / "keybox"))

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("khukuri.api")

app = Flask(__name__)
CORS(app)

# ---------------------------------------------------------------------------
# Lazy-loaded Khukuri components
# ---------------------------------------------------------------------------

_components = {}

def get_wm():
    if "wm" not in _components:
        from src.world_model.manager import WorldModelManager
        _components["wm"] = WorldModelManager(project=PROJECT)
    return _components["wm"]

def get_sentinel():
    if "sentinel" not in _components:
        from src.watchdog.drug_sentinel import DrugResistanceSentinel
        from src.resistance.threat_fitness import ThreatAwareFitnessFunction
        s = DrugResistanceSentinel(project=PROJECT,
                                   fitness_fn=ThreatAwareFitnessFunction())
        s.load_mrsa_panel()
        _components["sentinel"] = s
    return _components["sentinel"]

def get_health_watchdog():
    if "health" not in _components:
        from src.watchdog.system_health import SystemHealthWatchdog
        _components["health"] = SystemHealthWatchdog()
    return _components["health"]

def get_numen():
    if "numen" not in _components:
        from src.integrations.numen_index import KhukuriNumen
        _components["numen"] = KhukuriNumen(project=PROJECT)
    return _components["numen"]

def get_controller():
    if "controller" not in _components:
        from src.autonomous.controller import AutonomousLoopController, Mode
        from src.agents.llm_provider import create_provider
        _components["controller"] = AutonomousLoopController(
            mode=Mode.SUPERVISED,
            project=PROJECT,
            on_decision_needed=lambda dp: "approve",
        )
    return _components["controller"]

PROJECT = os.environ.get("KHUKURI_PROJECT", "default")

# ---------------------------------------------------------------------------
# Routes
# ---------------------------------------------------------------------------

@app.route("/api/status")
def api_status():
    try:
        ctrl = get_controller()
        wm = get_wm()
        s = wm.get_summary()
        return jsonify({
            "status": "running",
            "mode": ctrl.mode.value,
            "project": PROJECT,
            "timestamp": datetime.now().isoformat(),
            "runs_completed": s.get("runs_completed", 0),
            "compounds": s.get("state", {}).get("compounds", 0),
            "targets": s.get("state", {}).get("targets", 0),
            "hypotheses": s.get("state", {}).get("hypotheses", {}).get("total", 0),
            "consecutive_stagnant": ctrl._consecutive_stagnant,
        })
    except Exception as e:
        return jsonify({"status": "error", "error": str(e)}), 500


@app.route("/api/world")
def api_world():
    try:
        wm = get_wm()
        summary = wm.get_summary()
        top = wm.get_top_compounds(n=10)
        hyps = wm.get_actionable_hypotheses()
        return jsonify({
            "summary": summary,
            "top_compounds": top,
            "hypotheses": hyps[:5],
        })
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/sentinel")
def api_sentinel():
    try:
        sentinel = get_sentinel()
        # Load threats
        from src.resistance.pincer_engine import MutationSpaceMapper
        from src.resistance.hgt_mapper import HGTMapper
        from src.world_model.knowledge_graph import KnowledgeGraph
        mapper = MutationSpaceMapper()
        threats = mapper.map_binding_pocket("AMILVCFYWHDEKRST", list(range(12)))
        kg = KnowledgeGraph()
        hgt = HGTMapper(kg)
        cluster = hgt.map_threats("S. aureus", selective_pressure="high")
        sentinel.load_threats(threats, cluster.to_viable_threats())
        alerts = sentinel.run_surveillance()
        status = sentinel.get_status()
        critical = [d.name for d in sentinel.get_critical_drugs()]
        return jsonify({
            "status": status,
            "alerts": [a.to_dict() for a in alerts],
            "critical_drugs": critical,
            "threat_counts": {
                "mutational": len(threats),
                "hgt": len(cluster.alleles),
                "total": len(threats) + len(cluster.alleles),
            },
        })
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/compounds")
def api_compounds():
    try:
        wm = get_wm()
        metric = request.args.get("metric", "pincer_score")
        n = int(request.args.get("n", 20))
        compounds = wm.get_top_compounds(n=n, metric=metric)
        return jsonify({"compounds": compounds, "metric": metric, "count": len(compounds)})
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/hypotheses")
def api_hypotheses():
    try:
        wm = get_wm()
        hyps = wm.get_actionable_hypotheses(min_confidence=0.5)
        return jsonify({"hypotheses": hyps, "count": len(hyps)})
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/health")
def api_health():
    try:
        watchdog = get_health_watchdog()
        report = watchdog.check()
        return jsonify(report.to_dict())
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/runs")
def api_runs():
    try:
        state_dir = Path.home() / ".khukuri" / "memory" / PROJECT / "autonomous"
        cycles = []
        if state_dir.exists():
            for f in sorted(state_dir.glob("cycle_*.json"))[-10:]:
                try:
                    cycles.append(json.loads(f.read_text()))
                except Exception:
                    pass
        return jsonify({"runs": list(reversed(cycles)), "count": len(cycles)})
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/mode", methods=["POST"])
def api_set_mode():
    try:
        from src.autonomous.controller import Mode
        data = request.json
        mode_str = data.get("mode", "supervised")
        mode = Mode.SUPERVISED if mode_str == "supervised" else Mode.LIGHTS_OUT
        get_controller().set_mode(mode)
        return jsonify({"mode": mode.value, "ok": True})
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/run", methods=["POST"])
def api_trigger_run():
    """Trigger a discovery cycle (runs in background thread)."""
    try:
        import threading
        data = request.json or {}
        pathogen = data.get("pathogen", "S. aureus MRSA")
        mode = data.get("mode", "supervised")

        def run_bg():
            try:
                from src.workflows.amr_discovery import AMRDiscoveryWorkflow
                wf = AMRDiscoveryWorkflow(project=PROJECT)
                wf.run_discovery(pathogen=pathogen, n_iterations=2)
                logger.info(f"Background run complete for {pathogen}")
            except Exception as ex:
                logger.error(f"Background run failed: {ex}")

        thread = threading.Thread(target=run_bg, daemon=True)
        thread.start()
        return jsonify({
            "ok": True,
            "message": f"Discovery run started for {pathogen}",
            "pathogen": pathogen,
        })
    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/chat", methods=["POST"])
def api_chat():
    """LLM chat endpoint — routes to Anthropic or OpenAI."""
    try:
        data = request.json or {}
        model = data.get("model", "claude-sonnet-4-20250514")
        system = data.get("system", "")
        messages = data.get("messages", [])
        api_key = data.get("api_key") or None

        if model.startswith("claude"):
            import anthropic
            key = api_key or os.environ.get("ANTHROPIC_API_KEY")
            if not key:
                return jsonify({"error": "No Anthropic API key. Set ANTHROPIC_API_KEY or enter in GUI."}), 400
            client = anthropic.Anthropic(api_key=key)
            resp = client.messages.create(
                model=model,
                max_tokens=1000,
                system=system,
                messages=messages[-10:],  # last 10 turns
            )
            return jsonify({"reply": resp.content[0].text})

        else:
            import openai
            key = api_key or os.environ.get("OPENAI_API_KEY")
            if not key:
                return jsonify({"error": "No OpenAI API key. Set OPENAI_API_KEY or enter in GUI."}), 400
            client = openai.OpenAI(api_key=key)
            full_messages = [{"role": "system", "content": system}] + messages[-10:]
            resp = client.chat.completions.create(
                model=model, messages=full_messages, max_tokens=1000
            )
            return jsonify({"reply": resp.choices[0].message.content})

    except Exception as e:
        return jsonify({"error": str(e)}), 500


@app.route("/api/search")
def api_search():
    """Literature search via Numen index + optional PubMed."""
    try:
        query = request.args.get("q", "")
        if not query:
            return jsonify({"results": []})

        numen = get_numen()
        # Search local index
        local = numen.search_literature(query, top_k=5)
        results = []
        for r in local:
            meta = r.get("meta", {})
            if meta.get("title"):
                results.append({
                    "title": meta.get("title", ""),
                    "snippet": meta.get("title", "")[:200],
                    "url": meta.get("url", ""),
                    "score": r.get("score", 0),
                })

        # If few local results, try PubMed
        if len(results) < 3:
            n_indexed = numen.literature.index_pubmed_query(query, max_results=10)
            if n_indexed > 0:
                fresh = numen.search_literature(query, top_k=5)
                for r in fresh:
                    meta = r.get("meta", {})
                    if meta.get("title") and meta not in [x.get("meta") for x in local]:
                        results.append({
                            "title": meta.get("title", ""),
                            "snippet": meta.get("title", "")[:200],
                            "url": meta.get("url", ""),
                            "score": r.get("score", 0),
                        })

        return jsonify({"results": results[:8], "query": query})
    except Exception as e:
        return jsonify({"error": str(e), "results": []}), 500
def api_threats():
    try:
        from src.resistance.pincer_engine import MutationSpaceMapper
        from src.resistance.hgt_mapper import HGTMapper
        from src.world_model.knowledge_graph import KnowledgeGraph
        mapper = MutationSpaceMapper()
        threats = mapper.map_binding_pocket("AMILVCFYWHDEKRST", list(range(12)))
        kg = KnowledgeGraph()
        hgt = HGTMapper(kg)
        cluster = hgt.map_threats("S. aureus", selective_pressure="high")
        top_hgt = hgt.get_high_risk_genes(top_n=5)
        return jsonify({
            "mutational_count": len(threats),
            "hgt_count": len(cluster.alleles),
            "total": len(threats) + len(cluster.alleles),
            "top_hgt_genes": top_hgt,
            "network_stats": cluster.network_stats,
        })
    except Exception as e:
        return jsonify({"error": str(e)}), 500


# ---------------------------------------------------------------------------
# Serve static frontend
# ---------------------------------------------------------------------------

@app.route("/")
def index():
    frontend = Path(__file__).parent / "gui" / "index.html"
    if frontend.exists():
        return frontend.read_text()
    return "<h1>Khukuri API running. Frontend not found at gui/index.html</h1>"


if __name__ == "__main__":
    print(f"""
╔══════════════════════════════════════╗
║   Khukuri API Server                 ║
║   http://localhost:5000              ║
║   Project: {PROJECT:<26}║
╚══════════════════════════════════════╝
""")
    app.run(debug=False, host="0.0.0.0", port=5000)
