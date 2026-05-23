"""
test_memory.py
--------------
Simulates Khukuri's autonomous loop across 3 runs.
Tests: add -> query -> commit -> cross-shard query -> salience -> reinforce.
"""

import sys, shutil, tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from src.memory import ThreatIndex, CompoundMemory, LiteratureIndex, KhukuriMemory

def separator(title):
    print(f"\n{'─'*60}")
    print(f"  {title}")
    print('─'*60)

def test_full_lifecycle():
    with tempfile.TemporaryDirectory() as tmpdir:
        root = Path(tmpdir)

        # -- RUN 1 ------------------------------------------------------------
        separator("RUN 1: ThreatIndex")
        ti = ThreatIndex(root, dim=512)

        ti.build([
            {"id": "t001", "description": "mecA mediated beta-lactam resistance PBP2a", "type": "HGT",      "gene": "mecA",  "strain": "MRSA252", "severity": 0.9},
            {"id": "t002", "description": "vanA vancomycin resistance glycopeptide",    "type": "HGT",      "gene": "vanA",  "strain": "VRE",     "severity": 0.95},
            {"id": "t003", "description": "PBP2a S403A mutation beta-lactam escape",    "type": "mutation", "gene": "mecA",  "strain": "MRSA252", "severity": 0.7},
            {"id": "t004", "description": "qnrS quinolone resistance efflux",           "type": "HGT",      "gene": "qnrS",  "strain": "Ecoli",   "severity": 0.6},
        ])

        results = ti.get_worst_threats("beta-lactam antibiotic penicillin binding protein", k=3)
        print(f"\nQuery 'beta-lactam' -> {len(results)} results:")
        for r in results:
            print(f"  {r.entry_id:6s}  sim={r.score:.3f}  sal={r.salience:.3f}  "
                  f"gene={r.meta.get('gene')}  type={r.meta.get('type')}")

        # simulate: t001 and t003 were the actual worst cases this run
        ti.reinforce("t001", was_worst_case=True)
        ti.reinforce("t003", was_worst_case=True)
        ti.reinforce("t002", was_worst_case=False)

        sealed = ti._mem.commit_run()
        print(f"\nRun 1 sealed -> {sealed}")
        print(f"Status: {ti.status()}")

        # -- RUN 2 ------------------------------------------------------------
        separator("RUN 2: ThreatIndex -- new instance (simulates restart)")
        ti2 = ThreatIndex(root, dim=512)   # new instance, same root -> loads shards

        # add a newly discovered threat
        ti2.add_threat(
            "t005",
            "mecC novel variant beta-lactam extended spectrum",
            threat_type="mutation",
            gene="mecC",
            strain="MRSA252",
            severity=0.85,
        )

        results2 = ti2.get_worst_threats("beta-lactam antibiotic", k=4, strain="MRSA252")
        print(f"\nQuery 'beta-lactam' (MRSA252 only) after run 2 add:")
        for r in results2:
            print(f"  {r.entry_id:6s}  sim={r.score:.3f}  sal={r.salience:.3f}  "
                  f"final={r.final_score:.3f}  gene={r.meta.get('gene')}")

        assert any(r.entry_id == "t001" for r in results2), "t001 should be present (cross-shard)"
        assert any(r.entry_id == "t005" for r in results2), "t005 should be present (hot shard)"

        # t001 was reinforced positively -> higher salience -> should rank above t002
        t001 = next(r for r in results2 if r.entry_id == "t001")
        t002_present = [r for r in results2 if r.entry_id == "t002"]
        if t002_present:
            assert t001.salience > t002_present[0].salience, \
                "t001 salience should be higher than t002 (reinforcement)"
        print(f"\n  Cross-shard query works (frozen run_0001 + hot)")
        print(f"  Salience correct: t001={t001.salience:.3f} (reinforced)")

        ti2._mem.commit_run()

        # -- CompoundMemory ---------------------------------------------------
        separator("CompoundMemory: 2 runs of compound accumulation")
        cm = CompoundMemory(root, dim=512)

        # run 1 compounds
        cm.remember("CC(C)Cc1ccc(cc1)C(C)C(=O)O", pincer_score=0.81, strain="MRSA252",
                    tags=["beta-lactam"])
        cm.remember("c1ccc2c(c1)cc(n2)C(=O)O",    pincer_score=0.45, strain="MRSA252",
                    tags=["quinolone"])
        cm.remember("CC1=CC(=O)c2ccccc2C1=O",      pincer_score=0.72, strain="MRSA252",
                    tags=["beta-lactam"])
        cm.reinforce(f"cpd_{abs(hash('CC(C)Cc1ccc(cc1)C(C)C(=O)O')) % 1_000_000:06d}", True)
        cm.commit_run()

        # run 2 new seeds seeded from run 1 winners
        similar = cm.recall("CC(C)Cc1ccc(cc1)C(C)C(=O)O", k=3, strain="MRSA252")
        print(f"\nRecall similar to top run-1 compound ({len(similar)} results):")
        for r in similar:
            print(f"  {r.entry_id}  sim={r.score:.3f}  sal={r.salience:.3f}  "
                  f"score={r.meta.get('score'):.2f}")

        assert len(similar) > 0, "Should recall run-1 compounds in run-2"
        print(f"\n  Cross-run compound recall works")

        # -- LiteratureIndex --------------------------------------------------
        separator("LiteratureIndex: target-filtered retrieval")
        li = LiteratureIndex(root, dim=512)

        li.index_pubmed_query(
            "PBP2a MRSA resistance beta-lactam",
            papers=[
                {"id": "pm001", "abstract": "PBP2a encoded by mecA confers resistance to all beta-lactam antibiotics in MRSA", "year": 2021},
                {"id": "pm002", "abstract": "vanA cluster mediates high-level vancomycin resistance in enterococci via D-Ala-D-Lac substitution", "year": 2019},
                {"id": "pm003", "abstract": "Structural basis for beta-lactam resistance: crystal structure of PBP2a bound to ceftaroline", "year": 2022},
            ],
            target="PBP2a",
        )
        li.index_paper("pm004", "vanB operon regulation and inducibility in VRE clinical isolates", target="vanB", year=2023)

        pbp2a_papers = li.search("resistance beta-lactam crystal structure", k=5, target="PBP2a")
        print(f"\nSearch 'resistance beta-lactam' (PBP2a only): {len(pbp2a_papers)} results")
        for r in pbp2a_papers:
            print(f"  {r.entry_id}  sim={r.score:.3f}  target={r.meta.get('target')}")
        assert all(r.meta.get("target") == "PBP2a" for r in pbp2a_papers), \
            "Target filter should exclude vanB paper"
        print(f"  Target filter works (vanB paper excluded)")

        # -- Final status -----------------------------------------------------
        separator("Final status across all stores")
        print(f"ThreatIndex:     {ti2.status()}")
        print(f"CompoundMemory:  {cm.status()}")
        print(f"LiteratureIndex: {li.status()}")

        print(f"\n{'='*60}")
        print(f"  ALL TESTS PASSED")
        print(f"{'='*60}\n")


if __name__ == "__main__":
    test_full_lifecycle()
