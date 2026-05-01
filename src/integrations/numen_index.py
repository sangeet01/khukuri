"""
NumenIndex — Numen retrieval integrated into Khukuri.

Three use cases, one unified index class:

    1. ThreatIndex    — index viable threat cluster, retrieve worst-case
                        threats per drug (replaces O(n) loop in PINCER)

    2. LiteratureIndex — index PubMed/local abstracts for target discovery,
                         query with biological questions, feed NetworkAnalyzer

    3. CompoundMemory  — index world model compound history, retrieve similar
                         previously-evaluated compounds + their outcomes,
                         warm the LearningLoop

All three share the same NumenRetriever core (n-gram hashing, CRC32,
log-saturation, cosine similarity) extracted from limitnumen.

No training. No API. Arbitrary dimension. Works offline.
"""

import json
import logging
import zlib
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger('khukuri')


# ---------------------------------------------------------------------------
# Core: NumenRetriever (extracted from limitnumen/limit/numen.ipynb)
# ---------------------------------------------------------------------------

class NumenRetriever:
    """
    Character N-Gram Hashing retriever.

    Encodes text into a dense vector of arbitrary dimension using
    3/4/5-gram CRC32 hashing with log-saturation and cosine similarity.

    Beats BM25 on the LIMIT benchmark at 32768d (93.9% Recall@100).
    Works at any dimension — 512d is enough for Khukuri's use cases.
    """

    def __init__(self, dim: int = 512):
        self.dim = dim
        self._index: Dict[str, np.ndarray] = {}   # id -> unit vector
        self._meta: Dict[str, Any] = {}            # id -> metadata

    def encode(self, text: str) -> np.ndarray:
        """Encode text to a unit vector via n-gram hashing."""
        text = text.lower()
        words = text.split()
        if not words:
            return np.zeros(self.dim, dtype=np.float32)

        vector = np.zeros(self.dim, dtype=np.float32)
        for word in words:
            marked = f"^{word}$"
            for n in [3, 4, 5]:
                if len(marked) >= n:
                    for i in range(len(marked) - n + 1):
                        gram = marked[i:i+n]
                        idx = zlib.crc32(gram.encode()) % self.dim
                        weight = 10.0 if n >= 5 else (5.0 if n == 4 else 1.0)
                        vector[idx] += weight

        vector = np.log1p(vector)
        norm = np.linalg.norm(vector)
        return vector / norm if norm > 0 else vector

    def add(self, item_id: str, text: str, meta: Optional[Dict] = None):
        """Add one item to the index."""
        self._index[item_id] = self.encode(text)
        self._meta[item_id] = meta or {}

    def add_batch(self, items: Dict[str, str],
                  metas: Optional[Dict[str, Dict]] = None):
        """Add multiple items. items = {id: text}"""
        metas = metas or {}
        for item_id, text in items.items():
            self.add(item_id, text, metas.get(item_id))

    def search(self, query: str, top_k: int = 5) -> List[Dict]:
        """
        Search the index. Returns top_k results sorted by cosine similarity.
        Each result: {id, score, meta}
        """
        if not self._index:
            return []

        qvec = self.encode(query)
        ids = list(self._index.keys())
        matrix = np.stack([self._index[i] for i in ids])  # (N, dim)
        scores = matrix @ qvec                              # cosine sim

        top_idx = np.argsort(scores)[::-1][:top_k]
        return [
            {
                "id": ids[i],
                "score": float(scores[i]),
                "meta": self._meta[ids[i]],
            }
            for i in top_idx
        ]

    def search_vector(self, vec: np.ndarray, top_k: int = 5) -> List[Dict]:
        """Search by vector instead of text (for threat embeddings)."""
        if not self._index:
            return []

        # Align dimensions
        if len(vec) != self.dim:
            padded = np.zeros(self.dim, dtype=np.float32)
            n = min(len(vec), self.dim)
            padded[:n] = vec[:n]
            vec = padded
        norm = np.linalg.norm(vec)
        vec = vec / norm if norm > 0 else vec

        ids = list(self._index.keys())
        matrix = np.stack([self._index[i] for i in ids])
        scores = matrix @ vec

        top_idx = np.argsort(scores)[::-1][:top_k]
        return [
            {
                "id": ids[i],
                "score": float(scores[i]),
                "meta": self._meta[ids[i]],
            }
            for i in top_idx
        ]

    def __len__(self):
        return len(self._index)

    def save(self, path: str):
        """Persist index to disk."""
        data = {
            "dim": self.dim,
            "ids": list(self._index.keys()),
            "vectors": np.stack(list(self._index.values())).tolist(),
            "meta": self._meta,
        }
        Path(path).write_text(json.dumps(data))

    def load(self, path: str):
        """Restore index from disk."""
        data = json.loads(Path(path).read_text())
        self.dim = data["dim"]
        vecs = np.array(data["vectors"], dtype=np.float32)
        for i, item_id in enumerate(data["ids"]):
            self._index[item_id] = vecs[i]
            self._meta[item_id] = data["meta"].get(item_id, {})


# ---------------------------------------------------------------------------
# Use Case 1: ThreatIndex
# ---------------------------------------------------------------------------

class ThreatIndex:
    """
    Numen index over PINCER's viable threat cluster.

    Instead of evaluating every drug against all N threats (O(N) per drug),
    retrieve the k most dangerous threats per drug via cosine similarity
    then evaluate only those. Reduces fitness evaluation cost dramatically
    for large threat clusters (400+ threats).

    Also enables the README's architecture: 'LimitNumen Vector Index'
    on the Red Team side, producing 'Apex Skeleton Key' candidates.

    Usage:
        ti = ThreatIndex(dim=256)
        ti.build(viable_threats)                    # index the threat cluster
        worst = ti.get_worst_threats(drug_smiles, k=20)  # retrieve per drug
        score = fitness_fn(drug_smiles, worst)      # evaluate only top-k
    """

    def __init__(self, dim: int = 256):
        self._retriever = NumenRetriever(dim=dim)
        self._threats: Dict[str, Any] = {}   # id -> threat object

    def build(self, threats: List) -> int:
        """
        Index a list of ViableMutant or ResistanceGeneAllele threats.
        Returns number of threats indexed.
        """
        self._retriever = NumenRetriever(dim=self._retriever.dim)
        self._threats = {}

        for i, threat in enumerate(threats):
            threat_id = f"threat_{i}"
            text = self._threat_to_text(threat)
            meta = {
                "fitness": getattr(threat, "fitness", 0.5),
                "fitness_cost": getattr(threat, "fitness_cost",
                                        getattr(threat, "fitness_cost", 0.0)),
                "index": i,
            }
            self._retriever.add(threat_id, text, meta)
            self._threats[threat_id] = threat

        logger.info(f"ThreatIndex: indexed {len(threats)} threats (dim={self._retriever.dim})")
        return len(threats)

    def get_worst_threats(self, smiles: str, k: int = 20) -> List:
        """
        Retrieve the k threats most relevant to this drug SMILES.
        Returns actual threat objects (ViableMutant / ResistanceGeneAllele).
        """
        if not self._threats:
            return []

        results = self._retriever.search(smiles, top_k=min(k, len(self._threats)))

        # Sort by threat fitness (most dangerous first) among retrieved
        retrieved = [self._threats[r["id"]] for r in results]
        retrieved.sort(
            key=lambda t: getattr(t, "fitness", 0.5),
            reverse=True,
        )
        return retrieved

    def get_worst_by_vector(self, drug_vec: np.ndarray, k: int = 20) -> List:
        """Retrieve worst threats by drug embedding vector."""
        results = self._retriever.search_vector(drug_vec, top_k=k)
        return [self._threats[r["id"]] for r in results]

    def __len__(self):
        return len(self._threats)

    @staticmethod
    def _threat_to_text(threat) -> str:
        """Convert a threat object to a searchable text representation."""
        parts = []
        mutations = getattr(threat, "mutations", [])
        for m in mutations:
            wt = getattr(m, "wild_type", "")
            pos = getattr(m, "position", "")
            mut = getattr(m, "mutant", "")
            parts.append(f"{wt}{pos}{mut}")

        # HGT allele fields
        gene_id = getattr(threat, "gene_id", "")
        allele = getattr(threat, "allele", "")
        res_class = getattr(threat, "resistance_class", "")
        mge = getattr(threat, "mge_type", "")

        text = " ".join(filter(None, [
            gene_id, allele, res_class, mge,
            " ".join(parts),
        ]))
        return text or "unknown_threat"


# ---------------------------------------------------------------------------
# Use Case 2: LiteratureIndex
# ---------------------------------------------------------------------------

class LiteratureIndex:
    """
    Numen index for scientific literature mining.

    Indexes PubMed abstracts, local paper collections, or any text corpus.
    Powers Khukuri's target discovery with semantic retrieval — no API,
    no embeddings model, no internet required after initial indexing.

    Feeds directly into NetworkAnalyzer's literature_mining stub.

    Usage:
        li = LiteratureIndex()
        li.index_pubmed_query("PBP2a MRSA resistance")    # live PubMed fetch
        li.index_texts({"paper1": abstract1, ...})        # local corpus
        results = li.search("mecA horizontal transfer MRSA")
        targets = li.extract_targets(results)
    """

    # Known AMR-relevant terms for target extraction
    TARGET_KEYWORDS = [
        "pbp2a", "pbp2", "meca", "mecA", "katg", "inha", "gyra", "parc",
        "vanA", "vana", "ddl", "ftsl", "ftsz", "mura", "murb", "murc",
        "murd", "mure", "murf", "baiJ", "fabI", "fabH", "accD",
        "norA", "norB", "mepA", "lmrS", "arlRS",
    ]

    def __init__(self, dim: int = 1024, cache_dir: Optional[str] = None):
        self._retriever = NumenRetriever(dim=dim)
        self.cache_dir = Path(cache_dir or Path.home() / ".khukuri" / "literature")
        self.cache_dir.mkdir(parents=True, exist_ok=True)
        self._paper_count = 0

    def index_texts(self, texts: Dict[str, str],
                    metas: Optional[Dict[str, Dict]] = None):
        """Index a dict of {id: abstract_text}."""
        self._retriever.add_batch(texts, metas)
        self._paper_count += len(texts)
        logger.info(f"LiteratureIndex: {len(texts)} texts indexed "
                    f"(total={self._paper_count})")

    def index_pubmed_query(
        self,
        query: str,
        max_results: int = 50,
    ) -> int:
        """
        Fetch PubMed abstracts for query and index them.
        Uses NCBI E-utilities (no API key needed for <3 req/s).
        Returns number of abstracts indexed.
        """
        try:
            import urllib.request, urllib.parse, xml.etree.ElementTree as ET

            # 1. Search
            search_url = (
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?"
                + urllib.parse.urlencode({
                    "db": "pubmed", "term": query,
                    "retmax": max_results, "retmode": "json",
                })
            )
            with urllib.request.urlopen(search_url, timeout=15) as r:
                ids = json.loads(r.read())["esearchresult"]["idlist"]

            if not ids:
                logger.warning(f"LiteratureIndex: no PubMed results for '{query}'")
                return 0

            # 2. Fetch abstracts
            fetch_url = (
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?"
                + urllib.parse.urlencode({
                    "db": "pubmed", "id": ",".join(ids),
                    "rettype": "abstract", "retmode": "xml",
                })
            )
            with urllib.request.urlopen(fetch_url, timeout=30) as r:
                xml_data = r.read()

            # 3. Parse and index
            root = ET.fromstring(xml_data)
            texts = {}
            metas = {}
            for article in root.findall(".//PubmedArticle"):
                pmid_el = article.find(".//PMID")
                title_el = article.find(".//ArticleTitle")
                abstract_el = article.find(".//AbstractText")
                pmid = pmid_el.text if pmid_el is not None else "unknown"
                title = title_el.text if title_el is not None else ""
                abstract = abstract_el.text if abstract_el is not None else ""
                if abstract:
                    texts[f"pmid_{pmid}"] = f"{title} {abstract}"
                    metas[f"pmid_{pmid}"] = {
                        "pmid": pmid, "title": title,
                        "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
                    }

            self.index_texts(texts, metas)
            return len(texts)

        except Exception as exc:
            logger.warning(f"LiteratureIndex: PubMed fetch failed: {exc}")
            return 0

    def search(self, query: str, top_k: int = 10) -> List[Dict]:
        """Semantic search over indexed literature."""
        return self._retriever.search(query, top_k=top_k)

    def extract_targets(self, results: List[Dict]) -> List[str]:
        """
        Extract candidate drug targets from search results.
        Looks for known AMR target keywords in retrieved text.
        """
        found = set()
        for r in results:
            meta = r.get("meta", {})
            text = (meta.get("title", "") + " " +
                    str(meta.get("abstract", ""))).lower()
            for kw in self.TARGET_KEYWORDS:
                if kw.lower() in text:
                    found.add(kw)
        return list(found)

    def mine_for_targets(
        self,
        pathogen: str,
        resistance_genes: Optional[List[str]] = None,
        max_papers: int = 30,
    ) -> Dict[str, Any]:
        """
        Full literature mining pipeline for a pathogen.
        Fetches papers, indexes them, extracts targets.
        Returns dict ready to feed into NetworkAnalyzer.
        """
        queries = [
            f"{pathogen} drug target",
            f"{pathogen} resistance mechanism",
            f"{pathogen} essential gene antibiotic",
        ]
        if resistance_genes:
            queries.append(f"{' '.join(resistance_genes[:3])} resistance {pathogen}")

        n_indexed = 0
        for q in queries:
            n_indexed += self.index_pubmed_query(
                q, max_results=max_papers // len(queries)
            )

        results = self.search(
            f"{pathogen} drug target essential gene",
            top_k=20,
        )
        targets = self.extract_targets(results)

        return {
            "pathogen": pathogen,
            "papers_indexed": n_indexed,
            "candidate_targets": targets,
            "top_papers": [r["meta"] for r in results[:5]],
            "search_queries": queries,
        }

    def save(self, filename: Optional[str] = None):
        path = self.cache_dir / (filename or "literature_index.json")
        self._retriever.save(str(path))
        logger.info(f"LiteratureIndex: saved to {path}")

    def load(self, filename: Optional[str] = None):
        path = self.cache_dir / (filename or "literature_index.json")
        if path.exists():
            self._retriever.load(str(path))
            logger.info(f"LiteratureIndex: loaded {len(self._retriever)} papers")


# ---------------------------------------------------------------------------
# Use Case 3: CompoundMemory
# ---------------------------------------------------------------------------

class CompoundMemory:
    """
    Numen index over Khukuri's compound history.

    Every evaluated compound is indexed by SMILES.
    When a new candidate appears, retrieve the most similar
    previously-evaluated compounds and their outcomes —
    warm the LearningLoop without re-evaluation.

    Also enables the world model to answer:
        "Have we seen something like this before?"
        "What happened when we tried similar scaffolds?"

    Usage:
        cm = CompoundMemory()
        cm.remember("c1ccccc1", pincer_score=0.72, admet={...})
        similar = cm.recall("c1ccncc1", k=5)
        # returns [{smiles, pincer_score, admet, similarity}, ...]
    """

    def __init__(self, dim: int = 512):
        self._retriever = NumenRetriever(dim=dim)
        self._outcomes: Dict[str, Dict] = {}

    def remember(
        self,
        smiles: str,
        compound_id: Optional[str] = None,
        pincer_score: Optional[float] = None,
        docking_score: Optional[float] = None,
        admet: Optional[Dict] = None,
        script_str: Optional[str] = None,
    ):
        """Store a compound and its evaluation outcomes."""
        cid = compound_id or f"comp_{len(self._outcomes)}"
        outcome = {
            "smiles": smiles,
            "compound_id": cid,
            "pincer_score": pincer_score,
            "docking_score": docking_score,
            "admet": admet or {},
            "script": script_str,
        }
        # Index by SMILES — Numen naturally handles chemical substructure
        # similarity because fragment n-grams overlap for similar scaffolds
        text = smiles
        if script_str:
            text = f"{smiles} {script_str}"   # SCRIPT adds structural tokens
        self._retriever.add(cid, text, outcome)
        self._outcomes[cid] = outcome

    def remember_batch(self, compounds: List[Dict]):
        """Bulk remember. Each dict needs at least 'smiles'."""
        for c in compounds:
            self.remember(**{k: v for k, v in c.items()
                            if k in ['smiles','compound_id','pincer_score',
                                     'docking_score','admet','script_str']})

    def recall(self, smiles: str, k: int = 5) -> List[Dict]:
        """
        Retrieve k most similar previously-evaluated compounds.
        Returns list of outcome dicts with added 'similarity' field.
        """
        results = self._retriever.search(smiles, top_k=k)
        recalled = []
        for r in results:
            outcome = dict(self._outcomes.get(r["id"], {}))
            outcome["similarity"] = r["score"]
            recalled.append(outcome)
        return recalled

    def get_best(self, metric: str = "pincer_score", top_k: int = 10) -> List[Dict]:
        """Return top compounds by a stored metric."""
        scored = [
            v for v in self._outcomes.values()
            if v.get(metric) is not None
        ]
        scored.sort(key=lambda x: x.get(metric, 0.0), reverse=True)
        return scored[:top_k]

    def warm_learning_loop(self, learning_loop) -> int:
        """
        Feed stored compound outcomes into a LearningLoop instance.
        Returns number of observations added.
        """
        added = 0
        for cid, outcome in self._outcomes.items():
            smiles = outcome.get("smiles", "")
            if not smiles:
                continue
            try:
                from rdkit import Chem
                from rdkit.Chem import Descriptors
                mol = Chem.MolFromSmiles(smiles)
                if mol:
                    features = {
                        "mw": Descriptors.MolWt(mol),
                        "logp": Descriptors.MolLogP(mol),
                    }
                    activity = outcome.get("pincer_score") or 0.0
                    learning_loop.add_observation(
                        cid, features, {"activity": activity}
                    )
                    added += 1
            except Exception:
                pass
        logger.info(f"CompoundMemory: warmed LearningLoop with {added} observations")
        return added

    def __len__(self):
        return len(self._outcomes)

    def save(self, path: str):
        data = {
            "outcomes": self._outcomes,
            "index_dim": self._retriever.dim,
        }
        Path(path).write_text(json.dumps(data, default=str))
        # Also save the vector index
        self._retriever.save(path + ".vectors")

    def load(self, path: str):
        if Path(path).exists():
            data = json.loads(Path(path).read_text())
            self._outcomes = data.get("outcomes", {})
            if Path(path + ".vectors").exists():
                self._retriever.load(path + ".vectors")


# ---------------------------------------------------------------------------
# KhukuriNumen — unified interface wiring all 3 into Khukuri
# ---------------------------------------------------------------------------

class KhukuriNumen:
    """
    Unified Numen integration for Khukuri.

    Owns all three indexes and wires them into the right places.

    Usage in AMRDiscoveryWorkflow:
        self.numen = KhukuriNumen(project="mrsa")
        self.numen.index_threats(viable_threats)
        self.numen.mine_literature("S. aureus MRSA")
        worst = self.numen.threats.get_worst_threats(smiles, k=20)
        similar = self.numen.compounds.recall(smiles, k=5)
    """

    def __init__(
        self,
        project: str = "default",
        memory_dir: Optional[str] = None,
        threat_dim: int = 256,
        literature_dim: int = 1024,
        compound_dim: int = 512,
    ):
        self.project = project
        self.base_dir = Path(memory_dir or Path.home() / ".khukuri" / "numen") / project
        self.base_dir.mkdir(parents=True, exist_ok=True)

        self.threats   = ThreatIndex(dim=threat_dim)
        self.literature = LiteratureIndex(
            dim=literature_dim,
            cache_dir=str(self.base_dir / "literature"),
        )
        self.compounds = CompoundMemory(dim=compound_dim)

        self._restore()
        logger.info(
            f"KhukuriNumen: threats={len(self.threats)} "
            f"papers={len(self.literature._retriever)} "
            f"compounds={len(self.compounds)}"
        )

    def index_threats(self, threats: List) -> int:
        n = self.threats.build(threats)
        self._save_threats()
        return n

    def mine_literature(
        self,
        pathogen: str,
        resistance_genes: Optional[List[str]] = None,
    ) -> Dict:
        result = self.literature.mine_for_targets(pathogen, resistance_genes)
        self.literature.save()
        return result

    def remember_compound(self, smiles: str, **kwargs):
        self.compounds.remember(smiles, **kwargs)
        self._save_compounds()

    def get_worst_threats(self, smiles: str, k: int = 20) -> List:
        return self.threats.get_worst_threats(smiles, k=k)

    def recall_similar(self, smiles: str, k: int = 5) -> List[Dict]:
        return self.compounds.recall(smiles, k=k)

    def search_literature(self, query: str, top_k: int = 10) -> List[Dict]:
        return self.literature.search(query, top_k=top_k)

    def _save_threats(self):
        pass  # ThreatIndex is rebuilt each run from live threats

    def _save_compounds(self):
        self.compounds.save(str(self.base_dir / "compounds.json"))

    def _restore(self):
        # Literature index persists across sessions
        self.literature.load()
        # Compound memory persists across sessions
        path = str(self.base_dir / "compounds.json")
        self.compounds.load(path)
