"""
HGTMapper — Horizontal Gene Transfer threat modelling for PINCER.

Models HGT as a network process on the KnowledgeGraph:
    - Nodes: resistance genes / alleles, bacterial species
    - Edges: transfer events weighted by probability

Produces a TransferThreatCluster — the set of resistance gene variants
that could plausibly *arrive* in the target strain via MGE acquisition,
analogous to MutationSpaceMapper's ViableMutant cluster for vertical evolution.

Both threat sources feed PincerEngine as parallel red teams:
    vertical threats  → MutationSpaceMapper  (mutation)
    horizontal threats → HGTMapper           (gene transfer)

Data sources supported:
    - CARD (Comprehensive Antibiotic Resistance Database)
    - PATRIC / BV-BRC
    - User-supplied JSON (offline / curated)

The graph lives in KnowledgeGraph so it persists across sessions via
PersistentMemory and is accessible to KosmosEngine for reasoning.
"""

import json
import logging
import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

import numpy as np

logger = logging.getLogger('khukuri')


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class ResistanceGeneAllele:
    """
    A resistance gene variant that could be acquired via HGT.

    Analogous to ViableMutant — same interface so PincerEngine
    can treat both threat types uniformly.
    """
    gene_id: str                     # e.g. "mecA", "vanA", "blaZ"
    allele: str                      # e.g. "mecA_1", "mecA_MRSA252"
    donor_species: str               # e.g. "S. epidermidis"
    mge_type: str                    # plasmid | transposon | integron | phage | ICE
    resistance_class: str            # beta-lactam | glycopeptide | etc.
    target_gene: str                 # gene it disrupts / encodes resistance to
    transfer_probability: float      # P(transfer | contact) in [0, 1]
    fitness_cost: float              # fitness burden on recipient in [0, 1]
    fitness: float = 0.0             # 1 - fitness_cost (mirrors ViableMutant)
    protein_sequence: Optional[str] = None
    vector: Optional[np.ndarray] = None    # embedding for ThreatAwareFitnessFunction

    def __post_init__(self):
        self.fitness = max(0.0, 1.0 - self.fitness_cost)

    # Expose same attribute as ViableMutant for duck-typing in PincerEngine
    @property
    def mutations(self):
        """Duck-type compatibility with ViableMutant.mutations."""
        return [_FakeMutation(self.gene_id, self.fitness_cost)]


@dataclass
class _FakeMutation:
    """Shim so ResistanceGeneAllele.mutations is iterable like ViableMutant."""
    wild_type: str
    fitness_cost: float
    position: int = 0
    mutant: str = "HGT"


@dataclass
class TransferThreatCluster:
    """All HGT threats relevant to a target strain."""
    target_strain: str
    alleles: List[ResistanceGeneAllele] = field(default_factory=list)
    network_stats: Dict[str, Any] = field(default_factory=dict)

    def to_viable_threats(self) -> List[ResistanceGeneAllele]:
        """Return alleles sorted by transfer probability × fitness."""
        return sorted(
            self.alleles,
            key=lambda a: a.transfer_probability * a.fitness,
            reverse=True,
        )


# ---------------------------------------------------------------------------
# MGE transfer probability model
# ---------------------------------------------------------------------------

# Base transfer rates by MGE type (per cell per generation, approximate)
# Based on Baquero et al. 2013, Lipsitch et al. 2002
MGE_BASE_RATES = {
    "conjugative_plasmid": 1e-3,
    "plasmid":             1e-4,
    "transposon":          1e-5,
    "integron":            1e-5,
    "ICE":                 1e-4,   # Integrative Conjugative Elements
    "phage":               1e-6,
    "unknown":             1e-5,
}

# Phylogenetic distance penalty — transfer drops off with distance
# Donor/recipient in same genus: multiplier ~1.0
# Different genus, same family: ~0.1
# Different family: ~0.01
PHYLO_DISTANCE_PENALTY = {
    "same_species":  1.0,
    "same_genus":    0.8,
    "same_family":   0.1,
    "same_order":    0.01,
    "different":     0.001,
}

# Antibiotic selective pressure amplifies transfer
SELECTIVE_PRESSURE_MULTIPLIER = {
    "high":   10.0,
    "medium":  3.0,
    "low":     1.0,
    "none":    0.5,
}


def estimate_transfer_probability(
    mge_type: str,
    phylo_distance: str = "same_genus",
    selective_pressure: str = "medium",
    known_rate: Optional[float] = None,
) -> float:
    """
    Estimate P(transfer) from MGE type, phylogenetic distance,
    and antibiotic selective pressure.
    """
    if known_rate is not None:
        return float(np.clip(known_rate, 0.0, 1.0))

    base = MGE_BASE_RATES.get(mge_type, MGE_BASE_RATES["unknown"])
    phylo = PHYLO_DISTANCE_PENALTY.get(phylo_distance, 0.01)
    pressure = SELECTIVE_PRESSURE_MULTIPLIER.get(selective_pressure, 1.0)

    prob = base * phylo * pressure
    return float(np.clip(prob, 0.0, 1.0))


# ---------------------------------------------------------------------------
# HGT Knowledge Graph builder
# ---------------------------------------------------------------------------

class HGTKnowledgeGraphBuilder:
    """
    Populates a KnowledgeGraph with HGT network data.

    Relationship types added:
        species   --[transfers_to]-->  species      (MGE transfer event)
        species   --[carries]-->       gene_allele  (gene presence)
        gene_allele --[confers_resistance_to]--> drug_class
        gene_allele --[disrupts]-->    target_gene
    """

    # Curated MRSA-relevant HGT data (built-in, no API needed)
    # Sources: CARD, NCBI AMRFinderPlus, literature
    MRSA_BUILTIN_DATA = {
        "genes": [
            # Beta-lactam resistance
            {
                "gene_id": "mecA",
                "alleles": ["mecA_1", "mecA_MRSA252", "mecA_MW2"],
                "resistance_class": "beta-lactam",
                "target_gene": "pbp2a",
                "mge_type": "transposon",
                "mge_name": "SCCmec",
                "donor_species": ["S. epidermidis", "S. haemolyticus", "S. sciuri"],
                "fitness_cost": 0.05,
            },
            {
                "gene_id": "mecC",
                "alleles": ["mecC_1"],
                "resistance_class": "beta-lactam",
                "target_gene": "pbp2a",
                "mge_type": "transposon",
                "mge_name": "SCCmec_XI",
                "donor_species": ["S. xylosus", "M. sciuri"],
                "fitness_cost": 0.08,
            },
            {
                "gene_id": "blaZ",
                "alleles": ["blaZ_1", "blaZ_2", "blaZ_staph"],
                "resistance_class": "beta-lactam",
                "target_gene": "penicillin_binding_protein",
                "mge_type": "plasmid",
                "mge_name": "pI524",
                "donor_species": ["S. aureus", "S. epidermidis"],
                "fitness_cost": 0.03,
            },
            # Glycopeptide resistance
            {
                "gene_id": "vanA",
                "alleles": ["vanA_1", "vanA_Tn1546"],
                "resistance_class": "glycopeptide",
                "target_gene": "ddl",
                "mge_type": "transposon",
                "mge_name": "Tn1546",
                "donor_species": ["E. faecalis", "E. faecium"],
                "fitness_cost": 0.12,
            },
            {
                "gene_id": "vanB",
                "alleles": ["vanB_1", "vanB_2"],
                "resistance_class": "glycopeptide",
                "target_gene": "ddl",
                "mge_type": "transposon",
                "mge_name": "Tn5382",
                "donor_species": ["E. faecium", "E. faecalis"],
                "fitness_cost": 0.10,
            },
            # Aminoglycoside resistance
            {
                "gene_id": "aac6-aph2",
                "alleles": ["aac6-aph2_1"],
                "resistance_class": "aminoglycoside",
                "target_gene": "16S_rRNA",
                "mge_type": "transposon",
                "mge_name": "Tn4001",
                "donor_species": ["S. epidermidis", "E. faecalis"],
                "fitness_cost": 0.07,
            },
            # Fluoroquinolone resistance (less common via HGT in MRSA but present)
            {
                "gene_id": "qnrS",
                "alleles": ["qnrS_1"],
                "resistance_class": "fluoroquinolone",
                "target_gene": "gyrA",
                "mge_type": "plasmid",
                "mge_name": "pINF5",
                "donor_species": ["E. coli", "K. pneumoniae"],
                "fitness_cost": 0.15,
            },
            # Efflux pump (multi-drug)
            {
                "gene_id": "norA",
                "alleles": ["norA_1", "norA_overexp"],
                "resistance_class": "multidrug",
                "target_gene": "efflux_pump",
                "mge_type": "integron",
                "mge_name": "In_norA",
                "donor_species": ["S. aureus", "S. haemolyticus"],
                "fitness_cost": 0.04,
            },
            # Tetracycline
            {
                "gene_id": "tetM",
                "alleles": ["tetM_1", "tetM_Tn916"],
                "resistance_class": "tetracycline",
                "target_gene": "16S_rRNA",
                "mge_type": "ICE",
                "mge_name": "Tn916",
                "donor_species": ["E. faecalis", "S. agalactiae"],
                "fitness_cost": 0.06,
            },
        ],
        "transfer_events": [
            # Known MRSA acquisition events from literature
            {"donor": "S. epidermidis", "recipient": "S. aureus",
             "gene": "mecA", "evidence": "genomic", "selective_pressure": "high"},
            {"donor": "E. faecalis", "recipient": "S. aureus",
             "gene": "vanA", "evidence": "clinical", "selective_pressure": "high"},
            {"donor": "E. faecium", "recipient": "S. aureus",
             "gene": "vanB", "evidence": "experimental", "selective_pressure": "medium"},
            {"donor": "S. haemolyticus", "recipient": "S. aureus",
             "gene": "aac6-aph2", "evidence": "genomic", "selective_pressure": "medium"},
        ],
    }

    def __init__(self, knowledge_graph):
        self.kg = knowledge_graph

    def load_builtin_mrsa(self, selective_pressure: str = "high") -> int:
        """
        Load built-in MRSA HGT data into the KnowledgeGraph.
        Returns number of alleles loaded.
        """
        data = self.MRSA_BUILTIN_DATA
        count = 0

        for gene_entry in data["genes"]:
            gene_id = gene_entry["gene_id"]
            mge_type = gene_entry["mge_type"]

            for allele_name in gene_entry["alleles"]:
                node_id = f"hgt_allele_{allele_name}"
                self.kg.add_gene(node_id, {
                    "gene_id": gene_id,
                    "allele": allele_name,
                    "resistance_class": gene_entry["resistance_class"],
                    "target_gene": gene_entry["target_gene"],
                    "mge_type": mge_type,
                    "mge_name": gene_entry.get("mge_name", "unknown"),
                    "fitness_cost": gene_entry["fitness_cost"],
                    "node_subtype": "hgt_allele",
                })

                # Resistance class edge
                self.kg.add_relationship(
                    node_id,
                    gene_entry["resistance_class"],
                    "confers_resistance_to",
                    {"mechanism": mge_type},
                )
                # Target gene edge
                self.kg.add_relationship(
                    node_id,
                    gene_entry["target_gene"],
                    "disrupts",
                )

                # Donor species → allele transfer edges
                for donor in gene_entry["donor_species"]:
                    donor_node = f"species_{donor.replace(' ', '_')}"
                    self.kg.add_strain(donor_node, {
                        "species": donor,
                        "node_subtype": "donor_species",
                    })
                    phylo = _estimate_phylo_distance(donor, "S. aureus")
                    prob = estimate_transfer_probability(
                        mge_type, phylo, selective_pressure
                    )
                    self.kg.add_relationship(
                        donor_node,
                        node_id,
                        "carries",
                        {
                            "transfer_probability": prob,
                            "mge_type": mge_type,
                            "phylo_distance": phylo,
                            "selective_pressure": selective_pressure,
                        },
                    )
                count += 1

        # Add known transfer events as direct edges
        for event in data["transfer_events"]:
            donor_node = f"species_{event['donor'].replace(' ', '_')}"
            recipient_node = f"species_{event['recipient'].replace(' ', '_')}"
            self.kg.add_relationship(
                donor_node,
                recipient_node,
                "transfers_to",
                {
                    "gene": event["gene"],
                    "evidence": event["evidence"],
                    "selective_pressure": event["selective_pressure"],
                },
            )

        logger.info(f"HGTKnowledgeGraphBuilder: loaded {count} MRSA alleles")
        return count

    def load_from_json(self, path: str, selective_pressure: str = "medium") -> int:
        """
        Load HGT data from a user-supplied JSON file.

        Expected format (same as MRSA_BUILTIN_DATA):
        {
            "genes": [...],
            "transfer_events": [...]
        }
        """
        with open(path) as f:
            data = json.load(f)

        # Temporarily swap builtin data and reuse loader
        original = self.MRSA_BUILTIN_DATA
        self.MRSA_BUILTIN_DATA = data
        count = self.load_builtin_mrsa(selective_pressure)
        self.MRSA_BUILTIN_DATA = original
        return count


# ---------------------------------------------------------------------------
# HGT Mapper — the red team for horizontal threats
# ---------------------------------------------------------------------------

class HGTMapper:
    """
    Map the horizontal gene transfer threat landscape for a target strain.

    Queries the KnowledgeGraph to extract which resistance gene alleles
    could plausibly be acquired by the target strain, weighted by:
        - Transfer probability (MGE type × phylogenetic distance × pressure)
        - Fitness of the allele in the recipient
        - Evidence quality (genomic > clinical > experimental > predicted)

    Produces a TransferThreatCluster compatible with PincerEngine.
    """

    EVIDENCE_WEIGHTS = {
        "genomic":      1.0,
        "clinical":     0.9,
        "experimental": 0.7,
        "predicted":    0.4,
    }

    def __init__(
        self,
        knowledge_graph,
        min_transfer_probability: float = 1e-6,
        fitness_threshold: float = 0.2,
    ):
        self.kg = knowledge_graph
        self.min_prob = min_transfer_probability
        self.fitness_threshold = fitness_threshold
        self._builder = HGTKnowledgeGraphBuilder(knowledge_graph)

    def build_mrsa_network(self, selective_pressure: str = "high") -> int:
        """Convenience: load built-in MRSA data if graph is empty."""
        hgt_nodes = [
            n for n, d in self.kg.graph.nodes(data=True)
            if d.get("node_subtype") == "hgt_allele"
        ]
        if not hgt_nodes:
            return self._builder.load_builtin_mrsa(selective_pressure)
        logger.info(f"HGTMapper: using {len(hgt_nodes)} existing alleles in graph")
        return len(hgt_nodes)

    def map_threats(
        self,
        target_strain: str = "S. aureus",
        resistance_classes: Optional[List[str]] = None,
        selective_pressure: str = "high",
        auto_load: bool = True,
    ) -> TransferThreatCluster:
        """
        Map all viable HGT threats for the target strain.

        Args:
            target_strain: recipient species name
            resistance_classes: filter to specific classes (None = all)
            selective_pressure: "high" | "medium" | "low"
            auto_load: load built-in MRSA data if graph empty

        Returns:
            TransferThreatCluster
        """
        if auto_load:
            self.build_mrsa_network(selective_pressure)

        alleles = []

        # Walk all HGT allele nodes in the graph
        for node_id, node_data in self.kg.graph.nodes(data=True):
            if node_data.get("node_subtype") != "hgt_allele":
                continue

            res_class = node_data.get("resistance_class", "")
            if resistance_classes and res_class not in resistance_classes:
                continue

            # Find all donor species that carry this allele
            max_prob = 0.0
            for pred in self.kg.graph.predecessors(node_id):
                edge_data = self.kg.graph.get_edge_data(pred, node_id)
                if not edge_data:
                    continue
                for _, edata in edge_data.items():
                    if edata.get("relationship_type") == "carries":
                        prob = edata.get("transfer_probability", 0.0)
                        max_prob = max(max_prob, prob)

            if max_prob < self.min_prob:
                continue

            fitness_cost = node_data.get("fitness_cost", 0.5)
            fitness = 1.0 - fitness_cost

            if fitness < self.fitness_threshold:
                continue

            # Build the allele object
            allele = ResistanceGeneAllele(
                gene_id=node_data.get("gene_id", node_id),
                allele=node_data.get("allele", node_id),
                donor_species=self._get_donors(node_id),
                mge_type=node_data.get("mge_type", "unknown"),
                resistance_class=res_class,
                target_gene=node_data.get("target_gene", "unknown"),
                transfer_probability=max_prob,
                fitness_cost=fitness_cost,
            )
            allele.vector = self._embed_allele(allele)
            alleles.append(allele)

        cluster = TransferThreatCluster(
            target_strain=target_strain,
            alleles=alleles,
            network_stats=self._network_stats(),
        )

        logger.info(
            f"HGTMapper: {len(alleles)} transfer threats for {target_strain} "
            f"(pressure={selective_pressure})"
        )
        return cluster

    def get_high_risk_genes(self, top_n: int = 5) -> List[Dict]:
        """
        Return the highest-risk genes ranked by transfer_probability × fitness.
        Useful for agent meeting context.
        """
        cluster = self.map_threats()
        ranked = sorted(
            cluster.alleles,
            key=lambda a: a.transfer_probability * a.fitness,
            reverse=True,
        )
        return [
            {
                "gene_id": a.gene_id,
                "allele": a.allele,
                "mge_type": a.mge_type,
                "resistance_class": a.resistance_class,
                "transfer_probability": a.transfer_probability,
                "fitness": a.fitness,
                "donor_species": a.donor_species,
            }
            for a in ranked[:top_n]
        ]

    def _get_donors(self, allele_node_id: str) -> str:
        """Return comma-separated donor species for an allele node."""
        donors = []
        for pred in self.kg.graph.predecessors(allele_node_id):
            pred_data = self.kg.graph.nodes[pred]
            if pred_data.get("node_subtype") == "donor_species":
                donors.append(pred_data.get("species", pred))
        return ", ".join(donors) if donors else "unknown"

    def _embed_allele(self, allele: ResistanceGeneAllele) -> np.ndarray:
        """
        Embed a resistance gene allele into the same 64-dim space
        as ViableMutant.vector and drug embeddings.
        """
        # Resistance class one-hot (8 classes)
        classes = [
            "beta-lactam", "glycopeptide", "aminoglycoside",
            "fluoroquinolone", "tetracycline", "multidrug",
            "macrolide", "other",
        ]
        one_hot = [1.0 if allele.resistance_class == c else 0.0 for c in classes]

        # MGE mobility score (conjugative = highest mobility)
        mge_scores = {
            "conjugative_plasmid": 1.0,
            "plasmid": 0.8,
            "ICE": 0.7,
            "transposon": 0.5,
            "integron": 0.4,
            "phage": 0.3,
            "unknown": 0.2,
        }
        mobility = mge_scores.get(allele.mge_type, 0.2)

        feats = (
            one_hot
            + [
                allele.transfer_probability * 100,   # scale up
                allele.fitness_cost,
                allele.fitness,
                mobility,
            ]
        )
        dim = 64
        tiles = math.ceil(dim / len(feats))
        full = np.tile(np.array(feats, dtype=np.float64), tiles)[:dim]

        norm = np.linalg.norm(full)
        if norm > 0:
            full /= norm
        return full

    def _network_stats(self) -> Dict:
        try:
            import networkx as nx
            g = self.kg.graph
            hgt_nodes = [
                n for n, d in g.nodes(data=True)
                if d.get("node_subtype") == "hgt_allele"
            ]
            species_nodes = [
                n for n, d in g.nodes(data=True)
                if d.get("node_subtype") == "donor_species"
            ]
            return {
                "total_alleles": len(hgt_nodes),
                "donor_species": len(species_nodes),
                "transfer_edges": sum(
                    1 for _, _, d in g.edges(data=True)
                    if d.get("relationship_type") == "transfers_to"
                ),
            }
        except Exception:
            return {}


# ---------------------------------------------------------------------------
# PincerEngine extension — dual red team
# ---------------------------------------------------------------------------

class DualRedTeamMixin:
    """
    Mixin for PincerEngine that adds HGT as a second red team.

    Usage:
        from src.resistance.hgt_mapper import DualRedTeamMixin
        from src.resistance.pincer_engine import PincerEngine

        class PincerEngineWithHGT(DualRedTeamMixin, PincerEngine):
            pass

        engine = PincerEngineWithHGT(population_size=100, n_generations=50)
        engine.map_threats(wild_type_seq, active_site)          # vertical
        engine.map_hgt_threats(knowledge_graph)                  # horizontal
        apex = engine.evolve(seed_smiles, fitness_fn=fitness_fn)
    """

    def map_hgt_threats(
        self,
        knowledge_graph,
        target_strain: str = "S. aureus",
        selective_pressure: str = "high",
        resistance_classes: Optional[List[str]] = None,
    ) -> TransferThreatCluster:
        """Map HGT threats and merge into viable_threats."""
        mapper = HGTMapper(knowledge_graph)
        cluster = mapper.map_threats(
            target_strain=target_strain,
            selective_pressure=selective_pressure,
            resistance_classes=resistance_classes,
        )

        hgt_as_threats = cluster.to_viable_threats()

        # Merge with existing vertical threats
        if not hasattr(self, "viable_threats"):
            self.viable_threats = []

        self.viable_threats.extend(hgt_as_threats)
        self._hgt_cluster = cluster

        logger.info(
            f"DualRedTeam: {len(self.viable_threats)} total threats "
            f"({len(hgt_as_threats)} HGT + "
            f"{len(self.viable_threats) - len(hgt_as_threats)} mutational)"
        )
        return cluster

    def get_threat_breakdown(self) -> Dict:
        """Break down threats by type."""
        hgt = getattr(self, "_hgt_cluster", None)
        hgt_count = len(hgt.alleles) if hgt else 0
        total = len(getattr(self, "viable_threats", []))
        return {
            "total": total,
            "mutational": total - hgt_count,
            "hgt": hgt_count,
            "hgt_genes": (
                [a.gene_id for a in hgt.alleles] if hgt else []
            ),
        }


# ---------------------------------------------------------------------------
# Convenience: pre-wired PincerEngine with HGT
# ---------------------------------------------------------------------------

def make_pincer_with_hgt(population_size=100, n_generations=50, **kwargs):
    """
    Factory that returns a PincerEngine subclass with HGT support.

    from src.resistance.hgt_mapper import make_pincer_with_hgt
    engine = make_pincer_with_hgt()
    """
    from .pincer_engine import PincerEngine

    class PincerEngineHGT(DualRedTeamMixin, PincerEngine):
        pass

    return PincerEngineHGT(
        population_size=population_size,
        n_generations=n_generations,
        **kwargs,
    )


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _estimate_phylo_distance(donor: str, recipient: str) -> str:
    """
    Crude phylogenetic distance estimate from species names.
    Real implementation would use a taxonomy API or pre-built tree.
    """
    d_genus = donor.split()[0] if donor else ""
    r_genus = recipient.split()[0] if recipient else ""

    if donor == recipient:
        return "same_species"
    if d_genus == r_genus:
        return "same_genus"
    # Staphylococcus → Enterococcus is same family (Bacillota / Firmicutes)
    staph_family = {"Staphylococcus", "Macrococcus", "Salinicoccus"}
    entero_family = {"Enterococcus", "Streptococcus", "Lactococcus"}
    if d_genus in staph_family and r_genus in staph_family:
        return "same_family"
    if d_genus in (staph_family | entero_family) and r_genus in (staph_family | entero_family):
        return "same_order"
    return "different"
