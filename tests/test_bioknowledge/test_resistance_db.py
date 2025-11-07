"""Tests for ResistanceDatabase"""

import pytest
from src.bioknowledge import ResistanceDatabase


def test_resistance_db_initialization():
    """Test database initialization"""
    db = ResistanceDatabase()
    assert len(db.genes) > 0
    assert len(db.mechanisms) > 0


def test_query_gene():
    """Test gene query"""
    db = ResistanceDatabase()
    gene_info = db.query_gene('mecA')
    assert gene_info is not None
    assert gene_info['type'] == 'beta_lactam'
    assert gene_info['organism'] == 'S. aureus'


def test_get_genes_by_organism():
    """Test organism-specific gene query"""
    db = ResistanceDatabase()
    genes = db.get_genes_by_organism('E. coli')
    assert len(genes) > 0


def test_get_genes_by_drug_class():
    """Test drug class query"""
    db = ResistanceDatabase()
    genes = db.get_genes_by_drug_class('beta_lactam')
    assert len(genes) > 0


def test_add_gene():
    """Test adding new gene"""
    db = ResistanceDatabase()
    db.add_gene('testGene', {'type': 'test', 'organism': 'Test', 'mechanism': 'test'})
    assert 'testGene' in db.genes
