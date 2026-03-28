import pytest
from genedesign.checkers.internal_rbs_checker import InternalRBSChecker

@pytest.fixture
def checker():
    checker = InternalRBSChecker()
    checker.initiate()
    return checker

def test_detects_internal_rbs(checker):
    # Sequences with SD-like motifs 4-12 bp upstream of ATG
    rbs_seqs = [
        "TTTAAGGAGATCGATGCCC",     # AAGGAG 6bp before ATG
        "TTTTAGGAGTTTTTATGCCC",    # AGGAG 8bp before ATG
        "CCCAGGAGGTTTTATGAAA",     # AGGAGG before ATG
        "AAAGAAGGAGCCCCCATGTTT",   # AAGGAG with 6bp spacing
        "TTTGGAGGCCCATGCCC",       # GGAGG 5bp before ATG
    ]

    for seq in rbs_seqs:
        result, detail = checker.run(seq)
        assert result == False, f"Expected False for {seq}, got {result}"

def test_passes_no_internal_rbs(checker):
    # Sequences without SD motifs near ATG
    clean_seqs = [
        "TTTCCCGGGAAATTTCCC",          # no ATG at all
        "CCCCCCCCCATGCCC",              # ATG but no SD upstream
        "AAGGAGATGCCC",                 # SD too close to ATG (< 4bp spacing)
        "AAGGAGTTTTTTTTTTTTTTTATGCCC",  # SD too far from ATG (> 12bp)
        "TTTCCCTTTGGGATTTCGATCG",       # no ATG, no SD
    ]

    for seq in clean_seqs:
        result, detail = checker.run(seq)
        assert result == True, f"Expected True for {seq}, got {result} ({detail})"
