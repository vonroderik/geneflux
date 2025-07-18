import pytest
from utils import dna_tools as gene
from main import validate_sequence


# Should raise a ValueError in case of an invalid DNA/RNA sequence
def test_validate_sequence():

    assert validate_sequence("CGTA") == ("CGTA", "dna")
    assert validate_sequence("CGUA") == ("CGUA", "rna")

    with pytest.raises(ValueError):
        validate_sequence("FOOBAR")
    with pytest.raises(ValueError):
        validate_sequence("CGTAU")


# Should return the complimentary DNA of the user´s inputed sequence
def test_complementary_dna():

    assert gene.complementary_dna("ATCG", False) == "TAGC"
    assert gene.complementary_dna("AUCG", False) == "TAGC"


# Should return the the RNA transcript of the user´s inputed DNA
def test_transcription():

    assert gene.transcription("CGTA", False) == "CGUA"
    assert gene.transcription("AUGA", False) == "AUGA"


# Should return the aminoacid sequence of the user inputed DNA or RNA sequence
# if it has at least one Start Codong and one Stop Codon
def test_translation(monkeypatch):

    monkeypatch.setattr("builtins.input", lambda _: "1")
    assert gene.translation("AUGCGUUGA", False) == ["Met", "Arg", "Ter"]
    assert gene.translation("AA", False) == None


# Should return the original RNA sequence if it doesn´t have a Start and/or Stop codons
def test_orf_finder():

    assert gene.orf_finder("AUGCGACGAUGACGA", False, False) == ["AUGCGACGA"]
    assert gene.orf_finder("GGGGGGGGG", False, False) == []


def test_reverse_transcription():

    assert gene.reverse_transcription("UUUU") == "TTTT"
