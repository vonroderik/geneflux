from Bio.Seq import Seq
from Bio.SeqUtils import seq3, molecular_weight, gc_fraction
import re
from typing import List, Tuple, Literal

def validate_sequence(seq: str) -> Tuple[str, Literal["dna", "rna"]]:
    """Valida se a entrada é uma sequência de DNA ou RNA."""
    if not seq:
        raise ValueError("Entrada vazia. Por favor, digite uma sequência de DNA ou RNA")

    match = re.fullmatch(r"^([ATCGU]+)$", seq.upper())
    if not match:
        raise ValueError(
            "Sequência inválida. Deve conter apenas os caracteres A, T, C, G ou U."
        )

    validated_seq = match.group(1)
    if len(validated_seq) < 3:
        raise ValueError("Sua sequência deve ter pelo menos 3 nucleotídeos")

    if "U" in validated_seq and "T" in validated_seq:
        raise ValueError(
            "Sequência de DNA/RNA inválida: contém tanto timina (T) quanto uracila (U)"
        )

    return validated_seq, "dna" if "U" not in validated_seq else "rna"

def reverse_transcription(seq: str) -> str:
    """Converte uma sequência de RNA para DNA."""
    return seq if "U" not in seq else seq.replace("U", "T")

def complementary_dna(seq: str) -> str:
    """Retorna o filamento complementar de DNA."""
    dna_seq = reverse_transcription(seq)
    return dna_seq.translate(str.maketrans("ATCG", "TAGC"))

def transcription(seq: str) -> str:
    """Retorna um transcrito de RNA."""
    return seq.replace("T", "U")

def translation(seq: str, to_stop: bool = False) -> List[str]:
    """Traduz a sequência de RNA para aminoácidos."""
    # Ensure it's a multiple of 3
    trimmed_seq = seq[: len(seq) - len(seq) % 3]
    if not trimmed_seq:
        return []
    
    protein = str(Seq(trimmed_seq).translate(to_stop=to_stop))
    return [seq3(aa).title() for aa in protein]

def orf_finder(seq: str) -> List[str]:
    """Encontra ORFs em uma sequência de RNA."""
    rna_seq = transcription(seq)
    if len(rna_seq) < 3:
        return []

    orfs = []
    stop_codons = {"UAA", "UAG", "UGA"}
    for frame in range(3):
        i = frame
        while i + 3 <= len(rna_seq):
            codon = rna_seq[i : i + 3]
            if codon == "AUG":
                start = i
                j = i
                found_stop = False
                while j + 3 <= len(rna_seq):
                    next_codon = rna_seq[j : j + 3]
                    if next_codon in stop_codons:
                        orf = rna_seq[start:j]
                        orfs.append(orf)
                        found_stop = True
                        break
                    j += 3
                if found_stop:
                    i = j + 3
                else:
                    i += 3
            else:
                i += 3
    return orfs

def calculate_gc_content(seq: str) -> float:
    """Calcula o conteúdo GC da sequência."""
    dna_seq = reverse_transcription(seq)
    return gc_fraction(dna_seq) * 100

def calculate_molecular_weight(seq: str, seq_type: Literal["dna", "rna", "protein"]) -> float:
    """Calcula o peso molecular da sequência."""
    if seq_type == "dna":
        return molecular_weight(reverse_transcription(seq), seq_type="DNA")
    elif seq_type == "rna":
        return molecular_weight(transcription(seq), seq_type="RNA")
    elif seq_type == "protein":
        # Seq expects single letter amino acids for protein MW
        # If input is a list of 3-letter codes, we need to convert back
        # But here we assume it's the raw RNA/DNA sequence to be translated or already a protein string
        return molecular_weight(seq, seq_type="protein")
    return 0.0
