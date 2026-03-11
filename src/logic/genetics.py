from collections import Counter
from itertools import product
from typing import List, Dict, Tuple

def validate_genotype(genotype: str, mode: str = "monohybrid"):
    """Verifica se o genótipo entrado pelo usuário é válido."""
    if not genotype:
        raise ValueError("Entrada vazia. Por favor, digite um genótipo válido.")
    
    genotype = genotype.strip()
    
    if mode == "monohybrid":
        if len(genotype) != 2:
            raise ValueError("Um genótipo monohíbrido deve ter 2 alelos.")
    elif mode == "dihybrid":
        if len(genotype) != 4:
            raise ValueError("Um genótipo dihíbrido deve ter 4 alelos (ex: AaBb).")
    
    # Simple validation for now, can be expanded
    if not genotype.isalpha():
        raise ValueError("O genótipo deve conter apenas letras.")

def calculate_punnett(genotype1: str, genotype2: str, mode: str = "monohybrid") -> Dict:
    """Calcula o quadro de Punnett e frequências."""
    if mode == "monohybrid":
        gametes1 = list(genotype1)
        gametes2 = list(genotype2)
    else: # dihybrid
        # Gametes are combinations of the two genes
        # Ex: AaBb -> AB, Ab, aB, ab
        gametes1 = ["".join(p) for p in product(genotype1[0:2], genotype1[2:4])]
        gametes2 = ["".join(p) for p in product(genotype2[0:2], genotype2[2:4])]
    
    offspring = []
    for g1 in gametes1:
        for g2 in gametes2:
            # Sort alleles within each gene to normalize (A with a, B with b)
            if mode == "monohybrid":
                offspring.append("".join(sorted(g1 + g2)))
            else:
                # Dihybrid: AABB + aabb -> AaBb
                # Gene 1
                gene1 = "".join(sorted(g1[0] + g2[0]))
                # Gene 2
                gene2 = "".join(sorted(g1[1] + g2[1]))
                offspring.append(gene1 + gene2)
    
    return {
        "gametes1": gametes1,
        "gametes2": gametes2,
        "offspring": offspring,
        "genotype_freq": Counter(offspring)
    }

def get_abo_phenotype(genotype: str) -> str:
    """Determina o fenótipo ABO a partir do genótipo."""
    gen = genotype.upper()
    if "A" in gen and "B" in gen:
        return "AB"
    if "A" in gen:
        return "A"
    if "B" in gen:
        return "B"
    return "O"

def calculate_abo_frequencies(offspring: List[str]) -> Dict[str, int]:
    """Calcula frequências fenotípicas do sistema ABO."""
    phenotypes = [get_abo_phenotype(gen) for gen in offspring]
    return Counter(phenotypes)
