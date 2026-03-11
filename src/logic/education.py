import random
from typing import List, Dict, Tuple
from .dna_rna import transcription, translation

def generate_random_dna(length: int = 18) -> str:
    """Gera uma sequência de DNA aleatória começando com ATG e terminando com um códon de parada."""
    STOP_CODONS = ["TGA", "TAA", "TAG"]
    middle_length = length - 6 # 3 for start, 3 for stop
    if middle_length < 0: middle_length = 0
    
    random_middle = "".join(random.choices(["C", "G", "T", "A"], k=middle_length))
    return "ATG" + random_middle + random.choice(STOP_CODONS)

def get_mutation_impact(original_dna: str, mutated_dna: str) -> Dict:
    """Determina o tipo de mutação e seu impacto na proteína."""
    original_rna = transcription(original_dna)
    mutated_rna = transcription(mutated_dna)
    
    original_protein = translation(original_rna)
    mutated_protein = translation(mutated_rna)
    
    # Simple mutation type identification
    if len(original_dna) == len(mutated_dna):
        dna_mutation = "Substituição"
    elif len(mutated_dna) > len(original_dna):
        dna_mutation = "Inserção"
    else:
        dna_mutation = "Deleção"
        
    # Protein impact
    if dna_mutation == "Substituição":
        if original_protein == mutated_protein:
            protein_mutation = "Silenciosa"
        elif "Stop" in mutated_protein and len(mutated_protein) < len(original_protein):
             protein_mutation = "Nonsense (Sem Sentido)"
        else:
            protein_mutation = "Missense (Sentido Trocado)"
    else: # Inserção/Deleção
        if (len(mutated_dna) - len(original_dna)) % 3 != 0:
            protein_mutation = "Frameshift (Alteração de Quadro de Leitura)"
        else:
            protein_mutation = "In-frame (Manutenção de Quadro de Leitura)"
            
    return {
        "dna_type": dna_mutation,
        "protein_type": protein_mutation,
        "original_protein": original_protein,
        "mutated_protein": mutated_protein
    }

def generate_pcr_challenge() -> Dict:
    """Gera uma sequência e os primers corretos para o desafio de PCR."""
    template = generate_random_dna(length=40)
    # Forward primer: first 10 bases
    fwd = template[:10]
    # Reverse primer: reverse complement of the last 10 bases
    from .dna_rna import complementary_dna
    last_10 = template[-10:]
    rev = complementary_dna(last_10)[::-1]
    
    return {
        "template": template,
        "fwd_correct": fwd,
        "rev_correct": rev
    }

def generate_splicing_challenge() -> Dict:
    """Gera uma sequência com éxons (maiúsculas) e íntrons (minúsculas)."""
    exon1 = "ATGCGTAC"
    intron = "gtagctagctag"
    exon2 = "GTACTGAA"
    return {
        "pre_mrna": exon1 + intron + exon2,
        "exon1": exon1,
        "intron": intron,
        "exon2": exon2,
        "mature_mrna": exon1 + exon2
    }

BIOCHEM_QUIZ = [
    {
        "question": "Quais bases nitrogenadas são classificadas como Purinas (possuem dois anéis)?",
        "options": ["Adenina e Guanina", "Citosina e Timina", "Adenina e Uracila", "Guanina e Citosina"],
        "answer": "Adenina e Guanina",
        "explanation": "Purinas (A e G) têm estrutura de dois anéis fundidos. Pirimidinas (C, T e U) têm apenas um anel."
    },
    {
        "question": "Quantas pontes de hidrogênio estabilizam o par Guanina-Citosina (G-C)?",
        "options": ["1", "2", "3", "4"],
        "answer": "3",
        "explanation": "O par G-C possui 3 pontes de hidrogênio, tornando-o mais estável termicamente que o par A-T (2 pontes)."
    },
    {
        "question": "Qual base é exclusiva do RNA, substituindo a Timina do DNA?",
        "options": ["Adenina", "Citosina", "Guanina", "Uracila"],
        "answer": "Uracila",
        "explanation": "A Uracila (U) é a base pirimídica que se pareia com a Adenina no RNA."
    }
]

CODON_TABLE = {
    'T': {'T': {'T': 'Phe', 'C': 'Phe', 'A': 'Leu', 'G': 'Leu'},
          'C': {'T': 'Ser', 'C': 'Ser', 'A': 'Ser', 'G': 'Ser'},
          'A': {'T': 'Tyr', 'C': 'Tyr', 'A': 'STOP', 'G': 'STOP'},
          'G': {'T': 'Cys', 'C': 'Cys', 'A': 'STOP', 'G': 'Trp'}},
    'C': {'T': {'T': 'Leu', 'C': 'Leu', 'A': 'Leu', 'G': 'Leu'},
          'C': {'T': 'Pro', 'C': 'Pro', 'A': 'Pro', 'G': 'Pro'},
          'A': {'T': 'His', 'C': 'His', 'A': 'Gln', 'G': 'Gln'},
          'G': {'T': 'Arg', 'C': 'Arg', 'A': 'Arg', 'G': 'Arg'}},
    'A': {'T': {'T': 'Ile', 'C': 'Ile', 'A': 'Ile', 'G': 'Met'},
          'C': {'T': 'Thr', 'C': 'Thr', 'A': 'Thr', 'G': 'Thr'},
          'A': {'T': 'Asn', 'C': 'Asn', 'A': 'Lys', 'G': 'Lys'},
          'G': {'T': 'Ser', 'C': 'Ser', 'A': 'Arg', 'G': 'Arg'}},
    'G': {'T': {'T': 'Val', 'C': 'Val', 'A': 'Val', 'G': 'Val'},
          'C': {'T': 'Ala', 'C': 'Ala', 'A': 'Ala', 'G': 'Ala'},
          'A': {'T': 'Asp', 'C': 'Asp', 'A': 'Glu', 'G': 'Glu'},
          'G': {'T': 'Gly', 'C': 'Gly', 'A': 'Gly', 'G': 'Gly'}}
}
