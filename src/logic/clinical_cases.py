from typing import Dict, List

# Banco de dados de casos clínicos didáticos
CLINICAL_CASES = {
    "Anemia Falciforme (HBB)": {
        "description": """Paciente de 8 anos apresenta cansaço extremo e dores articulares. 
        O hemograma revela hemácias em formato de foice (drepanócitos). 
        Sua tarefa é analisar o fragmento do gene da Beta-Globina (HBB) para confirmar a mutação.""",
        "reference_dna": "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCC",
        "patient_dna": "ATGGTGCACCTGACTCCTGTGGAGAAGTCTGCC",
        "clinical_context": """A troca de um A por um T na posição 20 do DNA resulta na substituição de um 
        Ácido Glutâmico (Glu) por uma Valina (Val) na proteína. O Glu é hidrofílico (polar), enquanto a 
        Valina é hidrofóbica (apolar), o que causa a polimerização da hemoglobina em baixas tensões de oxigênio.""",
        "mutation_pos": 20,
        "ref_base": "A",
        "pat_base": "T"
    },
    "Talassemia Beta (Mutação de Ponto)": {
        "description": """Paciente apresenta anemia microcítica e hipocrômica grave. 
        A análise do gene HBB sugere uma mutação que cria um códon de parada prematuro (Nonsense).""",
        "reference_dna": "ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCC",
        "patient_dna": "ATGGTGCACCTGACTCCTGATGAGAAGTCTGCC",
        "clinical_context": """A mutação G -> T no códon 7 transforma o códon GAG (Glu) no códon de parada TAG (Stop). 
        Isso interrompe a síntese da proteína, resultando em uma globina incompleta e funcionalmente nula.""",
        "mutation_pos": 21,
        "ref_base": "G",
        "pat_base": "T"
    }
}

def compare_clinical_sequences(ref: str, pat: str) -> List[Dict]:
    """Identifica as diferenças exatas entre as sequências."""
    diffs = []
    min_len = min(len(ref), len(pat))
    for i in range(min_len):
        if ref[i] != pat[i]:
            diffs.append({
                "pos": i + 1,
                "ref": ref[i],
                "pat": pat[i]
            })
    return diffs
