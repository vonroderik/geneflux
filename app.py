import streamlit as st
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqUtils import seq3
import re
import random
from typing import List, Tuple, Literal
from collections import Counter
from itertools import product
from io import StringIO

# --- Funções de Análise do DNA, adaptadas para Streamlit ---


def validate_sequence(seq: str) -> Tuple[str, Literal["dna", "rna"]]:
    """Valida se a entrada é uma sequência de DNA ou RNA."""
    if not seq:
        raise ValueError("Entrada vazia. Por favor, digite uma sequência de DNA ou RNA")

    match = re.fullmatch(r"^([ATCGU]+)$", seq.upper())
    if not match:
        raise ValueError(
            """
            Sequência inválida.
            Por favor, digite uma sequência de ácido nucleico válida.
            Deve conter apenas os caracteres A, T, C, G ou U.
            """
        )

    validated_seq = match.group(1)
    if len(validated_seq) < 3:
        raise ValueError("Sua sequência deve ter pelo menos 3 nucleotídeos")

    if "U" in validated_seq and "T" in validated_seq:
        raise ValueError(
            """
            Sequência de DNA/RNA inválida: não é possível misturar bases de RNA e DNA.
            Sua sequência contém tanto timina (T) quanto uracila (U)
            """
        )

    return validated_seq, "dna" if "U" not in validated_seq else "rna"


def validate_genotype(genotype: str):
    """Verifica se o genótipo entrado pelo usuário é válido."""
    valid_genotypes = {
        "AA",
        "Aa",
        "aA",
        "aa",
        "BB",
        "Bb",
        "bB",
        "bb",
        "OO",
        "Oo",
        "oO",
        "oo",
        "AB",
        "BA",
        "Ao",
        "AO",
        "OA",
        "oA",
        "aO",
        "Oa",
        "Bo",
        "BO",
        "OBoB",
        "bO",
        "Ob",
    }
    genotype_normalized = genotype.strip()
    if not genotype_normalized:
        raise ValueError("Entrada vazia. Por favor, digite um genótipo válido.")
    if genotype_normalized not in valid_genotypes:
        raise ValueError(f"'{genotype}' não é um genótipo válido.")


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


def translation(seq: str) -> List[str]:
    """Traduz a sequência de RNA para aminoácidos."""
    rna_seq = seq
    if len(rna_seq) < 3:
        st.error("Sua sequência deve ter pelo menos 3 códons para ser traduzida.")
        return []

    trimmed_seq = rna_seq[: len(rna_seq) - len(rna_seq) % 3]
    protein = str(Seq(trimmed_seq).translate(to_stop=False))
    return [seq3(aa).title() for aa in protein]


def orf_finder(seq: str) -> List[str] | None:
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


def gene_identifier(seq: str):
    """Envia uma sequência para o NCBI BLAST para identificação de genes."""
    dna_seq = reverse_transcription(seq)
    if len(seq) < 11:
        st.error("Sua sequência deve ter pelo menos 11 nucleotídeos para o BLAST.")
        return None

    result_handle = None
    try:
        with st.spinner(
            "Enviando sequência para o NCBI BLAST. Isso pode levar até um minuto, por favor, aguarde..."
        ):
            result_handle = NCBIWWW.qblast(
                program="blastn",
                database="nt",
                sequence=dna_seq,
            )

        blast_record = NCBIXML.read(result_handle)

        if not blast_record.alignments:
            st.error("Nenhum gene associado à sequência de DNA foi encontrado")
            return None

        return blast_record
    except Exception as e:
        st.error(f"A solicitação BLAST falhou: {e}")
        st.warning("Verifique sua conexão ou tente novamente mais tarde.")
        return None
    finally:
        if result_handle is not None:
            result_handle.close()


# --- Funções do Módulo Educativo, adaptadas para Streamlit ---


def generate_sequence() -> str:
    """Gera uma sequência de DNA aleatória que começa com ATG e termina com um códon de parada."""
    STOP_CODONS = ["TGA", "TAA", "TAG"]
    random_middle = "".join(random.choices(["C", "G", "T", "A"], k=12))
    return "ATG" + random_middle + random.choice(STOP_CODONS)


def get_colored_feedback(correct_seq: str, user_seq: str) -> str:
    """
    Compara duas sequências e retorna uma string formatada com cores para destacar erros.
    - Verde para bases corretas.
    - Vermelho para bases incorretas.
    """
    feedback_str = ""
    min_len = min(len(correct_seq), len(user_seq))

    for i in range(min_len):
        if correct_seq[i] == user_seq[i]:
            feedback_str += (
                f"<span style='color:green; font-weight:bold;'>{user_seq[i]}</span>"
            )
        else:
            feedback_str += (
                f"<span style='color:red; font-weight:bold;'>{user_seq[i]}</span>"
            )

    if len(user_seq) > min_len:
        feedback_str += (
            f"<span style='color:red; font-weight:bold;'>{user_seq[min_len:]}</span>"
        )

    return feedback_str


def get_mutation_type(original_dna: str, mutated_dna: str) -> str:
    """Determina o tipo de mutação com base nas sequências de DNA original e mutada."""
    if len(original_dna) == len(mutated_dna):
        mismatches = sum(1 for a, b in zip(original_dna, mutated_dna) if a != b)
        if mismatches == 1:
            return "Substituição"

    if len(mutated_dna) > len(original_dna):
        return "Inserção"

    if len(mutated_dna) < len(original_dna):
        return "Deleção"

    return "N/A"


def punnet_square_display(genotype1: str, genotype2: str):
    """Cria e exibe um quadro de Punnett de forma mais visualmente agradável."""
    try:
        validate_genotype(genotype1)
        validate_genotype(genotype2)
    except ValueError as e:
        st.error(f"Erro: {e}")
        return

    gamete1 = list(genotype1)
    gamete2 = list(genotype2)
    offspring = [a + b for a, b in product(gamete1, gamete2)]

    is_mendelian_cross = not any(
        allele.lower() in {"b", "o"} for allele in genotype1 + genotype2
    )

    normalized_offspring = []
    for gen in offspring:
        if is_mendelian_cross:
            normalized_offspring.append("".join(sorted(gen)))
        else:
            normalized_gen_upper = "".join(sorted(gen.upper()))
            if normalized_gen_upper == "AO":
                normalized_offspring.append("Ao")
            elif normalized_gen_upper == "BO":
                normalized_offspring.append("Bo")
            elif normalized_gen_upper == "OO":
                normalized_offspring.append("oo")
            else:
                normalized_offspring.append(normalized_gen_upper)

    count = Counter(normalized_offspring)

    st.subheader("Quadro de Punnett")
    st.markdown("---")

    table_html = "<table style='width:100%; border-collapse: collapse;'>"
    table_html += "<tr>"
    table_html += "<th style='border: 1px solid black; padding: 8px;'>♂ / ♀</th>"
    for g in gamete2:
        table_html += (
            f"<th style='border: 1px solid black; padding: 8px;'><b>{g}</b></th>"
        )
    table_html += "</tr>"

    offspring_index = 0
    for g1 in gamete1:
        table_html += "<tr>"
        table_html += (
            f"<td style='border: 1px black solid; padding: 8px;'><b>{g1}</b></td>"
        )
        for _ in gamete2:
            child = normalized_offspring[offspring_index]

            if "aa" in child.lower() or "oo" in child.lower():
                color = "red"
            else:
                color = "green"

            table_html += f"<td style='border: 1px solid black; padding: 8px;'><span style='color:{color}; font-weight:bold;'>{child}</span></td>"
            offspring_index += 1
        table_html += "</tr>"

    table_html += "</table>"
    st.markdown(table_html, unsafe_allow_html=True)
    st.markdown("---")

    st.markdown("### Frequência Genotípica")
    for genotype, freq in count.items():
        perc = 100 * freq / 4
        st.write(f"- **{genotype}**: {freq}/{4} ({perc:.1f}%)")

    if is_mendelian_cross:
        st.markdown("### Frequência Fenotípica (Autossômica)")
        dominant_count = sum(freq for gen, freq in count.items() if "A" in gen)
        recessive_count = count.get("aa", 0)
        total = dominant_count + recessive_count

        if total > 0:
            perc_dom = 100 * dominant_count / total
            perc_rec = 100 * recessive_count / total
            st.write(f"- **Dominante**: {dominant_count}/{total} ({perc_dom:.1f}%)")
            st.write(f"- **Recessivo**: {recessive_count}/{total} ({perc_rec:.1f}%)")

    else:
        st.markdown("### Frequência Fenotípica (Sistema ABO)")
        phenotype_counts = Counter()
        for genotype, freq in count.items():
            if "AB" in genotype.upper():
                phenotype_counts["AB"] += freq
            elif "A" in genotype.upper():
                phenotype_counts["A"] += freq
            elif "B" in genotype.upper():
                phenotype_counts["B"] += freq
            elif "OO" in genotype.upper():
                phenotype_counts["O"] += freq

        total = sum(phenotype_counts.values())
        if total > 0:
            for phenotype, freq in phenotype_counts.items():
                perc = 100 * freq / total
                st.write(f"- **Tipo {phenotype}**: {freq}/{total} ({perc:.1f}%)")


def display_formatted_sequence(seq_type: str, seq: str, line_size: int = 12):
    """Exibe uma sequência de forma formatada."""
    st.subheader(f"{seq_type.upper()}")
    for i in range(0, len(seq), line_size):
        fragment = " ".join(seq[i : i + line_size])
        st.code(f"{i:04}  {fragment}")


def display_complementary(dna_seq: str, comp_dna_seq: str):
    """Exibe o filamento de DNA e seu complementar."""
    st.subheader("SEQUÊNCIA DE DNA COMPLEMENTAR")
    for i in range(0, len(dna_seq), 12):
        dna_fragment = " ".join(dna_seq[i : i + 12])
        comp_dna_fragment = " ".join(comp_dna_seq[i : i + 12])
        st.code(f"{i:04}  5'  {dna_fragment}  3'\n{i:04}  3'  {comp_dna_fragment}  5'")


def display_amino_acids(protein_list: List[str], line_size: int = 12):
    """Exibe a sequência de aminoácidos."""
    st.subheader("SEQUÊNCIA DE AMINOÁCIDOS")
    output = StringIO()
    for i in range(0, len(protein_list), line_size):
        aa_fragment = " - ".join(protein_list[i : i + line_size])
        output.write(f"{i:04} - {aa_fragment}\n")
    st.code(output.getvalue())


def get_colored_codons(original_seq: str, mutated_seq: str) -> str:
    feedback_str = ""
    min_len = min(len(original_seq), len(mutated_seq))

    first_diff = -1
    for i in range(min_len):
        if original_seq[i] != mutated_seq[i]:
            first_diff = i
            break
    if first_diff == -1 and len(original_seq) != len(mutated_seq):
        first_diff = min_len

    for i in range(0, len(mutated_seq), 3):
        codon = mutated_seq[i : i + 3]
        if first_diff != -1 and i >= first_diff:
            feedback_str += (
                f"<span style='color:red; font-weight:bold;'>{codon}</span> "
            )
        else:
            feedback_str += (
                f"<span style='color:black; font-weight:normal;'>{codon}</span> "
            )

    return feedback_str.strip()


def display_mutation_results(original_seq: str, mutated_seq: str, mutation_type: str):
    original_rna = transcription(original_seq)
    mutated_rna = transcription(mutated_seq)

    original_protein = translation(original_rna)
    mutated_protein = translation(mutated_rna)

    protein_mutation_type = ""
    STOP_CODONS_RNA = ["UAA", "UAG", "UGA"]

    if mutation_type == "Substituição":
        mutated_codon_index = -1
        for i in range(len(original_rna) // 3):
            if original_rna[i * 3 : i * 3 + 3] != mutated_rna[i * 3 : i * 3 + 3]:
                mutated_codon_index = i
                break

        mutated_codon = mutated_rna[
            mutated_codon_index * 3 : mutated_codon_index * 3 + 3
        ]

        if mutated_codon in STOP_CODONS_RNA:
            protein_mutation_type = "Nonsense (Sem Sentido)"
            mutated_protein = original_protein[:mutated_codon_index] + ["STOP"]
        elif "".join(original_protein) == "".join(mutated_protein):
            protein_mutation_type = "Silenciosa"
        else:
            protein_mutation_type = "Missense (Sentido Trocado)"

    elif mutation_type in ["Inserção", "Deleção"]:
        if len(mutated_seq) % 3 != 0:
            protein_mutation_type = "Frameshift (Alteração de Quadro de Leitura)"
        else:
            if len(original_protein) == len(mutated_protein):
                protein_mutation_type = "Silenciosa"
            else:
                protein_mutation_type = "Deleção de Códons (não frameshift)"

    st.subheader("Resultados da Mutação")
    st.markdown("---")

    st.write(f"**Tipo de Mutação no DNA:** **{mutation_type}**")
    st.write(f"**Tipo de Mutação na Proteína:** **{protein_mutation_type}**")

    st.markdown("### Sequência de DNA")
    original_codons_dna = " ".join(
        [original_seq[i : i + 3] for i in range(0, len(original_seq), 3)]
    )
    st.code(f"5'  {original_codons_dna}  3'")
    mutated_codons_dna = get_colored_codons(original_seq, mutated_seq)
    st.markdown(
        f"<p style='font-family:monospace; font-size:14px;'>5'  {mutated_codons_dna}  3'</p>",
        unsafe_allow_html=True,
    )
    st.markdown("---")

    st.markdown("### Sequência de RNA Mensageiro (mRNA)")
    original_codons_rna = " ".join(
        [original_rna[i : i + 3] for i in range(0, len(original_rna), 3)]
    )
    st.code(f"5'  {original_codons_rna}  3'")
    mutated_codons_rna = get_colored_codons(original_rna, mutated_rna)
    st.markdown(
        f"<p style='font-family:monospace; font-size:14px;'>5'  {mutated_codons_rna}  3'</p>",
        unsafe_allow_html=True,
    )
    st.markdown("---")

    st.markdown("### Sequência de Proteína")
    st.write(f"**Original:** `{' - '.join(original_protein)}`")
    st.write(f"**Mutada:** `{' - '.join(mutated_protein)}`")
    st.markdown("---")


# --- Estrutura da Aplicação Streamlit ---

st.set_page_config(page_title="GeneFlux: DNA Sequence Analyzer", layout="wide")

st.title("GeneFlux")
st.markdown("""
    #### DNA → RNA → ORFs → Proteína
    ###### por RMello
""")
st.markdown("---")

module_choice = st.sidebar.radio(
    "Escolha um módulo:", ("Módulo de Análises", "Módulo Educativo")
)

if module_choice == "Módulo de Análises":
    st.subheader("Módulo de Análises")

    if "seq" not in st.session_state:
        st.session_state.seq = None
    if "mol_type" not in st.session_state:
        st.session_state.mol_type = None

    user_sequence = (
        st.text_area(
            "Digite sua sequência de DNA ou RNA aqui:",
            height=150,
            placeholder="Ex: ATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT",
        )
        .strip()
        .upper()
        .replace(" ", "")
    )

    if st.button("Processar Sequência"):
        if user_sequence:
            try:
                seq, mol_type = validate_sequence(user_sequence)
                st.session_state.seq = seq
                st.session_state.mol_type = mol_type
                st.rerun()
            except ValueError as e:
                st.error(f"Erro: {e}")
        else:
            st.info("Por favor, digite uma sequência antes de processar.")

    if st.session_state.seq and st.session_state.mol_type:
        st.success(f"Sequência de **{st.session_state.mol_type.upper()}** detectada!")
        st.markdown("---")

        st.sidebar.markdown("### Ferramentas de Análise")
        option = st.sidebar.radio(
            "Escolha uma operação para realizar:",
            (
                "DNA Complementar",
                "Transcrever para RNA",
                "Traduzir para Aminoácidos",
                "Localizador de ORF",
                "Identificador de Gene",
            ),
        )

        st.subheader(f"Resultado da Análise: {option}")

        if option == "DNA Complementar":
            comp_dna = complementary_dna(st.session_state.seq)
            display_complementary(st.session_state.seq, comp_dna)
        elif option == "Transcrever para RNA":
            rna_seq = transcription(st.session_state.seq)
            display_formatted_sequence("Sequência de RNA", rna_seq)
        elif option == "Traduzir para Aminoácidos":
            rna_seq = transcription(st.session_state.seq)
            orfs = orf_finder(st.session_state.seq)

            if orfs:
                orf_options = [f"ORF {i + 1}" for i in range(len(orfs))]
                orf_options.insert(0, "Sequência Completa")

                selected_option = st.selectbox(
                    "Escolha a sequência para tradução:", options=orf_options, index=0
                )

                if selected_option == "Sequência Completa":
                    seq_to_translate = rna_seq
                else:
                    orf_index = orf_options.index(selected_option) - 1
                    seq_to_translate = orfs[orf_index]

            else:
                st.warning(
                    "Nenhuma ORF encontrada. Apenas a tradução da sequência completa está disponível."
                )
                seq_to_translate = rna_seq

            protein_list = translation(seq_to_translate)
            if protein_list:
                display_amino_acids(protein_list)

        elif option == "Localizador de ORF":
            orfs = orf_finder(st.session_state.seq)
            st.subheader("Open Reading Frames (ORFs) Encontradas")
            if not orfs:
                st.warning("Nenhuma ORF encontrada.")
            else:
                for idx, orf in enumerate(orfs, 1):
                    st.code(f"ORF - {idx:02}: {orf}")

        elif option == "Identificador de Gene":
            blast_record = gene_identifier(st.session_state.seq)
            if blast_record:
                st.subheader("Possíveis Genes Encontrados (Top 3)")
                for alignment in blast_record.alignments[:3]:
                    for hsp in alignment.hsps:
                        st.markdown("---")
                        st.success(f"**Gene Encontrado:** {alignment.title}")
                        st.info(f"**Identidade:** {hsp.identities}")
                        st.code(f"Sequência Correspondente: {hsp.sbjct}")

elif module_choice == "Módulo Educativo":
    st.subheader("Módulo Educativo")

    edu_option = st.sidebar.radio(
        "Escolha uma opção:",
        (
            "Herança Mendeliana",
            "Simulador de Mutação Gênica",
            "Exercício - DNA Complementar",
            "Exercício - Transcrição",
        ),
    )

    if "punnett_genotype1" not in st.session_state:
        st.session_state.punnett_genotype1 = ""
    if "punnett_genotype2" not in st.session_state:
        st.session_state.punnett_genotype2 = ""
    if "show_punnett_results" not in st.session_state:
        st.session_state.show_punnett_results = False
    if "mutation_sequence" not in st.session_state:
        st.session_state.mutation_sequence = ""
    if "exercise_state" not in st.session_state:
        st.session_state.exercise_state = "initial"

    if edu_option == "Herança Mendeliana":
        st.subheader("Herança Mendeliana")
        st.markdown("---")
        st.write(
            """Por favor, digite os genótipos de dois progenitores para gerar um Quadro de Punnett. 
            No caso de sistema ABO, utilize letras maiusculas"""
        )

        col1, col2 = st.columns(2)
        with col1:
            genotype1_input = st.text_input(
                "Genótipo Progenitor 1:",
                value=st.session_state.punnett_genotype1,
                placeholder="Ex: Aa, Ao, AB",
                key="g1_input",
            )
        with col2:
            genotype2_input = st.text_input(
                "Genótipo Progenitor 2:",
                value=st.session_state.punnett_genotype2,
                placeholder="Ex: aa, Bo, oo",
                key="g2_input",
            )

        if st.button("Gerar Quadro de Punnett"):
            try:
                validate_genotype(genotype1_input)
                validate_genotype(genotype2_input)
                st.session_state.punnett_genotype1 = genotype1_input
                st.session_state.punnett_genotype2 = genotype2_input
                st.session_state.show_punnett_results = True
            except ValueError as e:
                st.error(f"Erro: {e}")
                st.session_state.show_punnett_results = False

        if st.session_state.show_punnett_results:
            punnet_square_display(
                st.session_state.punnett_genotype1, st.session_state.punnett_genotype2
            )

    elif edu_option == "Simulador de Mutação Gênica":
        st.subheader("Simulador de Mutação Gênica")
        st.markdown("---")
        st.write("Crie uma mutação e veja o impacto em nível de DNA, RNA e Proteína.")

        if (
            st.button("Gerar Sequência de DNA")
            or not st.session_state.mutation_sequence
        ):
            st.session_state.mutation_sequence = generate_sequence()

        st.markdown(
            f"**Sequência de DNA Original:** `{st.session_state.mutation_sequence}`"
        )

        mutation_type_choice = st.selectbox(
            "Escolha o tipo de mutação:", ("Substituição", "Inserção", "Deleção")
        )

        st.markdown("---")

        mutated_seq = ""

        if mutation_type_choice == "Substituição":
            col_pos, col_base = st.columns(2)
            with col_pos:
                position = (
                    st.number_input(
                        "Posição da Mutação (1 a 18):",
                        min_value=1,
                        max_value=18,
                        value=1,
                    )
                    - 1
                )
            with col_base:
                new_base = st.selectbox("Nova Base:", ["A", "T", "C", "G"])

            dna_list = list(st.session_state.mutation_sequence)
            if 0 <= position < len(dna_list):
                dna_list[position] = new_base
            mutated_seq = "".join(dna_list)

        elif mutation_type_choice == "Inserção":
            col_pos, col_base = st.columns(2)
            with col_pos:
                position = (
                    st.number_input(
                        "Posição da Inserção (1 a 18):",
                        min_value=1,
                        max_value=18,
                        value=1,
                    )
                    - 1
                )
            with col_base:
                new_base = st.selectbox("Base para Inserir:", ["A", "T", "C", "G"])

            dna_list = list(st.session_state.mutation_sequence)
            dna_list.insert(position, new_base)
            mutated_seq = "".join(dna_list)

        elif mutation_type_choice == "Deleção":
            position = (
                st.number_input(
                    "Posição da Deleção (1 a 18):", min_value=1, max_value=18, value=1
                )
                - 1
            )

            dna_list = list(st.session_state.mutation_sequence)
            if 0 <= position < len(dna_list):
                del dna_list[position]
            mutated_seq = "".join(dna_list)

        if st.button("Aplicar Mutação"):
            mutation_type = get_mutation_type(
                st.session_state.mutation_sequence, mutated_seq
            )
            display_mutation_results(
                st.session_state.mutation_sequence, mutated_seq, mutation_type
            )

    elif edu_option.startswith("Exercício"):
        st.subheader(edu_option)
        st.markdown("---")

        exercise_type_map = {
            "Exercício - DNA Complementar": "DNA complementar",
            "Exercício - Transcrição": "RNA",
        }
        exercise_type = exercise_type_map[edu_option]

        st.write(
            f"Instruções: Digite a sequência de **{exercise_type}** para cada uma das sequências exibidas."
        )

        if (
            st.button("Iniciar Novo Exercício")
            or st.session_state.exercise_state == "initial"
        ):
            st.session_state.exercise_state = "running"
            st.session_state.show_exercise_results = False
            st.session_state.exercises = []

            for _ in range(3):
                exercise_sequence = generate_sequence()

                if exercise_type == "DNA complementar":
                    expected_sequence = complementary_dna(exercise_sequence)
                    valid_bases = "CGTA"
                elif exercise_type == "RNA":
                    expected_sequence = transcription(exercise_sequence)
                    valid_bases = "CGUA"

                st.session_state.exercises.append(
                    {
                        "original_sequence": exercise_sequence,
                        "expected_sequence": expected_sequence,
                        "valid_bases": valid_bases,
                    }
                )

            st.session_state.user_inputs = [""] * 3

        if st.session_state.exercise_state == "running":
            st.markdown("---")
            user_inputs = []
            for n in range(3):
                current_ex = st.session_state.exercises[n]
                st.subheader(f"Sequência {n + 1}")
                st.code(f"5' - {current_ex['original_sequence']} - 3'")

                user_input = st.text_input(
                    f"Sua resposta para {exercise_type} (use apenas {current_ex['valid_bases']}):",
                    value=st.session_state.user_inputs[n],
                    key=f"user_input_{n}",
                    placeholder=f"Digite a sequência de {exercise_type}",
                )
                st.session_state.user_inputs[n] = user_input

            if st.button("Verificar Respostas"):
                score = 0
                st.session_state.results = []
                for n in range(3):
                    current_ex = st.session_state.exercises[n]
                    student_sequence = (
                        st.session_state.user_inputs[n].strip().upper().replace(" ", "")
                    )
                    is_correct = False

                    if all(
                        base in current_ex["valid_bases"] for base in student_sequence
                    ) and len(student_sequence) == len(current_ex["expected_sequence"]):
                        if student_sequence == current_ex["expected_sequence"]:
                            is_correct = True

                    st.session_state.results.append(
                        {
                            "original_sequence": current_ex["original_sequence"],
                            "expected_sequence": current_ex["expected_sequence"],
                            "user_answer": student_sequence,
                            "correct": is_correct,
                        }
                    )
                    if is_correct:
                        score += 1

                st.session_state.final_score = score
                st.session_state.exercise_state = "finished"
                st.rerun()

        if st.session_state.exercise_state == "finished":
            st.subheader("Resultados Finais")
            total_score = st.session_state.final_score
            st.markdown(f"### Pontuação total: **{total_score}** de **3**")
            st.markdown("---")
            for i, result in enumerate(st.session_state.results, 1):
                st.write(f"**Exercício {i}**")

                st.markdown(f"**Resposta Correta:** `{result['expected_sequence']}`")

                colored_feedback = get_colored_feedback(
                    result["expected_sequence"], result["user_answer"]
                )
                st.markdown(
                    f"**Sua Resposta:** {colored_feedback}", unsafe_allow_html=True
                )

                if result["correct"]:
                    st.success("Correto")
                else:
                    st.error("Incorreto.")
                st.markdown("---")
