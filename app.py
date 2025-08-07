import streamlit as st
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqUtils import seq3
import re
from typing import List, Tuple, Literal

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
    
    protein = str(Seq(rna_seq).translate(to_stop=False))
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
        with st.spinner("Enviando sequência para o NCBI BLAST. Isso pode levar até um minuto, por favor, aguarde..."):
            result_handle = NCBIWWW.qblast(
                program="blastn",
                database="nt",
                sequence=dna_seq,
            )
            # time.sleep(2) # Simula um delay para a resposta
        
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

# --- Funções de Exibição Formatada ---

LINE_SIZE = 12

def display_formatted_sequence(seq_type: str, seq: str, line_size: int = LINE_SIZE):
    """Exibe uma sequência de forma formatada."""
    st.subheader(f"{seq_type.upper()}")
    for i in range(0, len(seq), line_size):
        fragment = " ".join(seq[i : i + line_size])
        st.code(f"{i:04}  {fragment}")

def display_complementary(dna_seq: str, comp_dna_seq: str, line_size: int = LINE_SIZE):
    """Exibe o filamento de DNA e seu complementar."""
    st.subheader("SEQUÊNCIA DE DNA COMPLEMENTAR")
    for i in range(0, len(dna_seq), line_size):
        dna_fragment = " ".join(dna_seq[i : i + line_size])
        comp_dna_fragment = " ".join(comp_dna_seq[i : i + line_size])
        st.code(f"{i:04}  5'  {dna_fragment}  3'\n"
                f"{i:04}  3'  {comp_dna_fragment}  5'")

def display_amino_acids(protein_list: List[str], line_size: int = LINE_SIZE):
    """Exibe a sequência de aminoácidos."""
    st.subheader("SEQUÊNCIA DE AMINOÁCIDOS")
    for i in range(0, len(protein_list), line_size):
        aa_fragment = " - ".join(protein_list[i : i + line_size])
        st.code(f"{i:04} - {aa_fragment}")

# --- Estrutura da Aplicação Streamlit ---

# Configuração da página e título
st.set_page_config(
    page_title="GeneFlux: DNA Sequence Analyzer",
    layout="wide"
)

# Inicializa st.session_state
if "seq" not in st.session_state:
    st.session_state.seq = None
if "mol_type" not in st.session_state:
    st.session_state.mol_type = None

# Título do aplicativo
st.title("GeneFlux")
st.markdown("""

    #### DNA → RNA → ORFs → Proteína
    ###### by RMello
""")
st.markdown("---")

# Seção de entrada de dados
st.subheader("Entrada de Sequência")
user_sequence = st.text_area(
    "Digite sua sequência de DNA ou RNA aqui:",
    height=150,
    placeholder="Ex: ATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT"
).strip().upper().replace(" ", "")

# Adição de um botão para processar a sequência
if st.button("Processar Sequência"):
    if user_sequence:
        try:
            seq, mol_type = validate_sequence(user_sequence)
            st.session_state.seq = seq
            st.session_state.mol_type = mol_type
            st.rerun()  # Força a página a recarregar para exibir os resultados
        except ValueError as e:
            st.error(f"Erro: {e}")
    else:
        st.info("Por favor, digite uma sequência antes de processar.")

# Bloco de análise que só é exibido se houver uma sequência válida na sessão
if st.session_state.seq and st.session_state.mol_type:
    st.success(f"Sequência de **{st.session_state.mol_type.upper()}** detectada!")
    st.markdown("---")

    # Menu lateral para as ferramentas
    st.sidebar.title("Opções de Análise")
    option = st.sidebar.radio(
        "Escolha uma operação para realizar:",
        ("DNA Complementar", "Transcrever para RNA", "Traduzir para Aminoácidos", "ORF Finder", "Identificador de Gene")
    )

    # Exibição de resultados com base na opção do menu
    st.subheader(f"Resultado da Análise: {option}")

    if option == "DNA Complementar":
        comp_dna = complementary_dna(st.session_state.seq)
        display_complementary(st.session_state.seq, comp_dna)
        
    elif option == "Transcrever para RNA":
        rna_seq = transcription(st.session_state.seq)
        display_formatted_sequence("Sequência de RNA", rna_seq)

    elif option == "Traduzir para Aminoácidos":
        st.write("---")
        rna_seq = transcription(st.session_state.seq)
        orfs = orf_finder(st.session_state.seq)

        # Se houver ORFs, dá a opção de escolher qual traduzir
        if orfs:
            orf_options = [f"ORF {i+1} (Início: {rna_seq.find(orf)} -> Fim: {rna_seq.find(orf) + len(orf)})" for i, orf in enumerate(orfs)]
            orf_options.insert(0, "Traduzir a sequência completa")
            
            selected_option = st.selectbox(
                "Escolha a sequência para tradução:",
                options=orf_options,
                index=0
            )

            if selected_option == "Traduzir a sequência completa":
                seq_to_translate = rna_seq
            else:
                orf_index = orf_options.index(selected_option) - 1
                seq_to_translate = orfs[orf_index]

        else:
            st.warning("Nenhuma ORF encontrada. Apenas a tradução da sequência completa está disponível.")
            seq_to_translate = rna_seq

        protein_list = translation(seq_to_translate)
        if protein_list:
            display_amino_acids(protein_list)
            
    elif option == "ORF Finder":
        orfs = orf_finder(st.session_state.seq)
        if orfs:
            st.write("### Open Reading Frames (ORFs) Encontradas")
            for idx, orf in enumerate(orfs, 1):
                st.code(f"ORF - {idx:02}: {orf}")
        else:
            st.warning("Nenhuma ORF encontrada.")

    elif option == "Identificador de Gene":
        blast_record = gene_identifier(st.session_state.seq)
        if blast_record:
            st.write("### Possíveis Genes Encontrados (Top 3)")
            for alignment in blast_record.alignments[:3]:
                for hsp in alignment.hsps:
                    st.write("---")
                    st.success(f"**Gene Encontrado:** {alignment.title}")
                    st.info(f"**Identidade:** {hsp.identities}")
                    st.code(f"Sequência Correspondente: {hsp.sbjct}")

