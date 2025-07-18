from typing import List
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqUtils import seq3
import streamlit as st

####### DNA Tools Functions #######

# Retorna a fita de DNA complementar (DNA -> DNA Complementar)
def complementary_dna_streamlit(seq: str) -> str:
    dna_seq = reverse_transcription(seq)
    comp_seq = dna_seq.translate(str.maketrans("ATCG", "TAGC"))
    st.session_state['comp_seq'] = comp_seq
    st.session_state['dna_seq_original'] = dna_seq
    st.session_state['show_complementary'] = True
    return comp_seq

# Retorna um transcrito de RNA (DNA -> RNA)
def transcription_streamlit(seq: str) -> str:
    rna_seq = seq.replace("T", "U")
    st.session_state['rna_seq'] = rna_seq
    st.session_state['show_transcription'] = True
    return rna_seq

# Usa o Seq do Biopython que converte a sequência em um objeto
# e retorna uma lista de aminoácidos (aa) de uma letra, e o Seq3 do Biopython.
# que converte o aa de uma letra para um formato de três letras.
# Chama sequence_to_translate() a sequência a ser traduzida: sequência completa ou
# ORF específica.
def translation_streamlit(seq: str):
    rna_seq = transcription(seq, show_results=False)
    if len(rna_seq) < 3:
        st.error("Sua sequência deve ter pelo menos 3 códons. Por favor, digite uma nova sequência.")
        return

    translation_option = st.radio(
        "Selecione uma opção para tradução:",
        ("Traduzir a sequência completa (ignorando Start e Stop códons)", "Selecionar uma ORF para traduzir")
    )

    selected_seq = None
    if translation_option == "Selecionar uma ORF para traduzir":
        orfs_list = orf_finder(rna_seq, show_results=False, wait=False)
        if orfs_list:
            orf_options = [f"ORF {i+1}: {orf}" for i, orf in enumerate(orfs_list)]
            selected_orf_str = st.selectbox("Selecione uma ORF para traduzir:", orf_options)
            chosen_orf_index = orf_options.index(selected_orf_str)
            selected_seq = orfs_list[chosen_orf_index]
        else:
            st.info("Nenhum ORF encontrado para traduzir.")
            return

    if selected_seq is None:
        rna_seq_to_translate = rna_seq[: len(rna_seq) - len(rna_seq) % 3]
    else:
        rna_seq_to_translate = selected_seq[: len(selected_seq) - len(selected_seq) % 3]

    protein = str(Seq(rna_seq_to_translate).translate(to_stop=False))
    aa_list = [seq3(aa).title() for aa in protein]

    st.session_state['aa_list'] = aa_list
    st.session_state['show_translation'] = True
    return aa_list


# Usa laços while e for para procurar por ORFs.
# Ele pesquisa os três quadros de RNA (+1, +2, +3) por um códon de início (AUG).
# Qualquer coisa entre um AUG e um códon de parada (UGA, UAG, UAA) é considerada uma ORF.
# A última ORF não precisa ter um códon de parada.
def orf_finder_streamlit(seq: str) -> List[str] | None:
    rna_seq = transcription(seq, show_results=False)
    if len(rna_seq) < 3:
        st.error("Sua sequência deve ter pelo menos 3 nucleotídeos. Por favor, digite uma nova sequência.")
        return

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
    st.session_state['orfs'] = orfs
    st.session_state['show_orfs'] = True
    return orfs


# Usa o NCBIWWW do Biopython para enviar uma consulta de pesquisa ao NCBI para BLAST,
# e recebe um arquivo XML que é convertido em um objeto de Classe por NCBIXML
def gene_identifier_streamlit(seq: str):
    dna_seq = reverse_transcription(seq)

    if len(seq) < 11:
        st.error("Sua sequência deve ter pelo menos 11 nucleotídeos para o BLAST. Por favor, digite uma nova sequência.")
        return

    st.info("Enviando sequência para NCBI BLAST. Isto pode levar até dois minutos, por favor aguarde...")

    try:
        result_handle = NCBIWWW.qblast(
            program="blastn",
            database="nt",
            sequence=dna_seq,
        )

        blast_record = NCBIXML.read(result_handle)

        if not blast_record.alignments:
            st.info("Nenhum gene associado à sequência foi encontrado.")
            return None

        st.session_state['blast_record'] = blast_record
        st.session_state['show_gene_identifier'] = True
        return blast_record

    except Exception as e:
        st.error(f"Requisição de BLAST falhou: {e}")
        st.error("Verifique sua conexão e/ou tente novamente mais tarde.")
        return

    finally:
        if 'result_handle' in locals() and result_handle is not None:
            result_handle.close()


# Retorna uma sequência de DNA a partir de uma sequência de RNA (RNA -> DNA)
def reverse_transcription(seq: str) -> str:
    return seq if "U" not in seq else seq.replace("U", "T")

def transcription(seq: str, show_results: bool) -> str:
    rna_seq = seq.replace("T", "U")
    return rna_seq

def orf_finder(seq: str, show_results: bool, wait: bool) -> List[str] | None:
    rna_seq = transcription(seq, show_results=False)
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
