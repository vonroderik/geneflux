import streamlit as st
import re
from typing import Tuple, Literal
from utils import dna_tools as gene
from utils import interface as gui

def validate_sequence(seq: str) -> Tuple[str, Literal["dna", "rna"]]:
    """
    Usa uma Regex para verificar se a entrada do usuário é uma sequência de DNA ou RNA.
    Também verifica se a sequência não é uma mistura de DNA e RNA.
    Sequências sem T/U são consideradas DNA.
    Retorna a sequência validada e uma string especificando o tipo de molécula.
    """
    if not seq:
        st.error("Entrada inválida. A sequência não pode ser vazia.")
        st.stop()

    match = re.fullmatch(r"^([ATCGU]+)$", seq)

    if not match:
        st.error("""
        Sequência inválida.
        Por favor, insira uma sequência de ácido nucleico válida.
        Ela deve conter apenas os caracteres A, T, C, G ou U.
        """)
        st.stop()

    validated_seq = match.group(1)

    if len(validated_seq) < 3:
        st.error("Sua sequência deve ter pelo menos 3 nucleotídeos.")
        st.stop()

    if "U" in validated_seq and "T" in validated_seq:
        st.error("""
        Sequência inválida: não é possível misturar bases de RNA e DNA.
        Sua sequência contém tanto timina (T) quanto uracila (U).
        """)
        st.stop()

    return validated_seq, "dna" if "U" not in validated_seq else "rna"

def main():
    gui.welcome_streamlit()

    if 'user_sequence_input' not in st.session_state:
        st.session_state['user_sequence_input'] = ""
    if 'current_sequence' not in st.session_state:
        st.session_state['current_sequence'] = ""
    if 'mol_type' not in st.session_state:
        st.session_state['mol_type'] = ""

    raw_user_input = st.text_input(
        "Por favor, digite uma sequência de DNA ou RNA:",
        value=st.session_state['user_sequence_input'],
        help="Ex: ATCGTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC",
        key="seq_input" # Adiciona uma chave para o widget
    )
    user_sequence_input = (raw_user_input or "").strip().upper().replace(" ", "")


    if user_sequence_input and user_sequence_input != st.session_state['user_sequence_input']:
        st.session_state['user_sequence_input'] = user_sequence_input
        try:
            seq, mol_type = validate_sequence(user_sequence_input)
            st.session_state['current_sequence'] = seq
            st.session_state['mol_type'] = mol_type
            st.success(f"Sequência de {mol_type.upper()} detectada!")
            for key in ['show_complementary', 'show_transcription', 'show_translation', 'show_orfs', 'show_gene_identifier']:
                if key in st.session_state:
                    st.session_state[key] = False
        except ValueError as e:
            st.error(f"Erro: {e}")
            st.session_state['current_sequence'] = ""
            st.session_state['mol_type'] = ""
            return
        except Exception as e:
            st.error(f"Ocorreu um erro inesperado: {e}")
            st.session_state['current_sequence'] = ""
            st.session_state['mol_type'] = ""
            return

    if st.session_state['current_sequence']:
        st.sidebar.header("Opções")
        choice = st.sidebar.radio(
            "Escolha uma operação:",
            (
                "Gerar fita de DNA complementar",
                "Transcrever",
                "Traduzir",
                "Encontrar ORF",
                "Identificar Gene",
                "Digitar nova sequência"
            ),
            key="operation_choice"
        )

        seq_to_process = st.session_state['current_sequence']

        for key in ['show_complementary', 'show_transcription', 'show_translation', 'show_orfs', 'show_gene_identifier']:
            if key in st.session_state:
                st.session_state[key] = False

        if choice == "Gerar fita de DNA complementar":
            gene.complementary_dna_streamlit(seq_to_process)
        elif choice == "Transcrever":
            gene.transcription_streamlit(seq_to_process)
        elif choice == "Traduzir":
            gene.translation_streamlit(seq_to_process)
        elif choice == "Encontrar ORF":
            gene.orf_finder_streamlit(seq_to_process)
        elif choice == "Identificar Gene":
            gene.gene_identifier_streamlit(seq_to_process)
        elif choice == "Digitar nova sequência":
            st.session_state.clear() 
            st.rerun()

    if st.session_state.get('show_complementary'):
        gui.show_complementary_streamlit(st.session_state['comp_seq'], st.session_state['dna_seq_original'])
        st.session_state['show_complementary'] = False

    if st.session_state.get('show_transcription'):
        gui.show_transcription_streamlit(st.session_state['rna_seq'])
        st.session_state['show_transcription'] = False

    if st.session_state.get('show_translation'):
        gui.show_translation_streamlit(st.session_state['aa_list'])
        st.session_state['show_translation'] = False

    if st.session_state.get('show_orfs'):
        gui.show_orfs_streamlit(st.session_state['orfs'])
        st.session_state['show_orfs'] = False

    if st.session_state.get('show_gene_identifier'):
        gui.show_gene_identifier_streamlit(st.session_state['blast_record'])
        st.session_state['show_gene_identifier'] = False

if __name__ == "__main__":
    main()
