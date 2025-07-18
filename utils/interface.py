import streamlit as st

LINE_SIZE = 12

def show_complementary_streamlit(comp_dna_seq: str, dna_seq: str) -> None:
    """
    Exibe a fita de DNA complementar (complementary()) de forma formatada no Streamlit.
    """
    st.subheader("SEQUÊNCIA DE DNA COMPLEMENTAR")

    for i in range(0, len(comp_dna_seq), LINE_SIZE):
        dna_seq_fragment = " ".join(dna_seq[i : i + LINE_SIZE])
        comp_dna_seq_fragment = " ".join(comp_dna_seq[i : i + LINE_SIZE])

        st.markdown(f"**{i:04}** &nbsp; 5' &nbsp; `{dna_seq_fragment}` &nbsp; 3'")
        st.markdown(f"**{i:04}** &nbsp; 3' &nbsp; `{comp_dna_seq_fragment}` &nbsp; 5'")
    st.markdown("---")


def show_transcription_streamlit(rna: str) -> None:
    """
    Exibe o transcrito de RNA (transcription()) de forma formatada no Streamlit.
    """
    st.subheader("SEQUÊNCIA DE RNA")

    for i in range(0, len(rna), LINE_SIZE):
        rna_seq = " ".join(rna[i : i + LINE_SIZE])
        st.markdown(f"**{i:04}** - 5' `{rna_seq}` 3'")
    st.markdown("---")


def show_translation_streamlit(protein: list) -> None:
    """
    Exibe a sequência de aminoácidos (translation()) de forma formatada no Streamlit.
    """
    st.subheader("SEQUÊNCIA DE AMINOÁCIDOS")

    colored_aa = []
    for aa in protein:

        if aa == "Met":
            colored_aa.append(f"<span style='color: yellow;'>{aa}</span>")
        else:
            colored_aa.append(f"<span style='color: magenta;'>{aa}</span>")

    for i in range(0, len(colored_aa), LINE_SIZE):
        aa_seq = " - ".join(colored_aa[i : i + LINE_SIZE])
        st.markdown(f"**{i:04}** - {aa_seq}", unsafe_allow_html=True)
    st.markdown("---")


def show_orfs_streamlit(orfs: list) -> None:
    """
    Exibe as ORFs identificadas (orf_finder()) de forma formatada no Streamlit.
    """
    st.subheader("OPEN READING FRAMES")

    if not orfs:
        st.info("Nenhum ORF encontrado.")
        return

    for index, orf in enumerate(orfs, 1):
        st.markdown(f"**ORF - {index:02}** - 5' `{orf}` 3'")
    st.markdown("---")


def show_gene_identifier_streamlit(gene_record) -> None:
    """
    Exibe o gene identificado (gene_identifier()) de forma formatada no Streamlit.
    """
    st.subheader("IDENTIFICADOR DE GENE")

    blast_record = gene_record

    if not blast_record.alignments:
        st.info("Nenhum gene associado à sequência foi encontrado.")
        return

    for alignment in blast_record.alignments[:3]:
        for hsp in alignment.hsps:
            st.markdown("---")
            st.write("**Possível gene encontrado:**")
            st.write(f"**Título:** {alignment.title}")
            st.write(f"**Identidade:** {hsp.identities}")
            st.write(f"**Sequência correspondente:** `{hsp.sbjct}`")
    st.markdown("---")


def welcome_streamlit() -> None:
    """
    Exibe o nome do programa a primeira vez que é inicializado no Streamlit.
    """

    st.title("GeneFlux")
    st.markdown("""

        #### DNA → RNA → ORFs → Proteína
        ###### by RMello
    """)
    st.markdown("---")
