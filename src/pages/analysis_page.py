import streamlit as st
from ..logic.dna_rna import (
    validate_sequence, 
    complementary_dna, 
    transcription, 
    translation, 
    orf_finder,
    calculate_gc_content,
    calculate_molecular_weight
)
from ..logic.ncbi import gene_identifier
from ..components.visualizations import (
    display_formatted_sequence, 
    display_complementary, 
    display_amino_acids
)

def render_analysis_page():
    st.subheader("Módulo de Análises")

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
                st.success(f"Sequência de **{st.session_state.mol_type.upper()}** detectada!")
            except ValueError as e:
                st.error(f"Erro: {e}")
        else:
            st.info("Por favor, digite uma sequência antes de processar.")

    if "seq" in st.session_state and st.session_state.seq:
        st.markdown("---")
        
        # New Metrics Row
        col1, col2, col3 = st.columns(3)
        with col1:
            gc = calculate_gc_content(st.session_state.seq)
            st.metric("Conteúdo GC", f"{gc:.2f}%")
        with col2:
            mw = calculate_molecular_weight(st.session_state.seq, st.session_state.mol_type)
            st.metric("Peso Molecular", f"{mw:.2f} Da")
        with col3:
            st.metric("Comprimento", f"{len(st.session_state.seq)} bp")

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
                    orf_index = int(selected_option.split()[-1]) - 1
                    seq_to_translate = orfs[orf_index]
            else:
                st.warning("Nenhuma ORF encontrada. Traduzindo sequência completa.")
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
                    st.code(f"ORF - {idx:02} (Size: {len(orf)}nt): {orf}")

        elif option == "Identificador de Gene":
            try:
                with st.spinner("Consultando NCBI BLAST..."):
                    blast_record = gene_identifier(st.session_state.seq)
                
                if blast_record and blast_record.alignments:
                    st.subheader("Possíveis Genes Encontrados (Top 3)")
                    for alignment in blast_record.alignments[:3]:
                        for hsp in alignment.hsps:
                            st.markdown("---")
                            st.success(f"**Gene Encontrado:** {alignment.title}")
                            st.info(f"**Identidade:** {hsp.identities}")
                            st.code(f"Sequência Correspondente: {hsp.sbjct}")
                else:
                    st.error("Nenhum gene encontrado no BLAST.")
            except Exception as e:
                st.error(str(e))
