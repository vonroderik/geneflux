import streamlit as st
from ..logic.education import (
    generate_random_dna, 
    get_mutation_impact, 
    generate_pcr_challenge, 
    generate_splicing_challenge,
    BIOCHEM_QUIZ
)
from ..logic.genetics import validate_genotype, calculate_punnett, calculate_abo_frequencies
from ..components.genetics_ui import render_punnett_square, render_abo_frequencies
from ..components.visualizations import render_codon_table

def render_education_page():
    st.subheader("Módulo Educativo")

    edu_option = st.sidebar.radio(
        "Escolha uma opção:",
        (
            "Herança Mendeliana",
            "Simulador de Mutação Gênica",
            "Desafio PCR: Primers",
            "Laboratório de Splicing",
            "Quiz de Bioquímica",
            "Tabela de Códons",
            "Glossário Bio Cel & Mol"
        ),
    )

    if edu_option == "Herança Mendeliana":
        st.subheader("Herança Mendeliana")
        st.write("Estude a segregação de alelos (1ª e 2ª Lei de Mendel).")
        
        cross_type = st.radio("Tipo de Cruzamento:", ["Monohíbrido / ABO", "Dihíbrido (2 genes)"])
        
        col1, col2 = st.columns(2)
        if cross_type == "Monohíbrido / ABO":
            with col1:
                g1 = st.text_input("Genótipo Progenitor 1:", "Aa", key="mon1").strip()
            with col2:
                g2 = st.text_input("Genótipo Progenitor 2:", "aa", key="mon2").strip()
            
            if st.button("Cruzamento Monohíbrido"):
                try:
                    validate_genotype(g1, "monohybrid")
                    validate_genotype(g2, "monohybrid")
                    results = calculate_punnett(g1, g2, "monohybrid")
                    render_punnett_square(results, "monohybrid")
                    
                    if any(c in "oO" for c in g1+g2) or "AB" in g1.upper() or "AB" in g2.upper():
                         st.markdown("---")
                         render_abo_frequencies(results["offspring"])
                except ValueError as e:
                    st.error(str(e))
        else:
            st.info("Formato: 4 letras, ex: AaBb")
            with col1:
                g1 = st.text_input("Genótipo Progenitor 1:", "AaBb", key="di1").strip()
            with col2:
                g2 = st.text_input("Genótipo Progenitor 2:", "AaBb", key="di2").strip()
                
            if st.button("Cruzamento Dihíbrido"):
                try:
                    validate_genotype(g1, "dihybrid")
                    validate_genotype(g2, "dihybrid")
                    results = calculate_punnett(g1, g2, "dihybrid")
                    render_punnett_square(results, "dihybrid")
                except ValueError as e:
                    st.error(str(e))

    elif edu_option == "Simulador de Mutação Gênica":
        st.subheader("Simulador de Mutação Gênica")
        if "mut_dna" not in st.session_state:
            st.session_state.mut_dna = generate_random_dna()
            
        if st.button("Gerar Nova Sequência"):
            st.session_state.mut_dna = generate_random_dna()
            
        st.info(f"Sequência Original: {st.session_state.mut_dna}")
        
        mut_type = st.selectbox("Tipo de Mutação:", ["Substituição", "Inserção", "Deleção"])
        pos = st.number_input("Posição (1-indexed):", 1, len(st.session_state.mut_dna), 1) - 1
        
        dna_list = list(st.session_state.mut_dna)
        if mut_type == "Substituição":
            base = st.selectbox("Nova Base:", ["A", "T", "C", "G"])
            dna_list[pos] = base
        elif mut_type == "Inserção":
            base = st.selectbox("Base Inserida:", ["A", "T", "C", "G"])
            dna_list.insert(pos, base)
        elif mut_type == "Deleção":
            dna_list.pop(pos)
            
        mutated_dna = "".join(dna_list)
        
        if st.button("Ver Impacto"):
            st.markdown(f"**Sequência Mutada:** `{mutated_dna}`")
            impact = get_mutation_impact(st.session_state.mut_dna, mutated_dna)
            st.write(f"**DNA:** {impact['dna_type']}")
            st.write(f"**Proteína:** {impact['protein_type']}")
            st.write(f"**Proteína Original:** {'-'.join(impact['original_protein'])}")
            st.write(f"**Proteína Mutada:** {'-'.join(impact['mutated_protein'])}")

    elif edu_option == "Desafio PCR: Primers":
        st.subheader("Desafio de Primers (PCR)")
        st.write("""O sucesso de uma PCR depende da escolha correta dos primers. 
        Lembre-se: o DNA polimeriza na direção **5' → 3'**.""")
        
        if "pcr_data" not in st.session_state or st.button("Novo Template de DNA"):
            st.session_state.pcr_data = generate_pcr_challenge()
        
        st.info(f"Sequência Template: **5' - {st.session_state.pcr_data['template']} - 3'**")
        
        fwd_user = st.text_input("Digite o Forward Primer (10nt):").upper().strip()
        rev_user = st.text_input("Digite o Reverse Primer (10nt):").upper().strip()
        
        if st.button("Validar Primers"):
            col1, col2 = st.columns(2)
            with col1:
                if fwd_user == st.session_state.pcr_data['fwd_correct']:
                    st.success("Forward Primer: Correto! (Pareia com o início da fita 3'->5')")
                else:
                    st.error(f"Forward incorreto. Esperado: {st.session_state.pcr_data['fwd_correct']}")
            
            with col2:
                if rev_user == st.session_state.pcr_data['rev_correct']:
                    st.success("Reverse Primer: Correto! (Lembrou de inverter e complementar!)")
                else:
                    st.error(f"Reverse incorreto. Dica: Deve ser o Complementar Reverso do final do template.")

    elif edu_option == "Laboratório de Splicing":
        st.subheader("Laboratório de Splicing")
        st.write("Identifique os Éxons (Maiúsculas) e Íntrons (Minúsculas) e processe o mRNA.")
        
        s_data = generate_splicing_challenge()
        st.markdown(f"**Pré-mRNA:** `{s_data['pre_mrna']}`")
        
        user_splice = st.text_input("Digite o mRNA Maduro (apenas os éxons):").strip()
        
        if st.button("Processar RNA"):
            if user_splice == s_data['mature_mrna']:
                st.success("Splicing perfeito! Os íntrons foram removidos e os éxons unidos.")
                st.info(f"O mRNA Maduro agora pode seguir para os Ribossomos para Tradução.")
            else:
                st.error("Erro no splicing. Você deixou resíduos de íntrons ou removeu partes dos éxons.")

    elif edu_option == "Quiz de Bioquímica":
        st.subheader("Quiz de Bioquímica: Ácidos Nucleicos")
        for i, q in enumerate(BIOCHEM_QUIZ):
            st.markdown(f"**Pergunta {i+1}:** {q['question']}")
            ans = st.radio(f"Escolha uma opção para Q{i+1}:", q['options'], key=f"q{i}")
            if st.button(f"Verificar Q{i+1}"):
                if ans == q['answer']:
                    st.success(f"Correto! {q['explanation']}")
                else:
                    st.error(f"Incorreto. A resposta certa era '{q['answer']}'.")

    elif edu_option == "Tabela de Códons":
        render_codon_table()
        
    elif edu_option == "Glossário Bio Cel & Mol":
        st.subheader("Glossário Rápido")
        glossary = {
            "Ácido Desoxirribonucleico (DNA)": "Molécula que carrega as instruções genéticas para o desenvolvimento e funcionamento dos seres vivos.",
            "Ácido Ribonucleico (RNA)": "Molécula responsável por traduzir a informação do DNA em proteínas.",
            "Código Genético": "Conjunto de regras que define como a sequência de nucleotídeos em um gene é traduzida em aminoácidos.",
            "Códon": "Sequência de três nucleotídeos que codifica um aminoácido específico.",
            "Éxon": "Segmento de um gene que contém a informação necessária para codificar uma proteína.",
            "Íntron": "Segmento de um gene que não codifica proteínas e é removido durante o splicing.",
            "ORF (Open Reading Frame)": "Sequência de DNA que começa com um códon de início e termina com um códon de parada, tendo o potencial de ser traduzida.",
            "Fenótipo": "As características observáveis de um organismo, resultantes da interação do seu genótipo com o ambiente.",
            "Genótipo": "A constituição genética de um indivíduo."
        }
        for term, desc in glossary.items():
            with st.expander(term):
                st.write(desc)
