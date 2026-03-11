import streamlit as st
from ..logic.clinical_cases import CLINICAL_CASES, compare_clinical_sequences
from ..logic.dna_rna import transcription, translation

def render_clinical_page():
    st.subheader("Casos Clínicos: Diagnóstico Molecular")
    st.markdown("Analise casos reais e descubra a origem molecular das patologias.")

    case_name = st.selectbox("Selecione um caso clínico:", list(CLINICAL_CASES.keys()))
    case_data = CLINICAL_CASES[case_name]

    st.info(f"**Histórico Clínico:** {case_data['description']}")
    
    col1, col2 = st.columns(2)
    with col1:
        st.markdown("### Sequência de Referência (Normal)")
        st.code(case_data['reference_dna'])
    with col2:
        st.markdown("### Sequência do Paciente")
        st.code(case_data['patient_dna'])

    st.markdown("---")
    st.subheader("🔬 Laboratório de Diagnóstico")
    
    st.write("Sua tarefa: Identifique a posição e a mudança de base no DNA do paciente.")
    
    col_pos, col_base = st.columns(2)
    with col_pos:
        user_pos = st.number_input("Posição da Mutação (1, 2, 3...):", 1, len(case_data['reference_dna']), 1)
    with col_base:
        user_base = st.selectbox("Qual base você encontrou no Paciente?", ["A", "T", "C", "G"])

    if st.button("Validar Diagnóstico"):
        diffs = compare_clinical_sequences(case_data['reference_dna'], case_data['patient_dna'])
        found_match = False
        for d in diffs:
            if d['pos'] == user_pos and d['pat'] == user_base:
                found_match = True
                break
        
        if found_match:
            st.success(f"**Correto!** Você identificou a mutação na posição {user_pos}.")
            
            # Show molecular impact
            st.markdown("---")
            st.subheader("Análise Molecular Detalhada")
            
            ref_rna = transcription(case_data['reference_dna'])
            pat_rna = transcription(case_data['patient_dna'])
            
            ref_prot = translation(ref_rna)
            pat_prot = translation(pat_rna)
            
            # Find which amino acid changed
            codon_idx = (user_pos - 1) // 3
            ref_aa = ref_prot[codon_idx]
            pat_aa = pat_prot[codon_idx]
            
            col_a, col_b = st.columns(2)
            with col_a:
                st.write("**Proteína Referência:**")
                st.write(f"`... {' - '.join(ref_prot[max(0, codon_idx-2):codon_idx+3])} ...`")
                st.markdown(f"**Aminoácido {codon_idx+1}:** {ref_aa}")
            with col_b:
                st.write("**Proteína Paciente:**")
                st.write(f"`... {' - '.join(pat_prot[max(0, codon_idx-2):codon_idx+3])} ...`")
                st.markdown(f"**Aminoácido {codon_idx+1}:** {pat_aa}")
            
            st.warning(f"**Conclusão Clínica:** {case_data['clinical_context']}")
            
        else:
            st.error("Diagnóstico incorreto. Compare as sequências cuidadosamente base por base.")
            if st.checkbox("Mostrar Dica"):
                st.write(f"Dica: Olhe para a região próxima à posição {case_data['mutation_pos']}.")
