import streamlit as st
from src.pages.analysis_page import render_analysis_page
from src.pages.education_page import render_education_page
from src.pages.clinical_page import render_clinical_page

# Configuração da Página
st.set_page_config(
    page_title="GeneFlux: Análise de Sequências e Genética", 
    layout="wide",
    initial_sidebar_state="expanded"
)

def main():
    st.sidebar.title("GeneFlux 🧬")
    st.sidebar.markdown("""
        **Auxiliar de Biologia Celular e Molecular**
        ---
    """)

    module_choice = st.sidebar.radio(
        "Escolha um módulo:", 
        ("Módulo de Análises", "Módulo Educativo", "Casos Clínicos")
    )

    st.title("GeneFlux")
    st.markdown("""
        #### DNA → RNA → Proteína | Genética Mendeliana
        *Ferramenta didática para estudantes de Biomedicina*
    """)
    st.markdown("---")

    if module_choice == "Módulo de Análises":
        render_analysis_page()
    elif module_choice == "Módulo Educativo":
        render_education_page()
    else:
        render_clinical_page()


    # Footer
    st.sidebar.markdown("---")
    st.sidebar.info("Desenvolvido por Rodrigo Mello para fins educacionais.")


if __name__ == "__main__":
    main()
