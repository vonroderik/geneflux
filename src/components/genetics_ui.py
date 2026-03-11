import streamlit as st
from typing import Dict, List
from ..logic.genetics import get_abo_phenotype

def render_punnett_square(punnett_results: Dict, mode: str = "monohybrid"):
    """Renderiza o quadro de Punnett visualmente."""
    gametes1 = punnett_results["gametes1"]
    gametes2 = punnett_results["gametes2"]
    offspring = punnett_results["offspring"]
    genotype_freq = punnett_results["genotype_freq"]

    st.subheader(f"Quadro de Punnett ({'Mono' if mode == 'monohybrid' else 'Di'}híbrido)")
    
    # CSS for the table
    st.markdown("""
        <style>
        .punnett-table { width: 100%; border-collapse: collapse; margin-top: 10px; }
        .punnett-table th, .punnett-table td { border: 1px solid #4B4B4B; padding: 10px; text-align: center; }
        .punnett-table th { background-color: #262730; font-weight: bold; }
        .punnett-gen { font-family: monospace; font-size: 1.1em; font-weight: bold; color: #ff4b4b; }
        </style>
    """, unsafe_allow_html=True)

    table_html = "<table class='punnett-table'>"
    table_html += "<tr><th>♂ / ♀</th>"
    for g2 in gametes2:
        table_html += f"<th>{g2}</th>"
    table_html += "</tr>"

    idx = 0
    for g1 in gametes1:
        table_html += f"<tr><th>{g1}</th>"
        for _ in gametes2:
            child = offspring[idx]
            table_html += f"<td><span class='punnett-gen'>{child}</span></td>"
            idx += 1
        table_html += "</tr>"
    table_html += "</table>"
    
    st.markdown(table_html, unsafe_allow_html=True)
    
    st.markdown("---")
    
    col1, col2 = st.columns(2)
    with col1:
        st.write("### Frequência Genotípica")
        total = len(offspring)
        for gen, count in genotype_freq.items():
            perc = (count / total) * 100
            st.write(f"- **{gen}**: {count}/{total} ({perc:.1f}%)")
            
    with col2:
        st.write("### Frequência Fenotípica")
        # Mendelian phenotype (simplified: A is dominant)
        if mode == "monohybrid":
             phenotypes = []
             for gen in offspring:
                 if any(c.isupper() for c in gen):
                     phenotypes.append("Dominante")
                 else:
                     phenotypes.append("Recessivo")
             
             from collections import Counter
             p_freq = Counter(phenotypes)
             for phen, count in p_freq.items():
                 perc = (count / total) * 100
                 st.write(f"- **{phen}**: {count}/{total} ({perc:.1f}%)")
        else:
            st.info("Frequência fenotípica complexa em dihibridismo.")

def render_abo_frequencies(offspring: List[str]):
    """Renderiza frequências ABO."""
    st.write("### Frequência Fenotípica (Sistema ABO)")
    from collections import Counter
    phenotypes = [get_abo_phenotype(gen) for gen in offspring]
    p_freq = Counter(phenotypes)
    total = len(offspring)
    
    for phen, count in p_freq.items():
        perc = (count / total) * 100
        st.write(f"- **Tipo {phen}**: {count}/{total} ({perc:.1f}%)")
