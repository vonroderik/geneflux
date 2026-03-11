import streamlit as st
from typing import List
from ..logic.education import CODON_TABLE


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
    output = ""
    for i in range(0, len(protein_list), line_size):
        aa_fragment = " - ".join(protein_list[i : i + line_size])
        output += f"{i:04} - {aa_fragment}\n"
    st.code(output)


def render_codon_table():
    """Renderiza a tabela de códons com cores adaptáveis e categorização bioquímica."""
    st.subheader("Tabela de Códons (DNA)")

    # Mapeamento de propriedades para cores (tons suaves que funcionam em ambos os temas)
    # Cores: Apolares (Amarelo), Polares (Verde), Ácidos (Vermelho), Básicos (Azul), Parada (Cinza/Preto)
    aa_colors = {
        "Phe": "#FFF9C4",
        "Leu": "#FFF9C4",
        "Ile": "#FFF9C4",
        "Met": "#FFF9C4",
        "Val": "#FFF9C4",
        "Ala": "#FFF9C4",
        "Trp": "#FFF9C4",
        "Pro": "#FFF9C4",
        "Gly": "#FFF9C4",
        "Ser": "#C8E6C9",
        "Thr": "#C8E6C9",
        "Tyr": "#C8E6C9",
        "Cys": "#C8E6C9",
        "Asn": "#C8E6C9",
        "Gln": "#C8E6C9",
        "Asp": "#FFCDD2",
        "Glu": "#FFCDD2",
        "Lys": "#BBDEFB",
        "Arg": "#BBDEFB",
        "His": "#BBDEFB",
        "STOP": "#B83737",
    }

    bases = ["T", "C", "A", "G"]

    html = """
    <style>
        .codon-table { 
            width: 100%; border-collapse: collapse; text-align: center; 
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            color: black !important; /* Forçar texto escuro dentro das células coloridas para contraste */
        }
        .codon-table th, .codon-table td { border: 1px solid #ddd; padding: 4px; }
        .header-base { background-color: #f0f2f6; font-weight: bold; color: #31333F; }
        .side-base { background-color: #f0f2f6; font-weight: bold; width: 40px; color: #31333F; }
        .codon-cell { font-size: 0.85em; height: 80px; }
        .codon-item { margin: 2px 0; padding: 2px; border-radius: 3px; }
        .stop-aa { color: white !important; font-weight: bold; }
    </style>
    """

    html += "<table class='codon-table'>"
    html += "<tr><th rowspan='2' class='header-base'>1ª</th><th colspan='4' class='header-base'>2ª Base</th><th rowspan='2' class='header-base'>3ª</th></tr>"
    html += "<tr>"
    for b in bases:
        html += f"<th class='header-base'>{b}</th>"
    html += "</tr>"

    for b1 in bases:
        html += "<tr>"
        html += f"<td class='side-base'>{b1}</td>"
        for b2 in bases:
            html += "<td class='codon-cell'>"
            for b3 in bases:
                codon = b1 + b2 + b3
                aa = CODON_TABLE[b1][b2][b3]
                bg_color = aa_colors.get(aa, "#ffffff")
                text_class = "stop-aa" if aa == "STOP" else ""
                html += (
                    f"<div class='codon-item' style='background-color: {bg_color};'>"
                )
                html += f"<span style='font-weight: normal;'>{codon}:</span> <span class='{text_class}' style='font-weight: bold;'>{aa}</span>"
                html += "</div>"
            html += "</td>"

        # 3rd base column
        html += "<td class='side-base'>"
        for b3 in bases:
            html += f"<div style='margin: 10px 0;'>{b3}</div>"
        html += "</td>"
        html += "</tr>"

    html += "</table>"

    # Legenda Pedagógica
    st.markdown(html, unsafe_allow_html=True)
    st.markdown(
        """
        <div style='display: flex; gap: 10px; font-size: 0.8em; margin-top: 10px; justify-content: center;'>
            <span style='background-color: #FFF9C4; padding: 2px 5px; border-radius: 3px; color: black;'>Apolar</span>
            <span style='background-color: #C8E6C9; padding: 2px 5px; border-radius: 3px; color: black;'>Polar</span>
            <span style='background-color: #FFCDD2; padding: 2px 5px; border-radius: 3px; color: black;'>Ácido</span>
            <span style='background-color: #BBDEFB; padding: 2px 5px; border-radius: 3px; color: black;'>Básico</span>
            <span style='background-color: #212121; padding: 2px 5px; border-radius: 3px; color: white;'>Parada</span>
        </div>
    """,
        unsafe_allow_html=True,
    )
