from Bio.Blast import NCBIWWW, NCBIXML
import streamlit as st
from .dna_rna import reverse_transcription

def gene_identifier(seq: str):
    """Envia uma sequência para o NCBI BLAST para identificação de genes."""
    dna_seq = reverse_transcription(seq)
    if len(seq) < 11:
        raise ValueError("Sua sequência deve ter pelo menos 11 nucleotídeos para o BLAST.")

    try:
        result_handle = NCBIWWW.qblast(
            program="blastn",
            database="nt",
            sequence=dna_seq,
        )
        blast_record = NCBIXML.read(result_handle)
        result_handle.close()
        return blast_record
    except Exception as e:
        raise RuntimeError(f"A solicitação BLAST falhou: {e}")
