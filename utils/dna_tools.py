from typing import List, Tuple, Literal, List
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.Seq import Seq
from Bio.SeqUtils import seq3
from rich import print as rprint
from . import interface as cli
import re

####### Validação de Entrada #######


# Usa uma Regex para verificar se a entrada do usuário é uma sequência de DNA ou RNA.
# Também verifica se a sequência não é uma mistura de DNA e RNA.
# Sequências sem 'U' são consideradas DNA.
# Retorna a sequência validada e uma string especificando o tipo de molécula.
def validate_sequence(seq) -> Tuple[str, Literal["dna", "rna"]]:

    if not seq:
        raise ValueError(
            "[bold red]\nEntrada vazia. Por favor, digite uma sequência de DNA ou RNA[/bold red]"
        )

    match = re.fullmatch(r"^([ATCGU]+)$", seq)

    if not match:
        raise ValueError(
            """[bold red]\n
    Sequência inválida.
    Por favor, digite uma sequência de ácido nucleico válida.
    Deve conter apenas os caracteres A, T, C, G ou U.
    [/bold red]"""
        )

    validated_seq = match.group(1)

    if len(validated_seq) < 3:
        raise ValueError(
            "[bold red]\nSua sequência deve ter pelo menos 3 nucleotídeos[/bold red]"
        )

    if "U" in validated_seq and "T" in validated_seq:
        raise ValueError(
            """[bold red]\n
    Sequência de DNA/RNA inválida: não é possível misturar bases de RNA e DNA.
    Sua sequência contém tanto timina (T) quanto uracila (U)
            [/bold red]"""
        )

    return validated_seq, "dna" if "U" not in validated_seq else "rna"


####### Funções de Ferramentas de DNA #######


# Retorna o filamento complementar de DNA (DNA -> DNA Complementar)
def complementary_dna(seq: str, show_results: bool) -> str:

    dna_seq = reverse_transcription(seq)
    comp_seq = dna_seq.translate(str.maketrans("ATCG", "TAGC"))

    if show_results:
        cli.show_complementary(comp_seq, dna_seq)
        cli.press_to_continue()
    return comp_seq


# Retorna um transcrito de RNA (DNA -> RNA)
def transcription(seq: str, show_results: bool) -> str:

    rna_seq = seq.replace("T", "U")

    if show_results:
        cli.show_transcription(rna_seq)
        cli.press_to_continue()
    return rna_seq


# Usa Biopython's Seq que converte a sequência em um objeto
# e retorna uma lista de aminoácidos (aa) de uma letra, e Biopython's Seq3
# que converte o aa de uma letra em um formato de três letras.
# Chama sequence_to_translate() para obter a sequência a ser traduzida: sequência inteira ou
# ORF específica.
def translation(seq: str, show_results: bool):

    rna_seq = transcription(seq, show_results=False)
    if len(rna_seq) < 3:
        rprint("[bold red]\nSua sequência deve ter pelo menos 3 códons. Por favor, digite uma nova sequência.[/bold red]")
        return
    selected_seq = sequence_to_translate(rna_seq)
    if selected_seq is None:
        return
    elif not selected_seq:
        rna_seq = rna_seq[: len(rna_seq) - len(rna_seq) % 3]
    else:
        rna_seq = rna_seq[: len(selected_seq) - len(selected_seq) % 3]
    protein = str(Seq(rna_seq).translate(to_stop=False))
    aa_list = [seq3(aa).title() for aa in protein]

    if show_results:
        cli.show_translation(aa_list)
        cli.press_to_continue()
    return aa_list


# Pede ao usuário para escolher entre traduzir toda a sequência de DNA/RNA,
# nesse caso retorna a sequência não modificada, ou escolher de uma ORF,
# nesse caso chama orf_finder(), e retorna uma lista de ORFs, se houver,
# para o usuário escolher.
def sequence_to_translate(seq: str):

    while True:
        try:
            option = int(
                input(
                    "[bold blue]Por favor, selecione uma opção:[/bold blue]\n"
                    "[bold blue]1.[/bold blue] Traduzir a sequência de RNA completa, ignorando códons de início e parada\n"
                    "[bold blue]2.[/bold blue] Selecionar uma ORF para traduzir\n"
                    "[bold blue]3.[/bold blue] Voltar ao menu principal\n"
                    "[yellow]Opção: [/yellow]"
                )
            )

            if option == 1:
                return False
            elif option == 2:
                orfs_choice = orf_finder(seq, show_results=True, wait=False)
                if orfs_choice:
                    while True:
                        try:
                            chosen_orf = int(
                                input(
                                    "[yellow]\nPor favor, selecione uma ORF, por número, para traduzir: [/yellow]"
                                )
                            )
                        except ValueError:
                            rprint("[bold red]\nPor favor, digite um número válido[/bold red]")
                            continue
                        return orfs_choice[chosen_orf - 1]
            elif option == 3:
                return None
            else:
                rprint("[bold red]\nPor favor, escolha uma opção válida[/bold red]")
        except ValueError:
            rprint("[bold red]\nEntrada inválida. Por favor, digite um número.[/bold red]")


# Usa loops while e for para procurar ORFs.
# Procura nos três quadros de RNA (+1, +2, +3) por um códon de início (AUG).
# Qualquer coisa entre um AUG e um códon de parada (UGA, UAG, UAA) é considerada uma ORF.
# A última ORF não precisa ter um códon de parada.
def orf_finder(seq: str, show_results: bool, wait: bool) -> List[str] | None:

    rna_seq = transcription(seq, show_results=False)
    if len(rna_seq) < 3:
        rprint("[bold red]\nSua sequência deve ter pelo menos 3 nucleotídeos. Por favor, digite uma nova sequência.[/bold red]")
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
    if show_results:
        if not orfs:
            rprint("[bold red]\nORF não encontrada[/bold red]")
        else:
            cli.show_orfs(orfs)
        if wait:
            cli.press_to_continue()

    return orfs


# Usa o NCBIWWW do Biopython para enviar uma consulta de pesquisa ao NCBI para BLAST,
# e recebe um arquivo XML que é convertido em um objeto de Classe por NCBIXML
def gene_identifier(seq: str, show_results: bool):

    dna_seq = reverse_transcription(seq)

    result_handle = None

    if len(seq) < 11:
        rprint(
            "[bold red]\nSua sequência deve ter pelo menos 11 nucleotídeos. Por favor, digite uma nova sequência[/bold red]"
        )
        return
    rprint("[cyan]\nEnviando sequência para o NCBI BLAST. Isso pode levar até um minuto, por favor, aguarde...[/cyan]")

    try:
        result_handle = NCBIWWW.qblast(
            program="blastn",
            database="nt",
            sequence=dna_seq,
        )

        blast_record = NCBIXML.read(result_handle)

        if not blast_record.alignments:
            rprint("[bold red]\nNenhum gene associado à sequência de DNA foi encontrado[/bold red]")
            return None

        if show_results:
            cli.show_gene_identifier(blast_record)
            cli.press_to_continue()
        return blast_record

    except Exception as e:
        rprint(f"[bold red]\nA solicitação BLAST falhou: {e}[/bold red]")
        rprint("[bold red]Verifique sua conexão ou tente novamente mais tarde.[/bold red]")
        return

    finally:
        if result_handle is not None:
            result_handle.close()


# Retorna uma sequência de DNA a partir de uma sequência de RNA (RNA -> DNA)
def reverse_transcription(seq: str) -> str:

    return seq if "U" not in seq else seq.replace("U", "T")

