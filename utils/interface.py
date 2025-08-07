import pyfiglet
from rich.console import Console
from rich.panel import Panel
from rich.text import Text
from rich import print as rprint
from rich.prompt import Prompt

console = Console()

LINE_SIZE = 12


####### Funções CLI #######


# Exibe o DNA complementar (complementary()) de forma formatada
def show_complementary(comp_dna_seq: str, dna_seq: str) -> None:

    rprint(Panel("[cyan]SEQUÊNCIA DE DNA COMPLEMENTAR[/cyan]"))

    for i in range(0, len(comp_dna_seq), LINE_SIZE):
        dna_seq_fragment = " ".join(dna_seq[i : i + LINE_SIZE])
        comp_dna_seq_fragment = " ".join(comp_dna_seq[i : i + LINE_SIZE])

        rprint(f"[green]{i:04}  5'  {dna_seq_fragment}  3'[/green]")
        rprint(f"[blue]{i:04}  3'  {comp_dna_seq_fragment}  5'[/blue]\n")


# Exibe o transcrito de RNA (transcription()) de forma formatada
def show_transcription(rna: str) -> None:

    rprint(Panel("[cyan]SEQUÊNCIA DE RNA[/cyan]"))

    for i in range(0, len(rna), LINE_SIZE):
        rna_seq = " ".join(rna[i : i + LINE_SIZE])

        rprint(
            f"[bold green]{i:04} - "
            f"[/bold green][yellow]5'  {rna_seq}  3'[/yellow]"
        )
        rprint()


# Exibe a sequência de aminoácidos (translation()) de forma formatada
def show_translation(protein: list) -> None:

    rprint(Panel("[cyan]SEQUÊNCIA DE AMINOÁCIDOS[/cyan]"))

    colored_aa = []
    for aa in protein:
        if aa == "Met":
            colored_aa.append(f"[yellow]{aa}[/yellow]")
        else:
            colored_aa.append(f"[magenta]{aa}[/magenta]")

    for i in range(0, len(colored_aa), LINE_SIZE):
        aa_seq = " - ".join(colored_aa[i : i + LINE_SIZE])

        rprint()
        rprint(
            f"[bold green]{i:04} - [/bold green]{aa_seq}"
        )


# Exibe as ORFs identificadas (orf_finder()) de forma formatada
def show_orfs(orfs: list) -> None:

    rprint(Panel("[cyan]OPEN READING FRAMES[/cyan]"))

    for index, orf in enumerate(orfs, 1):
        rprint(
            f"\n[bold green]ORF - {index:02} [/bold green] - "
            f"[magenta]5'  {orf}  3'[/magenta]"
        )


# Exibe o gene identificado (gene_identifier()) de forma formatada
def show_gene_identifier(gene) -> None:

    rprint(Panel("[cyan]IDENTIFICADOR DE GENE[/cyan]"))

    blast_record = gene

    for alignment in blast_record.alignments[:3]:
        for hsp in alignment.hsps:
            rprint("[green]\nGene Possível Encontrado:[/green]")
            rprint(f"[magenta]Título:[/magenta] {alignment.title}")
            rprint(f"[yellow]Identidade:[/yellow] {hsp.identities}")
            rprint(f"[green]Sequência Correspondente:[/green] {hsp.sbjct}")
            rprint("[bold cyan]\n--------[/bold cyan]")


# Pede ao usuário para pressionar ENTER antes
# que menu() seja chamado novamente por outras funções
def press_to_continue() -> str:
    return input("\n[red]Pressione ENTER para continuar...[/red]\n")


####### MENU #######


# Exibe o nome do programa na primeira vez que é inicializado
def welcome() -> None:
    # `pyfiglet.figlet_format` retorna uma string, então usamos Text.from_ansi
    # para renderizar as cores e estilos do Rich.
    figlet_text = pyfiglet.figlet_format("GeneFlux", font="slant")
    rich_text = Text.from_ansi(f"[bold magenta]{figlet_text}[/bold magenta]")
    rprint(rich_text, end="")
    rprint(
        Text.from_markup(
            """
          [cyan]GeneFlux v1.0
    DNA → RNA → ORFs → Protein
        [/cyan]"""
        )
    )


# Responsável por mostrar o menu principal e solicitar a escolha do usuário
def show_menu() -> int:

    rprint(Panel("[bold blue]=== GeneFlux: DNA Sequence Analyzer ===", style="blue"))

    rprint(
        Text.from_markup(
            """
[bold blue]Escolha uma operação para realizar:[/bold blue]

[bold blue]1.[/bold blue] Gerar filamento de DNA complementar       (DNA → DNA)
[bold blue]2.[/bold blue] Transcrever para RNA                       (DNA → RNA)
[bold blue]3.[/bold blue] Traduzir sequência para aminoácidos       (RNA → Proteína)
[bold blue]4.[/bold blue] Localizador de ORF                              (Códon de Início → Códon de Parada)
[bold blue]5.[/bold blue] Identificador de Gene                         (DNA → Gene)
[bold blue]6.[/bold blue] Digitar uma nova sequência de DNA/RNA

[bold blue]7.[/bold blue] Voltar
"""
        )
    )

    # O input não é formatado pelo Rich, então permanece como estava
    return int(input("[yellow]Digite sua opção: [/yellow]").strip())


def show_edu_menu() -> int:

    rprint(Panel("[bold blue]=== GeneFlux: DNA Sequence Analyzer ===", style="blue"))
    rprint(Panel("[bold blue]      === Módulo Educacional ===", style="blue"))
    rprint(
        Text.from_markup(
            """
[bold blue]Escolha uma opção:[/bold blue]

[bold blue]1.[/bold blue] Exercícios com DNA Complementar
[bold blue]2.[/bold blue] Exercícios de Transcrição
[bold blue]3.[/bold blue] Exercícios de Tradução
[bold blue]4.[/bold blue] Encontre uma ORF

[bold blue]0.[/bold blue] Voltar
"""
        )
    )
    
    # O input não é formatado pelo Rich, então permanece como estava
    return int(input("[yellow]Digite sua opção: [/yellow]").strip())
