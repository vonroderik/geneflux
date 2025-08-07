from rich import print as rprint
from . import dna_tools as gene
from . import interface as cli
from edu import educational as edu


# Responsável por controlar o Menu relacionado
# ao módulo de Ferramentas de DNA.
def tools_menu() -> None:

    while True:
        rprint(
            """
            [yellow]"Digite uma sequência de DNA ou RNA,
            ou digite [red]'sair'[/red][yellow] para voltar:[/yellow]"""
        )

        user_sequence = (
            input("[bold red]\n>> [/bold red]")
            .strip()
            .upper()
            .replace(" ", "")
        )

        if user_sequence == "SAIR":
            break

        try:
            seq, mol_type = gene.validate_sequence(user_sequence)

            rprint(
                f"[cyan]\n{mol_type.upper()} sequence detected\n[/cyan]"
            )

        except ValueError as e:
            rprint(f"[bold red] \n\tError: {e} [/bold red]")
            continue

        while True:
            try:
                choice = cli.show_menu()
                if choice == 1:
                    gene.complementary_dna(seq, show_results=True)
                elif choice == 2:
                    gene.transcription(seq, show_results=True)
                elif choice == 3:
                    gene.translation(seq, show_results=True)
                elif choice == 4:
                    gene.orf_finder(seq, show_results=True, wait=True)
                elif choice == 5:
                    gene.gene_identifier(seq, show_results=True)
                elif choice == 6:
                    break
                elif choice == 7:
                    return None
                else:
                    rprint("[bold red]Por favor, selecione uma opção válida do menu (1-7).[/bold red]")

            except ValueError:
                rprint("[bold red]Entrada inválida. Por favor, digite um número[/bold red]")
                continue

# Responsável por controlar o Menu relacionado
# ao módulo Educacional.
def edu_menu() -> None:
    
    while True:
        rprint(
            """
            [yellow]"Escolha uma opção,
            ou digite [red]'sair'[/red][yellow] para voltar:[/yellow]"""
        )
        
        while True:
            try:
                edu_menu_choice = cli.show_edu_menu()
                if edu_menu_choice == 1:
                    edu.edu_complementary_dna()
                elif edu_menu_choice == 2:
                    edu.edu_transcription()
                elif edu_menu_choice == 3:
                    edu.edu_translation()
                elif edu_menu_choice == 4:
                    edu.edu_orf()
                elif edu_menu_choice == 0:
                    return None
                else:
                    rprint("[bold red]Por favor, selecione uma opção válida do menu (0-4).[/bold red]")

            except ValueError:
                rprint("[bold red]Entrada inválida. Por favor, digite um número[/bold red]")
                continue
