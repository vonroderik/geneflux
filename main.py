import sys
from rich.console import Console
from rich.text import Text
from rich.panel import Panel
from rich import print as rprint
from utils import interface as cli
from utils import menus as menu

# Cria uma instância do console para usar Rich
console = Console()

def main() -> None:

    while True:
        try:
            cli.welcome()
            rprint(Panel("[bold blue]=== GeneFlux: DNA Sequence Analyzer ===", style="blue"))
            
            # Usando rprint para exibir o menu com estilo Rich
            rprint(
                Text.from_markup(
                    """
        [bold blue]Escolha uma opção:[/bold blue]

        [bold blue]1.[/bold blue] Módulo de Análises
        [bold blue]2.[/bold blue] Módulo Educativo
        [bold blue]3.[/bold blue] Sair
        """
                )
            )

            # O input não é formatado pelo Rich, então permanece como estava
            module_choice = int(
                input("[yellow]Digite sua opção: [/yellow]").strip()
            )

            if module_choice == 1:
                menu.tools_menu()
            elif module_choice == 2:
                menu.edu_menu()
            elif module_choice == 3:
                sys.exit("[green]Encerrando programa[/green]")
            else:
                rprint("[bold red]Por favor selecione uma opção válida (1-3).[/bold red]")
        except ValueError:
            rprint("[bold red]Opção inválida. Por favor digite um número[/bold red]")
            continue


if __name__ == "__main__":

    try:
        main()
    except (KeyboardInterrupt, EOFError):
        rprint("[red]\nEncerrando programa.[/red]")
    finally:
        rprint("[/]")
