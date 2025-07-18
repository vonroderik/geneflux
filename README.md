
-----

# GeneFlux

GeneFlux é um **aplicativo web interativo**, construído com **Streamlit**, que oferece ferramentas básicas de biologia molecular, como Transcrição de DNA, Tradução de RNA, localizador de ORF (Open Reading Frame) e localizador de genes, usando NCBI BLASTn.

-----

## Estrutura do Projeto:

```
GeneFlux/
├── main.py                   # Ponto de entrada do aplicativo Streamlit
├── requirements.txt          # Lista as dependências do projeto
├── README.md                 # Documentação principal
├── utils/                    # Módulos auxiliares e lógica de negócio
│   ├── __init__.py           # Inicializa o pacote utils
│   ├── dna_tools.py          # Funções para ferramentas de biologia molecular
│   └── interface.py          # Funções para exibição de resultados na interface web
└── tests/                    # Testes automatizados do projeto
    ├── __init__.py           # Inicializa o pacote tests
    └── test_dna_tools.py     # Testes unitários para as funções de dna_tools.py
```

-----

# Bibliotecas Importadas

As seguintes bibliotecas foram utilizadas neste projeto:

  * **Streamlit**: Usada para construir a interface web interativa do aplicativo.

  * **Biopython**: Três módulos foram utilizados:

      * `Seq` e `seq3` para traduzir RNA em sequências de aminoácidos.
      * `NCBIWWW` para enviar sequências ao servidor NCBI BLAST.
      * `NCBIXML` para analisar os resultados XML retornados.

  * **pytest**: Executa testes de unidade (**test\_dna\_tools.py**).

Duas bibliotecas internas também foram usadas: `typing` para dicas de tipo e `re` (regex) para validar a sequência de DNA/RNA inserida pelo usuário.

Todas as bibliotecas necessárias estão listadas no arquivo **requirements.txt**. Para instalá-las, basta digitar no seu terminal:
`pip install -r requirements.txt`

-----

# Funções do Código

Este programa contém um **test\_dna\_tools.py** que usa `pytest` para executar testes de unidade. Apenas as funções *DNA Tools* e `validate_sequence()` foram testadas. A função `gene_identifier()` não foi testada, pois a resposta dos servidores do NCBI pode levar até dois minutos.

Várias funções foram criadas para este programa. Elas foram divididas em categorias e estão organizadas nos módulos **dna\_tools.py** e **interface.py** dentro da pasta `utils/`, além de **main.py**.

-----

## Funções de Interface (Streamlit)

Estas funções são responsáveis por exibir o título do aplicativo e formatar os resultados das análises de biologia molecular na interface web do Streamlit. Elas estão localizadas em **utils/interface.py**.

**welcome\_streamlit()**: Exibe o título principal do aplicativo "GeneFlux" e informações da versão.

**show\_complementary\_streamlit()**: Exibe a sequência de DNA complementar de forma formatada.

**show\_transcription\_streamlit()**: Exibe o transcrito de RNA de forma formatada.

**show\_translation\_streamlit()**: Exibe a sequência de aminoácidos de forma formatada, destacando o códon de início (Metionina).

**show\_orfs\_streamlit()**: Exibe as Open Reading Frames (ORFs) encontradas.

**show\_gene\_identifier\_streamlit()**: Exibe os resultados do BLASTn para identificação de genes.

-----

## Funções de Ferramentas de DNA

O coração do programa. Responsáveis pelas ferramentas de biologia molecular, localizadas em **utils/dna\_tools.py**.

Todas as funções de ferramentas de DNA (com sufixo `_streamlit`) interagem com o `st.session_state` para armazenar os resultados e sinalizar à interface principal para exibi-los. As funções auxiliares internas (sem sufixo `_streamlit`) não interagem diretamente com o Streamlit.

**complementary\_dna\_streamlit()**: Recebe uma sequência de DNA e retorna a versão complementar da mesma sequência: *ATCG -\> TAGC*.

**transcription\_streamlit()**: Converte um DNA em um RNA substituindo (*T*)imina por (*U*)racila. *CGTA -\> CGUA*.

**translation\_streamlit()**: Verifica se a sequência de DNA/RNA tem pelo menos 3 códons. Se Verdadeiro, a sequência é convertida para RNA. Utiliza a biblioteca *Biopython* para traduzir a sequência em aminoácidos. Permite ao usuário escolher entre traduzir a sequência completa ou uma ORF específica.

**orf\_finder\_streamlit()**: Procura por Open Reading Frames (ORFs) em uma sequência de RNA. Para este projeto, as ORFs são definidas como qualquer sequência de RNA entre um *códon de início* (AUG - *incluído*) e um *códon de parada* (UAA, UAG, UGA - *excluído*), lidas em três quadros (+1, +2, +3). A função retorna uma lista de ORFs encontradas.

**gene\_identifier\_streamlit()**: Faz uso do *Biopython* para submeter uma sequência de DNA ao NCBI BLASTn. A sequência deve ter pelo menos 11 nucleotídeos. Os resultados são analisados e, se genes correspondentes forem encontrados, são exibidos. Esta função depende da resposta do servidor NCBI BLAST e pode levar alguns minutos.

**reverse\_transcription()**: Uma função auxiliar interna que converte uma sequência de RNA de volta para DNA substituindo (*U*)racila por (*T*)imina. *CGUA -\> CGTA*.

-----

## Outras Funções

**main()** é a função principal do aplicativo, localizada em **main.py**. Ela orquestra a interface do usuário no Streamlit, gerencia a entrada da sequência, valida-a e, com base na escolha do usuário, chama as funções apropriadas de ferramentas de DNA. Também é responsável por exibir os resultados utilizando as funções de interface.

**validate\_sequence()**: Esta função usa *expressões regulares* para verificar se uma sequência é válida. Para isso, uma sequência válida deve atender a todos os seguintes critérios:

  * Contém apenas caracteres válidos: *A, C, G, T, U*.
  * Ter pelo menos 3 nucleotídeos de comprimento.
  * Conter **ou** uma (*U*)racila ou uma (*T*)imina, mas **não ambas**.

-----

# COMO USAR E EXEMPLO DE SAÍDA:

Para executar o aplicativo GeneFlux, certifique-se de ter o Streamlit e as outras dependências instaladas (via `pip install -r requirements.txt`). Em seguida, navegue até o diretório raiz do projeto no seu terminal e execute:

```bash
streamlit run main.py
```

Isso abrirá o aplicativo no seu navegador padrão.

  * **Inserir Sequência**: Na barra lateral, digite sua sequência de DNA ou RNA.
  * **Validar e Detectar Tipo**: O aplicativo validará a sequência e indicará se é DNA ou RNA.
  * **Escolher Operação**: Selecione uma das operações disponíveis no menu lateral (Complementar DNA, Transcrever, Traduzir, Encontrar ORF, Identificar Gene).
  * **Visualizar Resultados**: Os resultados da operação escolhida serão exibidos na área principal do aplicativo.

-----

## Exemplo de Análise:

Com a sequência `ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAG`:

**Entrada da Sequência**: `ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAG`

**Detecção**: "Sequência de DNA detectada\!"

**Gerar Fita de DNA Complementar**:

```
SEQUÊNCIA DE DNA COMPLEMENTAR
0000 5' A T G G T G C A C C T G A C T C C T G A G G A G A A G T C T G C C G T T A C T G C C C T G T G G G G C A A G 3'
0000 3' T A C C A C G T G G A C T G A G G A C T C C T T C A G A C G G C A A T G A C G G G A C A C C C C G T T C 5'
```

**Transcrever**:

```
SEQUÊNCIA DE RNA
0000 - 5' A U G G U G C A C C U G A C U C C U G A G G A G A A G U C U G C C G U U A C U G C C C U G U G G G G C A A G 3'
```

**Traduzir**:

```
SEQUÊNCIA DE AMINOÁCIDOS
0000 - Met - Val - His - Leu - Thr - Pro - Glu - Glu - Lys - Ser - Ala - Val - Thr - Ala - Leu - Trp - Gly - Lys
```

**Encontrar ORF**:

```
OPEN READING FRAMES
Nenhum ORF encontrado. (Isso acontece porque há um códon de início, mas nenhum códon de parada na sequência fornecida, dentro dos critérios definidos para este projeto.)
```

**Identificar Gene**:

```
IDENTIFICADOR DE GENE
Possível gene encontrado:
Título: gi|1584133638|gb|MK476301.1| Homo sapiens voucher Nzime_44_0 hemoglobin subunit beta (HBB) gene, complete cds
Identidade: 54
Sequência Correspondente: ATGGTGCACCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAG
```

-----

# SOBRE O AUTOR

**Rodrigo Mello** é professor universitário com graduação em **Ciências Biomédicas** e mestrado em **Genética e Toxicologia**. Atualmente, ele está quase terminando seu doutorado em **Biologia Celular e Molecular**.

Recentemente, ele descobriu uma paixão por **Python** e programação, e tem aprendido a codificar para enriquecer sua carreira e as Ciências Biomédicas.