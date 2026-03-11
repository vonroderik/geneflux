# GeneFlux 🧬

GeneFlux é um ecossistema educacional interativo para estudantes de **Biomedicina**, com foco em **Biologia Celular e Molecular** e **Genética**. Desenvolvido com **Streamlit** e **Biopython**, o projeto visa transformar conceitos complexos em ferramentas práticas e visuais.

---

## ✨ Novidades da Versão 2.0

- **Arquitetura Modular:** Código refatorado seguindo padrões de separação de responsabilidades (Logic, Components, Pages).
- **Módulo Educativo Expandido:**
    - **Herança Mendeliana:** Suporte a cruzamentos Monohíbridos (incluindo Sistema ABO) e Dihíbridos (9:3:3:1).
    - **Simulador de Mutação:** Visualize o impacto de substituições, inserções e deleções no DNA, RNA e na proteína resultante.
    - **Tabela de Códons Interativa:** Uma referência visual dinâmica para tradução.
    - **Glossário Bio Cel & Mol:** Termos fundamentais sempre à mão.
- **Métricas de Bioinformática:** Cálculo automático de **Conteúdo GC**, **Peso Molecular** e comprimento da sequência.
- **Integração NCBI BLAST:** Identificação de genes reais diretamente da base de dados do NCBI.

---

## 🏗️ Estrutura do Projeto

```text
GeneFlux/
├── app.py                # Ponto de entrada (Router)
├── pyproject.toml        # Configuração do projeto e dependências
├── src/
│   ├── logic/            # O "Coração" do sistema
│   │   ├── dna_rna.py    # Bioinformática (Transcrição, Tradução, Métricas)
│   │   ├── genetics.py   # Genética (Mendel, ABO, Punnett)
│   │   ├── education.py  # Geradores de exercícios e simulações
│   │   └── ncbi.py       # Integração com API BLAST
│   ├── components/       # UI Reutilizável
│   │   ├── visualizations.py # Tabelas de Códons e Sequências
│   │   └── genetics_ui.py   # Renderização de Quadros de Punnett
│   └── pages/            # Layouts das Páginas
│       ├── analysis_page.py  # Ferramentas de análise profissional
│       └── education_page.py # Laboratório virtual de aprendizagem
```

---

## 🛠️ Tecnologias Utilizadas

- **Streamlit:** Interface web responsiva.
- **Biopython:** Processamento de sequências biológicas e acesso ao NCBI.
- **Python 3.12+:** Tipagem estática e lógica moderna.

---

## 🚀 Como Executar

1.  **Instale as dependências:**
    ```bash
    pip install streamlit biopython
    ```

2.  **Inicie a aplicação:**
    ```bash
    streamlit run app.py
    ```

---

## 📖 Funcionalidades

### 1. Módulo de Análises
- **DNA Complementar:** Gera a fita complementar 3'->5'.
- **Transcrição:** Converte DNA em mRNA.
- **Tradução:** Converte sequências em aminoácidos (3 letras).
- **Localizador de ORF:** Identifica Quadros de Leitura Aberta.
- **BLASTn:** Identifica espécies e genes correspondentes.

### 2. Módulo Educativo
- **Quadro de Punnett:** Simula cruzamentos genéticos com cálculo de frequências genotípicas e fenotípicas.
- **Mutação em Tempo Real:** Altere uma base no DNA e veja instantaneamente se a mutação é Silenciosa, Missense, Nonsense ou Frameshift.

---

## 👤 Sobre o Autor

**Rodrigo Mello**
*Biomédico, Mestre em Genética e Toxicologia, Doutorando em Biologia Celular e Molecular.*
Unindo a ciência da vida com o poder da programação Python.
