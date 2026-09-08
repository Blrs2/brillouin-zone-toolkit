# Brillouin-zone toolkit — Fortran 90

Implementação em **Fortran 90** dos cálculos de rede recíproca e primeira zona
de Brillouin usados no projeto de graduação. O núcleo numérico não depende de
Python, NumPy, bibliotecas de álgebra linear ou pacotes de cristalografia.

![Galeria das zonas de Brillouin](docs/figures/gallery.png)

## O que o código calcula

Para uma base direta `a_i`, o módulo `brillouin90` calcula a base recíproca
`b_i`, seguindo `a_i · b_j = 2π δ_ij`, e constrói a célula de Wigner–Seitz da
rede recíproca a partir das meias-espaços

\[
  \mathbf G\cdot\mathbf k \leq \frac{|\mathbf G|^2}{2}.
\]

O algoritmo enumera interseções de trios de planos, testa cada ponto contra
todos os planos, identifica as faces coplanares e calcula o volume por
tetraedrização em torno da origem. As operações de simetria cúbica e a redução
de pontos por órbitas também estão no módulo Fortran.

Os detalhes numéricos e as escolhas de portabilidade estão em
[`docs/algorithm.md`](docs/algorithm.md).

As quatro redes de referência são construídas diretamente no código:

| rede direta | poliedro da primeira zona | vértices | faces |
| --- | --- | ---: | ---: |
| SC | cubo | 8 | 6 |
| FCC | octaedro truncado | 24 | 14 |
| BCC | dodecaedro rômbico | 14 | 12 |
| hexagonal | prisma hexagonal | 12 | 8 |

O volume da zona coincide, dentro da tolerância numérica, com
`(2π)^3 / V_célula` em todos os casos.

## Compilar e testar

É necessário um compilador Fortran com suporte a Fortran 90; o fluxo abaixo usa
`gfortran`:

```bash
make test
```

O teste compila com verificações de limites (`-fcheck=all`) e valida dualidade
direta–recíproca, topologia, volume, contenção da origem, grupo pontual cúbico
de 48 operações e redução de uma órbita de pontos.

Para executar o programa demonstrativo:

```bash
make build
./build/brillouin_demo
```

A saída inclui uma linha por rede com número de vértices, faces, volume
calculado e volume da célula recíproca. A CI do GitHub recompila e executa os
dois binários em cada alteração.

## Estrutura

```text
fortran/brillouin90.f90       módulo Fortran 90 (núcleo numérico)
programs/brillouin_demo.f90   exemplo executável e regressão de volumes
tests_fortran/test_brillouin.f90  testes executáveis sem dependências externas
docs/figures/                 figuras e resumo visual do projeto
legacy/python_reference/      reconstrução Python anterior, somente referência
```

O diretório `legacy/python_reference` foi mantido para preservar a história do
repositório e facilitar comparação numérica; ele não é a implementação
principal, não é compilado pela CI e não deve ser citado como o código original
do projeto. A implementação mantida é a Fortran 90 em `fortran/`.

## Escopo e proveniência

Este repositório organiza uma reimplementação Fortran 90 reproduzível do método
e dos resultados geométricos associados ao projeto de graduação. O PDF do
relatório acadêmico permanece fora do repositório público até confirmação das
condições de redistribuição. As figuras em `docs/figures/` são artefatos
visuais selecionados do estudo; o código fonte atual é independente delas.

O código está sob a licença MIT. O relatório e quaisquer materiais acadêmicos
originais continuam sujeitos às suas próprias condições de autoria e uso.

## English summary

The maintained implementation is standard Fortran 90: reciprocal-lattice
construction, half-space first-Brillouin-zone geometry, polyhedron volumes,
cubic point-group operations and orbit reduction. The former Python version is
archived under `legacy/python_reference`; it is not the project’s primary
implementation.
