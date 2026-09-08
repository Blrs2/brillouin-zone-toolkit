# Algoritmo geométrico em Fortran 90

## 1. Rede recíproca

As três linhas da matriz `direct` são os vetores primitivos da rede direta.
Usando o produto vetorial, o módulo constrói

\[
 \mathbf b_1=2\pi\frac{\mathbf a_2\times\mathbf a_3}{V},\quad
 \mathbf b_2=2\pi\frac{\mathbf a_3\times\mathbf a_1}{V},\quad
 \mathbf b_3=2\pi\frac{\mathbf a_1\times\mathbf a_2}{V},
\]

onde `V = a1 · (a2 × a3)`. O teste Fortran verifica
`A Bᵀ = 2π I` antes de construir cada zona.

## 2. Planos de Wigner–Seitz

Para cada combinação inteira `(h,k,l)` em `[-search,search]^3`, exceto a
origem, calcula-se `G = h b1 + k b2 + l b3` e registra-se o plano

```text
dot(G, k) = dot(G, G) / 2
```

A região permitida é o lado que contém `k = 0`. Para `search = 2` são apenas
124 planos candidatos; isso é suficiente para as quatro redes de referência.

## 3. Vértices, faces e volume

Cada trio de planos define, quando linearmente independente, uma solução de um
sistema 3×3. O ponto é aceito como vértice se satisfaz todos os outros
semi-espaços. Vértices dentro da tolerância são fundidos.

Para cada plano que contém pelo menos três vértices, os pontos são projetados
em uma base ortonormal do plano e ordenados por `atan2`. A face poligonal é
triangulada com a origem como ponto interno; a soma dos volumes dos tetraedros

\[
  V_{tetra} = \frac{|v_0\cdot(v_i\times v_{i+1})|}{6}
\]

fornece o volume da zona.

## 4. Representação compatível com Fortran 90

Fortran 90 não possui componentes alocáveis em tipos derivados. Para manter o
código compilável por um compilador F90 simples, `zone_t` usa limites fixos
(`max_planes`, `max_vertices`, `max_faces`). Os limites são validados durante a
construção e retornam um código de erro se forem excedidos. Isso deixa explícito
o compromisso entre portabilidade e generalidade; para células muito
anisotrópicas, aumente `search` e, se necessário, os parâmetros máximos.

## 5. Simetria

`cubic_operations` gera as 48 matrizes de permutação de eixos com sinais
independentes, incluindo operações próprias e impróprias. `reduce_points`
aplica cada operação, escolhe o ponto lexicograficamente mínimo da órbita e
agrupa representantes dentro da tolerância. O teste executável confirma que
`(1,0,0)`, `(-1,0,0)` e `(0,1,0)` pertencem a uma única órbita cúbica.

## 6. Validação

O binário de teste compara, para SC/FCC/BCC/hexagonal:

- dualidade direta–recíproca;
- número de vértices e faces conhecidos;
- volume da zona contra `(2π)^3 / |det(A)|`;
- contenção da origem;
- ordem do grupo pontual cúbico;
- redução de órbitas.

Esses testes são numéricos e não substituem uma biblioteca geral de grupos
espaciais. Eles cobrem o escopo do projeto de graduação e mantêm o cálculo
legível para inspeção em Fortran 90.
