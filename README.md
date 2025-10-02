# Ordenação Paralela de Dados Genômicos com MPI

Implementação paralela de algoritmo de ordenação para grandes volumes de dados genômicos usando MPI, para comparação com código sequencial de ordenação quicksort.

**Gerar exemplos**  

```sh
gcc -Iincludes src/gerar_exemplos.c src/funcoes.c -o gerar_exemplos
```
```sh
# 100 mil sequências
./gerar_exemplos 100000 testes/entradas/ex_pequeno.txt
 
# 1 milhão de sequências
./gerar_exemplos 1000000 testes/entradas/ex_medio.txt

# 10 milhões de sequências
./gerar_exemplos 10000000 testes/entradas/ex_grande.txt
```
**Algoritmo sequencial**
```sh
mpicc -Iincludes src/sequencial.c src/funcoes.c -o codigo_sequencial

# Ordenar 10 milhões
./codigo_sequencial testes/entradas/ex_grande.txt testes/saidas/saida_grande.txt.txt

```
**Algoritmo paralelo**
```sh
mpicc -O3 -Iincludes src/paralelo.c src/funcoes.c -o codigo_paralelo

# Ordenar 10 milhões com 12 processadores
mpirun -np 12 ./codigo_paralelo testes/entradas/ex_grande.txt testes/saidas/saida_grande.txt


```
