// gcc -Iincludes src/gerar_exemplos.c src/funcoes.c -o gerar_exemplos
// ./gerar_exemplos 100000 testes/entradas/ex_pequeno.txt
// ./gerar_exemplos 1000000 testes/entradas/ex_medio.txt
// ./gerar_exemplos 10000000 testes/entradas/ex_grande.txt


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "../include/funcoes.h"

#define MAX_SEQ_LENGTH 100
#define DNA_CHARS "ACGT"

int main(int argc, char *argv[]) {
    
    if (argc != 3) {
        printf("Comando errado.");
        return 1;
    }
    
    int num_seqs = atoi(argv[1]);
    
    // Geração aleatória de keys com base na hora
    //srand(time(0));
    
    // Abertura de arquivo txt
    FILE* arquivo = fopen(argv[2], "w");
    if (!arquivo) {
        printf("Erro ao criar o arquivo txt.");
        return 1;
    }
    
    for (int i = 0; i < num_seqs; i++) {
        int tamanho = 10 + (rand() % 91); // gera valores de 10 a 100
        char sequencia[MAX_SEQ_LENGTH + 1];
        generate_dna_sequence(sequencia, tamanho);
        fprintf(arquivo, "%s\n", sequencia);
    }

    fclose(arquivo);
    return 0;
}