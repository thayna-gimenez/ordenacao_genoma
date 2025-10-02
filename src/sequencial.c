// mpicc -Iincludes src/sequencial.c src/funcoes.c -o codigo_sequencial
// ./codigo_sequencial entrada.txt saida.txt

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "../include/funcoes.h"

#define MAX_SEQ_LENGTH 100
#define DNA_CHARS "ACGT"

int main(int argc, char *argv[]) {
    clock_t inicio, fim;
    
    if (argc != 3) {
        printf("Erro de comando.");
        return 1;
    }
    
    char* arq_entrada = argv[1];
    char* arq_saida = argv[2];

    // Leitura de Arquivo
    int total_seqs;
    char** dna_sequencias = ler_arquivo(arq_entrada, &total_seqs);
    
    inicio = clock();
    sequential_sort(dna_sequencias, total_seqs);
    fim = clock();

    double tempo = (double)(fim - inicio)/CLOCKS_PER_SEC;

    printf("%f", tempo);

    // Arquivo de saída
    escrever_arquivo(arq_saida, dna_sequencias, total_seqs);

    free(dna_sequencias);
    return 0;
}
