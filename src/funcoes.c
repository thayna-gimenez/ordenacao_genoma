#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "../include/funcoes.h"

char** ler_arquivo(char* nome_arquivo, int* total_seqs){

    // Abertura do arquivo
    FILE* arquivo = fopen(nome_arquivo, "r");
    if (!arquivo) {
        printf("Erro ao abrir o arquivo.\n");
        return NULL;
    }
    
    // Contador de linhas
    *total_seqs = 0;
    char buffer[MAX_SEQ_LENGTH + 2];
    
    while (fgets(buffer, sizeof(buffer), arquivo)) {
        if (buffer[0] == '\n' || buffer[0] == '\0') {
            continue;
        }
        (*total_seqs)++;
    }

    // Volta pro início do arquivo (complementa o código acima)
    rewind(arquivo);

    // Alocacação de memória para armazenar as sequências
    char** sequencia = (char**)malloc(*total_seqs * sizeof(char*));
    if (!sequencia) {
        fclose(arquivo);
        printf("Erro ao alocar memória\n");
        return NULL;
    }
    
    // Leitura e gravação das sequências
    int cont = 0;
    while (fgets(buffer, sizeof(buffer), arquivo) && cont < *total_seqs) {
        // Ignora linhas vazias
        if (buffer[0] == '\n' || buffer[0] == '\0') {
            continue;
        }

        // Removendo quebra de linha que o fgets coloca
        buffer[strcspn(buffer, "\n")] = '\0';

        // Alocação de memória para as sequências
        sequencia[cont] = (char*)malloc((strlen(buffer) + 1) * sizeof(char));
        if (!sequencia[cont]) { // se houver erro
            printf("Erro de alocação de memória.");
            for (int i = 0; i < cont; i++) {
                free(sequencia[i]);
            }
        
            free(sequencia);
            fclose(arquivo);
            return NULL;
        }
    
        strcpy(sequencia[cont], buffer);
        cont++;
    }

    return sequencia;
}

void escrever_arquivo(char* nome_arquivo, char** dna_sequencia, int total_seqs){
    FILE* arquivo = fopen(nome_arquivo, "w");

    for (int i = 0; i < total_seqs; i++) {
        fprintf(arquivo, "%s\n", dna_sequencia[i]);
    }
    
    fclose(arquivo);
}

void generate_dna_sequence(char* seq, int length) {
    for (int i = 0; i < length; i++) {
        seq[i] = DNA_CHARS[rand() % 4];
    }
    seq[length] = '\0';
}

int compare_dna(const void* a, const void* b) {
    return strcmp(*(char**)a, *(char**)b);
}

void swap_dna(char **a, char **b) {
    char *temp = *a;
    *a = *b;
    *b = temp;
}

void sequential_sort(char** data, int n) {
    qsort(data, n, sizeof(char*), compare_dna);
}