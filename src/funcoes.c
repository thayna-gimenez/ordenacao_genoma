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
    
    // Volta pro início do arquivo (complementa o código acima)
    rewind(arquivo);

    // Alocacação de memória para armazenar as sequências
    char** sequencia = (char**)malloc(*total_seqs * sizeof(char*));
    if (!sequencia) {
        fclose(arquivo);
        printf("Erro ao alocar memória\n");
        return NULL;
    }
    
    // Leitura das sequências
    int cont = 0;
    while (fgets(buffer, sizeof(buffer), arquivo) && cont < *total_seqs) {
        // Ignora linhas vazias
        if (buffer[0] == '\n' || buffer[0] == '\0') {
            continue;
        }
    }
    
    // fgets adiciona uma quebra de linha, removendo...
    buffer[strcspn(buffer, "\n")] = '\0';

    // Gravação das sequências
    sequencia[cont] = (char*)malloc(((strlen(buffer) + 1) * sizeof(char)));
    if (!sequencia[cont]) { // se houver erro
        printf("Erro de alocação de memória");
        for (int i = 0; i < cont; i++) {
            free(sequencia[cont]);
        }
        
        free(sequencia);
        fclose(arquivo);
        return NULL;
    }
    
    strcpy(sequencia[cont], buffer);
    return sequencia;
}

void generate_dna_sequence(char* seq, int length) {
    for (int i = 0; i < length; i++) {
        seq[i] = DNA_CHARS[rand() % 4];
    }
    seq[length] = '\0';
}

int compare_dna(const void* a, const void* b) {
    return strcmp((char*)a, (char*)b);
}

void sequential_sort(char** data, int n) {
    qsort(data, n, sizeof(char*), compare_dna);
}
