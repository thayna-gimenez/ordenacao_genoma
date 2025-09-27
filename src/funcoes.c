#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include "../include/funcoes.h"
#define MAX_SEQ_LENGTH 100

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

void distribuir_sequencias(char** dna_sequencias, char*** vetor_local, char** buffer_local, int* counter_local, int total_seqs, int rank, int size){
    int base = total_seqs / size;
    int resto = total_seqs % size;
    *counter_local = base + (rank < resto ? 1 : 0);

    *vetor_local = malloc((*counter_local) * sizeof(char*));
    *buffer_local = malloc((*counter_local) * MAX_SEQ_LENGTH);

    for (int i = 0; i < *counter_local; i++) {
        (*vetor_local)[i] = &(*buffer_local)[i * MAX_SEQ_LENGTH];
    }

    if (rank == 0) {
        int *sendcounts = malloc(size * sizeof(int));
        int *displs = malloc(size * sizeof(int));
        int deslocamento = 0;
        
        for (int i = 0; i < size; i++) {
            int count_i = base + (i < resto ? 1 : 0);
            sendcounts[i] = count_i * MAX_SEQ_LENGTH;
            displs[i] = deslocamento * MAX_SEQ_LENGTH;
            deslocamento += count_i;
        }
        
        char *temp_buffer = malloc(total_seqs * MAX_SEQ_LENGTH);
        for (int i = 0; i < total_seqs; i++) {
            strncpy(&temp_buffer[i * MAX_SEQ_LENGTH], dna_sequencias[i], MAX_SEQ_LENGTH);
        }
        
        MPI_Scatterv(temp_buffer, sendcounts, displs, MPI_CHAR,
                    *buffer_local, (*counter_local) * MAX_SEQ_LENGTH, MPI_CHAR, 0, MPI_COMM_WORLD);
        
        free(temp_buffer);
        free(sendcounts);
        free(displs);
    } else {
        MPI_Scatterv(NULL, NULL, NULL, MPI_CHAR, *buffer_local, (*counter_local) * MAX_SEQ_LENGTH, MPI_CHAR, 0, MPI_COMM_WORLD);

    }
}