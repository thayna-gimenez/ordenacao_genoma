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

char** escolher_samples(char** vetor_local, int counter_local, int size, int* tamanho_samples){
    char** samples_locais; 
    
    if (counter_local > (size-1)){
        *tamanho_samples = size - 1; // vai ter m-1 samples
        samples_locais = malloc(*tamanho_samples * sizeof(char*));
        for (int i = 0; i < size-1; i++){
            samples_locais[i] = vetor_local[((i+1) * counter_local) / size];
        }
    } else {
        *tamanho_samples = counter_local; // vai pegar todos os elementos
        samples_locais = malloc(*tamanho_samples * sizeof(char*));
        for (int i = 0; i < counter_local; i++){
            samples_locais[i] = vetor_local[i];
        }
    }

    return samples_locais;
}

int* coletar_n_amostras_no_root(int meu_n_amostras, int rank, int size) {
    int *n_amostras_por_processo = NULL;
    
    if (rank == 0) {
        n_amostras_por_processo = malloc(size * sizeof(int));
    }
    
    // Cada processo envia seu valor para o root
    MPI_Gather(&meu_n_amostras, 1, MPI_INT,
               n_amostras_por_processo, 1, MPI_INT,
               0, MPI_COMM_WORLD);
    
    return n_amostras_por_processo;  // Root retorna o array, outros retornam NULL
}

char** enviar_samples(char** samples_locais, int* n_amostras_por_processo, int tamanho_samples, int* qtde_amostras, int rank, int size){
        // Enviando os samples para o processador 0
    int* recvcounts = NULL;
    int* displs = NULL;
    char** samples = NULL; // todas as amostras
    char *recv_buffer = NULL;
    *qtde_amostras = 0;

    //printf("%d\n", n_amostras_por_processo[rank]);

    // Definindo tamanho do que é necessário para o Gatherv
    if (rank == 0){
        recvcounts = malloc(size * sizeof(int));
        displs = malloc(size * sizeof(int));
        
        (*qtde_amostras) = 0;
        for (int i = 0; i < size; i++){
            recvcounts[i] = n_amostras_por_processo[i] * MAX_SEQ_LENGTH;
            displs[i] = (*qtde_amostras) * MAX_SEQ_LENGTH;
            (*qtde_amostras) += n_amostras_por_processo[i];
        }
        
        samples = malloc((*qtde_amostras) * sizeof(char*));
        for (int i = 0; i < (*qtde_amostras); i++) {
            samples[i] = malloc(MAX_SEQ_LENGTH * sizeof(char));
        }

        recv_buffer = malloc((*qtde_amostras) * MAX_SEQ_LENGTH * sizeof(char));
    } else {
        // Outros processos apenas alocam memória para receber
        recvcounts = malloc(size * sizeof(int));
        displs = malloc(size * sizeof(int));
    }

    MPI_Bcast(recvcounts, size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(displs, size, MPI_INT, 0, MPI_COMM_WORLD);

    // Cria buffer contíguo para envio
    char *send_buffer = malloc(tamanho_samples * MAX_SEQ_LENGTH * sizeof(char));
    for (int i = 0; i < tamanho_samples; i++) {
        strncpy(&send_buffer[i * MAX_SEQ_LENGTH], samples_locais[i], MAX_SEQ_LENGTH);
    }

    // Coleta todas as amostras no root usando Gatherv
    MPI_Gatherv(send_buffer, tamanho_samples * MAX_SEQ_LENGTH, MPI_CHAR, (rank == 0) ? recv_buffer : NULL, recvcounts, displs, MPI_CHAR, 0, MPI_COMM_WORLD);

    // printf("Processo %d: %d strings\n", rank, counter_local);
    // for (int i = 0; i < size-1; i++) {
    //     printf("%s\n", samples_locais[i]);
    // }

    if (rank == 0) {
        int current_index = 0;
        for (int proc = 0; proc < size; proc++) {
            int n_proc = (n_amostras_por_processo != NULL) ? 
                         n_amostras_por_processo[proc] : 
                         recvcounts[proc] / MAX_SEQ_LENGTH;
            
            for (int i = 0; i < n_proc; i++) {
                strncpy(samples[current_index], &recv_buffer[displs[proc] + i * MAX_SEQ_LENGTH], MAX_SEQ_LENGTH);
                samples[current_index][MAX_SEQ_LENGTH - 1] = '\0';
                current_index++;
            }
        }
    }

    free(send_buffer);
    free(recvcounts);
    free(displs);

    if (rank == 0) {
        free(recv_buffer);
    }

    return samples;
}