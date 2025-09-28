#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "../include/funcoes.h"

#define MAX_SEQ_LENGTH 100

int main(int argc, char *argv[]) {
    int rank, size, total_seqs;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    char* arq_entrada = argv[1];
    char* arq_saida = argv[2];
    char** dna_sequencias = NULL;
    
    if (rank == 0) {
        dna_sequencias = ler_arquivo(arq_entrada, &total_seqs);
    }

    MPI_Bcast(&total_seqs, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Enviando elementos para cada processador
    char** vetor_local = NULL;
    char* buffer_local = NULL;
    int counter_local;
    int tamanho_samples;
    distribuir_sequencias(dna_sequencias, &vetor_local, &buffer_local, &counter_local, total_seqs, rank, size);

    //Ordenação local
    sequential_sort(vetor_local, counter_local);

    // Escolhendo samples em cada processador
    char** samples_locais;

    if (counter_local > (size-1)){
        tamanho_samples = size - 1;
        samples_locais = malloc((tamanho_samples) * sizeof(char)); // vai ter m-1 samples
        for (int i = 0; i < size-1; i++){
            samples_locais[i] = vetor_local[((i+1) * counter_local) / size];
        }
    } else {
        tamanho_samples = counter_local;
        samples_locais = malloc(tamanho_samples * sizeof(char)); // vai pegar todos os elementos
        for (int i = 0; i < counter_local; i++){
            samples_locais[i] = vetor_local[i];
        }
    }

    // Enviando os samples para o processador 0
    int* recvcounts = NULL;
    int* displs = NULL;
    int qtde_amostras = tamanho_samples * size;
    char** samples = NULL;

    // Definindo tamanho do que é necessário para o Gatherv
    if (rank == 0){
        recvcounts = malloc(size * sizeof(int));
        displs = malloc(size * sizeof(int));
        
        for (int i = 0; i < size; i++){
            recvcounts[i] = tamanho_samples * MAX_SEQ_LENGTH;
            displs[i] = i * tamanho_samples * MAX_SEQ_LENGTH;
        }
        
        samples = malloc(qtde_amostras * sizeof(char*));
        for (int i = 0; i < qtde_amostras; i++) {
            samples[i] = malloc(MAX_SEQ_LENGTH * sizeof(char));
        }
    }

    // buffer temporário pra enviar no gather
    char* send_buffer = malloc(tamanho_samples * MAX_SEQ_LENGTH * sizeof(char));
    for (int i = 0; i < tamanho_samples; i++) {
        strncpy(&send_buffer[i * MAX_SEQ_LENGTH], samples_locais[i], MAX_SEQ_LENGTH);
    }

    // buffer temporário pra receber no gather (apenas processador 0)
    char *recv_buffer = NULL;
    if (rank == 0) {
        recv_buffer = malloc(qtde_amostras * MAX_SEQ_LENGTH * sizeof(char));
    }

    MPI_Gatherv(send_buffer, 
            tamanho_samples * MAX_SEQ_LENGTH, 
            MPI_CHAR, 
            recv_buffer, 
            recvcounts, 
            displs, 
            MPI_CHAR, 
            0, 
            MPI_COMM_WORLD);

    if (rank == 0) {
        for (int proc = 0; proc < size; proc++) {
            for (int i = 0; i < tamanho_samples; i++) {
                int index_global = proc * tamanho_samples + i;
                strncpy(samples[index_global], 
                    &recv_buffer[displs[proc] + i * MAX_SEQ_LENGTH], 
                    MAX_SEQ_LENGTH);
                samples[index_global][MAX_SEQ_LENGTH - 1] = '\0';
            }
        }
        printf("%d\n", tamanho_samples);

        // printf("Processo %d: %d strings\n", rank, counter_local);
        for (int i = 0; i < qtde_amostras; i++) {
            printf("%s\n", samples[i]);
        }
    }

    free(send_buffer);
    if (rank == 0) {
        free(recv_buffer);
        free(recvcounts);
        free(displs);
    }


    // printf("Processo %d: %d strings\n", rank, counter_local);
    // for (int i = 0; i < size-1; i++) {
    //     printf("%s\n", samples_locais[i]);
    // }


    // Liberando memória
    free(vetor_local); 
    free(buffer_local);
    free(samples_locais);
    
    if (rank == 0) {
        for (int i = 0; i < total_seqs; i++) free(dna_sequencias[i]);
        free(dna_sequencias);
    }
    
    MPI_Finalize();
    return 0;
}   