// mpicc -Iincludes src/paralelo.c src/funcoes.c -o codigo_paralelo
// mpirun -np 4 ./codigo_paralelo testes/entradas/ex_minusculo.txt testes/saidas/saida_minusculo.txt

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
    char** samples_locais = escolher_samples(vetor_local, counter_local, size, &tamanho_samples);
    int* n_amostras_por_processo = coletar_n_amostras_no_root(tamanho_samples, rank, size);

    // Enviando os samples para o processador 0
    int* recvcounts = NULL;
    int* displs = NULL;
    char** samples = NULL; // todas as amostras
    int qtde_amostras = 0;
    char *recv_buffer = NULL;

    //printf("%d\n", n_amostras_por_processo[rank]);

    // Definindo tamanho do que é necessário para o Gatherv
    if (rank == 0){
        recvcounts = malloc(size * sizeof(int));
        displs = malloc(size * sizeof(int));
        
        qtde_amostras = 0;
        for (int i = 0; i < size; i++){
            recvcounts[i] = n_amostras_por_processo[i] * MAX_SEQ_LENGTH;
            displs[i] = qtde_amostras * MAX_SEQ_LENGTH;
            qtde_amostras += n_amostras_por_processo[i];
        }
        
        samples = malloc(qtde_amostras * sizeof(char*));
        for (int i = 0; i < qtde_amostras; i++) {
            samples[i] = malloc(MAX_SEQ_LENGTH * sizeof(char));
        }

        recv_buffer = malloc(qtde_amostras * MAX_SEQ_LENGTH * sizeof(char));
    } else {
        // Outros processos apenas alocam memória para receber
        recvcounts = malloc(size * sizeof(int));
        displs = malloc(size * sizeof(int));
    }

    MPI_Bcast(recvcounts, size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(displs, size, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&qtde_amostras, 1, MPI_INT, 0, MPI_COMM_WORLD);

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

    // Liberando memória
    free(vetor_local); 
    free(buffer_local);
    
    free(samples_locais);
    free(send_buffer);
    free(recvcounts);
    free(displs);

    if (rank == 0) {
        free(recv_buffer);
    }
    
    if (rank == 0) {
        for (int i = 0; i < total_seqs; i++) free(dna_sequencias[i]);
        free(dna_sequencias);
    }
    
    MPI_Finalize();
    return 0;
}   