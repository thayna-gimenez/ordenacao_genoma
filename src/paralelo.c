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
    distribuir_sequencias(dna_sequencias, &vetor_local, &buffer_local, &counter_local, total_seqs, rank, size);

    //Ordenação local
    sequential_sort(vetor_local, counter_local);

    // Escolhendo samples em cada processador
    char** samples_locais;

    if (counter_local > (size-1)){
        samples_locais = malloc((size - 1) * sizeof(char)); // vai ter m-1 samples
        for (int i = 0; i < size-1; i++){
            samples_locais[i] = vetor_local[((i+1) * counter_local) / size];
        }
    } else {
        samples_locais = malloc(counter_local * sizeof(char)); // vai pegar todos os elementos
        for (int i = 0; i < counter_local; i++){
            samples_locais[i] = vetor_local[((i+1) * counter_local) / size];
        }
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