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
    char** vetor_local = NULL;
    char* buffer_local = NULL;
    
    if (rank == 0) {
        dna_sequencias = ler_arquivo(arq_entrada, &total_seqs);
    }

    MPI_Bcast(&total_seqs, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Enviando elementos para cada processador
    int base = total_seqs / size;
    int resto = total_seqs % size;
    int counter_local = base + (rank < resto ? 1 : 0);

    vetor_local = malloc(counter_local * sizeof(char*));
    buffer_local = malloc(counter_local * MAX_SEQ_LENGTH);

    for (int i = 0; i < counter_local; i++) {
        vetor_local[i] = &buffer_local[i * MAX_SEQ_LENGTH];
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
                    buffer_local, counter_local * MAX_SEQ_LENGTH, MPI_CHAR, 0, MPI_COMM_WORLD);
        
        free(temp_buffer);
        free(sendcounts);
        free(displs);
    } else {
        MPI_Scatterv(NULL, NULL, NULL, MPI_CHAR,
                    buffer_local, counter_local * MAX_SEQ_LENGTH, MPI_CHAR, 0, MPI_COMM_WORLD);
    }
    
    //Ordenação local
    sequential_sort(vetor_local, counter_local);

    // Escolhendo samples em cada processador
    char** samples_locais;

    printf("Processo %d: %d strings\n", rank, counter_local);
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
    
    
    
    for (int i = 0; i < counter_local; i++) {
        printf("  [P%d] %s\n", rank, vetor_local[i]);
    }
    // for (int i = 0; i < size-1; i++) {
    //     printf("%s\n", samples_locais[i]);
    // }

    // Liberando memória
    free(samples_locais); 
    free(vetor_local); 
    free(buffer_local); 
    
    if (rank == 0) {
        for (int i = 0; i < total_seqs; i++) free(dna_sequencias[i]);
        free(dna_sequencias);
    }
    
    MPI_Finalize();
    return 0;
}   