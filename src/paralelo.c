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

    int qtde_amostras;
    char** samples = enviar_samples(samples_locais, n_amostras_por_processo, tamanho_samples, &qtde_amostras, rank, size);
    MPI_Bcast(&qtde_amostras, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    

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