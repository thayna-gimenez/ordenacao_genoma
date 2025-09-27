#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include "../include/funcoes.h"

// mpicc -Iincludes src/paralelo.c src/funcoes.c -o codigo_paralelo
// gcc -o codigo_sequencial src/sequencial.c
// ./codigo_paralelo testes/entradas/ex_minusculo.txt testes/saidas/saida_minusculo.txt

#define MAX_SEQ_LENGTH 100
#define DNA_CHARS "ACGT"

int main(int argc, char *argv[]) {
    int rank, size;

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    char* arq_entrada = argv[1];
    char* arq_saida = argv[2];

    int total_seqs;
    char** dna_sequencias = NULL;
    if (rank == 0) {
        // Leitura de Arquivo no processador 0    
        dna_sequencias = ler_arquivo(arq_entrada, &total_seqs);
    }

    MPI_Bcast(&total_seqs, 1, MPI_INT, 0, MPI_COMM_WORLD);

    // Cálculo dos tamanhos locais
    int base = total_seqs / size;
    int resto_local = total_seqs;
    int tamanho_local = base + (rank < resto_local ? 1 : 0);

    char (*vetor_local)[MAX_SEQ_LENGTH] = malloc(tamanho_local * MAX_SEQ_LENGTH);

    if (rank == 0){
        // Vetores de tamanho p para especificar o número de elementos de cada vetor e o displacement 
        int *sendcounts = (int*)malloc(size * sizeof(int));
        int *displs = (int*)malloc(size * sizeof(int));

        // Cálculo dos vetores
        int deslocamento = 0;
        for (int i = 0; i < size; i++){
            sendcounts[i] = base;
            
            if (i < resto_local){ 
                sendcounts[i] += 1 ;
            }
            
            displs[i] = deslocamento;
            deslocamento += sendcounts[i];
        }

        // Distribuição dos dados entre os processadores
        MPI_Scatterv(dna_sequencias, sendcounts, displs, MPI_CHAR, vetor_local, tamanho_local * MAX_SEQ_LENGTH, MPI_CHAR, 0, MPI_COMM_WORLD);

        free(dna_sequencias);
        free(sendcounts);
        free(displs);
    } else {
        MPI_Scatterv(NULL, NULL, NULL, MPI_CHAR,
                    vetor_local, tamanho_local * MAX_SEQ_LENGTH, MPI_CHAR,
                    0, MPI_COMM_WORLD);
    }

    printf("Processo %d recebeu %d elementos: ", rank, tamanho_local);
    for (int i = 0; i < tamanho_local; i++) {
        printf("%s ", vetor_local[i]);
    }
    printf("\n");

    // Arquivo de saída
    //escrever_arquivo(arq_saida, dna_sequencias, total_seqs);

    free(vetor_local);
    MPI_Finalize();
    return 0;
}