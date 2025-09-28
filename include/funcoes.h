#ifndef FUNCOES_H
#define FUNCOES_H

#define MAX_SEQ_LENGTH 100
#define DNA_CHARS "ACGT"

#include <time.h>

// Arquivo
char** ler_arquivo(char* filename, int* total_seqs);
void escrever_arquivo(char* nome_arquivo, char** dna_sequencia, int total_seqs);

// Gerar DNA
void generate_dna_sequence(char* seq, int length);
int compare_dna(const void* a, const void* b);

// Função sequencial
void sequential_sort(char** data, int n);

// Funções paralelas
void distribuir_sequencias(char** dna_sequencias, char*** vetor_local, char** buffer_local, int* counter_local, int total_seqs, int rank, int size);
char** escolher_samples(char** vetor_local, int counter_local, int size, int* tamanho_samples);
int* coletar_n_amostras_no_root(int meu_n_amostras, int rank, int size);
#endif