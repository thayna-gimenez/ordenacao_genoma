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
char** enviar_samples(char** samples_locais, int* n_amostras_por_processo, int tamanho_samples, int* qtde_amostras, int rank, int size);
char** selecionar_samples_globais(char** samples, int qtde_amostras, int* qtde_amostras_globais, int rank, int size);
int* particionar_dados(char** vetor_local, int counter_local, char** samples_globais, int qtde_samples_globais, int size, int rank);
int trocar_dados(char** vetor_local, int counter_local, int* sendcounts, char*** dados_recebidos, char** samples_globais, int num_pivos, int rank, int size);
char** coletar_resultados(char** dados_recebidos, int total_recebido, int total_seqs, int rank, int size);

#endif