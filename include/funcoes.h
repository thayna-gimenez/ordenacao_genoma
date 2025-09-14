#ifndef FUNCOES_H
#define FUNCOES_H

#define MAX_SEQ_LENGTH 100
#define DNA_CHARS "ACGT"

#include <time.h>

// Leitura de Arquivo
char** ler_arquivo(char* filename, int* total_seqs);

// Gerar DNA
void generate_dna_sequence(char* seq, int length);
int compare_dna(const void* a, const void* b);
void sequential_sort(char** data, int n);

#endif