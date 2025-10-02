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

void sequential_sort_parallel(char** data, int n) {
    #pragma omp parallel
    {
        #pragma omp single nowait
        qsort(data, n, sizeof(char*), compare_dna);
    }
}

void distribuir_sequencias(char** dna_sequencias,
                           char*** vetor_local,
                           char** buffer_local,
                           int* counter_local,
                           int total_seqs,
                           int rank,
                           int size) {
    int base = total_seqs / size;
    int resto = total_seqs % size;
    *counter_local = base + (rank < resto ? 1 : 0);

    // Cada processo aloca o espaço local (contíguo)
    *buffer_local = malloc((*counter_local) * MAX_SEQ_LENGTH);
    *vetor_local = malloc((*counter_local) * sizeof(char*));
    for (int i = 0; i < *counter_local; i++) {
        (*vetor_local)[i] = &(*buffer_local)[i * MAX_SEQ_LENGTH];
    }

    // Criar tipo MPI para string fixa
    MPI_Datatype dna_type;
    MPI_Type_contiguous(MAX_SEQ_LENGTH, MPI_CHAR, &dna_type);
    MPI_Type_commit(&dna_type);

    if (rank == 0) {
        int *sendcounts = malloc(size * sizeof(int));
        int *displs = malloc(size * sizeof(int));
        int deslocamento = 0;

        for (int i = 0; i < size; i++) {
            int count_i = base + (i < resto ? 1 : 0);
            sendcounts[i] = count_i;
            displs[i] = deslocamento;
            deslocamento += count_i;
        }

        // Criar buffer contíguo apenas uma vez
        char *temp_buffer = malloc(total_seqs * MAX_SEQ_LENGTH);
        for (int i = 0; i < total_seqs; i++) {
            strncpy(&temp_buffer[i * MAX_SEQ_LENGTH], dna_sequencias[i], MAX_SEQ_LENGTH);
        }

        // Scatter usando o novo tipo
        MPI_Scatterv(temp_buffer, sendcounts, displs, dna_type,
                     *buffer_local, *counter_local, dna_type,
                     0, MPI_COMM_WORLD);

        free(temp_buffer);
        free(sendcounts);
        free(displs);
    } else {
        MPI_Scatterv(NULL, NULL, NULL, dna_type,
                     *buffer_local, *counter_local, dna_type,
                     0, MPI_COMM_WORLD);
    }

    MPI_Type_free(&dna_type);
}

char** escolher_samples(char** vetor_local, int counter_local, int size, int* tamanho_samples){
    char** samples_locais;
    
    if (counter_local > (size-1)){
        *tamanho_samples = size * (size - 1); // vai ter m-1 samples
    } else {
        *tamanho_samples = counter_local; // vai pegar todos os elementos
    }

    samples_locais = malloc(*tamanho_samples * sizeof(char*));
        for (int i = 0; i < *tamanho_samples; i++){
            samples_locais[i] = vetor_local[((i+1) * counter_local) / (*tamanho_samples + 1)];
        }

    return samples_locais;
}

int* coletar_n_amostras_no_root(int tamanho_samples, int rank, int size) {
    int *n_amostras_por_processo = NULL;
    
    if (rank == 0) {
        n_amostras_por_processo = malloc(size * sizeof(int));
    }
    
    // Cada processo envia seu valor para o root
    MPI_Gather(&tamanho_samples, 1, MPI_INT,
               n_amostras_por_processo, 1, MPI_INT,
               0, MPI_COMM_WORLD);
    
    return n_amostras_por_processo;  // Root retorna o array, outros retornam NULL
}

char** enviar_samples(char** samples_locais, int* n_amostras_por_processo, int tamanho_samples, int* qtde_amostras, int rank, int size)
{
    char **samples = NULL;
    char *recv_buffer = NULL;
    *qtde_amostras = 0;

    int *recvcounts = NULL;
    int *displs = NULL;

    if (rank == 0) {
        recvcounts = malloc(size * sizeof(int));
        displs = malloc(size * sizeof(int));

        // total de amostras e vetores de counts/displs em BYTES
        int offset = 0;
        for (int i = 0; i < size; i++) {
            int n_i = n_amostras_por_processo[i];
            recvcounts[i] = n_i * MAX_SEQ_LENGTH; // BYTES
            displs[i] = offset;                   // BYTES
            offset += recvcounts[i];
            *qtde_amostras += n_i;
        }

        samples = malloc((*qtde_amostras) * sizeof(char*));
        for (int i = 0; i < *qtde_amostras; i++) {
            samples[i] = malloc(MAX_SEQ_LENGTH * sizeof(char));
        }
        recv_buffer = malloc((*qtde_amostras) * MAX_SEQ_LENGTH);
    }

    // buffer contíguo de envio (BYTES)
    char *send_buffer = malloc((tamanho_samples) * MAX_SEQ_LENGTH);
    for (int i = 0; i < tamanho_samples; i++) {
        strncpy(&send_buffer[i * MAX_SEQ_LENGTH], samples_locais[i], MAX_SEQ_LENGTH);
    }

    // root coleta tudo via Gatherv; não-root não precisa conhecer counts/displs
    MPI_Gatherv(send_buffer, tamanho_samples * MAX_SEQ_LENGTH, MPI_CHAR,
                (rank == 0 ? recv_buffer : NULL),
                recvcounts, displs, MPI_CHAR, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        // desembrulha pro vetor de strings
        int current = 0;
        for (int proc = 0; proc < size; proc++) {
            int n_proc = n_amostras_por_processo[proc];
            int base = (displs ? displs[proc] : current * MAX_SEQ_LENGTH);
            for (int i = 0; i < n_proc; i++) {
                strncpy(samples[current], &recv_buffer[base + i * MAX_SEQ_LENGTH], MAX_SEQ_LENGTH);
                samples[current][MAX_SEQ_LENGTH - 1] = '\0';
                current++;
            }
        }
    }

    free(send_buffer);
    if (rank == 0) {
        free(recv_buffer);
        free(recvcounts);
        free(displs);
    }
    return samples;
}

char** selecionar_samples_globais(char** samples, int qtde_amostras, int* qtde_amostras_globais, int rank, int size){
    char** samples_globais = NULL;
    
    if (rank == 0) {
        sequential_sort_parallel(samples, qtde_amostras);
        samples_globais = malloc((*qtde_amostras_globais) * sizeof(char*));

        // Escolhendo m-1 samples globais igualmente espaçados dentro de cada processador
        for (int i = 0; i < (*qtde_amostras_globais); i++) {
            samples_globais[i] = malloc(MAX_SEQ_LENGTH * sizeof(char));
            int index = (i * (qtde_amostras - 1)) / (*qtde_amostras_globais);
            strncpy(samples_globais[i], samples[index], MAX_SEQ_LENGTH);
        }
    } else {
        samples_globais = malloc((*qtde_amostras_globais) * sizeof(char*));
        for (int i = 0; i < (*qtde_amostras_globais); i++) {
            samples_globais[i] = malloc(MAX_SEQ_LENGTH * sizeof(char));
        }
    }
    
    if (rank == 0) {
    char* global_buffer = malloc((*qtde_amostras_globais) * MAX_SEQ_LENGTH);
        for (int i = 0; i < (*qtde_amostras_globais); i++) {
            strncpy(&global_buffer[i * MAX_SEQ_LENGTH], samples_globais[i], MAX_SEQ_LENGTH);
        }
        MPI_Bcast(global_buffer, (*qtde_amostras_globais) * MAX_SEQ_LENGTH, MPI_CHAR, 0, MPI_COMM_WORLD);
        free(global_buffer);
    } else {
        char* global_buffer = malloc((*qtde_amostras_globais) * MAX_SEQ_LENGTH);
        MPI_Bcast(global_buffer, (*qtde_amostras_globais) * MAX_SEQ_LENGTH, MPI_CHAR, 0, MPI_COMM_WORLD);
        for (int i = 0; i < (*qtde_amostras_globais); i++) {
            strncpy(samples_globais[i], &global_buffer[i * MAX_SEQ_LENGTH], MAX_SEQ_LENGTH);
        }
        free(global_buffer);
    }

    return samples_globais;
}

int* particionar_dados(char** vetor_local, int counter_local, char** samples_globais, int qtde_amostras_globais, int size, int rank) {
    int* sendcounts = calloc(size, sizeof(int));
    
    // Ordenar localmente primeiro para busca binária
    sequential_sort_parallel(vetor_local, counter_local);
    
    // Usar busca binária para encontrar buckets
    for (int i = 0; i < counter_local; i++) {
        int left = 0, right = qtde_amostras_globais - 1;
        int bucket = size - 1; // default último bucket
        
        while (left <= right) {
            int mid = left + (right - left) / 2;
            int cmp = strcmp(vetor_local[i], samples_globais[mid]);
            
            if (cmp <= 0) {
                bucket = mid;
                right = mid - 1;
            } else {
                left = mid + 1;
            }
        }
        sendcounts[bucket]++;
    }
    
    return sendcounts;
}

int trocar_dados(char** vetor_local, int counter_local,
                 int* sendcounts, char*** dados_recebidos,
                 char** samples_globais, int num_pivos,
                 int rank, int size) {
    
    // Tipo MPI fixo
    MPI_Datatype dna_type;
    MPI_Type_contiguous(MAX_SEQ_LENGTH, MPI_CHAR, &dna_type);
    MPI_Type_commit(&dna_type);

    // Quantos cada processo vai receber
    int *recvcounts = malloc(size * sizeof(int));
    MPI_Alltoall(sendcounts, 1, MPI_INT, recvcounts, 1, MPI_INT, MPI_COMM_WORLD);

    // Deslocamentos (em elementos)
    int *sdispls = malloc(size * sizeof(int));
    int *rdispls = malloc(size * sizeof(int));
    sdispls[0] = 0; rdispls[0] = 0;
    for (int i = 1; i < size; i++) {
        sdispls[i] = sdispls[i-1] + sendcounts[i-1];
        rdispls[i] = rdispls[i-1] + recvcounts[i-1];
    }
    int total_recv = rdispls[size-1] + recvcounts[size-1];

    // Buffers
    char *send_buffer = malloc(counter_local * MAX_SEQ_LENGTH);
    char *recv_buffer = malloc(total_recv * MAX_SEQ_LENGTH);

    // Preencher send_buffer já no formato certo
    int *bucket_index = calloc(size, sizeof(int));
    for (int i = 0; i < counter_local; i++) {
        int b = 0;
        while (b < num_pivos && strcmp(vetor_local[i], samples_globais[b]) > 0) b++;
        int pos = sdispls[b] + bucket_index[b]++;
        strncpy(&send_buffer[pos * MAX_SEQ_LENGTH], vetor_local[i], MAX_SEQ_LENGTH);
    }

    // Alltoallv com tipo definido
    MPI_Alltoallv(send_buffer, sendcounts, sdispls, dna_type,
                  recv_buffer, recvcounts, rdispls, dna_type,
                  MPI_COMM_WORLD);

    // Reconstruir ponteiros
    *dados_recebidos = malloc(total_recv * sizeof(char*));
    for (int i = 0; i < total_recv; i++) {
        (*dados_recebidos)[i] = malloc(MAX_SEQ_LENGTH);
        strncpy((*dados_recebidos)[i], &recv_buffer[i * MAX_SEQ_LENGTH], MAX_SEQ_LENGTH);
        (*dados_recebidos)[i][MAX_SEQ_LENGTH-1] = '\0';
    }

    free(bucket_index);
    free(send_buffer);
    free(recv_buffer);
    free(sdispls);
    free(rdispls);
    free(recvcounts);

    MPI_Type_free(&dna_type);
    return total_recv;
}



char** coletar_resultados(char** dados_recebidos, int total_recebido,
                          int total_seqs, int rank, int size)
{
    char **resultado_final = NULL;
    int *recvcounts = NULL, *displs = NULL;


    if (rank == 0) {
        recvcounts = malloc(size * sizeof(int));
        displs = malloc(size * sizeof(int));
        resultado_final = malloc(total_seqs * sizeof(char*));
        for (int i = 0; i < total_seqs; i++) {
            resultado_final[i] = malloc(MAX_SEQ_LENGTH);
        }
    }

    MPI_Gather(&total_recebido, 1, MPI_INT,
               (rank == 0 ? recvcounts : NULL), 1, MPI_INT,
               0, MPI_COMM_WORLD);

    char *recv_buffer = NULL;
    if (rank == 0) {
        displs[0] = 0;
        for (int i = 1; i < size; i++) {
            displs[i] = displs[i-1] + recvcounts[i-1] * MAX_SEQ_LENGTH; // BYTES
        }
        for (int i = 0; i < size; i++) {
            recvcounts[i] = recvcounts[i] * MAX_SEQ_LENGTH;             // BYTES
        }
        recv_buffer = malloc((size_t)total_seqs * MAX_SEQ_LENGTH);
    }

    char *send_buffer = malloc((size_t)total_recebido * MAX_SEQ_LENGTH);
    for (int i = 0; i < total_recebido; i++) {
        strncpy(&send_buffer[(size_t)i * MAX_SEQ_LENGTH], dados_recebidos[i], MAX_SEQ_LENGTH);
    }

    MPI_Gatherv(send_buffer, total_recebido * MAX_SEQ_LENGTH, MPI_CHAR,
                recv_buffer, recvcounts, displs, MPI_CHAR,
                0, MPI_COMM_WORLD);

    if (rank == 0) {
        int current = 0;
        for (int proc = 0; proc < size; proc++) {
            int n_el = recvcounts[proc] / MAX_SEQ_LENGTH; // volta a ELEMENTOS
            int base = displs[proc];
            for (int i = 0; i < n_el; i++) {
                strncpy(resultado_final[current],
                        &recv_buffer[base + i * MAX_SEQ_LENGTH],
                        MAX_SEQ_LENGTH);
                resultado_final[current][MAX_SEQ_LENGTH - 1] = '\0';
                current++;
            }
        }
        free(recv_buffer);
        free(recvcounts);
        free(displs);
    }

    free(send_buffer);
    return resultado_final;
}
