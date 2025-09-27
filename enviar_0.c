// Enviando os samples para o processador 0
    int* recvcounts = NULL;
    int* displs = NULL;
    int qtde_amostras = tamanho_samples * size;
    char** samples;

    // Definindo tamanho do que é necessário para o Gatherv
    if (rank == 0){
        recvcounts = malloc(size * sizeof(char));
        displs = malloc(size * sizeof(char));
        for (int i = 0; i < size; i++){
            recvcounts[i] = tamanho_samples * MAX_SEQ_LENGTH; // vetor em que cada posição é a quantidade de amostras que cada processador envia 
            displs[i] = i * tamanho_samples * MAX_SEQ_LENGTH; // deslocamento no buffer
        }
        
        samples = malloc(qtde_amostras * sizeof(char*)); // recebe todas as samples de cada processador
        for (int i = 0; i < qtde_amostras; i++) {
            samples[i] = malloc(MAX_SEQ_LENGTH * sizeof(char));
        }
    }
    
    // buffer temporário pra enviar no gather
    char* send_buffer = malloc(tamanho_samples * MAX_SEQ_LENGTH * sizeof(char));
    for (int i = 0; i < tamanho_samples; i++) {
        strncpy(&send_buffer[i * MAX_SEQ_LENGTH], samples[i], MAX_SEQ_LENGTH);
    }

    // buffer temporário pra receber no gather (apenas processador 0)
    char *recv_buffer = NULL;
    if (rank == 0) {
        recv_buffer = malloc(qtde_amostras * MAX_SEQ_LENGTH * sizeof(char));
    }
    
    MPI_Gatherv(send_buffer, tamanho_samples * MAX_SEQ_LENGTH, MPI_CHAR, recv_buffer, recvcounts, displs, MPI_CHAR, 0, MPI_COMM_WORLD);
    
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
    }
