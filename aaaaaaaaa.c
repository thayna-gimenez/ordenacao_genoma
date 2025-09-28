char** coletar_amostras_gatherv(char **local_strings, int local_count, 
                               int rank, int size, int n_amostras) {
    char **amostras = NULL;
    char **todas_amostras = NULL;
    
    // Cada processo prepara suas amostras
    amostras = malloc(n_amostras * sizeof(char*));
    for (int i = 0; i < n_amostras; i++) {
        amostras[i] = malloc(MAX_STR_LEN * sizeof(char));
        
        // Seleciona amostras igualmente espaçadas
        int index = (i + 1) * local_count / (n_amostras + 1);
        if (index >= local_count) index = local_count - 1;
        
        if (local_count > 0) {
            strncpy(amostras[i], local_strings[index], MAX_STR_LEN);
        } else {
            strcpy(amostras[i], "");  // Processo sem dados
        }
        amostras[i][MAX_STR_LEN - 1] = '\0';
    }
    
    // Prepara arrays para Gatherv
    int *recvcounts = NULL;
    int *displs = NULL;
    int total_amostras = size * n_amostras;
    
    if (rank == 0) {
        recvcounts = malloc(size * sizeof(int));
        displs = malloc(size * sizeof(int));
        todas_amostras = malloc(total_amostras * sizeof(char*));
        
        for (int i = 0; i < size; i++) {
            recvcounts[i] = n_amostras * MAX_STR_LEN;  // Cada processo envia n_amostras strings
            displs[i] = i * n_amostras * MAX_STR_LEN;  // Deslocamento no buffer
        }
        
        // Aloca buffer contíguo para receber todas as amostras
        for (int i = 0; i < total_amostras; i++) {
            todas_amostras[i] = malloc(MAX_STR_LEN * sizeof(char));
        }
    }
    
    // Cria buffer contíguo para envio
    char *send_buffer = malloc(n_amostras * MAX_STR_LEN * sizeof(char));
    for (int i = 0; i < n_amostras; i++) {
        strncpy(&send_buffer[i * MAX_STR_LEN], amostras[i], MAX_STR_LEN);
    }
    
    // Buffer de recepção (apenas root usa)
    char *recv_buffer = NULL;
    if (rank == 0) {
        recv_buffer = malloc(total_amostras * MAX_STR_LEN * sizeof(char));
    }
    
    // Coleta todas as amostras no root usando Gatherv
    MPI_Gatherv(send_buffer, n_amostras * MAX_STR_LEN, MPI_CHAR,
               recv_buffer, recvcounts, displs, MPI_CHAR,
               0, MPI_COMM_WORLD);
    
    // Root converte buffer de volta para array de strings
    if (rank == 0) {
        for (int proc = 0; proc < size; proc++) {
            for (int i = 0; i < n_amostras; i++) {
                int index_global = proc * n_amostras + i;
                strncpy(todas_amostras[index_global], 
                       &recv_buffer[displs[proc] + i * MAX_STR_LEN], 
                       MAX_STR_LEN);
                todas_amostras[index_global][MAX_STR_LEN - 1] = '\0';
            }
        }
    }
    
    // Limpeza
    for (int i = 0; i < n_amostras; i++) {
        free(amostras[i]);
    }
    free(amostras);
    free(send_buffer);
    
    if (rank == 0) {
        free(recvcounts);
        free(displs);
        free(recv_buffer);
    }
    
    return todas_amostras;  // Apenas root retorna algo útil, outros retornam NULL
}