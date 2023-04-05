#include "mthread.h"
#include "map.h"

#ifdef MULTITHREADING

#ifndef MODE_COMPRESSED
    bp_t chunk[MAX_GENOME_CHUNK_SIZE];
#else
    bp32_t chunk[MAX_GENOME_CHUNK_SIZE/32];
#endif

pthread_mutex_t TF_mutex;           //Used by the the mapping function
pthread_mutex_t RF_mutex;           //Used here to fill the thread file


double filtering_time_vector[THREADSNUM];
double alignment_time_vector[THREADSNUM];

void mutex_group_create(){
    pthread_mutex_init(&TF_mutex, NULL);
    pthread_mutex_init(&RF_mutex, NULL);
}

void mutex_group_destroy(){
    pthread_mutex_destroy(&TF_mutex);
    pthread_mutex_destroy(&RF_mutex);
}

void * mapper_thread(void* in){

    int cpu;
    int node;
    struct mapper_ctx_t * thread_ctx = (struct mapper_ctx_t *)in;
    char thread_reads_path[256] = READS_FILE_PATH;

    #ifdef MODE_COMPRESSED

        createThreadReadFile_bp32(thread_ctx, thread_reads_path);

        pthread_mutex_lock(&RF_mutex);
        fillThreadReadFile_bp32(thread_ctx, thread_reads_path);
        pthread_mutex_unlock(&RF_mutex);

        #ifdef VERBOSE
            printf("[Main][Thread<%03u>][Cpu<%03u>] Running full genome mapping\n", thread_ctx->id, sched_getcpu());
        #endif

        MapReadsToGenome_bp32(thread_ctx->TF_global, &(thread_ctx->RF_local), NULL);

    #else
        createThreadReadFile(thread_ctx, thread_reads_path);
        printf("[Main][Thread<%03u>][Cpu<%03u>] Performing on %u reads\n", thread_ctx->id, sched_getcpu(), thread_ctx->reads_num);

        pthread_mutex_lock(&RF_mutex);
        fillThreadReadFile(thread_ctx, thread_reads_path);
        pthread_mutex_unlock(&RF_mutex);

        #ifdef VERBOSE
            printf("[Main][Thread<%03u>][Cpu<%03u>] Running full genome mapping\n", thread_ctx->id, sched_getcpu());
        #endif

        MapReadsToGenome(thread_ctx->TF_global, &(thread_ctx->RF_local), NULL);

        filtering_time_vector[thread_ctx->id] = filtering_time;
        alignment_time_vector[thread_ctx->id] = alignment_time;

        #ifdef VERBOSE
            printf("[Main][Thread<%03u>][Cpu<%03u>] Thread running\n", thread_ctx->id, sched_getcpu());
            printf("[Main][Thread<%03u>][Cpu<%03u>] Creating thread local read files\n", thread_ctx->id, sched_getcpu());
        #endif


    #endif

    fclose(thread_ctx->RF_local.file);
    //printf("[Main][Thread<%03u>] exiting ... \n",ctx->id);
        
}

#ifndef MODE_COMPRESSED

void createThreadReadFile(struct mapper_ctx_t * thread_ctx, char * thread_reads_path){

    /* Create all temporary files */

    strcat(thread_reads_path,"thread/tmp_000.fasta");

    thread_reads_path[strlen(thread_reads_path) - strlen("000.fasta") + 0] = thread_ctx->id/100+48;
    thread_reads_path[strlen(thread_reads_path) - strlen("000.fasta") + 1] = (thread_ctx->id%100)/10+48;
    thread_reads_path[strlen(thread_reads_path) - strlen("000.fasta") + 2] = (thread_ctx->id%10)+48;

    thread_ctx->RF_local.file = fopen(thread_reads_path, "w+");
    if(thread_ctx->RF_local.file == NULL) { printf("\033[1m\033[31m[Main] Error: thread[%u] read file was not created\033[37m\n", thread_ctx->id); exit; }
}

void fillThreadReadFile(struct mapper_ctx_t * thread_ctx, char * thread_reads_path){

    uint32_t chunk_size = MAX_GENOME_CHUNK_SIZE;
    uint32_t read_len = 0;
    //uint64_t chunk_size;
    uint32_t chunk_num;
    uint32_t size;

    /* Copy original read file into temporary thread files */

    for( uint32_t i = 0; i < thread_ctx->reads_num; i++ ){

		getNextRead(thread_ctx->RF_global);
        read_len = thread_ctx->RF_global->seqlen;

        //printf("1) read_len %u \n", read_len);

        fprintf(thread_ctx->RF_local.file, ">MICROMAP_MTHREAD %u length=%u\n", i+1, read_len);

       // printf("2) read_len %u \n", read_len);
        chunk_size = (read_len < chunk_size) ? read_len : chunk_size;
        chunk_num = ceil((double)  read_len / chunk_size);

        //if(i%thread_num == 0)printf("Writing a seq of len %05u, with a chunk of %05u and %01u chunks to file %03u-%016x\n", read_len, chunk_size, chunk_num, i%thread_num, ftell(thread_RF[i%thread_num].file));

        for( uint32_t j = 0; j < chunk_num; j++ ){
                    
            if( j == chunk_num-1 && (read_len % chunk_size) != 0 ){
                size = read_len % chunk_size;
            }else{
                size = chunk_size;
            } 

            getReadChunk(thread_ctx->RF_global, chunk, (j * chunk_size), size);
            
            for( uint32_t k = 0; k < size; k++ ){
                fputc(chunk[k], thread_ctx->RF_local.file);
            }

            fputc('\n', thread_ctx->RF_local.file);
        }
                
    }

    // reset thread file pointers

    fseek(thread_ctx->RF_local.file, 0, SEEK_SET);
}

#else

void createThreadReadFile_bp32(struct mapper_ctx_t * thread_ctx, char * thread_reads_path){

    /* Create all temporary files */

    strcat(thread_reads_path,"thread/tmp_000.bp32");

    thread_reads_path[strlen(thread_reads_path) - strlen("000.bp32") + 0] = thread_ctx->id/100+48;
    thread_reads_path[strlen(thread_reads_path) - strlen("000.bp32") + 1] = (thread_ctx->id%100)/10+48;
    thread_reads_path[strlen(thread_reads_path) - strlen("000.bp32") + 2] = (thread_ctx->id%10)+48;

    thread_ctx->RF_local.file = fopen(thread_reads_path, "wb+");
    if(thread_ctx->RF_local.file == NULL) { printf("\033[1m\033[31m[Main] Error: thread[%u] read file was not created\033[37m\n", thread_ctx->id); exit; }

    fwrite(&thread_ctx->reads_num, sizeof(uint64_t), 1, thread_ctx->RF_local.file);

}

void fillThreadReadFile_bp32(struct mapper_ctx_t * thread_ctx, char * thread_reads_path){

    uint64_t read_len = 0;
    uint64_t chunk_size;
    uint64_t chunk_num;
    uint32_t size;

    /* Copy original read file into temporary thread files */

    for( uint64_t i = 0; i < thread_ctx->reads_num; i++ ){

		getNextRead_bp32(thread_ctx->RF_global);
        read_len = (uint64_t)thread_ctx->RF_global->seqlen;
  
        fwrite(&read_len, sizeof(uint64_t), 1, thread_ctx->RF_local.file);
        chunk_size = (read_len < MAX_GENOME_CHUNK_SIZE) ? read_len : MAX_GENOME_CHUNK_SIZE;
        chunk_num = ceil(read_len / chunk_size);

        //if(i%thread_num == 0)printf("Writing a seq of len %05u, with a chunk of %05u and %01u chunks to file %03u-%016x\n", read_len, chunk_size, chunk_num, i%thread_num, ftell(thread_RF[i%thread_num].file));

        for( uint32_t j = 0; j < chunk_num; j++ ){
                    
            getReadChunk_bp32(thread_ctx->RF_global, chunk, j * chunk_size, chunk_size);
            size = (chunk_size%32 == 0) ? chunk_size/32 : chunk_size/32 + 1;
            for( uint32_t k = 0; k < size; k++ ){
                fwrite(&chunk[k], sizeof(uint64_t), 1, thread_ctx->RF_local.file);
            }
        }
                
    }

    // reset thread file pointers

    fseek(thread_ctx->RF_local.file, 0, SEEK_SET);
}

/*void splitReadsFile_bp32(struct seqfile_t * RF, struct seqfile_t * TF, char * thread_reads_path, const uint16_t thread_num, struct mapper_ctx_t * thread_ctx){
        
    uint64_t reads_num;
    uint64_t read_len = 0;
    uint64_t thread_reads_num = 0;
    uint64_t chunk_size;
    uint64_t chunk_num;

    // First, create all temporary files 

    fread(&reads_num, sizeof(uint64_t), 1, RF->file);

    thread_reads_num = reads_num / thread_num;

    for( uint16_t i = 0; i < thread_num; i++ ){
                
        thread_reads_path[strlen(thread_reads_path) - strlen("000.bp32") + 0] = i/100+48;
        thread_reads_path[strlen(thread_reads_path) - strlen("000.bp32") + 1] = (i%100)/10+48;
        thread_reads_path[strlen(thread_reads_path) - strlen("000.bp32") + 2] = (i%10)+48;
        thread_ctx[i].RF.file = fopen(thread_reads_path, "wb+");
        if(thread_ctx[i].RF.file == NULL) { printf("\033[1m\033[31m[Main] Error: thread[%u] read file was not created\033[37m\n", i); return -1; }

        thread_ctx[i].TF = *TF;
        fwrite(&thread_reads_num, sizeof(uint64_t), 1, thread_ctx[i].RF.file);
    }

    // Copy original read file into temporary thread files 

    uint32_t size;

    for( uint64_t i = 0; i < reads_num; i++ ){

		getNextRead_bp32(&RF);
        read_len = (uint64_t)RF.seqlen;
  
        fwrite(&read_len, sizeof(uint64_t), 1, thread_ctx[i%thread_num].RF.file);
        chunk_size = (RF.seqlen < MAX_GENOME_CHUNK_SIZE) ? RF.seqlen : MAX_GENOME_CHUNK_SIZE;
        chunk_num = ceil(RF.seqlen / chunk_size);

        //if(i%thread_num == 0)printf("Writing a seq of len %05u, with a chunk of %05u and %01u chunks to file %03u-%016x\n", read_len, chunk_size, chunk_num, i%thread_num, ftell(thread_RF[i%thread_num].file));

        for( uint32_t j = 0; j < chunk_num; j++ ){
                    
            getReadChunk_bp32(&RF, chunk, j * chunk_size, chunk_size);
            size = (chunk_size%32 == 0) ? chunk_size/32 : chunk_size/32 + 1;
            for( uint32_t k = 0; k < size; k++ ){
                fwrite(&chunk[k], sizeof(uint64_t), 1, thread_ctx[i%thread_num].RF.file);
            }
        }
                
    }

    // reset thread file pointers

    for( uint16_t i = 0; i < thread_num; i++ ) fseek(thread_ctx[i%thread_num].RF.file, 0, SEEK_SET);
}*/

#endif

#endif