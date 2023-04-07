#define _GNU_SOURCE
#include "types.h"
#include "param.h"
#include "time.h"

#include "seeding.h"
#include "map.h"

#ifdef MULTITHREADING

    #include "mthread.h"

    extern double filtering_time_vector[THREADSNUM];
    extern double alignment_time_vector[THREADSNUM];
    extern uint32_t total_reads_vector[THREADSNUM];

#endif

/* Time measuring global variables */

#ifdef TAKE_TIME

    struct timespec start, end;
    double time_elapsed;

    extern double filtering_time, alignment_time;
    double  seeding_time;

#endif

int main(){

    /* Reference Genome and Reads files structures*/

    char target_path[256] = TARGET_FILE_PATH;
    char reads_path[256] = READS_FILE_PATH;

	struct seqfile_t RF;
	struct seqfile_t TF;

    RF.seqid = 0;

    /* Printf Introduction */

    #ifdef VERBOSE
        printLogo();
        printf("[Main] Running version 0.6\n");

        #ifdef MODE_COMPRESSED
        printf("[Main] Compressed format enabled\n");
        #else
            printf("[Main] Compressed format disabled\n");
        #endif
    #endif

    /* Define the Index Structure: a hashband (look at types.h for more info) */

    TF.file = fopen(strcat(target_path,TARGET_FILE_NAME), "r");
    if(TF.file == NULL){ printf("\033[1m\033[31m[Main] Error: reference genome file not found\033[37m\n"); return -1; }
	#ifdef MODE_COMPRESSED
        if(strstr(target_path, ".bp32") == NULL) { printf("\033[1m\033[31m[Main] Error: Wrong format, file .bp32 not found\033[37m\n"); return -1; }
    #else
        if(strstr(target_path, ".fasta") == NULL) { printf("\033[1m\033[31m[Main] Error: Wrong format, file .fasta not found\033[37m\n"); return -1; }
    #endif

    /* Create the Index Structure; this part can be skipped if the index has already been built */

    #ifdef BUILD_GENOME_INDEX
        #ifdef TAKE_TIME
            clock_gettime(CLOCK_MONOTONIC, &start);
        #endif

        #ifdef VERBOSE
            printf("[Main] Creating genome Index\n");
        #endif

        #ifdef MODE_COMPRESSED
            BuildIndex_bp32(&TF);
        #else
            BuildIndex(&TF);
        #endif
        
        #ifdef TAKE_TIME 
            clock_gettime(CLOCK_MONOTONIC, &end);
            time_elapsed = (end.tv_sec - start.tv_sec);
            time_elapsed += (end.tv_nsec - start.tv_nsec) / 1000000000.0;
            #ifdef VERBOSE
                printf("[Main] Index built in %08.8f s\n", time_elapsed);
            #endif
        #endif
    #endif

    RF.file = fopen(strcat(reads_path,READS_FILE_NAME), "r");
    if(RF.file == NULL) { printf("\033[1m\033[31m[Main] Error: reads file not found\033[37m\n"); return -1; }

    #ifdef MODE_COMPRESSED
        if(strstr(target_path, ".bp32") == NULL) { printf("\033[1m\033[31m[Main] Error: Wrong format, file .bp32 not found\033[37m\n"); return -1; }
    #else
        if(strstr(target_path, ".fasta") == NULL) { printf("\033[1m\033[31m[Main] Error: Wrong format, file .fasta not found\033[37m\n"); return -1; }
    #endif

    #ifdef TAKE_TIME
        clock_gettime(CLOCK_MONOTONIC, &start);
    #endif

    #ifndef MULTITHREADING
        
        #ifdef VERBOSE
            printf("\n[Main] Full-Genome Mapping begins ... \n");
        #endif

        #ifdef MODE_COMPRESSED
            MapReadsToGenome_bp32(&TF, &RF, NULL);
        #else

        #ifdef TAKE_TIME
            uint32_t total_read_length = totalReadLength(&RF);
        #endif

        MapReadsToGenome(&TF, &RF, NULL, NULL);
            
        #endif



    #else
        
        uint64_t reads_num = 0;
        uint16_t thread_num;
        struct mapper_ctx_t * thread_ctx;  

        thread_num = (get_nprocs() <= THREADSNUM) ? get_nprocs() : THREADSNUM;
        thread_ctx = (struct mapper_ctx_t *) malloc(sizeof(struct mapper_ctx_t)*thread_num);
        
        /* Creating a read file for each thread */

        #ifdef VERBOSE
            printf("[Main] multithreading mode activated - %02u Threads\n", thread_num);
        #endif

        #ifdef MODE_COMPRESSED

            if(strstr(target_path, ".bp32") == NULL) { printf("\033[1m\033[31m[Main] Error: Wrong format, file .bp32 not found\033[37m\n"); return -1; }

            //Get the number of total reads to know how many reads each thread is going to handle
            fread(&reads_num, sizeof(uint64_t), 1, RF.file);

        #else

            if(strstr(target_path, ".fasta") == NULL) { printf("\033[1m\033[31m[Main] Error: Wrong format, file .fasta not found\033[37m\n"); return -1; }

            //Get the number of total reads to know how many reads each thread is going to handle
            while(!feof(RF.file)){
                if(fgetc(RF.file) == '\n') reads_num++;
            }
            reads_num = reads_num >> 1; //for each read 2 '\n' will be counted, therefore we must divide by two
            fseek(RF.file, 0, SEEK_SET);

        #endif

        /* Init all required mutexex for the threads to share the Index and the global Read file */
        mutex_group_create();
        
       // int[thread_num];

        for (int i = 0; i < reads_num%thread_num; i++) {

            thread_ctx[i].id = i;
            thread_ctx[i].reads_num = reads_num/thread_num+1;
            thread_ctx[i].TF_global = &TF;
            thread_ctx[i].RF_global = &RF;
            if(pthread_create(&(thread_ctx[i].thread_handler), NULL, mapper_thread, &thread_ctx[i]) != 0) {printf("\033[1m\033[31m[Main] Error: Thread number <%d> failed\033[37m\n", i); return -1; }; 
        }

        for (int i = reads_num%thread_num ; i < thread_num  ; i++){
            thread_ctx[i].id = i;
            thread_ctx[i].reads_num = reads_num/thread_num;
            thread_ctx[i].TF_global = &TF;
            thread_ctx[i].RF_global = &RF;
            if(pthread_create(&(thread_ctx[i].thread_handler), NULL, mapper_thread, &thread_ctx[i]) != 0) {printf("\033[1m\033[31m[Main] Error: Thread number <%d> failed\033[37m\n", i); return -1; }; 
        }


        
        for(int i = 0; i < thread_num; i++){
            if(pthread_join(thread_ctx[i].thread_handler, NULL)!= 0) {printf("\033[1m\033[31m[Main] Error: Thread number <%d> failed\033[37m\n", i); return -1; };
        }

        mutex_group_destroy();
        free(thread_ctx);

    #ifdef TAKE_TIME
        double sum_alignment = 0;
        double sum_filtering = 0;
     
        for(uint32_t i=0; i<thread_num; i++){
            sum_alignment += alignment_time_vector[i];
            sum_filtering += filtering_time_vector[i];
            
           // printf("throughput per thread %f \n", throughput_per_thread[i]);

        }

        filtering_time = sum_filtering/thread_num;
        alignment_time = sum_alignment/thread_num;
    #endif
        
    #endif

    #ifdef TAKE_TIME
        clock_gettime(CLOCK_MONOTONIC, &end);
        time_elapsed = (end.tv_sec - start.tv_sec);
        time_elapsed += (end.tv_nsec - start.tv_nsec) / 1000000000.0;
        
        #ifdef VERBOSE
            printf("[Main] Full-Genome Mapping completed in %08.8f s; %08u reads mapped\n\n", time_elapsed, RF.seqid/* +1 */ );
        #endif

        #ifdef MULTITHREADING

            double sum_throughput = 0;
            double throughput_per_thread[thread_num];
            double throughput;
        
            for(uint32_t i=0; i<thread_num; i++){
            throughput_per_thread[i] = total_reads_vector[i]/(pow(10,3)*time_elapsed);
            sum_throughput += throughput_per_thread[i];
            
            }

            throughput = sum_throughput;
        #else 

            double throughput = total_read_length/(pow(10,3)*time_elapsed);

        #endif

        seeding_time = time_elapsed - filtering_time - alignment_time;

        #ifdef VERBOSE

            printf("[Main] Throughput: %08.8f reads/s\n", throughput);
            printf("[Main] Seeding time %08.8f s\n", seeding_time);
            printf("[ChunkElaboration] Filtering time %08.8f s\n", filtering_time);
            printf("[ChunkElaboration] Alignment time %08.8f s\n\n", alignment_time);

        #endif  
    #endif
    



    freeHashTable(&TF.index);
    fclose(RF.file);
    fclose(TF.file);

   

}
