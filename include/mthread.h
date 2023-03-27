#ifndef MTHREAD
#define MTHREAD

#include "types.h"
#include "param.h"

    #ifdef MULTITHREADING

        void mutex_group_create();
        void mutex_group_destroy();
        void * mapper_thread(void* in);

        #ifndef MODE_COMPRESSED

        void createThreadReadFile(struct mapper_ctx_t * thread_ctx, char * thread_reads_path);
        void fillThreadReadFile(struct mapper_ctx_t * thread_ctx, char * thread_reads_path);

        #else

        void createThreadReadFile_bp32(struct mapper_ctx_t * thread_ctx, char * thread_reads_path);
        void fillThreadReadFile_bp32(struct mapper_ctx_t * thread_ctx, char * thread_reads_path);

        #endif

    #endif
    
#endif