#include "throughput.h"

int main(){

    uint32_t total_length = 0;
    struct seqfile_t RF;
    RF.seqid = 0;
    char reads_path[256] = READS_FILE_PATH;
    RF.file = fopen(strcat(reads_path,READS_FILE_NAME), "r");

    if(RF.file == NULL) { printf("\033[1m\033[31m[Main] Error: reads file not found\033[37m\n"); return -1; }
    
     	while (fgetc(RF.file) != EOF){
        getNextRead(&RF);
        total_length += RF.seqlen;
        printf("totalLength %lu \n", total_length);

        if(fgetc(RF.file) == EOF){
            break;
        } 
    } 

    printf("totalLength %lu \n", total_length);

return 0;

}