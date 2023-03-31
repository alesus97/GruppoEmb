
#include "compress.h"

int main( int argc, char *argv[] ){

    /* Reference Genome and Reads files structures*/

    FILE * readFile;
    FILE * genFile;
    FILE * cmpReadFile;
    FILE * cmpGenFile;

    uint8_t flag;

    char genomePath[1024];
    char readsPath[1024];

     if( argc == 1 ){
        printf("\nNo arguments provided\n\t");
        printf("Use -r 'path/to/reads' to load the reads fasta file\n\t");
        printf("Use -g 'path/to/genome' to load the genome fasta file\n\t");
        printf("Use -rg 'path/to/reads' 'path/to/genome' to load both reads and genome fasta files\n\n");
    }  
    else if( argc == 3 | argc == 4){
        if(!strcmp("-r",argv[1])){
            strcpy(readsPath,argv[2]);
            flag = 1;
        }
        else if(!strcmp("-g",argv[1])){
            strcpy(genomePath,argv[2]);
            flag = 2;
        }
        else if(!strcmp("-rg",argv[1])){
            strcpy(readsPath,argv[2]);
            strcpy(genomePath,argv[3]);
            flag = 3;
        }
        else{
            printf("\nWrong format\n\t");
            printf("Use -r 'path/to/reads' to load the reads fasta file\n\t");
            printf("Use -g 'path/to/genome' to load the genome fasta file\n\t");
            printf("Use -rg 'path/to/reads' 'path/to/genome' to load both reads and genome fasta files\n\n");
            return -1;
        }

        // Build bp32 file for the Genome

        if(flag == 2 | flag == 3){

           // genFile = readFile = fopen("../data/genomes/salmonellaEntericaShrtN.fasta", "r");
            genFile = fopen(genomePath, "r");
            if(genFile == NULL){ printf("[Main] Error : genome file to be compressed (.FASTA) not found\n"); return -1; }
            cmpGenFile = fopen("compressedGenome.bp32","wb+");

            printf("Compressing Reference Genome :: FASTA 2 BP32\t");
            compressReferenceGenome(genFile, cmpGenFile);
            printf("Done\n");

            fclose(genFile);
            fclose(cmpGenFile); 

        }

        // Build bp32 file for reads

         if(flag == 1 | flag == 3){

            //readFile  = fopen("../data/reads/salmonellaEntericaReadsShrt.fasta", "r");
            readFile = fopen(readsPath, "r");
            if(readFile == NULL){ printf("[Main] Error : reads file to be compressed (.FASTA) not found\n"); return -1; }
            cmpReadFile = fopen("compressedReads.bp32","wb+");

            printf("Compressing Read file :: FASTA 2 BP32\t");
            compressReads(readFile, cmpReadFile);
            printf("Done\n");

            fclose(readFile);  
            fclose(cmpReadFile);

        } 

    }
     else{
        printf("\nWrong format\n\t");
        printf("Use -r 'path/to/reads' to load the reads fasta file\n\t");
        printf("Use -g 'path/to/genome' to load the genome fasta file\n\t");
        printf("Use -rg 'path/to/reads' 'path/to/genome' to load both reads and genome fasta files\n\n");
    } 

    return 0;
}
