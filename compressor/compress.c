#include "compress.h"

void compressReferenceGenome(FILE * genFile, FILE * cmpGenFile){

    char tmp;
    uint64_t cnt = 0;
    uint64_t bp[32];
    uint64_t bp32 = 0;
    uint64_t genome_start;
    uint64_t sequence_size = 0;



    while( fgetc(genFile) != '\n' ); 	            //Skip first line: it only contains metadata
     genome_start = ftell(genFile);                  //Record genome starting point
    while( fgetc(genFile)!= EOF ) sequence_size++;    //Compute Genome Size
    fseek(genFile,genome_start,SEEK_SET);

    fwrite(&sequence_size, sizeof(uint64_t), 1, cmpGenFile);     //Write Genome size on file
    tmp = fgetc(genFile);
    while(tmp != EOF){

        switch(tmp){
			case 'A': bp[cnt] = 0x0; break;
			case 'T': bp[cnt] = 0x3; break;
			case 'C': bp[cnt] = 0x1; break;
			case 'G': bp[cnt] = 0x2; break;
		}

        cnt++;
        tmp = fgetc(genFile);

        if(cnt == 32){
            for(int i = 0; i < 32; i++) bp32 |= (bp[i] << 2*i); 
            fwrite(&bp32, sizeof(uint64_t), 1, cmpGenFile);
            bp32 = 0;
            cnt = 0;
        }

    } 
}

void compressReads(FILE * readFile, FILE * cmpReadFile){

    char *line_buf = NULL;
	char *subline_buf = NULL;
	size_t line_size;
    uint64_t readNum = 0;

    char tmp;
    uint64_t cnt = 0;
    uint64_t bp[32];
    uint64_t bp32 = 0;
    uint64_t sequence_size = 0;

    while(!feof(readFile)){
        if(fgetc(readFile) == '\n') readNum++;
    }

    readNum = readNum/2;
    fwrite(&readNum, sizeof(uint64_t), 1, cmpReadFile);
    fseek(readFile, 0, SEEK_SET);

    while(!feof(readFile)){

        getline(&line_buf, &line_size, readFile);
        subline_buf = strstr(line_buf, "length=");

        if(subline_buf != NULL){
            sequence_size = atoi(subline_buf + sizeof("length=") - 1);
            fwrite(&sequence_size, sizeof(uint64_t), 1, cmpReadFile);
        }

        for(int i = 0; i < sequence_size; i++){

            tmp = fgetc(readFile);

            switch(tmp){
                case 'A': bp[cnt] = 0x0; break;
                case 'T': bp[cnt] = 0x3; break;
                case 'C': bp[cnt] = 0x1; break;
                case 'G': bp[cnt] = 0x2; break;
		    }

            cnt++;
            
            if(cnt == 32){
                for(int i = 0; i < 32; i++) bp32 |= (bp[i] << 2*i); 
                fwrite(&bp32, sizeof(uint64_t), 1, cmpReadFile);
                
                bp32 = 0;
                cnt = 0;
            }
        }

        if(cnt != 0){
            
            for(int i = 0; i < cnt; i++) bp32 |= (bp[i] << 2*i); 
            fwrite(&bp32, sizeof(uint64_t), 1, cmpReadFile);

        }
        
        fgetc(readFile); //Read the 'n'

        bp32 = 0;
        cnt = 0;

    }

}