#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define READS_FILE_NAME         "salmonellaEntericaReadsShrt.fasta"
#define READS_FILE_PATH  	"../data/reads/"

struct seqfile_t{

	FILE * file;
	uint32_t seqid;
	uint32_t seqlen;
	uint32_t seqstart;
};

void getNextRead(struct seqfile_t * RF);
