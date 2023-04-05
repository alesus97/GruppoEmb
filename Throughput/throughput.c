#include "throughput.h"

void getNextRead(struct seqfile_t * RF){
	{
		char *line_buf = NULL;
		char *subline_buf = NULL;
		size_t line_size;

		
		fseek(RF->file, RF->seqstart + RF->seqlen + 1, SEEK_SET);

		getline(&line_buf, &line_size, RF->file);

		subline_buf = strstr(line_buf, "length=");
		if (subline_buf != NULL)
		{
			RF->seqlen = atoi(subline_buf + sizeof("length=") - 1);
			RF->seqid++;
		}
		printf("Lunghezza read : %d \n", RF->seqlen);
		RF->seqstart = ftell(RF->file);
	}
}