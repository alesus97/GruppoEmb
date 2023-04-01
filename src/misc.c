#include "misc.h"

/* hash-band related functions */

void printHashmizer( const struct hashmizer_t * mini ){

	printf(" -> { # 0x%016lx ,", mini->hash); printf(" %03u ,", mini->len ); printf(" %019lu ", mini->start );
	if( mini->strand == 0 ) printf(", + }");
	else printf(", - }");
}

void printHashTable( const struct hashband_t * ht, uint8_t param ){

	
	for( int i = 0; i < ht->size; i++ ){
		if( param == HASHTABLE_PRINT_FULL ){
			printf(" %04u ||", i); 
			printList( ht, i );
		}
		else if( param == HASHTABLE_PRINT_DENSE ){
			if( getListLength( ht, i ) ){
				printf(" %04u ||", i); 
				printList( ht, i );
			}
		}
	}
}

int getListLength( const struct hashband_t * ht, const uint32_t index ){

	int length = 0;
	struct hashnode_t *current;
	
	for(current = ht->entry[index]; current != NULL; current = current->next) {
     		length++;
	}

	return --length;
}

void printList( const struct hashband_t * ht, const uint32_t index ) {

	struct hashnode_t *ptr;
   	int len = getListLength( ht, index );
	ptr = ht->entry[index];
   	
	while(ptr != NULL) {
	
		if( ptr->next != NULL || !len ){
			if(len) printf("\033[0;32m");
			printf(" -> { # 0x%016lx ", ptr->hashmizer.hash); printf(", %019lu ", ptr->hashmizer.start );
			if( ptr->hashmizer.strand == 0 ) printf(", + }\033[0;37m");
			else printf(", - }\033[0;37m");
		}
		ptr = ptr->next;

	}

	printf("\n");
}

/* sequences related functions */

uint64_t min64( uint64_t a, uint64_t b){

	uint64_t min;
	
	min = ( a > b ) ? b : a ;
	
	return min;
}

#ifndef MODE_COMPRESSED

void printSequence( const char * seq, uint64_t start, uint64_t end, uint8_t newline ){

	printf("> ");
	for( uint32_t i = start; i < end; i++ ){
		switch(seq[i]){
			case 'A': printf("\033[0;33m"); break;
			case 'T': printf("\033[0;32m"); break;
			case 'C': printf("\033[0;35m"); break;
			case 'G': printf("\033[0;36m"); break;
			default : printf("\033[0;37m"); break;
			
		}
		printf("%c",seq[i]);
	}
	if( newline ) printf("\033[0;37m\n");
	else printf("\033[0;37m ");
}

void reverseComplement( char * A, char * B, const uint32_t L ){

	for( int i = 0; i < L; i++ ){
	
		switch( A[L-i-1] ){
			case 'A': B[i] = 'T'; break;
			case 'T': B[i] = 'A'; break;
			case 'C': B[i] = 'G'; break;
			case 'G': B[i] = 'C'; break;
			default : B[i] = 'N'; break;
		}
	}
}


void getNextRead(struct seqfile_t * RF){

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

	RF->seqstart = ftell(RF->file);
	
}

void getReadChunk(struct seqfile_t * RF, bp_t * chunk, uint32_t start, uint32_t size){

	uint32_t i = 0;
	char tmp;

	fseek(RF->file, RF->seqstart + start, SEEK_SET);
	fread(chunk, sizeof(uint8_t), size, RF->file);
	fseek(RF->file, RF->seqstart, SEEK_SET);

}


#else

void printSequence_bp32( const bp32_t * seq, uint64_t start, uint64_t end, uint8_t newline ){

	uint64_t start_in = start % 32;
	uint64_t end_in = end % 32;
	uint8_t tmp;
	uint8_t j;

	printf("> ");
	for( uint32_t i = start/32; i <= end/32; i++ ){

		if(i == start/32) 	start_in = start % 32;
		else				start_in = 0;
		if(i == end/32)		end_in = end % 32;
		else				end_in = 32;

		for( uint64_t j = start_in; j < end_in; j++ ){
			tmp = (seq[i] >> 2*(j)) & 0x3;
			switch(tmp){
				case 0: printf("\033[0;33mA"); break;
				case 3: printf("\033[0;32mT"); break;
				case 1: printf("\033[0;35mC"); break;
				case 2: printf("\033[0;36mG"); break;
				default : printf("\033[0;37mN"); break;
			}
		}
	}
	if( newline ) printf("\033[0;37m\n");
	else printf("\033[0;37m ");
}


void reverseComplement_bp32( bp32_t * seq, const uint32_t seq_size ){

	uint64_t word_num;
	bp32_t seqrev_a = 0;
	bp32_t seqrev_b = 0;
	bp32_t word_a;
	bp32_t word_b;
	uint64_t mask_b;
	uint8_t l = seq_size % 32;
	uint8_t h = 32 - l;

	word_num = ceil((float)seq_size/32);

	for( int j = 0; j < word_num/2; j++ ){

		//Select Reverse words
		word_a = seq[j];
		word_b = seq[word_num - j - 1];

		//Reverse
		seqrev_a = 0;
		seqrev_b = 0;

		for( int k = 0; k < 32; k++ ){
			seqrev_a |= (((word_a >> 2*k) & 0x3) << (62-2*k));
			seqrev_b |= (((word_b >> 2*k) & 0x3) << (62-2*k));
		}

		//Complement
		mask_b = pow(2,2*l) - 1;
		seqrev_a ^= 0xffffffffffffffff ;
		seqrev_b = (!j && l) ? (mask_b)^(seqrev_b >> 2*h) : 0xffffffffffffffff ^seqrev_b;

		//Reverse Assign
		seq[j] = seqrev_b;
		seq[word_num - j - 1] = seqrev_a;
	}


	//If word_num is odd reverse and complement word number word_num/2+1

	if(word_num % 2 != 0){
		seqrev_a = 0;
		for( int k = 0; k < 32; k++ ) seqrev_a |= (((seq[word_num/2] >> 2*k) & 0x3) << (62-2*k));
		seq[word_num/2] = 0xffffffffffffffff ^ seqrev_a;
	}


	//If seq_size is not a multiple of 32, we must perform a final shift

	if(l){
		seq[0] |= (seq[1] << 2*l);
		for(int i = 1; i < word_num-1; i++){
			seq[i] = (seq[i] >> 2*h) | (seq[i+1] << 2*l);
		}
		seq[word_num-1] = (seq[word_num - 1] >> 2*h);
	}


}

void getNextRead_bp32(struct seqfile_t * RF){

	if(RF->seqlen != 0){
		if(RF->seqlen % 32) fseek(RF->file, RF->seqstart + (RF->seqlen - RF->seqlen%32)/4 + 8, SEEK_SET);
		else fseek(RF->file, RF->seqstart + (RF->seqlen/4), SEEK_SET);
	}

	fread(&RF->seqlen, sizeof(uint64_t), 1, RF->file); 
	RF->seqid++;
	RF->seqstart = ftell(RF->file);
}

void getReadChunk_bp32(struct seqfile_t * RF, bp32_t * chunk, uint32_t start, uint32_t size){

	uint64_t start_in;

	start_in = start % 32;
	//RF->seqstart will always be aligned to 64 bits. So we only need to know if start is % 32 (start is in bps)

	if(!start_in){
		fseek(RF->file, RF->seqstart + start/4, SEEK_SET);
		fread(chunk, sizeof(uint64_t), ceil(size/32)+1, RF->file);
	}else{
		fseek(RF->file, RF->seqstart + (start/32)*8, SEEK_SET);

		fread(&chunk[0], sizeof(uint64_t), 1, RF->file);
		chunk[0] = chunk[0] >> 2*start_in;

		for(uint64_t i = 1; i < ceil(size/32)+1; i++){
			fread(&chunk[i], sizeof(uint64_t), 1, RF->file);
			chunk[i-1] |= (chunk[i] << 2*(32-start_in));
			chunk[i] = chunk[i] >> 2*start_in;
		}

	}

	if(((start + size) % 32)){																//if the end is not 32 bit aligned, clean the remaining bits
		chunk[size/32] &= (uint64_t)(pow(2,2*((start + size) % 32))-1);		
	}
	
	fseek(RF->file, RF->seqstart, SEEK_SET);
}

#endif

void printLogo(){

	/***
	 *            _                                   
	 *      /\/\ (_) ___ _ __ ___   /\/\   __ _ _ __  
	 *     /    \| |/ __| '__/ _ \ /    \ / _` | '_ \ 
	 *    / /\/\ | | (__| | | (_) / /\/\ | (_| | |_) |
	 *    \/    \|_|\___|_|  \___/\/    \/\__,_| .__/ 
	 *                                         |_|    
	 */

	printf("\n");
	printf("      /\\/\\ (_) ___ _ __ ___   /\\/\\   __ _ _ __  \n");
	printf("     /    \\| |/ __| '__/ _ \\ /    \\ / _` | '_ \\ \n");
	printf("    / /\\/\\ | | (__| | | (_) / /\\/\\ | (_| | |_) |\n");
	printf("    \\/    \\|_|\\___|_|  \\___/\\/    \\/\\__,_| .__/ \n");
	printf("                                         |_|    \n");
	printf("\n");

}