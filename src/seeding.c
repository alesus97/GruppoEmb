#include "seeding.h"

#ifndef MODE_COMPRESSED

void BuildIndex( struct seqfile_t * TF ){

    // Phase 1 variables: Read the FASTA file and take genome metadata

    uint64_t genome_start = 0;
    uint64_t genome_size = 0;
	uint64_t genome_end;
	uint64_t genome_line_size;

    // Phase 2 variables: Read genome chunks according to MAX_GENOME_CHUNK_SIZE

    bp_t * ref_seqA;
    bp_t * ref_seqB;
    uint64_t seq_size;
    uint64_t last_seq_size = 0;
    uint64_t genome_start_ptr;
    uint64_t genome_end_ptr;
    char tmp;

    // Phase 3 variables: Compute minimizers into the chunks

    uint8_t w,k;
    uint64_t fin_seq_size;
    struct hashmizer_t mini;
    struct hashmizer_t mini_cpy;

    uint64_t chunkpos = 0;

    /*  
        ########################################################################  
        ##### 1 - Read the reference genome file and get the genome length #####
        ######################################################################## 
    */

    #ifdef VERBOSE
        printf("[BuildIndex] Reading ref genome size ... ");
    #endif
        
    while( fgetc(TF->file) != '\n' ); 	//Skip first line: it only contains metadata
    genome_start = ftell(TF->file);

    while( !feof(TF->file) ) if( fgetc(TF->file) != '\n' ) genome_size++;
    genome_end = ftell(TF->file)-1;	//Last base pair is before the eof

    fseek(TF->file, genome_start, SEEK_SET);
    TF->seqlen = genome_size;
    TF->seqstart = genome_start;

    #ifdef VERBOSE
        printf(" %010lu bps\n", genome_size);
    #endif

    /*  
        #######################################################################################  
        ##### 2 - Read the reference genome chunks by chunks (reverse complement as well) #####
        #######################################################################################
    */

    #ifdef VERBOSE
        printf("[BuildIndex] Initializing the Index\t");
    #endif

    genome_start_ptr = genome_start;
    genome_end_ptr = genome_end + 1;

    seq_size = ( genome_size <= MAX_GENOME_CHUNK_SIZE ) ? genome_size : MAX_GENOME_CHUNK_SIZE ;

    ref_seqA = (bp_t *)malloc(sizeof(bp_t)*seq_size);
    ref_seqB = (bp_t *)malloc(sizeof(bp_t)*seq_size);

    initHashTable(&(TF->index));

    #ifdef VERBOSE
        printf("Done\n[BuildIndex] Reading genome chunks of size %010lu\n", seq_size);
    #endif

    // For each chunk to be read

    for( int i = 0; i < ceil((double)(genome_size)/seq_size); i++ ){

        #ifdef VERBOSE
            printf("\t[BuildIndex] Working on chunks %04d/%.0f\t",i,ceil((double)(genome_size)/seq_size));
        #endif

        // read positive strand chunk and save it into SeqA

        if( genome_start_ptr <= genome_end ) fseek(TF->file, genome_start_ptr, SEEK_SET);

        for( int j = 0; j < seq_size; j++ ){
            tmp = fgetc(TF->file);
            if( tmp == '\0' || i*seq_size + j >= genome_size ) { last_seq_size = j - 1; break; }
            //if( tmp == '\n' ) j = j - 1;			
            else ref_seqA[j] = tmp;			        
        }

        #ifdef VERBOSE
            printf("\t\tChunk Size %08lu\tGSptr = %010lu", fin_seq_size, genome_start_ptr - genome_start);
        #endif

        genome_start_ptr += ( seq_size );  

        // read negative strand chunk and save it into SeqB
        
        if( genome_end_ptr >= seq_size && genome_end_ptr >= genome_start ) genome_end_ptr -= ( seq_size );
        else genome_end_ptr = genome_start;    
        fseek(TF->file, genome_end_ptr, SEEK_SET);

        #ifdef VERBOSE
            printf("\tGEptr = %010lu\t\t", genome_end_ptr - genome_start);
        #endif

        for( int j = 0; j < seq_size; j++ ){
            tmp = fgetc(TF->file);
            if( last_seq_size && j == last_seq_size ) break;			
            else ref_seqB[j] = tmp;			        
        }

        /*  
            #####################################################################  
            ##### 3 - For each couple of chunks (+ & -) find the minimizers #####
            #####################################################################
        */

        
        fin_seq_size = ( last_seq_size ) ? last_seq_size : seq_size; 
        w = MINI_W;
        k = MINI_K;

        for( int j = 0; j < fin_seq_size - w; j++ ){

            findMinimizer(ref_seqA + j, ref_seqB + fin_seq_size-w-j, j, fin_seq_size-w-j, MAX_64BIT, k, w, &mini);

            //Offset minimizer position into the full genome
            if( mini.strand == STRAND_MINUS ) mini.start = genome_size - i*seq_size - mini.start;
            else mini.start = i*seq_size + mini.start;

             //If the minimizer is not a duplicate add it to the index
            if( mini_cpy.hash != mini.hash && mini.start != mini_cpy.start ) {
                mini_cpy = mini;
                addToHashTable( &(TF->index), &mini );
            }

        }
        chunkpos += fin_seq_size;
        #ifdef VERBOSE
            printf("Done \n");
        #endif
        
    }

    #ifdef VERBOSE
        printf("[BuildIndex] Index succesfully built\n");
    #endif

    // Free resources
    free(ref_seqA);
    free(ref_seqB);

}

void findMinimizer( const bp_t * seqA, const bp_t * seqB, const uint64_t seqOffsetA, const uint64_t seqOffsetB, const uint64_t max, const uint8_t k, const uint8_t w, struct hashmizer_t * mini ){

    // For a detailed explanation of this algorithm check minimap paper.
    
    uint64_t m;
    uint64_t u,v;
	bp_t * kmerA;
    bp_t * kmerB;
	bp_t * rev_kmer;
	
	kmerA = (bp_t *)malloc(k*sizeof(bp_t));
	rev_kmer = (bp_t *)malloc(k*sizeof(bp_t));
	kmerB = (bp_t *)malloc(k*sizeof(bp_t));

    m = max;

	for( int j = 0; j < w-k+1; j++ ){

		strncpy(kmerA, seqA + j, k);
		strncpy(kmerB, seqB + w - k - j, k);
		reverseComplement( kmerB, rev_kmer, k);

		u = hashSequence(kmerA, k);
		v = hashSequence(rev_kmer,k);
		if( u != v ) m = min64( m, min64(u,v) );
    }

	for( int j = 0; j < w-k+1; j++ ){
		
		strncpy(kmerA, seqA + j, k);
		strncpy(kmerB, seqB + w - k - j, k);
		reverseComplement( kmerB, rev_kmer, k);
			
		u = hashSequence(kmerA, k);
		v = hashSequence(rev_kmer,k);
			
		if( u < v && u == m ) {
			
			mini->hash = u;
			mini->len = k;
			mini->start = seqOffsetA+j;
			mini->strand = STRAND_PLUS;
		}

		else if( u > v && v == m ) {
			
			mini->hash = u;
			mini->len = k;
			mini->start = seqOffsetB + w - k - j;
			mini->strand = STRAND_MINUS;
		}
	}

    free(kmerA);
	free(kmerB);
	free(rev_kmer);
}

#else

void BuildIndex_bp32( struct seqfile_t * TF ){

    // Phase 1 variables: Read the FASTA file and take genome metadata

    uint64_t genome_start = 0;
    uint64_t genome_size = 0;
	uint64_t genome_end;
	uint64_t genome_line_size;

    // Phase 2 variables: Read genome chunks according to MAX_GENOME_CHUNK_SIZE

    bp32_t * ref_seqA;
    bp32_t * ref_seqB;
    uint64_t seq_size;
    uint64_t last_seq_size = 0;
    uint64_t genome_start_ptr;
    uint64_t genome_end_ptr;
    
    uint64_t chunk_num;
    uint64_t last_chunk_size;
    uint64_t ref_seqRev = 0x0;

    // Phase 3 variables: Compute minimizers into the chunks

    uint8_t w,k;
    uint64_t fin_seq_size;
    struct hashmizer_t mini;
    struct hashmizer_t mini_cpy;

    /*  
        ########################################################################  
        ##### 1 - Read the reference genome file and get the genome length #####
        ######################################################################## 
    */

    #ifdef VERBOSE 
        printf("[BuildIndex] Reading ref genome size ... ");
    #endif
        
    fread(&genome_size, sizeof(uint64_t), 1, TF->file); 
    genome_start = 0x8;                                     //First 8 bytes are for the genome size

    TF->seqlen = genome_size;
    TF->seqstart = genome_start;

    #ifdef VERBOSE
        printf(" %010lu bps\n", genome_size);
    #endif

    /*  
        #######################################################################################  
        ##### 2 - Read the reference genome chunks by chunks (reverse complement as well) #####
        #######################################################################################
    */

    #ifdef VERBOSE
        printf("[BuildIndex] Initializing the Index\t");
    #endif

    genome_start_ptr = genome_start;
    genome_end_ptr = genome_start + genome_size/4;
    #ifdef VERBOSE
        printf(" start = %lu\tgeptr = %lu\n", genome_start, genome_end_ptr - genome_start);
    #endif

    seq_size = ( genome_size <= MAX_GENOME_CHUNK_SIZE ) ? genome_size : MAX_GENOME_CHUNK_SIZE ;

    ref_seqA = (bp32_t *)malloc(sizeof(bp32_t)*ceil(seq_size/32));
    ref_seqB = (bp32_t *)malloc(sizeof(bp32_t)*ceil(seq_size/32));

    initHashTable(&(TF->index));

    #ifdef VERBOSE
        printf("Done\n[BuildIndex] Reading genome chunks of size %010lu, with a module of %ld\n", seq_size, genome_size%seq_size);
    #endif

    // For each chunk to be read

    chunk_num = ceil((double)(genome_size)/seq_size);
    last_chunk_size = genome_size%seq_size;

    for( int i = 0; i < chunk_num; i++ ){

        fin_seq_size = ( i == chunk_num-1 && last_chunk_size ) ? last_chunk_size : seq_size;

        //printf("\t\t%d\t%d\n", chunk_num, fin_seq_size);

        #ifdef VERBOSE
            printf("\t[BuildIndex] Working on chunks %04d/%.0f\t",i+1,ceil((double)(genome_size)/seq_size));
        #endif

        // read positive strand chunk and save it into SeqA

        if( genome_start_ptr <= genome_end ) fseek(TF->file, genome_start_ptr, SEEK_SET);

        for( int j = 0; j < ceil(fin_seq_size/32); j++ ){
            fread(&ref_seqA[j], sizeof(uint64_t), 1, TF->file); 	        
        }

        #ifdef VERBOSE
            printf("\t\tChunk Size %08lu\tGSptr = %010lu", fin_seq_size, genome_start_ptr - genome_start);
        #endif

        genome_start_ptr += ( fin_seq_size/4 );  

        // read negative strand chunk and save it into SeqB
        
        genome_end_ptr -= ( fin_seq_size/4 );
 
        fseek(TF->file, genome_end_ptr, SEEK_SET);
        #ifdef VERBOSE
            printf("\tGEptr = %010lu\t\t", genome_end_ptr - genome_start);
        #endif

        for( int j = 0; j < ceil(fin_seq_size/32); j++ ){
            fread(&ref_seqB[j], sizeof(uint64_t), 1, TF->file);       
        }
        reverseComplement_bp32(ref_seqB, fin_seq_size);
        
        /*
            e.g. :          chunk00. chunk01. chunk02. chunk03. chunk04. chunk05.

                ref_seqA :  CGTCGACT GATCGATG CTCGATGC TGATCGAT CGATGCTA GCTGCA##
                ref_seqB :  TAGCTGCA ATCGATGC GCTGATCG TGCTCGAT CTGATCGA CGTCGA##
                Reverse & Complement each chunk
                ref_seqB :  TGCAGCTA GCATCGAT CGATCAGC ATCGAGCA TCGATCAG TCGACG## 
        
        */

        /*  
            #####################################################################  
            ##### 3 - For each couple of chunks (+ & -) find the minimizers #####
            #####################################################################
        */


        w = MINI_W;
        k = MINI_K;

        bp32_t tmpA, tmpB;

        for( int j = 0; j < fin_seq_size - w; j++ ){

            findMinimizer_bp32(ref_seqA, ref_seqB, j, MAX_64BIT, k, w, &mini);

            //Offset minimizer position into the full genome
            if( mini.strand == STRAND_MINUS ) mini.start = genome_size - (i*seq_size + j + mini.start);
            else mini.start = i*seq_size + j + mini.start;

             //If the minimizer is not a duplicate add it to the index
            if( mini_cpy.hash != mini.hash && mini.start != mini_cpy.start ) {
                mini_cpy = mini;
                addToHashTable( &(TF->index), &mini );
            }
          
        }

        #ifdef VERBOSE
            printf("Done \n");
        #endif
        
    }

    #ifdef VERBOSE
        printf("[BuildIndex] Index succesfully built\n");
    #endif

    // Free resources
    free(ref_seqA);
    free(ref_seqB);

}

void findMinimizer_bp32( const bp32_t * seqA, const bp32_t * seqB, uint64_t seqOffset, const uint64_t max, const uint8_t k, const uint8_t w, struct hashmizer_t * mini ){

    // For a detailed explanation of this algorithm check minimap paper.

    // We will assume that W and K won't be greater than 32 bps; we will not need dynamic allocation
    
    uint64_t m;
    uint64_t u,v;
	bp32_t kmerA;
    bp32_t kmerB;
    uint8_t kmerStartOffset;
    uint8_t kmerEndOffset;

    m = max;

    for( int j = 0; j < w-k+1; j++ ){

        kmerA = 0;
        kmerB = 0;
        kmerStartOffset = j + seqOffset%32;

        if( k + kmerStartOffset < 32 ){         //Kmer fits in this 32bp word

            kmerA = (seqA[seqOffset/32] >> 2*kmerStartOffset) & (bp32_t)(pow(2,2*k)-1);
            kmerB = (seqB[seqOffset/32] >> 2*kmerStartOffset) & (bp32_t)(pow(2,2*k)-1);

        }
        else if( kmerStartOffset < 32 ){        //Kmer ovlps between 2 32bp words

            kmerA = (seqA[seqOffset/32] >> 2*kmerStartOffset);
            kmerA |= ((seqA[seqOffset/32 + 1] & (bp32_t)(pow(2,2*(k - (32 - kmerStartOffset)))-1)) << 2*(32 - kmerStartOffset));
            kmerA &= (bp32_t)(pow(2,2*k)-1);

            kmerB = (seqB[seqOffset/32] >> 2*kmerStartOffset);
            kmerB |= ((seqB[seqOffset/32 + 1] & (bp32_t)(pow(2,2*(k - (32 - kmerStartOffset)))-1)) << 2*(32 - kmerStartOffset));
            kmerB &= (bp32_t)(pow(2,2*k)-1);    

        }
        else{                                   //Kmer fits in the next 32bp word
            
            kmerA = (seqA[seqOffset/32+1] >> 2*(kmerStartOffset%32)) & (bp32_t)(pow(2,2*k)-1);
            kmerB = (seqB[seqOffset/32+1] >> 2*(kmerStartOffset%32)) & (bp32_t)(pow(2,2*k)-1);

        }

        u = hashSequence_bp32(kmerA, k);
		v = hashSequence_bp32(kmerB, k);
		if( u != v ) m = min64( m, min64(u,v) );

    }

	for( int j = 0; j < w-k+1; j++ ){
		
		kmerA = 0;
        kmerB = 0;
        kmerStartOffset = j + seqOffset%32;

        if( k + kmerStartOffset < 32 ){         //Kmer fits in this 32bp word

            kmerA = (seqA[seqOffset/32] >> 2*kmerStartOffset) & (bp32_t)(pow(2,2*k)-1);
            kmerB = (seqB[seqOffset/32] >> 2*kmerStartOffset) & (bp32_t)(pow(2,2*k)-1);

        }
        else if( kmerStartOffset < 32 ){        //Kmer ovlps between 2 32bp words

            kmerA = (seqA[seqOffset/32] >> 2*kmerStartOffset);
            kmerA |= ((seqA[seqOffset/32 + 1] & (bp32_t)(pow(2,2*(k - (32 - kmerStartOffset)))-1)) << 2*(32 - kmerStartOffset));
            kmerA &= (bp32_t)(pow(2,2*k)-1);

            kmerB = (seqB[seqOffset/32] >> 2*kmerStartOffset);
            kmerB |= ((seqB[seqOffset/32 + 1] & (bp32_t)(pow(2,2*(k - (32 - kmerStartOffset)))-1)) << 2*(32 - kmerStartOffset));
            kmerB &= (bp32_t)(pow(2,2*k)-1);         

        }
        else{                                   //Kmer fits in the next 32bp word
            
            kmerA = (seqA[seqOffset/32+1] >> 2*(kmerStartOffset%32)) & (bp32_t)(pow(2,2*k)-1);
            kmerB = (seqB[seqOffset/32+1] >> 2*(kmerStartOffset%32)) & (bp32_t)(pow(2,2*k)-1);
        }

        u = hashSequence_bp32(kmerA, k);
		v = hashSequence_bp32(kmerB, k);
        
		mini->len = k;
		mini->start = j;
			
		if( u < v && u == m ) {
            mini->strand = STRAND_PLUS;
            mini->hash = u;
        }
		else if( u > v && v == m ){
            mini->strand = STRAND_MINUS;
            mini->hash = v;
        }

	}

   

}

#endif