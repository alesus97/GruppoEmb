#include "map.h"
#include <inttypes.h>


uint64_t totalAccepted = 0;

#ifdef MULTITHREADING

	extern pthread_mutex_t TF_mutex;
	
#endif

#ifndef MODE_COMPRESSED

	bp_t chunkQ[MAX_GENOME_CHUNK_SIZE]; // Chunk for the Query  (i.e. the read)
	bp_t chunkT[MAX_GENOME_CHUNK_SIZE]; // Chunk for the Target (i.e. the reference)
	bp_t temp_chunk[MAX_GENOME_CHUNK_SIZE]; 

void MapReadsToGenome(struct seqfile_t *TF, struct seqfile_t *RF, FILE *SAMfile)
{

	bp_t *chunkA;
	bp_t *chunkB;
	uint32_t chunk_size;
	uint32_t prev_chunk_size = 0;
	uint32_t chunk_num;

	uint8_t w, k;

	cvector_vector_type(struct minimatch_t) matches = NULL;

	w = MINI_W;
	k = MINI_K;

	// For all reads in the file

	while (fgetc(RF->file) != EOF)
	{

		/*
			########################################
			##### 1 - Extract read chunks      #####
			########################################
		*/

		getNextRead(RF);

		/*#ifdef VERBOSE
			printf("[MapReadsToGenome] Mapping read <%010u>\n", RF->seqid);
		#endif*/

		chunk_size = (RF->seqlen <= MAX_READ_CHUNK_SIZE) ? RF->seqlen : MAX_READ_CHUNK_SIZE;
		chunk_num = ceil((double)(RF->seqlen) / chunk_size) ;


		// Check if chunk size contains at least one window
		if (RF->seqlen % chunk_size < w && RF->seqlen % chunk_size != 0)
		{
			chunk_size += w + 1;
			chunk_num -= 1;
		}

		// Allocate chunks whenever the size changes
		if (prev_chunk_size < chunk_size)
		{

			if (RF->seqid != 1)
			{
				free(chunkA);
				free(chunkB);
			}

			chunkA = (bp_t *)malloc(sizeof(bp_t) * chunk_size);
			chunkB = (bp_t *)malloc(sizeof(bp_t) * chunk_size);
			prev_chunk_size = chunk_size;
		}

		/*
			########################################################
			##### 2 - For each read, perform mapping on chunks #####
			########################################################
		*/

		MapChunksToGenome(chunkA, chunkB, chunk_num, chunk_size, k, w, &matches, RF, TF);

		if (fgetc(RF->file) == EOF)
			break;
	}

	printf("totalAccepted %d \n", totalAccepted);
	// free area
	free(chunkA);
	free(chunkB);
	cvector_free(matches);
}

void MapChunksToGenome(bp_t *chunkA, bp_t *chunkB, const uint32_t chunk_num, const uint32_t chunk_size, const uint8_t k, const uint8_t w, cvector_vector_type(struct minimatch_t) * matches, struct seqfile_t *RF, struct seqfile_t *TF)
{

	uint32_t last_chunk_size = 0;
	uint32_t effective_chunk_size;

	// for all chunks in the read

	for (int i = 0; i < chunk_num; i++)
	{

		// Extract read chunks

		if (i == chunk_num - 1)
		{
			last_chunk_size = RF->seqlen - (i)*chunk_size;
		}

		if (!last_chunk_size)
		{
			getReadChunk(RF, chunkA, i * chunk_size, chunk_size);
			getReadChunk(RF, chunkB, RF->seqlen - (i + 1) * chunk_size, chunk_size);
		}
		else
		{
			getReadChunk(RF, chunkB, 0, last_chunk_size);
			getReadChunk(RF, chunkA, RF->seqlen - last_chunk_size, last_chunk_size);
		}

		/*
			###############################################
			##### 3 - For each chunk, find minimizers #####
			###############################################
		*/

		effective_chunk_size = (last_chunk_size) ? last_chunk_size : chunk_size;

		// find chunked read (RF) minimizers and match them against the reference genome (TF). Matching minimizers are saved into the matches vector
		findAndMatchChunkMinimizers(chunkA, chunkB, chunk_size, i, k, w, matches, RF, TF);
		// printf("Read[%08u::%04u] matched [%08lu] times\n", RF->seqid, i, cvector_size(*matches));

		/*
			###############################################################
			##### 4 - For each match, perform filtering and alignment #####
			###############################################################
		*/

		// filter matches using the SneakySnake algorithm
		chunkElaboration(matches, RF, TF);

		/*
			##########################################################
			##### 5 - For each filtered match, perform alignment #####
			##########################################################
		*/

		last_chunk_size = 0;
	}
}

void findAndMatchChunkMinimizers(bp_t *chunkA, bp_t *chunkB, const uint32_t chunk_size, uint32_t chunkID, const uint8_t k, const uint8_t w, cvector_vector_type(struct minimatch_t) * matches, struct seqfile_t *RF, struct seqfile_t *TF){

	struct hashmizer_t mini;
	struct hashmizer_t mini_cpy;

	for (int j = 0; j < chunk_size - w; j++)
	{

		findMinimizer(chunkA + j, chunkB + chunk_size - w - j, j, chunk_size - w - j, MAX_64BIT, k, w, &mini);

		// Offset minimizer position into the full read
		if (mini.strand == STRAND_MINUS)
			mini.start = RF->seqlen - chunkID * chunk_size - mini.start;
		else
			mini.start = chunkID * chunk_size + mini.start;

		// If the minimizer is not a duplicate then proceed
		if (mini_cpy.hash != mini.hash && mini.start != mini_cpy.start)
		{
			mini_cpy = mini;

			// Query the index : returns all the locations in the ref genome where a match occurs
			// for cvector_vector_type look in types.h

			#ifdef MULTITHREADING
				pthread_mutex_lock(&TF_mutex);
				queryHashTable(&(TF->index), &mini, matches);
				pthread_mutex_unlock(&TF_mutex);
			#else
				queryHashTable(&(TF->index), &mini, matches);
			#endif
		}
	}
}


void strandTNegative(struct seqfile_t *RF, struct seqfile_t *TF, uint32_t chunk_num, struct minimatch_t minmatch, uint32_t leftOffsetT, uint64_t chunk_size, uint64_t last_chunk_size, int i){

	minmatch.startT = TF->seqlen - minmatch.startT;
	leftOffsetT = minmatch.startT - minmatch.startQ;
	getReadChunk(TF, temp_chunk, leftOffsetT + (chunk_size * (chunk_num - 1 - i)), last_chunk_size);
	reverseComplement(temp_chunk, chunkT, last_chunk_size);	
	
}

void strandTPositive(struct seqfile_t *RF, struct seqfile_t *TF, uint32_t chunk_num, struct minimatch_t minmatch, uint32_t leftOffsetT, uint64_t chunk_size, uint64_t last_chunk_size, int i){
	
	leftOffsetT = minmatch.startT - minmatch.startQ;
	getReadChunk(TF, chunkT, leftOffsetT + (chunk_size * i), last_chunk_size);

}


int chunkFilter(struct seqfile_t *RF, struct seqfile_t *TF, uint32_t chunk_num, struct minimatch_t minmatch, uint32_t leftOffsetT, uint32_t last_chunk_size, uint32_t chunk_size){

		uint32_t accepted = 0;

	 for (uint32_t i = 0; i < chunk_num; i++){

				// If strandQ is Negative
				if(minmatch.strandQ == 1){

					if (i == 0 && (RF->seqlen % chunk_size) != 0){
						last_chunk_size = (RF->seqlen) % chunk_size;

					}else{
						last_chunk_size = chunk_size;
					}

					getReadChunk(RF, temp_chunk, (chunk_size * (chunk_num - 1 - i)), last_chunk_size);
					//printSequence(chunkQ, 0, last_chunk_size, 1);
					reverseComplement(temp_chunk, chunkQ, last_chunk_size);

				} // If strandQ is Positive
				else if(minmatch.strandQ == 0){

					if (i == (chunk_num - 1) && (RF->seqlen % chunk_size) != 0){
						last_chunk_size = (RF->seqlen) % chunk_size;
					}else{
						last_chunk_size = chunk_size;
					}

					getReadChunk(RF, chunkQ, (chunk_size * i), last_chunk_size);
					//printSequence(chunkQ, 0, last_chunk_size, 1);

				}

				// If strandT is Negative
				if (minmatch.strandT == 1){

					if (i == 0 && (RF->seqlen % chunk_size) != 0){
						last_chunk_size = (RF->seqlen) % chunk_size;
					}else{
						last_chunk_size = chunk_size;
					}
					
					strandTNegative(RF,TF,chunk_num,minmatch,leftOffsetT,chunk_size, last_chunk_size, i);	
 
				} // If strandQ is Positive
				else if (minmatch.strandT == 0){

					if (i == (chunk_num - 1) && (RF->seqlen % chunk_size) != 0){
						last_chunk_size = (RF->seqlen) % chunk_size;
					}else{
						last_chunk_size = chunk_size;
					}

					strandTPositive(RF, TF, chunk_num, minmatch, leftOffsetT, chunk_size , last_chunk_size, i);
					
				}


				/* printf("strandQ %d \n", minmatch.strandQ);
				printf("strandT %d \n", minmatch.strandT);
				printf("last_chunk_size %d \n", last_chunk_size); */
				/* printSequence(chunkQ, 0, last_chunk_size, 1);
				printSequence(chunkT, 0, last_chunk_size, 1); */
 

				printf("lastChunkSize %d \n", last_chunk_size);
				 if(SneakySnake(last_chunk_size, chunkT, chunkQ, last_chunk_size/100, ceil(last_chunk_size/10.0), 0, 1)){
					accepted++;
				}else {
					printf("Chunk Accettati %d/%d \n", accepted, chunk_num);
					//printf("lastChunkSize %d \n", last_chunk_size);
					return 0;
				}
 
			}

	printf("Chunk Accettati %d/%d \n", accepted, chunk_num);

	return 1;

}


/* 	  wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
	wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default; */

void chunkAlign(struct seqfile_t *RF, struct seqfile_t *TF, uint32_t chunk_num, struct minimatch_t minmatch, uint32_t leftOffsetT, uint32_t last_chunk_size, uint32_t chunk_size, wavefront_aligner_attr_t attributes, wavefront_aligner_t* const wf_aligner){

		 	for(int i=0; i<chunk_num; i++){

				if(minmatch.strandQ == 1){

					if (i == 0){
						last_chunk_size = (RF->seqlen) % chunk_size;
					}else{
						last_chunk_size = chunk_size;
					}

					getReadChunk(RF, temp_chunk, (chunk_size * (chunk_num - 1 - i)), last_chunk_size);
					reverseComplement(temp_chunk, chunkQ, last_chunk_size);

				} // If strandQ is Positive
				else if(minmatch.strandQ == 0){

					if (i == (chunk_num - 1)){
						last_chunk_size = (RF->seqlen) % chunk_size;
					}else{
						last_chunk_size = chunk_size;
					}

					getReadChunk(RF, chunkQ, (chunk_size * i), last_chunk_size);

				}

					for(int j=0; j<chunk_num; j++){ 

								// If strandT is Negative
						if (minmatch.strandT == 1){

							if (j == 0){
								last_chunk_size = (RF->seqlen) % chunk_size;
							}else{
								last_chunk_size = chunk_size;
							}
							
							strandTNegative(RF,TF,chunk_num,minmatch,leftOffsetT,chunk_size, last_chunk_size, j);
		
						} // If strandQ is Positive
						else if (minmatch.strandT == 0){

							if (i == (chunk_num - 1)){
								last_chunk_size = (RF->seqlen) % chunk_size;
							}else{
								last_chunk_size = chunk_size;
							}

							strandTPositive(RF, TF, chunk_num, minmatch, leftOffsetT, chunk_size , last_chunk_size, j);
							
						}

						printf("ChunkT elaborati %d/%d \n",j, chunk_num - 1);

						//Allineo	
						wavefront_align(wf_aligner,chunkQ,last_chunk_size,chunkT,last_chunk_size);
					   	fprintf(stderr,"WFA-Alignment returns score %d\n",wf_aligner->cigar->score);
						fprintf(stderr,"  SCORE (RE)COMPUTED %d\n",
						cigar_score_gap_affine(wf_aligner->cigar,&attributes.affine_penalties));
					
					}		
			 } 
}



void chunkElaboration(cvector_vector_type(struct minimatch_t) * matches, struct seqfile_t *RF, struct seqfile_t *TF){

	struct minimatch_t minmatch;
	uint32_t leftOffsetT;
	uint32_t chunk_num;
	
	
	wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.distance_metric = gap_affine;
 	attributes.affine_penalties.match = 0;
 	attributes.affine_penalties.mismatch = 4;
    attributes.affine_penalties.gap_opening = 6;
    attributes.affine_penalties.gap_extension = 2;

   // Initialize Wavefront Aligner
    wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
	
	while (cvector_size(*matches))
	{

		uint32_t last_chunk_size = 0;
		uint32_t chunk_size = MAX_GENOME_CHUNK_SIZE;
		chunk_num = (RF->seqlen <= MAX_GENOME_CHUNK_SIZE) ? 1 : ceil(RF->seqlen / (MAX_GENOME_CHUNK_SIZE));

		//printf("chunkNum: %ld \n", chunk_num);
			
		minmatch = *(*matches + (cvector_size(*matches) - 1));
		
		if(chunkFilter(RF, TF, chunk_num, minmatch, leftOffsetT, last_chunk_size, chunk_size)){
			chunkAlign(RF, TF, chunk_num, minmatch, leftOffsetT, last_chunk_size, chunk_size, attributes, wf_aligner);
		}

		cvector_pop_back(*matches);
				
	}

	// Free
	wavefront_aligner_delete(wf_aligner);   
		
}

#else

	bp32_t chunkQ[MAX_GENOME_CHUNK_SIZE/32]; // Chunk for the Query  (i.e. the read)
	bp32_t chunkT[MAX_GENOME_CHUNK_SIZE/32]; // Chunk for the Target (i.e. the reference)
	bp32_t temp_chunk[MAX_GENOME_CHUNK_SIZE/32]; 
	


void MapReadsToGenome_bp32(struct seqfile_t *TF, struct seqfile_t *RF, FILE *SAMfile)
{

	bp32_t *chunkA;
	bp32_t *chunkB;
	uint32_t chunk_size;
	uint32_t prev_chunk_size = 0;
	uint32_t chunk_num;
	uint64_t read_num;

	uint8_t w, k;

	cvector_vector_type(struct minimatch_t) matches = NULL;

	w = MINI_W;
	k = MINI_K;

	// For all reads in the file

	fread(&read_num, sizeof(uint64_t), 1, RF->file);

	#ifdef VERBOSE
		printf("[MapReadsToGenome] Starting mapping process for <%lu> reads\n", read_num);
	#endif

	for( uint64_t i = 0; i < read_num; i++ )
	{

		/*
			########################################
			##### 1 - Extract read chunks      #####
			########################################
		*/

		getNextRead_bp32(RF);
		chunk_size = (RF->seqlen <= MAX_READ_CHUNK_SIZE) ? RF->seqlen : MAX_READ_CHUNK_SIZE;	
	
		// Check if chunk size contains at least one window
		/*if (RF->seqlen % chunk_size < w && RF->seqlen % chunk_size != 0)
		{
			chunk_size += w + 1;
		}*/

		chunk_num = ceil((double)(RF->seqlen) / chunk_size);

		// Allocate chunks whenever the size changes
		if (prev_chunk_size < chunk_size)
		{

			if (RF->seqid != 1)
			{
				free(chunkA);
				free(chunkB);
			}

			chunkA = (bp32_t *)malloc(sizeof(bp32_t) * ceil(chunk_size/32));
			chunkB = (bp32_t *)malloc(sizeof(bp32_t) * ceil(chunk_size/32));
			prev_chunk_size = chunk_size;
		}
		
		/*
			########################################################
			##### 2 - For each read, perform mapping on chunks #####
			########################################################
		*/
		
		MapChunksToGenome_bp32(chunkA, chunkB, chunk_num, chunk_size, k, w, &matches, RF, TF);
	}	

	printf("totalAccepted %d \n", totalAccepted);
	// free area
	free(chunkA);
	free(chunkB);
	cvector_free(matches);
}

void MapChunksToGenome_bp32(bp32_t *chunkA, bp32_t *chunkB, const uint32_t chunk_num, const uint32_t chunk_size, const uint8_t k, const uint8_t w, cvector_vector_type(struct minimatch_t) * matches, struct seqfile_t *RF, struct seqfile_t *TF)
{

	uint32_t last_chunk_size = 0;
	uint32_t effective_chunk_size;

	// for all chunks in the read

	for (int i = 0; i < chunk_num; i++)
	{

		// Extract read chunks

		if (i == chunk_num - 1){
			last_chunk_size = RF->seqlen - (i)*chunk_size;
		}

		if (!last_chunk_size){
			getReadChunk_bp32(RF, chunkA, i * chunk_size, chunk_size);
			getReadChunk_bp32(RF, chunkB, RF->seqlen - (i + 1) * chunk_size, chunk_size);
			reverseComplement_bp32(chunkB, chunk_size);
		}
		else{
			getReadChunk_bp32(RF, chunkB, 0, last_chunk_size);
			getReadChunk_bp32(RF, chunkA, RF->seqlen - last_chunk_size, last_chunk_size);
			reverseComplement_bp32(chunkB, chunk_size);
		}

		/*
			###############################################
			##### 3 - For each chunk, find minimizers #####
			###############################################
		*/

		effective_chunk_size = (last_chunk_size) ? last_chunk_size : chunk_size;

		// find chunked read (RF) minimizers and match them against the reference genome (TF). Matching minimizers are saved into the matches vector

		findAndMatchChunkMinimizers_bp32(chunkA, chunkB, chunk_size, i, k, w, matches, RF, TF);

		/*
			###############################################################
			##### 4 - For each match, perform filtering and alignment #####
			###############################################################
		*/


		// filter matches using the SneakySnake algorithm
		chunkElaboration_bp32(matches, RF, TF);

		/*
			##########################################################
			##### 5 - For each filtered match, perform alignment #####
			##########################################################
		*/

		last_chunk_size = 0;
	}
}

void findAndMatchChunkMinimizers_bp32(bp32_t *chunkA, bp32_t *chunkB, const uint32_t chunk_size, uint32_t chunkID, const uint8_t k, const uint8_t w, cvector_vector_type(struct minimatch_t) * matches, struct seqfile_t *RF, struct seqfile_t *TF)
{
	
	struct hashmizer_t mini;
	struct hashmizer_t mini_cpy;

	for (int j = 0; j < chunk_size - w; j++)
	{

		findMinimizer_bp32(chunkA, chunkB, j, MAX_64BIT, k, w, &mini);
		
		// Offset minimizer position into the full read
		if( mini.strand == STRAND_MINUS ) mini.start = RF->seqlen - (chunkID * chunk_size + j + mini.start);
        else mini.start = chunkID * chunk_size + j + mini.start;

		// If the minimizer is not a duplicate then proceed
		if (mini_cpy.hash != mini.hash && mini.start != mini_cpy.start)
		{
			mini_cpy = mini;

			// Query the index : returns all the locations in the ref genome where a match occurs
			// for cvector_vector_type look in types.h

			// TF is a thread-shared object, hence we must define a critical section for it if multithreading is active 

			#ifdef MULTITHREADING
				pthread_mutex_lock(&TF_mutex);
				queryHashTable(&(TF->index), &mini, matches);
				pthread_mutex_unlock(&TF_mutex);
			
			#else
				queryHashTable(&(TF->index), &mini, matches);
			#endif
			
		}
	}
}




void strandTPositive_bp32(struct seqfile_t *RF, struct seqfile_t *TF, uint32_t chunk_num, struct minimatch_t minmatch, uint32_t leftOffsetT, uint64_t chunk_size, uint64_t last_chunk_size, int i){
	
	leftOffsetT = minmatch.startT - minmatch.startQ;

		/* #ifdef MULTITHREADING
					pthread_mutex_lock(&TF_mutex);
					getReadChunk_bp32(TF, chunkT, leftOffsetT + (chunk_size * i), last_chunk_size);
					pthread_mutex_unlock(&TF_mutex);
		#else */
					getReadChunk_bp32(TF, chunkT, leftOffsetT + (chunk_size * i), last_chunk_size);
	//	#endif
}

void strandTNegative_bp32(struct seqfile_t *RF, struct seqfile_t *TF, uint32_t chunk_num, struct minimatch_t minmatch, uint32_t leftOffsetT, uint64_t chunk_size, uint64_t last_chunk_size, int i){

	minmatch.startT = TF->seqlen - minmatch.startT;
	leftOffsetT = minmatch.startT - minmatch.startQ;

		/* #ifdef MULTITHREADING
					pthread_mutex_lock(&TF_mutex);
					getReadChunk_bp32(TF, chunkT, leftOffsetT + (chunk_size * (chunk_num - 1 - i)), last_chunk_size);
					pthread_mutex_unlock(&TF_mutex);
		#else */
					getReadChunk_bp32(TF, chunkT, leftOffsetT + (chunk_size * (chunk_num - 1 - i)), last_chunk_size);
		//#endif
	
	reverseComplement_bp32(chunkT, last_chunk_size);	

}

int chunkFilter_bp32(struct seqfile_t *RF, struct seqfile_t *TF, uint32_t chunk_num, struct minimatch_t minmatch, uint32_t leftOffsetT, uint32_t last_chunk_size, uint32_t chunk_size){

		uint32_t accepted = 0;

	 for (uint32_t i = 0; i < chunk_num; i++){

				// If strandQ is Negative
				if(minmatch.strandQ == 1){

					if (i == 0 && (RF->seqlen % chunk_size) != 0){
						last_chunk_size = (RF->seqlen) % chunk_size;
					}else{
						last_chunk_size = chunk_size;
					}

				} // If strandQ is Positive
				else if(minmatch.strandQ == 0){

					if (i == (chunk_num - 1) && (RF->seqlen % chunk_size) != 0){
						last_chunk_size = (RF->seqlen) % chunk_size;
					}else{
						last_chunk_size = chunk_size;
					}

					getReadChunk_bp32(RF, chunkQ, (chunk_size * i), last_chunk_size);

				}

				// If strandT is Negative
				if (minmatch.strandT == 1){

					if (i == 0 && (RF->seqlen % chunk_size) != 0){
						last_chunk_size = (RF->seqlen) % chunk_size;
					}else{
						last_chunk_size = chunk_size;
					}

					strandTNegative_bp32(RF,TF,chunk_num,minmatch,leftOffsetT,chunk_size, last_chunk_size, i);	
 
				} // If strandQ is Positive
				else if (minmatch.strandT == 0){

					if (i == (chunk_num - 1) && (RF->seqlen % chunk_size) != 0){
						last_chunk_size = (RF->seqlen) % chunk_size;
					}else{
						last_chunk_size = chunk_size;
					}
					
					strandTPositive_bp32(RF, TF, chunk_num, minmatch, leftOffsetT, chunk_size , last_chunk_size, i);

				}
 
				/* printf("strandQ %d \n", minmatch.strandQ);
				printf("strandT %d \n", minmatch.strandT);
				printf("startQ %d \n", minmatch.startQ);
				printf("startT %d \n", minmatch.startT);
				printf("last_chunk_size %d \n", last_chunk_size); */
				/* printSequence_bp32(chunkQ, 0, last_chunk_size, 1);
				printSequence_bp32(chunkT, 0, last_chunk_size, 1);  */

				 if(SneakySnake_bp32(last_chunk_size, chunkT, chunkQ, last_chunk_size/100, ceil(last_chunk_size/10.0), 0, 1)){
					accepted++;
				}else {
					printf("Chunk Accettati %d/%d \n", accepted, chunk_num);
					return 0;
				}   

			}

	printf("Chunk Accettati %d/%d \n", accepted, chunk_num);

	return 1;

}


void chunkElaboration_bp32(cvector_vector_type(struct minimatch_t) * matches, struct seqfile_t *RF, struct seqfile_t *TF){

	
	uint32_t chunk_num;
	struct minimatch_t minmatch;
	uint32_t leftOffsetT;
	uint32_t last_chunk_size = 0;
	while (cvector_size(*matches))
	{
		uint32_t last_chunk_size = 0;
		uint32_t chunk_size = MAX_GENOME_CHUNK_SIZE;
		chunk_num = (RF->seqlen <= MAX_GENOME_CHUNK_SIZE) ? 1 : ceil((double) RF->seqlen / (MAX_GENOME_CHUNK_SIZE));
		minmatch = *(*matches + (cvector_size(*matches) - 1));

		
		/* #ifdef MULTITHREADING
					pthread_mutex_lock(&TF_mutex);
					if(chunkFilter_bp32(RF, TF, chunk_num, minmatch, leftOffsetT, last_chunk_size, chunk_size)){
						printf("Allineamento \n");
					}
					pthread_mutex_unlock(&TF_mutex);
		#else */


		if(chunkFilter_bp32(RF, TF, chunk_num, minmatch, leftOffsetT, last_chunk_size, chunk_size)){
			printf("Allineamento \n");
		} 
	/* 	#endif */

		cvector_pop_back(*matches);
	}

	// Free
	//wavefront_aligner_delete(wf_aligner); 
}



#endif