
#ifndef MAP
#define MAP

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "types.h"
#include "hash.h"
#include "misc.h"
#include "seeding.h"

#ifdef MULTITHREADING

    #include "mthread.h"

#endif

#ifndef MODE_COMPRESSED

void MapReadsToGenome( struct seqfile_t * TF, struct seqfile_t * RF, FILE * SAMfile );
void findAndMatchChunkMinimizers(bp_t * chunkA, bp_t * chunkB, const uint32_t chunk_size, uint32_t chunkID, const uint8_t k, const uint8_t w,  cvector_vector_type(struct minimatch_t) * matches, struct seqfile_t * RF, struct seqfile_t * TF);
void MapChunksToGenome(bp_t * chunkA, bp_t * chunkB, const uint32_t chunk_num, const uint32_t chunk_size, const uint8_t k, const uint8_t w, cvector_vector_type(struct minimatch_t) * matches, struct seqfile_t * RF, struct seqfile_t * TF);
void chunkFilter(cvector_vector_type(struct minimatch_t) * matches, struct seqfile_t * RF, struct seqfile_t * TF);

#else

void MapReadsToGenome_bp32( struct seqfile_t * TF, struct seqfile_t * RF, FILE * SAMfile );
void findAndMatchChunkMinimizers_bp32(bp32_t *chunkA, bp32_t *chunkB, const uint32_t chunk_size, uint32_t chunkID, const uint8_t k, const uint8_t w, cvector_vector_type(struct minimatch_t) * matches, struct seqfile_t *RF, struct seqfile_t *TF);
void MapChunksToGenome_bp32(bp32_t * chunkA, bp32_t * chunkB, const uint32_t chunk_num, const uint32_t chunk_size, const uint8_t k, const uint8_t w, cvector_vector_type(struct minimatch_t) * matches, struct seqfile_t * RF, struct seqfile_t * TF);
void chunkFilter_bp32(cvector_vector_type(struct minimatch_t) * matches, struct seqfile_t * RF, struct seqfile_t * TF);

#endif

#endif
