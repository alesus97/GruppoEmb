
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
#include "mthread.h"

#include "../WFA2/utils/commons.h"
#include "../WFA2/wavefront/wavefront_align.h"
#include "../include/SneakySnake.h"

#ifdef MULTITHREADING

    //#include "mthread.h"

#endif

#ifndef MODE_COMPRESSED

void MapReadsToGenome(struct seqfile_t *TF, struct seqfile_t *RF, FILE *SAMfile, struct mapper_ctx_t * thread_ctx);
void findAndMatchChunkMinimizers(bp_t * chunkA, bp_t * chunkB, const uint32_t chunk_size, uint32_t chunkID, const uint8_t k, const uint8_t w,  cvector_vector_type(struct minimatch_t) * matches, struct seqfile_t * RF, struct seqfile_t * TF);
void MapChunksToGenome(bp_t * chunkA, bp_t * chunkB, const uint32_t chunk_num, const uint32_t chunk_size, const uint8_t k, const uint8_t w, cvector_vector_type(struct minimatch_t) * matches, struct seqfile_t * RF, struct seqfile_t * TF,  struct mapper_ctx_t * thread_ctx);
void chunkElaboration(cvector_vector_type(struct minimatch_t) * matches, struct seqfile_t * RF, struct seqfile_t * TF,  struct mapper_ctx_t * thread_ctx);
int chunkFilter(struct seqfile_t *RF, struct seqfile_t *TF, uint32_t chunk_num, struct minimatch_t minmatch,  uint32_t leftOffsetT, uint32_t last_chunk_size, uint32_t chunk_size,  struct mapper_ctx_t * thread_ctx);
void chunkAlign(struct seqfile_t *RF, struct seqfile_t *TF, uint32_t chunk_num, struct minimatch_t minmatch, uint32_t leftOffsetT, uint32_t last_chunk_size, uint32_t chunk_size, wavefront_aligner_attr_t attributes, wavefront_aligner_t* const wf_aligner,struct mapper_ctx_t * thread_ctx);

#else

void MapReadsToGenome_bp32( struct seqfile_t * TF, struct seqfile_t * RF, FILE * SAMfile );
void findAndMatchChunkMinimizers_bp32(bp32_t *chunkA, bp32_t *chunkB, const uint32_t chunk_size, uint32_t chunkID, const uint8_t k, const uint8_t w, cvector_vector_type(struct minimatch_t) * matches, struct seqfile_t *RF, struct seqfile_t *TF);
void MapChunksToGenome_bp32(bp32_t * chunkA, bp32_t * chunkB, const uint32_t chunk_num, const uint32_t chunk_size, const uint8_t k, const uint8_t w, cvector_vector_type(struct minimatch_t) * matches, struct seqfile_t * RF, struct seqfile_t * TF);
void chunkElaboration_bp32(cvector_vector_type(struct minimatch_t) * matches, struct seqfile_t * RF, struct seqfile_t * TF);
int chunkFilter_bp32(struct seqfile_t *RF, struct seqfile_t *TF, uint32_t chunk_num, struct minimatch_t minmatch,  uint32_t leftOffsetT, uint32_t last_chunk_size, uint32_t chunk_size);
#endif

#endif
