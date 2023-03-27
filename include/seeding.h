#include "hash.h"

#ifndef MODE_COMPRESSED

void BuildIndex( struct seqfile_t * TF );
void findMinimizer( const bp_t * seqA, const bp_t * seqB, const uint64_t seqOffsetA, const uint64_t seqOffsetB, const uint64_t max, const uint8_t k, const uint8_t w, struct hashmizer_t * mini );

#else

void BuildIndex_bp32( struct seqfile_t * TF );
void findMinimizer_bp32( const bp32_t * seqA, const bp32_t * seqB, uint64_t seqOffset, const uint64_t max, const uint8_t k, const uint8_t w, struct hashmizer_t * mini );

#endif
