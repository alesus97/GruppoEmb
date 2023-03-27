
#ifndef HASH
#define HASH

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "types.h"
#include "misc.h"

//Hash functions

void initHashTable( struct hashband_t * ht );
void addToHashTable( struct hashband_t * ht, const struct hashmizer_t * hm );
void queryHashTable( const struct hashband_t * ht, const struct hashmizer_t * hm, cvector_vector_type(struct minimatch_t) * matches );
void freeHashTable( struct hashband_t * ht );

#ifndef MODE_COMPRESSED
    uint64_t hashSequence( const char * read, const uint32_t seqlen );
#else
    uint64_t hashSequence_bp32( const bp32_t read, const uint32_t seqlen );
#endif

uint64_t hashInv( uint64_t key );

//List functions

void insertFirstNode( struct hashband_t * ht, const struct hashmizer_t * hm, const uint64_t index );
struct hashnode_t* findNode( const struct hashband_t * ht, const struct hashmizer_t * hm, const uint32_t index, uint64_t key );



#endif
