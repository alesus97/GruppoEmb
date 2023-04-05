
#ifndef MISC
#define MISC

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "types.h"

/* Hash-Band related functions */

void printHashmizer( const struct hashmizer_t * mini );
void printHashTable( const struct hashband_t * ht, uint8_t param );
void printList( const struct hashband_t * ht, const uint32_t index );
int getListLength( const struct hashband_t * ht, const uint32_t index );

/* sequences related functions */

uint64_t min64( uint64_t a, uint64_t b);

#ifndef MODE_COMPRESSED

    void printSequence( const char * seq, uint64_t start, uint64_t end, uint8_t newline );
    void reverseComplement( char * A, char * B, const uint32_t L );

    /* File sequence reading - FASTA */

    void getNextRead(struct seqfile_t * RF);
    void getReadChunk(struct seqfile_t * RF, bp_t * chunk, uint32_t start, uint32_t size);
    uint32_t totalReadLength(struct seqfile_t * RF);

#else

    void printSequence_bp32( const bp32_t * seq, uint64_t start, uint64_t end, uint8_t newline );
    void reverseComplement_bp32( bp32_t * seq, const uint32_t seq_size );

    /* File sequence reading - bp32 */

    void getNextRead_bp32(struct seqfile_t * RF);
    void getReadChunk_bp32(struct seqfile_t * RF, bp32_t * chunk, uint32_t start, uint32_t size);

#endif


void printLogo();


#endif
