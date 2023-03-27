#ifndef TYPE
#define TYPE

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "param.h"
#include "cvector.h"

#ifdef MULTITHREADING

	#include <pthread.h>
    #include <sched.h>
    #include <sys/sysinfo.h>
	
#endif


/* 	Notes on cvector_vector_type 

	Definition:	cvector_vector_type(TYPENAME) VECTORNAME = NULL;
	Append:		cvector_push_back(VECTORNAME, ELEMENT);
	pop:		cvector_pop_back(VECTORNAME)
	Free:		cvector_free(VECTORNAME);

	It is also possible to Insert in any position, copy the vector etc. For more details look at cvector.h
*/

#ifdef MODE_COMPRESSED
	typedef uint64_t bp32_t;	// {A,C,G,T} = {00,01,10,11}, 32 bases per word
#else
	typedef char bp_t;
#endif


/* Hashmizer: it is a minimizer structure. See minimap2 paper for formal modeling and definition */

struct hashmizer_t{

	uint64_t hash;		// 64-bit hash of the k-mer
	uint8_t len;		// k-mer length (i.e. MINI_K)
	uint64_t start;		// Starting point into the sequence where the minimizer was extracted from
	uint8_t strand;		// Positive strand or Negative strand
};

/* Hashnode: it is a node in the hashband linked lists */

struct hashnode_t {
   struct hashmizer_t hashmizer;	// Hashmizer for the node
   struct hashnode_t * next;		// Pointer to the next node
   uint32_t list_offset;			// Position of node into the list
   uint64_t file_offset;			// Position of node onto the file

};

/* 	Hashband: it is a linked list hash table with a fixed number of entries and a maximum number of nodes.
	When a list is full, further elements will be added to file.
*/

struct hashband_t{

	struct hashnode_t * entry[MAX_HASHBAND_SIZE];	// HashBand entries
	uint64_t size;									// HashBand Number of entries
	FILE * index_file;								// HashBand pointer to the file for its extension
};

/* Minimatch: it contains two minimizers references. The two minimizers matched. T is for the target reference and Q is for the reads query*/
struct minimatch_t{

	uint64_t startT;		// Starting point into the Target Reference genome
	uint8_t strandT;		// Positive strand or Negative strand for the target minimizer
	uint64_t startQ;		//Starting point into the Query read
	uint8_t strandQ;		// Positive strand or Negative strand for the query minimizer
	uint8_t len;			// k-mer length (i.e. MINI_K)
};

/* Sequence Files Abstraction */

struct seqfile_t{

	FILE * file;
	struct hashband_t index;
	uint32_t seqid;
	uint32_t seqlen;
	uint32_t seqstart;
};

/* Multithreading related structure */

#ifdef MULTITHREADING

struct mapper_ctx_t{
	struct seqfile_t RF_local;
	struct seqfile_t * TF_global;
	struct seqfile_t * RF_global;
	pthread_t thread_handler;
	uint64_t reads_num;
	uint8_t id;
	
};

#endif

#endif
