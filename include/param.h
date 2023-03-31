#ifndef PARAM
#define PARAM

#include <math.h>

/*  ######### CONFIGURATION PARAMETERS #########  */

#define TARGET_FILE_PATH        "data/genomes/"
#define READS_FILE_PATH         "data/reads/"
#define INDEX_FILE_PATH         "data/index/"

#define MODE_COMPRESSED         0x1                 // Comment this to disable compressed mode

#ifdef MODE_COMPRESSED
    #define TARGET_FILE_NAME        "salmonellaEntericaShrt.bp32"
    #define READS_FILE_NAME         "salmonellaEntericaReadsShrt.bp32"
#else
    #define TARGET_FILE_NAME        "salmonellaEntericaShrt.fasta"
    #define READS_FILE_NAME         "salmonellaEntericaReadsShrt.fasta"
#endif

    #define INDEX_FILE_NAME         "salmonellaEntericaReadIndex.hb"

#define BUILD_GENOME_INDEX      0x1                 // Comment this to disable Genome Index Building
#define TAKE_TIME               0x1                 // Comment this to disable Time Measuring 
#define VERBOSE                 0x1                 // Comment this to disable terminal messages

/*  ######### MULTITHREADING PARAMETERS #########  */

//#define MULTITHREADING          0x1                 // Comment this to disable multithreading
#define THREADSNUM              16                  // Number of threads used if MULTITHREADING is defined

/*  ######### TUNING PARAMETERS #########  */

#define MAX_GENOME_CHUNK_SIZE   512//4*1024             // Maximum Chunk Size allowed to store Genome locally
#define MAX_READ_CHUNK_SIZE     1*1024              // Maximum Chunk Size allowed to store Read locally
#define MAX_HASHBAND_SIZE       64*1024             // Maximum HashBand number of entries
#define MAX_HASHNODES           32                  // Maximum number of Linked List Nodes per HashBand entry 

#define MINI_W                  15                  // Window Size for Minimizer Computations ( this has to be the same for both index and seeds ) ( its choice has performance implications )
#define MINI_K                  11                  // k-mer Size for Minimizer Computations ( this has to be the same for both index and seeds ) ( its choice has performance implications )

/*  ######### RUNTIME PARAMETERS #########  */

#define STRAND_PLUS             0
#define STRAND_MINUS            1

#define MAX_32BIT               4294967295
#define MAX_64BIT               9223372036854775807

#define HASHTABLE_PRINT_FULL    0                   // Print locally stored hashband; empty entries are printed as well
#define HASHTABLE_PRINT_DENSE   1                   // Print locally stored hashband; empty entries are ignored

#define IN_LINE                 0                   // Print a sequence with no '\n'
#define NEW_LINE                1                   // Print a sequence with '\n'

#endif
