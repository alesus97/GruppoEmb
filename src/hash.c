#include "hash.h"

struct hashnode_t *head = NULL;
struct hashnode_t *current = NULL;

void initHashTable( struct hashband_t * ht ){

	uint8_t k = MINI_K;
	uint8_t w = MINI_W;

	char index_path[256] = INDEX_FILE_PATH;

	ht->index_file = fopen(strcat(index_path,INDEX_FILE_NAME),"wb+");
	ht->size = MAX_HASHBAND_SIZE;

	if( ht->index_file != NULL ){

		for( int i = 0; i < MAX_HASHBAND_SIZE; i++ ){

			ht->entry[i] = (struct hashnode_t*) malloc(sizeof(struct hashnode_t));
		
			ht->entry[i]->next = NULL;
			ht->entry[i]->file_offset = 0x0;
			ht->entry[i]->list_offset = 0x0;
			ht->entry[i]->hashmizer.hash = 0;
			ht->entry[i]->hashmizer.len = 0;
			ht->entry[i]->hashmizer.strand = 0;
			ht->entry[i]->hashmizer.start = 0;
		}
	}
	else{
		printf("\033[31m[initHashTable] : Error in opening hash band file\n\033[37m");
	}
}

void addToHashTable( struct hashband_t * ht, const struct hashmizer_t * hm ){
	
	insertFirstNode( ht, hm, hm->hash % ht->size );
	
}

void queryHashTable( const struct hashband_t * ht, const struct hashmizer_t * hm, cvector_vector_type(struct minimatch_t) * matches ){
	
	uint32_t index = hm->hash % ht->size;
	uint64_t current_file_offset = ( ht->entry[index]->next != NULL ) ? ht->entry[index]->next->file_offset : 0;

	struct hashnode_t* current = ht->entry[index];
	struct hashmizer_t mini;
	struct minimatch_t minmatch;

	//navigate through in memory list
	while( current->next != NULL ) {
		if(current->hashmizer.hash == hm->hash) {

			//Save matching minimizers information
			minmatch.len = hm->len;
			minmatch.startQ = hm->start;
			minmatch.startT = current->hashmizer.start;
			minmatch.strandQ = hm->strand;
			minmatch.strandT = current->hashmizer.strand;
			cvector_push_back(*matches, minmatch);

			//printf("[4]matches is %16x\n", matches);
		}
		current = current->next;
	}

	//navigate through in file list

	while( current_file_offset != 0 ){
		fseek(ht->index_file, current_file_offset + 8, SEEK_SET); 			//Skip the list ID
		fread(&mini.hash, sizeof(uint64_t), 1, ht->index_file); 
		fread(&mini.start, sizeof(uint64_t), 1, ht->index_file);
		fread(&mini.strand, sizeof(uint8_t), 1, ht->index_file);		
		fread(&current_file_offset, sizeof(uint64_t), 1, ht->index_file); 	//Jump to the next minimizer on file

		if( mini.hash == hm->hash ){

			//Save matching minimizers information
			minmatch.len = hm->len;
			minmatch.startQ = hm->start;
			minmatch.startT = mini.start;
			minmatch.strandQ = hm->strand;
			minmatch.strandT = mini.strand;
			cvector_push_back(*matches, minmatch);
			//printf("[4]matches is %16x\n", matches);
			
		}
	}

	//if( cvector_size(matches) ) printf("found <%ld> matches on list [%d]\n", cvector_size(matches), index);

}

void freeHashTable( struct hashband_t * ht ){

	uint32_t len;
	struct hashnode_t * v;
	struct hashnode_t * u;
	
	for( int i = 0; i < MAX_HASHBAND_SIZE; i++ ){
	
		v = ht->entry[i];
		while( v != NULL ) {
			u = v;
			v = v->next;
			free(u);
		}
	}

}

#ifndef MODE_COMPRESSED

uint64_t hashSequence( const char * read, const uint32_t seqlen ){			//Supporting up to 31-mers

	uint64_t hash = 0;
	uint32_t val = 0;

	for( int i = 0; i < seqlen; i++ ){
	
		switch(read[i]){
			case 'A': val = 0; break;
			case 'T': val = 3; break;
			case 'C': val = 1; break;
			case 'G': val = 2; break;
			
		}
		hash += (val << 2*(seqlen-i-1));
	}
	
	return hashInv(hash);
}

#else

uint64_t hashSequence_bp32( const bp32_t read, const uint32_t seqlen ){

	uint64_t hash = 0;
	uint32_t val = 0;
	uint8_t tmp = 0;
	
	for( int i = 0; i < seqlen; i++ ){

		tmp = (read >> 2*(i)) & 0x3;
		hash += (tmp << 2*(seqlen-i-1));
		
	}

	return hashInv(hash);
}

#endif

uint64_t hashInv( uint64_t key ){

	key = (~key) + (key << 21); // key = (key << 21) - key - 1;
	key = key ^ (key >> 24);
	key = (key + (key << 3)) + (key << 8); // key * 265
	key = key ^ (key >> 14);
	key = (key + (key << 2)) + (key << 4); // key * 21
	key = key ^ (key >> 28);
	key = key + (key << 31);
	return key;
}

//List functions

void insertFirstNode( struct hashband_t * ht, const struct hashmizer_t * hm, const uint64_t index ){

	uint64_t cur_file_offset = 0x0;
	uint64_t prv_file_offset = 0x0;
	uint64_t last_item_offset;

	// End of file is where we will add our new item

	fseek(ht->index_file, 0, SEEK_END);
	last_item_offset = ftell(ht->index_file);
	fwrite(&index, sizeof(uint64_t), 1, ht->index_file); 
	fwrite(&hm->hash, sizeof(uint64_t), 1, ht->index_file); 
	fwrite(&hm->start, sizeof(uint64_t), 1, ht->index_file); 
	fwrite(&hm->strand, sizeof(uint8_t), 1, ht->index_file); 
	fwrite(&cur_file_offset, sizeof(uint64_t), 1, ht->index_file); 

	// If this wasn't the first item, find the previous last item in the list and replace the offset

	if( ht->entry[index]->list_offset != 0){
		fseek(ht->index_file, ht->entry[index]->next->file_offset + 25, SEEK_SET);
		fread(&cur_file_offset, sizeof(uint64_t), 1, ht->index_file); 

		prv_file_offset = ht->entry[index]->next->file_offset;

		while( cur_file_offset != 0){
			prv_file_offset = cur_file_offset;
			fseek(ht->index_file, cur_file_offset + 25, SEEK_SET);
			fread(&cur_file_offset, sizeof(uint64_t), 1, ht->index_file); 
		}

		fseek(ht->index_file, prv_file_offset + 25, SEEK_SET);
		fwrite(&last_item_offset, sizeof(uint64_t), 1, ht->index_file);

	}
	
	// Add the node to the Memory Cached copy ( if possible )

	if( ht->entry[index]->list_offset < MAX_HASHNODES ){
		//create a link
		struct hashnode_t *link = (struct hashnode_t*) malloc(sizeof(struct hashnode_t));
		
		ht->entry[index]->file_offset = last_item_offset;
		link->hashmizer = *hm;
		link->list_offset = ht->entry[index]->list_offset + 1;
		link->file_offset = 0x0;
		

		//point it to old first node
		link->next = ht->entry[index];
		
		//point first to new first node
		ht->entry[index] = link;
	}

}


struct hashnode_t* findNode( const struct hashband_t * ht, const struct hashmizer_t * hm, const uint32_t index, uint64_t key ){

	// Check if the node is Memory Cached, otherwise check into the indexFile

	//start from the first link
	struct hashnode_t* current = ht->entry[index];
	uint8_t onfile_flag = 0;
	
	if( current->list_offset >= 0 ){
		
		onfile_flag = ( !current->list_offset ) ? 1 : 0 ;
		//if list is empty
		if( getListLength( ht, index ) == 0 )  return NULL;

		//navigate through list
		while(current->hashmizer.hash != key) {
			
				//if it is last node
				if(current->next == NULL) {
					break;//return NULL; 
				} else {
				//go to next link
					current = current->next;
				}
		}      
	} 
	
	if( onfile_flag ) {
		printf("[findNode] Current node might be on file\n");
	}
	else printf("[findNode] Node found\n");
	
	//if data found, return the current Link
	return current;
}
