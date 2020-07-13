#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <ctype.h>
#include "radix_sa.h"
#include "bloom.h"


#define CHAR_LOWER_BOUND 64
#define ALPHABET_SIZE 6
#define ERR_OUT_OF_MEMORY "Out of memory"
#define FALSE_POSITIVE_PROB 0.01
/*
char *_word;
size_t _size;
size_t _window_size;
size_t *_suffix_array;*/
char alphabet[] = {0,1,2,3,4,5};


char symbolAt(size_t index, struct read_s *read){
	//if index is out of the bound of the word that we return 0 that is $
	if(index >= read->size)
		return 0;

	//if we are not using the bloom filters on this read we just check in memory
	if(!read->use_blooms){

		return read->word[index];
	}else{

		int bloom_results[ALPHABET_SIZE - 1] = {0};

		int checked = -1;
		// we query each bloom filter
		for(int i = 0; i < ALPHABET_SIZE - 1; ++i){
			if((bloom_results[i] = bloom_check(&read->blooms[i], &index, sizeof(size_t))) == 1 &&
					checked != -1){
				//if more than one filter report a positve result for the index
				//we check in memory
				return read->word[i];

			}else if(bloom_results[i] == 1){
				checked = i;
			}
		}
		//else if only one bloom filter report the index as a positve we are sure that it's the correct one
		return checked + 1;

	}

}



struct bloom *init_blooms(char *word){
	struct bloom *blooms = malloc(sizeof(struct bloom) * (ALPHABET_SIZE - 1));
	if(blooms == NULL){
		perror(ERR_OUT_OF_MEMORY);
		return NULL;
	}
	for(int i = 0; i < ALPHABET_SIZE - 1; ++i){
		if(bloom_init(&blooms[i], strlen(word) + 1, FALSE_POSITIVE_PROB) != 0){
			free(blooms);
			perror("Bloom error");
			return NULL;
		}
	}

	for(size_t i = 0; i < strlen(word); ++i){
		bloom_add(&blooms[word[i] - 1], &i, sizeof(size_t));
	}

	return blooms;
}

void radix_step(size_t i, struct read_s *read, size_t *old_suffix_array, size_t *new_suffix_array){
	char symbol;
	size_t index;
	offset_hits *offset;
	size_t next_index;
	for(size_t j = 0; j < read->size; ++j){
		index = old_suffix_array[j];
		//symbol = symbolAt(index + read->window_size - i, read);
		next_index = index + read->window_size - i;
		symbol = next_index < read->size ? read->word[index + read->window_size - i] : 0;
		offset = read->offsets + (size_t) symbol;
		new_suffix_array[offset->counting + offset->offset + (offset->hits++)] = index;
	}
}


void updateOffsets(struct read_s *read, size_t index_to_add){
	//char newSymbol = symbolAt(read->window_size - index_to_add, read);
	char newSymbol = read->word[read->window_size - index_to_add];
	char i = 0;
	read->offsets[0].hits = 0;
	do {
		++i;
		if(i < ALPHABET_SIZE){
			--(read->offsets[i].offset);
			read->offsets[i].hits = 0;
		}
	} while (alphabet[i] != newSymbol);

	for(; i < ALPHABET_SIZE; ++i){
		read->offsets[i].hits = 0;
	}
}

void down_convert(char *word, size_t size){
	for(size_t i = 0; i < size; ++i){
		char converted;
		switch(tolower(word[i])){
		case 'a':
			converted = 1;
			break;
		case 'c':
			converted = 2;
			break;
		case 'g':
			converted = 3;
			break;
		case 't':
			converted = 4;
			break;
		case 'n':
			converted = 5;
			break;
		default:
			converted = 0;
			break;
		}

		word[i] = converted;
	}
}

void init(struct read_s *read){
	char symbol;
	char *word = read->word;
	size_t size = read->size;
	offset_hits *offsets = read->offsets;

	for(int i = 0; i <= read->window_size - 1; ++i){
		//this is a sequencial read no need to use bloom filters
		symbol = read->word[i];
		//if symbol is $ no need to set it to -1
		if(symbol < ALPHABET_SIZE - 1){

			offsets[symbol + 1].offset = offsets[symbol + 1].offset - 1;

		}
	}
	for(size_t i = 0; i < size; ++i){
		symbol = read->word[i];
		offsets[symbol].counting = offsets[symbol].counting + 1;
	}

	size_t count = read->window_size;
	size_t current;
	offsets[0].counting = 0;
	for(int i = 1; i < ALPHABET_SIZE; ++i){
		current = offsets[i].counting;
		offsets[i].counting = count + offsets[i - 1].counting;
		offsets[i].offset = offsets[i].offset + offsets[i - 1].offset;
		count = current;
	}

}

size_t *init_suffix_array(size_t *suffix_array, size_t size){
	for(size_t i = 0; i < size; ++i){
		suffix_array[i] = i;
	}
	return suffix_array;
}


struct read_s *init_read(struct read_s *read){
	//convert char  A C G T N to 1 2 3 4 5
	down_convert(read->word, read->size);

	if(read->use_blooms){
		if((read->blooms = init_blooms(read->word)) == NULL){
			read->use_blooms = 0;
		}
	}
	else{
		read->blooms = NULL;
	}

	read->offsets = calloc(ALPHABET_SIZE, sizeof(offset_hits));
	//init the offsets
	init(read);

	return read;

}

struct read_s *init_read_old(char *word, size_t size, size_t window_size, int use_blooms){
	//convert char  A C G T N to 1 2 3 4 5
	down_convert(word, size);

	struct read_s *read = malloc(sizeof(struct read_s));
	read->word = word;
	read->size = size;
	read->window_size = window_size;

	read->use_blooms = use_blooms;
	if(use_blooms){
		read->blooms = init_blooms(read->word);
	}

	read->offsets = calloc(ALPHABET_SIZE, sizeof(offset_hits));
	//init the offsets
	init(read);

	return read;

}



/*
size_t *compute_suffix_array(char *word, size_t word_size, size_t window_size, int use_blooms){
	//if lydon word w_size = size + 1, we need to compute on w_size = size
	if(window_size > word_size) window_size = word_size;

	struct read_s *read = init_read_old(word, word_size, window_size, use_blooms);


	size_t *old_suffix_array = malloc(sizeof(size_t) * word_size);
	size_t *new_suffix_array = malloc(sizeof(size_t) * word_size);
	size_t *swap_support;
	init_suffix_array(new_suffix_array, word_size);


	//offset_hits offsets[alphabet_size];
	/*for(int i = 0; i < alphabet_size; ++i){
		offsets[i].counting = 0;
		offsets[i].hits = 0;
		offsets[i].offset = 0;
	}*//*
	for(size_t j = 1; j <= window_size; ++j){
		//swap gli array
		swap_support = new_suffix_array;
		new_suffix_array = old_suffix_array;
		old_suffix_array = swap_support;

		updateOffsets(read, j);
		radix_step(j, read, old_suffix_array, new_suffix_array);
	}

	free(old_suffix_array);
	clean_up(read);

	return new_suffix_array;
}*/


struct read_s *compute_suffix_array_read(struct read_s *read){
	if(read->window_size > read->size) read->window_size = read->size;

	init_read(read);

	size_t *old_suffix_array = malloc(sizeof(size_t) * read->size);
	if(old_suffix_array == NULL){
		free(read->offsets);
		perror(ERR_OUT_OF_MEMORY);
		return NULL;
	}

	size_t *new_suffix_array = malloc(sizeof(size_t) * read->size);
	if(new_suffix_array == NULL){
		free(read->offsets);
		free(old_suffix_array);
		perror(ERR_OUT_OF_MEMORY);
		return NULL;
	}

	size_t *swap_support;
	init_suffix_array(new_suffix_array, read->size);

	for(size_t j = 1; j <= read->window_size; ++j){
		//swap
		swap_support = new_suffix_array;
		new_suffix_array = old_suffix_array;
		old_suffix_array = swap_support;

		updateOffsets(read, j);
		radix_step(j, read, old_suffix_array, new_suffix_array);
	}

	free(old_suffix_array);
	free(read->offsets);

	read->suffix_array = new_suffix_array;

	return read;
}


size_t *compute_suffix_array(char *word, size_t word_size, size_t window_size, int use_blooms){
	struct read_s read;

	read.word = word;
	read.size = word_size;
	read.window_size = window_size;
	read.use_blooms = use_blooms;

	compute_suffix_array_read(&read);

	return read.suffix_array;
}

size_t *compute_suffix_array_no_bloom(char *word, size_t word_size, size_t window_size){
	return compute_suffix_array(word, word_size, window_size, 0);
}


