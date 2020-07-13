/*
 * radix_sa.h
 *
 *  Created on: Sep 15, 2019
 *      Author: vittorio
 */

#ifndef RADIX_SA_H_
#define RADIX_SA_H_


typedef struct {
	long counting;
	long offset;
	long hits;
} offset_hits;


struct read_s{
	char *word;
	size_t size;
	size_t window_size;
	offset_hits *offsets;
	struct bloom *blooms;
	size_t *suffix_array;
	int use_blooms;
};

struct read_s *compute_suffix_array_read(struct read_s *read);
size_t *compute_suffix_array_no_bloom(char *word, size_t word_size, size_t window_size);
size_t *compute_suffix_array(char *word, size_t word_size, size_t window_size, int use_blooms);

#endif /* RADIX_SA_H_ */
