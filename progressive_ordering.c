#include <string.h>
#include "suffix_array.h"
#include <math.h>
#include <stdlib.h>
#include "confrotation.h"
#include "progressive_order.h"
int calculate_position(int delta, int delta_m, int delta_M, int no_positions){
	return  floor( (delta_M - delta) * (no_positions - 1)  / (delta_M - delta_m) );
}

int partial_oreder(int *values, int *prev_ordered, node_t *current_fact, int old_boundary, int *min, int *max, int *no_collisions, int isFirstLap){
	int fact_len = strlen(current_fact->factor), boundary = old_boundary + fact_len;
	int stimed_position;

	//i can do that bc we jump old_boundary and we carry over this result
	int real_max;
	if(isFirstLap){
		*min = real_max = values[0];
		*max = values[1];

		for(int i = 1; i < boundary; ++i){
			if(real_max < values[i]){
				*max = real_max;
				real_max = values[i];
			}
		}
	}else if(current_fact->next != NULL){
		real_max = values[boundary];
	}else{
		real_max = values[old_boundary];
	}

	//findig min max(-) for the current boundary
	for(int i = old_boundary + 1; i < boundary; ++i){
		if(*min > values[i]) *min = values[i];
		 // not needed "values[i] - values[old_boundary]" bc we start at old_boundary + 1
		else if(values[i] != real_max && *max < values[i]) *max = values[i];
	}


	//printf("\n-\nbnd:%d\n", boundary);

	int delta_M, delta_m, last;
	last = values[old_boundary];
	delta_M = real_max - *min;
	delta_m = real_max - *max;

	int i;
	if(isFirstLap) i = 0;
	else i = 1;
	for(; i < fact_len; ++i){
		if(isFirstLap && real_max == values[old_boundary + i]) continue;
		stimed_position = calculate_position(real_max - values[i + old_boundary],  real_max - *max, real_max - *min, boundary);
		//TODO debug
	//	printf("%d: %d\n", i, stimed_position);
		if(stimed_position >= boundary){
			stimed_position = stimed_position + 1 - 1;
		}
		if(boundary >= 12){
			boundary = boundary +1 -1;
		}
		if(prev_ordered[stimed_position] >= 0){
			if(values[prev_ordered[stimed_position]] < values[old_boundary + i]){
				stimed_position += 1;
			}
			shift_positions(prev_ordered, stimed_position, boundary);
			*no_collisions += 1;
		}
		prev_ordered[stimed_position] = i + old_boundary;
	}
	//bc carry over
	*max = real_max;
	prev_ordered[boundary - 1] = old_boundary;

	return boundary;
}


int *create_values_array(char *str, int str_len, int window_size){
	int *values = malloc(sizeof(int) * str_len);
	if(values == NULL) return NULL;
	for(int i = 0; i < str_len; ++i){
		values[i] = strToNumber(str, i, window_size, str_len);
	}
	return values;
}


int *progressive_order(int *values, node_t *icfl_fact, int word_len, int *no_collisions){
	node_t *current_factor = icfl_fact;
	int *ordered = malloc(sizeof(int) * word_len);
	for(int i = 0; i < word_len; ++i){
		ordered[i] = -1;
	}
	int current_boundary = 0, min, max;
	*no_collisions = 0;
	min = max = -1;
	int isFirstLap = 1;
	while(current_factor != NULL){
		current_boundary = partial_oreder(values, ordered, current_factor, current_boundary, &min, &max, no_collisions, isFirstLap);
		isFirstLap = 0;
		current_factor = current_factor->next;
	}

	return ordered;
}

int *progressive_sa(char *word, char *alphabet, int *no_collisions, int *values, int is_prog){
	int **window_size = malloc(sizeof(int*));
	node_t *root = NULL;
	if(is_prog){
		root = ICFL_recursive_for_alphabet(word, alphabet, window_size);
	}else{
		root = malloc(sizeof(node_t));
		root->factor = malloc(sizeof(char) * (strlen(word) + 1));
		strcpy(root->factor, word);
		root->next = NULL;
	}

	values = create_values_array(word, strlen(word), **window_size);
	int *sa = progressive_order(values, root, strlen(word), no_collisions);
	free(*window_size);
	free(window_size);
	free(values);
	free_list(root);
	return sa;
}
