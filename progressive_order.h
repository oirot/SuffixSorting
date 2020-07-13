#ifndef PREGESSIVE_H_
#define PREGESSIVE_H_
#include "factorizations.h"

int calculate_position(int delta, int delta_m, int delta_M, int no_positions);
int partial_oreder(int *values, int *prev_ordered, node_t *current_fact, int old_boundary, int *min, int *max, int *no_collisions, int is_first_lap);
int *create_values_array(char *str, int str_len, int window_size);
int *progressive_order(int *values, node_t *icfl_fact, int word_len, int *no_collisions);
int *progressive_sa(char *word, char *alphabet, int *no_collisions, int *values, int is_porg);

#endif /*PREGESSIVE_H:*/
