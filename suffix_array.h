/*
 * suffix_array.h
 *
 *  Created on: Apr 23, 2019
 *      Author: test-vm
 */
#include "factorizations.h"
#include "matrix.h"

#ifndef SUFFIX_ARRAY_H_
#define SUFFIX_ARRAY_H_


int *shift_positions(int *arr, int start, int size);
void add_elem(int suffix, int *sa, int size, matrix* matrix);
char *SA_to_string(int *sa, int size);
int ordine(int x);
int sum10(int x);
int fast_10_pow(int n);
int fast_sum10(int x);
int order_lenght(int ord, int r_extra);
int opt_str_len(int x, int r_extra);
int reminder_lenght(int n, int r_extra, int ord);
int string_length(int n, int r_extra);

void print_SA(int *sa, int size);
void print_SA_ld(size_t *sa, size_t size);
void print_SA_reverse(int *sa, int size);
int initSuffixArray(char *word, char *alphabet, int *sa, int *window_size);
int initSuffixArrayICFL(char *word, char *alphabet, int *sa, int *window_size, node_t **ICFL_node);
int **window_size(char *word, char *alphabet, node_t **ICFL_node);

#endif /* SUFFIX_ARRAY_H_ */
