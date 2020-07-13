#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "matrix.h"
#include "suffix_array.h"
#include "factorizations.h"

int pow_vect[] = {1, 10, 100, 1000, 10000, 100000, 1000000, 10000000, 100000000, 1000000000};


int *shift_positions(int *arr, int start, int size){
	int tmp;
	int old = arr[start];
	int i = start + 1;
	while(i < size && old >= 0){
		tmp = arr[i];
		arr[i] = old;
		old = tmp;
		++i;
	}
	return arr;
}

void add_elem(int suffix, int *sa, int size, matrix* matrix){
	int i;
	for(i = 0; (sa[i] >= 0 && strcmp(get_elem(matrix, sa[i]), get_elem(matrix, suffix)) > 0); ++i);

	if(sa[i] >= 0) shift_positions(sa, i, size);

	sa[i] = suffix;
}



char *SA_to_string(int *sa, int size){

	int string_size = string_length(size, 1) + 3;
	char *SA_string = malloc(sizeof(char) * (string_size + 1));
	SA_string[0] = '[';
	SA_string[1] = '\0';
	char *number_buffer = malloc(sizeof(char) * (ordine(size) + 1));
	for(int i = 0; i <= size; ++i){
		snprintf(number_buffer, 10, "%d", sa[i]);
		strncat(strncat(SA_string, number_buffer, string_size), ",", string_size);
	}

 	free(number_buffer);
	SA_string[string_size - 1] = ']';
	SA_string[string_size] = '\0';
	return SA_string;
}


int ordine(int x){
	int ord = 0;
	while(x > 0){
		x /= 10;
		++ord;
	}
	return ord;
}

int sum10(int x){
	if(x < 0) return 0;
	if(x == 0) return 1;
	return 10 * sum10(x - 1) + 1;
}



int fast_10_pow(int n){
	if(n > 9)
		return -1;
	else
		return pow_vect[n];
}

int fast_sum10(int x){
	int res = 0;
	while(x >= 0){
		res += fast_10_pow(x--);
	}
	return res;
}



/*
int string_length(int n, int r_extra){
	int m_order = ordine(n) - 1;
	int n_s = (m_order * (fast_10_pow(m_order) - fast_10_pow(m_order - 1) + n));
	int a_s = (r_extra * (fast_10_pow(m_order) + n - 1));
	int sum = (sum10(m_order - 1));
	return	n_s + a_s - sum;
}*/

int order_lenght(int ord, int r_extra){
	--ord;
	return (fast_10_pow(ord) * ord) +
			(r_extra * (fast_10_pow(ord) - 1)) -
			sum10(ord - 1);
}


int opt_str_len(int x, int r_extra){
	int ord = ordine(x);
	return ((ord + r_extra) *
			(fast_10_pow(ord) - fast_10_pow(ord - 1) + x + 1))
			- (r_extra + fast_sum10(ord - 1));
}


int reminder_lenght(int n, int r_extra, int ord){
	return ((ord + r_extra) *
			(n - fast_10_pow(ord - 1) + 1));
}

int string_length(int n, int r_extra){
	int ord = ordine(n);
	int ord_len = order_lenght(ord, r_extra);
	int rem_len = reminder_lenght(n, r_extra, ord);
	return ord_len + rem_len;
}


void print_SA(int *sa, int size){
	printf("Begin INSERT SA\n------------------\n");
	for(int i = 0; i < size; ++i){
		printf("%d, ", sa[i]);
	}
	printf("\n------------------\nEND INSERT SA\n");
}

void print_SA_ld(size_t *sa, size_t size){

	printf("Begin INSERT SA\n------------------\n");
	for(size_t i = 0; i < size; ++i){
		printf("%ld, ", sa[i]);
	}
	printf("\n------------------\nEND INSERT SA\n");
}


void print_SA_reverse(int *sa, int size){
	printf("Begin INSERT SA\n------------------\n");
	for(int i = size - 1; i >= 0; --i){
		printf("%d, ", sa[i]);
	}
	printf("\n------------------\nEND INSERT SA\n");
}


//create the suffix array of word using the induced order by alphabet
//sa needs to be preallocated and it's in reverse order (>)
//w_size can be NULL otherwise there will be the window size calculated for the string word
//returns int => 0 if ok < 0 otherwise
int initSuffixArray(char *word, char *alphabet, int *sa, int *w_size){
	return initSuffixArrayICFL(word, alphabet, sa, w_size, NULL);
}

int initSuffixArrayICFL(char *word, char *alphabet, int *sa, int *w_size, node_t **ICFL_node){
	int **window_size = malloc(sizeof(int*));
	if(window_size == NULL) return -1;

	//Calculating ICFL(word)
	node_t *node = ICFL_recursive_for_alphabet(word, alphabet, window_size);
	if(ICFL_node == NULL){
		free_list(node);
	}else{
		*ICFL_node = node;
	}
	if(node == NULL) return -1;



	//calculating SA
	matrix *m = window_matrix(word, **window_size);
	if(m == NULL) return -1;

	int word_len = strlen(word);
	for(int i = 0; i < word_len; ++i){
		sa[i] = -1;
	}
	for(int i = 0; i < word_len; ++i){
		add_elem(i, sa, word_len, m);
	}
	free_matrix(m);

	if(w_size != NULL) *w_size = **window_size;

	free(*window_size);
	free(window_size);

	return 1;
}


int **window_size(char *word, char *alphabet, node_t **ICFL_node){
	int **window_size = malloc(sizeof(int*));
	if(window_size == NULL) return -1;

	//Calculating ICFL(word)
	node_t *node = ICFL_recursive_for_alphabet(word, alphabet, window_size);
	if(ICFL_node == NULL){
		free_list(node);
	}else{
		*ICFL_node = node;
	}
	if(node == NULL) return -1;

	return window_size;

}


