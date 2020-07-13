#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include "matrix.h"

char *get_elem(matrix *mat, int position){
	if(position >= mat->size) return NULL;
	return mat->start + (position * (mat->window_size + 1));
}

char *set_elem(matrix *mat, int position, char *elem){
	if(position >= mat->size) return NULL;
	*(mat->start + ((position + 1)* (mat ->window_size + 1)) - 1) = '\0';
	return strncpy(mat->start + (position * (mat->window_size + 1)), elem, mat->window_size);
}

matrix *create_matrix(int size, int window_size){
	matrix *m = malloc(sizeof(matrix));
	if(m == NULL){
		perror("Matrix creation struct");
		return NULL;
	}
	m->size = size;
	m->window_size = window_size;
	m->start = malloc(sizeof(char) * size * (window_size + 1));
	if(m->start == NULL){
		perror("Matrix creation string space");
		free(m);
		return NULL;
	}
	return m;
}

void free_matrix(matrix *m){
	free(m->start);
	free(m);
}


//Create a matrix m[word_len][window_size] where m[i] = word[i...i+window_size]
matrix* window_matrix(char *word, int window_size){
	size_t word_size = strlen(word);
	matrix *mat = create_matrix(strlen(word), window_size);
	if(mat == NULL) return NULL;

	for(int i = 0; i < word_size; i++){
		set_elem(mat, i, word + i);
	}

	return mat;
}

void print_matrix(matrix *matrix){
	printf("Start matrix\n----------------------------------\n");
	for(int i = 0; i < matrix->size; i++){
		printf("%s,\n", get_elem(matrix, i));
	}
	printf("\n----------------------------------\nEnd matrix\n");
}
