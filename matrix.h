#ifndef MATRIX_H_
#define MATRIX_H_

typedef struct matrix{
	int size;
	int window_size;
	char *start;
} matrix;

char *get_elem(matrix *mat, int position);
char *set_elem(matrix *mat, int position, char *elem);
matrix *create_matrix(int size, int window_size);
void free_matrix(matrix *m);
void print_matrix(matrix *matrix);
matrix* window_matrix(char *word, int window_size);

#endif /* MATRIX_H_ */
