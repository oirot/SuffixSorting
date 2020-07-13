#ifndef CONFROTATION_H_
#define CONFROTATION_H_


unsigned int strToNumber(char *str, size_t beg, size_t len, size_t full);;
void findMinMax(char *str, size_t size, int *min, int *max, size_t window);
//int *createLargeArray(char *str, size_t size, size_t window, size_t *ret_size);
int *createLargeArray(char *str, size_t size, size_t window, size_t *ret_size, int *delta_min, int *delta_max);

double *createLargeArrayD(char *str, size_t size, size_t window, size_t *ret_size, int *supp);
#endif /* CONFROTATION_H_ */
