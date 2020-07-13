#include <math.h>
#include "matrix.h"
#include <stdlib.h>
#include "confrotation.h"


unsigned int strToNumber(char *str, size_t beg, size_t len, size_t full){
	int val = 0;
	int mult = 0;
	int curr;
	if(beg + len > full){
		mult = beg + len - full;
		len = full - beg;
	}
	len += beg;
	for(int i = beg; i < len; ++i){
		switch(str[i]){
		case 'N':
			curr = 0;
			break;
		case 'A':
			curr = 1;
			break;
		case 'C':
			curr = 2;
			break;
		case 'G':
			curr = 3;
			break;
		case 'T':
			curr = 4;
			break;
		}
		val = val * 5 + curr;
	}
	if(mult != 0) val = val *  pow(5, mult);
	return val;
}

void findMinMax(char *str, size_t size, int *min, int *max, size_t window){
	*min = strToNumber(str, 0, window, size);
	*max = *min;
	int curr_elem;
	for(size_t i = 1; i < size; i++){
		curr_elem = strToNumber(str, i, window, size);
		if(*min > curr_elem) *min = curr_elem;
		else if(*max < curr_elem) *max = curr_elem;
	}
}


int *createLargeArray(char *str, size_t size, size_t window, size_t *ret_size, int *delta_min, int *delta_max){
	int min, max;
	findMinMax(str, size, &min, &max, window);
	/*size_t deltaSize = (max - min + 1);
	*ret_size = deltaSize;
	int *deltaVect = calloc(deltaSize, sizeof(int));
	 */
	size_t deltaSize = max + 1;
	*ret_size = deltaSize;
	int *deltaVect = calloc(deltaSize, sizeof(int));

	int last_before = strToNumber(str, 0, window, size);
	for(size_t i = 0; i < size; i++){

		//deltaVect[strToNumber(str, i, window, size) - min] = i + 1;
		int numb = strToNumber(str, i, window, size);
		if(last_before == max || (last_before != max && last_before < numb)) last_before = numb;
		//deltaVect[numb] = (int) ceil( (max - numb) * size / max);
		deltaVect[numb] = i + 1;
		numb = 0; //TODO SOLO PER DEBUG
	}
	*delta_min = max - last_before;
	*delta_max = max - min;

	return deltaVect;
}

double *createLargeArrayD(char *str, size_t size, size_t window, size_t *ret_size, int *supp){
	int min, max;
	findMinMax(str, size, &min, &max, window);
	/*size_t deltaSize = (max - min + 1);
	*ret_size = deltaSize;
	int *deltaVect = calloc(deltaSize, sizeof(int));
	 */
	size_t deltaSize = max + 1;
	*ret_size = deltaSize;
	double *deltaVect = calloc(deltaSize, sizeof(int));

	for(size_t i = 0; i < size; i++){
//		deltaVect[strToNumber(str, i, window, size) - min] = i + 1;
		int numb = strToNumber(str, i, window, size);
	//deltaVect[numb] = (max - numb) * size / (double) max;
		deltaVect[strToNumber(str, i, window, size)] = i + 1;

		supp[i] = numb;
		numb = 0; //TODO SOLO PER DEBUG
	}

	return deltaVect;
}
