/*
 * handlers.h
 *
 *  Created on: Nov 23, 2019
 *      Author: vittorio
 */

#ifndef HANDLERS_H_
#define HANDLERS_H_

#define _GNU_SOURCE
#include <sys/resource.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
//#include <divsufsort.h>
#include <math.h>
#include <sys/mman.h>
#include <pthread.h>
#include <unistd.h>
#include <getopt.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


#include "factorizations.h"
#include "utils.h"
#include "matrix.h"
#include "confrotation.h"
#include "progressive_order.h"
#include "radix_sa.h"
#include "bloom.h"

#include "suffix_array.h"
#include "sais.h"


#define BOOL int
#define FALSE 0
#define TRUE 1

typedef struct _icfl_report{
	char *word;
	node_t *icfl;

	size_t word_size;
	size_t window_size;
	size_t no_factors;
	size_t max_fact_len;
	size_t min_fact_len;
	double med_fact_len;
	double avg_fact_len;

	int no_thread;
	int seq_id;

	struct rusage rusage;

	BOOL *finished;

}icfl_report;

void print_icfl_report(icfl_report *report);
void icfl_report_print_stats(void *void_report);
void print_icfl_report(icfl_report *report);
void icfl_report_print_stats(void *void_report);
void icfl_report_clean_up(void * void_report);
void icfl_report_init(void** arg, BOOL* flag_finished, char* word, int no_thread, size_t seq_id, BOOL compare, BOOL use_bloom);
void *icfl_report_calculate(void *arg);


typedef struct {
	struct timeval real_time;
	struct timeval user_time;
	long page_faults;
	long hard_page_faults;
}thread_stat;


typedef struct sa_arg{
	char *word;
	size_t size;
	size_t window_size;
	size_t *sa;
	pthread_t thread;
	int index;
	int status;
	thread_stat thread_stat;
	struct rusage suffix_rusage;
	BOOL compare;
	BOOL use_bloom;
	size_t id;
	size_t no_fact;
	node_t *node;
	BOOL *finished;

}sa_info;

BOOL main_threads(FILE *fp, int out_file, int no_threads, BOOL compare, BOOL use_blooms, void *calculate_funxtion(void *));
BOOL single_thread(FILE *fp, int out_file, BOOL compare, BOOL use_blooms, void*calculate_function(void *));
void *only_memory(void * args);

void copy_info(sa_info *arg);

void *thread_worker(void *args);
void *thread_worker_only_sa(void *args);

void print_stat(void *args);

//void *sais_icfl(void *args);

//BOOL generic_multi_thread(char *file_name, int out_file, int no_threads, BOOL compare, BOOL use_blooms, void init_function(void*, BOOL*, char*, int, int, BOOL, BOOL),
	//	void *calculate_function(void *), void *stat_function(void*), void clean_up(void*));

#endif /* HANDLERS_H_ */
