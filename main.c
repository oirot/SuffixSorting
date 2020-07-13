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
#include "handlers.h"

#define BUFF_DIM 1024

#define NPAR 2

#define BOOL int
#define TRUE 1
#define FALSE 0


char *last_map_end = NULL;

char *read_fasta_map(char *map){

	if(last_map_end == NULL) last_map_end = map;
	else if(*last_map_end == '\0') return NULL;

	char *start = last_map_end;

	while(*start != '\n'){
		++start;
	}
	++start;

	size_t shift = 0;
	size_t i = 0;

	while(*(start + i) != '>' && *(start + i) != '\0'){
		if(*(start + i)== '\n'){
			++shift;
		}else{
			start[i - shift] = start[i];
		}
		++i;
	}

	last_map_end = start + i;

	start[i - shift + 1] = '\0';
	//_word_fasta[strlen(_word_fasta) - shift] = '\0';
	return start;
}


void no_merge_print_stat(void *_args){
	sa_info *args = (sa_info*) _args;

	printf("%d,%d,%ld,%ld,%ld,%ld.%06ld,%ld.%06ld,%ld,%ld\n",
			args->status,
			args->id,
			args->size,
			args->no_fact,
			args->window_size,
			args->suffix_rusage.ru_utime.tv_sec,
			args->suffix_rusage.ru_utime.tv_usec,
			args->suffix_rusage.ru_stime.tv_sec,
			args->suffix_rusage.ru_stime.tv_usec,
			args->suffix_rusage.ru_minflt,
			args->suffix_rusage.ru_majflt);
	fflush(stdout);
}

size_t no_merge_icfl_sais(char *word, size_t word_len, int *w_size, BOOL use_bloom, node_t **icfl){
	//TODO change internal logic to use size_t
	int **w_size_int;

	w_size_int =  window_size(word, "ACGTN", icfl);
	node_t *iter_node = *icfl;
	int i = 0;
	size_t *sa = malloc(sizeof(int) * word_len);

	while(iter_node != NULL){
		sais(iter_node->factor, sa, strlen(iter_node->factor));
		iter_node = iter_node->next;
		i++;
	}


	//size_t *sa = compute_suffix_array(word, word_len, **w_size_int, use_bloom);

	*w_size = **w_size_int;
	free(*w_size_int);
	free(w_size_int);
	free(sa);
	//print_SA_ld(sa, word_len);
	return i;
}


void *sais_icfl(void *args){
	sa_info *struct_args = ((sa_info*)args);
	int window_size = NULL;

	if(strcmp(struct_args->word,"") == 0 ||
			struct_args->size <= 0){

		fprintf(stderr, "Error on word: %s\tsize: %ld\tw_size: %ld\n", struct_args->word, struct_args->size, struct_args->window_size);
		fflush(stderr);

		struct_args->status = 0;
		*struct_args->finished = 1;
		return args;
	}

	struct_args->no_fact = no_merge_icfl_sais(struct_args->word, struct_args->size, &window_size, struct_args->use_bloom, &struct_args->node);

	struct_args->window_size = window_size;

	int res = getrusage(RUSAGE_THREAD, &struct_args->suffix_rusage); // @suppress("Symbol is not resolved")


	if(struct_args->compare){
		compare_to_libs(struct_args);
	}else{
		struct_args->status = 2;
	}


	struct_args->status = 1;
	*struct_args->finished = 1;
	return args;
}



void init_suffix_arg(sa_info *args, char *word, int thread_no, BOOL compare, BOOL use_blooms, size_t id){
	args->size = strlen(word);

	char *arg_word = malloc(sizeof(char) * args->size + 1);
	strncpy(arg_word, word, args->size + 1);

	args->word = arg_word;
	args->sa = NULL;
	args->index = thread_no;
	//todo controlla se ci vuole anche qua *finished
	args->finished = 0;
	args->compare = compare;
	args->window_size = -1;
	args->use_bloom = use_blooms;
	args->id = id;
}

void init_sa_no_merge(void **_args, BOOL *finished, char *word, int thread_no, size_t id, BOOL compare, BOOL use_blooms){
	sa_info *args = malloc(sizeof(sa_info));


	args->size = strlen(word);

	args->word = word;
	args->sa = NULL;
	args->index = thread_no;

	args->compare = 0;
	args->window_size = 0;
	args->use_bloom = 0;
	args->id = id;

	args->finished = finished;
	*finished = 0;


	*_args = args;
}

void no_clean_up(void *args){
	sa_info *info = (sa_info *) args;
	free_list(info->node);

	free((sa_info *) args);
}


/*
void * sa_no_merge(void *args){
	sa_info *struct_args = ((sa_info*)args);
	int window_size = NUL	sa_info *info = (sa_info *) args);
	free(info->)L;

	if(strcmp(struct_args->word,"") == 0 ||
			struct_args->size <= 0 ||
			struct_args->window_size <= 0){

		fprintf(stderr, "Error on word: %s\tsize: %ld\tw_size: %ld\n", struct_args->word, struct_args->size, struct_args->window_size);
		fflush(stderr);

		struct_args->status = 0;
		struct_args->finished = 1;
		return args;
	}

	struct_args->no_fact = radix_suffix_array(struct_args->word, struct_args->size, &window_size, struct_args->use_bloom);

	struct_args->window_size = window_size;

	int res = getrusage(RUSAGE_THREAD, &struct_args->suffix_rusage); // @suppress("Symbol is not resolved")


	if(struct_args->compare){
		compare_to_libs(struct_args);
	}else{
		struct_args->status = 2;
	}

	struct_args->finished = 1;
	return args;
}
*/

typedef struct _sais_info{
	char *word;
	int length;
	int thread_no;
	int id;
	struct rusage rusage;
	int *sa;
	BOOL *finished;
}sais_info;


void sais_init(void **_args, BOOL *finished, char *word, int thread_no, size_t id, BOOL use_blooms, BOOL compare){
	sais_info *info = malloc(sizeof(sais_info));

	info->thread_no = thread_no;
	info->id = id;
	info->length = strlen(word);
	info->word = word;
	info->sa = NULL;
	info->finished = finished;
	*finished = 0;

	*_args = info;
}

void *sais_compute(void *_info){
	sais_info *info = (sais_info *) _info;

	info->sa = malloc(sizeof(int) * info->length);

	sais(info->word,info->sa, info->length);

	getrusage(RUSAGE_THREAD, &info->rusage );// @suppress("Symbol is not resolved")

	*info->finished = 1;
}

void sais_clean(void *_info){
	sais_info *info = (sais_info *) _info;

	free(info->sa);
	free(info);
}

void sais_stat(void *_info){
	sais_info *info = (sais_info *) _info;

	printf("%d,%d,%d,%ld.%06ld,%ld.%06ld,%ld,%ld,%ld\n",
			info->id,
			info->thread_no,
			info->length,
			info->rusage.ru_utime.tv_sec,
			info->rusage.ru_utime.tv_usec,
			info->rusage.ru_stime.tv_sec,
			info->rusage.ru_stime.tv_usec,
			info->rusage.ru_majflt,
			info->rusage.ru_majflt,
			info->rusage.ru_maxrss);

}


BOOL generic_multi_thread(char *file_name, int out_file, int no_threads, BOOL compare, BOOL use_blooms, void init_function(void**, BOOL*, char*, int, size_t, BOOL, BOOL),
		void *calculate_function(void *), void stat_function(void*), void clean_up(void*)){

	char *fasta_line;
	void *args[no_threads];
	pthread_t threads[no_threads];
	BOOL finished_flags[no_threads];
	int t_res;
	int current_threads;
	size_t id = 0;

	struct stat stats;
	if(stat(file_name, &stats) < 0){
		perror("stat");
		return FALSE;
	}

	int fd = 0;
	if( (fd = open(file_name, O_RDONLY)) == NULL){
		perror("stat");
		return FALSE;
	}

	void *map;
	fasta_line = map = mmap(NULL, stats.st_size, PROT_READ | PROT_WRITE, MAP_PRIVATE, fd, 0);
	if(fasta_line == MAP_FAILED){
		perror("mmap");
		return FALSE;
	}



	for(current_threads = 0; current_threads < no_threads; ++current_threads){

		if((fasta_line = read_fasta_map(fasta_line)) != NULL){
			init_function(&(args[current_threads]), &finished_flags[current_threads], fasta_line, current_threads, id++, compare, use_blooms);
			t_res = pthread_create(&(threads[current_threads]), NULL, calculate_function, args[current_threads]);

			if(t_res < 0) perror("ngul o riavl");

		}else{
			break;
		}
	}


	int no_started = current_threads;
	while(no_started){
		for(int i = 0; i < current_threads; ++i){
			if(finished_flags[i]){
				pthread_join(threads[i], NULL);

				finished_flags[i] = 0;

				stat_function(args[i]);

				clean_up(args[i]);

				threads[i] = 0;
				fasta_line = read_fasta_map(fasta_line);
				if(fasta_line != NULL){
					init_function(&(args[i]), &finished_flags[i], fasta_line, i, id++, compare, use_blooms);
					t_res = pthread_create(&(threads[i]), NULL, calculate_function, args[i]);

					if(t_res < 0) perror("ngul o riavl");
				}else{
					--no_started;
				}
			}
		}

		if(out_file > 0){
			fprintf(stderr, "Output file not yet supported in multi threaded mode\n");
		}

		sleep(0.2);
	}

	munmap(map, stats.st_size);

	return TRUE;
}


void init_standard(void **_args, BOOL *finished, char *word, int thread_no, size_t id, BOOL use_blooms, BOOL compare){
	sa_info *args = malloc(sizeof(sa_info));

	args->size = strlen(word);

	char *arg_word = malloc(sizeof(char) * args->size + 1);
	strncpy(arg_word, word, args->size + 1);

	args->word = arg_word;
	args->sa = NULL;
	args->index = thread_no;
	//todo controlla se ci vuole anche qua *finished
	args->finished = finished;
	*(args->finished) = 0;

	args->compare = compare;
	args->window_size = -1;
	args->use_bloom = use_blooms;
	args->id = id;

	args->node = NULL;

	*_args = args;
}


void cleaup_standard(void *_args){
	sa_info *args = (sa_info*) _args;
	if(args->node != NULL){
		free_list(args->node);
	}

	free(args->sa);
	free(args->word);
	free(args);
}


void *compute_standard(void *_args){
	sa_info *args = (sa_info *) _args;

	int **w_size_int;

	w_size_int =  window_size(args->word, "ACGTN", &(args->node));


	size_t *sa = compute_suffix_array(args->word, args->size, **w_size_int, FALSE);
	args->window_size = **w_size_int;
	getrusage(RUSAGE_THREAD, &args->suffix_rusage); // @suppress("Symbol is not resolved")

	free(*w_size_int);
	free(w_size_int);

	*(args->finished) = TRUE;
	return NULL;
}

int main(int nargs, char **vargs){


	BOOL is_compare = FALSE;
	BOOL use_blooms = FALSE;
	BOOL test_only = FALSE;
	struct option long_options[] = {
			{"filename", required_argument, NULL, 'f'},
			{"nothread", required_argument, NULL, 't'},
			{"output", required_argument, NULL, 'o'},
			{"compare",  no_argument, &is_compare, TRUE},
			{"use-blooms", no_argument, &use_blooms, TRUE},
			{0,0,0,0}
	};

	char opt;
	int opt_index = 1;
	char file_name[255] = {0};
	char output_file_name[255] = {0};

	BOOL icfl_report = FALSE;
	BOOL huge = FALSE;
	BOOL no_merge = FALSE;
	BOOL sais_only = FALSE;
	int no_thread = 1;
	while(TRUE){
		opt = getopt_long(nargs, vargs, "f:t:o:cbTINS", long_options, &opt_index);
		if(opt == -1) break;
		switch(opt){
		case 'f':
			strncpy(file_name, optarg, 255);
			break;
		case 't':
			no_thread = atoi(optarg);
			break;
		case 'o':
			strncpy(output_file_name, optarg, 255);
			break;
		case 'c':
			is_compare = TRUE;
			break;
		case 'b':
			use_blooms = TRUE;
			break;
		case 'T':
			test_only = TRUE;
			break;
		case 'I':
			icfl_report = TRUE;
			break;
		case 'N':
			no_merge = TRUE;
			break;
		case 'S':
			sais_only = TRUE;
			break;
		default:
			break;
		}
	}





	void* (*function)(void *) = thread_worker;
	if(test_only){
		function = thread_worker_only_sa;
	}
	if(sais_only){

		printf("id,#thread,word len,user time, sys time, min pg, maj pf, max res\n");
		if(generic_multi_thread(file_name, 0, no_thread, FALSE, FALSE, sais_init,
				sais_compute, sais_stat, sais_clean)){
			exit(EXIT_SUCCESS);
		}
		exit(EXIT_FAILURE);
	}



	//FILE *fp = fopen(file_name, "r");


	if(strlen(file_name) <= 0){
		fprintf(stderr, "Input filename is required\n");
		exit(EXIT_FAILURE);
	}





	//single_thread(fp, out_file, is_compare, use_blooms, only_memory);
	//return 0;



	if(icfl_report){
		printf("Word size, #fact, min fact len, max fact len, avg fact len, med fact len, #thread, seq id, user time, real time, Maj pf, Min pf, max rss\n");
		if(generic_multi_thread(file_name, 0, no_thread, FALSE, FALSE, icfl_report_init,
				icfl_report_calculate, icfl_report_print_stats, icfl_report_clean_up)){
			exit(EXIT_SUCCESS);
		}
		exit(EXIT_FAILURE);
	}

	if(no_merge){
		printf("Status,Thread,Size,#fact,Window size,Time user,Time real,Page fault,Hard Page Fault\n");

		if(generic_multi_thread(file_name, 0, no_thread, FALSE, FALSE, init_sa_no_merge,
				sais_icfl, no_merge_print_stat, no_clean_up)){
			exit(EXIT_SUCCESS);
		}
		exit(EXIT_FAILURE);
	}
	if(sais_only){
		if(generic_multi_thread(file_name, 0, no_thread, FALSE, FALSE, sais_init,
				sais_compute, sais_stat, sais_clean)){
			exit(EXIT_SUCCESS);
		}
		exit(EXIT_FAILURE);
	}



	FILE *fp = fopen(file_name, "r");
	if(fp <= 0){
		perror("Cannot open file");
		exit(EXIT_FAILURE);
	}

	int out_file = 0;
	if(strlen(output_file_name) > 0){
		out_file = open(output_file_name, O_WRONLY | O_CREAT | O_TRUNC, 0666);
		if(out_file <= 0){
			perror("Cannot open file");
			exit(EXIT_FAILURE);
		}
	}

	/*if(no_thread == 1){
		if(single_thread(fp, out_file, is_compare, use_blooms, function))
			exit(EXIT_SUCCESS);
		else
			exit(EXIT_FAILURE);
	/*}else{
		if(main_threads(fp, out_file, no_thread, is_compare, use_blooms, function)){
			exit(EXIT_SUCCESS);
		}else{
			exit(EXIT_FAILURE);
		}
	}*/
	printf("Status,Id,Thread,Size,#fact,Window size,Time user,Time real,Page fault,Hard Page Fault\n");
	if(generic_multi_thread(file_name, 0, no_thread, is_compare, FALSE, init_standard, compute_standard, print_stat, cleaup_standard)){
		exit(EXIT_SUCCESS);
	}else{
		exit(EXIT_FAILURE);

	}

	return 0;
}


/*
size_t suffix_worker(char *word, size_t size){

}


*/

/*

int main2(int nargs, char **vargs){
	int protection = PROT_READ | PROT_WRITE;
	int visibility = MAP_ANONYMOUS | MAP_SHARED;
	void *map = mmap(NULL, PROCS * sizeof())

	FILE *fp = openFileFromParams(nargs, vargs);
	char *fasta_line = NULL;
	int processes = 0;
	while( (fasta_line = read_fasta(fp)) != NULL){
		if(processes >= PROCS){
			wait(NULL);
			--processes;
		}
		++processes;
		int pid;
		if( (pid = fork()) < 0){
			perror("FORK ERROR");
			exit(0);
		}else if(pid == 0){
			suffix_worker();
			return 0;
		}


	}
	return 0;


}
*/
