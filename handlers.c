/*
 * handlers.c
 *
 *  Created on: Nov 23, 2019
 *      Author: vittorio
 */
#include "handlers.h"

//
//BOOL generic_multi_thread(char *file_name, int out_file, int no_threads, BOOL compare, BOOL use_blooms, void init_function(void*, BOOL*, char*, int, int, BOOL, BOOL),
//		void *calculate_function(void *), void *stat_function(void*), void clean_up(void*)){
//
//	char *fasta_line;
//	void *args[no_threads];
//	pthread_t threads[no_threads];
//	BOOL finished_flags[no_threads];
//	int t_res;
//	int current_threads;
//	size_t id = 0;
//
//	struct stat stats;
//	if(stat(file_name, &stats) < 0){
//		perror("stat");
//		return FALSE;
//	}
//
//	int fd = 0;
//	if( (fd = open(file_name, O_RDONLY)) == NULL){
//		perror("stat");
//		return FALSE;
//	}
//
//	void *map;
//	fasta_line = map = mmap(NULL, stats.st_size, PROT_READ | PROT_WRITE, MAP_PRIVATE, fd, 0);
//	if(fasta_line == MAP_FAILED){
//		perror("mmap");
//		return FALSE;
//	}
//
//
//
//	for(current_threads = 0; current_threads < no_threads; ++current_threads){
//
//		if((fasta_line = read_fasta_map(fasta_line)) != NULL){
//			init_function(&(args[current_threads]), &finished_flags[current_threads], fasta_line, current_threads, id++, compare, use_blooms);
//			t_res = pthread_create(&(threads[current_threads]), NULL, calculate_function, args[current_threads]);
//
//			if(t_res < 0) perror("ngul o riavl");
//
//		}else{
//			break;
//		}
//	}
//
//
//	int no_started = current_threads;
//	while(no_started){
//		for(int i = 0; i < current_threads; ++i){
//			if(finished_flags[i]){
//				pthread_join(threads[i], NULL);
//
//				finished_flags[i] = 0;
//
//				stat_function(args[i]);
//
//				clean_up(args[i]);
//
//				threads[i] = 0;
//				fasta_line = read_fasta_map(fasta_line);
//				if(fasta_line != NULL){
//					init_function(&(args[i]), &finished_flags[i], fasta_line, i, id++, compare, use_blooms);
//					t_res = pthread_create(&(threads[i]), NULL, calculate_function, args[i]);
//
//					if(t_res < 0) perror("ngul o riavl");
//				}else{
//					--no_started;
//				}
//			}
//		}
//
//		if(out_file > 0){
//			fprintf(stderr, "Output file not yet supported in multi threaded mode\n");
//		}
//
//		sleep(0.2);
//	}
//
//	munmap(map, stats.st_size);
//
//	return TRUE;
//}


//***************************
//		ICFL REPORT			*
//***************************
void print_icfl_report(icfl_report *report){
	printf("%ld,%ld,%ld,%ld,%f,%f,%d,%d,%ld.%06ld,%ld.%06ld,%ld,%ld,%ld\n",
			report->word_size,
			report->no_factors,
			report->min_fact_len,
			report->max_fact_len,
			report->avg_fact_len,
			report->med_fact_len,
			report->no_thread,
			report->seq_id,
			report->rusage.ru_utime.tv_sec,
			report->rusage.ru_utime.tv_usec,
			report->rusage.ru_stime.tv_sec,
			report->rusage.ru_stime.tv_usec,
			report->rusage.ru_majflt,
			report->rusage.ru_minflt,
			report->rusage.ru_maxrss);


	fflush(stdout);
}

void icfl_report_print_stats(void *void_report){
	icfl_report *report = (icfl_report *) void_report;

	print_icfl_report(report);
}

void icfl_report_clean_up(void * void_report){
	icfl_report *report = (icfl_report *) void_report;
	free(report->word);
	free_list(report->icfl);
	free(report);
}

void icfl_report_init(void** arg, BOOL* flag_finished, char* word, int no_thread, size_t seq_id, BOOL compare, BOOL use_bloom){
	icfl_report *report = malloc(sizeof(icfl_report));

	report->word_size = strlen(word);

	report->word = malloc(sizeof(char) * (report->word_size + 1));
	strcpy(report->word, word);

	report->seq_id = seq_id;
	report->no_thread = no_thread;

	report->finished = flag_finished;
	*(report->finished) = FALSE;

	*arg = report;
}


void *icfl_report_calculate(void *arg){
	icfl_report *report = (icfl_report *) arg;

	size_t w_size =  window_size(report->word, "ACGTN", &report->icfl);

	if(getrusage(RUSAGE_THREAD, &report->rusage) < 0){  // @suppress("Symbol is not resolved")
		perror("rusage");
	}


	node_t *fact = report->icfl;

	if(fact == NULL){
		perror("Negg");
		return NULL;
	}
	size_t current_length = strlen(fact->factor);
	size_t sum = current_length;
	report->no_factors = 1;
	report->max_fact_len = current_length;
	report->min_fact_len = current_length;
	report->window_size = w_size;
	report->med_fact_len = 0;

	while((fact = fact->next) != NULL){
		current_length = strlen(fact->factor);
		sum += current_length;
		++report->no_factors;

		if(report->max_fact_len < current_length)
			report->max_fact_len = current_length;
		else if(report->min_fact_len > current_length)
			report->min_fact_len = current_length;
	}

	report->avg_fact_len = sum / (double) report->no_factors;

	*(report->finished) = TRUE;

	return report;

}

//*****************************************************************************+
//					NON GENERIC
//****************************************************************************



#define FASTA_BUFF_DIM 2048
#define MAX_WORD_LEN 1024
#define ALPHA "NACGTnacgtS"

char fasta_buff[FASTA_BUFF_DIM + 1];
int fasta_buff_position = 0;


char _word_fasta[BUFSIZ * 1024 * 4];


char *read_fasta(FILE *fp){
	_word_fasta[0] = '\0';
	//">%*[^\n]\n%m[" ALPHA "]\n" ;
	char * pattern = ">%*[^\n]\n%[" ALPHA "\n]";
	if(fscanf(fp, pattern , _word_fasta) == EOF){
		return NULL;
	}
	size_t shift = 0;
	for(int i = 0; i < strlen(_word_fasta); ++i){
		if(_word_fasta[i] == '\n') {
			++shift;
		}else{
			_word_fasta[i - shift] = _word_fasta[i];
		}
	}
	_word_fasta[strlen(_word_fasta) - shift + 1] = '\0';
	return _word_fasta;
}



BOOL main_threads(FILE *fp, int out_file, int no_threads, BOOL compare, BOOL use_blooms, void *calculate_funxtion(void *)){

	char *fasta_line;
	sa_info args[no_threads];

	int t_res;
	int current_threads;
	size_t id = 0;
	printf("Status,Thread,Size,No fact,Window size,Time user,Time real,Page fault,Hard Page Fault\n");
	for(current_threads = 0; current_threads < no_threads; ++current_threads){
		args[current_threads].word = NULL;
		if((fasta_line = read_fasta(fp)) != NULL){
			init_suffix_arg(&(args[current_threads]), fasta_line, current_threads, compare, use_blooms, id++);
			t_res = pthread_create(&(args[current_threads].thread), NULL, calculate_funxtion, &args[current_threads]);
		}else{
			break;
		}
	}


	int no_started = current_threads;
	while(no_started){
		for(int i = 0; i < current_threads; ++i){
			if(args[i].finished){
				pthread_join(args[i].thread, NULL);

				args[i].finished = 0;

				//print_stat(args[i]);

				free(args[i].word);
				free(args[i].sa);
				args[i].thread = 0;
				fasta_line = read_fasta(fp);
				if(fasta_line != NULL){
					init_suffix_arg(&(args[i]), fasta_line, i, compare, use_blooms, id++);
					t_res = pthread_create(&(args[i].thread), NULL, thread_worker, &args[i]);
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
	return TRUE;
}

BOOL single_thread(FILE *fp, int out_file, BOOL compare, BOOL use_blooms, void*calculate_function(void *)){

	printf("Status,Thread,Size,Window size,Time user,Time real,Page fault,Hard Page Fault\n");
	char *fasta_line = NULL;
	size_t id = 0;
	while((fasta_line = read_fasta(fp)) != NULL){
		sa_info arg;
		size_t size = strlen(fasta_line);
		if(size <= 0) continue;

		init_suffix_arg(&arg, fasta_line, 1, compare, use_blooms, id++);
		calculate_function(&arg);

		//print_stat(arg);

		free(arg.word);
		free(arg.sa);

		if(out_file > 0){
			write(out_file, arg.sa, sizeof(size_t));
			write(out_file, arg.sa, arg.size);
		}
	}
	return TRUE;

}




void *only_memory(void * args){
	sa_info *st_arg = (sa_info*) args;

	size_t arr[st_arg->size];
	size_t arr2[st_arg->size];
	for(int i = 0; i < st_arg->size; ++i){
		arr[i] = i;
	}
	for(int i = 0; i < st_arg->size; ++i){
		for(int j = 0; j < st_arg->size; ++j){
			arr2[j] = arr[j];
		}
	}
}

size_t radix_suffix_array(char *word, size_t word_len, int *w_size, BOOL use_bloom){
	//TODO change internal logic to use size_t
	int **w_size_int;
	node_t *node;
	w_size_int =  window_size(word, "ACGTN", &node);

	node_t *iter_node = node;
	int i = 0;
	/*while(iter_node != NULL){
		sa = compute_suffix_array(iter_node->factor, strlen(iter_node->factor), strlen(iter_node->factor), use_bloom);
		iter_node = iter_node->next;
		i++;
	}*/


	size_t *sa = compute_suffix_array(word, word_len, **w_size_int, use_bloom);

	*w_size = **w_size_int;

	free(*w_size_int);
	free(w_size_int);
	//print_SA_ld(sa, word_len);
	return i;
}



int *lib_sa(sa_info *arg){
	int *sa_lib = malloc(sizeof(int) * arg->size);
	divsufsort((unsigned char*)arg->word, sa_lib, arg->size);
	return sa_lib;
}

int compare_to_libs(sa_info *arg){
	int *sa_lib = lib_sa(arg);
	int i = 0;
	for(; i < arg->size; ++i){
		if(sa_lib[i] != arg->sa[i]) break;
	}
	free(sa_lib);

	if(i == arg->size) arg->status = 1;
	else arg->status = 0;

	return arg->status;
}


void copy_info(sa_info *arg){
	struct rusage thread_rusage;
	int res = getrusage(RUSAGE_THREAD, &thread_rusage); // @suppress("Symbol is not resolved")
		if(res >= 0){
			arg->thread_stat.hard_page_faults = thread_rusage.ru_majflt;
			arg->thread_stat.page_faults = thread_rusage.ru_minflt;
			arg->thread_stat.real_time = thread_rusage.ru_stime;
			arg->thread_stat.user_time = thread_rusage.ru_utime;
			arg->status = 1;
		}
}

void *thread_worker(void *args){
	sa_info *struct_args = ((sa_info*)args);
	int window_size = NULL;

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


/*	if(struct_args->compare){
		compare_to_libs(struct_args);
	}else{
		struct_args->status = 2;
	}*/

	*(struct_args)->finished = 1;
	return args;
}



void *thread_worker_only_sa(void *args){
	sa_info *struct_args = ((sa_info*)args);

	if(strcmp(struct_args->word,"") == 0 ||
			struct_args->size <= 0 ||
			struct_args->window_size <= 0){

		fprintf(stderr, "Error on word: %s\tsize: %ld\tw_size: %ld\n", struct_args->word, struct_args->size, struct_args->window_size);
		fflush(stderr);

		struct_args->status = 0;
		struct_args->finished = 1;
		return args;
	}

	int *sa = lib_sa(struct_args);
	free(sa);

	int res = getrusage(RUSAGE_THREAD, &struct_args->suffix_rusage); // @suppress("Symbol is not resolved")

	struct_args->finished = 1;
	return args;
}



void print_stat(void *_args){
	sa_info *args = (sa_info *) _args;
	printf("%d,%ld,%d,%ld,%ld,%ld,%ld.%06ld,%ld.%06ld,%ld,%ld\n",
			args->status,
			args->id,
			args->index,
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



/*************************************************************
 * 	SAIS ICFL
 *
 ****************************************************************/



