#include "main.h"


// function declarations
void dbg_print_indiv(struct Individual);
int search(char **list, char *item, int lower, int upper);
CCDS* make_master_list(int&);
CCDS* remove_duplicates(CCDS*, int&);
//void analyze(CCDS*, int, Individual*, int);
//void analyze2(CCDS*, int, Individual*, int);
void em(CCDS*, int, Individual*, Parameters*);
void single_em_iteration_no_hwe(Parameters *para, int *R, int *N, int *indices, int valid_individuals, double input[], double output[]);
void cov_matrix_no_hwe(Parameters *para, int *R, int *N, int *indices, int valid_individuals, double **para_estimates, int iterations, double sum_h, double sum_i, FILE *fptr);

double absf(const double &val);
int largest_of_three(double a, double b, double c);
//void print_result(CCDS*, int, Individual*, int);
void print_result(CCDS *master_list, int master_size, struct Individual *data, struct Parameters *para);


// Global variables
struct Main_list master;
struct Individual *data;
struct Parameters *para;
struct Sync_sorter *sort_syncher;
struct Sync_searcher *search_syncher;
struct Sync_ana *ana_syncher;

// Thread related
struct Help_em{
	int num_threads;						// number of threads
	int id;									// thread Id, 0..(numthreads-1)
	int start;								// 1st SNP to process
	int stop;								// last SNP to process
	Analyzer ana_instance;
};

struct Help_reader{
	Read_data reader_instance;
	int min_depth;
};

struct Help_sorter{
	Sorter_t sorter_instance;
	int thread_id;							// identifies thread. More a debugging feature than anything else

};

struct Help_searcher{
	Searcher_t searcher_instance;
	int thread_id;
};

void* reader_thread_startup(void* arg);
void* em_thread_startup(void*);
void* sorter_thread_startup(void*);
void* searcher_thread_startup(void*);

// main

int main(int argc, char *argv[]){
	char inctrl[100];
//#ifdef PROFILE
//	PROFILE_START("main");
//#endif


	if(argc != 2){
		fprintf(stderr,"%s\n", version);
		fprintf(stderr,"Useage: name_of_SEQEM_executable <SEQEM control file>\n");
		exit(0);
	}
	int i;
	strcpy(inctrl, argv[1]);
	Readparam preader;
	para = preader.master_reader(inctrl);
	// start sync data acquisition
printf("reading & sorting\n");
	data = NULL;
	data = new struct Individual[para->num_header_files];
	sort_syncher = new struct Sync_sorter;
	sort_syncher->items = para->num_header_files;
	sort_syncher->done = 0;
	sort_syncher->ready = 0;
	sort_syncher->search_start = 0;
	sort_syncher->list = new int[sort_syncher->items];
//#ifdef PROFILE
//	PROFILE_START("reading_sorting");
//#endif
	int rc;
	int tc = pthread_mutex_init(&sort_syncher->sync_sorter_mutex , NULL);
	if(tc){
		fprintf(stderr, "main::Failed to initialize sort_syncher->sync_sorter_mutex. Bailing out.\n");
		exit(0);
	}
	tc = pthread_cond_init(&sort_syncher->cond_waiting_consumer, NULL);
	if(tc){
		fprintf(stderr, "main::Failed to initialize sort_syncher->cond_waiting_consumer. Bailing out.\n");
		exit(0);
	}
	memset((void*)sort_syncher->list, 0, sort_syncher->items*sizeof(int));
	pthread_t reader_t;
	Help_reader *arg_r = new Help_reader;
	arg_r->min_depth = para->min_depth;
	rc = pthread_create(&reader_t, NULL, reader_thread_startup, (void*)arg_r);
	int num_sorter_threads = para->threads - 1;
	if(num_sorter_threads <= 0){
		num_sorter_threads = 1;
	}
	Help_sorter *arg_s = new Help_sorter[num_sorter_threads];
	pthread_t *thread_sorter;
	thread_sorter = new pthread_t[num_sorter_threads];
	for(i = 0; i < num_sorter_threads; i++){
		arg_s[i].thread_id = i;
		rc = pthread_create(&thread_sorter[i], NULL, sorter_thread_startup, (void*)&arg_s[i]);
	}
// now wait for threads to return
// first the reader thread as it must finish first
	rc = pthread_join(reader_t, NULL);
	if(rc){
		fprintf(stderr,"main:: Failed to join reader-thread. Error code %i. Bailing out.\n", rc);
		exit(0);
	}
	printf("done reading\n");
// now wait for sorter-threads
	for(i = 0; i < num_sorter_threads; i++){
		rc = pthread_join(thread_sorter[i], NULL);
		if(rc){
			fprintf(stderr,"main:: Failed to join sorter-thread number %i. Error code %i. Bailing out.\n", i, rc);
			exit(0);
		}
	}
	delete [] sort_syncher->list;
	sort_syncher->list = NULL;
	delete sort_syncher;
	sort_syncher = NULL;
	delete [] thread_sorter;
	delete [] arg_s;
printf("done sorting\n");
	// end sync data acquisition
//#ifdef PROFILE
//    PROFILE_STOP("reading_sorting");
//#endif
//#ifdef PROFILE
//	PROFILE_START("master_list");
//#endif
	int master_size = 0;
	// Now create a master_list that contains all SNPs from al individuals
	// Strategy:
	// We use data[0].ccds_data as reference and search data[n].ccds_data against it, where n = 1..(individuals-1)
	// Every time we find a SNP/item that's not in the reference list we add it to a linked list.
	// Eventually we will transfer the linked list into an array, sort it, remove duplicates and combine it
	// with the reference list. We will call the resulting list the master_list as it contains all SNPs.
	search_syncher = new struct Sync_searcher;
	search_syncher->next = 1;
	search_syncher->items = para->num_header_files;
	search_syncher->list_size = 0;
	search_syncher->head = NULL;
	search_syncher->current = NULL;
	search_syncher->previous = NULL;
	tc = pthread_mutex_init(&search_syncher->sync_searcher_mutex, NULL);
	if(tc){
		fprintf(stderr, "main::Failed to initialize search_syncher->sync_searcher_mutex. Bailing out.\n");
		exit(0);
	}
printf("start masterlist\n");
	pthread_t *thread_searcher;
	thread_searcher = new pthread_t[para->threads];
	Help_searcher *arg_ser = new Help_searcher[para->threads];
	// start all searcher threads
	for(i = 0; i < para->threads; i++){
		arg_ser[i].thread_id = i;
		rc = pthread_create(&thread_searcher[i], NULL, searcher_thread_startup, (void*)&arg_ser[i]);
	}

	// wait for threads to return
	for(i = 0; i < para->threads; i++){
		rc = pthread_join(thread_searcher[i], NULL);
		if(rc){
			fprintf(stderr,"main:: Failed to join sorter-thread number %i. Error code %i. Bailing out.\n", i, rc);
			exit(0);
		}
	}
// now merge the linked list with data[0]
	CCDS *master_list = make_master_list(master_size);
// clean up
	delete search_syncher;
	delete [] thread_searcher;
	delete [] arg_ser;
	// end searcher threads
	
	master.size = master_size;
	master.list = master_list;
	if(para->debug & 2){
		printf("### Master_list, reads = %i\n", master_size);
		for(i = 0; i < master_size; i++){
			printf("%i  %i start= %i stop= %i chr= %i\n",i, master_list[i].ccds_name, master_list[i].start, master_list[i].stop, master_list[i].chr);
		}
	}
printf("end masterlist\n");
//#ifdef PROFILE
//    PROFILE_STOP("master_list");
//#endif
	// create EM-threads
//#ifdef PROFILE
//	PROFILE_START("EM");
//#endif
printf("start EM\n");
	ana_syncher = new struct Sync_ana;
	ana_syncher->current = 0;
	ana_syncher->items = master.size;
	ana_syncher->threads = para->threads;
	tc = pthread_mutex_init(&ana_syncher->sync_ana_mutex , NULL);
	if(tc){
		fprintf(stderr, "main::Failed to initialize ana_syncher->sync_ana_mutex. Bailing out.\n");
		exit(0);
	}

	pthread_t *thread_em;
	thread_em = new pthread_t[para->threads];
	struct Help_em *arg2;
	arg2 = new struct Help_em[para->threads];
	int start, stop;										// first & last SNP to process for each thread
	stop = -1;
	for(i = 0; i < para->threads; i++){
		arg2[i].id = i;
		arg2[i].num_threads = para->threads;
		start = stop + 1;
		stop = start + (int)(1.7*((master_size-1-start) / (2*(para->threads-i))));
		arg2[i].start = start;
		arg2[i].stop = stop;
		ana_syncher->current = stop + 1;
		rc = pthread_create(&thread_em[i], NULL, em_thread_startup, (void*)&arg2[i]); 
		if(rc){
			fprintf(stderr,"main:: Failed to create thread_em[%i]. Bailing out\n",i);
		}	
	}
	// now wait for EM-threads to exit
	
	for(i = 0; i < para->threads; i++){
		rc = pthread_join(thread_em[i], NULL);
		if(rc){
			fprintf(stderr," main:: Failed to join thread_em[%i]. Bailing out\n", i);
			exit(0);
		}
	}
printf("end EM\n");
//#ifdef PROFILE
//    PROFILE_STOP("EM");
//#endif
// done, print result
	print_result(master_list, master_size, data, para);
//#ifdef PROFILE
//    PROFILE_STOP("main");
//#endif
//#ifdef PROFILE
//	ProfilePrint();
//#endif
	return(0);
}

//*************************************

void dbg_print_indiv(struct Individual indi)
{
	int i;
	int pos;
	Binsearch <CCDSI> searcher;
	printf("### Invividual reads= %i\n", indi.total_reads);
	for(i = 0; i < indi.total_reads; i++){
		printf("%i  %i\n",i,indi.ccds_data[i].location);
		if(0 == i%20){
			pos = searcher.search(indi.ccds_data, indi.ccds_data[i],0,indi.total_reads-1);
			printf("### position of %i is %i\n", indi.ccds_data[i].location, pos);
		}
	}
}


//*************************************

int search(char **list, char *item, int lower, int upper)
{
	int middle;
	int position = -1;
	int temp;
	while(upper >= lower){
		middle = (lower + upper) / 2;
		temp = strcmp(list[middle], item);
		if(0 == temp){
			position = middle;
			break;
		}
		else if(temp < 0){
			lower = middle + 1;
		}
		else{
			upper = middle - 1;
		}
	}
	return(position);
}

//*************************************

CCDS* make_master_list(int &size)
// 1. transform the linked list from the previous step into a sortable array
// 2. sort the array
// 3. remove duplicates
// 4. merge it with the list from data[0] into a master-list that contains ALL SNPs
{
	CCDS *result = NULL;
	Binsearch <CCDS> searcher;
	Sort <CCDS> sorter;
	int i,j;


	struct Node *cur = NULL;
	struct Node *prev = NULL;
	int position;
	size = search_syncher->list_size;

	if(0 < size){
		CCDS *new_ccds_list = new CCDS[size];								// now make a sortable array out of linked list
		cur = search_syncher->head;
		for(i = 0; i < size; i++){
			new_ccds_list[i].ccds_name = cur->item.location;
			prev = cur;
			cur = cur->next;
			delete prev;
			prev = NULL;
		}
		sorter.init(size, new_ccds_list);									// sort array made from linked list
		sorter.qsort();
		new_ccds_list = remove_duplicates(new_ccds_list, size);				// remove duplicates from sorted array
		result = new CCDS[data[0].total_reads+size];						// allocate space fr master list, size_of_reference + size_of_no_duplicates
		for(i = 0; i < data[0].total_reads; i++){							// add item from reference
			result[i].ccds_name = data[0].ccds_data[i].location;
		}
		size += data[0].total_reads;
		for(i = data[0].total_reads, j = 0; i < size; i++, j++){			// add item from no duplicates
			result[i].ccds_name = new_ccds_list[j].ccds_name;
		}
		sorter.init(size, result);											// sort master list
		sorter.qsort();
	}
	else{																	// the reference list contains ALL elements
		result = new CCDS[data[0].total_reads];
		for(i = 0; i < data[0].total_reads; i++){
			result[i].ccds_name = data[0].ccds_data[i].location;
		}
		size = data[0].total_reads;
	}
	return(result);
}

//*************************************

CCDS* remove_duplicates(CCDS* list, int &size)
{
	int i;
	int newsize = 1;
	struct Node *head = new struct Node;										// head of linked list
	head->item.location = list[0].ccds_name;
	head->next = NULL;
	struct Node *current = head;
	struct Node *previous = NULL;
	int prev;
	for(i = 1, prev = 0; i < size; i++, prev++){
		if(list[i] != list[prev]){
			current->next = new Node;
			current = current->next;
			current->item.location = list[i].ccds_name;
			current->next = NULL;
			newsize++;
		}
	}
	delete [] list;
	list = new CCDS[newsize];
	current = head;
	for(i = 0; i < newsize; i++){
		list[i].ccds_name = current->item.location;
		previous = current;
		current = current->next;
		delete previous;
		previous = NULL;
	}
	size = newsize;
	return(list);
}

//*************************************

double absf(const double &val)
{
	double result = val;
	if(result < 0) result *= -1;
	return(result);
}

//*************************************

int largest_of_three(double a, double b, double c)
{
	int result = 0;

	if(a > b){
		result = 0;
		if(c > a) result = 2;
	}
	else{
		result = 1;
		if(c > b) result = 2;
	}
	return(result);
}

//*************************************

void* em_thread_startup(void *arg)
{
	struct Help_em *arg1 = (struct Help_em*)arg;
	struct EM_info *em_arg;
	em_arg = new struct EM_info;
	em_arg->id = arg1->id;
	em_arg->num_threads = arg1->num_threads;
	em_arg->start = arg1->start;
	em_arg->stop = arg1->stop;
	arg1->ana_instance.run(em_arg);
	return(NULL);
}

//*************************************

void print_result(CCDS *masterlist, int mastersize, struct Individual *data, struct Parameters *para)
{
	int i,j;
	int temp;
	char rr[3] = {'r','r','\0'};
	char rv[3] = {'r','v','\0'};
	char vv[3] = {'v','v','\0'};
	Sort <CCDS> sorter;
	Binsearch <CCDSI> searcher;
	CCDSI item;
	sorter.init(mastersize, masterlist);
	sorter.qsort();
	FILE *fptr_seq = fopen(para->outfile, "w");
	fprintf(fptr_seq,"### Genotypes lexographically sorted by name, start and stop position. There are a total of %i records ###\n", mastersize);
	fprintf(fptr_seq,"genotypes: r = reference, v = variant, na = genotype not determined\n");
	for(i = 0; i < mastersize; i++){
		fprintf(fptr_seq,"%10i. %i  iter= %i  error= %7.5f  p= %7.5f hwe-chisq= %8.6f\n",i,masterlist[i].ccds_name,masterlist[i].iterations,masterlist[i].error,masterlist[i].p, masterlist[i].stat);
		if(masterlist[i].prr >= 0){
			fprintf(fptr_seq,"P(rr)= %7.5f p(rv)= %7.5f p(vv)= %7.5f\n", masterlist[i].prr,masterlist[i].prv,masterlist[i].pvv);
		}
		else if(para->hwe){
			fprintf(fptr_seq,"P(rr) according to assumed HWE\n");
		}
		else{
			fprintf(fptr_seq,"P(rr)= -1.0 p(rv)= -1.0 p(vv)= -1.0 ...general error condition\n");
		}
		for(j = 0; j < para->num_header_files; j++){
			item.location = masterlist[i].ccds_name;
			temp = searcher.search(data[j].ccds_data, item, 0, data[j].total_reads-1);
			if(temp >= 0){
				if(0 == masterlist[i].genotype[j]){
					fprintf(fptr_seq,"           genotype[%i]= %s props(rr,rv,vv)= %6.3f %6.3f %6.3f    reads= %6i , variants= %6i\n", j, rr,masterlist[i].geno_prop[0][j], masterlist[i].geno_prop[1][j], masterlist[i].geno_prop[2][j], data[j].ccds_data[temp].depth, data[j].ccds_data[temp].variants);
				}
				else if(1 == masterlist[i].genotype[j]){
					fprintf(fptr_seq,"           genotype[%i]= %s props(rr,rv,vv)= %6.3f %6.3f %6.3f    reads= %6i , variants= %6i\n", j, rv,masterlist[i].geno_prop[0][j], masterlist[i].geno_prop[1][j], masterlist[i].geno_prop[2][j], data[j].ccds_data[temp].depth, data[j].ccds_data[temp].variants);
				}
				else if(2 == masterlist[i].genotype[j]){
					fprintf(fptr_seq,"           genotype[%i]= %s props(rr,rv,vv)= %6.3f %6.3f %6.3f    reads= %6i , variants= %6i\n", j, vv,masterlist[i].geno_prop[0][j], masterlist[i].geno_prop[1][j], masterlist[i].geno_prop[2][j], data[j].ccds_data[temp].depth, data[j].ccds_data[temp].variants);
				}
				else{
					fprintf(fptr_seq,"           genotype[%i]= na \n", j);
				}
			}
			else{
				fprintf(fptr_seq,"           genotype[%i]= na \n", j);
			}
		}
		fprintf(fptr_seq,"\n");
	}
	fclose(fptr_seq);
}

//*************************************

void* reader_thread_startup(void* arg)
{
	Help_reader *arg1 = (Help_reader *)arg;
	arg1->reader_instance.main_reader(arg1->min_depth);
	return(NULL);
}

//*************************************

void* sorter_thread_startup(void* arg)
{
	struct Help_sorter *arg1 = (struct Help_sorter *)arg;
	arg1->sorter_instance.run(arg1->thread_id);
	return(NULL);
}

//*************************************

void* searcher_thread_startup(void* arg)
{
	struct Help_searcher *arg1 = (struct Help_searcher *)arg;
	arg1->searcher_instance.run(arg1->thread_id);
	return(NULL);
}

//*************************************
