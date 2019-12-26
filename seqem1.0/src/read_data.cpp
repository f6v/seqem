#include "read_data.h"

//*********************************

Read_data::Read_data()
{
}

//*********************************

Read_data::~Read_data()
{

}

//*********************************

void Read_data::main_reader(int mindepth)
{
	const int startstopsize = 10;
	ifstream infile;
	char space[2] = {' ', '\0'};
	char chr[5] = {'>', 'c', 'h', 'r', '\0'};
	char *line = new char[max_rdln];
	char *token;
	char *variant;
	char *reference;

	char *ascii = NULL;
	unsigned char u8hash;
	unsigned short u16hash;
	unsigned int u32hash;


	char ccdsname[40];
	char temp[10];
	char number[startstopsize];
	int eat_header = 0;

//	Sort <CCDS> sorter;
	int i, j, k, m, size;
	int count;
	int lines;
	unsigned int location;
	int depth;
	int variants;

	for(i = 0; i < sort_syncher->items ; i++){
		infile.open(para->header_file_names[i], ios::in);
		if(infile.is_open()){
			infile.clear();													// reset file to beginning for reading
			infile.seekg(0, ios::beg);										// normally this is not necessary after successfully opening a file, but for some bizarre reason this fails to work here

// we are trying to read something like:
// >chr1   1244704 1244704 C       G       3       100%    G       G       -3      CPSF3L  rs10907179
				lines = 0;
				while(infile.getline(line, max_rdln-1)){						// first count lines == CCDS labels
					token = strstr(line, chr);
					if(token != NULL){
						token = strtok(line, " ");
						token = strtok(NULL, " ");
						token = strtok(NULL, " ");
						token = strtok(NULL, " ");
						token = strtok(NULL, " ");
						token = strtok(NULL, " ");
						depth = atoi(token);
//						if((depth >= mindepth) && (depth < 256)){
						if(depth >= mindepth){
							lines++;
						}
					}
					else{
						eat_header++;
					}
				}
				infile.clear();													// reset file to beginning for reading
				infile.seekg(0, ios::beg);
				data[i].ccds_data = new CCDSI[lines];
				data[i].total_reads = lines;
				j = 0;
				while(eat_header && infile.getline(line, max_rdln-1)){			// eat throught whatever many header lines.
					eat_header--;
				}
				while(infile.getline(line, max_rdln-1)){						// now we actually start reading the data
					token = strtok(line, " ");									// this should look something like ">chr1" or ">chr15"
					token = strtok(NULL, " ");										// rs-name, irrelevant
					token = strtok(NULL, " ");										// location, unique identifier
//					idstr += token;

//					data[i].ccds_data[j].location = (unsigned int)atoi(token);		// might need to wite own routine cos there more base pairs than 2^31-1
					location = (unsigned int)atoi(token);
					reference = strtok(NULL, " ");										// reference
					size = (int)strlen(reference);

/*					idstr += token;
					idstr += "_";
*/
//					count = (int)strlen(token);
//					data[i].ccds_data[j].reference = new char[count+1];
//					strcpy(data[i].ccds_data[j].reference, token);

					variant = strtok(NULL, " ");										// variant
					size += (int)strlen(variant);
					size += 2;
					if(ascii != NULL){
						delete [] ascii;
					}
					ascii = new char[size];
					strcpy(ascii, reference);
					strcat(ascii, "_");
					strcat(ascii, variant);
					u8hash = hash_atou8(ascii, size-1);

	// test start
					/*
strcpy(ascii, "a_c");
u8hash = hash_atou8(ascii, 3);
strcpy(ascii, "a_g");
u8hash = hash_atou8(ascii, 3);
strcpy(ascii, "a_t");
u8hash = hash_atou8(ascii, 3);
strcpy(ascii, "a_-");
u8hash = hash_atou8(ascii, 3);
strcpy(ascii, "a_i");
u8hash = hash_atou8(ascii, 3);
strcpy(ascii, "c_a");
u8hash = hash_atou8(ascii, 3);
strcpy(ascii, "c_g");
u8hash = hash_atou8(ascii, 3);
strcpy(ascii, "c_t");
u8hash = hash_atou8(ascii, 3);
strcpy(ascii, "c_-");
u8hash = hash_atou8(ascii, 3);
strcpy(ascii, "c_i");
u8hash = hash_atou8(ascii, 3);
strcpy(ascii, "g_c");
u8hash = hash_atou8(ascii, 3);
strcpy(ascii, "g_a");
u8hash = hash_atou8(ascii, 3);
strcpy(ascii, "g_t");
u8hash = hash_atou8(ascii, 3);
strcpy(ascii, "g_-");
u8hash = hash_atou8(ascii, 3);
strcpy(ascii, "g_i");
u8hash = hash_atou8(ascii, 3);
strcpy(ascii, "t_c");
u8hash = hash_atou8(ascii, 3);
strcpy(ascii, "t_g");
u8hash = hash_atou8(ascii, 3);
strcpy(ascii, "t_a");
u8hash = hash_atou8(ascii, 3);
strcpy(ascii, "t_-");
u8hash = hash_atou8(ascii, 3);
strcpy(ascii, "t_i");
u8hash = hash_atou8(ascii, 3);
*/
	// test end
//					u16hash = hash_atou16(ascii, size-1);
//					u32hash = hash_atou32(ascii, size-1);
//					size = sizeof(CCDS);
/*
					idstr += token;
					md5code = md5instance.getHashFromString(idstr);
					strcpy(md5ptr, md5code.c_str());
					memcpy((void*)md5ptr_hi, (void*)md5ptr, 16*sizeof(char));
					memcpy((void*)md5ptr_lo, (void*)(md5ptr+16), 16*sizeof(char));
					md5ptr_hi[16] = '\0';
					md5ptr_lo[16] = '\0';
*/
					
					token = strtok(NULL, " ");										// read depth
//					data[i].ccds_data[j].depth = (unsigned char)atoi(token);
					depth = atoi(token);
//					if((depth >= mindepth) && (depth < 256)){
					if(depth >= mindepth){
						data[i].ccds_data[j].depth = (unsigned int)depth;
						data[i].ccds_data[j].location = location;
						data[i].ccds_data[j].variant_hash = u8hash;
						token = strtok(NULL, " ");										// variant%
/*						k = 0;
						while((token[k] >= '0') && (token[k] <= '9') && (k < (startstopsize-1))){
							number[k] = token[k];
							k++;
						}
						number[k] = '\0';
						data[i].ccds_data[j].variants = (unsigned char)atoi(number);
*/
						data[i].ccds_data[j].variants = (unsigned int)atoi(token);
						j++;
					}
				} // while(infile.getline(line, max_rdln-1))
//			} // else if(para->input == 2)
/*

	char *ccds_name;								// CCDS name/ID whatever you want to call it, always starts with ">"
	int start;										// where read starts
	int stop;										// where read stops...not important right now, but we might as well store the info
	int chr;										// chromosome. 1..22, X == 23, Y == 24
	char *reference;								// reference or expected read
	char *variant;									// what we read instead of the reference
	int num_reads;									// number of reads
	double proportion;								// of the variant in relation to "num_reads"
	unsigned char *presentin;						// 0|1, only for master record, present in individual or not
	int individuals;								// only for master record, size of presentin
	int iterations;									// how many iterations until it converges, only in masterlist 
	int *genotype;									// -1 = undetermined, 0 = "ref-ref", 1 = "ref-var", 2 = "var-var", only for master record
	double **geno_prop;								// instead of genotype is lists the probabilities for the 3 genotypes
	double error;									// only keep track in master_list;
	double p;										// prob. of reference, only keep trach in master_list
	double prr;
	double prv;
	double pvv;
*/
		}
		else{
			fprintf(stderr,"Read_data::main_reader. Failed to open %s. Bailing out.\n", para->header_file_names[i]);
			exit(0);			
		}
		infile.close();
		pthread_mutex_lock(&sort_syncher->sync_sorter_mutex);
		sort_syncher->list[i] = 1;
		sort_syncher->ready++;												// increment the item count that is ready to be processed
		pthread_mutex_unlock(&sort_syncher->sync_sorter_mutex);
		pthread_cond_broadcast(&sort_syncher->cond_waiting_consumer);		// let sorter-thread know we have some data to be processed
	} // for(i = 0; i < para->num_header_files; i++)
	delete [] line;
	if(ascii != NULL){
		delete [] ascii;
	}
	/*
	delete [] md5ptr;
	delete [] md5ptr_hi;
	delete [] md5ptr_lo;
	*/
}

//*********************************

unsigned char Read_data::hash_atou8(char *input, int size)
{
	unsigned char result = 0;
	int i;

	for(i = 0; i < size; i++){
		if(0 == (i % 2)){
			result ^= input[i] << 1;
		}
		else{
			 result ^= input[i];
		}		
	}
	return(result);
}

//*********************************

unsigned short Read_data::hash_atou16(char *input, int size)
{
	unsigned short result = 0;
	unsigned short temp;
	int i;
	int mod;

	for(i = 0; i < size; i++){
		mod = i % 4;
		switch(mod){
			case 0:
				result ^= input[i];
			break;
			case 1:
				temp = input[i];
				temp <<= 8;
				result ^= temp;
			break;
			case 2:
				result ^= input[i] << 1;
			break;
			case 3:
				temp = input[i];
				temp <<= 9;
				result ^= temp;
			break;
		}
	}
	return(result);
}

//*********************************

unsigned int Read_data::hash_atou32(char *input, int size)
{
	unsigned int result = 0;
	unsigned int temp;
	int i;
	int mod;

	for(i = 0; i < size; i++){
		mod = i % 8;
		switch(mod){
			case 0:
				result ^= input[i];
			break;
			case 1:
				temp = input[i];
				temp <<= 8;
				result ^= temp;
			break;
			case 2:
				temp = input[i];
				temp <<= 16;
				result ^= temp;
			break;
			case 3:
				temp = input[i];
				temp <<= 24;
				result ^= temp;
			break;
			case 4:
				result ^= input[i] << 1;
			break;
			case 5:
				temp = input[i];
				temp <<= 9;
				result ^= temp;
			break;
			case 6:
				temp = input[i];
				temp <<= 17;
				result ^= temp;
			break;
			case 7:
				temp = input[i];
				temp <<= 25;
				result ^= temp;
			break;

		}
	}
	return(result);
}

//*********************************

//*********************************

//*********************************

//*********************************

//*********************************

//*********************************

