#include "read_param.h"

//**************************************

Readparam::Readparam()
{

}

//**************************************

Readparam::~Readparam()
{

}

//**************************************

void Readparam::set_file_name(char name[])
{
	strcpy(filename, name);
}

//**************************************

int Readparam::open_infile(ifstream& infile)
{
	infile.open(filename, ios::in);
	if (infile.is_open()){
		return(1);
	}
	else{
		return(0);
	}
}

//**************************************

struct Parameters* Readparam::master_reader(char filename[])
{
	ifstream infile;
	set_file_name(filename);
	struct Parameters* seq_params = new struct Parameters;

	if(open_infile(infile)){
		seq_params->header_file_names = read_header_files(infile, seq_params->num_header_files);
		read_hwe(infile, seq_params->hwe);
		read_debug(infile, seq_params->debug);
		read_iter(infile, seq_params->max_iterations);
		read_emthres(infile, seq_params->em_threshold);
		read_outfile(infile, seq_params->outfile);
		read_threads(infile, seq_params->threads);
		read_forcep(infile, seq_params->forcepval);
		read_mindepth(infile, seq_params->min_depth);
		read_error(infile, seq_params->error);
		read_pickinit(infile, seq_params->pickinit);
		if(seq_params->pickinit){
			read_stabilizer(infile, seq_params->stab, seq_params->naive_thres);
			read_hweswitch(infile, seq_params->hwecutoff);
			read_minindiv(infile, seq_params->min_indiv);
			if(seq_params->pickinit == 8){
				read_naive_params(infile, seq_params->max_naive_error, seq_params->max_indiv, seq_params->max_rdepth, seq_params->and_or);
			}
			seq_params->fixerror = 0;
		}
		else{
			read_fixerror(infile, seq_params->fixerror);
		}
//		read_inputformat(infile, seq_params->input);
		if(0 == seq_params->hwe){
			read_cov(infile, seq_params->cov_matrix);
			if(1 == seq_params->cov_matrix){
				read_minden(infile, seq_params->min_denominator);
				 seq_params->threads = 1;							// we only run one thread otherwise we get into mutex hell deciding which thread gets to print the various matrices.
			}
			read_init_geno_pval(infile, seq_params->prr, seq_params->prv);
			if(seq_params->prr >= 0){
				seq_params->pvv = 1.0 - seq_params->prr - seq_params->prv;
			}
			else{
				seq_params->pvv = -1;
			}
		}
		else{
			seq_params->cov_matrix = 0;
		}
	}
	else{
		fprintf(stderr,"Readparam::master_reader. Can not open %s. Bailing out.\n", filename);
		exit(0);
	}
	return(seq_params);
}

//**************************************

char** Readparam::read_header_files(ifstream& infile, int& val)
{
	long position;
	int i;
	char **headers = NULL;
	position = find_file_pos(infile, header_files, 0);
	if(position >= 0){
		position += (long)strlen(header_files);
		infile.seekg(position, ios::beg);
		infile >> val;
		headers = new char*[val];
		for(i = 0; i < val; i++){
			headers[i] = new char[max_name];
		}
		for(i = 0; i < val; i++){
			infile >> headers[i];
		}
	}
	else{
		cerr <<"Readparam::read_header_files. Unable to find "<<header_files<<". Bailing out!"<<endl;
		exit(0);
	}
	return(headers);
}

//**************************************

void Readparam::read_hwe(ifstream &infile, int &val)
{
	long position;
	position = find_file_pos(infile, hwe, 0);
	if(position >= 0){
		position += (long)strlen(hwe);
		infile.seekg(position, ios::beg);
		infile >> val;
	}
	else{
		val = 0;
	}
}

//**************************************

void Readparam::read_debug(ifstream &infile, int &val)
{
	long position;
	position = find_file_pos(infile, debug, 0);
	if(position >= 0){
		position += (long)strlen(debug);
		infile.seekg(position, ios::beg);
		infile >> val;
	}
	else{
		val = 0;
	}
}

//**************************************

void Readparam::read_iter(ifstream &infile, int &val)
{
	long position;
	position = find_file_pos(infile, maxiter, 0);
	if(position >= 0){
		position += (long)strlen(maxiter);
		infile.seekg(position, ios::beg);
		infile >> val;
		if(val < 0){
			val = 100;
		}
		if(val < 100){
			printf("Readparam::read_iter. Warning, at least 100 EM-iterations are recommended.\n");
		}
	}
	else{
		val = 100;
	}
}

//**************************************

void Readparam::read_cov(ifstream &infile, int &val)
{
	long position;
	position = find_file_pos(infile, cov_matrix, 0);
	if(position >= 0){
		position += (long)strlen(cov_matrix);
		infile.seekg(position, ios::beg);
		infile >> val;
	}
	else{
		val = 0;
	}
}

//**************************************

void Readparam::read_threads(ifstream &infile, int &val)
{
	long position;
	position = find_file_pos(infile, threads, 0);
	if(position >= 0){
		position += (long)strlen(threads);
		infile.seekg(position, ios::beg);
		infile >> val;
		if(val < 1){
			val = 1;
		}
	}
	else{
		val = 1;
	}
}

//**************************************

void Readparam::read_mindepth(ifstream &infile, int &val)
{
	long position;
	position = find_file_pos(infile, min_depth, 0);
	if(position >= 0){
		position += (long)strlen(min_depth);
		infile.seekg(position, ios::beg);
		infile >> val;
		if(val < 0){
			val = 0;
		}
	}
	else{
		val = 0;
	}
}

//**************************************

void Readparam::read_pickinit(ifstream &infile, int &val)
{
	long position;
	position = find_file_pos(infile, initval, 0);
	if(position >= 0){
		position += (long)strlen(initval);
		infile.seekg(position, ios::beg);
		infile >> val;
		if(val < 0){
			val = 0;
		}
	}
	else{
		val = 0;
	}
}

//**************************************

void Readparam::read_fixerror(ifstream &infile, int &val)
{
	long position;
	position = find_file_pos(infile, fixerror, 0);
	if(position >= 0){
		position += (long)strlen(fixerror);
		infile.seekg(position, ios::beg);
		infile >> val;
		if((val != 0) && (val != 1)){
			val = 0;
		}
	}
	else{
		val = 0;
	}
}

//**************************************

void Readparam::read_minindiv(ifstream &infile, int &val)
{
	long position;
	position = find_file_pos(infile, minindiv, 0);
	if(position >= 0){
		position += (long)strlen(minindiv);
		infile.seekg(position, ios::beg);
		infile >> val;
		if(val < 1){
			val = 15;
		}
	}
	else{
		val = 15;
	}
}

//**************************************

void Readparam::read_stabilizer(ifstream &infile, int &naive_min_reads, double &thres)
{
	long position;
	position = find_file_pos(infile, stabilizer, 0);
	thres = 0.12;
	if(position >= 0){
		position += (long)strlen(stabilizer);
		infile.seekg(position, ios::beg);
		infile >> naive_min_reads;
		if(naive_min_reads < 0){
			naive_min_reads = 0;
			thres = 0.12;
		}
		if(naive_min_reads > 127){
			naive_min_reads = 127;
		}
		infile >> thres;
		if((thres < 0.001) || (thres > 0.5)){
				thres = 0.12;
		}
	}
	else{
		naive_min_reads = 10;
		thres = 0.12;
	}
}

//**************************************

void Readparam::read_hweswitch(ifstream &infile, double &val)
{
	long position;
	position = find_file_pos(infile, hwechisq, 0);
	if(position >= 0){
		position += (long)strlen(hwechisq);
		infile.seekg(position, ios::beg);
		infile >> val;
		if((val < 0) || (val > 0.5)){
			val = 0.01;
		}
	}
	else{
		val = -1;
	}
}

//**************************************

void Readparam::read_init_geno_pval(ifstream &infile, double &prr, double &prv)
{
	long position;
	position = find_file_pos(infile, prrs, 0);
	if(position >= 0){
		position += (long)strlen(prrs);
		infile.seekg(position, ios::beg);
		infile >> prr;
		if((prr < 0) || (prr >= 1.0)){
			prr = 0.25;
		}
	}
	else{
		prr = -1;
	}
	position = find_file_pos(infile, prvs, 0);
	if(position >= 0){
		position += (long)strlen(prvs);
		infile.seekg(position, ios::beg);
		infile >> prv;
		if((prv < 0) || (prv >= 1.0)){
			prv = 0.5;
		}
	}
	else{
		prv = -1;
	}
	double sum = prv + prr;
	if(sum > 1.0){
		prv = -1;
		prr = -1;
	}
	if((prv < 0) || (prr < 0)){
		prv = -1;
		prr = -1;
	}
 
}

//**************************************

void Readparam::read_naive_params(ifstream &infile, double &err, int &rdepth, int &maxindiv, int &andor)
{
	long position;
	position = find_file_pos(infile, naiveparams, 0);
	if(position >= 0){
		position += (long)strlen(naiveparams);
		infile.seekg(position, ios::beg);
		infile >> err;
		infile >> rdepth;
		infile >> maxindiv;
		infile >> andor;
	}
	else{
		err = 0.01;
		rdepth = 10;
		maxindiv = 15;
		andor = 1;				// 0 == OR, 1 = AND
	}
}

//**************************************
/*
void Readparam::read_inputformat(ifstream &infile, int &val)
{
	long position;
	val = -1;
	position = find_file_pos(infile, infiletype, 0);
	if(position >= 0){
		position += (long)strlen(infiletype);
		infile.seekg(position, ios::beg);
		infile >> val;
	}
	if((val < 1) || (val > 2)){
		fprintf(stderr,"Readparam::read_inputformat. No input format specified. See documentation for details. Bailing out.\n");
		exit(0);
	}
}
*/
//**************************************

void Readparam::read_emthres(ifstream &infile, double &val)
{
	long position;
	position = find_file_pos(infile, em_thres, 0);
	if(position >= 0){
		position += (long)strlen(em_thres);
		infile.seekg(position, ios::beg);
		infile >> val;
	}
	else{
		val = 0.00000001;
	}
}

//**************************************

void Readparam::read_forcep(ifstream &infile, double &val)
{
	long position;
	position = find_file_pos(infile, force_p, 0);
	if(position >= 0){
		position += (long)strlen(force_p);
		infile.seekg(position, ios::beg);
		infile >> val;
	}
	else{
		val = -1.0;
	}
}

//**************************************

void Readparam::read_error(ifstream &infile, double &val)
{
	long position;
	position = find_file_pos(infile, errorsetting, 0);
	if(position >= 0){
		position += (long)strlen(errorsetting);
		infile.seekg(position, ios::beg);
		infile >> val;
		if((val < 0) || (val > 1)){
			fprintf(stderr,"Readparam::read_error. Initial error estimate set to 0.01.\n");
			val = 0.01;
		}
	}
	else{
		val = 0.01;
	}
}

//**************************************

void Readparam::read_minden(ifstream &infile, double &val)
{
	long position;
	position = find_file_pos(infile, min_den, 0);
	if(position >= 0){
		position += (long)strlen(min_den);
		infile.seekg(position, ios::beg);
		infile >> val;
	}
	else{
		val = 0.00001;
	}
}

//**************************************

void Readparam::read_outfile(ifstream &infile, char name[])
{
	long position;
	position = find_file_pos(infile, outfilename, 0);
	if(position >= 0){
		position += (long)strlen(outfilename);
		infile.seekg(position, ios::beg);
		infile >> name;
	}
	else{
		strcpy((char*)name, "outfile.txt");
	}
}

//**************************************

void Readparam::ptr_to_buffer(char *ptr, char buffer[], int size)
{
	int i = 0;
	while((i < (size-1)) && ('\0' != *(ptr+i))){
		buffer[i] = *(ptr+i);
		i++;
	}
	buffer[i] = '\0';
}

//**************************************

long Readparam::find_file_pos(ifstream &search_file, const char find_me[], long start)
{
	long result = -1;			// return value if string not found
	static long current = 0;
	long temp;
	int not_found = 0;
	int loop_around;
	int passed_eof = 0;

	char line[max_rdln];
	char *substrptr;
	if(start >= 0){
		current = start;
	}
	else{
		start = current;
	}
	loop_around = start;
	if((search_file.is_open()) && (strlen(find_me) > 0)){
		search_file.clear();
		search_file.seekg(start, ios::beg);
		while((result == -1) && (not_found == 0) /*(!(search_file.eof()))*/){
			if(search_file.eof()){
				search_file.clear();
				search_file.seekg(0, ios::beg);
				passed_eof++;
				current = 0;
			}
		    search_file.getline(line, max_rdln);
			substrptr = NULL;
			if(line[0] != '#'){					// indicates a comment. We don't consider this line
				substrptr = strstr(line, find_me);
			}
			if(substrptr == NULL){
#if defined (WIN32) 
				current += (long)strlen(line) + 2;  // +2 for '\0' on DOS format
#endif
#if defined (SOLARIS) || defined (LINUX) || defined (OSX)
				current += strlen(line) + 1;  // +1 for '\0' 
#endif
			}
			else{
				temp = (substrptr - line) / sizeof(char);
				result = current + temp;
			}
			if((passed_eof == 1) && (current > loop_around)){      // we looped around 1 full time and still didn't find the string we're lookin for, time to quit
				not_found = 1;
			}
		}
		search_file.clear();
	}
	else{
		cerr << "Param::find_file_posCould not open infile or empty string, exiting." << endl;
	}
	return result;                     
}

//**************************************

long Readparam::find_file_pos(fstream &search_file, const char find_me[], long start)
{
	long result = -1;			// return value if string not found
	static long current = 0;
	long temp;
	int not_found = 0;
	int loop_around;
	int passed_eof = 0;

	char line[max_rdln];
	char *substrptr;
	if(start >= 0){
		current = start;
	}
	else{
		start = current;
	}
	loop_around = start;
	if((search_file.is_open()) && (strlen(find_me) > 0)){
		search_file.clear();
		search_file.seekg(start, ios::beg);
		while((result == -1) && (not_found == 0) /*(!(search_file.eof()))*/){
			if(search_file.eof()){
				search_file.clear();
				search_file.seekg(0, ios::beg);
				passed_eof++;
				current = 0;
			}
		    search_file.getline(line, max_rdln);
			substrptr = NULL;
			if(line[0] != '#'){
				substrptr = strstr(line, find_me);
			}
			if(substrptr == NULL){
#if defined (WIN32) 
				current += strlen(line) + 2;  // +2 for '\0' on DOS format
#endif
#if defined (SOLARIS) || defined (LINUX) || defined (OSX)
				current += strlen(line) + 1;  // +1 for '\0' 
#endif
			}
			else{
			  temp = (substrptr - line) / sizeof(char);
                    result = current + temp;
			}
			if((passed_eof == 1) && (current > loop_around)){      // we looped around 1 full time and still didn't find the string we're lookin for, time to quit
				not_found = 1;
			}
		}
		search_file.clear();
	}
	else{
		cerr << "Param::find_file_posCould not open infile or empty string, exiting." << endl;
	}
	return result;                             
}

//**************************************
