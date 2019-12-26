#include "ccds.h"


//*********************************

CCDS::CCDS()
{
	ccds_name = 0;
//	reference = NULL;
//	variant = NULL;
	presentin = NULL;
	genotype = NULL;
	geno_prop = NULL;
	individuals = 0;
	start = -1;
	stop = -1;
	chr = 0;
	num_reads = 0;
	num_var = 0;
	proportion = -1.0;
	prr = -1;
	prv = -1;
	pvv = -1;
}

//*********************************

CCDS::~CCDS()
{
	int i;
	/*
	if(NULL != reference){
		delete [] reference;
	}
	if(NULL != variant){
		delete [] variant;
	}
	*/
	if(NULL != presentin){
		delete [] presentin;
	}
	if(NULL != genotype){
		delete [] genotype;
	}
	if(NULL != geno_prop){
		for(i = 0; i < 3; i++){
			delete [] geno_prop[i];
		}
		delete [] geno_prop;
	}
}

//*********************************

CCDS & CCDS::operator=(const CCDS &rhs)
{
	int length;
	int i;
	if(&rhs != this){													// check if they're the same
		start = rhs.start;
		stop = rhs.stop;
		proportion = rhs.proportion;
		num_reads = rhs.num_reads;
		num_var = rhs.num_var;
		if(NULL != geno_prop){
			for(i = 0; i < 3; i++){
				delete [] geno_prop[i];
			}
			delete [] geno_prop;
			geno_prop = NULL;
		}
		individuals = rhs.individuals;
		error = rhs.error;
		chr = rhs.chr;
		iterations = rhs.iterations;
		p = rhs.p;
		prr = rhs.prr;
		prv = rhs.prv;
		pvv = rhs.pvv;
		stat = rhs.stat;
		ccds_name = rhs.ccds_name;
		/*
		if(reference != NULL){
			delete [] reference;
			reference = NULL;
		}
		if(variant != NULL){
			delete [] variant;
			variant = NULL;
		}
		*/
		if(presentin != NULL){
			delete [] presentin;
			presentin = NULL;
		}
		if(genotype != NULL){
			delete [] genotype;
			genotype = NULL;
		}
		/*
		if(NULL != rhs.reference){
			length = (int)(strlen(rhs.reference));
			if(0 < length){
				reference = new char[length+2];
				strcpy(reference, rhs.reference);
			}
		}
		if(NULL != rhs.variant){
			length = (int)(strlen(rhs.variant));
			if(0 < length){
				variant = new char[length+2];
				strcpy(variant, rhs.variant);
			}
		}
		*/
		if((0 < rhs.individuals) && (NULL != rhs.genotype)){
			genotype = new unsigned char[rhs.individuals];
			memcpy((void*)genotype, (void*)rhs.genotype, individuals*sizeof(unsigned char));
		}
		if((0 < rhs.individuals) && (NULL != rhs.genotype)){
			geno_prop = new float*[3];
			for(i = 0; i < 3; i++){
				geno_prop[i] = new float[rhs.individuals];
				memcpy((void*)geno_prop[i], rhs.geno_prop[i], rhs.individuals*sizeof(float));
			}
		}
	}
	return(*this);
}

//*********************************

int CCDS::operator==(const CCDS &other)
{
	/*
	int result = 0; 
	if(other.ccds_name == ccds_name){
		result = 1;
	}
	return(result);
	*/
	return(other.ccds_name == ccds_name);
}

//*********************************

int CCDS::operator!=(const CCDS &other)
{
	/*
	int result = 1;
	if(other.ccds_name == ccds_name){
		result = 0;
	}
	return(result);
	*/
	return(other.ccds_name != ccds_name);
}

//*********************************

int CCDS::operator>(const CCDS &other)
{
	/*
	int result = 0;
	if(ccds_name > other.ccds_name){
		result = 1;
	}
	return(result);
	*/
	return(ccds_name > other.ccds_name);
}

//*********************************

int CCDS::operator>=(const CCDS &other)
{
	/*
	int result;
	result = strcmp(ccds_name, other.ccds_name);
	if(result >= 0){
		result = 1;
	}
	else{
		result = 0;
	}	
	return(result);
	*/
	return(ccds_name >= other.ccds_name);
}

//*********************************

int CCDS::operator<(const CCDS &other)
{
	/*
	int result;
	result = strcmp(ccds_name, other.ccds_name);
	if(result < 0){
		result = 1;
	}
	else{
		result = 0;
	}
	return(result);
	*/
	return(ccds_name < other.ccds_name);
}

//*********************************

int CCDS::operator<=(const CCDS &other)
{
	/*
	int result;
	result = strcmp(ccds_name, other.ccds_name);
	if(result <= 0){
		result = 1;
	}
	else{
		result = 0;
	}
	return(result);
	*/
	return(ccds_name <= other.ccds_name);
}

//*********************************
