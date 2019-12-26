#ifndef __CCDSI_H__
#define __CCDSI_H__

#include <iostream>
#include <string>
#include <stdio.h>
#include "max.h"


class CCDSI{										// one line from a header file. Each CCDS represents a single read from one person
public:
	CCDSI();
	~CCDSI();
	CCDSI & operator=(const CCDSI &rhs);
	int operator==(const CCDSI&);
	int operator!=(const CCDSI&);
	int operator>(const CCDSI&);
	int operator>=(const CCDSI&);
	int operator<(const CCDSI&);
	int operator<=(const CCDSI&);

	unsigned int location;
	unsigned int depth;
	unsigned int variants;
	unsigned char variant_hash;

};


#endif
