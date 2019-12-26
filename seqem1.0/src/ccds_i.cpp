#include "ccds_i.h"

//****************************

CCDSI::CCDSI()
{
	location = 0;
	depth = 0;
	variants = 0;
	variant_hash = 0;
}

//****************************

CCDSI::~CCDSI()
{

}

//****************************

CCDSI & CCDSI::operator=(const CCDSI &rhs)
{
	location = rhs.location;
	depth = rhs.depth;
	variant_hash = rhs.variant_hash;
	variants = rhs.variants;
	return(*this);
}

//****************************

int CCDSI::operator==(const CCDSI &other)
{
	return(location == other.location);
}

//****************************

int CCDSI::operator!=(const CCDSI &other)
{
	return(location != other.location);
}

//****************************

int CCDSI::operator>(const CCDSI &other)
{
	return(location > other.location);
}

//****************************

int CCDSI::operator>=(const CCDSI &other)
{
	return(location >= other.location);
}

//****************************

int CCDSI::operator<(const CCDSI &other)
{
	return(location < other.location);
}

//****************************

int CCDSI::operator<=(const CCDSI &other)
{
	return(location <= other.location);
}

//****************************

//****************************

//****************************
