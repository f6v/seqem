#ifndef __SYNC_ANA_H__
#define __SYNC_ANA_H__

#ifdef WIN32
#include "pthread.h"
#include "sched.h"
#include "semaphore.h"
#include <windows.h>
#endif
#ifdef SOLARIS
#include <pthread.h>			// define _REENTRANT in makefile
#include <string.h>
#endif                          
#ifdef LINUX
#include <pthread.h>			// define _REENTRANT in makefile
#include <string.h>
#endif
#ifdef AIX
#include <pthread.h>			// define _REENTRANT in makefile
#include <string.h>
#endif


struct Sync_ana{
	int items;										// total number of SNPs
	int threads;									// number of threads
	int current;									// current = the next SNP to process
	pthread_mutex_t sync_ana_mutex;					// mutex to update current to next one to process. Once current >= item -> stop
};

#endif
