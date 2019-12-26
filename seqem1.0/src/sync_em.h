#ifndef __SYNCEM_H__
#define __SYNCEM_H__

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

struct Sync_em{
	int size;					// total items to be processed
	int next;					// next item to process
	pthread_mutex_t sync_em_mutex;
};

#endif
