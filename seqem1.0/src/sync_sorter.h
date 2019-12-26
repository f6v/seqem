#ifndef __SYNC_SORT_H__
#define __SYNC_SORT_H__

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


struct Sync_sorter{
	int items;
	int ready;									// number of items ready to sort
	int done;									// 0 if not done. 1 = done, time for sorter thread to exit
	int search_start;							// where to start search for new item
	pthread_mutex_t sync_sorter_mutex;			// 1 producer of data, reading from file(s)
	pthread_cond_t cond_waiting_consumer;		// consumer/sorter has to wait until data has been read from file
	int *list;									// list[items]. -1 = currently in use, 0 = not read yet, 1 = read, ready to be sorted, 2 = sorted data
};

#endif
