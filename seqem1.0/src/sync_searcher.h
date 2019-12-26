#ifndef __SYNC_SEARCH_H__
#define __SYNC_SEARCH_H__

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


struct Sync_searcher{
	int items;											// number of sorted lists that need to be searched against the one stored in data[0]
	struct Node *head;									// Linked list
	struct Node *current;
	struct Node *previous;
	int list_size;										// size of linked list that contains CCDSs not in list stored at data[0]
	int next;											// data[next] is the next one to be processed. next gets initialized to 1 (NOT 0) since [0] is the one we compare/search against
	pthread_mutex_t sync_searcher_mutex;				// 1 producer of data, reading from file(s)
};

#endif
