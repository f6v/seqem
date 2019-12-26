#include "sorter_t.h"

//*************************

Sorter_t::Sorter_t()
{

}

//*************************

Sorter_t::~Sorter_t()
{

}

//*************************

void Sorter_t::run(int id)
{
	Sort <CCDSI> sorter;
	int i;
	int finished;
	while(sort_syncher->done < sort_syncher->items){
		pthread_mutex_lock(&sort_syncher->sync_sorter_mutex);
		while((sort_syncher->ready == 0) && (sort_syncher->done < sort_syncher->items)){	// wait until some data becomes available
			pthread_cond_wait(&sort_syncher->cond_waiting_consumer, &sort_syncher->sync_sorter_mutex);
		}
		if(sort_syncher->done == sort_syncher->items){						// check if we might be done sorting
			pthread_mutex_unlock(&sort_syncher->sync_sorter_mutex);			// if YES, unlock mutex and return
			return;
		}
		i = 0;
		while((sort_syncher->list[i] != 1) && (i < sort_syncher->items)){	// to find an unsorted item (marked by 1) linearly search list from beginning
			i++;															// increment list position
		}
		if(i == sort_syncher->items){
			fprintf(stderr,"Sorter_t::run. Can not find item to be sorted. Bailing out.\n");
			exit(0);
		}
		sort_syncher->list[i] = -1;									// indicating this item is being worked on
		sort_syncher->ready--;
		pthread_mutex_unlock(&sort_syncher->sync_sorter_mutex);
		sorter.init(data[i].total_reads, data[i].ccds_data);		// sort item
		sorter.qsort();
		pthread_mutex_lock(&sort_syncher->sync_sorter_mutex);		// lock mutex
		sort_syncher->list[i] = 2;									// mark item as finished/sorted
		sort_syncher->done++;
		pthread_mutex_unlock(&sort_syncher->sync_sorter_mutex);		// unlock mutex
	}
	pthread_cond_broadcast(&sort_syncher->cond_waiting_consumer);	// very important, wake up any sleeping threads that might still be waiting
}

//*************************
