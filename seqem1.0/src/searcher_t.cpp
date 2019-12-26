#include "searcher_t.h"

//*************************

Searcher_t::Searcher_t()
{

}

//*************************

Searcher_t::~Searcher_t()
{

}

//*************************

void Searcher_t::run(int id)
{
	int i,j;
	Binsearch <CCDSI> searcher;
	int position;
	while(search_syncher->next < search_syncher->items){
		pthread_mutex_lock(&search_syncher->sync_searcher_mutex);
		if(search_syncher->next == search_syncher->items){
			pthread_mutex_unlock(&search_syncher->sync_searcher_mutex);
			return;
		}
		i = search_syncher->next;
		search_syncher->next++;
		pthread_mutex_unlock(&search_syncher->sync_searcher_mutex);
		for(j = 0; j < data[i].total_reads; j++){
			position = searcher.search(data[0].ccds_data, data[i].ccds_data[j], 0, data[0].total_reads-1);
			if(-1 == position){												// if item not in reference list, add to linked list
				pthread_mutex_lock(&search_syncher->sync_searcher_mutex);
				search_syncher->current = new struct Node;
				search_syncher->current->item = data[i].ccds_data[j];
				search_syncher->current->next = NULL;
				search_syncher->list_size++;
				if(NULL == search_syncher->head){											// empty list
					search_syncher->head = search_syncher->current;
					search_syncher->previous = search_syncher->current;
				}
				else{														// non-empty list
					search_syncher->previous->next = search_syncher->current;
					search_syncher->previous = search_syncher->current;
				}
				pthread_mutex_unlock(&search_syncher->sync_searcher_mutex);
			} // if(-1 == position)
		} // for(j = 0; j < data[i].total_reads; j++)
	} // while(search_syncher->next < search_syncher->items)
}

//*************************
