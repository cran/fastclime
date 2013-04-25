#include <stdio.h>
#include "heap.h"

void hfall(
        int     heapnum, 
        int     *key, 
        int     *iheap, 
        int     *heap, 
        int     cur
)
{
        int child, child1;
	int heapcur0 = heap[cur];
	int keyheapcur0 = key[heapcur0];

        child = 2*cur;
        while (child <= heapnum) {
		child1 = child+1;
                if (child < heapnum && key[heap[child1]] < key[heap[child]]) {
                        child=child1;
		}
                if (keyheapcur0 <= key[heap[child]]) { break; }
		heap[cur] = heap[child];
		iheap[heap[cur]] = cur;
		cur = child;
		child = 2*cur;
        }
	heap[cur] = heapcur0;
	iheap[heap[cur]] = cur;
}

void hrise(
        int     *key, 
        int     *iheap, 
        int     *heap, 
        int     cur
)
{
        int parent;
	int heapcur0 = heap[cur];
	int keyheapcur0 = key[heapcur0];

        parent = cur/2;
        while (parent > 0) {
                if (key[heap[parent]] <= keyheapcur0) { break; }
		heap[cur] = heap[parent];
		iheap[heap[cur]] = cur;
		cur = parent;
		parent = cur/2;
        }
	heap[cur] = heapcur0;
	iheap[heap[cur]] = cur;
}



