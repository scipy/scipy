/*
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 *
 */
/** Precision-independent memory-related routines.
    (Shared by [sdcz]memory.c) **/

#include "util.h"

/*
 * Set up pointers for integer working arrays.
 */
void
SetIWork(int m, int n, int panel_size, int *iworkptr, int **segrep,
	 int **parent, int **xplore, int **repfnz, int **panel_lsub,
	 int **xprune, int **marker)
{
    *segrep = iworkptr;
    *parent = iworkptr + m;
    *xplore = *parent + m;
    *repfnz = *xplore + m;
    *panel_lsub = *repfnz + panel_size * m;
    *xprune = *panel_lsub + panel_size * m;
    *marker = *xprune + n;
    ifill (*repfnz, m * panel_size, EMPTY);
    ifill (*panel_lsub, m * panel_size, EMPTY);
}


void
copy_mem_int(int howmany, void *old, void *new)
{
    register int i;
    int *iold = old;
    int *inew = new;
    for (i = 0; i < howmany; i++) inew[i] = iold[i];
}


void
user_bcopy(char *src, char *dest, int bytes)
{
    char *s_ptr, *d_ptr;

    s_ptr = src + bytes - 1;
    d_ptr = dest + bytes - 1;
    for (; d_ptr >= dest; --s_ptr, --d_ptr ) *d_ptr = *s_ptr;
}



int *intMalloc(int n)
{
    int *buf;
    buf = (int *) SUPERLU_MALLOC(n * sizeof(int));
    if ( !buf ) {
	ABORT("SUPERLU_MALLOC fails for buf in intMalloc()");
    }
    return (buf);
}

int *intCalloc(int n)
{
    int *buf;
    register int i;
    buf = (int *) SUPERLU_MALLOC(n * sizeof(int));
    if ( !buf ) {
	ABORT("SUPERLU_MALLOC fails for buf in intCalloc()");
    }
    for (i = 0; i < n; ++i) buf[i] = 0;
    return (buf);
}



#if 0
check_expanders()
{
    int p;
    printf("Check expanders:\n");
    for (p = 0; p < NO_MEMTYPE; p++) {
	printf("type %d, size %d, mem %d\n",
	       p, expanders[p].size, (int)expanders[p].mem);
    }

    return 0;
}


StackInfo()
{
    printf("Stack: size %d, used %d, top1 %d, top2 %d\n",
	   stack.size, stack.used, stack.top1, stack.top2);
    return 0;
}



PrintStack(char *msg, GlobalLU_t *Glu)
{
    int i;
    int *xlsub, *lsub, *xusub, *usub;

    xlsub = Glu->xlsub;
    lsub  = Glu->lsub;
    xusub = Glu->xusub;
    usub  = Glu->usub;

    printf("%s\n", msg);
    
/*    printf("\nUCOL: ");
    for (i = 0; i < xusub[ndim]; ++i)
	printf("%f  ", ucol[i]);

    printf("\nLSUB: ");
    for (i = 0; i < xlsub[ndim]; ++i)
	printf("%d  ", lsub[i]);

    printf("\nUSUB: ");
    for (i = 0; i < xusub[ndim]; ++i)
	printf("%d  ", usub[i]);

    printf("\n");*/
    return 0;
}   
#endif



