#
#include "defs.h"

int ntry, totalsearch;

/******************************************************************************/

ELinitialize()

/******************************************************************************/
{
int i;
	freeinit(&hfl, sizeof **ELhash);
	ELhashsize = 2 * sqrt_nsites;
	ELhash = (struct Halfedge **) myalloc ( sizeof *ELhash * ELhashsize);
	for(i=0; i<ELhashsize; i +=1) ELhash[i] = (struct Halfedge *)NULL;
	ELleftend = HEcreate( (struct Edge *)NULL, 0);
	ELrightend = HEcreate( (struct Edge *)NULL, 0);
	ELleftend -> ELleft = (struct Halfedge *)NULL;
	ELleftend -> ELright = ELrightend;
	ELrightend -> ELleft = ELleftend;
	ELrightend -> ELright = (struct Halfedge *)NULL;
	ELhash[0] = ELleftend;
	ELhash[ELhashsize-1] = ELrightend;
}
/******************************************************************************/

struct Halfedge *HEcreate(e, pm)
struct Edge *e;
int pm;

/******************************************************************************/
{
struct Halfedge *answer;
	answer = (struct Halfedge *) getfree(&hfl);
	answer -> ELedge = e;
	answer -> ELpm = pm;
	answer -> PQnext = (struct Halfedge *) NULL;
	answer -> vertex = (struct Site *) NULL;
	answer -> ELrefcnt = 0;
	return(answer);
}
/******************************************************************************/

ELinsert(lb, new)
struct	Halfedge *lb, *new;

/******************************************************************************/
{
	new -> ELleft = lb;
	new -> ELright = lb -> ELright;
	(lb -> ELright) -> ELleft = new;
	lb -> ELright = new;
}
/******************************************************************************/

struct Halfedge *ELgethash(b)
int b;
/******************************************************************************/
/* Get entry from hash table, pruning any deleted nodes */
{
struct Halfedge *he;

	if(b<0 || b>=ELhashsize) return((struct Halfedge *) NULL);
	he = ELhash[b]; 
	if (he == (struct Halfedge *) NULL || 
	    he -> ELedge != (struct Edge *) DELETED ) return (he);

/* Hash table points to deleted half edge.  Patch as necessary. */
	ELhash[b] = (struct Halfedge *) NULL;
	if ((he -> ELrefcnt -= 1) == 0) makefree(he, &hfl);
	return ((struct Halfedge *) NULL);
}	
/******************************************************************************/
struct Halfedge *ELleftbnd(p)
struct Point *p;
/******************************************************************************/
{
int i, bucket;
struct Halfedge *he;

/* Use hash table to get close to desired halfedge */
	bucket = (p->x - xmin)/deltax * ELhashsize;
	if(bucket<0) bucket =0;
	if(bucket>=ELhashsize) bucket = ELhashsize - 1;
	he = ELgethash(bucket);
	if(he == (struct Halfedge *) NULL)
	{   for(i=1; 1 ; i += 1)
	    {	if ((he=ELgethash(bucket-i)) != (struct Halfedge *) NULL) break;
		if ((he=ELgethash(bucket+i)) != (struct Halfedge *) NULL) break;
	    };
	totalsearch += i;
	};
	ntry += 1;
/* Now search linear list of halfedges for the corect one */
	if (he==ELleftend  || (he != ELrightend && right_of(he,p)))
	{do {he = he -> ELright;} while (he!=ELrightend && right_of(he,p));
	 he = he -> ELleft;
	}
	else 
	do {he = he -> ELleft;} while (he!=ELleftend && !right_of(he,p));

/* Update hash table and reference counts */
	if(bucket > 0 && bucket <ELhashsize-1)
	{	if(ELhash[bucket] != (struct Halfedge *) NULL) 
			ELhash[bucket] -> ELrefcnt -= 1;
		ELhash[bucket] = he;
		ELhash[bucket] -> ELrefcnt += 1;
	};
	return (he);
}
/******************************************************************************/
	
ELdelete(he)
struct Halfedge *he;

/******************************************************************************/
/* 

  This delete routine can't reclaim node, since pointers from hash table may be present.   
*/
{
	(he -> ELleft) -> ELright = he -> ELright;
	(he -> ELright) -> ELleft = he -> ELleft;
	he -> ELedge = (struct Edge *)DELETED;
}
/******************************************************************************/

struct Halfedge *ELright(he)
struct Halfedge *he;

/******************************************************************************/
{
	return (he -> ELright);
}
/******************************************************************************/

struct Halfedge *ELleft(he)
struct Halfedge *he;

/******************************************************************************/
{
	return (he -> ELleft);
}
/******************************************************************************/

struct Site *leftreg(he)
struct Halfedge *he;

/******************************************************************************/
{
	if(he -> ELedge == (struct Edge *)NULL) return(bottomsite);
	return( he -> ELpm == le ? 
		he -> ELedge -> reg[le] : he -> ELedge -> reg[re]);
}

struct Site *rightreg(he)
struct Halfedge *he;
{
	if(he -> ELedge == (struct Edge *)NULL) return(bottomsite);
	return( he -> ELpm == le ? 
		he -> ELedge -> reg[re] : he -> ELedge -> reg[le]);
}


