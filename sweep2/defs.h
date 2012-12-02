#ifndef NULL
#define NULL 0
#endif
#define DELETED -2

int triangulate, sorted, plot, debug;

struct	Freenode	{
struct	Freenode	*nextfree;
};
struct	Freelist	{
struct	Freenode	*head;
int			nodesize;
};
char *getfree();
char *malloc();
char *myalloc();

float xmin, xmax, ymin, ymax, deltax, deltay;


struct Point	{
float x,y;
};

/* structure used both for sites and for vertices */
struct Site	{
struct	Point	coord;
int		sitenbr;
int		refcnt;
};


struct	Site	*sites;
int		nsites;
int		siteidx;
int		sqrt_nsites;
int		nvertices;
struct 	Freelist sfl;
struct	Site	*bottomsite;


struct Edge	{
float		a,b,c;
struct	Site 	*ep[2];
struct	Site	*reg[2];
int		edgenbr;
};
#define le 0
#define re 1
int nedges;
struct	Freelist efl;

int has_endpoint(),right_of();
struct Site *intersect();
float dist();
struct Point PQ_min();
struct Halfedge *PQextractmin();
struct Edge *bisect();

struct Halfedge {
struct Halfedge	*ELleft, *ELright;
struct Edge	*ELedge;
int		ELrefcnt;
char		ELpm;
struct	Site	*vertex;
float		ystar;
struct	Halfedge *PQnext;
};

struct   Freelist	hfl;
struct	Halfedge *ELleftend, *ELrightend;
int 	ELhashsize;
struct	Halfedge **ELhash;
struct	Halfedge *HEcreate(), *ELleft(), *ELright(), *ELleftbnd();
struct	Site *leftreg(), *rightreg();


int PQhashsize;
struct	Halfedge *PQhash;
struct	Halfedge *PQfind();
int PQcount;
int PQmin;
int PQempty();


