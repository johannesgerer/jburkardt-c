# include <stdlib.h>
# include <limits.h>
# include <math.h>
# include <stdio.h>
# include <time.h>

# include "csparse.h"

typedef struct problem_struct
{
  cs *A ;
  cs *C ;
  int sym ;
  double *x ;
  double *b ;
  double *r ;
} problem ;

problem *get_problem ( FILE *f, double tol ) ;
int demo2 ( problem *Prob );
int demo3 ( problem *Prob );
problem *free_problem ( problem *Prob ) ;


/* 1 if A is square & upper tri., -1 if square & lower tri., 0 otherwise */
static int is_sym ( cs *A )
{
    int is_upper, is_lower, j, p, n = A->n, m = A->m, *Ap = A->p, *Ai = A->i ;
    if (m != n) return (0) ;
    is_upper = 1 ;
    is_lower = 1 ;
    for (j = 0 ; j < n ; j++)
    {
	for (p = Ap [j] ; p < Ap [j+1] ; p++)
	{
	    if (Ai [p] > j) is_upper = 0 ;
	    if (Ai [p] < j) is_lower = 0 ;
	}
    }
    return (is_upper ? 1 : (is_lower ? -1 : 0)) ;
}

/* true for off-diagonal entries */
static int dropdiag ( int i, int j, double aij, void *other ) 
{ 
  return (i != j);
}

/* C = A + triu(A,1)' */
static cs *make_sym ( cs *A )
{
    cs *AT, *C ;
    AT = cs_transpose (A, 1) ;		/* AT = A' */
    cs_fkeep (AT, &dropdiag, NULL) ;	/* drop diagonal entries from AT */
    C = cs_add (A, AT, 1, 1) ;		/* C = A+AT */
    cs_spfree (AT) ;
    return (C) ;
}

/* create a right-hand-side */
static void rhs (double *x, double *b, int m)
{
    int i ;
    for (i = 0 ; i < m ; i++) b [i] = 1 + ((double) i) / m ;
    for (i = 0 ; i < m ; i++) x [i] = b [i] ;
}

/* infinity-norm of x */
static double norm (double *x, int n)
{
    int i ;
    double normx = 0 ;
    for (i = 0 ; i < n ; i++) normx = CS_MAX (normx, fabs (x [i])) ;
    return (normx) ;
}

/* compute residual, norm(A*x-b,inf) / (norm(A,1)*norm(x,inf) + norm(b,inf)) */
static void resid (int ok, cs *A, double *x, double *b, double *r)
{
    int i, m, n ;
    if (!ok) { printf ("    (failed)\n") ; return ; }
    m = A->m ; n = A->n ;
    for (i = 0 ; i < m ; i++) r [i] = -b [i] ;	    /* r = -b */
    cs_gaxpy (A, x, r) ;			    /* r = r + A*x  */
    printf ("resid: %8.2e\n",
	norm (r,m) / ((n == 0) ? 1 : (cs_norm (A) * norm (x,n) + norm (b,m)))) ;
}

static double tic (void) 
{ 
  return (clock () / (double) CLOCKS_PER_SEC) ; 
}

static double toc (double t) 
{ 
  double s = tic (); 
  return (CS_MAX (0, s-t));
}

static void print_order (int order)
{
    switch (order)
    {
	case -1: printf ("natural    ") ; break ;
	case  0: printf ("amd(A+A')  ") ; break ;
	case  1: printf ("amd(S'*S)  ") ; break ;
	case  2: printf ("amd(A'*A)  ") ; break ;
    }
}

/* read a problem from a file */
problem *get_problem (FILE *f, double tol)
{
    cs *T, *A, *C ;
    int sym, m, n, mn, nz1, nz2 ;
    problem *Prob ;
    Prob = cs_calloc (1, sizeof (problem)) ;
    if (!Prob) return (NULL) ;
    T = cs_load (f) ;			/* load triplet matrix T from a file */
    Prob->A = A = cs_triplet (T) ;	/* A = compressed-column form of T */
    cs_spfree (T) ;			/* clear T */
    if (!cs_dupl (A)) return (free_problem (Prob)) ; /* sum up duplicates */
    Prob->sym = sym = is_sym (A) ;	/* determine if A is symmetric */
    m = A->m ; n = A->n ;
    mn = CS_MAX (m,n) ;
    nz1 = A->p [n] ;
    cs_dropzeros (A) ;			/* drop zero entries */
    nz2 = A->p [n] ;
    if (tol > 0) cs_droptol (A, tol) ;	/* drop tiny entries (just to test) */
    Prob->C = C = sym ? make_sym (A) : A ;  /* C = A + triu(A,1)', or C=A */
    if (!C) return (free_problem (Prob)) ;
    printf ("\n--- Matrix: %d-by-%d, nnz: %d (sym: %d: nnz %d), norm: %8.2e\n",
	    m, n, A->p [n], sym, sym ? C->p [n] : 0, cs_norm (C)) ;
    if (nz1 != nz2) printf ("zero entries dropped: %d\n", nz1 - nz2) ;
    if (nz2 != A->p [n]) printf ("tiny entries dropped: %d\n", nz2 - A->p [n]) ;
    Prob->b = cs_malloc (mn, sizeof (double)) ;
    Prob->x = cs_malloc (mn, sizeof (double)) ;
    Prob->r = cs_malloc (mn, sizeof (double)) ;
    return ((!Prob->b || !Prob->x || !Prob->r) ? free_problem (Prob) : Prob) ;
}

/* free a problem */
problem *free_problem (problem *Prob)
{
    if (!Prob) return (NULL) ;
    cs_spfree (Prob->A) ;
    if (Prob->sym) cs_spfree (Prob->C) ;
    cs_free (Prob->b) ;
    cs_free (Prob->x) ;
    cs_free (Prob->r) ;
    return (cs_free (Prob)) ;
}

/* solve a linear system using Cholesky, LU, and QR, with various orderings */
int demo2 (problem *Prob)
{
    cs *A, *C ;
    double *b, *x, *r,  t, tol ;
    int k, m, n, ok, order, nb, ns, *R, *S, *rr, sprank ;
    csd *D ;
    if (!Prob) return (0) ;
    A = Prob->A ; C = Prob->C ; b = Prob->b ; x = Prob->x ; r = Prob->r ;
    m = A->m ; n = A->n ;
    tol = Prob->sym ? 0.001 : 1 ;		/* partial pivoting tolerance */
    D = cs_dmperm (C) ;				/* dmperm analysis */
    if (!D) return (0) ;
    nb = D->nb ; R = D->R ; S = D->S ; rr = D->rr ;
    sprank = rr [3] ;
    for (ns = 0, k = 0 ; k < nb ; k++)
    {
	ns += ((R [k+1] == R [k]+1) && (S [k+1] == S [k]+1)) ;
    }
    printf ("blocks: %d singletons: %d structural rank: %d\n", nb, ns, sprank) ;
    cs_dfree (D) ;
    for (order = -1 ; order <= 2 ; order += 3)	/* natural and amd(A'*A) */
    {
	if (order == -1 && m > 1000) continue ;
	printf ("QR   ") ;
	print_order (order) ;
	rhs (x, b, m) ;				/* compute right-hand-side */
	t = tic () ;
	ok = cs_qrsol (C, x, order) ;		/* min norm(Ax-b) with QR */
	printf ("time: %8.2f ", toc (t)) ;
	resid (ok, C, x, b, r) ;		/* print residual */
    }
    if (m != n || sprank < n) return (1) ;	/* return if rect. or singular*/
    for (order = -1 ; order <= 2 ; order++)	/* try all orderings */
    {
	if (order == -1 && m > 1000) continue ;
	printf ("LU   ") ;
	print_order (order) ;
	rhs (x, b, m) ;				/* compute right-hand-side */
	t = tic () ;
	ok = cs_lusol (C, x, order, tol) ;	/* solve Ax=b with LU */
	printf ("time: %8.2f ", toc (t)) ;
	resid (ok, C, x, b, r) ;		/* print residual */
    }
    if (!Prob->sym) return (1) ;
    for (order = -1 ; order <= 0 ; order++)	/* natural and amd(A+A') */
    {
	if (order == -1 && m > 1000) continue ;
	printf ("Chol ") ;
	print_order (order) ;
	rhs (x, b, m) ;				/* compute right-hand-side */
	t = tic () ;
	ok = cs_cholsol (C, x, order) ;		/* solve Ax=b with Cholesky */
	printf ("time: %8.2f ", toc (t)) ;
	resid (ok, C, x, b, r) ;		/* print residual */
    }
    return (1) ;
} 

/* free workspace for demo3 */
static int done3 (int ok, css *S, csn *N, double *y, cs *W, cs *E, int *P)
{
    cs_sfree (S) ;
    cs_nfree (N) ;
    cs_free (y) ;
    cs_spfree (W) ;
    cs_spfree (E) ;
    cs_free (P) ;
    return (ok) ;
}

/* Cholesky update/downdate */
int demo3 (problem *Prob)
{
    cs *A, *C, *W = NULL, *WW, *WT, *E = NULL, *W2 ;
    int n, k, *Li, *Lp, *Wi, *Wp, p, p2, *P = NULL, ok ;
    double *b, *x, *r, *y = NULL, *Lx, *Wx, s,  t, t1 ;
    css *S = NULL ;
    csn *N = NULL ;
    if (!Prob || !Prob->sym || Prob->A->n == 0) return (0) ;
    A = Prob->A ; C = Prob->C ; b = Prob->b ; x = Prob->x ; r = Prob->r ;
    n = A->n ;
    if (!Prob->sym || n == 0) return (1) ;
    rhs (x, b, n) ;				/* compute right-hand-side */
    printf ("\nchol then update/downdate ") ;
    print_order (0) ;
    y = cs_malloc (n, sizeof (double)) ;
    t = tic () ;
    S = cs_schol (C, 0) ;			/* symbolic Cholesky */
    printf ("\nsymbolic chol time %8.2f\n", toc (t)) ;
    t = tic () ;
    N = cs_chol (C, S) ;			/* numeric Cholesky */
    printf ("numeric  chol time %8.2f\n", toc (t)) ;
    if (!S || !N || !y) return (done3 (0, S, N, y, W, E, P)) ;
    t = tic () ;
    cs_ipvec (n, S->Pinv, b, y) ;		/* y = P*b */
    cs_lsolve (N->L, y) ;			/* y = L\y */
    cs_ltsolve (N->L, y) ;			/* y = L'\y */
    cs_pvec (n, S->Pinv, y, x) ;		/* x = P'*y */
    printf ("solve    chol time %8.2f\n", toc (t)) ;
    printf ("original: ") ;
    resid (1, C, x, b, r) ;			/* print residual */
    k = n/2 ;					/* construct W  */
    W = cs_spalloc (n, 1, n, 1, 0) ;
    if (!W) return (done3 (0, S, N, y, W, E, P)) ;
    Lp = N->L->p ; Li = N->L->i ; Lx = N->L->x ;
    Wp = W->p ; Wi = W->i ; Wx = W->x ;
    Wp [0] = 0 ;
    p = Lp [k] ;
    Wp [1] = Lp [k+1] - p ;
    s = Lx [p] ;
    srand (1) ;
    for ( ; p < Lp [k+1] ; p++)
    {
	p2 = p - Lp [k] ;
	Wi [p2] = Li [p] ;
	Wx [p2] = s * rand () / ((double) RAND_MAX) ;
    }
    t = tic () ;
    ok = cs_updown (N->L, +1, W, S->parent) ;	/* update: L*L'+W*W' */
    t1 = toc (t) ;
    printf ("update:   time: %8.2f\n", t1) ;
    if (!ok) return (done3 (0, S, N, y, W, E, P)) ;
    t = tic () ;
    cs_ipvec (n, S->Pinv, b, y) ;		/* y = P*b */
    cs_lsolve (N->L, y) ;			/* y = L\y */
    cs_ltsolve (N->L, y) ;			/* y = L'\y */
    cs_pvec (n, S->Pinv, y, x) ;		/* x = P'*y */
    t = toc (t) ;
    P = cs_pinv (S->Pinv, n) ;
    W2 = cs_permute (W, P, NULL, 1) ;		/* E = C + (P'W)*(P'W)' */
    WT = cs_transpose (W2,1) ;
    WW = cs_multiply (W2, WT) ;
    cs_spfree (WT) ;
    cs_spfree (W2) ;
    E = cs_add (C, WW, 1, 1) ;
    cs_spfree (WW) ;
    if (!E || !P) return (done3 (0, S, N, y, W, E, P)) ;
    printf ("update:   time: %8.2f (incl solve) ", t1+t) ;
    resid (1, E, x, b, r) ;			/* print residual */
    cs_nfree (N) ;				/* clear N */
    t = tic () ;
    N = cs_chol (E, S) ;			/* numeric Cholesky */
    if (!N) return (done3 (0, S, N, y, W, E, P)) ;
    cs_ipvec (n, S->Pinv, b, y) ;		/* y = P*b */
    cs_lsolve (N->L, y) ;			/* y = L\y */
    cs_ltsolve (N->L, y) ;			/* y = L'\y */
    cs_pvec (n, S->Pinv, y, x) ;		/* x = P'*y */
    t = toc (t) ;
    printf ("rechol:   time: %8.2f (incl solve) ", t) ;
    resid (1, E, x, b, r) ;			/* print residual */
    t = tic () ;
    ok = cs_updown (N->L, -1, W, S->parent) ;	/* downdate: L*L'-W*W' */
    t1 = toc (t) ;
    if (!ok) return (done3 (0, S, N, y, W, E, P)) ;
    printf ("downdate: time: %8.2f\n", t1) ;
    t = tic () ;
    cs_ipvec (n, S->Pinv, b, y) ;		/* y = P*b */
    cs_lsolve (N->L, y) ;			/* y = L\y */
    cs_ltsolve (N->L, y) ;			/* y = L'\y */
    cs_pvec (n, S->Pinv, y, x) ;		/* x = P'*y */
    t = toc (t) ;
    printf ("downdate: time: %8.2f (incl solve) ", t1+t) ;
    resid (1, C, x, b, r) ;			/* print residual */
    return (done3 (1, S, N, y, W, E, P)) ;
} 
