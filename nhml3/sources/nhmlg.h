#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <limits.h>
#include <float.h>

#define TRUE 1
#define FALSE 0
#define MAXNSP 100
#define MAXDATASET 100
#define MAXLNAME 20
#define MAXLLINE 10000
#define lmot (8*sizeof(int))

#define seed 1234
#define deltax 1.e-4

#define nb_param_type 7
#define lmin 1.e-6
#define lmax 9999.
#define GCmin 1.e-6
#define GCmax 0.999999
#define titvmin 1.e-6
#define titvmax 9999.
#define frmin 1.e-6
#define frmax 0.999999
#define gammamin 0.005
#define gammamax 100.
#define covmin 1.e-6
#define covmax 100.
#define pimin 1.e-6
#define pimax 0.999999

#define GAMMA_RANGE_ITER 5
#define GAMMA_RANGE_SIZE 10
#define ERROR_THRESHOLD -1.e-9



typedef double proba[4][4];
typedef double dproba[5][4][4]; /* param type, nt1, nt2 */
typedef double like[4];

		/********* STRUCTURES *********/


typedef struct noeud{
	struct noeud *v1, *v2, *v3;	/* neighbors */
	struct noeud *remain;		/* unchanged neighbor */
	double l1, l2, l3;		/* branch lengths */
	double lremain;			/* unchanged length */
	int alremain;			/* unchanged status */
	double b1, b2, b3;		/* bootstrap values */
	int alive1, alive2, alive3;
	int numerobr;			
	double aeq, ceq, geq, teq;	/* eq. base comp. in parent branch */
	int nbpattern;
	int *patternson1, *patternson2;
	int *underlyingbr, *underlyinggc;
	int depth;			
	like **x;			/* likelihood (rate class, pattern) */
	like **dx;			/* lkh 1st derivative (rate class, pattern) */
	like **d2x;			/* lkh 2nd derivative (rate class, pattern) */
/* dx and d2x are for the parameter currently derived */
	proba** p1;			/* subst. proba matrix toward v1 */
	proba** p2;			/* subst. proba matrix toward v2 */
	dproba** dp1;			/* p1 first derivative */
	dproba** dp2;			/* p2 first derivative */
	dproba** d2p1;			/* p1 second derivative */
	dproba** d2p2;			/* p2 second derivative */
/* these probas have two additional dimensions with 
respct to nhml1: initial rate class and final rate class */
	double pdif1, pdif2;		/* proba switch rate class i to rate class j */
        double pide1, pide2;		/* proba no switch */
        double** rr1, **rr2;		/* rij */
	double** drr1, **drr2;		/* rij derivatives with respect */
        double** d2rr1, **d2rr2;		/* to x=l.covar */
	double ***x2sites;		/* lkh of the current two sites (rate class, state1, state2) */
	char nom[MAXLNAME+1];		/* name */
	char* seq;			/* data if tip */ 	
	double a, c, g, t;		/* base composition */	
} *noeud;


typedef struct branche{
	noeud bout1;
	noeud bout2;
	double length;
	double bootstrap;
} *branche;


typedef struct tree{
	noeud* node;			/* array of nodes */
	int nbseq;			/* number of compared sequences */
	int lgseq;			/* length of compared sequences (patterns) */
	char** seq;			/* sequences */
	char** names;			/* sequence names */
        int* sitetopatt; 		/* site to "pattern" correspondence */
	double** dist;			/* distances between sequences */
	double* weight;			/* weights at each site */
	double weightsum;		/* sum of weights */
	double titv;			/* transitions / transversions ratio */
	double froot;			/* root->l1 / root->l2 ratio */
	double GCanc;			/* ancestral G+C content */
	double gamma_shp;
	double covar;
	double pi;
	double omega1, omega2;
	branche* listbr;		/* list of branches */
	int* alivebr;			/* status of branches */
	int nbalivebr;
	double* lgbr_init;		/* list of initial branch lengths */
	double* listgc;			/* list of G+C contents */
	int nb_var_param;		/* number of optimized parameters */
	int which_var_param[7];
	int* param_nature;		/* 0:ti/tv, 1:fr, 2:gcanc, 3:length, 4:comp 5:cov 6:pi */		
	double** var_param;		/* array of addresses of optimized parameters */
	double** classlike_s;		/* class-dependant likelihood at each site */
	double* like_s;			/* likelihood at each site (site) */
        double* likeasrv_s;
	double*** classdlike_s;		/* class-dependant derivative at each site */
	double** dlike_s;		/* derivatives of likelihood at each site (param, site) */
	double** dlikeasrv_s;
	double*** classd2like_s;	/* class-dependant second derivative at each site */
	double** d2like_s;		/* second deriv. of like. at each site (param, site) */
	double** d2likeasrv_s;
	double lkh;			/* total log likelihood ln(L) */
        double lssrv, lasrv;
	double* deriv;			/* ln(L) derivatives (param) */
	double* deriv2;			/* ln(L) second and cross derivatives (param) */
} tree;


typedef struct upgmanoeud{
	struct upgmanoeud *v1, *v2, *v3;	/* neighbors */
	double l1, l2, l3;
        double b1, b2, b3;
        int num;
        int alive;
        double depth;
	char* nom;
} *upgmanoeud;


typedef struct print_option{
	int PRINT0;		/* more */
 	int PRINT1;		/* and  */
	int PRINT2;		/* more */
	int PRINT3;		/* output */
	int COV_SITE_PATTERN;	/* details about site-specific covariation */
	int OUTPUT_COMP;
	int EVAL_OUT;		/* detailed results in file detailed_out */
} *print_option;


typedef struct compute_option{
	char MODEL[10];
	int OPTIMIZE_LENGTH;
	int OPTIMIZE_COMP;
	int OPTIMIZE_TITV;
	int OPTIMIZE_ROOT;
	int OPTIMIZE_ANC;
	int OPTIMIZE_GAMMA;
	int OPTIMIZE_COV;
	int OPTIMIZE_PI;
} *compute_option;	


typedef struct init_option{
	char INIT_LENGTH[10];
	char INIT_COMP[10];
	double INIT_TITV;
	double INIT_ROOT;
	double INIT_ANC;
	double INIT_GAMMA;
	double INIT_COV;
        double INIT_PI;
	int GAMMA_NCL;
	int NBRANDOM;
	int NOISE;
}*init_option;


typedef struct converge_option{
	int PRECISION;
	int MAXITERATION;
	int SIMPLEX;
}* converge_option;


typedef struct options{
	init_option init;
	compute_option compute;
	converge_option converge;
	print_option print;
	int ALLCOUPLES;			/* for program EVAL  */
	int SH_G;			/* for program SHAKE */
	double SH_MAXLCROSSED;		/* for program SHAKE */
	int SH_MAXBOOTCROSSED;          /* for program SHAKE */
	int SH_RESTART;			/* for program SHAKE */
} *options;








