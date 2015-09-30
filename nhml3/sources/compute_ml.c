#include "nhmlg.h"

#define PointGamma(prob,alpha,beta) PointChi2(prob,2.0*(alpha))/(2.0*(beta))

#define pos(i,j,n)      ((i)*(n)+(j))

time_t tdeb, tfin, ttotdeb, ttotfin;
double probatime=0., dliketime=0., d2liketime=0., tottime;

FILE* probaoutcov;

extern tree* s_tree;
extern compute_option compute_o;
extern noeud root;

upgmanoeud upgmaroot;

double param[4*MAXNSP];
double min_value[7]={titvmin, frmin, GCmin, lmin, GCmin, covmin, pimin};
double max_value[7]={titvmax, frmax, GCmax, lmax, GCmax, covmax, pimax};

#define GAMMA_RANGE_SIZE 9

double gamma_range[GAMMA_RANGE_SIZE]={0.1, 0.2, 0.4, 0.6, 0.8, 1., 2., 5., 10.};

int paramcase[33][17];

double *gamma_clmean; 
int gamma_nbcl;
double covar;
double pi;
int exact_cov;
double **R, **P, **T92;
int noprint;

proba **p1prov, **p2prov;
dproba **dp1prov, **dp2prov; /* utilities */
dproba **d2p1prov, **d2p2prov;


double *vector(long nl, long nh);
void free_vector(double *v, long nl, long nh);
void amoeba_cont(double **p, double y[], int ndim, double ftol,
	double (*funk)(double []), int *nfunk, double** contraintes);
void like_node(noeud nd);

FILE* dbg, *dbgprov;


void printtime(double rtime){

  int d, h, m, s;

  d=(int)rtime/86400;
  rtime-=d*86400;
  if(d) printf("%d day", d);
  if(d>1) printf("s");
  if(d) printf(" ");

  h=(int)rtime/3600;
  rtime-=h*3600;
  if(h) printf("%d hour", h);
  if(h>1) printf("s");
  if(h) printf(" ");

  m=(int)rtime/60;
  rtime-=m*60;
  if(m) printf("%d minute", m);
  if(m>1) printf("s");
  if(m) printf(" ");

  s=(int)rtime;
  if(s) printf("%d second", s);
  if(s>1) printf("s");

  if(d==0 && h==0 && m==0 && s==0)
    printf("less than 1 second");

  printf("\n");
}



int eigen(int job, double A[], int n, double rr[], double ri[], double vr[], double vi[], double w[]);



/* init_cas */
/* Initialise le tableau paramcase qui sert au calcul des derivees secondes */
/* et croisees de la vraisemblance. */
/* premier indice (1->32): type de cas */
/* Un cas est un type de "relations" entre deux parametres vis-a-vis */
/* du processus de derivation. Il en existe 32, selon la nature des 2 */
/* parametres (branch length, Ts/Tv ration, GC equilibre, ...) et leur */
/* position dans l'arbre par rapport au noeud ou l'on calcule la derivee. */
/* deuxieme indice (1->16): rang des 16 termes de l'expression de la derivee */
/*			    seconde d'un produit de 4 facteurs: */
/* 										*/
/* d2(p1.L1.p2.L2)/dxdy =   d2p1/dxdy . L1 . p2 . L2   (terme 1)  */
/*			  + dp1/dx . dL1/dy . p2 . L2  (terme 2)  */
/*			  + dp1/dx . L1 . dp2/dy . L2  (terme 3)  */
/*			  + ...					  */
/*				...				  */
/*			  + p1 . L1 . p2 . d2L2/dxdy   (terme 16) */
/* 										*/
/* param_case explicite les termes de la somme qui valent 0 (valeur a 0) et ceux */
/* qu'il faut calculer (valeur a 1), pour chaque cas. */


void init_cas(){

  int i, j;

  for(i=0;i<33;i++)
    for(j=0;j<17;j++) 
      paramcase[i][j]=0;

  for(i=1;i<=16;i++) paramcase[i][i]=1;  
  for(i=1;i<=16;i++) paramcase[17][i]=1;
  for(i=1;i<=8;i++) paramcase[18][2*i-1]=1;
  for(i=1;i<=4;i++) paramcase[19][4*i-3]=1;
  for(i=1;i<=4;i++) paramcase[20][4*i-1]=1;
  for(i=1;i<=4;i++) paramcase[21][4*i-2]=1;
  for(i=1;i<=4;i++) paramcase[22][4*i]=1;
  paramcase[23][1]=1; paramcase[23][3]=1; paramcase[23][9]=1; paramcase[23][11]=1;
  paramcase[24][2]=1; paramcase[24][10]=1;
  paramcase[25][4]=1; paramcase[25][12]=1;
  paramcase[26][1]=1; paramcase[26][9]=1;
  paramcase[27][3]=1; paramcase[27][11]=1;
  paramcase[28][5]=1; paramcase[28][7]=1;
  paramcase[29][13]=1; paramcase[29][15]=1;
  paramcase[30][1]=1; paramcase[30][3]=1;
  paramcase[31][9]=1; paramcase[31][11]=1;

}



/* which_cas */
/* Determine dans quel cas (vis-a-vis du processus de derivation) */
/* se trouve une paire de parametres, decrits par leur nature */
/* np1 et np2 (0=Ts/Tv ratio, 1=root fraction, 2=GCanc, 3=branch length, */
/* 4=GC eq) et par leur position dans l'arbre par rapport au noeud */
/* dont on calcule la derivee, und1 et und2 (0 = au dessus, 1=en dessous */
/* a gauche mais pas branche fille, 2=en dessous a droite mais pas branche */
/* fille, 10=branche fille gauche, 20=branche fille droite, 12=branche */
/* racine). Si le noeud ou l'on calcule la derivee est la racine, rootnode */
/* est passe a 1 (cas particulier). */



int which_cas(int np1, int np2, int und1, int und2, int rootnode){

  if(np1==2 || np2==2) { printf("GCanc exception\n"); return -1; }
  if(np2<np1) {printf("Bad param order\n"); return -1; }

  if((np1==1 || np2==1) && !rootnode) return 32; 
  if(((np1==3 && und1==12) || (np2==3 && und2==12)) && !rootnode) return 32;
  if((np1==3 || np1==4) && und1==0) return 32;
  if((np2==3 || np2==4) && und2==0) return 32;
  
  if((np1==3 || np1==4) && (np2==3 || np2==4) && und1==10 && und2==10) return 1;
  if((np1==3 || np1==4) && (np2==3 || np2==4) && und1==10 && und2==1 ) return 2;
  if((np1==3 || np1==4) && (np2==3 || np2==4) && und1==10 && und2==20) return 3;
  if((np1==3 || np1==4) && (np2==3 || np2==4) && und1==10 && und2==2 ) return 4;
  if((np1==3 || np1==4) && (np2==3 || np2==4) && und1==1  && und2==10) return 5;
  if((np1==3 || np1==4) && (np2==3 || np2==4) && und1==1  && und2==1 ) return 6;
  if((np1==3 || np1==4) && (np2==3 || np2==4) && und1==1  && und2==20) return 7;
  if((np1==3 || np1==4) && (np2==3 || np2==4) && und1==1  && und2==2 ) return 8;
  if((np1==3 || np1==4) && (np2==3 || np2==4) && und1==20 && und2==10) return 9;
  if((np1==3 || np1==4) && (np2==3 || np2==4) && und1==20 && und2==1 ) return 10;
  if((np1==3 || np1==4) && (np2==3 || np2==4) && und1==20 && und2==20) return 11;
  if((np1==3 || np1==4) && (np2==3 || np2==4) && und1==20 && und2==2 ) return 12;
  if((np1==3 || np1==4) && (np2==3 || np2==4) && und1==2  && und2==10) return 13;
  if((np1==3 || np1==4) && (np2==3 || np2==4) && und1==2  && und2==1 ) return 14;
  if((np1==3 || np1==4) && (np2==3 || np2==4) && und1==2  && und2==20) return 15;
  if((np1==3 || np1==4) && (np2==3 || np2==4) && und1==2  && und2==2 ) return 16;

  if(np1==0 && np2==0) return 17;
  if(np1==5 && np2==5) return 17;

  if(np1==0 && (np2==1 || (np2==3 && und2==12))) return 18;
  
  if(np1==0 && (np2==3 || np2==4) && und2==10) return 19;
  if(np1==0 && (np2==3 || np2==4) && und2==20) return 20;
  if(np1==0 && (np2==3 || np2==4) && und2==1 ) return 21;
  if(np1==0 && (np2==3 || np2==4) && und2==2 ) return 22;
  if((np1==1 || (np1==3 && und1==12)) && (np2==1 || (np2==3 && und2==12))) return 23;
  if(np1==1 && (np2==3 || np2==4) && und2==1) return 24;
  if((np1==3 && und1==12) && ((np2==3 || np2==4) && und2==1)) return 24;
  if(np1==1 && (np2==3 || np2==4) && und2==2) return 25;
  if((np1==3 && und1==12) && ((np2==3 || np2==4) && und2==2)) return 25;
  if(np1==1 && np2==4 && und2==10) return 26;
  if((np1==3 && und1==12) && ((np2==3 || np2==4) && und2==10)) return 26;
  if(np1==1 && np2==4 && und2==20) return 27;
  if((np1==3 && und1==12) && ((np2==3 || np2==4) && und2==20)) return 27;
  if(np1==3 && np2==3 && und1==1 && und2==12) return 28;
  if(np1==3 && np2==3 && und1==2 && und2==12) return 29;
  if(np1==3 && np2==3 && und1==10 && und2==12) return 30;
  if(np1==3 && np2==3 && und1==20 && und2==12) return 31;

  
  printf("unexpected case : %d %d %d %d \n", np1, np2, und1, und2);
exit(EXIT_FAILURE);
  return 32;
  
}


void bit1(int *plist, int num){
  num--;
  plist+=(num/lmot);
  *plist |= (1<<(num%lmot));
}

void bit0(int *plist, int num){
  num--;
  plist+=(num/lmot);
  *plist &=  ~(1<<(num%lmot));
}


void et(int *listet, int *list1, int *list2, int len){
  int i;
  for(i=0; i< len; i++)
	listet[i]= list1[i] & list2[i];
}


void ou(int *listou, int *list1, int *list2, int len){
  int i;
  for(i=0; i< len; i++)
    listou[i]= list1[i] | list2[i];
}


void non(int *listnon, int *list, int len){
  int i;
  for(i=0; i< len; i++)
	listnon[i]= ~list[i];
}


int testbit(int *plist, int num){
  num--;
  plist += (num/lmot);
  return (*plist) & (1<<(num%lmot));
}



/* mini */

int mini(int a, int b){
if (a<b) return(a); else return(b);
}

/* maxi */

int maxi(int a, int b){
if (a>b) return(a); else return(b);
}


/* check_alloc */
/* Return a pointer with needed allocated memory or leave program */
/* if not enough memory. */

void *check_alloc(int nbrelt, int sizelt)
{
void *retval;
if( (retval=calloc(nbrelt,sizelt)) != NULL ) {
  return retval; 
}
printf("Not enough memory\n");
exit(EXIT_FAILURE);
}

void printmemory_init(){
  char *line, *ret;
  int maxlline=1000;


  dbg=fopen("debugNHML", "w");
  system("ps -aux > debugNHML.prov");
  dbgprov=fopen("debugNHML.prov", "r");
  if(dbgprov==NULL){printf("Warning: problem when debugging\n"); return;}
  line=check_alloc(maxlline+1, sizeof(char));
  ret=fgets(line, maxlline, dbgprov);
  if(ret!=NULL) fprintf(dbg, "             %s", line);

  fclose(dbgprov);
  fclose(dbg);
  free(line);
}


void printmemory(char title_arg[]){

  char* line, *ret, *syscall, *title, *prov;
  int maxlline=1000;
  int thispid, lgtit=10, lgarg, i;

  thispid=(int)getpid();
  lgarg=(int)strlen(title_arg);
  title=check_alloc(lgtit+3, sizeof(char));
  sprintf(title, "%.10s", title_arg);
  prov=title;
  while(*prov) prov++;
  for(i=0;i<lgtit-lgarg;i++){*prov=' '; prov++;}


  syscall=check_alloc(500, sizeof(char));
  sprintf(syscall, "ps -aux | grep eval_nhg | grep %d > debugNHML.prov", thispid);

  system(syscall);

  dbgprov=fopen("debugNHML.prov", "r");
  dbg=fopen("debugNHML", "a");
  line=check_alloc(maxlline+1, sizeof(char));

  while(1){
    ret=fgets(line, maxlline, dbgprov);
    if(ret==NULL) break;
    fprintf(dbg, "%s: %s", title, line);
  }

  fclose(dbgprov);
  fclose(dbg);
  system("rm -f debugNHML.prov");
  free(line);
  free(syscall);
  free(title);
}






/* ludcmp */
/* from Numerical Recipes in C */
/* Replace matrix a by a rowwise permutation of its LU decomposition. */

int ludcmp(double **a, int n, int *indx, double *d){
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;

	vv=vector(1,n);
	*d=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(a[i][j])) > big) big=temp;
		if (big == 0.0) return 0;
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=a[i][j];
			for (k=1;k<i;k++) sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=a[i][j];
			for (k=1;k<j;k++)
				sum -= a[i][k]*a[k][j];
			a[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=a[imax][k];
				a[imax][k]=a[j][k];
				a[j][k]=dum;
			}
			*d = -(*d);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (a[j][j] == 0.0) a[j][j]=1.0e-20;
		if (j != n) {
			dum=1.0/(a[j][j]);
			for (i=j+1;i<=n;i++) a[i][j] *= dum;
		}
	}
	free_vector(vv,1,n);
	return 1;
}



/* lubksb */
/* from Numerical Recipes in C */

void lubksb(double **a, int n, int *indx, double b[]){
	int i,ii=0,ip,j;
	double sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= a[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= a[i][j]*b[j];
		b[i]=sum/a[i][i];
	}
}



/* matvect */
/* Return the product vector of matrix m1(nbl1, nbc1) by column vector v2. */

int matvect(double** m1, int nbl1, int nbc1, double* v2, int nbl2, double* vprod){

  int i, k;

  if(nbc1!=nbl2) return 0;  

  for(i=0;i<nbl1;i++){
    vprod[i]=0;
    for(k=0;k<nbc1;k++) vprod[i]+=m1[i][k]*v2[k];
  }
  return 1;
}


/* matmat */
/* Product of matrices. */
/* Matrices m1(nbl1, nbc1) and m2(nbl2, nbc2) are multiplicated. */
/* Product is written into prod(nbl1, nbc2). */
/* prod allocation is not performed. */

int matmat(double** m1, int nbl1, int nbc1, double** m2, int nbl2, int nbc2, double** prod){

  int i, j, k;

  if(nbc1!=nbl2)  return 0;
  
  for(i=0;i<nbl1;i++){
    for(j=0;j<nbc2;j++){
      prod[i][j]=0.;
      for(k=0;k<nbc1;k++)
	prod[i][j]+=m1[i][k]*m2[k][j];
    }
  }
  return 1;
}



/* vectmat */
/* Return the product vector of row vector v1 by matrix m2(nbl1, nbc1). */

double* vectmat(double* v1, int nbc1, double** m2, int nbl2, int nbc2){
  
  int j, k;
  double* vprod;

  if(nbc1!=nbl2) return 0;
  vprod=(double*)check_alloc(nbc2+1, sizeof(double));
  for(j=0;j<nbc2;j++){
    vprod[j]=0;
    for(k=0;k<nbc1;k++) vprod[j]+=v1[k]*m2[k][j];
  }
  return vprod;
}





/* invmat */
/* Invert square matrix mat(n,n). Result in invmat.*/
/* If mat is singular, 0 is returned (1 otherwise) */

int invmat(double** mat, int n, double** invmat){

  int i, j, *indx;
  double **lu, d, *col;

  lu=check_alloc(n+1, sizeof(double*));
  for(i=1;i<n+1;i++){
    lu[i]=(double*)check_alloc(n+1, sizeof(double));
    for(j=1;j<n+1;j++) lu[i][j]=mat[i-1][j-1];
  }
  indx=(int*)check_alloc(n+1, sizeof(int));
  col=(double*)check_alloc(n+1, sizeof(double));

  if(!ludcmp(lu, n, indx, &d)) return 0;
  
  for(j=1;j<=n;j++){
    for(i=1;i<=n;i++) col[i]=0.;
    col[j]=1.;
    lubksb(lu, n, indx, col);
    for(i=1;i<=n;i++) invmat[i-1][j-1]=col[i];
  }
  for(i=1; i<n+1;i++) free(lu[i]);
  free(lu);
  free(indx); free(col);
  return 1;
}



/* functions concerning the CDF and percentage points of the gamma and
   Chi2 distribution
*/


double PointNormal (double prob)
{
/* returns z so that Prob{x<z}=prob where x ~ N(0,1) and (1e-12)<prob<1-(1e-12)
   returns (-9999) if in error
   Odeh RE & Evans JO (1974) The percentage points of the normal distribution.
   Applied Statistics 22: 96-97 (AS70)

   Newer methods:
     Wichura MJ (1988) Algorithm AS 241: the percentage points of the
       normal distribution.  37: 477-484.
     Beasley JD & Springer SG  (1977).  Algorithm AS 111: the percentage 
       points of the normal distribution.  26: 118-121.

*/
   double a0=-.322232431088, a1=-1, a2=-.342242088547, a3=-.0204231210245;
   double a4=-.453642210148e-4, b0=.0993484626060, b1=.588581570495;
   double b2=.531103462366, b3=.103537752850, b4=.0038560700634;
   double y, z=0, p=prob, p1;

   p1 = (p<0.5 ? p : 1-p);
   if (p1<1e-20) return (-9999);

   y = sqrt (log(1/(p1*p1)));   
   z = y + ((((y*a4+a3)*y+a2)*y+a1)*y+a0) / ((((y*b4+b3)*y+b2)*y+b1)*y+b0);
   return (p<0.5 ? -z : z);
}


double LnGamma (double alpha)
{
/* returns ln(gamma(alpha)) for alpha>0, accurate to 10 decimal places.  
   Stirling's formula is used for the central polynomial part of the procedure.
   Pike MC & Hill ID (1966) Algorithm 291: Logarithm of the gamma function.
   Communications of the Association for Computing Machinery, 9:684
*/
   double x=alpha, f=0, z;

   if (x<7) {
      f=1;  z=x-1;
      while (++z<7)  f*=z;
      x=z;   f=-log(f);
   }
   z = 1/(x*x);
   return  f + (x-0.5)*log(x) - x + .918938533204673 
	  + (((-.000595238095238*z+.000793650793651)*z-.002777777777778)*z
	       +.083333333333333)/x;  
}



double IncompleteGamma (double x, double alpha, double ln_gamma_alpha)
{
/* returns the incomplete gamma ratio I(x,alpha) where x is the upper 
	   limit of the integration and alpha is the shape parameter.
   returns (-1) if in error
   ln_gamma_alpha = ln(Gamma(alpha)), is almost redundant.
   (1) series expansion     if (alpha>x || x<=1)
   (2) continued fraction   otherwise
   RATNEST FORTRAN by
   Bhattacharjee GP (1970) The incomplete gamma integral.  Applied Statistics,
   19: 285-287 (AS32)
*/
   int i;
   double p=alpha, g=ln_gamma_alpha;
   double accurate=1e-8, overflow=1e30;
   double factor, gin=0, rn=0, a=0,b=0,an=0,dif=0, term=0, *pn;

   if (x==0) return (0);
   if (x<0 || p<=0) return (-1);

  pn=check_alloc(6, sizeof(double));

   factor=exp(p*log(x)-x-g);   
   if (x>1 && x>=p) goto l30;
   /* (1) series expansion */
   gin=1;  term=1;  rn=p;
 l20:
   rn++;
   term*=x/rn;   gin+=term;

   if (term > accurate) goto l20;
   gin*=factor/p;
   goto l50;
 l30:
   /* (2) continued fraction */
   a=1-p;   b=a+x+1;  term=0;
   pn[0]=1;  pn[1]=x;  pn[2]=x+1;  pn[3]=x*b;
   gin=pn[2]/pn[3];
 l32:
   a++;  b+=2;  term++;   an=a*term;
   for (i=0; i<2; i++) pn[i+4]=b*pn[i+2]-an*pn[i];
   if (pn[5] == 0) goto l35;
   rn=pn[4]/pn[5];   dif=fabs(gin-rn);
   if (dif>accurate) goto l34;
   if (dif<=accurate*rn) goto l42;
 l34:
   gin=rn;
 l35:
   for (i=0; i<4; i++) pn[i]=pn[i+2];
   if (fabs(pn[4]) < overflow) goto l32;
   for (i=0; i<4; i++) pn[i]/=overflow;
   goto l32;
 l42:
   gin=1-factor*gin;

 l50:
   free(pn);
   return (gin);
}



double PointChi2 (double prob, double v)
{
/* returns z so that Prob{x<z}=prob where x is Chi2 distributed with df=v
   returns -1 if in error.   0.000002<prob<0.999998
   RATNEST FORTRAN by
       Best DJ & Roberts DE (1975) The percentage points of the 
       Chi2 distribution.  Applied Statistics 24: 385-388.  (AS91)
   Converted into C by Ziheng Yang, Oct. 1993.
*/
   double e=.5e-6, aa=.6931471805, p=prob, g;
   double xx, c, ch, a=0,q=0,p1=0,p2=0,t=0,x=0,b=0,s1,s2,s3,s4,s5,s6;

   if (p<.000002 || p>.999998 || v<=0) return (-1);

   g = LnGamma (v/2);
   xx=v/2;   c=xx-1;
   if (v >= -1.24*log(p)) goto l1;

   ch=pow((p*xx*exp(g+xx*aa)), 1/xx);
   if (ch-e<0) return (ch);
   goto l4;
l1:
   if (v>.32) goto l3;
   ch=0.4;   a=log(1-p);
l2:
   q=ch;  p1=1+ch*(4.67+ch);  p2=ch*(6.73+ch*(6.66+ch));
   t=-0.5+(4.67+2*ch)/p1 - (6.73+ch*(13.32+3*ch))/p2;
   ch-=(1-exp(a+g+.5*ch+c*aa)*p2/p1)/t;
   if (fabs(q/ch-1)-.01 <= 0) goto l4;
   else                       goto l2;
  
l3: 
   x=PointNormal (p);
   p1=0.222222/v;   ch=v*pow((x*sqrt(p1)+1-p1), 3.0);
   if (ch>2.2*v+6)  ch=-2*(log(1-p)-c*log(.5*ch)+g);
l4:
   q=ch;   p1=.5*ch;
   if ((t=IncompleteGamma (p1, xx, g))<0) {
      printf ("\nerr IncompleteGamma");
      return (-1);
   }
   p2=p-t;
   t=p2*exp(xx*aa+g+p1-c*log(ch));   
   b=t/ch;  a=0.5*t-b*c;

   s1=(210+a*(140+a*(105+a*(84+a*(70+60*a))))) / 420;
   s2=(420+a*(735+a*(966+a*(1141+1278*a))))/2520;
   s3=(210+a*(462+a*(707+932*a)))/2520;
   s4=(252+a*(672+1182*a)+c*(294+a*(889+1740*a)))/5040;
   s5=(84+264*a+c*(175+606*a))/2520;
   s6=(120+c*(346+127*c))/5040;
   ch+=t*(1+0.5*t*s1-b*c*(s1-b*(s2-b*(s3-b*(s4-b*(s5-b*s6))))));
   if (fabs(q/ch-1) > e) goto l4;

   return (ch);
}



int DiscreteGamma (double freqK[], double rK[], 
    double alfa, double beta, int K, int median)
{
/* discretization of gamma distribution with equal proportions in each 
   category
*/
   int i;
   double gap05=1.0/(2.0*K), t, factor=alfa/beta*K, lnga1;

   if (median) {
      for (i=0; i<K; i++) rK[i]=PointGamma((i*2.0+1)*gap05, alfa, beta);
      for (i=0,t=0; i<K; i++) t+=rK[i];
      for (i=0; i<K; i++)     rK[i]*=factor/t;
   }
   else {
      lnga1=LnGamma(alfa+1);
      for (i=0; i<K-1; i++)
	 freqK[i]=PointGamma((i+1.0)/K, alfa, beta);
      for (i=0; i<K-1; i++)
	 freqK[i]=IncompleteGamma(freqK[i]*beta, alfa+1, lnga1);
      rK[0] = freqK[0]*factor;
      rK[K-1] = (1-freqK[K-2])*factor;
      for (i=1; i<K-1; i++)  rK[i] = (freqK[i]-freqK[i-1])*factor;
   }
   for (i=0; i<K; i++) freqK[i]=1.0/K;

   return (0);
}



/* samename */
/* TRUE if names name1 and name2 are equal, FALSE otherwise */

int samename(char* name1, char* name2){
 
  int i=0;
  char c1, c2;
  
  while(i<20){
    c1=name1[i]; c2=name2[i];
    if((c1==' ' || c1=='\n' || c1=='\t' || c1==0) && (c2==' ' || c2=='\n' || c2=='\t' || c2==0)) return TRUE;
    if(c1!=c2) return FALSE;
    i++;
  }
  return TRUE;
}





/* refresh */
/* Remove all sites with gaps only from sequence alignment seq. */
/* If option==0, all sites with at least one gap are removed. */

void refresh(char** seq, int nbseq, int option){

  int lgseq, l=-1, drapeau, i, j, k;
  char **seqref ;

  lgseq=(int)strlen(seq[0]);
  seqref=(char**)check_alloc(nbseq, sizeof(char*));
  for(i=0;i<nbseq;i++)
     seqref[i]=(char*)check_alloc(lgseq+1, sizeof(char));

  if (option==0){
    for(i=0;i<lgseq;i++){
      drapeau=0;
      for(j=0;j<nbseq;j++){
	if (seq[j][i]!='A' && seq[j][i]!='C' && seq[j][i]!='G' && seq[j][i]!='T' ){
          drapeau=1; 
          break;
        }
      }
      if (drapeau==0){
	l++;
	for(k=0;k<nbseq;k++){
    	  *(seqref[k]+l)=*(seq[k]+i);
	}
      }	
    }
  }
  else{
    for(i=0;i<lgseq;i++){
      drapeau=0;
      for(j=0;j<nbseq;j++){
	if (*(seq[j]+i)=='A' || *(seq[j]+i)=='C' || *(seq[j]+i)=='T' || *(seq[j]+i)=='G') { 
	  drapeau=1;
	  break;
	}
      }
      if (drapeau==1){
	l++;
	for(k=0;k<nbseq;k++)
	  if (seq[k][i]=='A' || seq[k][i]=='C' || seq[k][i]=='G' || seq[k][i]=='T')
	    *(seqref[k]+l)=*(seq[k]+i);
	  else
	    *(seqref[k]+l)=' ';
      }		
    }
  }
  for(i=0;i<nbseq;i++)
    for (j=l+1;j<lgseq;j++) 
      *(seqref[i]+j)='\0';
  for (i=0;i<nbseq;i++) 
    for (j=0;j<lgseq;j++)
      *(seq[i]+j)=*(seqref[i]+j);
  for(i=0;i<nbseq;i++) free(seqref[i]);
  free(seqref);
}


int samesite(char* site, tree* s_tree, int i){

  int j;

  for(j=0;j<s_tree->nbseq;j++)
    if(s_tree->node[j]->seq[i]!=site[j]) return 0;
  
  return 1;
}




void err_mess(char* option){
  printf("Missing or wrong line specifying option %s in option file\nBy-default option used\n", option);
}




void getoptions(options* opt_arg, FILE* in_opt){

  char* line, *text, *prov;
  int i, a, ret;
  options opt;

/* by-default options */
  
  opt=(struct options*)check_alloc(1, sizeof(struct options));
  opt->init=(struct init_option*)check_alloc(1, sizeof(struct init_option));
  opt->compute=(struct compute_option*)check_alloc(1, sizeof(struct compute_option));
  opt->converge=(struct converge_option*)check_alloc(1, sizeof(struct converge_option));
  opt->print=(struct print_option*)check_alloc(1, sizeof(struct print_option));

  sprintf(opt->compute->MODEL, "T3");
  opt->compute->OPTIMIZE_LENGTH=TRUE;
  opt->compute->OPTIMIZE_COMP=TRUE;
  opt->compute->OPTIMIZE_TITV=TRUE;
  opt->compute->OPTIMIZE_ROOT=TRUE;
  opt->compute->OPTIMIZE_ANC=TRUE;
  opt->compute->OPTIMIZE_GAMMA=TRUE;
  opt->compute->OPTIMIZE_COV=FALSE;
  opt->compute->OPTIMIZE_PI=FALSE;
  sprintf(opt->init->INIT_LENGTH, "REDO");
  sprintf(opt->init->INIT_COMP, "VAR");
  opt->init->INIT_TITV=-1.;
  opt->init->INIT_ANC=-1.;
  opt->init->INIT_ROOT=-1.;
  opt->init->NBRANDOM=0;
  opt->init->NOISE=0;
  opt->init->INIT_GAMMA=-1.;
  opt->init->GAMMA_NCL=4;
  opt->init->INIT_COV=-1.;
  opt->init->INIT_PI=-1.;
  opt->converge->PRECISION=4;
  opt->converge->SIMPLEX=TRUE;
  opt->print->PRINT1=TRUE;
  opt->print->PRINT2=FALSE;
  opt->print->PRINT3=FALSE;
  opt->print->EVAL_OUT=FALSE;
  opt->print->OUTPUT_COMP=FALSE;
  opt->print->COV_SITE_PATTERN=FALSE;
  opt->SH_G=-1;
  opt->SH_MAXLCROSSED=lmax;
  opt->SH_MAXBOOTCROSSED=INT_MAX;
  opt->SH_RESTART=0;
  opt->ALLCOUPLES=0;

  
  line=(char*)check_alloc(100, sizeof(char));
  text=(char*)check_alloc(5000, sizeof(char));

  while(fgets(line, 100, in_opt))
    strcat(text, line);





  prov=strstr(text, "ALLCOUPLES");
  if(prov)
    prov=strchr(prov, '=');
  else { err_mess("ALLCOUPLES"); goto finac; }
  if(prov)
    prov++;
  else { err_mess("ALLCOUPLES"); goto finac; }
  while(*prov==' ') prov++;
  sscanf(prov, "%d", &(opt->ALLCOUPLES));
  if(opt->ALLCOUPLES!=0 && opt->ALLCOUPLES!=1) {
    err_mess("ALLCOUPLES"); opt->ALLCOUPLES=1;
  }
  finac:


  prov=strstr(text, "INIT_GC");
  if(prov)
    prov=strchr(prov, '=');
  else { err_mess("INIT_GC"); goto fininitcomp; }
  if(prov)
    prov++;
  else { err_mess("INIT_GC"); goto fininitcomp; }
  while(*prov==' ') prov++;
  i=0;
  while(*prov!=' ' && *prov!='\t' && *prov!='\n'){
    opt->init->INIT_COMP[i]=*prov;
    i++; prov++;
  }
  opt->init->INIT_COMP[i]=0;
  if(strcmp(opt->init->INIT_COMP, "CONST") && strcmp(opt->init->INIT_COMP, "VAR") && strcmp(opt->init->INIT_COMP, "BALANCED")){
    err_mess("INIT_GC"); sprintf(opt->init->INIT_COMP, "VAR");
  }
  fininitcomp:




  prov=strstr(text, "INIT_LENGTH");
  if(prov)
    prov=strchr(prov, '=');
  else { err_mess("INIT_LENGTH"); goto fininitlg; }
  if(prov)
    prov++;
  else { err_mess("INIT_LENGTH"); goto fininitlg; }
  while(*prov==' ') prov++;
  i=0;
  while(*prov!=' ' && *prov!='\t' && *prov!='\n'){
    opt->init->INIT_LENGTH[i]=*prov;
    i++; prov++;
  }
  opt->init->INIT_LENGTH[i]=0;
  if(strcmp(opt->init->INIT_LENGTH, "KEEP") && strcmp(opt->init->INIT_LENGTH, "REDO")){
    err_mess("INIT_LENGTH"); sprintf(opt->init->INIT_LENGTH, "REDO");
  }
  fininitlg:



  prov=strstr(text, "OPTIMIZE_LENGTH");
  if(prov)
    prov=strchr(prov, '=');
  else { err_mess("OPTIMIZE_LENGTH"); goto finol; }
  if(prov)
    prov++;
  else { err_mess("OPTIMIZE_LENGTH"); goto finol; }
  while(*prov==' ') prov++;
  sscanf(prov, "%d", &(opt->compute->OPTIMIZE_LENGTH));
  if(opt->compute->OPTIMIZE_LENGTH!=0 && opt->compute->OPTIMIZE_LENGTH!=1) {
    err_mess("OPTIMIZE_LENGTH"); opt->compute->OPTIMIZE_LENGTH=1;
  }
  finol:

 
  prov=strstr(text, "OPTIMIZE_GC");
  if(prov)
    prov=strchr(prov, '=');
  else { err_mess("OPTIMIZE_GC"); goto finoc; }
  if(prov)
    prov++;
  else { err_mess("OPTIMIZE_GC"); goto finoc; }
  while(*prov==' ') prov++;
  sscanf(prov, "%d", &(opt->compute->OPTIMIZE_COMP));
  if(opt->compute->OPTIMIZE_COMP!=0 && opt->compute->OPTIMIZE_COMP!=1) {
    err_mess("OPTIMIZE_GC"); opt->compute->OPTIMIZE_COMP=1;
  }
  finoc:

  prov=strstr(text, "OPTIMIZE_TITV");
  if(prov)
    prov=strchr(prov, '=');
  else { err_mess("OPTIMIZE_TITV"); goto finoo; }
  if(prov)
    prov++;
  else { err_mess("OPTIMIZE_TITV"); goto finoo; }
  while(*prov==' ') prov++;
  sscanf(prov, "%d", &(opt->compute->OPTIMIZE_TITV));
  if(opt->compute->OPTIMIZE_TITV!=0 && opt->compute->OPTIMIZE_TITV!=1) {
    err_mess("OPTIMIZE_TITV"); opt->compute->OPTIMIZE_TITV=1;
  }
  finoo:

  prov=strstr(text, "OPTIMIZE_ROOT");
  if(prov)
    prov=strchr(prov, '=');
  else { err_mess("OPTIMIZE_ROOT"); goto finor; }
  if(prov)
    prov++;
  else { err_mess("OPTIMIZE_ROOT"); goto finor; }
  while(*prov==' ') prov++;
  sscanf(prov, "%d", &(opt->compute->OPTIMIZE_ROOT));
  if(opt->compute->OPTIMIZE_ROOT!=0 && opt->compute->OPTIMIZE_ROOT!=1) {
    err_mess("OPTIMIZE_ROOT"); opt->compute->OPTIMIZE_ROOT=1;
  }
  finor:

  prov=strstr(text, "OPTIMIZE_ANC");
  if(prov)
    prov=strchr(prov, '=');
  else { err_mess("OPTIMIZE_ANC"); goto finoa; }
  if(prov)
    prov++;
  else { err_mess("OPTIMIZE_ANC"); goto finoa; }
  while(*prov==' ') prov++;
  sscanf(prov, "%d", &(opt->compute->OPTIMIZE_ANC));
  if(opt->compute->OPTIMIZE_ANC!=0 && opt->compute->OPTIMIZE_ANC!=1) {
    err_mess("OPTIMIZE_ANC"); opt->compute->OPTIMIZE_ANC=1;
  }
  finoa:


  prov=strstr(text, "OPTIMIZE_GAMMA");
  if(prov)
    prov=strchr(prov, '=');
  else { err_mess("OPTIMIZE_GAMMA"); goto finog; }
  if(prov)
    prov++;
  else { err_mess("OPTIMIZE_GAMMA"); goto finog; }
  while(*prov==' ') prov++;
  sscanf(prov, "%d", &(opt->compute->OPTIMIZE_GAMMA));
  if(opt->compute->OPTIMIZE_GAMMA!=0 && opt->compute->OPTIMIZE_GAMMA!=1) {
    err_mess("OPTIMIZE_GAMMA"); opt->compute->OPTIMIZE_GAMMA=1;
  }
  finog:

  prov=strstr(text, "OPTIMIZE_COV");
  if(prov)
    prov=strchr(prov, '=');
  else { err_mess("OPTIMIZE_COV"); goto finoco; }
  if(prov)
    prov++;
  else { err_mess("OPTIMIZE_COV"); goto finoco; }
  while(*prov==' ') prov++;
  sscanf(prov, "%d", &(opt->compute->OPTIMIZE_COV));
  if(opt->compute->OPTIMIZE_COV!=0 && opt->compute->OPTIMIZE_COV!=1) {
    err_mess("OPTIMIZE_COV"); opt->compute->OPTIMIZE_COV=0;
  }
  finoco:

  prov=strstr(text, "OPTIMIZE_PI");
  if(prov)
    prov=strchr(prov, '=');
  else { err_mess("OPTIMIZE_PI"); goto finopi; }
  if(prov)
    prov++;
  else { err_mess("OPTIMIZE_PI"); goto finopi; }
  while(*prov==' ') prov++;
  sscanf(prov, "%d", &(opt->compute->OPTIMIZE_PI));
  if(opt->compute->OPTIMIZE_PI!=0 && opt->compute->OPTIMIZE_PI!=1) {
    err_mess("OPTIMIZE_PI"); opt->compute->OPTIMIZE_PI=0;
  }
  finopi:


  prov=strstr(text, "INIT_ROOT");
  if(prov)
    prov=strchr(prov, '=');
  else { err_mess("INIT_ROOT"); goto finkrl; }
  if(prov)
    prov++;
  else { err_mess("INIT_ROOT"); goto finkrl; }
  while(*prov==' ') prov++;
  a=sscanf(prov, "%le", &(opt->init->INIT_ROOT));
  if( opt->init->INIT_ROOT>1. || a<1){
    err_mess("INIT_ROOT"); opt->init->INIT_ROOT=-1.;
  }
  finkrl:



  prov=strstr(text, "INIT_TITV");
  if(prov)
    prov=strchr(prov, '=');
  else { err_mess("INIT_TITV"); goto finitt; }
  if(prov)
    prov++;
  else { err_mess("INIT_TITV"); goto finitt; }
  while(*prov==' ') prov++;
  a=sscanf(prov, "%le", &(opt->init->INIT_TITV));
  if(a<1){
    err_mess("INIT_TITV"); opt->init->INIT_TITV=-1.;
  }
  finitt:


  prov=strstr(text, "INIT_ANC");
  if(prov)
    prov=strchr(prov, '=');
  else { err_mess("INIT_ANC"); goto finian; }
  if(prov)
    prov++;
  else { err_mess("INIT_ANC"); goto finian; }
  while(*prov==' ') prov++;
  a=sscanf(prov, "%le", &(opt->init->INIT_ANC));
  if(a<1){
    err_mess("INIT_ANC"); opt->init->INIT_ANC=-1.;
  }
  finian:

  prov=strstr(text, "INIT_COV");
  if(prov)
    prov=strchr(prov, '=');
  else { err_mess("INIT_COV"); goto finico; }
  if(prov)
    prov++;
  else { err_mess("INIT_COV"); goto finico; }
  while(*prov==' ') prov++;
  a=sscanf(prov, "%le", &(opt->init->INIT_COV));
  if(a<1){
    err_mess("INIT_COV"); opt->init->INIT_COV=-1.;
  }
  finico:

  prov=strstr(text, "INIT_PI");
  if(prov)
    prov=strchr(prov, '=');
  else { err_mess("INIT_PI"); goto finipi; }
  if(prov)
    prov++;
  else { err_mess("INIT_PI"); goto finipi; }
  while(*prov==' ') prov++;
  a=sscanf(prov, "%le", &(opt->init->INIT_PI));
  if(a<1){
    err_mess("INIT_PI"); opt->init->INIT_PI=-1.;
  }
  finipi:




  prov=strstr(text, "PRECISION");
  if(prov)
    prov=strchr(prov, '=');
  else { err_mess("PRECISION"); goto finpr; }
  if(prov)
    prov++;
  else { err_mess("PRECISION"); goto finpr; }
  while(*prov==' ') prov++;
  sscanf(prov, "%d", &(opt->converge->PRECISION));
  if(opt->converge->PRECISION<=-5) {
    err_mess("PRECISION"); opt->converge->PRECISION=4;
  }
  finpr:


  prov=strstr(text, "SIMPLEX");
  if(prov)
    prov=strchr(prov, '=');
  else { err_mess("SIMPLEX"); goto finsp; }
  if(prov)
    prov++;
  else { err_mess("SIMPLEX"); goto finsp; }
  while(*prov==' ') prov++;
  sscanf(prov, "%d", &(opt->converge->SIMPLEX));
  if(opt->converge->SIMPLEX!=0 && opt->converge->SIMPLEX!=1) {
    err_mess("SIMPLEX"); opt->converge->SIMPLEX=1;
  }
  finsp:


  prov=strstr(text, "PRINT1");
  if(prov)
    prov=strchr(prov, '=');
  else { err_mess("PRINT1"); goto finp1; }
  if(prov)
    prov++;
  else { err_mess("PRINT1"); goto finp1; }
  while(*prov==' ') prov++;
  sscanf(prov, "%d", &(opt->print->PRINT1));
  if(opt->print->PRINT1!=0 && opt->print->PRINT1!=1) {
    err_mess("PRINT1"); opt->print->PRINT1=1;
  }
  finp1:

  prov=strstr(text, "PRINT2");
  if(prov)
    prov=strchr(prov, '=');
  else { err_mess("PRINT2"); goto finp2; }
  if(prov)
    prov++;
  else { err_mess("PRINT2"); goto finp2; }
  while(*prov==' ') prov++;
  sscanf(prov, "%d", &(opt->print->PRINT2));
  if(opt->print->PRINT2!=0 && opt->print->PRINT2!=1) {
    err_mess("PRINT2"); opt->print->PRINT2=1;
  }

  finp2:

  prov=strstr(text, "EVAL_OUT");
  if(prov)
    prov=strchr(prov, '=');
  else { err_mess("EVAL_OUT"); goto finevo; }
  if(prov)
    prov++;
  else { err_mess("EVAL_OUT"); goto finevo; }
  while(*prov==' ') prov++;
  sscanf(prov, "%d", &(opt->print->EVAL_OUT));
  if(opt->print->EVAL_OUT!=0 && opt->print->EVAL_OUT!=1) {
    err_mess("EVAL_OUT"); opt->print->EVAL_OUT=1;
  }

  finevo:

/*
  prov=strstr(text, "COV_SITE_PATTERN");
  if(prov)
    prov=strchr(prov, '=');
  else { err_mess("COV_SITE_PATTERN"); goto fincsp; }
  if(prov)
    prov++;
  else { err_mess("COV_SITE_PATTERN"); goto fincsp; }
  while(*prov==' ') prov++;
  sscanf(prov, "%d", &(opt->print->COV_SITE_PATTERN));
  if(opt->print->COV_SITE_PATTERN!=0 && opt->print->COV_SITE_PATTERN!=1) {
    err_mess("COV_SITE_PATTERN"); opt->print->COV_SITE_PATTERN=1;
  }

  fincsp:
*/



  prov=strstr(text, "SH_G");
  if(prov)
    prov=strchr(prov, '=');
  else { err_mess("SH_G"); goto finshg; }
  if(prov)
    prov++;
  else { err_mess("SH_G"); goto finshg; }
  while(*prov==' ') prov++;
  sscanf(prov, "%d", &(opt->SH_G));

  finshg:

  prov=strstr(text, "SH_MAXLCROSSED");
  if(prov)
    prov=strchr(prov, '=');
  else { err_mess("SH_MAXLCROSSED"); goto finshm; }
  if(prov)
    prov++;
  else { err_mess("SH_MAXLCROSSED"); goto finshm; }
  while(*prov==' ') prov++;
  sscanf(prov, "%le", &(opt->SH_MAXLCROSSED));
  if(opt->SH_MAXLCROSSED<=0) opt->SH_MAXLCROSSED=lmax;

  finshm:


  prov=strstr(text, "SH_MAXBOOTCROSSED");
  if(prov)
    prov=strchr(prov, '=');
  else { err_mess("SH_MAXBOOTCROSSED"); goto finshb; }
  if(prov)
    prov++;
  else { err_mess("SH_MAXBOOTCROSSED"); goto finshb; }
  while(*prov==' ') prov++;
  sscanf(prov, "%d", &(opt->SH_MAXBOOTCROSSED));

  if(opt->SH_MAXBOOTCROSSED<=0) opt->SH_MAXBOOTCROSSED=INT_MAX;

  finshb:


  prov=strstr(text, "SH_RESTART");
  if(prov)
    prov=strchr(prov, '=');
  else { err_mess("SH_RESTART"); goto finshr; }
  if(prov)
    prov++;
  else { err_mess("SH_RESTART"); goto finshr; }
  while(*prov==' ') prov++;
  sscanf(prov, "%d", &(opt->SH_RESTART));

  if(opt->SH_RESTART<=0) opt->SH_RESTART=0;

  finshr:


  prov=strstr(text, "INIT_GAMMA");
  if(prov)
    prov=strchr(prov, '=');
  else { err_mess("INIT_GAMMA"); goto fingam; }
  if(prov)
    prov++;
  else { err_mess("INIT_GAMMA"); goto fingam; }
  while(*prov==' ') prov++;
  ret=sscanf(prov, "%le", &(opt->init->INIT_GAMMA));
  if(ret==0) { err_mess("INIT_GAMMA"); opt->init->INIT_GAMMA=-1.;}

  fingam:


  prov=strstr(text, "GAMMA_NBCLASS");
  if(prov)
    prov=strchr(prov, '=');
  else { err_mess("GAMMA_NBCLASS"); goto fingac; }
  if(prov)
    prov++;
  else { err_mess("GAMMA_NBCLASS"); goto fingac; }
  while(*prov==' ') prov++;
  ret=sscanf(prov, "%d", &(opt->init->GAMMA_NCL));
  if(ret==0 || opt->init->GAMMA_NCL<2){
    err_mess("GAMMA_NBCLASS"); 
    opt->init->GAMMA_NCL=4;
  }

  fingac:

  if(opt->init->INIT_GAMMA<0. && opt->compute->OPTIMIZE_GAMMA==0)
    opt->init->GAMMA_NCL=1;


  *opt_arg=opt;

  free(line);
  free(text);

  return;

}




char* readtree(FILE* treef, int nbsp){

  char *line;
  char *ctree, *ctreefin, *prov;
  int i;

  line=(char*)check_alloc(MAXLLINE+1, sizeof(char));
  ctree=check_alloc(3*MAXLNAME*nbsp, sizeof(char));
  ctreefin=ctree;
  
  line[0]='\n';
  while(line[0]=='\n') fgets(line, MAXLLINE, treef);
  i=0;
  while(line[i]==' ' || line[i]=='\t') i++;
  if(line[i]!='(' && line[i]!='['){
    printf("Tree first character should be '[' or '('\n");
    return NULL;
  }
  if(line[i]=='['){
    while((prov=strchr(line, ']'))==NULL){
      if(!fgets(line, MAXLLINE, treef)){
	printf("Unterminated comment in tree ('[' missing)\n");
        return NULL;
      }
    }
    prov++;
    while(*prov==' ' || *prov=='\t') prov++;
    if(*prov=='\n'){
      line[0]='\n';
      while(line[0]=='\n' && fgets(line, MAXLLINE, treef));
      prov=line;
      while(*prov==' ' || *prov=='\t') prov++;
    }
    if(*prov!='('){
      printf("Comments should be followed by '('\n");
      return NULL;
    }
  }
  else prov=line+i;	/* prov = premiere parenthese */


  sprintf(ctree, "%s", prov);
  while(*ctreefin!='\n' && *ctreefin!=';' && *ctreefin) ctreefin++;
  if(*ctreefin==';') return ctree;
  while(fgets(line, MAXLLINE, treef)){
    sprintf(ctreefin, "%s", line);
    if(prov=strchr(line, ';')){
      *(prov+1)=0;
      return ctree;
    }
    while(*ctreefin!='\n' && *ctreefin) ctreefin++;  
  }
  *ctreefin=0;

  free(line);
  return ctree;
}

/* upgma_ctreewrite */

char* upgma_ctreewrite(upgmanoeud comesfrom, upgmanoeud node, char* ctree, int nodegc){
  upgmanoeud vd, vg;
  double l, b;

  if(comesfrom==node->v1) { vg=node->v2; vd=node->v3; l=node->l1; b=node->b1;}
  else if(comesfrom==node->v2) { vg=node->v1; vd=node->v3; l=node->l2; b=node->b2;}
  else { vg=node->v1; vd=node->v2; l=node->l3; b=node->b3;}


  if (vg==NULL){
    sprintf(ctree,"%s", node->nom);
    while(*ctree) ctree++;
    if (l>-0.5) sprintf(ctree,":%.4f",l);
    while(*ctree) ctree++;
    return ctree;
  }
  else{
    *(ctree++)='(';
    ctree=upgma_ctreewrite(node, vg, ctree, nodegc);
  }
  sprintf(ctree,", ");
  ctree+=2;
  ctree=upgma_ctreewrite(node, vd, ctree, nodegc);
  *(ctree++)=')';
  if (b>-0.1) sprintf(ctree,"%.2f",b);
  while(*ctree) ctree++;
  if (l>-0.9) sprintf(ctree,":%.6f", l);
  while(*ctree) ctree++;
  return ctree;
}



/* upgma_stoc */

void upgma_stoc(upgmanoeud* arbre_s, int racine, int notu, char* ctree, int nodegc){
  char* tring, *prov;

  tring=(char*)check_alloc((50+MAXLNAME)*notu,sizeof(char));

  if (racine) {
    tring=upgma_ctreewrite(NULL, upgmaroot, ctree, nodegc);
    prov=ctree+strlen(ctree)-1;
    while(*prov!=')') *(prov--)='\0';
    strcat(ctree,";");
  }
  else{
    tring=upgma_ctreewrite(arbre_s[2*notu-3]->v3, arbre_s[2*notu-3], ctree, nodegc);
    prov=ctree+strlen(ctree)-1;
    while(*prov!=')') *(prov--)='\0';
    *(prov++)=',';
    *(prov++)=' ';
    tring=upgma_ctreewrite(arbre_s[2*notu-3], arbre_s[2*notu-3]->v3, prov, nodegc);
    strcat(ctree,");");
  }
}




/* ctreewrite */
/* Write string describing the s_tree underlying to noeud node at address ctree. */

char* ctreewrite(noeud comesfrom, noeud node, char* ctree, int nodegc){
  noeud vd, vg;
  double l, b;

  if(comesfrom==node->v1) { vg=node->v2; vd=node->v3; l=node->l1;}
  else if(comesfrom==node->v2) { vg=node->v1; vd=node->v3; l=node->l2;}
  else { vg=node->v1; vd=node->v2; l=node->l3;}
  if(nodegc==1) b=node->c+node->g;
  else b=node->ceq+node->geq;
  if(nodegc<0) b=-1.;
  b*=100.;

  if (vg==NULL){
    sprintf(ctree,"%s", node->nom);
    while(*ctree) ctree++;
    if (b>-0.1) sprintf(ctree," %.2f",b);
    while(*ctree) ctree++;
    if (l>-0.5) sprintf(ctree,":%.4f",l);
    while(*ctree) ctree++;
    return ctree;
  }
  else{
    *(ctree++)='(';
    ctree=ctreewrite(node, vg, ctree, nodegc);
  }
  sprintf(ctree,", ");
  ctree+=2;
  ctree=ctreewrite(node, vd, ctree, nodegc);
  *(ctree++)=')';
  if (b>-0.1) sprintf(ctree,"%.2f",b);
  while(*ctree) ctree++;
  if (l>-0.9) sprintf(ctree,":%.6f", l);
  while(*ctree) ctree++;
  return ctree;
}



/* stoc */
/* StructTOChaine. Write string describing s_tree arbre_s at address ctree */
/* racine = unrooted/rooted (0/other) , notu = number of leaves. */

void stoc(noeud* arbre_s, int racine, int notu, char* ctree, int nodegc){
  char* tring, *prov;

  prov=ctree;
  while(*prov) {*prov=0; prov++;}
  tring=(char*)check_alloc(50*notu,sizeof(char));

  if (racine) {
    tring=ctreewrite(NULL, root, ctree, nodegc);
    prov=ctree+strlen(ctree)-1;
    while(*prov!=')') *(prov--)='\0';
    strcat(ctree,";");
  }
  else{
    tring=ctreewrite(arbre_s[2*notu-3]->v3, arbre_s[2*notu-3], ctree, nodegc);
    prov=ctree+strlen(ctree)-1;
    while(*prov!=')') *(prov--)='\0';
    *(prov++)=',';
    *(prov++)=' ';
    tring=ctreewrite(arbre_s[2*notu-3], arbre_s[2*notu-3]->v3, prov, nodegc);
    strcat(ctree,");");
  }
}


int retder(int *liste){
  int i=0, j;
  while (liste[i] != 0) i++;
  j = *(liste + i - 1);
  *(liste + i - 1) = 0;
  return j;
}

void aj(int *liste, int nb) {
  int  i=0;
  while (liste[i] != 0) i++;
  *(liste + i) = nb;
  return;
}


/* ctob */
/* ChaineTOBits. Read c_tree carbre and write b_tree barbre. */

int ctob(char *carbre, char *nom[], int *barbre[]){


  int   i=0, j, k, fin=0, nbpo=0,
	nbpf = 0, cptv1 = 0, nbotu, br_ouverte, *listecour, 
	otu = -1, pomoinspf, cpttree=0, t=1;
  char  c, cc, dejalu = '\0', readname[30];


  sscanf(carbre+(cpttree++), "%c", &c);
  if (c == '[') {
    while ((c != ']') && c) sscanf(carbre+(cpttree++), "%c", &c);
    if (c != ']'){
      printf("Unmatched '[' ']'\n");
      return -1;
    }
  } 
  else
    if (c == '(') cpttree=0;
    else{
      printf("Tree file 1st character must be '(' or '['\n");
      return -1;
    }
  while ((c != ';') && c ) {
    sscanf(carbre+(cpttree++), "%c", &c);
    if (c == '(') nbpo++;
    if (c == ')') nbpf++;
    if ((nbpo == nbpf + 1) && (c == ',')) cptv1++;
  }

  if (c != ';'){
    printf("';' missing at end of tree\n");
    return -1;
  }

  if (nbpo != nbpf){
    printf("Unmatched parenthesis\n");
    return -1;
  }

  if (cptv1 == 2) /* unrooted : ok */;
  else 
    if (cptv1 == 1){ /* rooted : problem */ 
      printf("Unexpected rooted tree.\n"); 
      return -1; 
    }
    else{ /* bad tree string */
      printf("Bad number of ',' in tree\n");
      return -1;
    }

  nbotu=nbpo+2;


  listecour=(int*) check_alloc(nbotu, sizeof(int));

  cpttree=0;
  sscanf(carbre+(cpttree++), "%c", &c);
  if (c == '['){
    while (c != ']') sscanf(carbre+(cpttree++), "%c", &c);
    while((c==']') || (c==' ') || (c=='\n') || (c=='\t')) sscanf(carbre+(cpttree++), "%c", &c);
    if (c!='(') return -1;
  }
  else
    while(c!='(') sscanf(carbre+(cpttree++), "%c", &c);

  pomoinspf=1;
  for (i = 0; i < nbotu; i++) listecour[i] = 0;
  br_ouverte = 0;


  for (k = 0; t==1; k++) {
    if (dejalu == '\0') sscanf(carbre+(cpttree++), "%c", &c);
    else {
      c = dejalu;
      dejalu = '\0';
    }
    switch (c) {
    case ';': fin = 1; break;
    case ',': case '\n': case '\t': case '\'': case ' ': break;
    case '(':
      pomoinspf ++;
      br_ouverte++;
      aj(listecour, br_ouverte); break;
    case ')':
      pomoinspf--;
      if (*(carbre+cpttree+1) == ';' || pomoinspf==0) {
	fin = 1; break;
      }
      j = retder(listecour);
      while(carbre[cpttree]!=',' && carbre[cpttree]!=')') cpttree++;
      break;
    default:
      otu=-1; cc = c; i = 0;
      while ((cc != ':') && (cc != ',') && (cc != ')') && (cc != '\n') && (cc != ' ')) {
	if (cc != '\'') { readname[i++] = cc; readname[i]='\0'; }
	sscanf(carbre+(cpttree++), "%c", &cc);
      }
      for(j=0;j<nbotu;j++) if(samename(readname, nom[j])) { otu=j; break; }
      if (otu==-1) {printf("Unknown name : %s\n", readname); return -1; }
      for(i=0;i<nbotu-3;i++) bit0(barbre[i], otu+1);
      for (i = 0; i < nbotu-3; i++)
	if (listecour[i] != 0)
	  bit1(barbre[listecour[i]-1], otu+1);
      cpttree--;
      while(carbre[cpttree]!=',' && carbre[cpttree]!=')') cpttree++; 
    }
    if (fin == 1) break;
  }

  free(listecour);
  return nbotu; 
}



/* ctot */
/* ChaineTOTable. Read c_tree input (string) and write t_tree arbre (int**), */
/* branch length lgbi (internal) and lgbp (terminal), bootstrap values, */
/* species names (nom) and rooted/unrooted (racine-> r (rooted) or n (not). */

int ctot(char *input, int *arbre[], double *lgbi, double *lgbp, double* bootstrap, char *nom[], char *racine, int* nbbi){


  int   i=0, j, k, fdf=0, fin=0, nbpo=0,
	nbpf = 0, cptv1 = 0, nbotu, br_ouverte, *listecour,
	otu = -1, pomoinspf, cpttree=0, t=1, nbnom=0;
  char  c, cc, cas, dejalu = '\0';
  double  f;

  sscanf(input+(cpttree++), "%c", &c);
  if (c == '[') {
    while ((c != ']') && (fdf != EOF)) fdf = sscanf(input+(cpttree++), "%c", &c);
    if (c != ']'){
      printf("Unmatched '[' ']'\n");
      return -1;
    }
  } 
  else
    if (c == '(') cpttree=0;
    else{
      printf("Tree file 1st character must be '(' or '['\n");
      return -1;
    }
  while ((c != ';') && (fdf != EOF)) {
    fdf = sscanf(input+(cpttree++), "%c", &c);
      if (c == '(') nbpo++;
      if (c == ')') nbpf++;
      if ((nbpo == nbpf + 1) && (c == ',')) cptv1++;
  }

  if (c != ';'){
    printf("';' missing at end of tree file\n");
    return -1;
  }

  if (nbpo != nbpf){
    printf("Unmatched parenthesis\n");
    return -1;
  }

  if (cptv1 == 1) cas = 'c';
  if (cptv1 == 2) cas = 'a';

  if ((cptv1!=1) && (cptv1!=2)){
    printf("Bad number of ',' in tree file\n");
    return -1;
  }

  nbotu = nbpo + 2;


  listecour=(int*) check_alloc(nbotu, sizeof(int));

  if (cas == 'a') *racine='n'; else *racine='r';

  if (cas=='c'){
    if (lgbp) lgbp[0] = 0.;
    strcpy(nom[0],"ROOT");
    for(i=0;i<nbotu-3;i++) arbre[0][i]=0;
    otu++;
  }


  cpttree=0;
  sscanf(input+(cpttree++), "%c", &c);
  if (c == '['){
    while (c != ']') sscanf(input+(cpttree++), "%c", &c);
    while((c==']') || (c==' ') || (c=='\n') || (c=='\t')) sscanf(input+(cpttree++), "%c", &c);
    if (c!='(') return -1;
  }
  else
    while(c!='(') sscanf(input+(cpttree++), "%c", &c);

  pomoinspf=1;
  for (i = 0; i < nbotu; i++) listecour[i] = 0;
  br_ouverte = 0;


  for (k = 0; t==1; k++) {
    if (dejalu == '\0') sscanf(input+(cpttree++), "%c", &c);
    else {
      c = dejalu;
      dejalu = '\0';
    }
    switch (c) {
    case ';': fin = 1; break;
    case ',': case '\n': case '\t': case '\'': case ' ': break;
    case '(':
      pomoinspf ++;
      br_ouverte++;
      aj(listecour, br_ouverte); break;
    case ')':
      pomoinspf--;
      sscanf(input+(cpttree++), "%c", &cc);
      if (cc == ';' || pomoinspf==0) {
	fin = 1; break;
      }
      j = retder(listecour);
      while (cc=='\n' || cc==' ' || cc=='\t') sscanf(input+(cpttree++),"%c",&cc);
      if (strpbrk(input+cpttree-1, "-0123456789.Ee")==input+cpttree-1){
        sscanf(input+cpttree-1, "%le", bootstrap+j-1);
        cpttree+=strspn(input+cpttree-1,"-.0123456789Ee");
        cc=*(input+cpttree-1);
        while (cc=='\n' || cc==' ' || cc=='\t') sscanf(input+(cpttree++),"%c",&cc);
      }
      if (cc == ':') {
	while(input[cpttree]==' ') cpttree++;
	sscanf(input+cpttree, "%le", &f);
	cpttree+=strspn(input+cpttree,"-0123456789.Ee");
	lgbi[j - 1] = f;
      }
      else dejalu=cc;
      break;
    default:
      nbnom++;
      otu++; cc = c; i = 0;
      while ((cc != ':') && (cc != ',') && (cc != ')') && (cc != '\n') && (cc!=' ')) {
	if (cc != '\'') { nom[otu][i] = cc; i++; nom[otu][i]='\0'; }
	sscanf(input+cpttree, "%c", &cc);
        cpttree++;
      }
      while ((cc != ':') && (cc != ',') && (cc != ')') && (cc != '\n')){
	sscanf(input+cpttree, "%c", &cc);
        cpttree++;
      }
      if (cc == ':') {
	while(input[cpttree]==' ') cpttree++;
	sscanf(input+(cpttree), "%le", &f);
	cpttree+=strspn(input+cpttree,"-0123456789.Ee");
	lgbp[otu] = f;
      } 
      else dejalu = cc;
      for (i = 0; i < nbotu - 3; i++) arbre[otu][i] = 0;
      for (i = 0; i < nbotu; i++)
	if (listecour[i] != 0)
	  arbre[otu][listecour[i] - 1] = 1;
    }
    if (fin == 1) break;
  }

  if(nbbi) *nbbi=br_ouverte;

  free(listecour);

  return nbnom;

/*  if (cas=='a') return nbotu; else return (nbotu-1);*/
}



/* create_node */
/* Return a node whose "parent", "childs" and values are */
/* v1, v2, v3, l1, l2, l3, nom, order. */

noeud create_node(noeud v1, noeud v2, noeud v3, double l1, double l2, double l3, double b1, double b2, double b3, char* nom) {

  noeud nd;

  nd = (struct noeud*) check_alloc(1, sizeof(struct noeud));
  nd->v1 = v1; nd->v2 = v2; nd->v3 = v3; 
  nd->l1 = l1; nd->l2 = l2; nd->l3 = l3;
  nd->alive1 = (int)b1; nd->alive2 = (int)b2; nd->alive3 = (int)b3;

  if (nom!=NULL) {strncpy(nd->nom, nom, MAXLNAME); nd->nom[MAXLNAME]='\0'; }
  return nd;
}


/* bottomnode */
/* Return the root node of rooted s_tree including node nd. */
/*    !!BEWARE!!    */
/* If the s_tree including node nd is unrooted, program will bug */
/* (infinite loop). */

noeud bottomnode(noeud nd){
  if(nd->v3 == NULL) return nd;
  return (bottomnode(nd->v3));
}


/* ttos */
/* TableTOStruct. Create s_tree arbre_s from : t_tree arbre_t, leaves number notu, */
/* branch lengths lgbi (internal branches) , lgbp (terminal branches). */
/* If tree has no branch length, arguments 3 and 4 must be NULL. */


int ttos(int** arbre_t, int notu, double* lgbi, double* lgbp, double* alive, char** nom, noeud* arbre_s){

  noeud   	nd, p1, p2, p3;
  int           notuv, i, j, k, sommebi, tax1, tax2, n1, n2, n3, cpt1=-1;
  int           *kill_tax, *kill_bi;
  double	arg6, arg9;

  kill_tax = (int *)check_alloc(notu, sizeof(int));
  kill_bi = (int *)check_alloc(notu, sizeof(int));

  notuv = notu;
  for (i = 0; i <notu; i++){
    kill_tax[i] = 0;
    kill_bi[i] = 0;
  }

	/* terminal nodes */

  for (i = 0; i < notu; i++){
    if(lgbp==NULL) arg6=-1.0; else arg6=lgbp[i];
    arbre_s[i] = create_node(NULL, NULL, NULL, 0., 0., arg6, -1., -1., 1., nom[i]);
  }


	/* internal nodes */

  for (i = 0; i < notu - 3; i++) {
	/* determination de la bi a creer */
    for (j = 0; j < notu - 3 ; j++) {
      if (kill_bi[j] == 0) {
	sommebi = 0;
	for (k = 0; k < notu; k++)
	  if (kill_tax[k] == 0) {
	    sommebi += arbre_t[k][j];
	  }
	if (sommebi == 2 || notuv - sommebi == 2)
	  break;
      }
    }

	/* determination des 2 otus/noeuds fils */
    if (sommebi == 2) {
      for (k = 0; k < notu; k++)
	if (arbre_t[k][j] == 1 && kill_tax[k] == 0) {
	  tax1 = k;
	  break;
	}
      for (k = tax1 + 1; k < notu; k++)
	if (arbre_t[k][j] == 1 && kill_tax[k] == 0) {
	  tax2 = k;
	  break;
	}
    } else {
      for (k = 0; k < notu; k++)
	if (arbre_t[k][j] == 0 && kill_tax[k] == 0) {
	  tax1 = k;
	  break;
	}
      for (k = tax1 + 1; k < notu; k++)
	if (arbre_t[k][j] == 0 && kill_tax[k] == 0) {
	  tax2 = k;
	  break;
	}
    }

    p1 = bottomnode(arbre_s[tax1]);
    p2 = bottomnode(arbre_s[tax2]);

    if (lgbp) arg6=lgbi[j]; else arg6=0.;
    if (alive) arg9=alive[j]; else arg9=-1.;

    arbre_s[notu + (++cpt1)] = create_node(p1, p2, NULL, p1->l3, p2->l3, arg6, (double)p1->alive3, (double)p2->alive3, arg9, NULL);


    p1->v3 = arbre_s[notu + cpt1];
    p2->v3 = arbre_s[notu + cpt1];

    kill_tax[tax1] = kill_bi[j] = 1;
    notuv--;
  }


	/* last node */


  for (i = 0; i < 2 * notu - 3; i++) 
    if (arbre_s[i]->v3 == NULL) {
      n1 = i;
      break;
    }

  for (i = n1 + 1; i < 2 * notu - 3; i++) 
    if (arbre_s[i]->v3 == NULL) {
      n2 = i;
      break;
    }
  
  for (i = n2 + 1; i < 2 * notu - 3; i++) 
    if (arbre_s[i]->v3 == NULL) {
      n3 = i;
      break;
    }
 

  arbre_s[2 * notu - 3] = create_node(arbre_s[n1], arbre_s[n2], arbre_s[n3], arbre_s[n1]->l3, arbre_s[n2]->l3, arbre_s[n3]->l3, arbre_s[n1]->alive3, arbre_s[n2]->alive3, arbre_s[n3]->alive3, NULL);
  arbre_s[n1]->v3 = arbre_s[2 * notu - 3];
  arbre_s[n2]->v3 = arbre_s[2 * notu - 3];
  arbre_s[n3]->v3 = arbre_s[2 * notu - 3];


  free(kill_tax);
  free(kill_bi);
  return 0;
}





/* mindepth */
/* Set node's depth in a rooted tree. */
/* mindepth(terminal node) = 0 */
/* mindepth(internal node) = number of branches between node and closest terminal node */
/* node comesfrom (one of node's neighbors) is root. */

int mindepth (noeud comesfrom, noeud node){
  noeud droit, gauche;

  if (node->v1==NULL) { node->depth=0; return 0; }
  if (node->v2==NULL) { node->depth=0; return 0; }
  
  if (comesfrom==node->v1) { 
    gauche=node->v2; 
    droit=node->v3; 
  }
  else if (comesfrom==node->v2) { 
    gauche=node->v1; 
    droit=node->v3; 
  }
  else { 
    gauche=node->v1; 
    droit=node->v2; 
  }

  node->depth=mini(mindepth(node, gauche), mindepth(node, droit))+1;
  return node->depth;
}



void makelistbr_unrooted(noeud from, noeud nd, branche* br, int* alivebr,int nbbranch){
  noeud to1, to2;
  double lto1, lto2;
  int i=0, alive1, alive2;



  if (nd==root){
    if(from==root->v1) makelistbr_unrooted(root, nd->v2, br, alivebr, nbbranch);
    if(from==root->v2) makelistbr_unrooted(root, nd->v1, br, alivebr, nbbranch);
    return;
  }


  if (from==nd->v1) {to1=nd->v2; to2=nd->v3; lto1=nd->l2; lto2=nd->l3; alive1=nd->alive2; alive2=nd->alive3;}
  else if (from==nd->v2) {to1=nd->v1; to2=nd->v3; lto1=nd->l1; lto2=nd->l3; alive1=nd->alive1; alive2=nd->alive3;}
  else if (from==nd->v3) {to1=nd->v1; to2=nd->v2; lto1=nd->l1; lto2=nd->l2; alive1=nd->alive1; alive2=nd->alive2;}
  else { 
    printf("erreur making branch list.\n");
  }


  while(i<nbbranch && br[i]->bout1) i++;

  if (to1 && to1==root) {
    if(nd==root->v1){
      lto1+=root->l2;
      br[i]->bout1=nd; br[i]->bout2=root->v2; br[i]->length=lto1; i++;
    }
    else{
      lto1+=root->l1;
      br[i]->bout1=nd; br[i]->bout2=root->v1; br[i]->length=lto1; i++;
    }
    if(alivebr) alivebr[i-1]=1;
  }


  if (to2 && to2==root) {
    if(nd==root->v1){
      lto2+=root->l2;
      br[i]->bout1=nd; br[i]->bout2=root->v2; br[i]->length=lto2; i++;
    }
    else{
      lto2+=root->l1;
      br[i]->bout1=nd; br[i]->bout2=root->v1; br[i]->length=lto2; i++;
    }
    if(alivebr) alivebr[i-1]=1;
  }

  

  if (to1 && to1!=root) { 
    br[i]->bout1=nd; 
    br[i]->bout2=to1; 
    br[i]->length=lto1; 
    if(alivebr) alivebr[i]=alive1; 
    i++;
  }

  if (to2 && to2!=root) { 
    br[i]->bout1=nd; 
    br[i]->bout2=to2; 
    br[i]->length=lto2; 
    if(alivebr) alivebr[i]=alive2;
  }

  if (to1) makelistbr_unrooted(nd, to1, br, alivebr, nbbranch); 
  if (to2) makelistbr_unrooted(nd, to2, br, alivebr, nbbranch);
}




void organize_tree(noeud from, noeud nd){
  noeud prov;
  double lprov;
  int alprov;
  static int numappel=0;

  if(!nd->v1) return;
  if(nd!=root){
    numappel++;
    sprintf(nd->nom, "int%d", numappel);
  }
  if(from==nd->v1) {
    prov=nd->v3; nd->v3=nd->v1; nd->v1=prov;
    lprov=nd->l3; nd->l3=nd->l1; nd->l1=lprov;
    alprov=nd->alive3; nd->alive3=nd->alive1; nd->alive1=alprov;
  }
  if(from==nd->v2) {
    prov=nd->v3; nd->v3=nd->v2; nd->v2=prov;
    lprov=nd->l3; nd->l3=nd->l2; nd->l2=lprov;
    alprov=nd->alive3; nd->alive3=nd->alive2; nd->alive2=alprov;
  }
  organize_tree(nd, nd->v1);
  organize_tree(nd, nd->v2);

  if(nd==root) numappel=0;
}



void invert_branch(branche br){
  noeud prov;
      
  prov=br->bout1; br->bout1=br->bout2; br->bout2=prov;

}



void organize_listbr_node(noeud nd, branche* listbr, int nbseq){
 
  int i, break1=0, break2=0;

  if(nd->v1==NULL) return;

  for(i=0;i<nbseq;i++){
    if(!break1 && (listbr[i]->bout1==nd->v1 && listbr[i]->bout2==nd)){
      invert_branch(listbr[i]);
      break1=1;
    }
    if(!break1 && (listbr[i]->bout1==nd && listbr[i]->bout2==nd->v1)){
      break1=1;
    }
    if(!break2 && (listbr[i]->bout1==nd->v2 && listbr[i]->bout2==nd)){
      invert_branch(listbr[i]);
      break2=1;
    }
    if(!break2 && (listbr[i]->bout1==nd && listbr[i]->bout2==nd->v2)){
      break2=1;
    }
    if(break1 && break2) break;
  }

  organize_listbr_node(nd->v1, listbr, nbseq);
  organize_listbr_node(nd->v2, listbr, nbseq);

}




void organize_listbr(tree* s_tree){

  int i, nbseq;
  branche* listbr;

  listbr=s_tree->listbr;

  nbseq=s_tree->nbseq;

  for(i=0;i<2*nbseq-3;i++){
    if(listbr[i]->bout1==root->v2 && listbr[i]->bout2==root->v1){
      invert_branch(listbr[i]);
      break;
    }
  }

  organize_listbr_node(root->v1, listbr, nbseq);
  organize_listbr_node(root->v2, listbr, nbseq);
}



int setpattern_site(noeud nd, int site){
  int i, son1, son2;

  if(nd->v1==NULL){
    nd->nbpattern=4;
    if(nd->seq[site]=='A') return 0;
    else if(nd->seq[site]=='C') return 1;
    else if(nd->seq[site]=='G') return 2;
    return 3;
  }

  son1=setpattern_site(nd->v1, site);
  son2=setpattern_site(nd->v2, site);

  for(i=0;i<nd->nbpattern;i++){
    if(nd->patternson1[i]==son1 && nd->patternson2[i]==son2)
      return i;
  }
  
  nd->patternson1[nd->nbpattern]=son1;
  nd->patternson2[nd->nbpattern]=son2;
  nd->nbpattern++;
  return nd->nbpattern-1;

}



void setpattern(tree* s_tree){
  int i, j, k, l, nbnd;
  static int nbappel=0;

  nbappel++;

  nbnd=2*s_tree->nbseq-2;
  for(i=0;i<nbnd;i++)
    s_tree->node[i]->nbpattern=0;
  root->nbpattern=0;
  
  for(i=0;i<s_tree->lgseq;i++)
    setpattern_site(root, i);
  
}



int setunderlyingbr(noeud nd, branche br, int i){

  int value=0;

  if(nd->v1==NULL){
    nd->underlyingbr[i]=0;
    return 0;
  }

  if(setunderlyingbr(nd->v1, br, i)) value=1;
  else if(setunderlyingbr(nd->v2, br, i)) value=2;
  else if(nd==br->bout1 && nd->v1==br->bout2) value=10;
  else if(nd==br->bout2 && nd->v1==br->bout1) value=10;
  else if(nd==br->bout1 && nd->v2==br->bout2) value=20;
  else if(nd==br->bout2 && nd->v2==br->bout1) value=20;
  else if(nd->v1==br->bout1 && nd->v2==br->bout2) value=12;  /* nd=root, br=root branch */
  else if(nd->v2==br->bout1 && nd->v1==br->bout2) value=12;  /* nd=root, br=root branch */
  
  nd->underlyingbr[i]=value;

  return value;
}



int setunderlyinggc(noeud nd, int nb){

  int i;

  for(i=0;i<2*nb-3;i++){
    nd->underlyinggc[i]=nd->underlyingbr[i];
    if(nd->underlyinggc[i]==12) nd->underlyinggc[i]=10;
  }
  if(nd==root)
    nd->underlyinggc[2*nb-3]=20;
  else
    nd->underlyinggc[2*nb-3]=0;

  if(nd->v1==NULL) return;

  setunderlyinggc(nd->v1, nb);
  setunderlyinggc(nd->v2, nb);

}




void set_numerobr(noeud nd, branche* listbr, int nbbr){

  int i;
  noeud br1, br2;

  if(nd==root) 
    nd->numerobr=-1;
  else{
    for(i=0;i<nbbr;i++){
      br1=listbr[i]->bout1;
      br2=listbr[i]->bout2;
      if(nd->v3==root){
	if((root->v1==br1 && root->v2==br2)||(root->v2==br1 && root->v1==br2)){
	  nd->numerobr=i;
	  break;
	}
      }
      else if((nd->v3==br1 && nd==br2)||(nd==br1 && nd->v3==br2)){
	nd->numerobr=i;
	break;
      } 
    }
  }
  
  if(!nd->v1) return;

  set_numerobr(nd->v1, listbr, nbbr);
  set_numerobr(nd->v2, listbr, nbbr);

}






void root_s(branche br, double frac1){
  noeud nd1, nd2;
  double l;

  nd1=br->bout1;
  nd2=br->bout2;
  l=br->length;

  root->v1=nd1; root->v2=nd2; root->l1=l*frac1; root->l2=l*(1-frac1);
  root->alive1=root->alive2=1; root->alive3=-1;
  nd1->remain=nd2; nd2->remain=nd1; 
  nd1->lremain=l; nd2->lremain=l;
  if (nd1->v1==nd2) { nd1->v1=root; nd1->l1=l*frac1; nd1->alremain=nd1->alive1; nd1->alive1=1;}
  else if (nd1->v2==nd2) { nd1->v2=root; nd1->l2=l*frac1; nd1->alremain=nd1->alive2; nd1->alive2=1; }
  else if (nd1->v3==nd2) { nd1->v3=root; nd1->l3=l*frac1; nd1->alremain=nd1->alive3; nd1->alive3=1; }
  else printf("Error rooting.\n");
  if (nd2->v1==nd1) { nd2->v1=root; nd2->l1=l*(1-frac1); nd2->alremain=nd2->alive1; nd2->alive1=1; }
  else if (nd2->v2==nd1) { nd2->v2=root; nd2->l2=l*(1-frac1); nd2->alremain=nd2->alive2; nd2->alive2=1; }
  else if (nd2->v3==nd1) { nd2->v3=root; nd2->l3=l*(1-frac1); nd2->alremain=nd2->alive3; nd2->alive3=1; }
  else printf("Error rooting.\n");
}





/* unroot */
/* Unroot t_tree treer. All parameters are modified (branch lengths, names, ...). */
/* Root info may be kept in list1-list2 l1-l2 fracroot1 if l1 is non NULL and list1 and list2  */
/* are allocated as nbtaxa-sized char* arrays. */

int unroot(int** treer, int notu, double* lgbi, double* lgbp, double* bootvals, char** nom, char** list1, char** list2, int* l1, int* l2, double* fracroot1){

int i, j, j1, j2, br_a_virer=-1, bp_a_modifier=-1, bi_a_modifier=-1, somme1, somme2, drapeau;
  double l;

  for(j=0;j<notu-2;j++) if(treer[0][j]!=0) return 1;
  for(j1=0;j1<notu-3;j1++){
    somme1=0;
    for(i=0;i<notu+1;i++) somme1 += treer[i][j1];
    if (somme1==notu-1){ 
      br_a_virer=j1; 
      if(lgbp){
        l=lgbi[j1]; 
        for(i=1;i<notu+1;i++) if(treer[i][j1]==0) {bp_a_modifier=i; break;} 
      }
      break;
    }
    for(j2=j1+1;j2<notu-2;j2++){
      somme2=0;
      for(i=0;i<notu+1;i++) somme2+=treer[i][j2];
      if (somme1==notu-somme2){
	drapeau=0;
	for(i=1;i<notu;i++) if(abs(treer[i][j1]-treer[i][j2])!=1) drapeau=1;
	if (drapeau==0){ 
          br_a_virer=j1; 
          if(lgbi) { bi_a_modifier=j2; l=lgbi[j1]; }
	  break; 
	}
      }

      if (br_a_virer!=-1) break;
    }
    if (br_a_virer!=-1) break;
  }


  if(l1){
    j1=j2=0;
    for(i=1;i<=notu;i++){
      if(treer[i][br_a_virer]==0) {list1[j1]=nom[i]; j1++;}
      else {list2[j2]=nom[i]; j2++;}
    }
    if(j1<1 || j2<1) printf("probleme\n");
    *l1=j1; *l2=j2;
    if (bp_a_modifier!=-1) *fracroot1=lgbp[bp_a_modifier]/(lgbp[bp_a_modifier]+l);
    if (bi_a_modifier!=-1) *fracroot1=lgbi[bi_a_modifier]/(lgbi[bi_a_modifier]+l);
  }
  

  if (bp_a_modifier!=-1) lgbp[bp_a_modifier]+=l;
  if (bi_a_modifier!=-1) lgbi[bi_a_modifier]+=l;
  for(i=br_a_virer;i<notu-3;i++){
    if (lgbi) lgbi[i]=lgbi[i+1];
    if (bootvals) bootvals[i]=bootvals[i+1];
  }
  if(lgbp) for(i=0;i<notu;i++) lgbp[i]=lgbp[i+1];
  if(nom){
    free(nom[0]);
    for(i=0;i<notu;i++) nom[i]=nom[i+1];
    nom[notu]=check_alloc(MAXLNAME+1, sizeof(char));
  }
		
  for(i=1;i<notu+1;i++)
    for(j=0;j<notu-2;j++)
      treer[i-1][j]=treer[i][j];

  for(j=br_a_virer+1;j<notu-2;j++)
    for(i=0;i<notu;i++)
      treer[i][j-1]=treer[i][j];

  return 0;
}




/* equallist */
/* Check the list of taxa down a branch. */

int equallist(noeud from, noeud nd, int nb, char** list){
  noeud n1, n2;
  int i, flag;
  

  if(from==nd->v1) { n1=nd->v2; n2=nd->v3; }
  else if(from==nd->v2) { n1=nd->v1; n2=nd->v3; }
  else if(from==nd->v3) { n1=nd->v1; n2=nd->v2; }
  else { printf("Erreur\n"); }

  if(!n1){
    for(i=0;i<nb;i++){
      flag=0; 
      if(samename(nd->nom, list[i])) {flag=1; break; }
    }
    return flag;
  }
  
  return(equallist(nd, n1, nb, list) && equallist(nd, n2, nb, list));
}





/* isbranch */
/* Check wether lists of taxa name list1 and list2 (respectively l1 */
/* and l2 taxa) exactly match the lists of taxa laying left-side  */
/* and right-side branch br. */

int isbranch(branche br, int l1, char** list1, int l2, char** list2){

  if( equallist(br->bout1, br->bout2, l1, list1) && equallist(br->bout2, br->bout1, l2, list2)) return 2;
  if( equallist(br->bout2, br->bout1, l1, list1) && equallist(br->bout1, br->bout2, l2, list2)) return 1;
  return 0;
}
  



/* taxa_count */
/* Return the number of taxa in c_tree ctree. */

int taxa_count(char* ctree){
  char* prov;
  int nbpo=0, nbpf=0, nbv=0;

  prov=ctree;
  while(*prov!='[' && *prov!='(' && *prov) prov++;
  if(*prov=='\0'){
    printf("Trees should begin with '[' or '('\n");
    return -1;
  }
  if(*prov=='['){
    while(*prov!=']') prov++;
    while(*prov!='(') prov++;
  }

/* debut arbre */

  while(*prov && *prov!=';'){
    if(*prov=='(') nbpo++;
    if(*prov==')') nbpf++;
    if(nbpo-nbpf==1 && *prov==',') nbv++;
    prov++;
  }
  
  if(nbv==1){ /* rooted */
    return (nbpo+1);
  }
  if(nbv==2){ /* unrooted */
    return (nbpo+2);
  }
  printf("Bad number of commas ',' in tree\n");
  return -1;
}




noeud* ctos(char* ctree, int* notu, int *tree_is_rooted){

  int i, nb, **ttree, l1, l2, nbbi;
  double *lgbp, *lgbi, *bootvals, fracroot1; 
  char **nom, racine, *list1[MAXNSP], *list2[MAXNSP];
  noeud* stree;
  branche* br_list;

  nb=taxa_count(ctree);

  lgbp=(double*)check_alloc(nb+1, sizeof(double));
  lgbi=(double*)check_alloc(nb+1, sizeof(double));
  bootvals=(double*)check_alloc(nb+1, sizeof(double));
  /*list1=check_alloc(2*nb+2, sizeof(char*));
  list2=check_alloc(2*nb+2, sizeof(char*));*/

  nom=check_alloc(nb+1, sizeof(char*));
  ttree=check_alloc(nb+1, sizeof(int**));
  for(i=0;i<nb+1;i++){
    nom[i]=(char*)check_alloc(MAXLNAME+1, sizeof(char));
    ttree[i]=(int*)check_alloc(nb, sizeof(int));
  }
  stree=(noeud*)check_alloc(2*nb, sizeof(noeud));
  br_list=(branche*)check_alloc(2*nb, sizeof(branche));
  for(i=0;i<2*nb-2;i++){
    br_list[i]=(struct branche*)check_alloc(1, sizeof(struct branche));
  }

  nb=ctot(ctree, ttree, lgbi, lgbp, bootvals, nom, &racine, &nbbi);

  if(racine=='r'){
    unroot(ttree, nb, lgbi, lgbp, bootvals, nom, list1, list2, &l1, &l2, &fracroot1);
  }
  ttos(ttree, nb, lgbi, lgbp, bootvals, nom, stree);
  makelistbr_unrooted(NULL, stree[0], br_list, NULL, 2*nb-3);
  root=create_node(NULL, NULL, NULL, -1., -1., -1., -1., -1., -1., "ROOT");

  if(racine=='r'){
    *tree_is_rooted=1;
    for(i=0;i<2*nb-3;i++){
      if(isbranch(br_list[i], l1, list1, l2, list2)==1){
 	root_s(br_list[i], fracroot1);
	break;
      }
      if(isbranch(br_list[i], l1, list1, l2, list2)==2){
 	root_s(br_list[i], 1.-fracroot1);
	break;
      }
    }
    if(i==2*nb-3) printf("Erreur racinage\n");
  }
  else{
    *tree_is_rooted=0;
    root_s(br_list[0], 0.5);
  }
  organize_tree(NULL, root);
  *notu=nb;

  free(lgbp);
  free(lgbi);
  free(bootvals);

  for(i=0;i<nb;i++){
    free(ttree[i]); free(nom[i]); 
  }
  free(ttree); free(nom);
  /* free(list1); free(list2); */

  for(i=0;i<2*nb-2;i++)
    free(br_list[i]);
  free(br_list);

  return stree;

}




void reajust_comp_node(noeud nd){

  double somme, *a, *c, *g, *t;
  
  if(nd==root){
    a=&(nd->a); c=&(nd->c); g=&(nd->g); t=&(nd->t); 
  }
  else{
    a=&(nd->aeq); c=&(nd->ceq); g=&(nd->geq); t=&(nd->teq); 
  }

  somme=*a+*c+*g+*t;
  *a/=somme;
  *c/=somme;
  *g/=somme;
  *t/=somme;
}



/* lslgbr */
/* Compute the least-square estimates of branch lengths of tree btree */
/* according to effective distances dist. Branch length are returned */
/* via a double array, ordered as in btree. sce returned value is the sum */
/* of squares of residuals between patristic and effective distances. */
/* totallength returned value is the sum of all branch lengths. */
/* Negative branch lengths are forced to 0 in sce and tot computing. */

double* lslgbr(int** btree, int n, double** dist, double* sce, double* totallength){
  int i, j, k, ii, nbdist, nbbr, *ind1, *ind2, b1, b2;
  double **ta, **B, **invB, **prod1, *lgbr, *cdist, *pdist;

  nbdist=n*(n-1)/2; nbbr=2*n-3;
  ta=(double**)check_alloc(nbbr, sizeof(double*));
  B=(double**)check_alloc(nbbr, sizeof(double*));
  invB=(double**)check_alloc(nbbr, sizeof(double*));
  prod1=(double**)check_alloc(nbbr, sizeof(double*));

  /* move double** dist to double* cdist */

  cdist=(double*)check_alloc(nbdist, sizeof(double));
  ind1=(int*)check_alloc(nbdist, sizeof(int));
  ind2=(int*)check_alloc(nbdist, sizeof(int));

  ii=0;
  for(i=0;i<n-1;i++){
    for(j=i+1;j<n;j++){
      cdist[ii]=dist[i][j];
      ind1[ii]=i; ind2[ii]=j;
      ii++;
    }
  }

  /* build ta */

  for(i=0;i<nbbr;i++) ta[i]=(double*)check_alloc(nbdist, sizeof(double));

  for(i=0;i<n;i++) {
    for(j=0;j<nbdist;j++){
      if(ind1[j]==i || ind2[j]==i) ta[i][j]=1.; else ta[i][j]=0.;
    }
  }
  for(i=n;i<nbbr;i++){
    for(j=0;j<nbdist;j++){
      b1=testbit(btree[i-n],ind1[j]+1); b2=testbit(btree[i-n],ind2[j]+1);
      if((b1 && !b2) || (b2 && !b1)) ta[i][j]=1.; else ta[i][j]=0.;
    }
  }

  /* compute B */

  for(i=0;i<nbbr;i++) B[i]=(double*)check_alloc(nbbr, sizeof(double));
  for(i=0;i<nbbr;i++) invB[i]=(double*)check_alloc(nbbr, sizeof(double));

  for(i=0;i<nbbr;i++) {
    for(j=0;j<nbbr;j++){
      B[i][j]=0.;
      for(k=0;k<nbdist;k++) B[i][j]+=ta[i][k]*ta[j][k];
    }
  }

  /* compute lgbr and pdist */

  for(i=0;i<nbbr;i++) prod1[i]=(double*)check_alloc(nbdist, sizeof(double));
  lgbr=(double*)check_alloc(nbbr, sizeof(double));

  if(!invmat(B, nbbr, invB)) return NULL;
  matmat(invB, nbbr, nbbr, ta, nbbr, nbdist, prod1);
  matvect(prod1, nbbr, nbdist, cdist, nbdist, lgbr);
  pdist=vectmat(lgbr, nbbr, ta, nbbr, nbdist);

  /* compute sce and total length */
 
  if(sce){
    *sce=0.;
    for(i=0;i<nbdist;i++) 
      if(pdist[i]>0.) *sce+=(pdist[i]-cdist[i])*(pdist[i]-cdist[i]);
      else *sce+=cdist[i]*cdist[i];
  }

  if(totallength){
    *totallength=0.;
    for(i=0;i<nbbr;i++) if(lgbr[i]>0.) *totallength+=lgbr[i];
  }

  free(cdist); free(ind1); free(ind2);
  for(i=0;i<nbbr;i++){
    free(ta[i]);
    free(B[i]);
    free(invB[i]);
    free(prod1[i]);
  }
  free(ta); free(B); free(invB); free(prod1);
  free(pdist);
  return lgbr;
}




int congruent(noeud nd, noeud from, int** btree, char** name, int nbseq, int j){

  noeud to1, to2;
  int i, c1, c2;


  if(nd==root){
    if(from==root->v1) return congruent(root->v2, root, btree, name, nbseq, j);
    else if(from==root->v2) return congruent(root->v1, root, btree, name, nbseq, j);
  }

  if(nd->v1==NULL){
    for(i=0;i<nbseq;i++){
      if(samename(nd->nom, name[i])){
	if(testbit(btree[j], i+1)) return 1; else return 0;
      }
    }
    if(i==nbseq){ printf("erreur cong\n"); exit(EXIT_FAILURE); }
  }

  if(from==nd->v1) {to1=nd->v2; to2=nd->v3;}
  else if(from==nd->v2) {to1=nd->v1; to2=nd->v3;}
  else if(from==nd->v3) {to1=nd->v1; to2=nd->v2;}
  else if(from==root->v1 && nd==root->v2) {to1=nd->v1; to2=nd->v2;}
  else if(from==root->v2 && nd==root->v1) {to1=nd->v1; to2=nd->v2;}
  else { printf("erreur cong2 : nd=%s , from=%s\n", nd->nom, from->nom); exit(EXIT_FAILURE);}

  c1=congruent(to1, nd, btree, name, nbseq, j);
  c2=congruent(to2, nd, btree, name, nbseq, j);

  if(c1==c2) return c1;
  return 10;

}



void reset_bl(tree* s_tree, noeud nd){
   
  double l;

  if(nd==NULL) return;


  if(nd!=root){
    l=s_tree->listbr[nd->numerobr]->length;
    if(nd==root->v1)
      nd->l3=root->l1=l*s_tree->froot;
    else if(nd==root->v2)
      nd->l3=root->l2=l*(1.-s_tree->froot);
    else{
      nd->l3=l;
      if(nd==nd->v3->v1) nd->v3->l1=l;
      if(nd==nd->v3->v2) nd->v3->l2=l;
    }
  }
  
  reset_bl(s_tree, nd->v1);
  reset_bl(s_tree, nd->v2);
}


void set_parambl(tree* s_tree, noeud nd){

  double l;

  if(nd==NULL) return;


  if(nd!=root){
    if(nd->v3==root) l=root->l1+root->l2;
    else l=nd->l3;
    s_tree->listbr[nd->numerobr]->length=l;

  }

  set_parambl(s_tree, nd->v1);
  set_parambl(s_tree, nd->v2);
}







int* btltostl(tree* stree, int** btree, char** name, double* treelg){

  int i, j, nbseq, left, right;
  branche br;
  noeud boutterm;

  nbseq=stree->nbseq;

  for(i=0;i<2*nbseq-3;i++){
    br=stree->listbr[i];
    if(br->bout1->v1==NULL || br->bout2->v1==NULL){ /* terminal */
      boutterm=(br->bout1->v1)?br->bout2:br->bout1;
      for(j=0;j<nbseq;j++){
	if(samename(boutterm->nom, name[j])){
	  br->length=treelg[j];
    	  break;
	}
      }
      if(j==nbseq) printf("erreur1\n");
    }
    else{					/* internal */
      for(j=0;j<nbseq-3;j++){
        left=congruent(br->bout1, br->bout2, btree, name, nbseq, j);
        right=congruent(br->bout2, br->bout1, btree, name, nbseq, j);
	if((left==0 && right==1) || (left==1 && right==0)){
	  br->length=treelg[nbseq+j];
	  break;
	}
      }
      if(j==nbseq-3) printf("erreur2\n");
    }
  }
  reset_bl(stree, root);
}



/* unroot_c */
/* Remove root from rooted c_tree carbre. */

int unroot_c(char* carbre){
  int i=0, diff;
  
  while((carbre[i]==' ' || carbre[i]=='\n' || carbre[i]=='\t') && carbre[i]) i++;
  if(carbre[i]!='[' && carbre[i]!='(') {
    printf("Tree first char must be ( or [\n"); 
    return -1; 
  }

  if(carbre[i]=='[') while(carbre[i]!=']' && carbre[i]) i++;
  if(!carbre[i]){ printf("Unmatched '[' ']'\n"); return -1; }
  while(carbre[i]!='(' && carbre[i]) i++;
  if(!carbre[i]) {printf("No initial parenthesis\n"); return -1; }
 
  i++;
  while(carbre[i]!='(' && carbre[i] && carbre[i]!=';') i++;
  if(!carbre[i] || carbre[i]==';') return 0;
  carbre[i]=' ';
  diff=0; 
  while(diff!=-1 && carbre[i] && carbre[i]!=';'){
    if(carbre[i]=='(')  diff++; 
    if(carbre[i]==')')  diff--; 
    i++;
  }

  if(!carbre[i] || carbre[i]==';') return 0;
  carbre[i-1]=' ';
  while(carbre[i]!=',' && carbre[i]!=')') {
    carbre[i]=' ';
    i++;
  }

  return 1;
}  
  


void init_bl(tree* stree){

  int i, j, nbseq, lliste, **btree;
  char* ctree;
  double* newbl, sumalive=0., sum=0., rap;

  nbseq=stree->nbseq;
  ctree=(char*)check_alloc(50*nbseq, sizeof(char));
  lliste=(nbseq+lmot-1)/lmot;
  btree=(int**)check_alloc(nbseq-3, sizeof(int*));
  for(i=0;i<nbseq-3;i++)
    btree[i]=(int*)check_alloc(lliste, sizeof(int));
  newbl=(double*)check_alloc(2*nbseq-3, sizeof(double));
  
  stoc(stree->node, 1, nbseq, ctree, 0);
  unroot_c(ctree);
  ctob(ctree, stree->names, btree);

  newbl=lslgbr(btree, nbseq, stree->dist, NULL, NULL);

  for(i=0;i<2*nbseq-3;i++){
    if(newbl[i]<lmin) newbl[i]=lmin;
    if(newbl[i]>lmax) newbl[i]=lmax;
  }

  for(i=0;i<2*nbseq-3;i++){
    sum+=newbl[i];
    if(stree->alivebr[i]) sumalive+=newbl[i];
  }
  rap=sum/sumalive;
  for(i=0;i<2*nbseq-3;i++)
    if(stree->alivebr[i])
      newbl[i]*=rap;
  

  btltostl(stree, btree, stree->names, newbl);


  for(i=0;i<nbseq-3;i++) free(btree[i]);
  free(btree);
  free(ctree);
  free(newbl);
}



/* freq_obs */
/* Write at address freq observed frequencies of 16 di-nucleotides XY */
/* (X= A,C,G,T  ,  Y=A,C,G,T) X and Y being homologous nucleotides of sequences */
/* seq1 and seq2. Alphabetic order is used : freq[0]=AA frequency, freq[1]=AC, */
/* ..., freq[15]=TT. */

int freq_obs(char* seq1, char* seq2, double* freq){

  int i, lgseq, lgseqvrai;

  lgseq=(int)strlen(seq1);
  if (lgseq!=(int)strlen(seq2)){
    printf ("Longueurs inegales.\n");
    return -1;
  }
  for(i=0;i<16;i++) freq[i]=0;
  lgseqvrai=lgseq;
  for(i=0;i<lgseq;i++){
    switch(seq1[i]){
    case 'A':
      switch(seq2[i]){
      case 'A' : freq[0]++; break;
      case 'C' : freq[1]++; break;
      case 'G' : freq[2]++; break;
      case 'T' : freq[3]++; break;
      default : lgseqvrai --; break;
      }
      break;
    case 'C':
      switch(seq2[i]){
      case 'A' : freq[4]++; break;
      case 'C' : freq[5]++; break;
      case 'G' : freq[6]++; break;
      case 'T' : freq[7]++; break;
      default : lgseqvrai --; break;
      }
      break;
    case 'G':
      switch(seq2[i]){
      case 'A' : freq[8]++; break;
      case 'C' : freq[9]++; break;
      case 'G' : freq[10]++; break;
      case 'T' : freq[11]++; break;
      default : lgseqvrai --; break;
      }
      break;
    case 'T':
      switch(seq2[i]){
      case 'A' : freq[12]++; break;
      case 'C' : freq[13]++; break;
      case 'G' : freq[14]++; break;
      case 'T' : freq[15]++; break;
      default : lgseqvrai --; break;
      }
      break;
    default :
      lgseqvrai --;
    }
  }
  if(lgseqvrai!=0){
    for(i=0;i<16;i++) freq[i]/=lgseqvrai;
    return 1;
  }
  else return 0;
}



/* alsurbet_kim */
/* Return the quotient between transition substitution rate (alpha) and */
/* transversion substitution rate (beta) according to Kimura 1980 */
/* 2-parameter model between sequence seq1 and seq2. */

double alsurbet_kim(char* seq1, char* seq2, double* freq){

  double alphakim, betakim, P, Q;

  if(!freq_obs(seq1, seq2, freq)) return -1.;

  Q=freq[1]+freq[3]+freq[4]+freq[6]+freq[9]+freq[11]+freq[12]+freq[14];
  P=1-Q-freq[0]-freq[5]-freq[10]-freq[15];

  if (2*Q>=1.) return -1.;
  if (2*P+Q>=1.) return -1.;

  betakim=-0.125*log(1.-2*Q);
  alphakim=-0.25*log(1.-2*P-Q)+0.125*log(1.-2*Q);

  if (betakim==0.) return -1.;
 	
  return(alphakim/betakim);
}




/* gg95 */
/* Galtier & Gouy distance */

double gg95(char* seq1, char* seq2, double a){

  double k11, k12, k21, k22, t0, rt, *freq, Q, t1, t2, e;
  int i;


  freq=check_alloc(17, sizeof(double));
  if(!freq_obs(seq1, seq2, freq)) return -1.;

  Q=freq[1]+freq[3]+freq[4]+freq[6]+freq[9]+freq[11]+freq[12]+freq[14];
  t1=0.;
  for(i=4;i<=11;i++) t1+=freq[i];
  t2=freq[1]+freq[2]+freq[5]+freq[6]+freq[9]+freq[10]+freq[13]+freq[14];
			
  if(2*Q>=1.) return -1.;
  rt=-0.5*log(1-Q*2);
	
  t0=(t1+t2)/2;

  e=1-exp(-rt*(a+1.)/2);
  k11=(0.5+a*t1*(1-t1))*rt ; 
  k12=(a/(a+1.)) * ((t0-t1)*(1-2*t1)) * e;	
  k21=(0.5+a*t2*(1-t2))*rt;
  k22=(a/(a+1.)) * ((t0-t2)*(1-2*t2)) * e;

  free(freq);

  return(k11+k12+k21+k22);

}





void init4bases_constant(tree* s_tree){

  int i, j, nbseq, lgseq, weightsum;
  double aglob=0., cglob=0., gglob=0., tglob=0.;

  nbseq=s_tree->nbseq;
  lgseq=s_tree->lgseq;
  weightsum=s_tree->weightsum;

  for(i=0;i<nbseq;i++){
    for(j=0;j<lgseq;j++){
      if(s_tree->node[i]->seq[j]=='A') aglob+=s_tree->weight[j]; 
      if(s_tree->node[i]->seq[j]=='C') cglob+=s_tree->weight[j];
      if(s_tree->node[i]->seq[j]=='G') gglob+=s_tree->weight[j];
      if(s_tree->node[i]->seq[j]=='T') tglob+=s_tree->weight[j];
    }
  }
  aglob/=(weightsum*nbseq); cglob/=(weightsum*nbseq); gglob/=(weightsum*nbseq); tglob/=(weightsum*nbseq);
  for(i=0;i<2*nbseq-2;i++){
    s_tree->node[i]->aeq=aglob;
    s_tree->node[i]->ceq=cglob;
    s_tree->node[i]->geq=gglob;
    s_tree->node[i]->teq=tglob;
    reajust_comp_node(s_tree->node[i]);
  }  
  root->aeq=aglob;
  root->ceq=cglob;
  root->geq=gglob;
  root->teq=tglob;
}




void init4bases_simplemean(noeud from, noeud nd, int lgseq, double weightsum, double* weight){

  noeud nd1, nd2;
  double a=0., c=0., g=0., t=0.;

  if (from==nd->v1) {nd1=nd->v2; nd2=nd->v3;}
  else if (from==nd->v2) {nd1=nd->v1; nd2=nd->v3;}
  else if (from==nd->v3) {nd1=nd->v1; nd2=nd->v2;}
  else { printf("error initiating gc.\n"); exit(EXIT_FAILURE); }

  if (nd1==NULL) {
    int i; 
    for(i=0;i<lgseq;i++){
      if(nd->seq[i]=='A') a+=weight[i];
      if(nd->seq[i]=='C') c+=weight[i];
      if(nd->seq[i]=='G') g+=weight[i];
      if(nd->seq[i]=='T') t+=weight[i];
    }
    a/=weightsum; c/=weightsum; g/=weightsum; t/=weightsum;
    nd->aeq=a;
    nd->ceq=c;
    nd->geq=g;
    nd->teq=t;
    return;
  }

  init4bases_simplemean(nd, nd1, lgseq, weightsum, weight);
  init4bases_simplemean(nd, nd2, lgseq, weightsum, weight);
  nd->aeq=(nd1->aeq+nd2->aeq)/2;
  nd->ceq=(nd1->ceq+nd2->ceq)/2;
  nd->geq=(nd1->geq+nd2->geq)/2;
  nd->teq=(nd1->teq+nd2->teq)/2;
  reajust_comp_node(nd);
  return ;
}




void init4bases_balanced(noeud from, noeud nd, int lgseq, double weightsum, double* weight){

  noeud nd1, nd2;
  double a=0., c=0., g=0., t=0.;

  if (from==nd->v1) {nd1=nd->v2; nd2=nd->v3;}
  else if (from==nd->v2) {nd1=nd->v1; nd2=nd->v3;}
  else if (from==nd->v3) {nd1=nd->v1; nd2=nd->v2;}
  else { printf("error initiating gc.\n"); exit(EXIT_FAILURE); }

  nd->aeq=nd->ceq=nd->geq=nd->teq=0.25;

  if (nd1==NULL) 
    return;
  
  init4bases_balanced(nd, nd1, lgseq, weightsum, weight);
  init4bases_balanced(nd, nd2, lgseq, weightsum, weight);
}





init4bases_random(noeud nd){

  double ran1=-1., ran2=-1., ran3=-1., ran4=-1., sum;

  while(ran1<=0.) ran1=drand48(); 
  while(ran2<=0.) ran2=drand48(); 
  while(ran3<=0.) ran3=drand48(); 
  while(ran4<=0.) ran4=drand48();
  sum=ran1+ran2+ran3+ran4;

  nd->aeq=ran1/sum;
  nd->ceq=ran2/sum;
  nd->geq=ran3/sum;
  nd->teq=ran4/sum;

  reajust_comp_node(nd);
  if(nd->v1==NULL) return;

  init4bases_random(nd->v1);
  init4bases_random(nd->v2);

}





void compute_proba_t3_g(noeud nd, double a, int OPTIMIZE_LENGTH, int OPTIMIZE_TITV, int OPTIMIZE_ROOT, int OPTIMIZE_COMP, int cl1, int cl2, double rr){

/* Parameters:    branch lengths  :	 nd->l1 , nd->l2 				   */
/* 		  ti/tv ratio 	  :	 a 						   */
/* 		  eq. GC contents :	 nd->v1->geq+nd->v1->ceq , nd->v2->geq+nd->v2->ceq */
/* 		  root location	  : 	 nd->l1/(nd->l1+nd->l2)				   */
/* 		  ancestral GC    :      useless 					   */
/*		  relative rate	  :      rr */

  double k;					/* (1+a)/2 */
  double ert1, ert2;				/* exp(l1) , exp(l2) */
  double ekrt1, ekrt2;				/* exp(k*l1) , exp(k*l2) */
  double upert1, upert2;			/* 1+exp(l1) , 1+exp(l2) */
  double umert1, umert2;			/* 1-exp(l1) , 1-exp(l2) */
  double gc1, gc2;				
  double gc1s2, gc2s2;				/* gc1/2 , gc2/2 */
  double umgc1, umgc2, umgc1s2, umgc2s2;	/* 1-gc1 , 1-gc2 , (1-gc1)/2 , (1-gc2)/2 */

  double d, umd;

  int i, j, kk, ll, mm, nn;



  if(nd->v1==NULL)
    return;



tdeb=time(NULL);

  gc1=nd->v1->geq+nd->v1->ceq; gc1s2=gc1/2.; umgc1=1.-gc1; umgc1s2=umgc1/2;
  gc2=nd->v2->geq+nd->v2->ceq; gc2s2=gc2/2.; umgc2=1.-gc2; umgc2s2=umgc2/2;
  k=(a+1)/2;
  ert1=exp(-nd->l1*rr);
  upert1=1.+ert1;
  umert1=1.-ert1;
  ekrt1=exp(-k*nd->l1*rr);
  ert2=exp(-nd->l2*rr);
  upert2=1.+ert2;
  umert2=1.-ert2;
  ekrt2=exp(-k*nd->l2*rr);

  if(nd==root){
    d=nd->l1/(nd->l1+nd->l2); 
    umd=1.-d;
  }
  else d=umd=1.;




		/* probabilities */

  if(nd->alive1){  
  nd->p1[cl1][cl2][0][0]=umgc1s2*upert1+gc1*ekrt1; 			
  nd->p1[cl1][cl2][0][1]=gc1s2*umert1;				
  nd->p1[cl1][cl2][0][2]=gc1s2*upert1-gc1*ekrt1;			
  nd->p1[cl1][cl2][0][3]=umgc1s2*umert1;				
  nd->p1[cl1][cl2][1][0]=umgc1s2*umert1;				
  nd->p1[cl1][cl2][1][1]=gc1s2*upert1+umgc1*ekrt1;		
  nd->p1[cl1][cl2][1][2]=gc1s2*umert1;				
  nd->p1[cl1][cl2][1][3]=umgc1s2*upert1-umgc1*ekrt1;		
  nd->p1[cl1][cl2][2][0]=umgc1s2*upert1-umgc1*ekrt1;		
  nd->p1[cl1][cl2][2][1]=gc1s2*umert1;				
  nd->p1[cl1][cl2][2][2]=gc1s2*upert1+umgc1*ekrt1;		
  nd->p1[cl1][cl2][2][3]=umgc1s2*umert1;				
  nd->p1[cl1][cl2][3][0]=umgc1s2*umert1;				
  nd->p1[cl1][cl2][3][1]=gc1s2*upert1-gc1*ekrt1;			
  nd->p1[cl1][cl2][3][2]=gc1s2*umert1;				
  nd->p1[cl1][cl2][3][3]=umgc1s2*upert1+gc1*ekrt1; 	
  }


if(nd->p1[cl1][cl2][0][1]<0.){
  printf("gc1s2=%f, umert1=%f\n", gc1s2, umert1);
  printf("ert1=%f, l1=%f, rr=%f\n", ert1, nd->l1, rr);
  printf("alive1: %d, alive2: %d\n", nd->alive1, nd->alive2);
}

  if(nd->alive2){
  nd->p2[cl1][cl2][0][0]=umgc2s2*upert2+gc2*ekrt2; 			
  nd->p2[cl1][cl2][0][1]=gc2s2*umert2;				
  nd->p2[cl1][cl2][0][2]=gc2s2*upert2-gc2*ekrt2;			
  nd->p2[cl1][cl2][0][3]=umgc2s2*umert2;				
  nd->p2[cl1][cl2][1][0]=umgc2s2*umert2;				
  nd->p2[cl1][cl2][1][1]=gc2s2*upert2+umgc2*ekrt2;		
  nd->p2[cl1][cl2][1][2]=gc2s2*umert2;				
  nd->p2[cl1][cl2][1][3]=umgc2s2*upert2-umgc2*ekrt2;		
  nd->p2[cl1][cl2][2][0]=umgc2s2*upert2-umgc2*ekrt2;		
  nd->p2[cl1][cl2][2][1]=gc2s2*umert2;				
  nd->p2[cl1][cl2][2][2]=gc2s2*upert2+umgc2*ekrt2;		
  nd->p2[cl1][cl2][2][3]=umgc2s2*umert2;				
  nd->p2[cl1][cl2][3][0]=umgc2s2*umert2;				
  nd->p2[cl1][cl2][3][1]=gc2s2*upert2-gc2*ekrt2;			
  nd->p2[cl1][cl2][3][2]=gc2s2*umert2;				
  nd->p2[cl1][cl2][3][3]=umgc2s2*upert2+gc2*ekrt2;
  } 	
 	



		/* first and second derivatives */


if(OPTIMIZE_TITV){
  double rt1ekrt1s2, rt2ekrt2s2, rt1rt1ekrt1s4, rt2rt2ekrt2s4;

  rt1ekrt1s2=nd->l1*rr*ekrt1/2.;
  rt2ekrt2s2=nd->l2*rr*ekrt2/2.;
  rt1rt1ekrt1s4=rr*rt1ekrt1s2*nd->l1/2.;
  rt2rt2ekrt2s4=rr*rt2ekrt2s2*nd->l2/2.;

  if(nd->alive1){
  nd->dp1[cl1][cl2][0][0][0]=-gc1*rt1ekrt1s2;		
  nd->dp1[cl1][cl2][0][0][1]=0.;				
  nd->dp1[cl1][cl2][0][0][2]=gc1*rt1ekrt1s2;			
  nd->dp1[cl1][cl2][0][0][3]=0.;				
  nd->dp1[cl1][cl2][0][1][0]=0.;				
  nd->dp1[cl1][cl2][0][1][1]=-umgc1*rt1ekrt1s2;		
  nd->dp1[cl1][cl2][0][1][2]=0.;				
  nd->dp1[cl1][cl2][0][1][3]=umgc1*rt1ekrt1s2;		
  nd->dp1[cl1][cl2][0][2][0]=umgc1*rt1ekrt1s2;		
  nd->dp1[cl1][cl2][0][2][1]=0.;				
  nd->dp1[cl1][cl2][0][2][2]=-umgc1*rt1ekrt1s2;		
  nd->dp1[cl1][cl2][0][2][3]=0.;				
  nd->dp1[cl1][cl2][0][3][0]=0.;				
  nd->dp1[cl1][cl2][0][3][1]=gc1*rt1ekrt1s2;			
  nd->dp1[cl1][cl2][0][3][2]=0.;				
  nd->dp1[cl1][cl2][0][3][3]=-gc1*rt1ekrt1s2; 
  }

  if(nd->alive2){
  nd->dp2[cl1][cl2][0][0][0]=-gc2*rt2ekrt2s2; 			
  nd->dp2[cl1][cl2][0][0][1]=0.;				
  nd->dp2[cl1][cl2][0][0][2]=gc2*rt2ekrt2s2;			
  nd->dp2[cl1][cl2][0][0][3]=0.;				
  nd->dp2[cl1][cl2][0][1][0]=0.;				
  nd->dp2[cl1][cl2][0][1][1]=-umgc2*rt2ekrt2s2;		
  nd->dp2[cl1][cl2][0][1][2]=0.;				
  nd->dp2[cl1][cl2][0][1][3]=umgc2*rt2ekrt2s2;		
  nd->dp2[cl1][cl2][0][2][0]=umgc2*rt2ekrt2s2;		
  nd->dp2[cl1][cl2][0][2][1]=0.;				
  nd->dp2[cl1][cl2][0][2][2]=-umgc2*rt2ekrt2s2;		
  nd->dp2[cl1][cl2][0][2][3]=0.;				
  nd->dp2[cl1][cl2][0][3][0]=0.;				
  nd->dp2[cl1][cl2][0][3][1]=gc2*rt2ekrt2s2;			
  nd->dp2[cl1][cl2][0][3][2]=0.;				
  nd->dp2[cl1][cl2][0][3][3]=-gc2*rt2ekrt2s2; 
  }

  if(nd->alive1){
  nd->d2p1[cl1][cl2][0][0][0]=gc1*rt1rt1ekrt1s4; 			
  nd->d2p1[cl1][cl2][0][0][1]=0.;				
  nd->d2p1[cl1][cl2][0][0][2]=-gc1*rt1rt1ekrt1s4;			
  nd->d2p1[cl1][cl2][0][0][3]=0.;				
  nd->d2p1[cl1][cl2][0][1][0]=0.;				
  nd->d2p1[cl1][cl2][0][1][1]=umgc1*rt1rt1ekrt1s4;		
  nd->d2p1[cl1][cl2][0][1][2]=0.;				
  nd->d2p1[cl1][cl2][0][1][3]=-umgc1*rt1rt1ekrt1s4;		
  nd->d2p1[cl1][cl2][0][2][0]=-umgc1*rt1rt1ekrt1s4;		
  nd->d2p1[cl1][cl2][0][2][1]=0.;				
  nd->d2p1[cl1][cl2][0][2][2]=umgc1*rt1rt1ekrt1s4;		
  nd->d2p1[cl1][cl2][0][2][3]=0.;				
  nd->d2p1[cl1][cl2][0][3][0]=0.;				
  nd->d2p1[cl1][cl2][0][3][1]=-gc1*rt1rt1ekrt1s4;			
  nd->d2p1[cl1][cl2][0][3][2]=0.;				
  nd->d2p1[cl1][cl2][0][3][3]=gc1*rt1rt1ekrt1s4; 
  }

  if(nd->alive2){
  nd->d2p2[cl1][cl2][0][0][0]=gc2*rt2rt2ekrt2s4; 			
  nd->d2p2[cl1][cl2][0][0][1]=0.;				
  nd->d2p2[cl1][cl2][0][0][2]=-gc2*rt2rt2ekrt2s4;			
  nd->d2p2[cl1][cl2][0][0][3]=0.;				
  nd->d2p2[cl1][cl2][0][1][0]=0.;				
  nd->d2p2[cl1][cl2][0][1][1]=umgc2*rt2rt2ekrt2s4;		
  nd->d2p2[cl1][cl2][0][1][2]=0.;				
  nd->d2p2[cl1][cl2][0][1][3]=-umgc2*rt2rt2ekrt2s4;		
  nd->d2p2[cl1][cl2][0][2][0]=-umgc2*rt2rt2ekrt2s4;		
  nd->d2p2[cl1][cl2][0][2][1]=0.;				
  nd->d2p2[cl1][cl2][0][2][2]=umgc2*rt2rt2ekrt2s4;		
  nd->d2p2[cl1][cl2][0][2][3]=0.;				
  nd->d2p2[cl1][cl2][0][3][0]=0.;				
  nd->d2p2[cl1][cl2][0][3][1]=-gc2*rt2rt2ekrt2s4;			
  nd->d2p2[cl1][cl2][0][3][2]=0.;				
  nd->d2p2[cl1][cl2][0][3][3]=gc2*rt2rt2ekrt2s4; 
  }
}
else{
    for(i=0;i<4;i++) for(j=0;j<4;j++){
      nd->dp1[cl1][cl2][0][i][j]=0.;
      nd->dp2[cl1][cl2][0][i][j]=0.;
      nd->d2p1[cl1][cl2][0][i][j]=0.;
      nd->d2p2[cl1][cl2][0][i][j]=0.;
    }
}


if(OPTIMIZE_ROOT && nd==root){
  double l, kekrt1, k2ekrt1, kekrt2, k2ekrt2, dert1, dert2, d2ert1, d2ert2;

  l=nd->l1+nd->l2;
  l*=rr;
  kekrt1=k*l*ekrt1;
  k2ekrt1=k*l*kekrt1;
  kekrt2=-k*l*ekrt2;
  k2ekrt2=-k*l*kekrt2;
  dert1=l*ert1; dert2=-l*ert2; d2ert1=l*dert1; d2ert2=-l*dert2;

  if(nd->alive1){
  nd->dp1[cl1][cl2][1][0][0]=-umgc1s2*dert1-gc1*kekrt1; 			
  nd->dp1[cl1][cl2][1][0][1]=gc1s2*dert1;				
  nd->dp1[cl1][cl2][1][0][2]=-gc1s2*dert1+gc1*kekrt1;			
  nd->dp1[cl1][cl2][1][0][3]=umgc1s2*dert1;				
  nd->dp1[cl1][cl2][1][1][0]=umgc1s2*dert1;				
  nd->dp1[cl1][cl2][1][1][1]=-gc1s2*dert1-umgc1*kekrt1;		
  nd->dp1[cl1][cl2][1][1][2]=gc1s2*dert1;				
  nd->dp1[cl1][cl2][1][1][3]=-umgc1s2*dert1+umgc1*kekrt1;		
  nd->dp1[cl1][cl2][1][2][0]=-umgc1s2*dert1+umgc1*kekrt1;		
  nd->dp1[cl1][cl2][1][2][1]=gc1s2*dert1;				
  nd->dp1[cl1][cl2][1][2][2]=-gc1s2*dert1-umgc1*kekrt1;		
  nd->dp1[cl1][cl2][1][2][3]=umgc1s2*dert1;				
  nd->dp1[cl1][cl2][1][3][0]=umgc1s2*dert1;				
  nd->dp1[cl1][cl2][1][3][1]=-gc1s2*dert1+gc1*kekrt1;			
  nd->dp1[cl1][cl2][1][3][2]=gc1s2*dert1;				
  nd->dp1[cl1][cl2][1][3][3]=-umgc1s2*dert1-gc1*kekrt1; 
  }

  if(nd->alive2){
  nd->dp2[cl1][cl2][1][0][0]=-umgc2s2*dert2-gc2*kekrt2; 			
  nd->dp2[cl1][cl2][1][0][1]=gc2s2*dert2;				
  nd->dp2[cl1][cl2][1][0][2]=-gc2s2*dert2+gc2*kekrt2;			
  nd->dp2[cl1][cl2][1][0][3]=umgc2s2*dert2;				
  nd->dp2[cl1][cl2][1][1][0]=umgc2s2*dert2;				
  nd->dp2[cl1][cl2][1][1][1]=-gc2s2*dert2-umgc2*kekrt2;		
  nd->dp2[cl1][cl2][1][1][2]=gc2s2*dert2;				
  nd->dp2[cl1][cl2][1][1][3]=-umgc2s2*dert2+umgc2*kekrt2;		
  nd->dp2[cl1][cl2][1][2][0]=-umgc2s2*dert2+umgc2*kekrt2;		
  nd->dp2[cl1][cl2][1][2][1]=gc2s2*dert2;				
  nd->dp2[cl1][cl2][1][2][2]=-gc2s2*dert2-umgc2*kekrt2;		
  nd->dp2[cl1][cl2][1][2][3]=umgc2s2*dert2;				
  nd->dp2[cl1][cl2][1][3][0]=umgc2s2*dert2;				
  nd->dp2[cl1][cl2][1][3][1]=-gc2s2*dert2+gc2*kekrt2;			
  nd->dp2[cl1][cl2][1][3][2]=gc2s2*dert2;				
  nd->dp2[cl1][cl2][1][3][3]=-umgc2s2*dert2-gc2*kekrt2; 
  }

  if(nd->alive1){
  nd->d2p1[cl1][cl2][1][0][0]=umgc1s2*d2ert1+gc1*k2ekrt1; 			
  nd->d2p1[cl1][cl2][1][0][1]=-gc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][1][0][2]=gc1s2*d2ert1-gc1*k2ekrt1;			
  nd->d2p1[cl1][cl2][1][0][3]=-umgc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][1][1][0]=-umgc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][1][1][1]=gc1s2*d2ert1+umgc1*k2ekrt1;		
  nd->d2p1[cl1][cl2][1][1][2]=-gc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][1][1][3]=umgc1s2*d2ert1-umgc1*k2ekrt1;		
  nd->d2p1[cl1][cl2][1][2][0]=umgc1s2*d2ert1-umgc1*k2ekrt1;		
  nd->d2p1[cl1][cl2][1][2][1]=-gc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][1][2][2]=gc1s2*d2ert1+umgc1*k2ekrt1;		
  nd->d2p1[cl1][cl2][1][2][3]=-umgc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][1][3][0]=-umgc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][1][3][1]=gc1s2*d2ert1-gc1*k2ekrt1;			
  nd->d2p1[cl1][cl2][1][3][2]=-gc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][1][3][3]=umgc1s2*d2ert1+gc1*k2ekrt1; 
  }

  if(nd->alive2){
  nd->d2p2[cl1][cl2][1][0][0]=umgc2s2*d2ert2+gc2*k2ekrt2; 			
  nd->d2p2[cl1][cl2][1][0][1]=-gc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][1][0][2]=gc2s2*d2ert2-gc2*k2ekrt2;			
  nd->d2p2[cl1][cl2][1][0][3]=-umgc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][1][1][0]=-umgc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][1][1][1]=gc2s2*d2ert2+umgc2*k2ekrt2;		
  nd->d2p2[cl1][cl2][1][1][2]=-gc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][1][1][3]=umgc2s2*d2ert2-umgc2*k2ekrt2;		
  nd->d2p2[cl1][cl2][1][2][0]=umgc2s2*d2ert2-umgc2*k2ekrt2;		
  nd->d2p2[cl1][cl2][1][2][1]=-gc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][1][2][2]=gc2s2*d2ert2+umgc2*k2ekrt2;		
  nd->d2p2[cl1][cl2][1][2][3]=-umgc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][1][3][0]=-umgc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][1][3][1]=gc2s2*d2ert2-gc2*k2ekrt2;			
  nd->d2p2[cl1][cl2][1][3][2]=-gc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][1][3][3]=umgc2s2*d2ert2+gc2*k2ekrt2;   
  }
}
else{
  for(i=0;i<4;i++) for(j=0;j<4;j++){ 
    nd->dp1[cl1][cl2][1][i][j]=0.;
    nd->dp2[cl1][cl2][1][i][j]=0.;
    nd->d2p1[cl1][cl2][1][i][j]=0.;
    nd->d2p2[cl1][cl2][1][i][j]=0.;
  }
}


if(OPTIMIZE_LENGTH){
  double kekrt1, k2ekrt1, kekrt2, k2ekrt2, dert1, dert2, d2ert1, d2ert2;
  kekrt1=k*d*rr*ekrt1;
  k2ekrt1=k*d*rr*kekrt1;
  kekrt2=k*umd*rr*ekrt2;
  k2ekrt2=k*umd*rr*kekrt2;
  dert1=d*rr*ert1; dert2=umd*rr*ert2; d2ert1=d*rr*dert1; d2ert2=umd*rr*dert2;

  if(nd->alive1){
  nd->dp1[cl1][cl2][2][0][0]=-umgc1s2*dert1-gc1*kekrt1; 			
  nd->dp1[cl1][cl2][2][0][1]=gc1s2*dert1;				
  nd->dp1[cl1][cl2][2][0][2]=-gc1s2*dert1+gc1*kekrt1;			
  nd->dp1[cl1][cl2][2][0][3]=umgc1s2*dert1;				
  nd->dp1[cl1][cl2][2][1][0]=umgc1s2*dert1;				
  nd->dp1[cl1][cl2][2][1][1]=-gc1s2*dert1-umgc1*kekrt1;		
  nd->dp1[cl1][cl2][2][1][2]=gc1s2*dert1;				
  nd->dp1[cl1][cl2][2][1][3]=-umgc1s2*dert1+umgc1*kekrt1;		
  nd->dp1[cl1][cl2][2][2][0]=-umgc1s2*dert1+umgc1*kekrt1;		
  nd->dp1[cl1][cl2][2][2][1]=gc1s2*dert1;				
  nd->dp1[cl1][cl2][2][2][2]=-gc1s2*dert1-umgc1*kekrt1;		
  nd->dp1[cl1][cl2][2][2][3]=umgc1s2*dert1;				
  nd->dp1[cl1][cl2][2][3][0]=umgc1s2*dert1;				
  nd->dp1[cl1][cl2][2][3][1]=-gc1s2*dert1+gc1*kekrt1;			
  nd->dp1[cl1][cl2][2][3][2]=gc1s2*dert1;				
  nd->dp1[cl1][cl2][2][3][3]=-umgc1s2*dert1-gc1*kekrt1; 
  }

  if(nd->alive2){
  nd->dp2[cl1][cl2][2][0][0]=-umgc2s2*dert2-gc2*kekrt2; 			
  nd->dp2[cl1][cl2][2][0][1]=gc2s2*dert2;				
  nd->dp2[cl1][cl2][2][0][2]=-gc2s2*dert2+gc2*kekrt2;			
  nd->dp2[cl1][cl2][2][0][3]=umgc2s2*dert2;				
  nd->dp2[cl1][cl2][2][1][0]=umgc2s2*dert2;				
  nd->dp2[cl1][cl2][2][1][1]=-gc2s2*dert2-umgc2*kekrt2;		
  nd->dp2[cl1][cl2][2][1][2]=gc2s2*dert2;				
  nd->dp2[cl1][cl2][2][1][3]=-umgc2s2*dert2+umgc2*kekrt2;		
  nd->dp2[cl1][cl2][2][2][0]=-umgc2s2*dert2+umgc2*kekrt2;		
  nd->dp2[cl1][cl2][2][2][1]=gc2s2*dert2;				
  nd->dp2[cl1][cl2][2][2][2]=-gc2s2*dert2-umgc2*kekrt2;		
  nd->dp2[cl1][cl2][2][2][3]=umgc2s2*dert2;				
  nd->dp2[cl1][cl2][2][3][0]=umgc2s2*dert2;				
  nd->dp2[cl1][cl2][2][3][1]=-gc2s2*dert2+gc2*kekrt2;			
  nd->dp2[cl1][cl2][2][3][2]=gc2s2*dert2;				
  nd->dp2[cl1][cl2][2][3][3]=-umgc2s2*dert2-gc2*kekrt2; 
  }

  if(nd->alive1){
  nd->d2p1[cl1][cl2][2][0][0]=umgc1s2*d2ert1+gc1*k2ekrt1; 			
  nd->d2p1[cl1][cl2][2][0][1]=-gc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][2][0][2]=gc1s2*d2ert1-gc1*k2ekrt1;			
  nd->d2p1[cl1][cl2][2][0][3]=-umgc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][2][1][0]=-umgc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][2][1][1]=gc1s2*d2ert1+umgc1*k2ekrt1;		
  nd->d2p1[cl1][cl2][2][1][2]=-gc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][2][1][3]=umgc1s2*d2ert1-umgc1*k2ekrt1;		
  nd->d2p1[cl1][cl2][2][2][0]=umgc1s2*d2ert1-umgc1*k2ekrt1;		
  nd->d2p1[cl1][cl2][2][2][1]=-gc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][2][2][2]=gc1s2*d2ert1+umgc1*k2ekrt1;		
  nd->d2p1[cl1][cl2][2][2][3]=-umgc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][2][3][0]=-umgc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][2][3][1]=gc1s2*d2ert1-gc1*k2ekrt1;			
  nd->d2p1[cl1][cl2][2][3][2]=-gc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][2][3][3]=umgc1s2*d2ert1+gc1*k2ekrt1; 
  }

  if(nd->alive2){
  nd->d2p2[cl1][cl2][2][0][0]=umgc2s2*d2ert2+gc2*k2ekrt2; 			
  nd->d2p2[cl1][cl2][2][0][1]=-gc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][2][0][2]=gc2s2*d2ert2-gc2*k2ekrt2;			
  nd->d2p2[cl1][cl2][2][0][3]=-umgc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][2][1][0]=-umgc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][2][1][1]=gc2s2*d2ert2+umgc2*k2ekrt2;		
  nd->d2p2[cl1][cl2][2][1][2]=-gc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][2][1][3]=umgc2s2*d2ert2-umgc2*k2ekrt2;		
  nd->d2p2[cl1][cl2][2][2][0]=umgc2s2*d2ert2-umgc2*k2ekrt2;		
  nd->d2p2[cl1][cl2][2][2][1]=-gc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][2][2][2]=gc2s2*d2ert2+umgc2*k2ekrt2;		
  nd->d2p2[cl1][cl2][2][2][3]=-umgc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][2][3][0]=-umgc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][2][3][1]=gc2s2*d2ert2-gc2*k2ekrt2;			
  nd->d2p2[cl1][cl2][2][3][2]=-gc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][2][3][3]=umgc2s2*d2ert2+gc2*k2ekrt2; 
  }
}
else{
  for(i=0;i<4;i++) for(j=0;j<4;j++){
    nd->dp1[cl1][cl2][2][i][j]=0.;
    nd->dp2[cl1][cl2][2][i][j]=0.;
    nd->d2p1[cl1][cl2][2][i][j]=0.;
    nd->d2p2[cl1][cl2][2][i][j]=0.;
  }
}



if(OPTIMIZE_COMP){
  double upert1s2, upert2s2, umert1s2, umert2s2;

  upert1s2=upert1/2.; upert2s2=upert2/2.; umert1s2=umert1/2.; umert2s2=umert2/2.;

  if(nd->alive1){
  nd->dp1[cl1][cl2][3][0][0]=-upert1s2+ekrt1; 			
  nd->dp1[cl1][cl2][3][0][1]=umert1s2;				
  nd->dp1[cl1][cl2][3][0][2]=upert1s2-ekrt1;			
  nd->dp1[cl1][cl2][3][0][3]=-umert1s2;				
  nd->dp1[cl1][cl2][3][1][0]=-umert1s2;				
  nd->dp1[cl1][cl2][3][1][1]=upert1s2-ekrt1;		
  nd->dp1[cl1][cl2][3][1][2]=umert1s2;				
  nd->dp1[cl1][cl2][3][1][3]=-upert1s2+ekrt1;		
  nd->dp1[cl1][cl2][3][2][0]=-upert1s2+ekrt1;		
  nd->dp1[cl1][cl2][3][2][1]=umert1s2;				
  nd->dp1[cl1][cl2][3][2][2]=upert1s2-ekrt1;		
  nd->dp1[cl1][cl2][3][2][3]=-umert1s2;				
  nd->dp1[cl1][cl2][3][3][0]=-umert1s2;				
  nd->dp1[cl1][cl2][3][3][1]=upert1s2-ekrt1;			
  nd->dp1[cl1][cl2][3][3][2]=umert1s2;				
  nd->dp1[cl1][cl2][3][3][3]=-upert1s2+ekrt1; 
  }

  if(nd->alive2){
  nd->dp2[cl1][cl2][3][0][0]=-upert2s2+ekrt2; 			
  nd->dp2[cl1][cl2][3][0][1]=umert2s2;				
  nd->dp2[cl1][cl2][3][0][2]=upert2s2-ekrt2;			
  nd->dp2[cl1][cl2][3][0][3]=-umert2s2;				
  nd->dp2[cl1][cl2][3][1][0]=-umert2s2;				
  nd->dp2[cl1][cl2][3][1][1]=upert2s2-ekrt2;		
  nd->dp2[cl1][cl2][3][1][2]=umert2s2;				
  nd->dp2[cl1][cl2][3][1][3]=-upert2s2+ekrt2;		
  nd->dp2[cl1][cl2][3][2][0]=-upert2s2+ekrt2;		
  nd->dp2[cl1][cl2][3][2][1]=umert2s2;				
  nd->dp2[cl1][cl2][3][2][2]=upert2s2-ekrt2;		
  nd->dp2[cl1][cl2][3][2][3]=-umert2s2;				
  nd->dp2[cl1][cl2][3][3][0]=-umert2s2;				
  nd->dp2[cl1][cl2][3][3][1]=upert2s2-ekrt2;			
  nd->dp2[cl1][cl2][3][3][2]=umert2s2;				
  nd->dp2[cl1][cl2][3][3][3]=-upert2s2+ekrt2; 
  }

  for(i=0;i<4;i++) for(j=0;j<4;j++) {
    nd->d2p1[cl1][cl2][3][i][j]=0.;
    nd->d2p2[cl1][cl2][3][i][j]=0.;
  }
}
else{
  for(i=0;i<4;i++) for(j=0;j<4;j++){
    nd->dp1[cl1][cl2][3][i][j]=0.;
    nd->dp2[cl1][cl2][3][i][j]=0.;
    nd->d2p1[cl1][cl2][3][i][j]=0.;
    nd->d2p2[cl1][cl2][3][i][j]=0.;
  }
}


	/* IF UNRESOLVED NODE */

  if(nd->alive1==0){
    for(i=0;i<4;i++) for(j=0;j<4;j++)
      if(i==j) 
	nd->p1[cl1][cl2][i][j]=1.; 
      else 
	nd->p1[cl1][cl2][i][j]=0.;

    for(i=0;i<4;i++) for(j=0;j<4;j++) for(kk=0;kk<4;kk++)
      nd->dp1[cl1][cl2][i][j][kk]=0;

    for(i=0;i<4;i++) for(j=0;j<4;j++) for(kk=0;kk<4;kk++)
      nd->d2p1[cl1][cl2][i][j][kk]=0;
  }

  if(nd->alive2==0){
    for(i=0;i<4;i++) for(j=0;j<4;j++)
      if(i==j) 
	nd->p2[cl1][cl2][i][j]=1.; 
      else 
	nd->p2[cl1][cl2][i][j]=0.;

    for(i=0;i<4;i++) for(j=0;j<4;j++) for(kk=0;kk<4;kk++)
      nd->dp2[cl1][cl2][i][j][kk]=0;

    for(i=0;i<4;i++) for(j=0;j<4;j++) for(kk=0;kk<4;kk++)
      nd->d2p2[cl1][cl2][i][j][kk]=0;
  }

  tfin=time(NULL);
  probatime+=difftime(tfin, tdeb);

  compute_proba_t3_g(nd->v1, a, OPTIMIZE_LENGTH, OPTIMIZE_TITV, OPTIMIZE_ROOT, OPTIMIZE_COMP, cl1, cl2, rr);	
  compute_proba_t3_g(nd->v2, a, OPTIMIZE_LENGTH, OPTIMIZE_TITV, OPTIMIZE_ROOT, OPTIMIZE_COMP, cl1, cl2, rr);	

}


/* Returns rij, ie the ratio of effective length to length */
/* for a branch of length lg, covar rate nu, initial and final */
/* rate classes cl1 and cl2 */
/* and derivatievs with respect to x=lg.nu */

void effective_length(double lg, int cl1, int cl2, double nu, double *rij, double *drij, double* d2rij){

  double ri, rj, r;
  double a, b, x, A, B, P1, P2, P3, P4, D, N;
  double emx, em2x, em3x, em4x;
  double x2, x3, x4;
  double gm1, gm12, gm13;
  int g;

double exact, dl;

  if(cl1==cl2) g=gamma_nbcl; else g=0;
  ri=gamma_clmean[cl1];
  rj=gamma_clmean[cl2];
  r=ri+rj;

  a=2.-r;
  b=1.+(g-2.)*r/2.;

  x=lg*nu;

  if(x<=0){
    *rij=r/2.;
    *drij=*d2rij=0.;
    return;
  }

  emx=exp(-x);
  em2x=emx*emx;
  em3x=em2x*emx;
  em4x=em3x*emx;
  x2=x*x; x3=x2*x; x4=x3*x;

  D=x*(1.+(g-1.)*emx);

	/* rij */

  if(x<=0.05){   /* devlp limite */
    if(cl1!=cl2)
      *rij=r/2. + (2.-r)*x/12. - (2.-r)*x*x*x/720.;
    else
      *rij=r/2. + (r/2.-1.)*(3.*g-2.)*x*x/(6.*g);

  }
  else{  /* exact calculation */
    N=(a+b*x)*emx+x-a;
    *rij=N/D;
  }

  if(*rij<0.) *rij=0.000001;


	/* drij/dx */

  if(x<=0.05){   /* devlp limite */
    if(cl1!=cl2)
      *drij=(2.-r)/12. - (2.-r)*x*x/240.;
    else
      *drij=(1.-r/2.)*x/(3.*g);
  }
  else{  /* exact calculation */
    A=-a*(g-1.);
    B=(g-1.-b)*x2 -a*g*x +a*(g-2.);
    N=A*em2x+B*emx+a;
    *drij=N/(D*D);
  }


	/* d2rij/dx2 */

  if(x<=0.05){   /* devlp limite */
    if(g==0)
      *d2rij=(r-2.)*x/120.;
    else
      *d2rij=(1.-r/2.)/(3.*g);
  }
  else{  /* exact calculation */
    gm1=g-1.; gm12=gm1*gm1; gm13=gm12*gm1;
    P4=2.*a*(gm13)*x;
    P3=gm12*(gm1-b)*x4 - g*gm12*a*x3 + 2.*g*gm12*a*x2 - 2.*(g-4.)*gm12*a*x;
    P2=4.*g*gm1*a*x2-6.*(g-2.)*gm1*a*x;
    P1=(b-gm1)*x4+g*a*x3+2.*g*a*x2-2.*(3.*g-4.)*a*x;
    *d2rij=P4*em4x+P3*em3x+P2*em2x+P1*emx-2.*a*x;
    *d2rij/=(D*D*D*D);
  }
}



void compute_proba_t3_cov(noeud nd, double a, int OPTIMIZE_LENGTH, int OPTIMIZE_TITV, int OPTIMIZE_ROOT, int OPTIMIZE_COMP, int OPTIMIZE_COV, int cl1, int cl2){

/* Parameters:    branch lengths  :	 nd->l1 , nd->l2 				   */
/* 		  ti/tv ratio 	  :	 a 						   */
/* 		  eq. GC contents :	 nd->v1->geq+nd->v1->ceq , nd->v2->geq+nd->v2->ceq */
/* 		  root location	  : 	 nd->l1/(nd->l1+nd->l2)				   */
/* 		  ancestral GC    :      useless 					   */


  double k;					/* (1+a)/2 */
  double ert1, ert2;				/* exp(l1) , exp(l2) */
  double ekrt1, ekrt2;				/* exp(k*l1) , exp(k*l2) */
  double upert1, upert2;			/* 1+exp(l1) , 1+exp(l2) */
  double umert1, umert2;			/* 1-exp(l1) , 1-exp(l2) */
  double gc1, gc2;				
  double gc1s2, gc2s2;				/* gc1/2 , gc2/2 */
  double umgc1, umgc2, umgc1s2, umgc2s2;	/* 1-gc1 , 1-gc2 , (1-gc1)/2 , (1-gc2)/2 */
  double rr1, rr2; 				/* relative rate */
  double nul1, nul2;
  double er1, er2;				/* exp(-nu.l1), exp(-nu.l2) */

  double d, umd;
  double ri, rj, aa, bb;
  int gg;

  int i, j, kk, ll, mm, nn;


  if(nd->v1==NULL)
    return;

tdeb=time(NULL);


	/* calculs des rr (=rij) */

  effective_length(nd->l1, cl1, cl2, covar, &(nd->rr1[cl1][cl2]), &(nd->drr1[cl1][cl2]), &(nd->d2rr1[cl1][cl2]));
  effective_length(nd->l2, cl1, cl2, covar, &(nd->rr2[cl1][cl2]), &(nd->drr2[cl1][cl2]), &(nd->d2rr2[cl1][cl2]));
  rr1=nd->rr1[cl1][cl2]; rr2=nd->rr2[cl1][cl2];


	
	/* initializations */

  gc1=nd->v1->geq+nd->v1->ceq; gc1s2=gc1/2.; umgc1=1.-gc1; umgc1s2=umgc1/2;
  gc2=nd->v2->geq+nd->v2->ceq; gc2s2=gc2/2.; umgc2=1.-gc2; umgc2s2=umgc2/2;
  k=(a+1.)/2.;
  ert1=exp(-nd->l1*rr1);
  upert1=1.+ert1;
  umert1=1.-ert1;
  ekrt1=exp(-k*nd->l1*rr1);
  ert2=exp(-nd->l2*rr2);
  upert2=1.+ert2;
  umert2=1.-ert2;
  ekrt2=exp(-k*nd->l2*rr2);

  if(nd==root){
    d=nd->l1/(nd->l1+nd->l2); 
    umd=1.-d;
  }
  else d=umd=1.;


		/* probabilities */

  if(nd->alive1){  
  nd->p1[cl1][cl2][0][0]=umgc1s2*upert1+gc1*ekrt1; 			
  nd->p1[cl1][cl2][0][1]=gc1s2*umert1;				
  nd->p1[cl1][cl2][0][2]=gc1s2*upert1-gc1*ekrt1;			
  nd->p1[cl1][cl2][0][3]=umgc1s2*umert1;				
  nd->p1[cl1][cl2][1][0]=umgc1s2*umert1;				
  nd->p1[cl1][cl2][1][1]=gc1s2*upert1+umgc1*ekrt1;		
  nd->p1[cl1][cl2][1][2]=gc1s2*umert1;				
  nd->p1[cl1][cl2][1][3]=umgc1s2*upert1-umgc1*ekrt1;		
  nd->p1[cl1][cl2][2][0]=umgc1s2*upert1-umgc1*ekrt1;		
  nd->p1[cl1][cl2][2][1]=gc1s2*umert1;				
  nd->p1[cl1][cl2][2][2]=gc1s2*upert1+umgc1*ekrt1;		
  nd->p1[cl1][cl2][2][3]=umgc1s2*umert1;				
  nd->p1[cl1][cl2][3][0]=umgc1s2*umert1;				
  nd->p1[cl1][cl2][3][1]=gc1s2*upert1-gc1*ekrt1;			
  nd->p1[cl1][cl2][3][2]=gc1s2*umert1;				
  nd->p1[cl1][cl2][3][3]=umgc1s2*upert1+gc1*ekrt1; 	
  }	
	

if(nd->p1[cl1][cl2][0][1]<0.){
  printf("gc1s2=%f, umert1=%f\n", gc1s2, umert1);
  printf("ert1=%f, l1=%f, rr1=%f\n", ert1, nd->l1, rr1);
  printf("cl1=%d, cl2=%d\n", cl1, cl2);
  printf("alive1: %d, alive2: %d\n", nd->alive1, nd->alive2);
}

  if(nd->alive2){
  nd->p2[cl1][cl2][0][0]=umgc2s2*upert2+gc2*ekrt2; 			
  nd->p2[cl1][cl2][0][1]=gc2s2*umert2;				
  nd->p2[cl1][cl2][0][2]=gc2s2*upert2-gc2*ekrt2;			
  nd->p2[cl1][cl2][0][3]=umgc2s2*umert2;				
  nd->p2[cl1][cl2][1][0]=umgc2s2*umert2;				
  nd->p2[cl1][cl2][1][1]=gc2s2*upert2+umgc2*ekrt2;		
  nd->p2[cl1][cl2][1][2]=gc2s2*umert2;				
  nd->p2[cl1][cl2][1][3]=umgc2s2*upert2-umgc2*ekrt2;		
  nd->p2[cl1][cl2][2][0]=umgc2s2*upert2-umgc2*ekrt2;		
  nd->p2[cl1][cl2][2][1]=gc2s2*umert2;				
  nd->p2[cl1][cl2][2][2]=gc2s2*upert2+umgc2*ekrt2;		
  nd->p2[cl1][cl2][2][3]=umgc2s2*umert2;				
  nd->p2[cl1][cl2][3][0]=umgc2s2*umert2;				
  nd->p2[cl1][cl2][3][1]=gc2s2*upert2-gc2*ekrt2;			
  nd->p2[cl1][cl2][3][2]=gc2s2*umert2;				
  nd->p2[cl1][cl2][3][3]=umgc2s2*upert2+gc2*ekrt2;
  } 	

/*	
if(nd==root){
 fprintf(probaoutcov, "class %d\n", cl1);
 for(i=0;i<4;i++)
   for(j=0;j<4;j++)
     fprintf(probaoutcov, "p1[%d][%d]=%f\n", i, j, nd->p1[cl1][cl1][i][j]);
 for(i=0;i<4;i++)
   for(j=0;j<4;j++)
     fprintf(probaoutcov, "p2[%d][%d]=%f\n", i, j, nd->p2[cl1][cl1][i][j]);

 if(cl1==gamma_nbcl-1) exit(EXIT_FAILURE);
}
 */

		/* first and second derivatives */


if(OPTIMIZE_TITV){
  double rt1ekrt1s2, rt2ekrt2s2, rt1rt1ekrt1s4, rt2rt2ekrt2s4;

  rt1ekrt1s2=nd->l1*rr1*ekrt1/2.;
  rt2ekrt2s2=nd->l2*rr2*ekrt2/2.;
  rt1rt1ekrt1s4=rr1*rt1ekrt1s2*nd->l1/2.;
  rt2rt2ekrt2s4=rr2*rt2ekrt2s2*nd->l2/2.;

  if(nd->alive1){
  nd->dp1[cl1][cl2][0][0][0]=-gc1*rt1ekrt1s2;		
  nd->dp1[cl1][cl2][0][0][1]=0.;				
  nd->dp1[cl1][cl2][0][0][2]=gc1*rt1ekrt1s2;			
  nd->dp1[cl1][cl2][0][0][3]=0.;				
  nd->dp1[cl1][cl2][0][1][0]=0.;				
  nd->dp1[cl1][cl2][0][1][1]=-umgc1*rt1ekrt1s2;		
  nd->dp1[cl1][cl2][0][1][2]=0.;				
  nd->dp1[cl1][cl2][0][1][3]=umgc1*rt1ekrt1s2;		
  nd->dp1[cl1][cl2][0][2][0]=umgc1*rt1ekrt1s2;		
  nd->dp1[cl1][cl2][0][2][1]=0.;				
  nd->dp1[cl1][cl2][0][2][2]=-umgc1*rt1ekrt1s2;		
  nd->dp1[cl1][cl2][0][2][3]=0.;				
  nd->dp1[cl1][cl2][0][3][0]=0.;				
  nd->dp1[cl1][cl2][0][3][1]=gc1*rt1ekrt1s2;			
  nd->dp1[cl1][cl2][0][3][2]=0.;				
  nd->dp1[cl1][cl2][0][3][3]=-gc1*rt1ekrt1s2; 
  }

  if(nd->alive2){
  nd->dp2[cl1][cl2][0][0][0]=-gc2*rt2ekrt2s2; 			
  nd->dp2[cl1][cl2][0][0][1]=0.;				
  nd->dp2[cl1][cl2][0][0][2]=gc2*rt2ekrt2s2;			
  nd->dp2[cl1][cl2][0][0][3]=0.;				
  nd->dp2[cl1][cl2][0][1][0]=0.;				
  nd->dp2[cl1][cl2][0][1][1]=-umgc2*rt2ekrt2s2;		
  nd->dp2[cl1][cl2][0][1][2]=0.;				
  nd->dp2[cl1][cl2][0][1][3]=umgc2*rt2ekrt2s2;		
  nd->dp2[cl1][cl2][0][2][0]=umgc2*rt2ekrt2s2;		
  nd->dp2[cl1][cl2][0][2][1]=0.;				
  nd->dp2[cl1][cl2][0][2][2]=-umgc2*rt2ekrt2s2;		
  nd->dp2[cl1][cl2][0][2][3]=0.;				
  nd->dp2[cl1][cl2][0][3][0]=0.;				
  nd->dp2[cl1][cl2][0][3][1]=gc2*rt2ekrt2s2;			
  nd->dp2[cl1][cl2][0][3][2]=0.;				
  nd->dp2[cl1][cl2][0][3][3]=-gc2*rt2ekrt2s2; 
  }

  if(nd->alive1){
  nd->d2p1[cl1][cl2][0][0][0]=gc1*rt1rt1ekrt1s4; 			
  nd->d2p1[cl1][cl2][0][0][1]=0.;				
  nd->d2p1[cl1][cl2][0][0][2]=-gc1*rt1rt1ekrt1s4;			
  nd->d2p1[cl1][cl2][0][0][3]=0.;				
  nd->d2p1[cl1][cl2][0][1][0]=0.;				
  nd->d2p1[cl1][cl2][0][1][1]=umgc1*rt1rt1ekrt1s4;		
  nd->d2p1[cl1][cl2][0][1][2]=0.;				
  nd->d2p1[cl1][cl2][0][1][3]=-umgc1*rt1rt1ekrt1s4;		
  nd->d2p1[cl1][cl2][0][2][0]=-umgc1*rt1rt1ekrt1s4;		
  nd->d2p1[cl1][cl2][0][2][1]=0.;				
  nd->d2p1[cl1][cl2][0][2][2]=umgc1*rt1rt1ekrt1s4;		
  nd->d2p1[cl1][cl2][0][2][3]=0.;				
  nd->d2p1[cl1][cl2][0][3][0]=0.;				
  nd->d2p1[cl1][cl2][0][3][1]=-gc1*rt1rt1ekrt1s4;			
  nd->d2p1[cl1][cl2][0][3][2]=0.;				
  nd->d2p1[cl1][cl2][0][3][3]=gc1*rt1rt1ekrt1s4; 
  }

  if(nd->alive2){
  nd->d2p2[cl1][cl2][0][0][0]=gc2*rt2rt2ekrt2s4; 			
  nd->d2p2[cl1][cl2][0][0][1]=0.;				
  nd->d2p2[cl1][cl2][0][0][2]=-gc2*rt2rt2ekrt2s4;			
  nd->d2p2[cl1][cl2][0][0][3]=0.;				
  nd->d2p2[cl1][cl2][0][1][0]=0.;				
  nd->d2p2[cl1][cl2][0][1][1]=umgc2*rt2rt2ekrt2s4;		
  nd->d2p2[cl1][cl2][0][1][2]=0.;				
  nd->d2p2[cl1][cl2][0][1][3]=-umgc2*rt2rt2ekrt2s4;		
  nd->d2p2[cl1][cl2][0][2][0]=-umgc2*rt2rt2ekrt2s4;		
  nd->d2p2[cl1][cl2][0][2][1]=0.;				
  nd->d2p2[cl1][cl2][0][2][2]=umgc2*rt2rt2ekrt2s4;		
  nd->d2p2[cl1][cl2][0][2][3]=0.;				
  nd->d2p2[cl1][cl2][0][3][0]=0.;				
  nd->d2p2[cl1][cl2][0][3][1]=-gc2*rt2rt2ekrt2s4;			
  nd->d2p2[cl1][cl2][0][3][2]=0.;				
  nd->d2p2[cl1][cl2][0][3][3]=gc2*rt2rt2ekrt2s4; 
  }
}
else{
    for(i=0;i<4;i++) for(j=0;j<4;j++){
      nd->dp1[cl1][cl2][0][i][j]=0.;
      nd->dp2[cl1][cl2][0][i][j]=0.;
      nd->d2p1[cl1][cl2][0][i][j]=0.;
      nd->d2p2[cl1][cl2][0][i][j]=0.;
    }
}


if(OPTIMIZE_ROOT && nd==root){
  double l, kekrt1, k2ekrt1, kekrt2, k2ekrt2, dert1, dert2, d2ert1, d2ert2;

  l=nd->l1+nd->l2;
  kekrt1=k*l*rr1*ekrt1;
  k2ekrt1=k*l*rr1*kekrt1;
  kekrt2=-k*l*rr2*ekrt2;
  k2ekrt2=-k*l*rr2*kekrt2;
  dert1=l*rr1*ert1; dert2=-l*rr2*ert2; d2ert1=l*rr1*dert1; d2ert2=-l*rr2*dert2;

  if(nd->alive1){
  nd->dp1[cl1][cl2][1][0][0]=-umgc1s2*dert1-gc1*kekrt1; 			
  nd->dp1[cl1][cl2][1][0][1]=gc1s2*dert1;				
  nd->dp1[cl1][cl2][1][0][2]=-gc1s2*dert1+gc1*kekrt1;			
  nd->dp1[cl1][cl2][1][0][3]=umgc1s2*dert1;				
  nd->dp1[cl1][cl2][1][1][0]=umgc1s2*dert1;				
  nd->dp1[cl1][cl2][1][1][1]=-gc1s2*dert1-umgc1*kekrt1;		
  nd->dp1[cl1][cl2][1][1][2]=gc1s2*dert1;				
  nd->dp1[cl1][cl2][1][1][3]=-umgc1s2*dert1+umgc1*kekrt1;		
  nd->dp1[cl1][cl2][1][2][0]=-umgc1s2*dert1+umgc1*kekrt1;		
  nd->dp1[cl1][cl2][1][2][1]=gc1s2*dert1;				
  nd->dp1[cl1][cl2][1][2][2]=-gc1s2*dert1-umgc1*kekrt1;		
  nd->dp1[cl1][cl2][1][2][3]=umgc1s2*dert1;				
  nd->dp1[cl1][cl2][1][3][0]=umgc1s2*dert1;				
  nd->dp1[cl1][cl2][1][3][1]=-gc1s2*dert1+gc1*kekrt1;			
  nd->dp1[cl1][cl2][1][3][2]=gc1s2*dert1;				
  nd->dp1[cl1][cl2][1][3][3]=-umgc1s2*dert1-gc1*kekrt1; 
  }

  if(nd->alive2){
  nd->dp2[cl1][cl2][1][0][0]=-umgc2s2*dert2-gc2*kekrt2; 			
  nd->dp2[cl1][cl2][1][0][1]=gc2s2*dert2;				
  nd->dp2[cl1][cl2][1][0][2]=-gc2s2*dert2+gc2*kekrt2;			
  nd->dp2[cl1][cl2][1][0][3]=umgc2s2*dert2;				
  nd->dp2[cl1][cl2][1][1][0]=umgc2s2*dert2;				
  nd->dp2[cl1][cl2][1][1][1]=-gc2s2*dert2-umgc2*kekrt2;		
  nd->dp2[cl1][cl2][1][1][2]=gc2s2*dert2;				
  nd->dp2[cl1][cl2][1][1][3]=-umgc2s2*dert2+umgc2*kekrt2;		
  nd->dp2[cl1][cl2][1][2][0]=-umgc2s2*dert2+umgc2*kekrt2;		
  nd->dp2[cl1][cl2][1][2][1]=gc2s2*dert2;				
  nd->dp2[cl1][cl2][1][2][2]=-gc2s2*dert2-umgc2*kekrt2;		
  nd->dp2[cl1][cl2][1][2][3]=umgc2s2*dert2;				
  nd->dp2[cl1][cl2][1][3][0]=umgc2s2*dert2;				
  nd->dp2[cl1][cl2][1][3][1]=-gc2s2*dert2+gc2*kekrt2;			
  nd->dp2[cl1][cl2][1][3][2]=gc2s2*dert2;				
  nd->dp2[cl1][cl2][1][3][3]=-umgc2s2*dert2-gc2*kekrt2; 
  }

  if(nd->alive1){
  nd->d2p1[cl1][cl2][1][0][0]=umgc1s2*d2ert1+gc1*k2ekrt1; 			
  nd->d2p1[cl1][cl2][1][0][1]=-gc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][1][0][2]=gc1s2*d2ert1-gc1*k2ekrt1;			
  nd->d2p1[cl1][cl2][1][0][3]=-umgc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][1][1][0]=-umgc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][1][1][1]=gc1s2*d2ert1+umgc1*k2ekrt1;		
  nd->d2p1[cl1][cl2][1][1][2]=-gc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][1][1][3]=umgc1s2*d2ert1-umgc1*k2ekrt1;		
  nd->d2p1[cl1][cl2][1][2][0]=umgc1s2*d2ert1-umgc1*k2ekrt1;		
  nd->d2p1[cl1][cl2][1][2][1]=-gc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][1][2][2]=gc1s2*d2ert1+umgc1*k2ekrt1;		
  nd->d2p1[cl1][cl2][1][2][3]=-umgc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][1][3][0]=-umgc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][1][3][1]=gc1s2*d2ert1-gc1*k2ekrt1;			
  nd->d2p1[cl1][cl2][1][3][2]=-gc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][1][3][3]=umgc1s2*d2ert1+gc1*k2ekrt1; 
  }

  if(nd->alive2){
  nd->d2p2[cl1][cl2][1][0][0]=umgc2s2*d2ert2+gc2*k2ekrt2; 			
  nd->d2p2[cl1][cl2][1][0][1]=-gc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][1][0][2]=gc2s2*d2ert2-gc2*k2ekrt2;			
  nd->d2p2[cl1][cl2][1][0][3]=-umgc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][1][1][0]=-umgc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][1][1][1]=gc2s2*d2ert2+umgc2*k2ekrt2;		
  nd->d2p2[cl1][cl2][1][1][2]=-gc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][1][1][3]=umgc2s2*d2ert2-umgc2*k2ekrt2;		
  nd->d2p2[cl1][cl2][1][2][0]=umgc2s2*d2ert2-umgc2*k2ekrt2;		
  nd->d2p2[cl1][cl2][1][2][1]=-gc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][1][2][2]=gc2s2*d2ert2+umgc2*k2ekrt2;		
  nd->d2p2[cl1][cl2][1][2][3]=-umgc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][1][3][0]=-umgc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][1][3][1]=gc2s2*d2ert2-gc2*k2ekrt2;			
  nd->d2p2[cl1][cl2][1][3][2]=-gc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][1][3][3]=umgc2s2*d2ert2+gc2*k2ekrt2;   
  }
}
else{
  for(i=0;i<4;i++) for(j=0;j<4;j++){ 
    nd->dp1[cl1][cl2][1][i][j]=0.;
    nd->dp2[cl1][cl2][1][i][j]=0.;
    nd->d2p1[cl1][cl2][1][i][j]=0.;
    nd->d2p2[cl1][cl2][1][i][j]=0.;
  }
}


if(OPTIMIZE_LENGTH){
  double kekrt1, k2ekrt1, kekrt2, k2ekrt2, dert1, dert2, d2ert1, d2ert2;
  kekrt1=k*d*rr1*ekrt1;
  k2ekrt1=k*d*rr1*kekrt1;
  kekrt2=k*umd*rr2*ekrt2;
  k2ekrt2=k*umd*rr2*kekrt2;
  dert1=d*rr1*ert1; dert2=umd*rr2*ert2; d2ert1=d*rr1*dert1; d2ert2=umd*rr2*dert2;

  if(nd->alive1){
  nd->dp1[cl1][cl2][2][0][0]=-umgc1s2*dert1-gc1*kekrt1; 			
  nd->dp1[cl1][cl2][2][0][1]=gc1s2*dert1;				
  nd->dp1[cl1][cl2][2][0][2]=-gc1s2*dert1+gc1*kekrt1;			
  nd->dp1[cl1][cl2][2][0][3]=umgc1s2*dert1;				
  nd->dp1[cl1][cl2][2][1][0]=umgc1s2*dert1;				
  nd->dp1[cl1][cl2][2][1][1]=-gc1s2*dert1-umgc1*kekrt1;		
  nd->dp1[cl1][cl2][2][1][2]=gc1s2*dert1;				
  nd->dp1[cl1][cl2][2][1][3]=-umgc1s2*dert1+umgc1*kekrt1;		
  nd->dp1[cl1][cl2][2][2][0]=-umgc1s2*dert1+umgc1*kekrt1;		
  nd->dp1[cl1][cl2][2][2][1]=gc1s2*dert1;				
  nd->dp1[cl1][cl2][2][2][2]=-gc1s2*dert1-umgc1*kekrt1;		
  nd->dp1[cl1][cl2][2][2][3]=umgc1s2*dert1;				
  nd->dp1[cl1][cl2][2][3][0]=umgc1s2*dert1;				
  nd->dp1[cl1][cl2][2][3][1]=-gc1s2*dert1+gc1*kekrt1;			
  nd->dp1[cl1][cl2][2][3][2]=gc1s2*dert1;				
  nd->dp1[cl1][cl2][2][3][3]=-umgc1s2*dert1-gc1*kekrt1; 
  }

  if(nd->alive2){
  nd->dp2[cl1][cl2][2][0][0]=-umgc2s2*dert2-gc2*kekrt2; 			
  nd->dp2[cl1][cl2][2][0][1]=gc2s2*dert2;				
  nd->dp2[cl1][cl2][2][0][2]=-gc2s2*dert2+gc2*kekrt2;			
  nd->dp2[cl1][cl2][2][0][3]=umgc2s2*dert2;				
  nd->dp2[cl1][cl2][2][1][0]=umgc2s2*dert2;				
  nd->dp2[cl1][cl2][2][1][1]=-gc2s2*dert2-umgc2*kekrt2;		
  nd->dp2[cl1][cl2][2][1][2]=gc2s2*dert2;				
  nd->dp2[cl1][cl2][2][1][3]=-umgc2s2*dert2+umgc2*kekrt2;		
  nd->dp2[cl1][cl2][2][2][0]=-umgc2s2*dert2+umgc2*kekrt2;		
  nd->dp2[cl1][cl2][2][2][1]=gc2s2*dert2;				
  nd->dp2[cl1][cl2][2][2][2]=-gc2s2*dert2-umgc2*kekrt2;		
  nd->dp2[cl1][cl2][2][2][3]=umgc2s2*dert2;				
  nd->dp2[cl1][cl2][2][3][0]=umgc2s2*dert2;				
  nd->dp2[cl1][cl2][2][3][1]=-gc2s2*dert2+gc2*kekrt2;			
  nd->dp2[cl1][cl2][2][3][2]=gc2s2*dert2;				
  nd->dp2[cl1][cl2][2][3][3]=-umgc2s2*dert2-gc2*kekrt2; 
  }

  if(nd->alive1){
  nd->d2p1[cl1][cl2][2][0][0]=umgc1s2*d2ert1+gc1*k2ekrt1; 			
  nd->d2p1[cl1][cl2][2][0][1]=-gc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][2][0][2]=gc1s2*d2ert1-gc1*k2ekrt1;			
  nd->d2p1[cl1][cl2][2][0][3]=-umgc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][2][1][0]=-umgc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][2][1][1]=gc1s2*d2ert1+umgc1*k2ekrt1;		
  nd->d2p1[cl1][cl2][2][1][2]=-gc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][2][1][3]=umgc1s2*d2ert1-umgc1*k2ekrt1;		
  nd->d2p1[cl1][cl2][2][2][0]=umgc1s2*d2ert1-umgc1*k2ekrt1;		
  nd->d2p1[cl1][cl2][2][2][1]=-gc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][2][2][2]=gc1s2*d2ert1+umgc1*k2ekrt1;		
  nd->d2p1[cl1][cl2][2][2][3]=-umgc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][2][3][0]=-umgc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][2][3][1]=gc1s2*d2ert1-gc1*k2ekrt1;			
  nd->d2p1[cl1][cl2][2][3][2]=-gc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][2][3][3]=umgc1s2*d2ert1+gc1*k2ekrt1; 
  }

  if(nd->alive2){
  nd->d2p2[cl1][cl2][2][0][0]=umgc2s2*d2ert2+gc2*k2ekrt2; 			
  nd->d2p2[cl1][cl2][2][0][1]=-gc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][2][0][2]=gc2s2*d2ert2-gc2*k2ekrt2;			
  nd->d2p2[cl1][cl2][2][0][3]=-umgc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][2][1][0]=-umgc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][2][1][1]=gc2s2*d2ert2+umgc2*k2ekrt2;		
  nd->d2p2[cl1][cl2][2][1][2]=-gc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][2][1][3]=umgc2s2*d2ert2-umgc2*k2ekrt2;		
  nd->d2p2[cl1][cl2][2][2][0]=umgc2s2*d2ert2-umgc2*k2ekrt2;		
  nd->d2p2[cl1][cl2][2][2][1]=-gc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][2][2][2]=gc2s2*d2ert2+umgc2*k2ekrt2;		
  nd->d2p2[cl1][cl2][2][2][3]=-umgc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][2][3][0]=-umgc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][2][3][1]=gc2s2*d2ert2-gc2*k2ekrt2;			
  nd->d2p2[cl1][cl2][2][3][2]=-gc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][2][3][3]=umgc2s2*d2ert2+gc2*k2ekrt2; 
  }
}
else{
  for(i=0;i<4;i++) for(j=0;j<4;j++){
    nd->dp1[cl1][cl2][2][i][j]=0.;
    nd->dp2[cl1][cl2][2][i][j]=0.;
    nd->d2p1[cl1][cl2][2][i][j]=0.;
    nd->d2p2[cl1][cl2][2][i][j]=0.;
  }
}



if(OPTIMIZE_COMP){
  double upert1s2, upert2s2, umert1s2, umert2s2;

  upert1s2=upert1/2.; upert2s2=upert2/2.; umert1s2=umert1/2.; umert2s2=umert2/2.;

  if(nd->alive1){
  nd->dp1[cl1][cl2][3][0][0]=-upert1s2+ekrt1; 			
  nd->dp1[cl1][cl2][3][0][1]=umert1s2;				
  nd->dp1[cl1][cl2][3][0][2]=upert1s2-ekrt1;			
  nd->dp1[cl1][cl2][3][0][3]=-umert1s2;				
  nd->dp1[cl1][cl2][3][1][0]=-umert1s2;				
  nd->dp1[cl1][cl2][3][1][1]=upert1s2-ekrt1;		
  nd->dp1[cl1][cl2][3][1][2]=umert1s2;				
  nd->dp1[cl1][cl2][3][1][3]=-upert1s2+ekrt1;		
  nd->dp1[cl1][cl2][3][2][0]=-upert1s2+ekrt1;		
  nd->dp1[cl1][cl2][3][2][1]=umert1s2;				
  nd->dp1[cl1][cl2][3][2][2]=upert1s2-ekrt1;		
  nd->dp1[cl1][cl2][3][2][3]=-umert1s2;				
  nd->dp1[cl1][cl2][3][3][0]=-umert1s2;				
  nd->dp1[cl1][cl2][3][3][1]=upert1s2-ekrt1;			
  nd->dp1[cl1][cl2][3][3][2]=umert1s2;				
  nd->dp1[cl1][cl2][3][3][3]=-upert1s2+ekrt1; 
  }

  if(nd->alive2){
  nd->dp2[cl1][cl2][3][0][0]=-upert2s2+ekrt2; 			
  nd->dp2[cl1][cl2][3][0][1]=umert2s2;				
  nd->dp2[cl1][cl2][3][0][2]=upert2s2-ekrt2;			
  nd->dp2[cl1][cl2][3][0][3]=-umert2s2;				
  nd->dp2[cl1][cl2][3][1][0]=-umert2s2;				
  nd->dp2[cl1][cl2][3][1][1]=upert2s2-ekrt2;		
  nd->dp2[cl1][cl2][3][1][2]=umert2s2;				
  nd->dp2[cl1][cl2][3][1][3]=-upert2s2+ekrt2;		
  nd->dp2[cl1][cl2][3][2][0]=-upert2s2+ekrt2;		
  nd->dp2[cl1][cl2][3][2][1]=umert2s2;				
  nd->dp2[cl1][cl2][3][2][2]=upert2s2-ekrt2;		
  nd->dp2[cl1][cl2][3][2][3]=-umert2s2;				
  nd->dp2[cl1][cl2][3][3][0]=-umert2s2;				
  nd->dp2[cl1][cl2][3][3][1]=upert2s2-ekrt2;			
  nd->dp2[cl1][cl2][3][3][2]=umert2s2;				
  nd->dp2[cl1][cl2][3][3][3]=-upert2s2+ekrt2; 
  }

  for(i=0;i<4;i++) for(j=0;j<4;j++) {
    nd->d2p1[cl1][cl2][3][i][j]=0.;
    nd->d2p2[cl1][cl2][3][i][j]=0.;
  }
}
else{
  for(i=0;i<4;i++) for(j=0;j<4;j++){
    nd->dp1[cl1][cl2][3][i][j]=0.;
    nd->dp2[cl1][cl2][3][i][j]=0.;
    nd->d2p1[cl1][cl2][3][i][j]=0.;
    nd->d2p2[cl1][cl2][3][i][j]=0.;
  }
}


if(OPTIMIZE_COV){
  double kekrt1, k2ekrt1, kekrt2, k2ekrt2, dert1, dert2, d2ert1, d2ert2;
  kekrt1=k*d*nd->l1*ekrt1;
  k2ekrt1=k*d*nd->l1*kekrt1;
  kekrt2=k*umd*nd->l2*ekrt2;
  k2ekrt2=k*umd*nd->l2*kekrt2;
  dert1=d*nd->l1*ert1; dert2=umd*nd->l2*ert2; d2ert1=d*nd->l1*dert1; d2ert2=umd*nd->l2*dert2;

  if(nd->alive1){
  nd->dp1[cl1][cl2][4][0][0]=-umgc1s2*dert1-gc1*kekrt1; 			
  nd->dp1[cl1][cl2][4][0][1]=gc1s2*dert1;				
  nd->dp1[cl1][cl2][4][0][2]=-gc1s2*dert1+gc1*kekrt1;			
  nd->dp1[cl1][cl2][4][0][3]=umgc1s2*dert1;				
  nd->dp1[cl1][cl2][4][1][0]=umgc1s2*dert1;				
  nd->dp1[cl1][cl2][4][1][1]=-gc1s2*dert1-umgc1*kekrt1;		
  nd->dp1[cl1][cl2][4][1][2]=gc1s2*dert1;				
  nd->dp1[cl1][cl2][4][1][3]=-umgc1s2*dert1+umgc1*kekrt1;		
  nd->dp1[cl1][cl2][4][2][0]=-umgc1s2*dert1+umgc1*kekrt1;		
  nd->dp1[cl1][cl2][4][2][1]=gc1s2*dert1;				
  nd->dp1[cl1][cl2][4][2][2]=-gc1s2*dert1-umgc1*kekrt1;		
  nd->dp1[cl1][cl2][4][2][3]=umgc1s2*dert1;				
  nd->dp1[cl1][cl2][4][3][0]=umgc1s2*dert1;				
  nd->dp1[cl1][cl2][4][3][1]=-gc1s2*dert1+gc1*kekrt1;			
  nd->dp1[cl1][cl2][4][3][2]=gc1s2*dert1;				
  nd->dp1[cl1][cl2][4][3][3]=-umgc1s2*dert1-gc1*kekrt1;
  }

  if(nd->alive2){
  nd->dp2[cl1][cl2][4][0][0]=-umgc2s2*dert2-gc2*kekrt2; 			
  nd->dp2[cl1][cl2][4][0][1]=gc2s2*dert2;				
  nd->dp2[cl1][cl2][4][0][2]=-gc2s2*dert2+gc2*kekrt2;			
  nd->dp2[cl1][cl2][4][0][3]=umgc2s2*dert2;				
  nd->dp2[cl1][cl2][4][1][0]=umgc2s2*dert2;				
  nd->dp2[cl1][cl2][4][1][1]=-gc2s2*dert2-umgc2*kekrt2;		
  nd->dp2[cl1][cl2][4][1][2]=gc2s2*dert2;				
  nd->dp2[cl1][cl2][4][1][3]=-umgc2s2*dert2+umgc2*kekrt2;		
  nd->dp2[cl1][cl2][4][2][0]=-umgc2s2*dert2+umgc2*kekrt2;		
  nd->dp2[cl1][cl2][4][2][1]=gc2s2*dert2;				
  nd->dp2[cl1][cl2][4][2][2]=-gc2s2*dert2-umgc2*kekrt2;		
  nd->dp2[cl1][cl2][4][2][3]=umgc2s2*dert2;				
  nd->dp2[cl1][cl2][4][3][0]=umgc2s2*dert2;				
  nd->dp2[cl1][cl2][4][3][1]=-gc2s2*dert2+gc2*kekrt2;			
  nd->dp2[cl1][cl2][4][3][2]=gc2s2*dert2;				
  nd->dp2[cl1][cl2][4][3][3]=-umgc2s2*dert2-gc2*kekrt2;
  }

  if(nd->alive1){
  nd->d2p1[cl1][cl2][4][0][0]=umgc1s2*d2ert1+gc1*k2ekrt1; 			
  nd->d2p1[cl1][cl2][4][0][1]=-gc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][4][0][2]=gc1s2*d2ert1-gc1*k2ekrt1;			
  nd->d2p1[cl1][cl2][4][0][3]=-umgc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][4][1][0]=-umgc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][4][1][1]=gc1s2*d2ert1+umgc1*k2ekrt1;		
  nd->d2p1[cl1][cl2][4][1][2]=-gc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][4][1][3]=umgc1s2*d2ert1-umgc1*k2ekrt1;		
  nd->d2p1[cl1][cl2][4][2][0]=umgc1s2*d2ert1-umgc1*k2ekrt1;		
  nd->d2p1[cl1][cl2][4][2][1]=-gc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][4][2][2]=gc1s2*d2ert1+umgc1*k2ekrt1;		
  nd->d2p1[cl1][cl2][4][2][3]=-umgc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][4][3][0]=-umgc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][4][3][1]=gc1s2*d2ert1-gc1*k2ekrt1;			
  nd->d2p1[cl1][cl2][4][3][2]=-gc1s2*d2ert1;				
  nd->d2p1[cl1][cl2][4][3][3]=umgc1s2*d2ert1+gc1*k2ekrt1;
  }

  if(nd->alive2){
  nd->d2p2[cl1][cl2][4][0][0]=umgc2s2*d2ert2+gc2*k2ekrt2; 			
  nd->d2p2[cl1][cl2][4][0][1]=-gc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][4][0][2]=gc2s2*d2ert2-gc2*k2ekrt2;			
  nd->d2p2[cl1][cl2][4][0][3]=-umgc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][4][1][0]=-umgc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][4][1][1]=gc2s2*d2ert2+umgc2*k2ekrt2;		
  nd->d2p2[cl1][cl2][4][1][2]=-gc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][4][1][3]=umgc2s2*d2ert2-umgc2*k2ekrt2;		
  nd->d2p2[cl1][cl2][4][2][0]=umgc2s2*d2ert2-umgc2*k2ekrt2;		
  nd->d2p2[cl1][cl2][4][2][1]=-gc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][4][2][2]=gc2s2*d2ert2+umgc2*k2ekrt2;		
  nd->d2p2[cl1][cl2][4][2][3]=-umgc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][4][3][0]=-umgc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][4][3][1]=gc2s2*d2ert2-gc2*k2ekrt2;			
  nd->d2p2[cl1][cl2][4][3][2]=-gc2s2*d2ert2;				
  nd->d2p2[cl1][cl2][4][3][3]=umgc2s2*d2ert2+gc2*k2ekrt2;
  }
}
else{
  for(i=0;i<4;i++) for(j=0;j<4;j++){
    nd->dp1[cl1][cl2][4][i][j]=0.;
    nd->dp2[cl1][cl2][4][i][j]=0.;
    nd->d2p1[cl1][cl2][4][i][j]=0.;
    nd->d2p2[cl1][cl2][4][i][j]=0.;
  }
}




	/* IF UNRESOLVED NODE */

  if(nd->alive1==0){
    for(i=0;i<4;i++) for(j=0;j<4;j++)
      if(i==j) 
	nd->p1[cl1][cl2][i][j]=1.; 
      else 
	nd->p1[cl1][cl2][i][j]=0.;

    for(i=0;i<4;i++) for(j=0;j<4;j++) for(kk=0;kk<4;kk++)
      nd->dp1[cl1][cl2][i][j][kk]=0;

    for(i=0;i<4;i++) for(j=0;j<4;j++) for(kk=0;kk<4;kk++)
      nd->d2p1[cl1][cl2][i][j][kk]=0;
  }

  if(nd->alive2==0){
    for(i=0;i<4;i++) for(j=0;j<4;j++)
      if(i==j) 
	nd->p2[cl1][cl2][i][j]=1.; 
      else 
	nd->p2[cl1][cl2][i][j]=0.;

    for(i=0;i<4;i++) for(j=0;j<4;j++) for(kk=0;kk<4;kk++)
      nd->dp2[cl1][cl2][i][j][kk]=0;

    for(i=0;i<4;i++) for(j=0;j<4;j++) for(kk=0;kk<4;kk++)
      nd->d2p2[cl1][cl2][i][j][kk]=0;
  }

  tfin=time(NULL);
  probatime+=difftime(tfin, tdeb);

  compute_proba_t3_cov(nd->v1, a, OPTIMIZE_LENGTH, OPTIMIZE_TITV, OPTIMIZE_ROOT, OPTIMIZE_COMP, OPTIMIZE_COV, cl1, cl2);	
  compute_proba_t3_cov(nd->v2, a, OPTIMIZE_LENGTH, OPTIMIZE_TITV, OPTIMIZE_ROOT, OPTIMIZE_COMP, OPTIMIZE_COV, cl1, cl2);	

}




void compute_probacov(noeud nd){

  int i, j, cl1, cl2, par;
  int dl1, dl2;
  double rt1, ert1, rt2, ert2, x1, x2;
  double pide1, pdif1, pide2, pdif2;
  double dpide1, dpide2, dpdif1, dpdif2;
  double d2pide1, d2pide2, d2pdif1, d2pdif2;
  double p1, p2, dp1, dp2, d2p1, d2p2;
  double** rr1, **rr2, **drr1dx, **drr2dx, **d2rr1dx, **d2rr2dx;
  double dlrr1, d2lrr1, dlrr2, d2lrr2;


  if(nd==NULL) return;

  rt1=nd->l1; rt2=nd->l2;
  x1=covar*rt1; x2=covar*rt2;
if(x1<0.05) dl1=1; else dl1=0;
if(x2<0.05) dl2=1; else dl2=0;
  ert1=exp(-x1); ert2=exp(-x2);

  nd->pide1=pide1=(1./gamma_nbcl)+(gamma_nbcl-1.)*ert1/gamma_nbcl;
  nd->pdif1=pdif1=(1.-ert1)/gamma_nbcl;
  nd->pide2=pide2=(1./gamma_nbcl)+(gamma_nbcl-1.)*ert2/gamma_nbcl;
  nd->pdif2=pdif2=(1.-ert2)/gamma_nbcl;

  rr1=nd->rr1; rr2=nd->rr2;
  drr1dx=nd->drr1; drr2dx=nd->drr2;
  d2rr1dx=nd->d2rr1; d2rr2dx=nd->d2rr2;


	/* STORE PIJ's */

  for(cl1=0;cl1<gamma_nbcl;cl1++){
    for(cl2=0;cl2<gamma_nbcl;cl2++){
      for(i=0;i<4;i++){
        for(j=0;j<4;j++){
          p1prov[cl1][cl2][i][j]=nd->p1[cl1][cl2][i][j];
          p2prov[cl1][cl2][i][j]=nd->p2[cl1][cl2][i][j];
          for(par=0;par<5;par++){
            dp1prov[cl1][cl2][par][i][j]=nd->dp1[cl1][cl2][par][i][j];
            dp2prov[cl1][cl2][par][i][j]=nd->dp2[cl1][cl2][par][i][j];
            d2p1prov[cl1][cl2][par][i][j]=nd->d2p1[cl1][cl2][par][i][j];
            d2p2prov[cl1][cl2][par][i][j]=nd->d2p2[cl1][cl2][par][i][j];
          }
        }
      }
    }
  }


  	/* MODIFY PROBAS */

  for(cl1=0;cl1<gamma_nbcl;cl1++){
    for(cl2=0;cl2<gamma_nbcl;cl2++){
      for(i=0;i<4;i++){
        for(j=0;j<4;j++){
          p1=p1prov[cl1][cl2][i][j];
          p2=p2prov[cl1][cl2][i][j];
          if(cl1==cl2){
            nd->p1[cl1][cl2][i][j]=p1*pide1;
            nd->p2[cl1][cl2][i][j]=p2*pide2;
          }
          else{
            nd->p1[cl1][cl2][i][j]=p1*pdif1;
            nd->p2[cl1][cl2][i][j]=p2*pdif2;
          }
        }
      }
    }
  }


  	/* MODIFY DERIVATIVES */
  
  /* all param but branch length and covar: */

  for(par=0;par<4;par++){
    if(par==2) continue;
    for(cl1=0;cl1<gamma_nbcl;cl1++){
      for(cl2=0;cl2<gamma_nbcl;cl2++){
        for(i=0;i<4;i++){
          for(j=0;j<4;j++){
            dp1=dp1prov[cl1][cl2][par][i][j];
            dp2=dp2prov[cl1][cl2][par][i][j];
            d2p1=d2p1prov[cl1][cl2][par][i][j];
            d2p2=d2p2prov[cl1][cl2][par][i][j];
            if(cl1==cl2){
              nd->dp1[cl1][cl2][par][i][j]=dp1*pide1;
              nd->dp2[cl1][cl2][par][i][j]=dp2*pide2;
              nd->d2p1[cl1][cl2][par][i][j]=d2p1*pide1;
              nd->d2p2[cl1][cl2][par][i][j]=d2p2*pide2;
            }
            else{
              nd->dp1[cl1][cl2][par][i][j]=dp1*pdif1;
              nd->dp2[cl1][cl2][par][i][j]=dp2*pdif2;
              nd->d2p1[cl1][cl2][par][i][j]=d2p1*pdif1;
              nd->d2p2[cl1][cl2][par][i][j]=d2p2*pdif2;
            }
          }
        }
      }
    }
  }


  /* branch length: */

  dpide1=-covar*(gamma_nbcl-1.)*ert1/gamma_nbcl;
  dpide2=-covar*(gamma_nbcl-1.)*ert2/gamma_nbcl;
  dpdif1=covar*ert1/gamma_nbcl;
  dpdif2=covar*ert2/gamma_nbcl;

  d2pide1=-covar*dpide1;
  d2pide2=-covar*dpide2;
  d2pdif1=-covar*dpdif1;
  d2pdif2=-covar*dpdif2;

  for(cl1=0;cl1<gamma_nbcl;cl1++){
    for(cl2=0;cl2<gamma_nbcl;cl2++){
      dlrr1=rr1[cl1][cl2]+rt1*covar*drr1dx[cl1][cl2];
      dlrr2=rr2[cl1][cl2]+rt2*covar*drr2dx[cl1][cl2];
      for(i=0;i<4;i++){
        for(j=0;j<4;j++){
          p1=p1prov[cl1][cl2][i][j];
          p2=p2prov[cl1][cl2][i][j];
          dp1=dp1prov[cl1][cl2][2][i][j]/rr1[cl1][cl2];
          dp2=dp2prov[cl1][cl2][2][i][j]/rr2[cl1][cl2];
          if(cl1==cl2){
            nd->dp1[cl1][cl2][2][i][j]=p1*dpide1+(dp1*pide1*dlrr1);
            nd->dp2[cl1][cl2][2][i][j]=p2*dpide2+(dp2*pide2*dlrr2);
          }
          else{
            nd->dp1[cl1][cl2][2][i][j]=p1*dpdif1+(dp1*pdif1*dlrr1);
            nd->dp2[cl1][cl2][2][i][j]=p2*dpdif2+(dp2*pdif2*dlrr2);
          }
        }
      }
    }
  }


  for(cl1=0;cl1<gamma_nbcl;cl1++){
    for(cl2=0;cl2<gamma_nbcl;cl2++){
      dlrr1=rr1[cl1][cl2]+rt1*covar*drr1dx[cl1][cl2];
      dlrr2=rr2[cl1][cl2]+rt2*covar*drr2dx[cl1][cl2];
      d2lrr1=2.*covar*drr1dx[cl1][cl2]+rt1*covar*covar*d2rr1dx[cl1][cl2];
      d2lrr2=2.*covar*drr2dx[cl1][cl2]+rt2*covar*covar*d2rr2dx[cl1][cl2];
      for(i=0;i<4;i++){
        for(j=0;j<4;j++){
          p1=p1prov[cl1][cl2][i][j];
          p2=p2prov[cl1][cl2][i][j];
          dp1=dp1prov[cl1][cl2][2][i][j]/rr1[cl1][cl2];
          dp2=dp2prov[cl1][cl2][2][i][j]/rr2[cl1][cl2];
          d2p1=d2p1prov[cl1][cl2][2][i][j]/(rr1[cl1][cl2]*rr1[cl1][cl2]);
          d2p2=d2p2prov[cl1][cl2][2][i][j]/(rr2[cl1][cl2]*rr2[cl1][cl2]);
          if(cl1==cl2){
            nd->d2p1[cl1][cl2][2][i][j]=d2pide1*p1;
            nd->d2p1[cl1][cl2][2][i][j]+=(2.*dpide1*dlrr1+pide1*d2lrr1)*dp1;
            nd->d2p1[cl1][cl2][2][i][j]+=(pide1*dlrr1*dlrr1)*d2p1;

            nd->d2p2[cl1][cl2][2][i][j]=d2pide2*p2;
            nd->d2p2[cl1][cl2][2][i][j]+=(2.*dpide2*dlrr2+pide2*d2lrr2)*dp2;
            nd->d2p2[cl1][cl2][2][i][j]+=(pide2*dlrr2*dlrr2)*d2p2;
          }
          else{
            nd->d2p1[cl1][cl2][2][i][j]=d2pdif1*p1;
            nd->d2p1[cl1][cl2][2][i][j]+=(2.*dpdif1*dlrr1+pdif1*d2lrr1)*dp1;
            nd->d2p1[cl1][cl2][2][i][j]+=(pdif1*dlrr1*dlrr1)*d2p1;

            nd->d2p2[cl1][cl2][2][i][j]=d2pdif2*p2;
            nd->d2p2[cl1][cl2][2][i][j]+=(2.*dpdif2*dlrr2+pdif2*d2lrr2)*dp2;
            nd->d2p2[cl1][cl2][2][i][j]+=(pdif2*dlrr2*dlrr2)*d2p2;
          }
        }
      }
    }
  }


  /* covar: */

  dpide1=-rt1*(gamma_nbcl-1.)*ert1/gamma_nbcl;
  dpide2=-rt2*(gamma_nbcl-1.)*ert2/gamma_nbcl;
  dpdif1=rt1*ert1/gamma_nbcl;
  dpdif2=rt2*ert2/gamma_nbcl;

  d2pide1=-rt1*dpide1;
  d2pide2=-rt2*dpide2;
  d2pdif1=-rt1*dpdif1;
  d2pdif2=-rt2*dpdif2;

  for(cl1=0;cl1<gamma_nbcl;cl1++){
    for(cl2=0;cl2<gamma_nbcl;cl2++){
      dlrr1=rt1*rt1*drr1dx[cl1][cl2];
      dlrr2=rt2*rt2*drr2dx[cl1][cl2];
      for(i=0;i<4;i++){
        for(j=0;j<4;j++){
          p1=p1prov[cl1][cl2][i][j];
          p2=p2prov[cl1][cl2][i][j];
          dp1=dp1prov[cl1][cl2][2][i][j]/rr1[cl1][cl2];
          dp2=dp2prov[cl1][cl2][2][i][j]/rr2[cl1][cl2];
          if(cl1==cl2){
            nd->dp1[cl1][cl2][4][i][j]=p1*dpide1+(dp1*pide1*dlrr1);
            nd->dp2[cl1][cl2][4][i][j]=p2*dpide2+(dp2*pide2*dlrr2);
          }
          else{
            nd->dp1[cl1][cl2][4][i][j]=p1*dpdif1+(dp1*pdif1*dlrr1);
            nd->dp2[cl1][cl2][4][i][j]=p2*dpdif2+(dp2*pdif2*dlrr2);
          }
        }
      }
    }
  }


  for(cl1=0;cl1<gamma_nbcl;cl1++){
    for(cl2=0;cl2<gamma_nbcl;cl2++){
      dlrr1=rt1*rt1*drr1dx[cl1][cl2];
      dlrr2=rt2*rt2*drr2dx[cl1][cl2];
      d2lrr1=rt1*rt1*rt1*d2rr1dx[cl1][cl2];
      d2lrr2=rt2*rt2*rt2*d2rr2dx[cl1][cl2];
      for(i=0;i<4;i++){
        for(j=0;j<4;j++){
          p1=p1prov[cl1][cl2][i][j];
          p2=p2prov[cl1][cl2][i][j];
          dp1=dp1prov[cl1][cl2][2][i][j]/rr1[cl1][cl2];
          dp2=dp2prov[cl1][cl2][2][i][j]/rr2[cl1][cl2];
          d2p1=d2p1prov[cl1][cl2][2][i][j]/(rr1[cl1][cl2]*rr1[cl1][cl2]);
          d2p2=d2p2prov[cl1][cl2][2][i][j]/(rr2[cl1][cl2]*rr2[cl1][cl2]);
          if(cl1==cl2){
            nd->d2p1[cl1][cl2][4][i][j]=d2pide1*p1;
            nd->d2p1[cl1][cl2][4][i][j]+=(2.*dpide1*dlrr1+pide1*d2lrr1)*dp1;
            nd->d2p1[cl1][cl2][4][i][j]+=(pide1*dlrr1*dlrr1)*d2p1;

            nd->d2p2[cl1][cl2][4][i][j]=d2pide2*p2;
            nd->d2p2[cl1][cl2][4][i][j]+=(2.*dpide2*dlrr2+pide2*d2lrr2)*dp2;
            nd->d2p2[cl1][cl2][4][i][j]+=(pide2*dlrr2*dlrr2)*d2p2;
          }
          else{
            nd->d2p1[cl1][cl2][4][i][j]=d2pdif1*p1;
            nd->d2p1[cl1][cl2][4][i][j]+=(2.*dpdif1*dlrr1+pdif1*d2lrr1)*dp1;
            nd->d2p1[cl1][cl2][4][i][j]+=(pdif1*dlrr1*dlrr1)*d2p1;

            nd->d2p2[cl1][cl2][4][i][j]=d2pdif2*p2;
            nd->d2p2[cl1][cl2][4][i][j]+=(2.*dpdif2*dlrr2+pdif2*d2lrr2)*dp2;
            nd->d2p2[cl1][cl2][4][i][j]+=(pdif2*dlrr2*dlrr2)*d2p2;
          }
        }
      }
    }
  }


  compute_probacov(nd->v1);
  compute_probacov(nd->v2);
}



int compute_exact_probacov(noeud nd, double titv){

  double th, ths2, umths2, matcov;
  int dim, ri, rf, i, j, m1, m2, ret1, ret2;

  if(nd->v1==NULL) return 1;

  dim=4*gamma_nbcl;
  matcov=covar/(double)gamma_nbcl;

  /* bl1: set T92 */

  th=nd->v1->geq+nd->v1->ceq; ths2=th/2.; umths2=0.5-ths2;
  T92[0][1]=T92[1][2]=T92[2][1]=T92[3][2]=ths2;
  T92[0][2]=T92[3][1]=titv*ths2;
  T92[0][3]=T92[1][0]=T92[2][3]=T92[3][0]=umths2;
  T92[1][3]=T92[2][0]=titv*umths2;

  /* bl1: set R */

  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      for(ri=0;ri<gamma_nbcl;ri++){
        for(rf=0;rf<gamma_nbcl;rf++){
          m1=4*ri+i;
          m2=4*rf+j;
          if(m1==m2) continue;
          if(i==j) {R[m1][m2]=matcov; continue;}
          if(ri==rf) {R[m1][m2]=T92[i][j]*gamma_clmean[ri]; continue;}
          R[m1][m2]=0.;
        }
      }
    }
  }

  for(m1=0;m1<dim;m1++){
    R[m1][m1]=0.;
    for(m2=0;m2<dim;m2++){
      if(m2==m1) continue;
      R[m1][m1]-=R[m1][m2];
    }
  }


  /* bl1: compute P */

  rtop(dim, R, nd->l1, P);

  /* bl1: change to nhml proba */

  for(m1=0;m1<dim;m1++){
    for(m2=0;m2<dim;m2++){
      ri=m1/4; i=m1%4;
      rf=m2/4; j=m2%4;
      nd->p1[ri][rf][i][j]=P[m1][m2];

if(P[m1][m2]<0.){
  if(P[m1][m2]<-0.0000001) return 0;
  else nd->p1[ri][rf][i][j]=0.;
}

    }
  }

  /* bl2: set T92 */

  th=nd->v2->geq+nd->v2->ceq; ths2=th/2.; umths2=0.5-ths2;
  T92[0][1]=T92[1][2]=T92[2][1]=T92[3][2]=ths2;
  T92[0][2]=T92[3][1]=titv*ths2;
  T92[0][3]=T92[1][0]=T92[2][3]=T92[3][0]=umths2;
  T92[1][3]=T92[2][0]=titv*umths2;

  /* bl2: set R */

  for(i=0;i<4;i++){
    for(j=0;j<4;j++){
      for(ri=0;ri<gamma_nbcl;ri++){
        for(rf=0;rf<gamma_nbcl;rf++){
          m1=4*ri+i;
          m2=4*rf+j;
          if(m1==m2) continue;
          if(i==j) {R[m1][m2]=matcov; continue;}
          if(ri==rf) {R[m1][m2]=T92[i][j]*gamma_clmean[ri]; continue;}
          R[m1][m2]=0.;
        }
      }
    }
  }

  for(m1=0;m1<dim;m1++){
    R[m1][m1]=0.;
    for(m2=0;m2<dim;m2++){
      if(m2==m1) continue;
      R[m1][m1]-=R[m1][m2];
    }
  }

  /* bl2: compute P */

  rtop(dim, R, nd->l2, P);

  /* bl2: change to nhml proba */

  for(m1=0;m1<dim;m1++){
    for(m2=0;m2<dim;m2++){
      ri=m1/4; i=m1%4;
      rf=m2/4; j=m2%4;
      nd->p2[ri][rf][i][j]=P[m1][m2];

if(P[m1][m2]<0.){
  if(P[m1][m2]<-0.0000001) return 0;
  else nd->p2[ri][rf][i][j]=0.;
}

    }
  }

  ret1=compute_exact_probacov(nd->v1, titv);
  ret2=compute_exact_probacov(nd->v2, titv);

  return ret1*ret2;

}




int sumproba(proba pr, char* s, int verifpositive){
  int i, j;
  double p;
  for(i=0;i<4;i++){
    p=0.;
    for(j=0;j<4;j++){
      if(verifpositive && pr[i][j]<0.){
        printf("%s[%d][%d] negatif (%.20f)\n",s, i, j, pr[i][j]);
        return 1;
      }
      p+=pr[i][j];
    }
    if(fabs(p-1.)>0.0001 && fabs(p)>0.0001){
      printf("%s[%d] faux : %e\n", s, i, p);
      return 1;
    }
  }

  return 0;

}



double set_node_comp(noeud nd, double gc, tree* s_tree){
  double a, c, g, t;

  a=(1.-gc)*(1.+s_tree->omega2)/2.;
  t=(1.-gc)*(1.-s_tree->omega2)/2.;
  g=gc*(1.+s_tree->omega1)/2.;
  c=gc*(1.-s_tree->omega1)/2.;

  if(nd==root){
    nd->a=a; nd->c=c; nd->g=g; nd->t=t;
  }
  else{
    nd->aeq=a; nd->ceq=c; nd->geq=g; nd->teq=t;
  }

  reajust_comp_node(nd);

  
  if(nd==root){
    return nd->c+nd->g;
  }
  else{
    return nd->ceq+nd->geq;
  }

}



void reset_gc(tree* s_tree, noeud nd){

  double gc, newgc;

  if(nd==NULL) return;


  if(nd!=root){
    if(nd==root->v2)
      gc=s_tree->listgc[2*s_tree->nbseq-3];
    else
      gc=s_tree->listgc[nd->numerobr];

    newgc=set_node_comp(nd, gc, s_tree);


    if(nd==root->v2)
      s_tree->listgc[2*s_tree->nbseq-3]=newgc;
    else
      s_tree->listgc[nd->numerobr]=newgc; 
  

  }

  
  reset_gc(s_tree, nd->v1);
  reset_gc(s_tree, nd->v2);

}





void reset_tree(tree* s_tree){
  double newgcanc;

  if(s_tree->which_var_param[3]>=0 || s_tree->which_var_param[1]>=0) reset_bl(s_tree, root);
  if(s_tree->which_var_param[4]>=0) reset_gc(s_tree, root);
  if(s_tree->which_var_param[2]>=0) set_node_comp(root, s_tree->GCanc, s_tree);
  

}



void boundaries(tree* s_tree, int* bounded){

  int i, limit;

  limit=s_tree->nbalivebr;

/*
  limit=2*s_tree->nbseq-3;
*/

  for(i=0;i<s_tree->nb_var_param;i++)
    bounded[i]=0;


  if(s_tree->which_var_param[0]>=0){
    if(*(s_tree->var_param[s_tree->which_var_param[0]])<titvmin){
      *(s_tree->var_param[s_tree->which_var_param[0]])=titvmin;
      bounded[s_tree->which_var_param[0]]=1;
    }
    if(*(s_tree->var_param[s_tree->which_var_param[0]])>titvmax){
      *(s_tree->var_param[s_tree->which_var_param[0]])=titvmin;
      bounded[s_tree->which_var_param[0]]=1;
    }
  }


  if(s_tree->which_var_param[1]>=0){
    if(*(s_tree->var_param[s_tree->which_var_param[1]])<frmin){
      *(s_tree->var_param[s_tree->which_var_param[1]])=frmin;
      bounded[s_tree->which_var_param[1]]=1;
    }
    if(*(s_tree->var_param[s_tree->which_var_param[1]])>frmax){
      *(s_tree->var_param[s_tree->which_var_param[1]])=frmax;
      bounded[s_tree->which_var_param[1]]=1;
    }
  }

  if(s_tree->which_var_param[2]>=0){
    if(*(s_tree->var_param[s_tree->which_var_param[2]])<GCmin){
      *(s_tree->var_param[s_tree->which_var_param[2]])=GCmin;
      bounded[s_tree->which_var_param[2]]=1;
    }
    if(*(s_tree->var_param[s_tree->which_var_param[2]])>GCmax){
      *(s_tree->var_param[s_tree->which_var_param[2]])=GCmax;
      bounded[s_tree->which_var_param[2]]=1;
    }
  }

  if(s_tree->which_var_param[3]>=0){
    for(i=0;i<limit;i++)
      if(*(s_tree->var_param[s_tree->which_var_param[3]+i])<lmin){
        *(s_tree->var_param[s_tree->which_var_param[3]+i])=lmin;   
        bounded[s_tree->which_var_param[3]+i]=1;
      } 
  }

  if(s_tree->which_var_param[4]>=0){
    for(i=0;i<limit+1;i++){
      if(*(s_tree->var_param[s_tree->which_var_param[4]+i])<GCmin){
        *(s_tree->var_param[s_tree->which_var_param[4]+i])=GCmin; 
        bounded[s_tree->which_var_param[4]+i]=1;
      }
      if(*(s_tree->var_param[s_tree->which_var_param[4]+i])>GCmax){
        *(s_tree->var_param[s_tree->which_var_param[4]+i])=GCmax; 
        bounded[s_tree->which_var_param[4]+i]=1;
      }
    }   
  }

  if(s_tree->which_var_param[5]>=0){
    if(*(s_tree->var_param[s_tree->which_var_param[5]])<covmin){
      *(s_tree->var_param[s_tree->which_var_param[5]])=covmin;
      bounded[s_tree->which_var_param[5]]=1;
    }
    if(*(s_tree->var_param[s_tree->which_var_param[5]])>covmax){
      *(s_tree->var_param[s_tree->which_var_param[5]])=covmax;
      bounded[s_tree->which_var_param[5]]=1;
    }
  }


  if(s_tree->which_var_param[6]>=0){
    if(*(s_tree->var_param[s_tree->which_var_param[6]])<pimin){
      *(s_tree->var_param[s_tree->which_var_param[6]])=pimin;
      bounded[s_tree->which_var_param[6]]=1;
    }
    if(*(s_tree->var_param[s_tree->which_var_param[6]])>pimax){
      *(s_tree->var_param[s_tree->which_var_param[6]])=pimax;
      bounded[s_tree->which_var_param[6]]=1;
    }
  }

}




void resetdx(noeud nd, int param){
  int k, i, cl;

  if(nd->v1==NULL) return;
  
  for(cl=0;cl<gamma_nbcl;cl++){
    for(i=0;i<nd->nbpattern;i++){
      for(k=0;k<4;k++){
        nd->dx[cl][i][k]=0.;
      }
    }
  }

  resetdx(nd->v1, param);
  resetdx(nd->v2, param);

}


void resetd2x(noeud nd){
  int k, i, cl;

  if(nd->v1==NULL) return;
  
  for(cl=0;cl<gamma_nbcl;cl++){
    for(i=0;i<nd->nbpattern;i++){
      for(k=0;k<4;k++){
        nd->d2x[cl][i][k]=0.;
      }
    }
  }

  resetd2x(nd->v1);
  resetd2x(nd->v2);

}






void set_node_deduced_comp(noeud nd){

  noeud nd1, nd2;
  int k, nbcl;

  nd1=nd->v1; nd2=nd->v2;
  if(nd1==NULL) return;
  nd1->a=nd1->c=nd1->g=nd1->t=0.;
  nd2->a=nd2->c=nd2->g=nd2->t=0.;
  nbcl=gamma_nbcl;

  for(k=0;k<nbcl;k++){
    nd1->a+=(nd->a*nd->p1[k][k][0][0]+nd->c*nd->p1[k][k][1][0]+nd->g*nd->p1[k][k][2][0]+nd->t*nd->p1[k][k][3][0])/nbcl;
    nd1->c+=(nd->a*nd->p1[k][k][0][1]+nd->c*nd->p1[k][k][1][1]+nd->g*nd->p1[k][k][2][1]+nd->t*nd->p1[k][k][3][1])/nbcl;
    nd1->g+=(nd->a*nd->p1[k][k][0][2]+nd->c*nd->p1[k][k][1][2]+nd->g*nd->p1[k][k][2][2]+nd->t*nd->p1[k][k][3][2])/nbcl;
    nd1->t+=(nd->a*nd->p1[k][k][0][3]+nd->c*nd->p1[k][k][1][3]+nd->g*nd->p1[k][k][2][3]+nd->t*nd->p1[k][k][3][3])/nbcl;
  
    nd2->a+=(nd->a*nd->p2[k][k][0][0]+nd->c*nd->p2[k][k][1][0]+nd->g*nd->p2[k][k][2][0]+nd->t*nd->p2[k][k][3][0])/nbcl;
    nd2->c+=(nd->a*nd->p2[k][k][0][1]+nd->c*nd->p2[k][k][1][1]+nd->g*nd->p2[k][k][2][1]+nd->t*nd->p2[k][k][3][1])/nbcl;
    nd2->g+=(nd->a*nd->p2[k][k][0][2]+nd->c*nd->p2[k][k][1][2]+nd->g*nd->p2[k][k][2][2]+nd->t*nd->p2[k][k][3][2])/nbcl;
    nd2->t+=(nd->a*nd->p2[k][k][0][3]+nd->c*nd->p2[k][k][1][3]+nd->g*nd->p2[k][k][2][3]+nd->t*nd->p2[k][k][3][3])/nbcl;
  }

  set_node_deduced_comp(nd->v1);
  set_node_deduced_comp(nd->v2);
  
}




void set_listgc(tree* s_tree, noeud nd){

  if(nd==NULL) return;

  if(nd!=root){
    if(nd==root->v2)
      s_tree->listgc[2*s_tree->nbseq-3]=nd->ceq + nd->geq;
    else
      s_tree->listgc[nd->numerobr]=nd->ceq + nd->geq;
  }
  
  set_listgc(s_tree, nd->v1);
  set_listgc(s_tree, nd->v2);

}



void set_var_param(tree* s_tree, compute_option comp){

  int i, j, k;

  s_tree->nb_var_param=0;

  if(comp->OPTIMIZE_TITV) s_tree->nb_var_param++;
  if(comp->OPTIMIZE_ROOT) s_tree->nb_var_param++;
  if(comp->OPTIMIZE_ANC) s_tree->nb_var_param++;
  if(comp->OPTIMIZE_LENGTH) s_tree->nb_var_param+=s_tree->nbalivebr;
  if(comp->OPTIMIZE_COMP) s_tree->nb_var_param+=s_tree->nbalivebr+1;
  if(comp->OPTIMIZE_COV) s_tree->nb_var_param++;
  if(comp->OPTIMIZE_PI) s_tree->nb_var_param++;


  s_tree->var_param=(double**)check_alloc(s_tree->nb_var_param, sizeof(double*));
  s_tree->param_nature=(int*)check_alloc(s_tree->nb_var_param, sizeof(int));
  for(i=0;i<nb_param_type;i++) s_tree->which_var_param[i]=-1;
 
  i=0;
  if(comp->OPTIMIZE_TITV){
    s_tree->var_param[i]=&(s_tree->titv);
    s_tree->param_nature[i]=0;
    s_tree->which_var_param[0]=i;
    i++;
  }

  if(comp->OPTIMIZE_ROOT){
    s_tree->var_param[i]=&(s_tree->froot);
    s_tree->param_nature[i]=1;
    s_tree->which_var_param[1]=i;
    i++;
  }

  if(comp->OPTIMIZE_ANC){
    s_tree->var_param[i]=&(s_tree->GCanc);
    s_tree->param_nature[i]=2;
    s_tree->which_var_param[2]=i;
    i++;
  }

  if(comp->OPTIMIZE_LENGTH){
    k=0;
    for(j=0;j<2*s_tree->nbseq-3;j++){
      if(s_tree->alivebr[j]){
	s_tree->var_param[i+k]=&(s_tree->listbr[j]->length);
	s_tree->param_nature[i+k]=3;
	k++;
      }
    }
    s_tree->which_var_param[3]=i;
    i+=s_tree->nbalivebr;
  }

  if(comp->OPTIMIZE_COMP){
    k=0;
    for(j=0;j<2*s_tree->nbseq-3;j++){
      if(s_tree->alivebr[j]){
        s_tree->var_param[i+k]=&(s_tree->listgc[j]);
        s_tree->param_nature[i+k]=4;
	k++;
      }
    }
    s_tree->var_param[i+k]=&(s_tree->listgc[2*s_tree->nbseq-3]);
    s_tree->param_nature[i+k]=4;
    s_tree->which_var_param[4]=i;
    i+=s_tree->nbalivebr; i++;
  }

  if(comp->OPTIMIZE_COV){
    s_tree->var_param[i]=&covar;
    s_tree->param_nature[i]=5;
    s_tree->which_var_param[5]=i;
    i++;
  }

  if(comp->OPTIMIZE_PI){
    s_tree->var_param[i]=&pi;
    s_tree->param_nature[i]=6;
    s_tree->which_var_param[6]=i;
  }

}



void reset_like(tree* s_tree){

  int i, j, k, l, nbpattern;

  for(i=s_tree->nbseq;i<2*s_tree->nbseq-2;i++){
    nbpattern=s_tree->node[i]->nbpattern;
    for(j=0;j<nbpattern;j++){
      for(l=0;l<gamma_nbcl;l++){
        for(k=0;k<4;k++) {
          s_tree->node[i]->x[l][j][k]=-1.;
        }
      }
    }
  }
}



dlikedfr_node(noeud nd, int param){

  int i, j, k, son1, son2, nbpattern, cl, cl2, deb, fin;
  double s11, sd1, s22, sd2;
  like *x1, *x2;
  noeud v1, v2;

  if(nd->v1==NULL)
     return;

  v1=nd->v1; v2=nd->v2;
  nbpattern=nd->nbpattern;

  for(i=0;i<nbpattern;i++){
    son1=nd->patternson1[i];
    son2=nd->patternson2[i];
    for(cl=0;cl<gamma_nbcl;cl++){
      if(covar>0.){deb=0; fin=gamma_nbcl;} else {deb=cl; fin=cl+1;}
      for(j=0;j<4;j++){
        s11=sd1=s22=sd2=0.;
        for(cl2=deb;cl2<fin;cl2++){
          x1=v1->x[cl2]; x2=v2->x[cl2];
          for(k=0;k<4;k++){
	    s11+=nd->p1[cl][cl2][j][k]*x1[son1][k];
	    sd1+=nd->dp1[cl][cl2][1][j][k]*x1[son1][k];
	    s22+=nd->p2[cl][cl2][j][k]*x2[son2][k];
	    sd2+=nd->dp2[cl][cl2][1][j][k]*x2[son2][k];
/*if(i==0 && j==3 && cl==0) printf("k=%d: %e %e %e %e\n", k, s11, sd1, s22, sd2);*/
          }
        }
        nd->dx[cl][i][j]=s11*sd2+s22*sd1;
/*if(i==0 && j==3 && cl==0) printf("dx 003=%e\n", nd->dx[cl][i][j]);*/
      }
    }
  }
}



dlikedoth_node(noeud nd, int param){

  int i, j, k, son1, son2, nbpattern, cl, cl2, deb, fin;
  double s11, sd1, s1d, s22, sd2, s2d;
  like *x1, *x2;
  noeud v1, v2;

  if(nd->v1==NULL)
     return;

  v1=nd->v1; v2=nd->v2;

  nbpattern=nd->nbpattern;


  dlikedoth_node(v1, param);
  dlikedoth_node(v2, param);

  
  for(i=0;i<nbpattern;i++){
    son1=nd->patternson1[i];
    son2=nd->patternson2[i];
    for(cl=0;cl<gamma_nbcl;cl++){
      if(covar>0.){deb=0; fin=gamma_nbcl;} else {deb=cl; fin=cl+1;}
      for(j=0;j<4;j++){
        s11=s1d=sd1=s22=s2d=sd2=0.;
        for(cl2=deb;cl2<fin;cl2++){
          x1=v1->x[cl2]; x2=v2->x[cl2];
          for(k=0;k<4;k++){
	    s11+=nd->p1[cl][cl2][j][k]*x1[son1][k];
	    sd1+=nd->dp1[cl][cl2][0][j][k]*x1[son1][k];
	    s1d+=nd->p1[cl][cl2][j][k]*v1->dx[cl2][son1][k];
	    s22+=nd->p2[cl][cl2][j][k]*x2[son2][k];
	    sd2+=nd->dp2[cl][cl2][0][j][k]*x2[son2][k];
	    s2d+=nd->p2[cl][cl2][j][k]*v2->dx[cl2][son2][k];
          }
        }
        nd->dx[cl][i][j]=s11*(sd2+s2d)+s22*(sd1+s1d);
      }
    }
  }

}


dlikedcov_node(noeud nd, int param){

  int i, j, k, son1, son2, nbpattern, cl, cl2, deb, fin;
  double s11, sd1, s1d, s22, sd2, s2d;
  like *x1, *x2;
  noeud v1, v2;

  if(nd->v1==NULL)
     return;

  v1=nd->v1; v2=nd->v2;

  nbpattern=nd->nbpattern;

  dlikedcov_node(v1, param);
  dlikedcov_node(v2, param);

  deb=0; fin=gamma_nbcl;
  
  for(i=0;i<nbpattern;i++){
    son1=nd->patternson1[i];
    son2=nd->patternson2[i];
    for(cl=0;cl<gamma_nbcl;cl++){
      for(j=0;j<4;j++){
        s11=s1d=sd1=s22=s2d=sd2=0.;
        for(cl2=deb;cl2<fin;cl2++){
          x1=v1->x[cl2]; x2=v2->x[cl2];
          for(k=0;k<4;k++){
	    s11+=nd->p1[cl][cl2][j][k]*x1[son1][k];
	    sd1+=nd->dp1[cl][cl2][4][j][k]*x1[son1][k];
	    s1d+=nd->p1[cl][cl2][j][k]*v1->dx[cl2][son1][k];
	    s22+=nd->p2[cl][cl2][j][k]*x2[son2][k];
	    sd2+=nd->dp2[cl][cl2][4][j][k]*x2[son2][k];
	    s2d+=nd->p2[cl][cl2][j][k]*v2->dx[cl2][son2][k];
          }
        }
        nd->dx[cl][i][j]=s11*(sd2+s2d)+s22*(sd1+s1d);
      }
    }
  }

}





dlikedbl_node(noeud nd, int numbr, int param){

  int i, j, k, son1, son2, under, nbpattern, cl, cl2, deb, fin;
  double s11, sd1, s1d, s22, sd2, s2d;
  like* x1, *x2;
  noeud v1, v2;

  if(nd->v1==NULL)
     return;

  v1=nd->v1; v2=nd->v2;
  nbpattern=nd->nbpattern;

  dlikedbl_node(v1, numbr, param);
  dlikedbl_node(v2, numbr, param);

  under=nd->underlyingbr[numbr];

  if(under==0) {
  } 
  else if(under==1){
    for(i=0;i<nbpattern;i++){
      son1=nd->patternson1[i];
      son2=nd->patternson2[i];
      for(cl=0;cl<gamma_nbcl;cl++){
        if(covar>0.){deb=0; fin=gamma_nbcl;} else {deb=cl; fin=cl+1;}
        for(j=0;j<4;j++){
	  s1d=s22=0.;
          for(cl2=deb;cl2<fin;cl2++){
            x1=v1->x[cl2]; x2=v2->x[cl2];
	    for(k=0;k<4;k++){
	      s1d+=nd->p1[cl][cl2][j][k]*v1->dx[cl2][son1][k];
	      s22+=nd->p2[cl][cl2][j][k]*x2[son2][k];
            }
          }
	  nd->dx[cl][i][j]=s1d*s22;
        }
      }
    }
  }    
  else if(under==2){
    for(i=0;i<nbpattern;i++){
      son1=nd->patternson1[i];
      son2=nd->patternson2[i];
      for(cl=0;cl<gamma_nbcl;cl++){
        if(covar>0.){deb=0; fin=gamma_nbcl;} else {deb=cl; fin=cl+1;}
	for(j=0;j<4;j++){
	  s11=s2d=0.;
          for(cl2=deb;cl2<fin;cl2++){
            x1=v1->x[cl2]; x2=v2->x[cl2];
	    for(k=0;k<4;k++){
	      s11+=nd->p1[cl][cl2][j][k]*x1[son1][k];
	      s2d+=nd->p2[cl][cl2][j][k]*v2->dx[cl2][son2][k];
	    }
          }
	  nd->dx[cl][i][j]=s11*s2d;
	}
      }
    }
  }
  else if(under==10){
    for(i=0;i<nbpattern;i++){
      son1=nd->patternson1[i];
      son2=nd->patternson2[i];
      for(cl=0;cl<gamma_nbcl;cl++){
        if(covar>0.){deb=0; fin=gamma_nbcl;} else {deb=cl; fin=cl+1;}
	for(j=0;j<4;j++){
	  sd1=s22=0.;
          for(cl2=deb;cl2<fin;cl2++){
            x1=v1->x[cl2]; x2=v2->x[cl2];
	    for(k=0;k<4;k++){
	      sd1+=nd->dp1[cl][cl2][2][j][k]*x1[son1][k];
	      s22+=nd->p2[cl][cl2][j][k]*x2[son2][k];
            }
	  }
	  nd->dx[cl][i][j]=sd1*s22;
	}
      }
    }
  }
  else if(under==20){
    for(i=0;i<nbpattern;i++){
      son1=nd->patternson1[i];
      son2=nd->patternson2[i];
      for(cl=0;cl<gamma_nbcl;cl++){
        if(covar>0.){deb=0; fin=gamma_nbcl;} else {deb=cl; fin=cl+1;}
        for(j=0;j<4;j++){
	  s11=sd2=0.;
          for(cl2=deb;cl2<fin;cl2++){
            x1=v1->x[cl2]; x2=v2->x[cl2];
            for(k=0;k<4;k++){
	      s11+=nd->p1[cl][cl2][j][k]*x1[son1][k];
	      sd2+=nd->dp2[cl][cl2][2][j][k]*x2[son2][k];
	    }
          }
	  nd->dx[cl][i][j]=s11*sd2;
	}
      }
    }
  }
  else if(under==12){
    for(i=0;i<nbpattern;i++){
      son1=nd->patternson1[i];
      son2=nd->patternson2[i];
      for(cl=0;cl<gamma_nbcl;cl++){
        if(covar>0.){deb=0; fin=gamma_nbcl;} else {deb=cl; fin=cl+1;}
	for(j=0;j<4;j++){
	  s11=sd2=s22=sd1=0.;
          for(cl2=deb;cl2<fin;cl2++){
            x1=v1->x[cl2]; x2=v2->x[cl2];
	    for(k=0;k<4;k++){
	      s11+=nd->p1[cl][cl2][j][k]*x1[son1][k];
	      sd2+=nd->dp2[cl][cl2][2][j][k]*x2[son2][k];
	      s22+=nd->p2[cl][cl2][j][k]*x2[son2][k];
	      sd1+=nd->dp1[cl][cl2][2][j][k]*x1[son1][k];
	    }
          }
	  nd->dx[cl][i][j]=s11*sd2+s22*sd1;
	}
      }
    }
  }

}



dlikedcomp_node(noeud nd, int numgc, int param){

  int i, j, k, son1, son2, under, nbpattern, cl, cl2, deb, fin;
  double s11, sd1, s1d, s22, sd2, s2d, x1x2;
  like *x1, *x2;
  noeud v1, v2;

  if(nd->v1==NULL)
     return;

  v1=nd->v1; v2=nd->v2;
  nbpattern=nd->nbpattern;

  dlikedcomp_node(v1, numgc, param);
  dlikedcomp_node(v2, numgc, param);

  under=nd->underlyinggc[numgc];

  if(under==0) {
  } 
  else if(under==1){
    for(i=0;i<nbpattern;i++){
      son1=nd->patternson1[i];
      son2=nd->patternson2[i];
      for(cl=0;cl<gamma_nbcl;cl++){
        if(covar>0.){deb=0; fin=gamma_nbcl;} else {deb=cl; fin=cl+1;}
	for(j=0;j<4;j++){
	  s1d=s22=0.;
          for(cl2=deb;cl2<fin;cl2++){
            x1=v1->x[cl2]; x2=v2->x[cl2];
	    for(k=0;k<4;k++){
	      s1d+=nd->p1[cl][cl2][j][k]*v1->dx[cl2][son1][k];
	      s22+=nd->p2[cl][cl2][j][k]*x2[son2][k];
	    }
          }
	  nd->dx[cl][i][j]=s1d*s22;
	}
      }
    }
  }
  else if(under==2){
    for(i=0;i<nbpattern;i++){
      son1=nd->patternson1[i];
      son2=nd->patternson2[i];
      for(cl=0;cl<gamma_nbcl;cl++){
        if(covar>0.){deb=0; fin=gamma_nbcl;} else {deb=cl; fin=cl+1;}
	for(j=0;j<4;j++){
	  s11=s2d=0.;
          for(cl2=deb;cl2<fin;cl2++){
            x1=v1->x[cl2]; x2=v2->x[cl2];
	    for(k=0;k<4;k++){
	      s11+=nd->p1[cl][cl2][j][k]*x1[son1][k];
	      s2d+=nd->p2[cl][cl2][j][k]*v2->dx[cl2][son2][k];
	    }
          }
	  nd->dx[cl][i][j]=s11*s2d;
	}
      }
    }
  }
  else if(under==10){
    for(i=0;i<nbpattern;i++){
      son1=nd->patternson1[i];
      son2=nd->patternson2[i];
      for(cl=0;cl<gamma_nbcl;cl++){
        if(covar>0.){deb=0; fin=gamma_nbcl;} else {deb=cl; fin=cl+1;}
	for(j=0;j<4;j++){
	  sd1=s22=0.;
          for(cl2=deb;cl2<fin;cl2++){
            x1=v1->x[cl2]; x2=v2->x[cl2];
  	    for(k=0;k<4;k++){
	      sd1+=nd->dp1[cl][cl2][3][j][k]*x1[son1][k];
	      s22+=nd->p2[cl][cl2][j][k]*x2[son2][k];
	    }
          }
	  nd->dx[cl][i][j]=sd1*s22;
	}
      }
    }
  }
  else if(under==20){ 
    for(i=0;i<nbpattern;i++){
      son1=nd->patternson1[i];
      son2=nd->patternson2[i];
      for(cl=0;cl<gamma_nbcl;cl++){
        if(covar>0.){deb=0; fin=gamma_nbcl;} else {deb=cl; fin=cl+1;}
	for(j=0;j<4;j++){
	  s11=sd2=0.;
          for(cl2=deb;cl2<fin;cl2++){
            x1=v1->x[cl2]; x2=v2->x[cl2];
	    for(k=0;k<4;k++){
	      s11+=nd->p1[cl][cl2][j][k]*x1[son1][k];
	      sd2+=nd->dp2[cl][cl2][3][j][k]*x2[son2][k];
	    }
          }
	  nd->dx[cl][i][j]=s11*sd2;
	}
      }
    }
  }
}




void dlike_dgc(tree* s_tree, int param, double* dlike_s, int cross){

  int i, k, nbseq, numgc, cl;


  nbseq=s_tree->nbseq;

  i=0;
  while(s_tree->param_nature[i]!=4)
    i++;
  k=param-i;
  i=numgc=-1;
  while(i!=k && numgc<2*nbseq-3){
    numgc++;
    if(s_tree->alivebr[numgc]) i++;
  }


  if(!cross) param=0;

  dlikedcomp_node(root, numgc, param);

  for(i=0;i<s_tree->lgseq;i++){
    dlike_s[i]=0.;
    for(cl=0;cl<gamma_nbcl;cl++){
      s_tree->classdlike_s[cl][param][i]=root->dx[cl][i][0]*root->a + root->dx[cl][i][1]*root->c + 
	    root->dx[cl][i][2]*root->g + root->dx[cl][i][3]*root->t;
      dlike_s[i]+=s_tree->classdlike_s[cl][param][i];
    }
    dlike_s[i]/=gamma_nbcl;
  }
}





void dlike_dbl(tree* s_tree, int param, double* dlike_s, int cross){

  int i, k, numbr, cl;

  i=0;
  while(s_tree->param_nature[i]!=3)
    i++;
  k=param-i;
  i=numbr=-1;
  while(i!=k){
    numbr++;
    if(s_tree->alivebr[numbr]) i++;
  }
  
  if(!cross) param=0;

  dlikedbl_node(root, numbr, param);

  for(i=0;i<s_tree->lgseq;i++){
    dlike_s[i]=0.;
    for(cl=0;cl<gamma_nbcl;cl++){
      s_tree->classdlike_s[cl][param][i]=root->dx[cl][i][0]*root->a + root->dx[cl][i][1]*root->c + 
	    root->dx[cl][i][2]*root->g + root->dx[cl][i][3]*root->t;
      dlike_s[i]+=s_tree->classdlike_s[cl][param][i];
    }
    dlike_s[i]/=gamma_nbcl;
  }
}




void dlike_dalpha(tree* s_tree, int param, double* dlike_s, int cross){

  double like;
  int i, nbseq, cl;

  nbseq=s_tree->nbseq;

  dlikedoth_node(root, param);

  for(i=0;i<s_tree->lgseq;i++){
    dlike_s[i]=0.;
    for(cl=0;cl<gamma_nbcl;cl++){
      s_tree->classdlike_s[cl][param][i]=root->dx[cl][i][0]*root->a + root->dx[cl][i][1]*root->c + 
	    root->dx[cl][i][2]*root->g + root->dx[cl][i][3]*root->t;
      dlike_s[i]+=s_tree->classdlike_s[cl][param][i];
    }
    dlike_s[i]/=gamma_nbcl;
  }
}




void dlike_dcov(tree* s_tree, int param, double* dlike_s, int cross){

  double like;
  int i, nbseq, cl;

  nbseq=s_tree->nbseq;

  dlikedcov_node(root, param);

  for(i=0;i<s_tree->lgseq;i++){
    dlike_s[i]=0.;
    for(cl=0;cl<gamma_nbcl;cl++){
      s_tree->classdlike_s[cl][param][i]=root->dx[cl][i][0]*root->a + root->dx[cl][i][1]*root->c + 
	    root->dx[cl][i][2]*root->g + root->dx[cl][i][3]*root->t;
      dlike_s[i]+=s_tree->classdlike_s[cl][param][i];
    }
    dlike_s[i]/=gamma_nbcl;
  }
}





void dlike_dfr(tree* s_tree, int param, double* dlike_s, int cross){

  int i, cl;
         /*
  if(!cross) param=0;
       */

  dlikedfr_node(root, param);

  for(i=0;i<s_tree->lgseq;i++){
    dlike_s[i]=0.;
    for(cl=0;cl<gamma_nbcl;cl++){
      s_tree->classdlike_s[cl][param][i]=root->dx[cl][i][0]*root->a;
      s_tree->classdlike_s[cl][param][i]+=root->dx[cl][i][1]*root->c;
      s_tree->classdlike_s[cl][param][i]+=root->dx[cl][i][2]*root->g;
      s_tree->classdlike_s[cl][param][i]+=root->dx[cl][i][3]*root->t;
      dlike_s[i]+=s_tree->classdlike_s[cl][param][i];
    }
    dlike_s[i]/=gamma_nbcl;
  }
}



void dlike_danc(tree* s_tree, int param, double* dlike_s, int cross){

  int i, j;
  double like, upo1s2, umo1s2, upo2s2, umo2s2;
  double rxa, rxc, rxg, rxt;

  upo1s2=(1. + s_tree->omega1)/2;
  umo1s2=(1. - s_tree->omega1)/2;
  upo2s2=(1. + s_tree->omega2)/2;
  umo2s2=(1. - s_tree->omega2)/2;


  for(i=0;i<s_tree->lgseq;i++){
    rxa=rxc=rxg=rxt=0.;
    for(j=0;j<gamma_nbcl;j++){
      rxa+=root->x[j][i][0];
      rxc+=root->x[j][i][1];
      rxg+=root->x[j][i][2];
      rxt+=root->x[j][i][3];
    }
    rxa/=gamma_nbcl; rxc/=gamma_nbcl;
    rxg/=gamma_nbcl; rxt/=gamma_nbcl;

    dlike_s[i]=-rxa*upo2s2 + rxc*umo1s2 + rxg*upo1s2 - rxt*umo2s2;
  }
}


void dlike_dpi(tree* s_tree, int param, double* dlike_s, int cross){

  int i;

  for(i=0;i<s_tree->lgseq;i++)
    dlike_s[i]=s_tree->like_s[i]-s_tree->likeasrv_s[i];

}




void d2like_node(noeud nd, int param, int pn, int urank1, int urank2){

  int i, j, k, nbpattern, und1, und2, son1, son2, cas, rootnode, cl, cl2, deb, fin;
  double p1, p2, dp1a, dp1b, dp2a, dp2b, d2p1, d2p2;
  double l1, l2, dl1a, dl1b, dl2a, dl2b, d2l1, d2l2;

/* 1,2 : neighbour    a,b : param */
  double s11, s1a, s1b, sa1, sb1, sab, sba, s12, s21;
  double t11, t1a, t1b, ta1, tb1, tab, tba, t12, t21;
  noeud v1, v2;
  
  if(nd->v1==NULL) return;
  if(pn==2) return;

tdeb=time(NULL);
  
  v1=nd->v1; v2=nd->v2;
  nbpattern=nd->nbpattern;
  
  switch(pn){
    case 0: case 1: und1=-1; break;
    case 3: und1=nd->underlyingbr[urank1]; break;
    case 4: und1=nd->underlyinggc[urank1]; break;
  }
  switch(pn){
    case 0: case 1: und2=-1; break;
    case 3: und2=nd->underlyingbr[urank2]; break;
    case 4: und2=nd->underlyinggc[urank2]; break;
  }

  if(nd==root) rootnode=1; else rootnode=0;
  cas=which_cas(pn, pn, und1, und2, rootnode);
  if(cas==32) return;
  if(cas<=0) exit(EXIT_FAILURE); 

  d2like_node(v1, param, pn, urank1, urank2);
  d2like_node(v2, param, pn, urank1, urank2);

  

  if(pn>2) pn--; 

  for(i=0;i<nbpattern;i++){
    son1=nd->patternson1[i];
    son2=nd->patternson2[i];

    for(cl=0;cl<gamma_nbcl;cl++){
      if(covar>0.){deb=0; fin=gamma_nbcl;} else {deb=cl; fin=cl+1;}

      for(j=0;j<4;j++){
        s11=sa1=s1a=sb1=s1b=sab=sba=s12=s21=0.;
        t11=ta1=t1a=tb1=t1b=tab=tba=t12=t21=0.;
        for(cl2=deb;cl2<fin;cl2++){
          for(k=0;k<4;k++){
            p1=nd->p1[cl][cl2][j][k];
	    dp1a=nd->dp1[cl][cl2][pn][j][k];
	    dp1b=nd->dp1[cl][cl2][pn][j][k];
	    d2p1=nd->d2p1[cl][cl2][pn][j][k];
	    l1=v1->x[cl2][son1][k]; 
    	    dl1a=v1->dx[cl2][son1][k];
    	    dl1b=v1->dx[cl2][son1][k];
    	    if(v1->v1) d2l1=v1->d2x[cl2][son1][k]; else d2l1=0.;
            p2=nd->p2[cl][cl2][j][k];
	    dp2a=nd->dp2[cl][cl2][pn][j][k];
	    dp2b=nd->dp2[cl][cl2][pn][j][k];
	    d2p2=nd->d2p2[cl][cl2][pn][j][k];
	    l2=v2->x[cl2][son2][k]; 
    	    dl2a=v2->dx[cl2][son2][k];
    	    dl2b=v2->dx[cl2][son2][k];
    	    if(v2->v2) d2l2=v2->d2x[cl2][son2][k]; else d2l2=0.;
	
            if(paramcase[cas][11] || paramcase[cas][12] || paramcase[cas][15] || paramcase[cas][16])
              s11+=p1*l1;  
            if(paramcase[cas][3] || paramcase[cas][4]) sa1+=dp1a*l1;  
            if(paramcase[cas][9] || paramcase[cas][13]) sb1+=dp1b*l1;  
            if(paramcase[cas][7] || paramcase[cas][8]) s1a+=p1*dl1a;  
            if(paramcase[cas][10] || paramcase[cas][14]) s1b+=p1*dl1b;  
            if(paramcase[cas][2]) sab+=dp1a*dl1b;  
            if(paramcase[cas][5]) sba+=dp1b*dl1a;  
            if(paramcase[cas][6]) s12+=p1*d2l1; 
	    if(paramcase[cas][1]) s21+=d2p1*l1; 
        
            if(paramcase[cas][1] || paramcase[cas][2] || paramcase[cas][5] || paramcase[cas][6])
              t11+=p2*l2;  
            if(paramcase[cas][9] || paramcase[cas][10]) ta1+=dp2a*l2;  
            if(paramcase[cas][3] || paramcase[cas][7]) tb1+=dp2b*l2;  
            if(paramcase[cas][13] || paramcase[cas][14]) t1a+=p2*dl2a;  
            if(paramcase[cas][4] || paramcase[cas][8]) t1b+=p2*dl2b;  
            if(paramcase[cas][12]) tab+=dp2a*dl2b;  
            if(paramcase[cas][15]) tba+=dp2b*dl2a;  
            if(paramcase[cas][16]) t12+=p2*d2l2;  
	    if(paramcase[cas][11]) t21+=d2p2*l2; 
          }
        }

        nd->d2x[cl][i][j]=0;
        if(paramcase[cas][1]) nd->d2x[cl][i][j]+=s21*t11;
        if(paramcase[cas][2]) nd->d2x[cl][i][j]+=sab*t11;
        if(paramcase[cas][3]) nd->d2x[cl][i][j]+=sa1*tb1;
        if(paramcase[cas][4]) nd->d2x[cl][i][j]+=sa1*t1b;
        if(paramcase[cas][5]) nd->d2x[cl][i][j]+=sba*t11;
        if(paramcase[cas][6]) nd->d2x[cl][i][j]+=s12*t11;
        if(paramcase[cas][7]) nd->d2x[cl][i][j]+=s1a*tb1;
        if(paramcase[cas][8]) nd->d2x[cl][i][j]+=s1a*t1b;
        if(paramcase[cas][9]) nd->d2x[cl][i][j]+=sb1*ta1;
        if(paramcase[cas][10]) nd->d2x[cl][i][j]+=s1b*ta1;
        if(paramcase[cas][11]) nd->d2x[cl][i][j]+=s11*t21;
        if(paramcase[cas][12]) nd->d2x[cl][i][j]+=s11*tab;
        if(paramcase[cas][13]) nd->d2x[cl][i][j]+=sb1*t1a;
        if(paramcase[cas][14]) nd->d2x[cl][i][j]+=s1b*t1a;
        if(paramcase[cas][15]) nd->d2x[cl][i][j]+=s11*tba;
        if(paramcase[cas][16]) nd->d2x[cl][i][j]+=s11*t12;

      }
    }
  }

  tfin=time(NULL);
  d2liketime+=difftime(tfin, tdeb);

}




void d2like_fct(tree* s_tree, int param, double* d2x_s){

  int i, pn, urank1, urank2, lgseq, cl;
  double upo1s2, umo1s2, upo2s2, umo2s2;

  lgseq=s_tree->lgseq;

  pn=s_tree->param_nature[param];

  switch(pn){
    case 0: case 1: case 2: case 5: case 6: urank1=-1; break;
    case 3: case 4:
    urank1=param; i=0;
    while(s_tree->param_nature[i]!=pn){
      i++; urank1--;
    }
    break;
    default: printf("erreur nature param\n"); return;
  }

  switch(pn){
    case 0: case 1: case 2: case 5: case 6: urank2=-1; break;
    case 3: case 4:
    urank2=param; i=0;
    while(s_tree->param_nature[i]!=pn){
      i++; urank2--;
    }
    break;
    default: printf("erreur nature param\n"); return;
  }

  if(pn==2){
    upo1s2=(1. + s_tree->omega1)/2;
    umo1s2=(1. - s_tree->omega1)/2;
    upo2s2=(1. + s_tree->omega2)/2;
    umo2s2=(1. - s_tree->omega2)/2;
  }



  if(pn!=2 && pn!=6){
    resetd2x(root);

    d2like_node(root, param, pn, urank1, urank2);

    for(i=0;i<lgseq;i++){
      d2x_s[i]=0.;
      for(cl=0;cl<gamma_nbcl;cl++){
        s_tree->classd2like_s[cl][param][i]= root->a*root->d2x[cl][i][0] + root->c*root->d2x[cl][i][1] 
	         +root->g*root->d2x[cl][i][2] + root->t*root->d2x[cl][i][3];
        d2x_s[i]+=s_tree->classd2like_s[cl][param][i];
      }
      d2x_s[i]/=gamma_nbcl;
    }
    return;
  }
    

/* GCanc and pi case */


  if(pn==2 || pn==6){
    for(i=0;i<lgseq;i++){
      d2x_s[i]=0.;
    }
  }

  return;

}




/* Computes first and second derivatives of log likelihood conditional */
/* to rate class (class). */

void derivatives(tree* s_tree){

  int i, j, k, k1, k2, nbparam, lgseq;
  double like, dlike, d2like, memcovar;

  nbparam=s_tree->nb_var_param;
  lgseq=s_tree->lgseq;

  if(pi>0. && covar>0.){
    memcovar=covar;
    covar=-1.;
    for(k1=0;k1<gamma_nbcl;k1++)
      compute_proba_t3_cov(root, s_tree->titv, compute_o->OPTIMIZE_LENGTH, compute_o->OPTIMIZE_TITV, compute_o->OPTIMIZE_ROOT, compute_o->OPTIMIZE_COMP, compute_o->OPTIMIZE_COV, k1, k1);

    like_node(root);
    for(i=0;i<nbparam;i++){
      resetdx(root, 0);
      if(s_tree->param_nature[i]==0) dlike_dalpha(s_tree, i, s_tree->dlikeasrv_s[i], 0);
      else if(s_tree->param_nature[i]==1) dlike_dfr(s_tree, i, s_tree->dlikeasrv_s[i], 0);
      else if(s_tree->param_nature[i]==2) dlike_danc(s_tree, i, s_tree->dlikeasrv_s[i], 0);
      else if(s_tree->param_nature[i]==3) dlike_dbl(s_tree, i, s_tree->dlikeasrv_s[i], 0);
      else if(s_tree->param_nature[i]==4) dlike_dgc(s_tree, i, s_tree->dlikeasrv_s[i], 0);
      else if(s_tree->param_nature[i]==5) for(k=0;k<lgseq;k++) s_tree->dlikeasrv_s[i][k]=0.;
      else if(s_tree->param_nature[i]==6) for(k=0;k<lgseq;k++) s_tree->dlikeasrv_s[i][k]=0.;

      if(s_tree->param_nature[i]==5 || s_tree->param_nature[i]==6)
        for(k=0;k<lgseq;k++) s_tree->d2likeasrv_s[i][k]=0.;
      else
        d2like_fct(s_tree, i, s_tree->d2likeasrv_s[i]);
    }

    covar=memcovar;
    for(k1=0;k1<gamma_nbcl;k1++){
      for(k2=0;k2<gamma_nbcl;k2++){
        compute_proba_t3_cov(root, s_tree->titv, compute_o->OPTIMIZE_LENGTH, compute_o->OPTIMIZE_TITV, compute_o->OPTIMIZE_ROOT, compute_o->OPTIMIZE_COMP, compute_o->OPTIMIZE_COV, k1, k2);
      }
    }
    compute_probacov(root);
    like_node(root);
  }

  for(i=0;i<nbparam;i++){
    resetdx(root, 0);
    if(s_tree->param_nature[i]==0) dlike_dalpha(s_tree, i, s_tree->dlike_s[i], 0);
    else if(s_tree->param_nature[i]==1) dlike_dfr(s_tree, i, s_tree->dlike_s[i], 0);
    else if(s_tree->param_nature[i]==2) dlike_danc(s_tree, i, s_tree->dlike_s[i], 0);
    else if(s_tree->param_nature[i]==3) dlike_dbl(s_tree, i, s_tree->dlike_s[i], 0);
    else if(s_tree->param_nature[i]==4) dlike_dgc(s_tree, i, s_tree->dlike_s[i], 0);
    else if(s_tree->param_nature[i]==5) dlike_dcov(s_tree, i, s_tree->dlike_s[i], 0);
    else if(s_tree->param_nature[i]==6) dlike_dpi(s_tree, i, s_tree->dlike_s[i], 0);

    d2like_fct(s_tree, i, s_tree->d2like_s[i]);
  }

  if(pi>0. && covar>0.){
    for(i=0;i<nbparam;i++){
      s_tree->deriv[i]=s_tree->deriv2[i]=0.;
      for(j=0;j<lgseq;j++){
        if(s_tree->param_nature[i]==6){
          like=pi*s_tree->like_s[j]+(1.-pi)*s_tree->likeasrv_s[j];
          dlike=s_tree->dlike_s[i][j];
          d2like=0.;
        }
        else{
          like=pi*s_tree->like_s[j]+(1.-pi)*s_tree->likeasrv_s[j];
          dlike=pi*s_tree->dlike_s[i][j]+(1.-pi)*s_tree->dlikeasrv_s[i][j];
          d2like=pi*s_tree->d2like_s[i][j]+(1.-pi)*s_tree->d2likeasrv_s[i][j];
        }
/*
if(j==0 && i>17) printf("i=%d: %e (%e %e) %e (%e %e) %e (%e %e)\n", i, like, s_tree->like_s[j], s_tree->likeasrv_s[j], dlike, s_tree->dlike_s[i][j], s_tree->dlikeasrv_s[i][j], d2like, s_tree->d2like_s[i][j], s_tree->d2likeasrv_s[i][j]);
*/
        s_tree->deriv[i]+=s_tree->weight[j]*dlike/like;
        s_tree->deriv2[i]+=s_tree->weight[j]*((d2like/like)-(dlike*dlike/(like*like)));
      }
    }
  }
  else{
    for(i=0;i<nbparam;i++){
      s_tree->deriv[i]=s_tree->deriv2[i]=0.;
      for(k=0;k<lgseq;k++){
        like=s_tree->like_s[k];
        dlike=s_tree->dlike_s[i][k];
        d2like=s_tree->d2like_s[i][k];
/*
if(k==0) printf("i=%d: %e %e %e\n", i, like, dlike, d2like);
*/
        s_tree->deriv[i]+=s_tree->weight[k]*dlike/like;
        s_tree->deriv2[i]+=s_tree->weight[k]*((d2like/like)-(dlike*dlike/(like*like)));
      }
    }
  }

}

  

/* like_node */
/* compute likelihood for node nd */

void like_node(noeud nd){

  int i, j, k, ret, cl, cl2, deb, fin;
  double s1, s2;

  if(nd->v1==NULL){
    return;
  }

  like_node(nd->v1);
  like_node(nd->v2);

  for(i=0;i<nd->nbpattern;i++){
    for(cl=0;cl<gamma_nbcl;cl++){
      if(covar>.0){deb=0; fin=gamma_nbcl;} else {deb=cl; fin=cl+1;}
      for(j=0;j<4;j++){
        s1=s2=0.;
        for(cl2=deb;cl2<fin;cl2++){
          for(k=0;k<4;k++){
	    s1+=nd->p1[cl][cl2][j][k]*nd->v1->x[cl2][nd->patternson1[i]][k];
	    s2+=nd->p2[cl][cl2][j][k]*nd->v2->x[cl2][nd->patternson2[i]][k];
          }
        }
        nd->x[cl][i][j]=s1*s2;
        if(nd->x[cl][i][j]<0.){
          if(nd->x[cl][i][j]<ERROR_THRESHOLD){
	    printf("i:%d , j:%d , nd->x[cl][i][j]:%.20f\n", i, j, nd->x[cl][i][j]);
            exit(EXIT_FAILURE);
          }
          else nd->x[cl][i][j]=0.;
        }
      }
    }
  }
}


/* like_node_2sites */
/* Compute joint likelihood of the current two sites for node nd */

void like_node_2sites(noeud nd){

  int i, j, i1, j1, ret, cl, cl1, deb, fin;
  double s, s1, s2, t1, t2;

  if(nd->v1==NULL){
    return;
  }

  like_node_2sites(nd->v1);
  like_node_2sites(nd->v2);

  for(cl=0;cl<gamma_nbcl;cl++){
    if(covar>0.){deb=0; fin=gamma_nbcl;} else {deb=cl; fin=cl+1;}
    for(i=0;i<4;i++){
      for(j=0;j<4;j++){
        s1=s2=0.;
        for(cl1=deb;cl1<fin;cl1++){
          for(i1=0;i1<4;i1++){
            for(j1=0;j1<4;j1++){
              t1=nd->v1->x2sites[cl1][i1][j1]*nd->p1[cl][cl1][i][i1]*nd->p1[cl][cl1][j][j1];
              t2=nd->v2->x2sites[cl1][i1][j1]*nd->p2[cl][cl1][i][i1]*nd->p2[cl][cl1][j][j1];
              s1+=t1;
              s2+=t2;
            }
          }
        }
        nd->x2sites[cl][i][j]=s1*s2;
        if(nd->x2sites[cl][i][j]<0.){
          printf("i:%d , j:%d , nd->x2sites[cl][i][j]:%.20f\n", i, j, nd->x2sites[cl][i][j]);
          exit(EXIT_FAILURE);
        }
      }
    }
  }
}



/* init_x2sites */
/* Initialize x2sites at tips for sites s1 and s2 */

void init_x2sites(noeud nd, int s1, int s2){

  int j, k, l, p1, p2;
  

  if(nd->v1){ init_x2sites(nd->v1, s1, s2); init_x2sites(nd->v2, s1 ,s2); return;}

  p1=setpattern_site(nd, s1);
  p2=setpattern_site(nd, s2);
 
  for(j=0;j<gamma_nbcl;j++){
    for(k=0;k<4;k++){ /* 4 possible states at site 1 */
      for(l=0;l<4;l++){ /* 4 possible states at site 2 */
	if(k==p1 && l==p2) nd->x2sites[j][k][l]=1.;
	else nd->x2sites[j][k][l]=0.;
      }
    }
  }
  
}


/* like_2sites */
/* Compute the joint likelihood of sites s1 and s2 assuming */
/* common changes of rate class */
/* note: patterns and patternsons are not used for this calculation */

double like_2sites(tree* s_tree, int s1, int s2){
  int i, j, k;
  double like;

  init_x2sites(root, s1, s2);

  like_node_2sites(root);

  like=0.;
  for(k=0;k<gamma_nbcl;k++){
    like+=root->x2sites[k][0][0]*root->a*root->a;
    like+=root->x2sites[k][0][1]*root->a*root->c;
    like+=root->x2sites[k][0][2]*root->a*root->g;
    like+=root->x2sites[k][0][3]*root->a*root->t;
    like+=root->x2sites[k][1][0]*root->c*root->a;
    like+=root->x2sites[k][1][1]*root->c*root->c;
    like+=root->x2sites[k][1][2]*root->c*root->g;
    like+=root->x2sites[k][1][3]*root->c*root->t;
    like+=root->x2sites[k][2][0]*root->g*root->a;
    like+=root->x2sites[k][2][1]*root->g*root->c;
    like+=root->x2sites[k][2][2]*root->g*root->g;
    like+=root->x2sites[k][2][3]*root->g*root->t;
    like+=root->x2sites[k][3][0]*root->t*root->a;
    like+=root->x2sites[k][3][1]*root->t*root->c;
    like+=root->x2sites[k][3][2]*root->t*root->g;
    like+=root->x2sites[k][3][3]*root->t*root->t;
  }

  if(like<0.){
    printf("x2sites %d %d: like <0\n", s1, s2);
    return 1.;
  }

  like/=(gamma_nbcl*gamma_nbcl);
  return like;
}




double log_like_total(tree* s_tree){
  int i, j, k, k1, k2, ret;
  double llike=0., memcovar;

/* calculation under ASRV for USSRV model */

  if(covar>0. && pi>0.){
    memcovar=covar;
    covar=-1;
    for(k1=0;k1<gamma_nbcl;k1++)
        compute_proba_t3_cov(root, s_tree->titv, compute_o->OPTIMIZE_LENGTH, compute_o->OPTIMIZE_TITV, compute_o->OPTIMIZE_ROOT, compute_o->OPTIMIZE_COMP, compute_o->OPTIMIZE_COV, k1, k1);

    like_node(root);

    for(k=0;k<gamma_nbcl;k++){
      for(i=0;i<s_tree->lgseq;i++){
        s_tree->classlike_s[k][i]=root->x[k][i][0]*root->a;
        s_tree->classlike_s[k][i]+=root->x[k][i][1]*root->c;
        s_tree->classlike_s[k][i]+=root->x[k][i][2]*root->g;
        s_tree->classlike_s[k][i]+=root->x[k][i][3]*root->t;
        if(s_tree->classlike_s[k][i]<0.){
          printf("site : %d , argument=%e\n", i, s_tree->classlike_s[k][i]);
          printf("%d : %e %e %e %e\n%e %e %e %e\n", i, root->x[k][i][0], root->x[k][i][1], root->x[k][i][2], root->x[k][i][3], root->a, root->c, root->g, root->t);
          return 1.;
        }
      }
    }

    for(i=0;i<s_tree->lgseq;i++){
      s_tree->likeasrv_s[i]=0.;
      for(k=0;k<gamma_nbcl;k++)
        s_tree->likeasrv_s[i]+=s_tree->classlike_s[k][i];
       s_tree->likeasrv_s[i]/=gamma_nbcl;
    }

    covar=memcovar;
  }

/* calculation under SSRV */


  if(exact_cov && covar>0.){
    ret=compute_exact_probacov(root, s_tree->titv);
    if(ret==0) return 1.;
  }
  else{
    for(k1=0;k1<gamma_nbcl;k1++){
      for(k2=0;k2<gamma_nbcl;k2++){
        if(k1!=k2 && covar<0.) continue;
        compute_proba_t3_cov(root, s_tree->titv, compute_o->OPTIMIZE_LENGTH, compute_o->OPTIMIZE_TITV, compute_o->OPTIMIZE_ROOT, compute_o->OPTIMIZE_COMP, compute_o->OPTIMIZE_COV, k1, k2);
      }
    }
    if(covar>0.)
      compute_probacov(root);
  }


  like_node(root);


  for(k=0;k<gamma_nbcl;k++){
    for(i=0;i<s_tree->lgseq;i++){
      s_tree->classlike_s[k][i]=root->x[k][i][0]*root->a;
      s_tree->classlike_s[k][i]+=root->x[k][i][1]*root->c;
      s_tree->classlike_s[k][i]+=root->x[k][i][2]*root->g;
      s_tree->classlike_s[k][i]+=root->x[k][i][3]*root->t;
      if(s_tree->classlike_s[k][i]<0.){
        printf("site : %d , argument=%e\n", i, s_tree->classlike_s[k][i]);
        printf("%d : %e %e %e %e\n%e %e %e %e\n", i, root->x[k][i][0], root->x[k][i][1], root->x[k][i][2], root->x[k][i][3], root->a, root->c, root->g, root->t);
        return 1.;
      }
    }
  }


  for(i=0;i<s_tree->lgseq;i++){
    s_tree->like_s[i]=0.;
    for(k=0;k<gamma_nbcl;k++)
      s_tree->like_s[i]+=s_tree->classlike_s[k][i];
     s_tree->like_s[i]/=gamma_nbcl;
  }

  if(pi>0.){
    for(i=0;i<s_tree->lgseq;i++)
      llike+=log(pi*s_tree->like_s[i]+(1.-pi)*s_tree->likeasrv_s[i])*s_tree->weight[i];
  }
  else{
    for(i=0;i<s_tree->lgseq;i++)
      llike+=log(s_tree->like_s[i])*s_tree->weight[i];
  }


  return llike;
}



void set_proba_zero(noeud nd){

  int i, j, cl1, cl2, p;
  if(nd==NULL) return;

  for(cl1=0;cl1<gamma_nbcl;cl1++){
    for(cl2=0;cl2<gamma_nbcl;cl2++){
      for(i=0;i<4;i++){
        for(j=0;j<4;j++){
          nd->p1[cl1][cl2][i][j]=0.;
          nd->p2[cl1][cl2][i][j]=0.;
          for(p=0;p<5;p++){
            nd->dp1[cl1][cl2][p][i][j]=0.;
            nd->d2p1[cl1][cl2][p][i][j]=0.;
            nd->dp2[cl1][cl2][p][i][j]=0.;
            nd->d2p2[cl1][cl2][p][i][j]=0.;
          }
        }
      }
    }
  }

  set_proba_zero(nd->v1);
  set_proba_zero(nd->v2);

}


/* minus_log_likelihood */

/* Compute minus the log-likelihood of tree s_tree (global) */
/* for parameters paramvprov and Gamma distribution given by */
/* gamma_nbcl and gamma_clmean (global) */
/* WARNING : array paramprov must start at index 1. */


double minus_log_likelihood(double *paramprov){

  int i, k1, k2, ret, nbparam;
  double l, memcovar, lssrv, lasrv;


  nbparam=s_tree->nb_var_param;

/* make param array starting at 0 */

  for(i=1;i<=nbparam;i++){
    param[i-1]=paramprov[i];
  }

/* reset tree */

  for(i=0;i<nbparam;i++)
    *(s_tree->var_param[i])=param[i];
  reset_tree(s_tree);


/* compute likelihood */

  l=log_like_total(s_tree);


  if(l>0.) {
/*
    printf("parameters :\n");
    for(i=0;i<nbparam;i++) 
      printf("%.10f\n", *(s_tree->var_param[i]));
    exit(EXIT_FAILURE);
*/
  }

  return -l;
}





double true_length(noeud nd, double alpha){

  double mylg, teta, lat, lcg, glat, glcg, coeff;
  noeud v3;
  int i;

  if(nd==root) return -1.;

  coeff=0.;
  v3=nd->v3;
  mylg=nd->l3;
  teta=nd->ceq+nd->geq;
  lat=(1+alpha*teta)/2.;
  lcg=(1+alpha*(1-teta))/2.;
 
  for(i=0;i<gamma_nbcl;i++){
    glat=gamma_clmean[i]*lat;
    glcg=gamma_clmean[i]*lcg;

    coeff+= v3->a*glat + v3->c*glcg + v3->g*glcg + v3->t*glat;
  }

  return mylg*coeff/gamma_nbcl;
}



double my_length(noeud nd, double alpha){

  double truelg, coeff, teta, lat, lcg, glat, glcg;
  noeud v3;
  int i;

  if(nd==root) return -1.;

  coeff=0.;
  v3=nd->v3;
  truelg=nd->l3;
  teta=nd->ceq+nd->geq;
  lat=(1+alpha*teta)/2.;
  lcg=(1+alpha*(1-teta))/2.;

  for(i=0;i<gamma_nbcl;i++){
    glat=gamma_clmean[i]*lat;
    glcg=gamma_clmean[i]*lcg;

    coeff+= nd->aeq*glat + nd->ceq*glcg + nd->geq*glcg + nd->teq*glat;
/* note: equil base comp are used here instead of node base comp */
/* because this (approximate) calculation occurs at the beginning */
/* of the prog, i.e. when node base comp are unknown */
  }

  return truelg*gamma_nbcl/coeff;
}



void mylg_to_truelg(noeud nd, double alpha){

  double tl;
  noeud v3;

  if(nd==NULL) return;

  tl=true_length(nd, alpha);

  if(tl>=0.){
    nd->l3=tl;
    v3=nd->v3;
    if(nd==v3->v1) v3->l1=tl; else v3->l2=tl;
  }

  mylg_to_truelg(nd->v1, alpha);
  mylg_to_truelg(nd->v2, alpha);
}



void truelg_to_mylg(noeud nd, double alpha){

  double ml;
  noeud v3;

  if(nd==NULL) return;

  ml=my_length(nd, alpha);

  if(ml>=0.){
    nd->l3=ml;
    v3=nd->v3;
    if(nd==v3->v1) v3->l1=ml; else v3->l2=ml;
  }


  truelg_to_mylg(nd->v1, alpha);
  truelg_to_mylg(nd->v2, alpha);

}


void output_tree_comp(noeud nd, FILE* out, tree* s_tree){

  static int numappel=0;
  int indbl, indgc;
  float m1=-1.;
  double truelg;


  if(nd==root){
    numappel++;
    if(out) fprintf(out, ";;data set %d\n", numappel);
    else return;
    set_node_deduced_comp(root);
  }

/* trouver le numero de parametre de la branche et du GC du noeud */

  indbl=s_tree->which_var_param[3]+nd->numerobr;
  indgc=s_tree->which_var_param[4]+nd->numerobr;
  if(nd==root->v2) indgc=s_tree->which_var_param[4]+2*s_tree->nbseq-3;
  

  truelg=true_length(nd, s_tree->titv);

  if(nd->v1==NULL){
    fprintf(out, "noeud %s (depth 0): bl=%f , aeq=%f , ceq=%f , geq=%f , teq=%f, a=%f, c=%f, g=%f, t=%f, truel=%f, dlkh/dbl=%f, d2lkh/dbl2=%f, dlkh/dgc=%f, d2lkh/dgc2=%f, gc=%f\n", nd->nom, nd->l3, nd->aeq, nd->ceq, nd->geq, nd->teq, nd->a, nd->c, nd->g, nd->t, truelg, (nd==root || indbl<0)?m1:s_tree->deriv[indbl], (nd==root|| indbl<0)?m1:s_tree->deriv2[indbl], (nd==root || indgc<0)?m1:s_tree->deriv[indgc], (nd==root || indgc<0)?m1:s_tree->deriv2[indgc], nd->c+nd->g); 
    return;
  }
  
  fprintf(out, "noeud %s (v1=%s , v2=%s) (depth %d) : bl=%f , aeq=%f , ceq=%f , geq=%f , teq=%f, a=%f, c=%f, g=%f, t=%f, truel=%f, dlkh/dbl=%f, d2lkh/dbl2=%f, dlkh/dgc=%f, d2lkh/dgc2=%f, gc=%f\n", nd->nom, nd->v1->nom, nd->v2->nom, nd->depth, nd->l3, nd->aeq, nd->ceq, nd->geq, nd->teq, nd->a, nd->c, nd->g, nd->t, truelg, (nd==root || indbl<0)?m1:s_tree->deriv[indbl], (nd==root || indbl<0)?m1:s_tree->deriv2[indbl], (nd==root || indgc<0)?m1:s_tree->deriv[indgc], (nd==root || indgc<0)?m1:s_tree->deriv2[indgc], nd->c+nd->g);
  
  

  output_tree_comp(nd->v1, out, s_tree);
  output_tree_comp(nd->v2, out, s_tree);

}



double simplex(double* initpoint, double bound, double* finalvect){

  int i, j, nbparam, nbiter, min;
  double ** point, *pointlike, ran, **contraintes, ret;

  nbparam=s_tree->nb_var_param;


/* set contraintes */

  contraintes=(double**)check_alloc(nbparam+1, sizeof(double*));
  for(i=1;i<=nbparam;i++)
    contraintes[i]=(double*)check_alloc(2, sizeof(double));

  for(i=1;i<=nbparam;i++){
    contraintes[i][0]=min_value[s_tree->param_nature[i-1]];
    contraintes[i][1]=max_value[s_tree->param_nature[i-1]];
  }

/* choose N+1 initial points around initial values given as argument */
/* WARNING : arrays start from 1 */

  point=(double**)check_alloc(nbparam+2, sizeof(double*));
  for(i=0;i<=nbparam+1;i++)
    point[i]=(double*)check_alloc(nbparam+1, sizeof(double));
  pointlike=(double*)check_alloc(nbparam+2, sizeof(double));

  for(i=2;i<=nbparam+1;i++){
    for(j=1;j<=nbparam;j++){
      ran=drand48();	/*  0<ran<1  */
      ran/=5;		/*  0<ran<0.2  */
      ran+=0.9;		/* 0.9<ran<1.1 */
      point[i][j]=initpoint[j-1]*ran;
      if(point[i][j]<contraintes[j][0]) point[i][j]=contraintes[j][0];
      else if(point[i][j]>contraintes[j][1]) point[i][j]=contraintes[j][1];
    }
  }
  for(j=1;j<=nbparam;j++)
    point[1][j]=initpoint[j-1];

/* compute the likelihoods of each initial point */
  
  for(i=1;i<=nbparam+1;i++){
    pointlike[i]=minus_log_likelihood(point[i]);
  }

/* simplex algorithm */
  amoeba_cont(point, pointlike, nbparam, bound, minus_log_likelihood, &nbiter, contraintes);

/* output best point */
  min=1;
  for(i=2;i<=nbparam+1;i++){
    if(pointlike[i]<pointlike[min])
      min=i;
  }
  for(i=0;i<nbparam;i++) finalvect[i]=point[min][i+1];

  ret=-pointlike[min];

  for(i=1;i<=nbparam;i++)
    free(contraintes[i]);
  free(contraintes);
  for(i=0;i<=nbparam+1;i++)
    free(point[i]);
  free(point);
  free(pointlike);

  return ret;
}



void set_alive_0(upgmanoeud nd){

  if(nd->v1==NULL){
    nd->alive=0;
    return;
  }

  set_alive_0(nd->v1);
  set_alive_0(nd->v2);
}


void set_tiplg_min(upgmanoeud nd, double tiplgmin){

  if(nd->v1==NULL){
    nd->l3-=tiplgmin;
    nd->v3->l1=nd->v3->l2=nd->l3;
    return;
  }

  set_tiplg_min(nd->v1, tiplgmin);
  set_tiplg_min(nd->v2, tiplgmin);
}



char* labelled_upgma(int nb, char** nom, double** dist, double *label){

  int i, j, j1, j2, k1, k2, iter, nbn, nbm1, ndnum;
  double dmin, d1, d2, tiplgmin;
  char* ctree;
  upgmanoeud* stree;


  /* alloc nodes */

  nbn=2*nb-1;
  stree=check_alloc(nbn, sizeof(upgmanoeud));

  /* set tip nodes */

  for(i=0;i<nb;i++){
    stree[i]=check_alloc(1, sizeof(struct upgmanoeud));
    stree[i]->v1=stree[i]->v2=NULL;
    stree[i]->l1=stree[i]->l2=-1;
    stree[i]->num=i;
    stree[i]->alive=1;
    stree[i]->depth=0.;
    stree[i]->nom=nom[i];
    stree[i]->b3=label[i];
  }

  /* main loop */

  nbm1=nb-1;
  ndnum=nb;
  for(iter=0;iter<nbm1;iter++){

	/* find min dist */
    dmin=DBL_MAX;
    for(i=0;i<nb;i++){
      if(stree[i]->alive==0) continue;
      for(j=i+1;j<nb;j++){
        if(stree[j]->alive==0) continue;
        if(dist[i][j]<dmin){
          dmin=dist[i][j];
          j1=i; j2=j;
        }
      }
    }

	/* find nodes to cluster */
    for(i=0;i<ndnum;i++)
      if(stree[i]->num==j1){k1=i; break;}
    for(i=0;i<ndnum;i++)
      if(stree[i]->num==j2){k2=i; break;}

	/* cluster */
    stree[ndnum]=check_alloc(1, sizeof(struct upgmanoeud));
    stree[ndnum]->v1=stree[k1];
    stree[ndnum]->v2=stree[k2];
    stree[k1]->v3=stree[ndnum];
    stree[k2]->v3=stree[ndnum];
    stree[ndnum]->l1=stree[k1]->l3=dmin/2.-stree[k1]->depth;
    stree[ndnum]->l2=stree[k2]->l3=dmin/2.-stree[k2]->depth;
    if(iter==0) tiplgmin=stree[ndnum]->l1;
    stree[ndnum]->depth=dmin/2.;
    stree[ndnum]->num=stree[k1]->num;
    stree[ndnum]->b3=(stree[k1]->b3+stree[k2]->b3)/2.;
    stree[k1]->num=-stree[k1]->num-1;
    stree[k2]->num=-stree[k2]->num-1;
    set_alive_0(stree[k2]);
    ndnum++;

	/* update distances */

    for(i=0;i<nb;i++){
      if(i<j1) d1=dist[i][j1];
      else d1=dist[j1][i];
      if(i<j2) d2=dist[i][j2];
      else d2=dist[j2][i];
      if(i<j1) dist[i][j1]=(d1+d2)/2.;
      if(i>j1) dist[j1][i]=(d1+d2)/2.;
    }

  }

  /* struct to string */

  ctree=check_alloc(nb*100, sizeof(char));
  upgmaroot=stree[ndnum-1];
  upgmaroot->v3=NULL;
  set_tiplg_min(upgmaroot, tiplgmin);
  upgma_stoc(stree, 1, nb, ctree, 0);

  return ctree;

}




/* output_covar_pattern */
/* Compute and print into file outcovfile : */
/* - the "variation index" of each site (contribution to the total likelihood */
/*   of scenarii involving at least one change of rate class ) */
/* - the "covariation index" of each pair of site (likelihood conditional */
/*   on simultaneous chanegs of rate class divided by unconstrained likelihood) */

void output_covar_pattern(tree* s_tree, double lcovar, double* bestvect, int optcov){

  int i, j, k, ii, jj, lgseq, truelgseq, *pattnumber;
  int maxsitepername=5;
  double l, lnoch, lcov, *varind, **covarind;
  char **name, **longname, **finlongname, *prov, *covartree;
  FILE* outcov;

  printf("computing variation and covariation indices...\n");
  
  outcov=fopen("covarfile", "w");
  if(outcov==NULL) return;
  lgseq=s_tree->lgseq;
  truelgseq=s_tree->weightsum;

  varind=check_alloc(lgseq, sizeof(double));
  covarind=check_alloc(lgseq, sizeof(double*));
  for(k=0;k<lgseq;k++) covarind[k]=check_alloc(lgseq, sizeof(double));
  pattnumber=check_alloc(lgseq, sizeof(int));
  longname=check_alloc(lgseq, sizeof(char*));
  finlongname=check_alloc(lgseq, sizeof(char*));
  name=check_alloc(lgseq, sizeof(char*));
  for(i=0;i<lgseq;i++) name[i]=check_alloc(maxsitepername*(log10(lgseq)+2), sizeof(char));


  minus_log_likelihood(bestvect-1);

  /* store current likelihoods */

  for(i=0;i<lgseq;i++){
    varind[i]=s_tree->like_s[i];
    for(j=i;j<lgseq;j++) covarind[i][j]=s_tree->like_s[i]*s_tree->like_s[j];
  }

  /* variation indices */

  covar=-1.;
  if(optcov) bestvect[s_tree->which_var_param[5]]=-1.;
  log_like_total(s_tree);

  for(i=0;i<lgseq;i++){
    lnoch=s_tree->like_s[i];
    varind[i]=1.-(lnoch/varind[i]);
  }

  covar=lcovar;
  if(optcov) bestvect[s_tree->nb_var_param]=lcovar;

  /* covariation indices */

  for(i=0;i<lgseq;i++){
    for(j=i;j<lgseq;j++){
      lcov=like_2sites(s_tree, i, j);
      covarind[i][j]=lcov/covarind[i][j];
    }
  }


  /* UPGMA tree for covariation indices */

	/* set names of site groups */
  ii=0;
  for(i=0;i<truelgseq;i++){
    j=s_tree->sitetopatt[i];
    k=pattnumber[j];
    pattnumber[j]++;
    if(k<maxsitepername){
      prov=name[j];
      while(*prov) prov++;
      if(prov==name[j]) sprintf(prov, "%d", i+1);
      else sprintf(prov, "_%d", i+1);
    }
    if(k==maxsitepername){
      longname[ii]=check_alloc((log10(lgseq)+1)*lgseq, sizeof(char));
      sprintf(longname[ii], "set%d: %s_%d", ii+1, name[j], i+1);
      sprintf(name[j], "set%d", ii+1);
      finlongname[ii]=longname[ii];
      while(*(finlongname[ii])) finlongname[ii]++;
      ii++;
    }
    if(k>maxsitepername){
      prov=name[j];
      sscanf(prov, "set%d", &jj);
      jj--;
      sprintf(finlongname[jj], "_%d", i+1);
      while(*(finlongname[jj])) finlongname[jj]++;
    }
  }

	/* compute distances */

  for(i=0;i<lgseq;i++)
    for(j=i+1;j<lgseq;j++)
       covarind[i][j]=1.-covarind[i][j];


	/* build tree */

  covartree=labelled_upgma(lgseq, name, covarind, varind);



	/* output tree */

  fprintf(outcov, "\nCovariation gaph:\n");
  fprintf(outcov, "%s\n", covartree);
  for(i=0;i<ii;i++){
    fprintf(outcov, "%s\n", longname[i]);
  }

  /* output variation indices */

  fprintf(outcov, "\nVariation index for each site:\n\n");
/*
  ii=0;
  while(1){
    for(j=0;j<10;j++){
      for(k=0;k<5;k++){
        i=ii*50+j*5+k;
        if(i>=truelgseq) break;
        fprintf(outcov, "%d:%.4e\t", i+1, varind[s_tree->sitetopatt[i]]);
      }
      fprintf(outcov, "\n");
      if(i>=truelgseq) break;
    }
    fprintf(outcov, "\n");
    if(i>=truelgseq) break;
    ii++;
  }
*/

  ii=0;
  while(1){
    for(j=0;j<10;j++){
      for(k=0;k<5;k++){
        i=ii*50+j*5+k;
        if(i>=truelgseq) break;
        fprintf(outcov, "%d\t%.4e\n", i+1, varind[s_tree->sitetopatt[i]]);
      }
      if(i>=truelgseq) break;
    }
    if(i>=truelgseq) break;
    ii++;
  }

  fclose(outcov);
}



/* my_eigen */
/* Diagonalize square matrix A (dim=n) */
/* job=0 -> compute eigen values only */
/* job=1 -> compute eigen values and eigen vectors */
/* eigen values returned in rr (real part) and ri (complex part) */
/* right eigen vectors returned as vr and vi */
/* A = v.diag(r).v-1 , where diag(r) is the diagonal */
/* matrix of eigenvalues, r is the matrix of right */
/* eigen vectors, and r-1 is the inverse of r */
/* rr, ri vr and vi must be passed allocated */

int my_eigen(int job, double** A, int n, double *rr, double *ri , double** vr, double** vi){

  int i,j, ret;
  double* newA;
  double* newvr, *newvi, *w;

  newA=check_alloc(n*n, sizeof(double));
  newvr=check_alloc(n*n, sizeof(double));
  newvi=check_alloc(n*n, sizeof(double));
  w=check_alloc(2*n, sizeof(double));

  for(i=0;i<n;i++)
    for(j=0;j<n;j++)
      newA[pos(i,j,n)]=A[i][j];

  ret=eigen(job, newA, n, rr, ri, newvr, newvi, w);

  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      vr[i][j]=newvr[pos(i,j,n)];
      vi[i][j]=newvi[pos(i,j,n)];
    }
  }

  free(w);
  free(newvr);
  free(newvi);
  free(newA);

  return ret;

}


/* rtop */
/* Compute transition probability matrix P from rate matrix R */
/* and branch length (time) t. */
/* lines (first dimension) of R should sum to zero:
for(i=0;i<dim;i++){
  R[i][i]=0.;
  for(j=0;j<dim;j++){
    if(j==i) continue;
    R[i][i]-=R[i][j];
  }
}
*/	

int rtop(int n, double** R, double t, double** P){

  double **vr, **vi, *rr, *ri, **invvr, **ediago, **provprod;
  int i, j, ret;

  vr=check_alloc(n, sizeof(double*));
  for(i=0;i<n;i++) vr[i]=check_alloc(n, sizeof(double));
  invvr=check_alloc(n, sizeof(double*));
  for(i=0;i<n;i++) invvr[i]=check_alloc(n, sizeof(double));
  ediago=check_alloc(n, sizeof(double*));
  for(i=0;i<n;i++) ediago[i]=check_alloc(n, sizeof(double));
  provprod=check_alloc(n, sizeof(double*));
  for(i=0;i<n;i++) provprod[i]=check_alloc(n, sizeof(double));
  vi=check_alloc(n, sizeof(double*));
  for(i=0;i<n;i++) vi[i]=check_alloc(n, sizeof(double));
  rr=check_alloc(n, sizeof(double));
  ri=check_alloc(n, sizeof(double));

  ret=my_eigen(1, R, n, rr, ri, vr, vi);
  if(ret==-1) {printf("Diagonalization problem: not converged\n"); return -1;}
  if(ret==1) {printf("Error: complex eigen values\n"); return -1;}

  ret=invmat(vr, n, invvr);
  if(ret==0){printf("Error: singular eigen vector matrix\n"); return -1;}

  for(i=0;i<n;i++){
    for(j=0;j<n;j++) ediago[i][j]=0.;
    ediago[i][i]=exp(rr[i]*t);
  }

  matmat(vr, n, n, ediago, n, n, provprod);
  matmat(provprod, n, n, invvr, n, n, P);

  for(i=0;i<n;i++){
    free(vr[i]); free(vi[i]); free(ediago[i]); free(provprod[i]); free(invvr[i]);
  }

  free(rr); free(ri); free(vr); free(vi); free(ediago); free(provprod); free(invvr);

  return;

}



double compute(tree* s_tree, init_option init_o, compute_option compute_opt, converge_option converge_o, print_option print_o){


  int i, j, k, ii, nbseq, lgseq, nbparam, firstcomp, nbreduction, best_gamma;
  int OPTIMIZE, nbiter, *bounded, do_simplex=0, do_simplex_now=0;
  double gc, rem, bound, *toadd, totallength, newtotallength, gamma_l, best_gamma_l;
  double llike, maxllike, **invmatparam, *paramvect, **covmat, *remvect, exact_llike;
  double prec, l, gamma_shp, *useless, lll, dlll;
  double* bestvect, *finalvect;
  char cov[10], type[20];
  static FILE* out=NULL, *outcov=NULL;

int den;
double l1, l2, l_1, l_2, delta, d_gamma, d2_gamma, toadd_gamma, old_gamma_shp;
char *carb;
double**updist, *uplab, d1, d2;
char**upnom;
int upnb;


ttotdeb=time(NULL);

/* printmemory("deb compute");*/

/* probaoutcov=fopen("probaoutcov", "w");*/

  if(print_o->PRINT3){
    if(!out) out=fopen("resfile", "w");
    if(!out) {
      printf("Cannot create file\n");
      exit(EXIT_FAILURE);
    }
  }

  if(print_o->PRINT1==0) noprint=1; else noprint=0;


  nbseq=s_tree->nbseq;
  lgseq=s_tree->lgseq;
  nbparam=s_tree->nb_var_param;

  if(!noprint)
    printf("%d species, %d sites, %d patterns\n", nbseq, (int)s_tree->weightsum, lgseq);

  prec=converge_o->PRECISION;


  compute_o=compute_opt;
  
  OPTIMIZE = compute_o->OPTIMIZE_LENGTH || compute_o->OPTIMIZE_COMP || compute_o->OPTIMIZE_TITV 
	                                || compute_o->OPTIMIZE_ROOT || compute_o->OPTIMIZE_ANC || compute_o->OPTIMIZE_COV;

  bound=pow(10., -1.*(prec));



			/* INIT PARAMETERS */

        /* titv */

  s_tree->titv=init_o->INIT_TITV;


	/* root */

  s_tree->froot=init_o->INIT_ROOT;


	/* Gamma distribution */

  gamma_shp=init_o->INIT_GAMMA;
  gamma_nbcl=init_o->GAMMA_NCL;

  if(gamma_nbcl>1 && !(gamma_shp<=0. && compute_o->OPTIMIZE_GAMMA==0)){
    if(gamma_shp<=0.) gamma_shp=gammamax;
    gamma_clmean=(double*)check_alloc(gamma_nbcl, sizeof(double));
    useless=(double*)check_alloc(gamma_nbcl, sizeof(double));
    DiscreteGamma (useless, gamma_clmean, gamma_shp, gamma_shp, gamma_nbcl, 0);
  }
  else{
    if(gamma_shp<=0.) gamma_shp=gammamax;
    gamma_nbcl=1;
    gamma_clmean=(double*)check_alloc(1, sizeof(double));
    gamma_clmean[0]=1.;
    compute_o->OPTIMIZE_GAMMA=0;
    useless=(double*)check_alloc(gamma_nbcl, sizeof(double));
  }


	/* compositional parameters */

  if(strcmp(init_o->INIT_COMP, "CONST")==0) init4bases_constant(s_tree);
  if(strcmp(init_o->INIT_COMP, "VAR")==0) init4bases_simplemean(NULL, root, lgseq, s_tree->weightsum, s_tree->weight);
  if(strcmp(init_o->INIT_COMP, "BALANCED")==0) init4bases_balanced(NULL, root, lgseq, s_tree->weightsum, s_tree->weight);

  s_tree->omega1=(root->geq-root->ceq)/(root->geq+root->ceq);
  s_tree->omega2=(root->aeq-root->teq)/(root->aeq+root->teq);
/*
   teta=G+C => C=teta/2 * (1-omega1) 
	       G=teta/2 * (1+omega1)
       	       T=(1-teta)/2 * (1-omega2)
       	       A=(1-teta)/2 * (1+omega2)
*/

/* uncomment the following 5 lines for true T92 model */

  s_tree->omega1=s_tree->omega2=0.;
  root->aeq=(root->aeq+root->teq)/2.;
  root->teq=root->aeq;
  root->ceq=(root->ceq+root->geq)/2.;
  root->geq=root->ceq;

  if(init_o->INIT_ANC>=0.){ 
    gc=init_o->INIT_ANC;
    root->aeq=(1.-gc)*(1.+s_tree->omega2)/2;
    root->teq=(1.-gc)*(1.-s_tree->omega2)/2.; 
    root->ceq=gc*(1.-s_tree->omega1)/2.;
    root->geq=gc*(1.+s_tree->omega1)/2.;
  }


  root->a=root->aeq; root->c=root->ceq; root->g=root->geq; root->t=root->teq;
  root->aeq=root->ceq=root->geq=root->teq=-1.;


  s_tree->GCanc=root->g+root->c;
  set_listgc(s_tree, root);


	/* lengths */

  if(strcmp(init_o->INIT_LENGTH, "REDO")==0)
    init_bl(s_tree);

/*
  truelg_to_mylg(root, s_tree->titv);
  set_parambl(s_tree, root);
*/

	/* covarions */

  covar=init_o->INIT_COV;
  pi=init_o->INIT_PI;
  if(covar<covmin && compute_o->OPTIMIZE_COV==1) covar=covmin;
  if(covar>covmax && compute_o->OPTIMIZE_COV==1) covar=covmax;
  if(pi<pimin && compute_o->OPTIMIZE_PI==1) pi=pimin;
  if(pi>pimax && compute_o->OPTIMIZE_PI==1) pi=pimax;
  p1prov=check_alloc(gamma_nbcl, sizeof(proba*));
  p2prov=check_alloc(gamma_nbcl, sizeof(proba*));
  for(i=0;i<gamma_nbcl;i++){
    p1prov[i]=check_alloc(gamma_nbcl, sizeof(proba));
    p2prov[i]=check_alloc(gamma_nbcl, sizeof(proba));
  }
  dp1prov=check_alloc(gamma_nbcl, sizeof(dproba*));
  dp2prov=check_alloc(gamma_nbcl, sizeof(dproba*));
  d2p1prov=check_alloc(gamma_nbcl, sizeof(dproba*));
  d2p2prov=check_alloc(gamma_nbcl, sizeof(dproba*));
  for(i=0;i<gamma_nbcl;i++){
    dp1prov[i]=check_alloc(gamma_nbcl, sizeof(dproba));
    dp2prov[i]=check_alloc(gamma_nbcl, sizeof(dproba));
    d2p1prov[i]=check_alloc(gamma_nbcl, sizeof(dproba));
    d2p2prov[i]=check_alloc(gamma_nbcl, sizeof(dproba));
  }


	/* for exact covarion calculation */
  j=4*gamma_nbcl;
  R=check_alloc(j, sizeof(double*));
  P=check_alloc(j, sizeof(double*));
  for(i=0;i<j;i++){
    R[i]=check_alloc(j, sizeof(double));
    P[i]=check_alloc(j, sizeof(double));
  }
  T92=check_alloc(4, sizeof(double*));
  for(i=0;i<4;i++) T92[i]=check_alloc(4, sizeof(double));


	/* random */

  if(init_o->NBRANDOM>0){
    if(init_o->NOISE<0){
      totallength=newtotallength=0.;
      for(i=0;i<2*nbseq-3;i++)
        totallength+=s_tree->listbr[i]->length;
    
      s_tree->titv=1.+9.*drand48();
      s_tree->froot=drand48();
      s_tree->GCanc=drand48();
      init4bases_random(root);
      for(i=0;i<2*nbseq-3;i++){
        s_tree->listbr[i]->length=drand48();
        newtotallength+=s_tree->listbr[i]->length;
      }
      for(i=0;i<2*nbseq-3;i++){
        s_tree->listbr[i]->length*=(totallength/newtotallength);
      }
    }
    else{
      for(i=0;i<nbparam;i++){
	*(s_tree->var_param[i])*= 1+ (drand48()-0.5)*2.*init_o->NOISE/100.;
      }
    }
  }
  

  toadd=(double*)check_alloc(nbparam, sizeof(double));  
  invmatparam=(double**)check_alloc(nbparam, sizeof(double*));
  for(i=0;i<nbparam;i++)
    invmatparam[i]=(double*)check_alloc(nbparam, sizeof(double));
  paramvect=check_alloc(nbparam+1, sizeof(double));
  bestvect=(double*)check_alloc(nbparam+1, sizeof(double));
  finalvect=(double*)check_alloc(nbparam+1, sizeof(double));
  remvect=(double*)check_alloc(nbparam+1, sizeof(double));
  bounded=(int*)check_alloc(nbparam, sizeof(int));




	/* COMPUTE MAX LIKELIHOOD */

  reset_like(s_tree);
  nbiter=0; rem=0.;
  for(i=0;i<nbparam;i++) bestvect[i]=*(s_tree->var_param[i]);
  exact_cov=0;
  if(covar>0. && print_o->PRINT2) printf("\nApproximate likelihood:\n");

  while(TRUE){ /* start of Newton-Raphson loop */

    nbiter++;

    firstcomp=1;
    nbreduction=0;
    for(i=0;i<nbparam;i++) remvect[i]=*(s_tree->var_param[i]);

/* compute likelihood, compare to previous, */
/* redo (up to 10 times) until it increases */

    while(firstcomp || llike<rem){

      for(i=0;i<nbparam;i++){
        *(s_tree->var_param[i])=remvect[i]+toadd[i];
      }


      for(i=0;i<nbparam;i++) bounded[i]=0;

      boundaries(s_tree, bounded);
      reset_tree(s_tree);

      for(i=0;i<nbparam;i++)
        paramvect[i+1]=*(s_tree->var_param[i]);

      llike=-minus_log_likelihood(paramvect);

      firstcomp=0;
      if(nbiter<GAMMA_RANGE_ITER+2) break;

      if(llike<rem){ 
        if(print_o->PRINT2)
	  printf("likelihood decrease (%f->%f): moving backward\n", rem, llike);
	for(i=0;i<nbparam;i++)
	  toadd[i]/=2.; 
	nbreduction++;
	if(nbreduction==10) {do_simplex_now=1; break;}
      }
    }
    

    if(do_simplex_now){
      do_simplex=1;
      break;
    }

    if(nbiter==1) maxllike=llike;
    else{
      if(maxllike<llike){
	maxllike=llike;
	for(i=0;i<nbparam;i++) bestvect[i]=paramvect[i+1];
      }
    }


/* check positivity */

    if(llike>0.) {
      printf("log negatif\n");
      if(print_o->PRINT3){
        fprintf(out, "Negative log\n");
        output_tree_comp(root, NULL, NULL);
      }
      return 1.;
    }

    if(print_o->PRINT2) printf("log likelihood : %f\n", llike); 

    if(!OPTIMIZE) break;


/* leave if converged */

    if(nbiter!=1 && llike-rem<bound && llike-rem>-bound){	
      for(i=0;i<nbparam;i++)
 	if(bounded[i]){ /* perform simplex if some param is at the*/
	  do_simplex=1;	/* boundary of its definition interval */
	}
      break; 
    }


/* remember likelihood */

    rem=llike;
    s_tree->lkh=llike;

/* first steps: try several gamma distributions */

    old_gamma_shp=gamma_shp;
    if(compute_o->OPTIMIZE_GAMMA && nbiter<GAMMA_RANGE_ITER){
      best_gamma_l=llike; best_gamma=-1;
      for(i=0;i<GAMMA_RANGE_SIZE;i++){
        gamma_shp=gamma_range[i];
        DiscreteGamma (useless, gamma_clmean, gamma_shp, gamma_shp, gamma_nbcl, 0);
        gamma_l=-minus_log_likelihood(paramvect);
        if(gamma_l>best_gamma_l){
          best_gamma=i;
          best_gamma_l=gamma_l;
        }
      }
      if(best_gamma!=-1){
        llike=best_gamma_l;
        rem=llike;
        s_tree->lkh=llike;
        gamma_shp=gamma_range[best_gamma];
      }
      else gamma_shp=old_gamma_shp;
      DiscreteGamma (useless, gamma_clmean, gamma_shp, gamma_shp, gamma_nbcl, 0);
      minus_log_likelihood(paramvect);
    }



/* analytical derivatives (1st + 2nd) */
 	
    derivatives(s_tree);


/* not first steps: */
/* numerical derivative with respect to shape param of gamma distribution, */
/* and modification of gamma shape parameter. */

  if(compute_o->OPTIMIZE_GAMMA && nbiter>GAMMA_RANGE_ITER){
    delta=0.000001;
    lll=s_tree->lkh;
    gamma_shp+=delta;
    DiscreteGamma (useless, gamma_clmean, gamma_shp, gamma_shp, gamma_nbcl, 0);
    l1=-minus_log_likelihood(paramvect);
    gamma_shp+=delta;
    DiscreteGamma (useless, gamma_clmean, gamma_shp, gamma_shp, gamma_nbcl, 0);
    l2=-minus_log_likelihood(paramvect);
    gamma_shp-=3*delta;    
    DiscreteGamma (useless, gamma_clmean, gamma_shp, gamma_shp, gamma_nbcl, 0);
    l_1=-minus_log_likelihood(paramvect);
    gamma_shp-=delta;    
    DiscreteGamma (useless, gamma_clmean, gamma_shp, gamma_shp, gamma_nbcl, 0);
    l_2=-minus_log_likelihood(paramvect);
    gamma_shp+=2*delta;    
    d_gamma=(l1-l_1)/(2*delta);
    d2_gamma=(l2+l_2-2*lll)/(4*delta*delta);
    den=1;
    old_gamma_shp=gamma_shp;
    if(d2_gamma<0.){ 
      toadd_gamma=-d_gamma/d2_gamma;
        while(den<1000){
        gamma_shp=old_gamma_shp+toadd_gamma/den;
 	if(gamma_shp<gammamin) gamma_shp=gammamin;
        DiscreteGamma (useless, gamma_clmean, gamma_shp, gamma_shp, gamma_nbcl, 0);
        lll=-minus_log_likelihood(paramvect);
        if(lll<llike) den*=2;
        else break;
      }
    }
    if(den>1000 || d2_gamma>=0.){
      gamma_shp=old_gamma_shp;
      DiscreteGamma (useless, gamma_clmean, gamma_shp, gamma_shp, gamma_nbcl, 0);
    }
    else llike=lll;
  }



/* modify parameters */

/*
delta=0.0001;
lll=s_tree->lkh;
for(i=0;i<nbparam;i++){
if(s_tree->param_nature[i]==3 || s_tree->param_nature[i]==4) continue;
  paramvect[i+1]+=delta;
  l1=-minus_log_likelihood(paramvect);
  paramvect[i+1]+=delta;
  l2=-minus_log_likelihood(paramvect);
  paramvect[i+1]-=(3*delta);
  l_1=-minus_log_likelihood(paramvect);
  paramvect[i+1]-=delta;
  l_2=-minus_log_likelihood(paramvect);
  paramvect[i+1]+=(2*delta);
  d1=(l1-l_1)/(2*delta);
  d2=(l2+l_2-2*lll)/(4*delta*delta);

  printf("%d, a%e n%e, a%e n%e\n", i, s_tree->deriv[i], d1,  s_tree->deriv2[i], d2);
}
*/

    for(ii=0;ii<s_tree->nb_var_param;ii++){
      toadd[ii]=-s_tree->deriv[ii]/s_tree->deriv2[ii];
    }


    for(i=0;i<nbparam;i++){
      if(s_tree->deriv2[i]>0.){
        toadd[i]=0.;
        if(!print_o->PRINT2) continue;
	switch(s_tree->param_nature[i]){
	  case 0: sprintf(type, "titv"); break;
	  case 1: sprintf(type, "root"); break;
	  case 2: sprintf(type, "GCanc"); break;
	  case 3: sprintf(type, "length"); break;
	  case 4: sprintf(type, "comp"); break;
	  case 5: sprintf(type, "covar"); break;
	  case 6: sprintf(type, "pi"); break;
	}
        printf("param %d (%s) diverging\n", i+1, type);
      }
    }





  } /* end of Newton-Raphson loop */


  if(covar>0.) exact_cov=1;
  exact_llike=-minus_log_likelihood(paramvect);
  if(exact_llike>0.) exact_llike=maxllike;
  if(covar>0. && print_o->PRINT2)
    printf("Exact log likelihood: %f\n\n", exact_llike);
  maxllike=exact_llike;

  if(covar>0.1) do_simplex=1;
  if(covar<=0.1) exact_cov=0;

	/* SIMPLEX */

  if((converge_o->SIMPLEX || covar>0.1) && do_simplex){
    if(print_o->PRINT2) 
      printf("Simplex...\n");
    llike=simplex(bestvect, 0.01, finalvect);
    if(llike>maxllike){
      for(i=0;i<nbparam;i++)
        *(s_tree->var_param[i])=finalvect[i];
      reset_tree(s_tree);
      maxllike=llike;
    }
    else{
      for(i=0;i<nbparam;i++)
        *(s_tree->var_param[i])=bestvect[i];
      reset_tree(s_tree);
      llike=maxllike; 
    }   
  }


	/* OUTPUT */
  
  s_tree->gamma_shp=gamma_shp;
  s_tree->covar=covar;
  s_tree->pi=pi;

  if (print_o->PRINT1){
    printf("\nln(L) = %.6f\n\n", maxllike, nbiter);

    if(compute_o->OPTIMIZE_TITV) printf("Estimated Ts/Tv ratio : %.4f\n", s_tree->titv);
    else printf("Fixed Ts/Tv ratio: %.4f\n", s_tree->titv);

    if(compute_o->OPTIMIZE_ANC) printf("Estimated ancestral GC%% : %.2f%%\n", 100*(root->c+root->g));
    else printf("Fixed ancestral GC% :%.2f%%\n", 100*(root->c+root->g));

    if(compute_o->OPTIMIZE_GAMMA) printf("Estimated Gamma shape parameter : %.4f\n", gamma_shp);
    else if(init_o->INIT_GAMMA>0. && gamma_nbcl>1) printf("Fixed Gamma shape parameter : %.4f\n", gamma_shp);
    else printf("Constant rates accross sites\n");

    if(compute_o->OPTIMIZE_COV) printf("Estimated covarion parameter : %.4f\n", covar);
    else if(init_o->INIT_COV>0.) printf("Fixed covarion parameter : %.4f\n", covar);
    else printf("No covarions\n");

    if(compute_o->OPTIMIZE_PI) printf("Estimated proportion of covarions : %.4f\n", pi);
    else if(init_o->INIT_PI>0.) printf("Fixed proportion of covarions : %.4f\n", pi);



    printf("\n");
  }


  if(print_o->PRINT1 && print_o->PRINT3){
    mindepth(NULL, root);
    output_tree_comp(root, out, s_tree);
    printf("Detailed results are written into file : resfile\n");   
  }

/*
  if((compute_o->OPTIMIZE_COV || init_o->INIT_COV>0.) && print_o->COV_SITE_PATTERN){
    if(covar>covmin){
      output_covar_pattern(s_tree, covar, bestvect, compute_o->OPTIMIZE_COV);
      printf("Site-specific covariation patterns are written into file : covarfile\n");
    }
  }
*/
    

	/* ESTIMATE ACTUAL BASE COMPOSITIONS AND BRANCH LENGTHS  */

  for(k=0;k<gamma_nbcl;k++)
    compute_proba_t3_cov(root, s_tree->titv, 0, 0, 0, 0, 0, k, k);

  set_node_deduced_comp(root);
  mylg_to_truelg(root, s_tree->titv);


	/* FREE */

  free(toadd);
  free(bestvect);
  free(paramvect);
  free(finalvect);
  free(remvect);
  for(i=0;i<s_tree->nb_var_param;i++)
    free(invmatparam[i]);
  free(invmatparam);

/*
  if(0 || outcov){
    for(i=0;i<nbparam;i++) free(covmat[i]);
    free(covmat);
  }
*/
  free(bounded);

  free(gamma_clmean);
  free(useless);

  for(i=0;i<gamma_nbcl;i++){
    free(p1prov[i]); free(p2prov[i]);
    free(dp1prov[i]); free(dp2prov[i]);
    free(d2p1prov[i]); free(d2p2prov[i]);
  }
  free(p1prov); free(p2prov);
  free(dp1prov); free(dp2prov);
  free(d2p1prov); free(d2p2prov);

  j=4*gamma_nbcl;
  for(i=0;i<j;i++) {free(R[i]); free(P[i]);}
  free(R); free(P);
  for(i=0;i<4;i++) free(T92[i]);
  free(T92);

/* printmemory("fin compute");*/

ttotfin=time(NULL);
tottime=difftime(ttotfin, ttotdeb);


if(!noprint){
  printf("Running time: ");
  printtime(tottime);
  printf("\n");
}

  return maxllike;
}








