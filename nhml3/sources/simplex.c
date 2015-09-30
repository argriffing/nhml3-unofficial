
#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#define NR_END 1
#define FREE_ARG char*
#include <math.h>
#define NRANSI
/* #include "nrutil.h" */
#define NMAX 10000
#define GET_PSUM \
					for (j=1;j<=ndim;j++) {\
					for (sum=0.0,i=1;i<=mpts;i++) sum += p[i][j];\
					psum[j]=sum;}
#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}


void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

double *vector(int nl, int nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}


double **matrix(int nrl, int nrh, int ncl, int nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	int i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}


void free_vector(double *v, int nl, int nh)
/* free a double vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}


void free_matrix(double **m, int nrl, int nrh, int ncl, int nch)
/* free a double matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}



 double amotry_cont( double **p,  double y[],  double psum[], int ndim,
	 double (*funk)( double []), int ihi,  double fac, double** contraintes)
{
	int j;
	double fac1,fac2,ytry,*ptry;

	ptry=vector(1,ndim);
	fac1=(1.0-fac)/ndim;
	fac2=fac1-fac;
	for (j=1;j<=ndim;j++) {
          ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;
          if(ptry[j]<contraintes[j][0]) ptry[j]=contraintes[j][0];
          else if(ptry[j]>contraintes[j][1]) ptry[j]=contraintes[j][1];
        }
	ytry=(*funk)(ptry);
	if (ytry < y[ihi]) {
/* amotry : ihi ameliore => modification du simplex */
		y[ihi]=ytry;
		for (j=1;j<=ndim;j++) {
			psum[j] += ptry[j]-p[ihi][j];
			p[ihi][j]=ptry[j];
		}
	}
/*else printf("amotry : ihi pas ameliore\n");*/
	free_vector(ptry,1,ndim);
	return ytry;
}


void amoeba_cont( double **p,  double *y, int ndim,  double ftol,
	 double (*funk)( double []),  int *nfunk, double** contraintes)
{
	 double amotry_cont( double **p,  double y[],  double psum[], int ndim,
		 double (*funk)( double []), int ihi,  double fac, double** contraintes);
	int i,ihi,ilo,inhi,j,mpts=ndim+1, step=0, centaine=0;
	 double rtol,sum,swap,ysave,ytry,*psum;
char chaine[10];


	psum=vector(1,ndim);
	*nfunk=0;

for (j=1;j<=ndim;j++) {
  for (sum=0.0,i=1;i<=mpts;i++) 
    sum += p[i][j];
  psum[j]=sum;
}

	/*GET_PSUM*/
	for (;;) {
	   	step++;
		ilo=1;
		ihi = y[1]>y[2] ? (inhi=2,1) : (inhi=1,2);
		for (i=1;i<=mpts;i++) {
			if (y[i] <= y[ilo]) ilo=i;
			if (y[i] > y[ihi]) {
				inhi=ihi;
				ihi=i;
			} else if (y[i] > y[inhi] && i != ihi) inhi=i;
		}
		rtol=fabs(y[ihi]-y[ilo]);

		if (rtol < ftol) {
			SWAP(y[1],y[ilo])
			for (i=1;i<=ndim;i++) SWAP(p[1][i],p[ilo][i])
			break;
		}
		if (*nfunk >= NMAX) 
                {
                 /* printf("\namoeba::attention NMAX exceeded"); */
                  free_vector(psum,1,ndim);
                  return;
                }
		*nfunk += 2;
/*if(PRINT && *nfunk>20*centaine){
  printf("%de appel: min=%f  max=%f\n", *nfunk, -y[ilo], -y[ihi]);
  centaine++;
}*/
		ytry=amotry_cont(p,y,psum,ndim,funk,ihi,-1.0, contraintes);
/* printf("1e essai (inversion) : ytry=%f\n", ytry); */
		if (ytry <= y[ilo]){
/* printf("ilo ameliore : on en rajoute dans cette direction\n"); */
			ytry=amotry_cont(p,y,psum,ndim,funk,ihi,2.0, contraintes);
		}
		else if (ytry >= y[inhi]) {
/* printf("ytry >= y[inhi] : reduction du simplex\n"); */
			ysave=y[ihi];
			ytry=amotry_cont(p,y,psum,ndim,funk,ihi,0.5, contraintes);
			if (ytry >= ysave) {
/* printf("toujours pas mieux : on rapproche tous les points du min\n"); */
				for (i=1;i<=mpts;i++) {
					if (i != ilo) {
						for (j=1;j<=ndim;j++){
						  p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);
						  if(p[i][j]<contraintes[j][0]) 
						    p[i][j]=psum[j]=contraintes[j][0];
          					  else if(p[i][j]>contraintes[j][1]) 
						    p[i][j]=psum[j]=contraintes[j][1];
						}
						y[i]=(*funk)(psum);
					}
				}
				*nfunk += ndim;
				GET_PSUM
			}
		} else --(*nfunk);
	}
	free_vector(psum,1,ndim);
}
#undef SWAP
#undef GET_PSUM
#undef NMAX
#undef NRANSI
