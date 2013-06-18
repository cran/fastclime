#include <stdlib.h>
#include "linalg.h"
#include "myalloc.h"
#include "macros.h"

#define EPS1 1.0e-8
/*---------------------------------------------------------------+
|  inner product between n-vectors x and y                      */

double dotprod( double *x, double *y, int n)         
{
        int i; 
        double dotprod=0.0e0;

        for (i=0; i<n; i++) dotprod += x[i]*y[i];

        return (dotprod);
}

/*---------------------------------------------------------------+
|  y = basis submatrix of (a,ka,ia) times x                     */

void bmx( int m, double *a, int *ka, int *ia, int *basis, 
          double *x, double *y)    
{
        int i,j,k;

        for (i=0; i<m; i++) y[i] = 0.0e0;
        for (i=0; i<m; i++) {
                j = basis[i];
                for (k=ka[j]; k<ka[j+1]; k++)
                        y[ia[k]] += a[k]*x[i];
        }
}

/*---------------------------------------------------------------+
|  y = basis submatrix of (a,ka,ia) transpose times x           */

void btmx( int m, double *a, int *ka, int *ia, int *basis, 
          double *x, double *y)    
{
        int i,j,k;

        for (i=0; i<m; i++) y[i] = 0.0e0;
        for (i=0; i<m; i++) {
                j = basis[i];
                for (k=ka[j]; k<ka[j+1]; k++)
                        y[i] += a[k]*x[ia[k]];
        }
}

/*---------------------------------------------------------------+
|  y = sparse matrix (a,ka,ia) times x                          */

void smx( int m, int n, double *a, int *ka, int *ia, double *x, double *y)    
{
        int i,j,k;

        for (i=0; i<m; i++) y[i] = 0.0e0;
        for (j=0; j<n; j++) 
                for (k=ka[j]; k<ka[j+1]; k++)
                        y[ia[k]] += a[k]*x[j];
}

/*---------------------------------------------------------------+
|  (kat,iat,at) = transpose of (ka,ia,a)                        */

void atnum( int m, int n, int *ka, int *ia, double *a,      
        int *kat,int *iat, double *at
        )  
{
        int i,j,k,row,addr;
        int *iwork;

        CALLOC( iwork, m, int );

        for (k=0; k<ka[n]; k++) {
                row = ia[k];
                iwork[row]++;
        }
        kat[0] = 0;
        for (i=0; i<m; i++) {
                kat[i+1] = kat[i] + iwork[i];
                iwork[i] = 0;
        }
        for (j=0; j<n; j++) {
                for (k=ka[j]; k<ka[j+1]; k++) {
                        row = ia[k];
                        addr = kat[row] +iwork[row];
                        iwork[row]++;
                        iat[addr] = j;
                        at[addr]  = a[k];
                }
        }
        FREE( iwork );
}

/*---------------------------------------------------------------+
|  compute componentwise maximum of n-vector x                  */

double maxv( double *x, int n)    
{
        int i;
        double maxv=0.0e0;

        for (i=0; i<n; i++) maxv = MAX(maxv, ABS(x[i]));

        return (maxv);
}


double sdotprod(double *c, double *x_B, int *basics, int m)
{
	int i;
	double prod = 0.0;
	for (i=0; i<m; i++) { prod += c[basics[i]]*x_B[i]; }
	return prod;
}


void Nt_times_y( 
    int n, 
    double *at, 
    int *iat, 
    int *kat, 
    int *basicflag,
    double *y, 
    int *iy, 
    int ny, 
    double *yN,
    int *iyN,
    int *pnyN
)
{
    int i,j,jj,k,kk;

    static double *a=NULL;
    static int  *tag=NULL;
    static int *link=NULL;
    static int  currtag=1;

    if (n == -1) {
	if (a != NULL) FREE(a);
	if (tag != NULL) FREE(tag);
	link--;
	if (link != NULL) FREE(link);
	return;
    }


    if (  a  == NULL) MALLOC(  a,  n,   double);
    if ( tag == NULL) CALLOC( tag, n,   int);
    if (link == NULL) {CALLOC(link, n+2, int); link++;}

    jj = -1;
    for (k=0; k<ny; k++) {
	i = iy[k];

	for (kk=kat[i]; kk<kat[i+1]; kk++) {
	    j = iat[kk];
	    if (basicflag[j] < 0) {
		if (tag[j] != currtag) {
		    a[j] = 0.0;
		    tag[j] = currtag;
		    link[jj] = j;
		    jj = j;

		}
		a[j] += y[k]*at[kk];

	    }

	}
    }

    link[jj] = n;
    currtag++;

    k = 0;

    for (jj=link[-1]; jj<n; jj=link[jj]) {
	if ( ABS(a[jj]) > EPS1 ) {
             yN[k] = a[jj];
            iyN[k] = -basicflag[jj]-1;
            k++;
	}

    }
 
    *pnyN = k;
}
