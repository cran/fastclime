/*****************************************************************************
Implementation of the Primal Dual (i.e. Self Dual) Simplex Method on Sparse Precision Matrix Estimation
H. Pang, H. Liu & R. Vanderbei, March 2013
******************************************************************************/
         
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <R.h>
#include "myalloc.h"
#include "lu.h"
#include "linalg.h"
#include "macros.h"

#define EPS1 1.0e-8
#define EPS2 1.0e-12
#define EPS3 1.0e-5

static int ColNum; 
static double **mu_mtx;
static int lambda;
static int status;
static double lambda_ratio;
static int **path_mtx;
static double **icov_mtx;
static int maxiter;
static int *max_row_iter;

double sdotprod(double *c, double *x_B, int *basics, int m);

void Nt_times_y( 
    int N, 
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
);

int ratio_test(
	double *dy, 
	int   *idy,
	int    ndy,
	double *y, 
	double mu
);

void solver2(
    int m,		
    int n,		
    int nz,		
    int *ia, 		
    int *ka, 		
    double *a,		
    double *b, 		
    double *c        
    );

void parametric(double *SigmaInput, int *m1, double *mu_input, double *ratio, int *nlambda, int *ipath, int *maxnlambda, double *iicov)
{

    int m;		/* number of constraints */
    int n;		/* number of variables */
    int nz;		/* number of nonzeros in sparse constraint matrix */
    int *ia; 		/* array row indices */
    int *ka; 		/* array of indices into ia and a */
    double *a;		/* array of nonzeros in the constraint matrix */
    double *b; 		/* right-hand side */
    double *c;          /* objective coefficients */
    double **LMATRIX;
    int i, j, k;
    int m0=*m1;
           
    lambda=*nlambda;
    lambda_ratio=*ratio;
    m=2*m0;
    n=m;
    nz=0;
    status=0;
    maxiter=0;

    MALLOC(max_row_iter,m0, int);

    MALLOC(mu_mtx, lambda,  double*);
    for (i=0; i<lambda; i++) {
        CALLOC(mu_mtx[i], m0,  double);
    }


    MALLOC(path_mtx, lambda,  int*);
    for (i=0; i<lambda; i++) {
        CALLOC(path_mtx[i], m0*m0,  int);
    }
    MALLOC(icov_mtx, lambda,  double*);
    for (i=0; i<lambda; i++) {
        CALLOC(icov_mtx[i], m0*m0,  double);
    }
    

    MALLOC(LMATRIX, m,  double*);
    for (i=0; i<m; i++) {
        MALLOC(LMATRIX[i], n,  double);
    }


    for (i=0; i<m0; i++){
         for (j=0; j<m0; j++){
	        LMATRIX[i][j]      =  SigmaInput[i*m0+j];
	        LMATRIX[i+m0][j]   = -SigmaInput[i*m0+j];
		LMATRIX[i][j+m0]   = -SigmaInput[i*m0+j];
		LMATRIX[i+m0][j+m0]=  SigmaInput[i*m0+j];
		if(SigmaInput[i*m0+j]!=0){
		      nz+=4;
		}
	}
    }

    MALLOC(        a, nz+m,  double );      
    MALLOC(       ia, nz+m,   int );      
    MALLOC(       ka, n+m+1,  int );       
    CALLOC(        c, n+m,   double );      

    for (i=0;i<n;i++){
	c[i]=-1.0;	
    }

    k=0;
    for (j=0; j<n; j++) {
	ka[j] = k;
	for (i=0; i<m; i++) {
	   if (LMATRIX[i][j]!=0) {
	     a[k] = LMATRIX[i][j];
             ia[k] = i;
             k++;
           }	    	    
	}
    }
    ka[n]=nz;
    

    for(i=0;i<m;i++){
        FREE(LMATRIX[i])
    }
    FREE(LMATRIX);

   

    for(ColNum=0;ColNum<m0;ColNum++){
          MALLOC(        b, m,   double );  
          for (i=0;i<m;i++){
	      b[i]=0.0;
	  }
	  b[ColNum]=1.0;
	  b[ColNum+m0]=-1.0;
          solver2(m,n,nz,ia,ka,a,b,c);
	  FREE(b);
        
                
    }
      
 
      for(i=0;i<lambda;i++){
	for(j=0;j<m0;j++){
	    mu_input[j*lambda+i]=mu_mtx[i][j];
          
        }
      }
      for(j=0;j<m0*m0;j++){
      for(i=1;i<lambda;i++){             
                if(i>max_row_iter[j/m0]){
                    path_mtx[i][j]=path_mtx[i-1][j];
                    icov_mtx[i][j]=icov_mtx[i-1][j];                
                }            
	    ipath[j*lambda+i]=path_mtx[i][j];
            iicov[j*lambda+i]=icov_mtx[i][j];
         }
      }
    
    *maxnlambda=maxiter;
    
    FREE(a);
    FREE(ia);
    FREE(ka);
    FREE(c);
    for(i=0;i<lambda;i++){
        FREE(mu_mtx[i])
    } 
    FREE(mu_mtx);
 
    for(i=0;i<lambda;i++){
        FREE(path_mtx[i])
    } 
    FREE(path_mtx);
    for(i=0;i<lambda;i++){
        FREE(icov_mtx[i])
    } 
    FREE(icov_mtx);
    FREE(max_row_iter);
    

}


void solver2(
    int m,		/* number of constraints */
    int n,		/* number of variables */
    int nz,		/* number of nonzeros in sparse constraint matrix */
    int *ia, 		/* array row indices */
    int *ka, 		/* array of indices into ia and a */
    double *a,		/* array of nonzeros in the constraint matrix */
    double *b, 		/* right-hand side */
    double *c          /* objective coefficients */
    )
{

    int *basics;
    int *nonbasics;
    int *basicflag;
    double  *x_B;	/* primal basics */
    double  *y_N;	/* dual nonbasics */
    double  *xbar_B;	/* primal basic perturbation */
    double  *dy_N;	/*  dual  basics step direction - values (sparse) */
    int    *idy_N;	/*  dual  basics step direction - row indices */
    int     ndy_N=0;	/* number of nonz in dy_N */
    double  *dx_B;	/* primal basics step direction - values (sparse) */
    int    *idx_B;	/* primal basics step direction - row indices */
    int     ndx_B;	/* number of nonz in dx_B */
    double  *at;	/* sparse data structure for a^t */
    int    *iat;
    int    *kat;
    int     col_in;	/* entering column; index in 'nonbasics' */
    int     col_out;	/* leaving column; index in 'basics' */
    int     iter = 0;	/* number of iterations */
    int     i,j,k,v=0;
    double  s, t, tbar, mu=HUGE_VAL, primal_obj;
    double  *vec;
    int    *ivec;
    int     nvec;
    int     from_scratch;
    int     N;
    double *output_vec;
    int     count;
    double ratio_convert;
  

    N=m+n;
    i = 0;
    k = ka[n];
    for (j=n; j<N; j++) {
	
	a[k] = 1.0;
	ia[k] = i;
	i++;
	k++;
	ka[j+1] = k;
    }
    nz = k;

    MALLOC(    x_B, m,   double );      
    MALLOC( xbar_B, m,   double );      
    MALLOC(   dx_B, m,   double );  
    MALLOC(    y_N, n,   double );      
    MALLOC(   dy_N, n,   double );  
    MALLOC(    vec, N,   double );
    MALLOC(   ivec, N,    int );
    MALLOC(  idx_B, m,    int );      
    MALLOC(  idy_N, n,    int );      
    MALLOC(     at, nz,  double );
    MALLOC(    iat, nz,   int );
    MALLOC(    kat, m+1,  int );
    MALLOC(   basics,    m,   int );      
    MALLOC(   nonbasics, n,   int );      
    MALLOC(   basicflag, N,   int );
    CALLOC(   output_vec, N, double );

    /**************************************************************** 
    *  initialization.              				    *
    ****************************************************************/

    atnum(m,N,ka,ia,a,kat,iat,at);	

    for (j=0; j<n; j++) {
	nonbasics[j] = j;
	basicflag[j] = -j-1;
	      y_N[j] = -c[j];  
    }

    for (i=0; i<m; i++) {
	    basics[i] = n+i;
       basicflag[n+i] = i;
	       x_B[i] = b[i];
	    xbar_B[i] = 1;
    }


    lufac( m, ka, ia, a, basics, 0 );

    for (iter=0; iter<lambda; iter++) {
       count=0;
       if(iter>maxiter){
           maxiter=iter;
       }

      /*************************************************************
      * step 1: find mu                                            *
      *************************************************************/
    
      mu = -HUGE_VAL;
      col_in  = -1;
      col_out = -1;
      for (i=0; i<m; i++) {
		if (xbar_B[i] > EPS2) { 
			if ( mu < -x_B[i]/xbar_B[i] ) {
			     mu = -x_B[i]/xbar_B[i];                                                    
			     col_out = i;
			     col_in  = -1;
			}
		}
      }  
      
      mu_mtx[iter][ColNum]=mu;
          
      for (i=0; i<m; i++){
         output_vec[basics[i]] = x_B[i];
      }
      for(i=0;i<m/2;i++){	      
         if((output_vec[i]-output_vec[i+m/2])>EPS3){
            count++;         
            path_mtx[iter][m/2*ColNum+i]=1;
            icov_mtx[iter][m/2*ColNum+i]=output_vec[i]-output_vec[i+m/2];           
         }
      }
       
    
      ratio_convert=(double)count/((m/2)-1);
    

      if(ratio_convert>=lambda_ratio){
       
          break;
          status=2;
      }

      if ( mu <= EPS3 ) {	/* optimal, we only need primal feasible */
          status=3;       
	  break;

      }

        /*************************************************************
	*                          -1  t                             *
	* step 2: compute dy  = -(b  n) e                            * 
	*                   n            i			     *
	*         where i = col_out                                  *
        *************************************************************/

	vec[0] = -1.0;
	ivec[0] = col_out;
	nvec = 1;
	btsolve( m, vec, ivec, &nvec ); 
	Nt_times_y( N, at, iat, kat, basicflag, vec, ivec, nvec, 
		     dy_N, idy_N, &ndy_N );

        /*************************************************************
	* step 3: ratio test to find entering column                 * 
        *************************************************************/

	col_in = ratio_test( dy_N, idy_N, ndy_N, y_N, mu );
	if (col_in == -1) { 	/* infeasible */
	    break;
	}

        /*************************************************************
	*                        -1                                  *
	* step 4: compute dx  = b  n e                               * 
	*                   b         j                              *
	*                                                            *
        *************************************************************/
      
	j = nonbasics[col_in];
	for (i=0, k=ka[j]; k<ka[j+1]; i++, k++) {
	     dx_B[i] =  a[k];
	    idx_B[i] = ia[k];
	}
	ndx_B = i;
	bsolve( m, dx_B, idx_B, &ndx_B );	 

      /*************************************************************
      *                                                            *
      * step 5: put       t = x /dx                                *
      *                        i   i                               *
      *                   _   _                                    *
      *                   t = x /dx                                *
      *                        i   i                               *
      *                   s = y /dy                                *
      *                        j   j                               *
      *                   _   _                                    *
      *                   s = y /dy                                *
      *                        j   j                               *
      *************************************************************/
     
      for (k=0; k<ndx_B; k++) if (idx_B[k] == col_out) break;

      t    =    x_B[col_out]/dx_B[k];
      tbar = xbar_B[col_out]/dx_B[k];


      for (k=0; k<ndy_N; k++) if (idy_N[k] == col_in) break;

      s    =    y_N[col_in]/dy_N[k];

      /*************************************************************
      *                                _    _    _                 *
      * step 6: set y  = y  - s dy     y  = y  - s dy              *
      *              n    n       n     n    n       n             *
      *                                _    _                      *
      *             y  = s             y  = s                      *
      *              i                  i                          *
      *             _    _    _                                    *
      *             x  = x  - t dx     x  = x  - t dx              *
      *              b    b       b     b    b       b             *
      *             _    _                                         *
      *             x  = t             x  = t                      *
      *              j                  j                          *
      *************************************************************/
     
      for (k=0; k<ndy_N; k++) {
		j = idy_N[k];
		y_N[j]    -= s   *dy_N[k];
        y_N[col_in]    = s;
      }

      for (k=0; k<ndx_B; k++) {
		i = idx_B[k];
		x_B[i]    -= t   *dx_B[k];
		xbar_B[i] -= tbar*dx_B[k];

      }

      x_B[col_out]     = t;
      xbar_B[col_out]  = tbar;


      /*************************************************************
      * step 7: update basis                                       * 
      *************************************************************/   
      i =    basics[col_out];
      j = nonbasics[col_in];
      basics[col_out]   = j;
      nonbasics[col_in] = i;
      basicflag[i] = -col_in-1;
      basicflag[j] = col_out;

      /*************************************************************
      * step 8: refactor basis                                     *
      *************************************************************/
      from_scratch = refactor( m, ka, ia, a, basics, col_out, v );
    
    }
   
    max_row_iter[ColNum]=iter;
    Nt_times_y( -1, at, iat, kat, basicflag, vec, ivec, nvec, 
		     dy_N, idy_N, &ndy_N );

    FREE(  vec );
    FREE( ivec );
    FREE(  x_B );
    FREE(  y_N );
    FREE(xbar_B);
    FREE( dx_B );
    FREE(idx_B );
    FREE( dy_N );
    FREE(idy_N );
    FREE( nonbasics );
    FREE( basics );
    FREE( output_vec );
    FREE(at);
    FREE(iat);
    FREE(basicflag);
    FREE(kat);
    lu_clo();
    btsolve(0, vec, ivec, &nvec);
    bsolve(0, vec, ivec, &nvec);

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

int ratio_test(
	double *dy, 
	int   *idy,
	int    ndy,
	double *y,  
	double mu
)
{
	int j, jj = -1, k;
	double min = HUGE_VAL;

	for (k=0; k<ndy; k++) {
	    if ( dy[k] > EPS1 ) {
	        j = idy[k];
		if ( y[j]/dy[k] < min ) {
			min = y[j]/dy[k];
			 jj = j;			
		}
	    }
	}

	return jj;
}

double sdotprod(double *c, double *x_B, int *basics, int m)
{
	int i;
	double prod = 0.0;
	for (i=0; i<m; i++) { prod += c[basics[i]]*x_B[i]; }
	return prod;
}







