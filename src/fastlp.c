/*****************************************************************************

                Implementation of the 
		Primal Dual (i.e. Self Dual) Simplex Method Linear Programming Solver for R
		R. Vanderbei & H. Pang, June 2013

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
#define MAX_ITER 1000000

static int status0;
static double lambda0;
static double *x;

int ratio_test0(
	double *dy, 
	int   *idy,
	int    ndy,
	double *y, 
        double *ybar,
	double mu
);


void solver20(
    int m,		/* number of constraints */
    int n,		/* number of variables */
    int nz,		/* number of nonzeros in sparse constraint matrix */
    int *ia, 		/* array row indices */
    int *ka, 		/* array of indices into ia and a */
    double *a,		/* array of nonzeros in the constraint matrix */
    double *b, 		/* right-hand side */
    double *c          /* objective coefficients */
    );



void fastlp(double *obj, double *mat, double *rhs, int *m0 , int *n0, double *opt, int *status, double *lambda)
{

    int m=*m0;		/* number of constraints */
    int n=*n0;		/* number of variables */
    int nz=0;		/* number of nonzeros in sparse constraint matrix */
    int *ia; 		/* array row indices */
    int *ka; 		/* array of indices into ia and a */
    double *a;		/* array of nonzeros in the constraint matrix */
    double *b; 		/* right-hand side */
    double *c;          /* objective coefficients */
    int i, j, k; 
    status0 = *status;
    lambda0 = *lambda;

    if(lambda0<=EPS3){
      lambda0=EPS3;
    }


    MALLOC(        a, m*n+m,  double );      
    MALLOC(       ia, m*n+m,   int );      
    MALLOC(       ka, n+m+1,  int );        
    MALLOC(        c, n,   double );
    MALLOC(        b, m,   double );    

    /*****************************************************************
     * Structure of the problem:

     ***************************************************************/

    for (i=0;i<n;i++){
	c[i]=obj[i];
    }

    for (i=0;i<m;i++){
	b[i]=rhs[i];
    }

    for (i=0;i<m;i++){
       for(j=0;j<n;j++){
	
       }	
    }


    k=0;
	//Sparse matrix representation
    for (j=0; j<n; j++) {
	    ka[j] = k; 
	    for (i=0; i<m; i++) {
		  if (mat[i*n+j]!=0)
		  {
	            a[k] = mat[i*n+j];
		    ia[k] = i;
                    k++;
                    nz++;  
                  
		  }    
	    
	    }
    }
    ka[n]=k;
    solver20(m,n,nz,ia,ka,a,b,c);
    *status=status0;

    for(i=0;i<n;i++){
        opt[i]=x[i];
    }

    FREE(b);
    FREE(a); 
    FREE(ia);
    FREE(ka);
    FREE(c);
    FREE(x);
}


void solver20(
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

	/*structure of the solver*/

    int *basics;
    int *nonbasics;
    int *basicflag;
    double  *x_B;	/* primal basics */
    double  *y_N;	/* dual nonbasics */
    double  *xbar_B;	/* primal basic perturbation */
    double  *ybar_N;    /* dual nonbasic perturbation*/
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
    double  s, t, sbar, tbar, mu=HUGE_VAL;
    double  *vec;
    int    *ivec;
    int     nvec;
    int     N;

    N=m+n;

	 /*******************************************************************
    * read in the data and initialize the common memory sites.
    *******************************************************************/

	//add the slack variables

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
    MALLOC( ybar_N, n,   double );           
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
    CALLOC(   x, N, double );

    /**************************************************************** 
    *  initialization.              				    *
    ****************************************************************/

    atnum(m,N,ka,ia,a,kat,iat,at);

    for (j=0; j<n; j++) {
	nonbasics[j] = j;
	basicflag[j] = -j-1;
	      y_N[j] = -c[j];
           ybar_N[j] = 1;
    }

    for (i=0; i<m; i++) {
	    basics[i] = n+i;
       basicflag[n+i] = i;
	       x_B[i] = b[i];
	    xbar_B[i] = 1;
    }

    lufac( m, ka, ia, a, basics, 0 );

 for (iter=0; iter<MAX_ITER; iter++) {

      /*************************************************************
      * step 1: find mu                                            *
      *************************************************************/
      mu = -HUGE_VAL;
      col_in  = -1;
      for (j=0; j<n; j++) {
		if (ybar_N[j] > EPS2) { 
			if ( mu < -y_N[j]/ybar_N[j] ) {
			     mu = -y_N[j]/ybar_N[j];
			     col_in  = j;
			}
		}
      }
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
     
       if ( mu <= lambda0 ) {	/* optimal */
          status0=0;       
	  break;

      }

        /*************************************************************
	*                          -1  t                             *
	* step 2: compute dy  = -(b  n) e                            * 
	*                   n            i			     *
	*         where i = col_out                                  *
        *************************************************************/
     if ( col_out >= 0 ) {
	vec[0] = -1.0;
	ivec[0] = col_out;
	nvec = 1;

	btsolve( m, vec, ivec, &nvec );  
	Nt_times_y( N, at, iat, kat, basicflag, vec, ivec, nvec, 
		     dy_N, idy_N, &ndy_N );

	col_in = ratio_test0( dy_N, idy_N, ndy_N, y_N, ybar_N,mu );

        /*************************************************************
	* STEP 3: Ratio test to find entering column                 * 
        *************************************************************/

	if (col_in == -1) { 	/* infeasible */
	    status0 = 1;
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

        }

        else {

        /*************************************************************
	*                        -1                                  *
	* STEP 2: Compute dx  = B  N e                               * 
	*                   B         j                              *
        *************************************************************/

	j = nonbasics[col_in];
	for (i=0, k=ka[j]; k<ka[j+1]; i++, k++) {
	     dx_B[i] =  a[k];
	    idx_B[i] = ia[k];
	}
	ndx_B = i;
	bsolve( m, dx_B, idx_B, &ndx_B );

        /*************************************************************
	* STEP 3: Ratio test to find leaving column                  * 
        *************************************************************/

	col_out = ratio_test0( dx_B, idx_B, ndx_B, x_B, xbar_B, mu );

	if (col_out == -1) {	/* UNBOUNDED */
	    status0 = 2;
	    break;
	}

        /*************************************************************
	*                          -1  T                             *
	* STEP 4: Compute dy  = -(B  N) e                            * 
	*                   N            i			     *
	*                                                            *
        *************************************************************/

	 vec[0] = -1.0;
	ivec[0] = col_out;
	nvec = 1;

	btsolve( m, vec, ivec, &nvec );  		
	Nt_times_y( N, at, iat, kat, basicflag, vec, ivec, nvec, 
		     dy_N, idy_N, &ndy_N );

      }

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
      sbar = ybar_N[col_in]/dy_N[k];


      /*************************************************************
      *                                _    _    _                 *
      * step 7: set y  = y  - s dy     y  = y  - s dy              *
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
                ybar_N[j] -= sbar*dy_N[k];
      }
      
      y_N[col_in]    = s;
      ybar_N[col_in] = sbar;

      for (k=0; k<ndx_B; k++) {
		i = idx_B[k];
		x_B[i]    -= t   *dx_B[k];
		xbar_B[i] -= tbar*dx_B[k];

      }

      x_B[col_out]     = t;
      xbar_B[col_out]  = tbar;

      /*************************************************************
      * step 8: update basis                                       * 
      *************************************************************/

      i =    basics[col_out];
      j = nonbasics[col_in];
      basics[col_out]   = j;
      nonbasics[col_in] = i;
      basicflag[i] = -col_in-1;
      basicflag[j] = col_out;


      /*************************************************************
      * step 9: refactor basis and print statistics                *
      *************************************************************/

      refactor( m, ka, ia, a, basics, col_out, v );

  } 

   

      for (i=0; i<m; i++) {
	  x[basics[i]] = x_B[i];
      }


      if(iter>=1){
          Nt_times_y( -1, at, iat, kat, basicflag, vec, ivec, nvec, 
		     dy_N, idy_N, &ndy_N );
      }


    /****************************************************************
    * 	free work space                                             *
    ****************************************************************/

    FREE(  vec );
    FREE( ivec );
    FREE(  x_B );
    FREE(  y_N );
    FREE( dx_B );
    FREE(idx_B );
    FREE( dy_N );
    FREE(idy_N );
    FREE(xbar_B);
    FREE(ybar_N);
    FREE( nonbasics );
    FREE( basics );
    FREE(at);
    FREE(iat);
    FREE(basicflag);
    FREE(kat);

    if(iter>=1){
       lu_clo();
       btsolve(0, vec, ivec, &nvec);
       bsolve(0, vec, ivec, &nvec);
    }


}



int ratio_test0(
	double *dy, 
	int   *idy,
	int    ndy,
	double *y, 
	double *ybar, 
	double mu
)
{
	int j, jj = -1, k;
	double min = HUGE_VAL;

	for (k=0; k<ndy; k++) {
	    if ( dy[k] > EPS1 ) {
	        j = idy[k];
		if ( (y[j] + mu*ybar[j])/dy[k] < min ) {
			min = (y[j] + mu*ybar[j])/dy[k];
			 jj = j;
		}
	    }
	}

	return jj;
}












