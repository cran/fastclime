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
static int N;
static int lambda;
static int *max_row_iter;


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
    double *c,
    double *lambdamin,
    int *maxnlambda,
    double *mu_input,
    double *iicov        
    );

void parametric(double *SigmaInput, int *m1, double *mu_input, double *lambdamin, int *nlambda, int *maxnlambda, double *iicov)
{
  //printf("hello6 \n");

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
  int m0 = *m1;


  /* lambdamin is the smallest lambda asked by the user */
  /* nlambda is the maximum iteration asked by the user */
  /* maxnlambda is the maxium number of iteration found among each column */
      
  lambda = *nlambda;

  /* Actual matrix dimension is the 2m1*2m1, so m=2*m1, n=2*m1 */
  m = 2*m0;
  n = 2*m0;
  nz = 0;
  N = m+n;


  MALLOC(max_row_iter,m0, int);
    

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
      if(SigmaInput[i*m0+j] != 0){
		    nz+=4;
	   	}
	   }
  }



  MALLOC(        a, nz+m,  double );      
  MALLOC(       ia, nz+m,   int );      
  MALLOC(       ka, n+m+1,  int );       
  CALLOC(        c, n+m,   double );      

/*Form vector c*/
  for ( i = 0; i < n; i++){
	  c[i] = -1.0;	
  }

/*Form sparse matrix A*/
  k=0;


  for (j = 0; j < n; j++) {
	  ka[j] = k;
	  for (i = 0; i < m; i++) {
	    if (LMATRIX[i][j]!=0.0) {
	      a[k] = LMATRIX[i][j];
        ia[k] = i;
        k++;
      }
	  }
  }
  ka[n] = k;

  for( i = 0; i < m; i++){
    FREE(LMATRIX[i])
  }
  FREE(LMATRIX);

  /* Add slack variable to the sparse matrix A */
  i = 0;
  for (j = n; j < N; j++) {
    a[k] = 1.0;
    ia[k] = i;
    i++;
    k++;
    ka[j+1] = k;
  }
  nz = k;


  for(ColNum = 0; ColNum < m0; ColNum++){

    /* Form vector b for each problem */
    MALLOC(        b, m,   double );  
    for (i = 0; i < m; i++){
	    b[i] = 0.0;
	  }

	  b[ColNum] = 1.0;
	  b[ColNum+m0] = -1.0;
    
    /* Call the parametric simplex method solver here */
    solver2(m,n,nz,ia,ka,a,b,c,lambdamin,maxnlambda,mu_input,iicov);
	  FREE(b);
        
  }
      

  for(j = 0; j < m0*m0; j++){
    for(i = 1; i < lambda; i++){

      if(i > max_row_iter[j/m0]){
        iicov[j*lambda+i] = iicov[j*lambda+i-1];                
      }            

    }
  }
    
    
  FREE(a);
  FREE(ia);
  FREE(ka);
  FREE(c);

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
    double *c,          /* objective coefficients */
    double *lambdamin,
    int *maxnlambda,
    double *mu_input,
    double *iicov
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
  double  s, t, tbar, mu=HUGE_VAL;
  double  *vec;
  int    *ivec;
  int     nvec;
  int     status;
  double *output_vec;

  

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
 

  for (iter = 0; iter < lambda; iter++) {

    CALLOC(   output_vec, N, double );
    if(iter > *maxnlambda){
      *maxnlambda = iter;
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
      
    mu_input[lambda*ColNum+iter] = mu;
    
    /*Find the current basic solution and store it to output_vec*/      
    for (i=0; i<m; i++){
      output_vec[basics[i]] = x_B[i] + mu*xbar_B[i];
    }


    for(i=0; i < m/2; i++){	      
      if(fabs(output_vec[i]-output_vec[i+m/2])>EPS3){
        iicov[ColNum * lambda * (m/2) + i * lambda + iter] = output_vec[i]-output_vec[i+m/2]; 
      }
    }
        

    if(mu <= *lambdamin){  /*mu becomes smaller than the required lambdamin, stop*/
      status = 0;
      break;
          
    }

    if ( mu <= EPS3 ) {	/* optimal, we only need primal feasible */ 
      status = 1;                    
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
	  if (col_in == -1) { 	/* infeasible for the current value of mu, stop*/
      status = 2;       
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
    refactor( m, ka, ia, a, basics, col_out, v );
    FREE( output_vec );
    
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









