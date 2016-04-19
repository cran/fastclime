/*****************************************************************************

                Implementation of the 
		Primal Dual (i.e. Self Dual) Simplex Method on Dantzig Selector
		R. Vanderbei & H. Pang, Novermeber 2012

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


int ratio_test2(
  double *dy, 
  int   *idy,
  int    ndy,
  double *y, 
  double mu
);


void dantzig(double *X2, double *Xy, double *BETA0, int *d0, 
  double *lambda, int *nlambda, double *lambdalist)
{


  int m;    /* number of constraints */
  int n;    /* number of variables */
  int nz;   /* number of nonzeros in sparse constraint matrix */
  int *ia;    /* array row indices */
  int *ka;    /* array of indices into ia and a */
  double *a;    /* array of nonzeros in the constraint matrix */
//  double **LMATRIX;
 

  int *basics;
  int *nonbasics;
  int *basicflag;
  double  *x_B; /* primal basics */
  double  *y_N; /* dual nonbasics */
  double  *xbar_B;  /* primal basic perturbation */
  double  *dy_N;  /*  dual  basics step direction - values (sparse) */
  int    *idy_N;  /*  dual  basics step direction - row indices */
  int     ndy_N=0;  /* number of nonz in dy_N */
  double  *dx_B;  /* primal basics step direction - values (sparse) */
  int    *idx_B;  /* primal basics step direction - row indices */
  int     ndx_B;  /* number of nonz in dx_B */
  double  *at;  /* sparse data structure for a^t */
  int    *iat;
  int    *kat;
  int     col_in; /* entering column; index in 'nonbasics' */
  int     col_out;  /* leaving column; index in 'basics' */
  int     iter = 0; /* number of iterations */
  int     i,j,k,v=0;
  double  s, t, tbar, mu=HUGE_VAL;
  double  *vec;
  int    *ivec;
  int     nvec;
  int     N;
  double *output_vec;
  int d;
  //double temp_sum;
  //double *temp_vec;

   
  d=*d0;         
  m=2*d;
  n=2*d;
  nz=m*n;
    
  MALLOC(        a, nz+m,  double );      
  MALLOC(       ia, nz+m,   int );      
  MALLOC(       ka, n+m+1,  int );    


   // MALLOC(LMATRIX, d,  double*);
   //  for (i=0; i<d; i++) {
   //      MALLOC(LMATRIX[i], n,  double);
   //  }


   //  for (i=0; i<d; i++){
   //       for (j=0; j<d; j++){
   //        LMATRIX[i][j]      =  X2[i*d+j];
   //        printf("LMATRIX[%d][%d]= %e \n",i, j, output_vec[i]);
   //    }
   //  }   
  
  //Form sparse matrix

  k=0;
  for (j=0; j<n; j++) {
    ka[j] = k;
    for (i=0; i<m; i++) {
      if (i<d && j<d) {
        a[k] = X2[i*d+j];
      }
      else if (i<d && j>=d) {
        a[k] = -X2[i*d+j-d];
      }
      else if(i>=d && j<d) {
        a[k] = -X2[(i-d)*d+j];
      }
      else{
        a[k] = X2[(i-d)*d+j-d];
      }
        ia[k] = i;
        k++;
    }
  }


  ka[n]=k;
    
//Add slack variables to matrix

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

  CALLOC(    x_B, m,   double );      
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
    *  initialization.                          *
    ****************************************************************/

  atnum(m,N,ka,ia,a,kat,iat,at);  

  for (j=0; j<n; j++) {
    nonbasics[j] = j;
    basicflag[j] = -j-1;
    //y_N[j] = -c[j];
    y_N[j] = 1.0; 
  }

  for (i=0; i<m; i++) {
    basics[i] = n+i;
    basicflag[n+i] = i;
    //x_B[i] = Xy[i];
    xbar_B[i] = 1.0;
  }

  for (i=0; i<d; i++) {
    x_B[i] = Xy[i];
    x_B[i+d]=-Xy[i];
    //printf("x_B[%d]= %e \n",i, x_B[i]);
  }

  lufac( m, ka, ia, a, basics, 0 );

  //Start the Main loop of PSM in the solver
  for (iter=0; iter<*nlambda; iter++) {
    //printf("iter=%d \n",iter);
      /*************************************************************
      * step 1: find mu                                            *
      *************************************************************/
    CALLOC(   output_vec, N, double );
    //MALLOC(   temp_vec, d,   double );

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
     
      
    lambdalist[iter]=mu;
    //printf("mu=%e \n",mu);  
    //printf("lambdalist[%d]=%e \n", iter, lambdalist[iter]);
          
    for (i=0; i<m; i++){
      output_vec[basics[i]] = x_B[i]+mu*xbar_B[i];

      //printf("x_B[%d]= %e \n",i, x_B[i]);
      //printf("output_vec[basics[%d]]= %e \n",i, output_vec[basics[i]]);
    }

    // for (i=0; i<m; i++){
    //   printf("output_vec[%d]= %e \n",i, output_vec[i]);
    // }

    for(i=0; i<m/2; i++){
        //printf("output_vec[%d]= %e \n",i, output_vec[i]);
        //printf("output_vec[%d]= %e \n",i+m/2, output_vec[i+m/2]);       
      if(fabs(output_vec[i]-output_vec[i+m/2])>EPS3){
        //printf("output_vec[%d]= %e \n",i, output_vec[i]);
        //printf("output_vec[%d]= %e \n",i+m/2, output_vec[i+m/2]);
        BETA0[m/2*iter+i]=output_vec[i]-output_vec[i+m/2];           
      }

    }

    // for (i = 0; i < d; i++) {
    // temp_sum = 0.0;
    // // requires a single access to the i-th row
    //   for (j = 0; j < d; j++) {
    //     temp_sum += LMATRIX[i][j] *(output_vec[j]-output_vec[j+m/2]);
    //   }
    //    temp_vec[i]=temp_sum;// requires a single access of element c[i]
    //    printf("result[%d]= %e \n",i, temp_sum-Xy[i]);
       
    // }
    // printf("lambda= %e \n",mu);
        

    if(mu <= *lambda){ /*Find the solution corresponding to the required lambda*/
      break;
    }

        /*************************************************************
  *                          -1  t                             *
  * step 2: compute dy  = -(b  n) e                            * 
  *                   n            i           *
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
    col_in = ratio_test2( dy_N, idy_N, ndy_N, y_N, mu );
    if (col_in == -1) {   /* infeasible */
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
    FREE( output_vec );
    //FREE(temp_vec);
    refactor( m, ka, ia, a, basics, col_out, v );
    
  }

  Nt_times_y( -1, at, iat, kat, basicflag, vec, ivec, nvec, 
          dy_N, idy_N, &ndy_N );
     

  FREE(a);
  FREE(ia);
  FREE(ka);
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


int ratio_test2(
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











