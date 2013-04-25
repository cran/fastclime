double dotprod( double *x, double *y, int n);
void bmx(int m, double *a, int *ka, int *ia, int *basis, double *x, double *y);
void btmx(int m, double *a, int *ka, int *ia, int *basis, double *x, double *y);
void smx( int m, int n, double *a, int *ka, int *ia, double *x, double *y);
void atnum( int m, int n, int *ka, int *ia, double *a,      
        int *kat,int *iat, double *at
        );
double maxv( double *x, int n);
double sdotprod(double *c, double *x_B, int *basics, int m);
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
);
