#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define F77(a) d ## a ## _
typedef int INT;
typedef double R;

void F77(steqr)(char *, int *, R *, R *, R *, int *, R *, int *);

/*  0<=i<=m-1, 0<=j<=n-1  without specifying

    i,                            j
    i,   j

    i,                            0            i<m-1
    i,   n 
    i+1, n+1
*/

void LegendMat(INT n, R **b, R **U)
{ INT i, j, k, info;
  R s, t; 
	
/*  basis  
    phi[k] = (L[k]-L[k-2])/sqrt(4*k-2)  2<=k<=N
           = (1-x)/2                      k = 0                            
           = (1+x)/2                      k = 1
*/                                                                         

   *U = (R *) malloc((n-1)*(n/2)*sizeof(R));
   *b = (R *) malloc((3*n+3)*sizeof(R));
   
   i  = n/2;
   for ( k=2; k<=n; k++ ) { /*  0,   n/2-1;   n/2,     n-2   */ 
   	   j       = k/2+(k%2?i-1:-1);
   	   	
   	   (*b)[j] = .5/((k+.5)*(k-1.5));
   	   
   	   if  ( k<=n-2 )       /*  n-1, n+n/2-2; n+n/2-1, 2*n-3 */
   	   (*b)[n-1+j] = -.25/(sqrt((k-.5)*(k+1.5))*(k+.5));
   }
   
   /* symmetic tridiagonal eigenvalues and eigenvectors
      calling the fortran interface
   */
   n--;
   if ( n>=1 ) F77(steqr)("I", &i, *b,   *b+n,   *U,   &n, *b+2*n, &info);
   if ( n>=2 ) F77(steqr)("I", &j, *b+i, *b+i+n, *U+i, &n, *b+2*n, &info);
   n++;

   *b += (n+1);
   
   s = -1./sqrt(6.);
   for ( k=0; k<i; k++ ) {     
   	   t           = s*(*U)[k*(n-1)];
   	   (*b)[k]     = t;          /* (L0+L1)/2, psie */
   	   (*b)[k+n+1] = t;          /* (L0-L1)/2, psie */
   }
   
   s = 1./sqrt(90.);
   for ( k=i; k<n-1; k++ ) {     
   	   t           = s*(*U)[i+(k-i)*(n-1)];
   	   (*b)[k]     = -t;        /* (L0+L1)/2, psio */
   	   (*b)[k+n+1] =  t;        /* (L0-L1)/2, psio */
   }
      
   (*b)[n-1]   =  2./3.;
   (*b)[n]     =  1./3.;
   (*b)[2*n]   =  1./3.;
   (*b)[2*n+1] =  2./3.;

   *b -= (n+1);
   
   return;                                                                     
}                                                                            

void assemble1d(int m, int n, R *h, R *b, int **iu, int **ju, R **vu, 
                int **iv, int **jv, R **vv)
{ int i, j, k, p, q, info;
  
  p   = m*n-1;
  q   = m*(n+2)-(m>1?5:3);
  *iu = (INT *) malloc((p+2)*sizeof(INT));
  *ju = (INT *) malloc(q*sizeof(INT));
  *vu = (R   *) malloc(q*sizeof(R));
  
  (*iu)[0] = p;
  (*iu)++;
  (*iu)[0] = 0;
  p = q = 0;
  
  for ( i=0; i<m;   i++ ) 
  for ( j=0; j<n-1; j++ ) {
  	
      (*ju)[q]   = p;
      (*vu)[q++] = 1./h[i];
      (*iu)[++p] = q;
  }
  
  for ( i=0; i<m-1; i++ ) {
      if  ( i>0 ) {
          (*ju)[q]   = m*(n-1)+i-1;
          (*vu)[q++] = -0.5/h[i];
      }
      
      (*ju)[q]   = m*(n-1)+i;
      (*vu)[q++] = 0.5/h[i]+0.5/h[i+1];
      
      if  ( i<m-2 ) {
          (*ju)[q]   = m*(n-1)+i+1;
          (*vu)[q++] = -0.5/h[i+1];
      }
      
      (*iu)[++p] = q;
  }
  (*iu)--;
   printf("nnz(S) = [%d, %d]\n",q,m*(n+2)-(m>1?5:3));
  
  p   = m*n-1;
  q   = 5*m*n-2*m-4*n+1-(m>1?2:0);
  *iv = (INT *) malloc((p+2)*sizeof(INT));
  *jv = (INT *) malloc(q*sizeof(INT));
  *vv = (R   *) malloc(q*sizeof(R));
  
  (*iv)[0] = p;
  (*iv)++;
  (*iv)[0] = 0;
  p = q = 0;
  
  for ( i=0; i<m;   i++ ) 
  for ( j=0; j<n-1; j++ ) {
  	
      (*jv)[q]   = p;
      (*vv)[q++] = b[j]*h[i];
      
      if  ( i>0 ) {
          (*jv)[q]   = m*(n-1)+(i-1);
          (*vv)[q++] = b[j+2*n+2]*h[i];
      }
      
      if  ( i<m-1 ) {
          (*jv)[q]   = m*(n-1)+i;
          (*vv)[q++] = b[j+n+1]*h[i];
      }
      
      (*iv)[++p] = q;
  }
  
  for ( i=0; i<m-1; i++ ) {
  	  for ( j=0; j<n-1; j++ ) {/* i */
          (*jv)[q]   = i*(n-1)+j;
          (*vv)[q++] = b[j+n+1]*h[i];
      }
      
  	  for ( j=0; j<n-1; j++ ) {/* i+1 */
          (*jv)[q]   = (i+1)*(n-1)+j;
          (*vv)[q++] = b[j+2*n+2]*h[i+1];
      }
      
      if ( i>0 ) {
          (*jv)[q]   = m*(n-1)+i-1;
          (*vv)[q++] = b[2*n+1]*h[i];
      }
      
      (*jv)[q]   = m*(n-1)+i;
      (*vv)[q++] = b[2*n]*h[i]+b[3*n+2]*h[i+1];
      
      if ( i<m-2 ) {
          (*jv)[q]   = m*(n-1)+i+1;
          (*vv)[q++] = b[3*n+1]*h[i+1];
      }
      
      (*iv)[++p] = q;
  }
  (*iv)--;
  printf("nnz(M) = [%d, %d]\n",q,5*m*n-2*m-4*n+1-(m>1?2:0));
  
}   


int main(int argc, char *argv[])
{
	int i, j, k, m=2, n=10;
	
	if (argc>1) m = atoi(argv[1]);
	if (argc>2) n = atoi(argv[2]);
    
    R *h = (R *) malloc(m*sizeof(R));
    
    for ( i=0; i<m; i++ ) h[i] = 2./m;
    
    for ( i=0; i<(argc-3<m?argc-3:m); i++ ) h[i] = atof(argv[i+3]);
	
	for ( i=0; i<m; i++ ) h[i] *= 0.5;
	
	R *b, *V;
	
	LegendMat(n, &b, &V);
	
	int *iu, *ju, *iv, *jv;
	R   *vu, *vv;
	
	assemble1d(m, n, h, b, &iu, &ju, &vu, &iv, &jv, &vv);
	
	FILE *fp=fopen("stiff","w+");
	for ( j=0; j<iu[0]; j++ )
	for ( k=iu[j+1]; k<iu[j+2]; k++ )
	    fprintf(fp,"%4d %4d %.16g\n",j,ju[k],vu[k]);
    fclose(fp);
	
	FILE *fd=fopen("mass","w+");
	for ( j=0; j<iv[0]; j++ ) 
	for ( k=iv[j+1]; k<iv[j+2]; k++ )
	    fprintf(fd,"%4d %4d %.16g\n",j,jv[k],vv[k]);
    fclose(fd);
    
    free(V);
    free(b);
    free(iu);
    free(ju);
    free(vu);
    free(iv);
    free(jv);
    free(vv);
}