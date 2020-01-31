#include <stdlib.h>
#include <math.h>

#define F77(a) d ## a ## _
typedef int INT;
typedef double R;

void F77(steqr)(char *, int *, R *, R *, R *, int *, R *, int *);

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
  
/*****************************************************************/
void matCcube(int n1, int n2, int n3, R h1, R h2, R h3, R *b1, R *b2, R *b3, R c, int **iw, int **jw, R **vw)
{
 
 int j1,j2,j3, k1,k2,k3, p, q, off[8];
 R u1,u2,u3, v1,v2,v3;
 
   
  off[0] = 0;
  off[1] = off[0] + 2*(n1-1)*(n2-1);
  off[2] = off[1] + 2*(n1-1)*(n3-1);
  off[3] = off[2] + 2*(n2-1)*(n3-1);
  off[4] = off[3] + 4*(n1-1);
  off[5] = off[4] + 4*(n2-1);
  off[6] = off[5] + 4*(n3-1);
  off[7] = off[6] + 8;
  
  p = off[7]-off[0];
  q = 4*(18*n1*n2*n3+7*n1*n2+7*n1*n3+7*n2*n3-12*n1-12*n2-12*n3+29);

  *iw = (int *) malloc(sizeof(int)*(p+2));
  *jw = (int *) malloc(sizeof(int)*q);
  *vw = (R   *) malloc(sizeof(R  )*q);
 
  /* number of rows
     = 2*(n1*n2+n2*n3+n3*n1+1)
     
     number of columns
     = 2*(n1*n2+n2*n3+n3*n1+1)
     
     number of nonzero entries 
     =[2*(n1-1)*(n2-1)]*[5*n3+10]
     +[2*(n1-1)*(n3-1)]*[5*n2+10]
     +[2*(n2-1)*(n3-1)]*[5*n1+10]
     +[4*(n1-1)]*[2*n2*n3+4*n2+4*n3+2]
     +[4*(n2-1)]*[2*n1*n3+4*n1+4*n3+2]
     +[4*(n3-1)]*[2*n1*n2+4*n1+4*n2+2]     
     +8*[2*n1*n2+2*n2*n3+2*n3*n1+10]
     =4*(18*n1*n2*n3+7*n1*n2+7*n1*n3+7*n2*n3-12*n1-12*n2-12*n3+29)
  */	
	
  (*iw)[0] = p;
  (*iw)[1] = p;
  
  
  (*iw)[2] = q = 0;
  p        = 2;
  
  for ( j3=n3-1; j3<n3+1; j3++ )
  for ( j1=0;    j1<n1-1; j1++ )
  for ( j2=0;    j2<n2-1; j2++ ) {
  	  
  	  /* k1=j1 and k2=j2 and k3=n-1,n */
  	  u1 = 1./h1;
  	  u2 = 1./h2;
  	  
  	  v1 = h1*b1[j1];
  	  v2 = h2*b2[j2];
  	  
  	  for ( k3=n3-1; k3<n3+1; k3++ ) {
  	  	
  	  	  u3 = k3-j3?-.5/h3:.5/h3;
  	  	  v3 = h3*b3[(j3-n3+2)*(n3+1)+k3];
  	  	  
          (*vw)[q]   = u1*v2*v3 + v1*u2*v3 + v1*v2*u3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[0] + ((k3-n3+1)*(n1-1)+j1)*(n2-1)+j2;
  	  }
  	  
  	  /* k1=j1 and k2=n-1,n and 0<=k3<=n-2 */
  	  u1 = 1./h1;
   /* u2 = 0.; */
   /* u3 = 0.; */
  	  
  	  v1 = h1*b1[j1];
  	  
  	  for ( k2=n2-1; k2<n2+1; k2++ ) 
  	  for ( k3=0;    k3<n3-1; k3++ ) {
  	  	
       	  v2 = h2*b2[(k2-n2+2)*(n2+1)+j2];
  	  	  v3 = h3*b3[(j3-n3+2)*(n3+1)+k3];
  	  	  
          (*vw)[q]   = u1*v2*v3 +  c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[1] + ((k2-n2+1)*(n1-1)+j1)*(n3-1)+k3;
  	  }
  	  
  	  /* k1=n-1,n and k2=j2 and 0<=k3<=n-2 */
   /* u1 = 0.; */
  	  u2 = 1./h2;
   /* u3 = 0.; */
  	  
  	  v2 = h2*b2[j2];
  	  
  	  for ( k1=n1-1; k1<n1+1; k1++ ) 
  	  for ( k3=0;    k3<n3-1; k3++ ) {
  	  	
       	  v1 = h1*b1[(k1-n1+2)*(n1+1)+j1];
  	  	  v3 = h3*b3[(j3-n3+2)*(n3+1)+k3];
  	  	  
          (*vw)[q]   =  v1*u2*v3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[2] + ((k1-n1+1)*(n2-1)+j2)*(n3-1)+k3;
  	  }
  	  
  	  /* k1=j1 and k2=n-1,n and k3=n-1,n */
  	  u1 = 1./h1;
   /* u2 = 0.; */
  	  
  	  v1 = h1*b1[j1];
  	  
  	  for ( k2=n2-1; k2<n2+1; k2++ ) 
  	  for ( k3=n3-1; k3<n3+1; k3++ ) {
  	  	
  	  	  u3 = k3-j3?-.5/h3:.5/h3;
  	  	  	
       	  v2 = h2*b2[(k2-n2+2)*(n2+1)+j2];
  	  	  v3 = h3*b3[(j3-n3+2)*(n3+1)+k3];
  	  	  
          (*vw)[q]   = u1*v2*v3 + v1*v2*u3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[3] + ((k2-n1+1)*2+k3-n3+1)*(n1-1)+j1;
  	  }

  	  /* k1=n-1,n and k2=j2 and k3=n-1,n */
   /* u1 = 0.; */
  	  u2 = 1./h2;
  	  
  	  v2 = h2*b2[j2];
  	  
  	  for ( k1=n1-1; k1<n1+1; k1++ ) 
  	  for ( k3=n3-1; k3<n3+1; k3++ ) {
  	  	
  	  	  u3 = k3-j3?-.5/h3:.5/h3;
  	  	  	
       	  v1 = h1*b1[(k1-n1+2)*(n1+1)+j1];
  	  	  v3 = h3*b3[(j3-n3+2)*(n3+1)+k3];
  	  	  
          (*vw)[q]   = v1*u2*v3 + v1*v2*u3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[4] + ((k1-n1+1)*2+k3-n3+1)*(n2-1)+j2;
  	  }
  	  
  	  /* k1=n-1,n and k2=n-1,n and 0<=k3<=n-2 */
   /* u1 = 0.; */
   /* u2 = 0.; */
   /* u3 = 0.; */
  	  
  	  for ( k1=n1-1; k1<n1+1; k1++ ) 
  	  for ( k2=n2-1; k2<n2+1; k2++ ) 
  	  for ( k3=0;    k3<n3-1; k3++ ) {
  	  	
       	  v1 = h1*b1[(k1-n1+2)*(n1+1)+j1];
       	  v2 = h2*b2[(k2-n2+2)*(n2+1)+j2];
  	  	  v3 = h3*b3[(j3-n3+2)*(n3+1)+k3];
  	  	  
          (*vw)[q]   = c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[5] + ((k1-n1+1)*2+j2-n2+1)*(n3-1)+k3;
  	  }
  	  
  	  /* k1=n-1,n and k2=n-1,n and k3=n-1,n */
   /* u1 = 0.; */
   /* u2 = 0.; */
  	  
  	  for ( k1=n1-1; k1<n1+1; k1++ ) 
  	  for ( k2=n2-1; k2<n2+1; k2++ ) 
  	  for ( k3=n3-1; k3<n3+1; k3++ ) {
  	  	  
  	  	  u3 = k3-j3?-.5/h3:.5/h3;
  	  	
       	  v1 = h1*b1[(k1-n1+2)*(n1+1)+j1];
       	  v2 = h2*b2[(k2-n2+2)*(n2+1)+j2];
  	  	  v3 = h3*b3[(j3-n3+2)*(n3+1)+k3];
  	  	  
          (*vw)[q]   = v1*v2*u3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[6] + ((k1-n1+1)*2+k2-n2+1)*2+k3-n3+1;
  	  }
  	
  	  (*iw)[++p] = q;
  }
  
  /************************/
  for ( j2=n2-1; j2<n2+1; j2++ )
  for ( j1=0;    j1<n1-1; j1++ )
  for ( j3=0;    j3<n3-1; j3++ ) {
  	  
  	  /* k1=j1 and 0<=k2<=n-2 and k3=n-1,n */
  	  u1 = 1./h1;
  	//u2 = 0.;
  	//u3 = 0.;
  	  
  	  v1 = h1*b1[j1];
  	  
  	  for ( k3=n3-1; k3<n3+1; k3++ ) 
  	  for ( k2=0;    k2<n2-1; k2++ ) {
  	  	
  	  	  v2      = h2*b2[(j2-n2+2)*(n2+1)+k2];
  	  	  v3      = h3*b3[(k3-n3+2)*(n3+1)+j3];
  	  	  
          (*vw)[q]   = u1*v2*v3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[0] + ((k3-n3+1)*(n1-1)+j1)*(n2-1)+k2;
  	  }
  	  
  	  /* k1=j1 and k2=n-1,n and k3=j3 */
  	  u1 = 1./h1;
  	  u3 = 1./h3;
  	  
  	  v1 = h1*b1[j1];
  	  v3 = h3*b3[j3];
  	  
  	  for ( k2=n2-1; k2<n2+1; k2++ ) {
  	  	
  	  	  u2 = k2-j2?-.5/h2:.5/h2;
  	  	  	
       	  v2 = h2*b2[(k2-n2+2)*(n2+1)+j2];
  	  	  
          (*vw)[q]   = u1*v2*v3 + v1*u2*v3 + v1*v2*u3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[1] + ((k2-n2+1)*(n1-1)+j1)*(n3-1)+j3;
  	  }
  	  
  	  /* k1=n-1,n and 0<=k2<=n-2 and k3=j3 */
  	//u1 = 0.;
  	//u2 = 0.;
  	  u3 = 1./h3;
  	  
  	  v3 = h3*b3[j3];
  	  
  	  for ( k1=n1-1; k1<n1+1; k1++ ) 
  	  for ( k2=0;    k2<n2-1; k2++ ) {
  	  	
       	  v1 = h1*b1[(k1-n1+2)*(n1+1)+j1];
  	  	  v2 = h2*b2[(j2-n2+2)*(n2+1)+k2];
  	  	  
          (*vw)[q]   = v1*v2*u3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[2] + ((k1-n1+1)*(n2-1)+k2)*(n3-1)+j3;
  	  }
  	  
  	  /* k1=j1 and k2=n-1,n and k3=n-1,n */
  	  u1 = 1./h1;
  	//u3 = 0.;
  	  
  	  v1 = h1*b1[j1];
  	  
  	  for ( k2=n2-1; k2<n2+1; k2++ ) 
  	  for ( k3=n3-1; k3<n3+1; k3++ ) {
  	  	
  	  	  u2 = k2-j2?-.5/h2:.5/h2;
  	  	  	
       	  v2 = h2*b2[(k2-n2+2)*(n2+1)+j2];
  	  	  v3 = h3*b3[(k3-n3+2)*(n3+1)+j3];
  	  	  
          (*vw)[q]   = u1*v2*v3 + v1*u2*v3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[3] + ((k2-n2+1)*2+k3-n3+1)*(n1-1)+j1;
  	  }
  	  
  	  /* k1=n-1,n and 0<=k2<=n-2 and k3=n-1,n */
  	  //u1 = 0.;
  	  //u2 = 0.;
  	  //u3 = 0.;
  	  
  	  for ( k1=n1-1; k1<n1+1; k1++ ) 
  	  for ( k3=n3-1; k3<n3+1; k3++ ) 
  	  for ( k2=0;    k2<n2-1; k2++ ) {
  	  	
       	  v1 = h1*b1[(k1-n1+2)*(n1+1)+j1];
       	  v2 = h2*b2[(k2-n2+2)*(n2+1)+j2];
  	  	  v3 = h3*b3[(k3-n3+2)*(n3+1)+j3];
  	  	  
          (*vw)[q]   = c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[4] + ((k1-n1+1)*2+k3-n3+1)*(n2-1)+k2;
  	  }
  	  
  	  /* k1=n-1,n and k2=n-1,n and k3=j3 */
    //u1 = 0;
  	  u3 = 1./h3;
  	  
  	  v3 = h3*b3[j3];
  	  
  	  for ( k1=n1-1; k1<n1+1; k1++ ) 
  	  for ( k2=n2-1; k2<n2+1; k2++ ) {
  	  	
  	  	  u2 = k2-j2?-.5/h2:.5/h2;
  	  	  
       	  v1 = h1*b1[(k1-n1+2)*(n1+1)+j1];
       	  v2 = h2*b2[(k2-n2+2)*(n2+1)+j2];
  	  	  
          (*vw)[q]   = v1*u2*v3 + v1*v2*u3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[5] + ((k1-n1+1)*2+k2-n2+1)*(n3-1)+j3;
  	  }
  	  
  	  /* k1=n-1,n and k2=n-1,n and k3=n-1,n */
  	  //u1 = 0;
  	  //u3 = 0;
  	  
  	  for ( k1=n1-1; k1<n1+1; k1++ ) 
  	  for ( k2=n2-1; k2<n2+1; k2++ ) 
  	  for ( k3=n3-1; k3<n3+1; k3++ ) {
  	  	  
  	  	  u2 = k2-j2?-.5/h2:.5/h2;
  	  	
       	  v1 = h1*b1[(k1-n1+2)*(n1+1)+j1];
       	  v2 = h2*b2[(k2-n2+2)*(n2+1)+j2];
  	  	  v3 = h3*b3[(k3-n3+2)*(n3+1)+j3];
  	  	  
          (*vw)[q]   = v1*u2*v3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[6] + ((k1-n1+1)*2+k2-n2+1)*2+k3-n3+1;
  	  }
  	
  	  (*iw)[++p] = q;
  }


  /************************/
  for ( j1=n1-1; j1<n1+1; j1++ )
  for ( j2=0;    j2<n2-1; j2++ )
  for ( j3=0;    j3<n3-1; j3++ ) {
  	  
  	  /* 0<=k1=n-2 and k2=j2 and k3=n-1,n */
  	//u1 = 0.;
  	  u2 = 1./h2;
  	//u3 = 0.;
  	  
  	  v2 = h2*b2[j1];
  	  
  	  for ( k3=n3-1; k3<n3+1; k3++ ) 
  	  for ( k1=0;    k1<n1-1; k1++ ) {
  	  	
  	  	  v1      = h1*b1[(j1-n1+2)*(n1+1)+k1];
  	  	  v3      = h3*b3[(k3-n3+2)*(n3+1)+j3];
  	  	  
          (*vw)[q]   = v1*u2*v3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[0] + ((k3-n3+1)*(n1-1)+k1)*(n2-1)+j2;
  	  }
  	  
  	  /* 0<=k1<=n-2 and k2=n-1,n and k3=j3 */
  	//u1 = 0.;
  	//u2 = 0.;
  	  u3 = 1./h3;
  	  
  	  v3 = h3*b3[j3];
  	  
  	  for ( k2=n2-1; k2<n2+1; k2++ ) 
  	  for ( k1=0;    k1<n1-1; k1++ ) {
  	  	
       	  v1 = h1*b1[(j1-n1+2)*(n1+1)+k1];
       	  v2 = h2*b2[(k2-n2+2)*(n2+1)+j2];
  	  	  
          (*vw)[q]   = v1*v2*u3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[1] + ((k2-n2+1)*(n1-1)+k1)*(n3-1)+j3;
  	  }
  	  
  	  /* k1=n-1,n and k2=j2 and k3=j3 */
  	  u2 = 1./h2;
  	  u3 = 1./h3;
  	  
  	  v2 = h2*b2[j2];
  	  v3 = h3*b3[j3];
  	  
  	  for ( k1=n1-1; k1<n1+1; k1++ ) {
  	  	
  	  	  u1 = k1-j1?-.5/h1:.5/h1;
       	  v1 = h1*b1[(k1-n1+2)*(n1+1)+j1];
  	  	  
          (*vw)[q]   = u1*v2*v3 + v1*u2*v3 + v1*v2*u3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[2] + ((k1-n1+1)*(n2-1)+j2)*(n3-1)+j3;
  	  }
  	  
  	  /* 0<=k1<n-2 and k2=n-1,n and k3=n-1,n */
  	  //u1 = 0.;
  	  //u2 = 0.;
  	  //u3 = 0.;
  	  
  	  for ( k2=n2-1; k2<n2+1; k2++ ) 
  	  for ( k3=n3-1; k3<n3+1; k3++ ) 
  	  for ( k1=0;    k1<n1-1; k1++ ) {
  	  	
       	  v1 = h1*b1[(j1-n1+2)*(n1+1)+k1];
       	  v2 = h2*b2[(k2-n2+2)*(n2+1)+j2];
  	  	  v3 = h3*b3[(k3-n3+2)*(n3+1)+j3];
  	  	  
          (*vw)[q]   = c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[3] + ((k2-n2+1)*2+k3-n3+1)*(n1-1)+k1;
  	  }
  	  
  	  /* k1=n-1,n and k2=j2 and k3=n-1,n */
  	  u2 = 1./h2;
  	//u3 = 0.;
  	  
  	  v2 = h2*b2[j2];
  	  
  	  for ( k1=n1-1; k1<n1+1; k1++ ) 
  	  for ( k3=n3-1; k3<n3+1; k3++ ) {
  	  	
  	  	  u1 = k1-j1?-.5/h1:.5/h1;
  	  	  	
       	  v1 = h1*b1[(k1-n1+2)*(n1+1)+j1];
  	  	  v3 = h3*b3[(k3-n3+2)*(n3+1)+j3];
  	  	  
          (*vw)[q]   = u1*v2*v3 + v1*u2*v3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[4] + ((k1-n1+1)*2+k3-n3+1)*(n2-1)+j2;
  	  }
  	  
  	  /* k1=n-1,n and k2=n-1,n and k3=j3 */
  	//u2 = 0.;
  	  u3 = 1./h3;
  	  
  	  v3 = h3*b3[j3];
  	  
  	  for ( k1=n1-1; k1<n1+1; k1++ ) 
  	  for ( k2=n2-1; k2<n2+1; k2++ ) {
  	  	  
  	  	  u1 = k1-j1?-.5/h1:.5/h1;
       	  v1 = h1*b1[(k1-n1+2)*(n1+1)+j1];
       	  v2 = h2*b2[(k2-n2+2)*(n2+1)+j2];
  	  	  
          (*vw)[q]   = u1*v2*v3 + v1*v2*u3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[5] + ((k1-n1+1)*2+k2-n2+1)*(n3-1)+j3;
  	  }
  	  
  	  /* k1=n-1,n and k2=n-1,n and k3=n-1,n */
  	//u2 = 0.;
  	//u3 = 0.;
  	  
  	  for ( k1=n1-1; k1<n1+1; k1++ ) 
  	  for ( k2=n2-1; k2<n2+1; k2++ ) 
  	  for ( k3=n3-1; k3<n3+1; k3++ ) {
  	  	  
  	  	  u1 = k1-j1?-.5/h1:.5/h1;
  	  	
       	  v1 = h1*b1[(k1-n1+2)*(n1+1)+j1];
       	  v2 = h2*b2[(k2-n2+2)*(n2+1)+j2];
  	  	  v3 = h3*b3[(k3-n3+2)*(n3+1)+j3];
  	  	  
          (*vw)[q]   = u1*v2*v3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[6] + ((k1-n1+1)*2+k2-n2+1)*2+k3-n3+1;
  	  }
  	  
  	  (*iw)[++p] = q;
  }
	
  /************************/
  for ( j2=n2-1; j2<n2+1; j2++ )
  for ( j3=n3-1; j3<n3+1; j3++ ) 
  for ( j1=0;    j1<n1-1; j1++ ) {
  	  
  	  /* k1=j1 and 0<=k2<=n-2 and k3=n-1,n */
  	  u1 = 1./h1;
  	//u2 = 0.;
  	  
  	  v1 = h1*b1[j1];
  	  
  	  for ( k3=n3-1; k3<n3+1; k3++ ) 
  	  for ( k2=0;    k2<n2-1; k2++ ) {
  	  	
  	  	  u3      = k3-j3?-.5/h3:.5/h3;
  	  	  	
  	  	  v2      = h2*b2[(j2-n2+2)*(n2+1)+k2];
  	  	  v3      = h3*b3[(k3-n3+2)*(n3+1)+j3];
  	  	  
          (*vw)[q]   = u1*v2*v3 + v1*v2*u3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[0] + ((k3-n3+1)*(n1-1)+j1)*(n2-1)+k2;
  	  }
  	  
  	  /* k1=j1 and k2=n-1,n and 0<=k3<=n-2 */
  	  u1 = 1./h1;
  	//u3 = 0.;
  	  
  	  v1 = h1*b1[j1];
  	  
  	  for ( k2=n2-1; k2<n2+1; k2++ ) 
  	  for ( k3=0;    k3<n3-1; k3++ ) {
  	  	
  	  	  u2 = k2-j2?-.5/h2:.5/h2;
  	  	  	
       	  v2 = h2*b2[(k2-n2+2)*(n2+1)+j2];
       	  v3 = h3*b3[(j3-n3+2)*(n3+1)+k3];
  	  	  
          (*vw)[q]   = u1*v2*v3 + v1*u2*v3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[1] + ((k2-n2+1)*(n1-1)+j1)*(n3-1)+k3;
  	  }
  	  
  	  /* k1=n-1,n and 0<=k2<=n-2 and 0<=k3<=n-2 */
  	//u1 = 0.;
  	//u2 = 0.;
  	//u3 = 0.;
  	  
  	  for ( k1=n1-1; k1<n1+1; k1++ ) 
  	  for ( k2=0;    k2<n2-1; k2++ ) 
  	  for ( k3=0;    k3<n3-1; k3++ ) 
  	  {
  	  	
       	  v1 = h1*b1[(k1-n1+2)*(n1+1)+j1];
       	  v2 = h2*b2[(j2-n2+2)*(n2+1)+k2];
       	  v3 = h3*b3[(j3-n3+2)*(n3+1)+k3];
  	  	  
          (*vw)[q]   = c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[2] + ((k1-n1+1)*(n2-1)+k2)*(n3-1)+k3;
  	  }

  	  /* k1=j1 and k2=n-1,n and k3=n-1,n */
  	  u1 = 1./h1;
  	  
  	  for ( k2=n2-1; k2<n2+1; k2++ ) 
  	  for ( k3=n3-1; k3<n3+1; k3++ ) {
  	  	
      	  u2 = k2-j2?-.5/h2:.5/h2;
  	      u3 = k3-j3?-.5/h3:.5/h3;
  	      	
       	  v2 = h2*b2[(k2-n2+2)*(n2+1)+j2];
  	  	  v3 = h3*b3[(k3-n3+2)*(n3+1)+j3];
  	  	  
          (*vw)[q]   = u1*v2*v3 + v1*u2*v3 + v1*v2*u3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[3] + ((k2-n2+1)*2+k3-n3+1)*(n1-1)+j1;
  	  }
  	  
  	  /* k1=n-1,n and 0<=k2<=n-2 and k3=n-1,n */
  	//u1 = 0.;
  	//u2 = 0.;
  	  
  	  for ( k1=n1-1; k1<n1+1; k1++ ) 
  	  for ( k3=n3-1; k3<n3+1; k3++ ) 
  	  for ( k2=0;    k2<n2-1; k2++ ) {
  	  	
  	  	  u3 = k3-j3?-.5/h3:.5/h3;
  	  	  	
       	  v1 = h1*b1[(k1-n1+2)*(n1+1)+j1];
       	  v2 = h2*b2[(j2-n2+2)*(n2+1)+k2];
  	  	  v3 = h3*b3[(k3-n3+2)*(n3+1)+j3];
  	  	  
          (*vw)[q]   = v1*v2*u3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[4] + ((k1-n1+1)*2+k3-n3+1)*(n2-1)+k2;
  	  }
  	  
  	  
  	  /* k1=n-1,n and k2=n-1,n and 0<=k3<=n-2 */
  	//u1 = 0;
  	//u3 = 0.;
  	  
  	  for ( k1=n1-1; k1<n1+1; k1++ ) 
  	  for ( k2=n2-1; k2<n2+1; k2++ ) 
  	  for ( k3=0;    k3<n3-1; k3++ ) {
  	  	  
  	  	  u2 = k2-j2?-.5/h2:.5/h2;
  	  	  	
       	  v1 = h1*b1[(k1-n1+2)*(n1+1)+j1];
       	  v2 = h2*b2[(k2-n2+2)*(n2+1)+j2];
       	  v3 = h3*b3[(j3-n3+2)*(n3+1)+k3];
  	  	  
          (*vw)[q]   = v1*u2*v3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[5] + ((k1-n1+1)*2+k2-n2+1)*(n3-1)+k3;
  	  }
  	  
  	  /* k1=n-1,n and k2=n-1,n and k3=n-1,n */
  	 //u1 = 0;
  	  
  	  for ( k1=n1-1; k1<n1+1; k1++ ) 
  	  for ( k2=n2-1; k2<n2+1; k2++ ) 
  	  for ( k3=n3-1; k3<n3+1; k3++ ) {
  	  	  
  	  	  u2 = k2-j2?-.5/h2:.5/h2;
  	  	  u3 = k3-j3?-.5/h3:.5/h3;
  	  	
       	  v1 = h1*b1[(k1-n1+2)*(n1+1)+j1];
       	  v2 = h2*b2[(k2-n2+2)*(n2+1)+j2];
  	  	  v3 = h3*b3[(k3-n3+2)*(n3+1)+j3];
  	  	  
          (*vw)[q]   = v1*u2*v3 + v1*v2*u3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[6] + ((k1-n1+1)*2+k2-n2+1)*2+k3-n3+1;
  	  }
  	
  	  (*iw)[++p] = q;
  }
  
  /************************/
  for ( j1=n1-1; j1<n1+1; j1++ )
  for ( j3=n3-1; j3<n3+1; j3++ ) 
  for ( j2=0;    j2<n2-1; j2++ ) {
  	  
  	  /* 0<=k1<=n-2 and k2=j2 and k3=n-1,n */
  	//u1 = 0.;
  	  u2 = 1./h2;
  	  
  	  v2 = h2*b2[j2];
  	  
  	  for ( k3=n3-1; k3<n3+1; k3++ ) 
  	  for ( k1=0;    k1<n1-1; k1++ ) {
  	  	
  	  	  u3      = k3-j3?-.5/h3:.5/h3;
  	  	  v1      = h1*b1[(j1-n1+2)*(n1+1)+k1];
  	  	  v3      = h3*b3[(k3-n3+2)*(n3+1)+j3];
  	  	  
          (*vw)[q]   = v1*u2*v3 + v1*v2*u3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[0] + ((k3-n3+1)*(n1-1)+k1)*(n2-1)+j2;
  	  }
  	  
  	  /* 0<=k1<=n-2 and k2=n-1,n and 0<=k3<=n-2 */
  	//u1 = 0.;
  	//u2 = 0.;
  	//u3 = 0.;
  	  
  	  for ( k2=n2-1; k2<n2+1; k2++ )
  	  for ( k1=0;    k1<n1-1; k1++ )  
  	  for ( k3=0;    k3<n3-1; k3++ ) {
  	  	  	  	 
       	  v1 = h1*b1[(j1-n1+2)*(n1+1)+k1];
       	  v2 = h2*b2[(k2-n2+2)*(n2+1)+j2];
       	  v3 = h3*b3[(j3-n3+2)*(n3+1)+k3];
  	  	  
          (*vw)[q]   = c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[1] + ((k2-n2+1)*(n1-1)+k1)*(n3-1)+k3;
  	  }
  	  
  	  /* k1=n-1,n and k2=j2 and 0<=k3<=n-2 */
  	  u2 = 1./h2;
  	//u3 = 0.;
  	  
  	  v2 = h2*b2[j2];
  	  
  	  for ( k1=n1-1; k1<n1+1; k1++ ) 
  	  for ( k3=0;    k3<n3-1; k3++ ) {
  	  	
  	  	  u1 = k1-j1?-.5/h1:.5/h1;
       	  v1 = h1*b1[(j1-n1+2)*(n1+1)+k1];
       	  v3 = h3*b3[(j3-n3+2)*(n3+1)+k3];
  	  	  
          (*vw)[q]   = u1*v2*v3 + v1*u2*v3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[2] + ((k1-n1+1)*(n2-1)+j2)*(n3-1)+k3;
  	  }
  	  
  	  /* 0<=k1<=n-2 and k2=n-1,n and k3=n-1,n */
  	//u1 = 0.;
  	//u2 = 0.;
  	  
  	  for ( k2=n2-1; k2<n2+1; k2++ ) 
  	  for ( k3=n3-1; k3<n3+1; k3++ ) 
  	  for ( k1=0;    k1<n1-1; k1++ ) {
  	  	
  	      u3 = k3-j3?-.5/h3:.5/h3;
  	      	
       	  v1 = h1*b1[(j1-n1+2)*(n1+1)+k1];
       	  v2 = h2*b2[(k2-n2+2)*(n2+1)+j2];
  	  	  v3 = h3*b3[(k3-n3+2)*(n3+1)+j3];
  	  	  
          (*vw)[q]   = v1*v2*u3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[3] + ((k2-n2+1)*2+k3-n3+1)*(n1-1)+k1;
  	  }

  	  /* k1=n-1,n and k2=j2 and k3=n-1,n */
  	  u2 = 1./h2;
  	  
  	  v2 = h2*b2[j2];
  	  
  	  for ( k1=n1-1; k1<n1+1; k1++ ) 
  	  for ( k3=n3-1; k3<n3+1; k3++ ) {
  	  	
  	  	  u1 = k1-j1?-.5/h1:.5/h1;
  	  	  u3 = k3-j3?-.5/h3:.5/h3;
  	  	  	
       	  v1 = h1*b1[(k1-n1+2)*(n1+1)+j1];
  	  	  v3 = h3*b3[(k3-n3+2)*(n3+1)+j3];
  	  	  
          (*vw)[q]   = u1*v2*v3 + v1*u2*v3 + v1*v2*u3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[4] + ((k1-n1+1)*2+k3-n3+1)*(n2-1)+j2;
  	  }
  	  
  	  /* k1=n-1,n and k2=n-1,n and 0<=k3<=n-2 */
  	//u2 = 0.;
  	//u3 = 0.;
  	  
  	  for ( k1=n1-1; k1<n1+1; k1++ ) 
  	  for ( k2=n2-1; k2<n2+1; k2++ ) 
  	  for ( k3=0;    k3<n3-1; k3++ ) {
  	  	  
  	  	  u1 = k1-j1?-.5/h1:.5/h1;
  	  	  	
       	  v1 = h1*b1[(k1-n1+2)*(n1+1)+j1];
       	  v2 = h2*b2[(k2-n2+2)*(n2+1)+j2];
       	  v3 = h3*b3[(j3-n3+2)*(n3+1)+k3];
  	  	  
          (*vw)[q]   = u1*v2*v3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[5] + ((k1-n1+1)*2+k2-n2+1)*(n3-1)+k3;
  	  }
  	  
  	  /* k1=n-1,n and k2=n-1,n and k3=n-1,n */
  	//u2 = 0.;
  	  
  	  for ( k1=n1-1; k1<n1+1; k1++ ) 
  	  for ( k2=n2-1; k2<n2+1; k2++ ) 
  	  for ( k3=n3-1; k3<n3+1; k3++ ) {
  	  	  
  	  	  u1 = k1-j1?-.5/h1:.5/h1;
  	  	  u3 = k3-j3?-.5/h3:.5/h3;
  	  	
       	  v1 = h1*b1[(k1-n1+2)*(n1+1)+j1];
       	  v2 = h2*b2[(k2-n2+2)*(n2+1)+j2];
  	  	  v3 = h3*b3[(k3-n3+2)*(n3+1)+j3];
  	  	  
          (*vw)[q]   = u1*v2*v3 + v1*v2*u3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[6] + ((k1-n1+1)*2+k2-n2+1)*2+k3-n3+1;
  	  }
  	  
  	  (*iw)[++p] = q;  	
  }
  
  /************************/
  for ( j1=n1-1; j1<n1+1; j1++ )
  for ( j2=n2-1; j2<n2+1; j2++ ) 
  for ( j3=0;    j3<n3-1; j3++ ) {
  	  
  	  /* 0<=k1<=n-2 and 0<=k2<=n-2 and k3=n-1,n */
  	//u1 = 0.;
  	//u2 = 0.;
  	//u3 = 0.;
  	  
  	  for ( k3=n3-1; k3<n3+1; k3++ ) 
  	  for ( k1=0;    k1<n1-1; k1++ )
  	  for ( k2=0;    k2<n2-1; k2++ ) {
  	  	
  	  	  v1      = h1*b1[(j1-n1+2)*(n1+1)+k1];
  	  	  v2      = h2*b2[(j2-n2+2)*(n2+1)+k2];
  	  	  v3      = h3*b3[(k3-n3+2)*(n3+1)+j3];
  	  	  
          (*vw)[q]   = c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[0] + ((k3-n3+1)*(n1-1)+k1)*(n2-1)+k2;
  	  }

  	  /* 0<=k1<=n-2 and k2=n-1,n and k3=j3 */
  	//u1 = 0.;
  	  u3 = 1./h3;
  	  
  	  v3 = h3*b3[j3];
  	  
  	  for ( k2=n2-1; k2<n2+1; k2++ )
  	  for ( k1=0;    k1<n1-1; k1++ ) {
  	  	  	  	 
  	  	  u2 = k2-j2?-.5/h2:.5/h2;
  	  	    	 
       	  v1 = h1*b1[(j1-n1+2)*(n1+1)+k1];
       	  v2 = h2*b2[(k2-n2+2)*(n2+1)+j2];
  	  	  
          (*vw)[q]   = v1*u2*v3 + v1*v2*u3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[1] + ((k2-n2+1)*(n1-1)+k1)*(n3-1)+j3;
  	  }
  	  
  	  /* k1=n-1,n and 0<=k2<n-2 and k3=j3 */
  	//u2 = 0.;
  	  u3 = 1./h3;
  	  
  	  v3 = h3*b3[j3];
  	  
  	  for ( k1=n1-1; k1<n1+1; k1++ ) 
  	  for ( k2=0;    k2<n2-1; k2++ ) {
  	  	
  	  	  u1 = k1-j1?-.5/h1:.5/h1;
  	  	  	
       	  v1 = h1*b1[(j1-n1+2)*(n1+1)+k1];
       	  v2 = h2*b2[(j2-n2+2)*(n2+1)+k2];
  	  	  
          (*vw)[q]   = u1*v2*v3 + v1*v2*u3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[2] + ((k1-n1+1)*(n2-1)+k2)*(n3-1)+j3;
  	  }
  	  
  	  /* 0<=k1<=n-2 and k2=n-1,n and k3=n-1,n */
  	//u1 = 0.;
  	//u3 = 0.;
  	  
  	  for ( k2=n2-1; k2<n2+1; k2++ ) 
  	  for ( k3=n3-1; k3<n3+1; k3++ ) 
  	  for ( k1=0;    k1<n1-1; k1++ ) {
  	  	
  	      u2 = k2-j2?-.5/h2:.5/h2;
  	      	
       	  v1 = h1*b1[(j1-n1+2)*(n1+1)+k1];
       	  v2 = h2*b2[(k2-n2+2)*(n2+1)+j2];
  	  	  v3 = h3*b3[(k3-n3+2)*(n3+1)+j3];
  	  	  
          (*vw)[q]   = v1*u2*v3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[3] + ((k2-n2+1)*2+k3-n3+1)*(n1-1)+k1;
  	  }

  	  /* k1=n-1,n and 0<=k2<=n-2 and k3=n-1,n */
  	//u2 = 0.;
  	//u3 = 0.;
  	  
  	  for ( k1=n1-1; k1<n1+1; k1++ )  
  	  for ( k3=n3-1; k3<n3+1; k3++ ) 
  	  for ( k2=0;    k2<n2-1; k2++ ) {
  	  	
  	  	  u1 = k1-j1?-.5/h1:.5/h1;
  	  	  	
       	  v1 = h1*b1[(k1-n1+2)*(n1+1)+j1];
  	  	  v2 = h2*b2[(j2-n2+2)*(n2+1)+k2];
  	  	  v3 = h3*b3[(k3-n3+2)*(n3+1)+j3];
  	  	  
          (*vw)[q]   = u1*v2*v3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[4] + ((k1-n1+1)*2+k3-n3+1)*(n2-1)+k2;
  	  }
  	  
  	  /* k1=n-1,n and k2=n-1,n and k3=j3 */
  	  u3 = 1./h3;
  	  
  	  v3 = h3*b3[j3];
  	  
  	  for ( k1=n1-1; k1<n1+1; k1++ ) 
  	  for ( k2=n2-1; k2<n2+1; k2++ ) {
  	  	  
  	  	  u1 = k1-j1?-.5/h1:.5/h1;
  	  	  u2 = k2-j2?-.5/h2:.5/h2;
  	  	  	
       	  v1 = h1*b1[(k1-n1+2)*(n1+1)+j1];
       	  v2 = h2*b2[(k2-n2+2)*(n2+1)+j2];
  	  	  
          (*vw)[q]   = u1*v2*v3 + v1*u2*v3 + v1*v2*u3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[5] + ((k1-n1+1)*2+k2-n2+1)*(n3-1)+j3;
  	  }

  	  /* k1=n-1,n and k2=n-1,n and k3=n-1,n */
  	//u3 = 0.;
  	  
  	  for ( k1=n1-1; k1<n1+1; k1++ ) 
  	  for ( k2=n2-1; k2<n2+1; k2++ ) 
  	  for ( k3=n3-1; k3<n3+1; k3++ ) {
  	  	  
  	  	  u1 = k1-j1?-.5/h1:.5/h1;
  	  	  u2 = k2-j2?-.5/h2:.5/h2;
  	  	
       	  v1 = h1*b1[(k1-n1+2)*(n1+1)+j1];
       	  v2 = h2*b2[(k2-n2+2)*(n2+1)+j2];
  	  	  v3 = h3*b3[(k3-n3+2)*(n3+1)+j3];
  	  	  
          (*vw)[q]   = u1*v2*v3 + v1*u2*v3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[6] + ((k1-n1+1)*2+k2-n2+1)*2+k3-n3+1;
  	  }
  	
  	  (*iw)[++p] = q;
  }	

  /************************/
  for ( j1=n1-1; j1<n1+1; j1++ )
  for ( j2=n2-1; j2<n2+1; j2++ ) 
  for ( j3=n3-1; j3<n3+1; j3++ ) {
  	  
  	  /* 0<=k1<=n-2 and 0<=k2<=n-2 and k3=n-1,n */
  	//u1 = 0.;
  	//u2 = 0.;
  	  
  	  for ( k3=n3-1; k3<n3+1; k3++ ) 
  	  for ( k1=0;    k1<n1-1; k1++ )
  	  for ( k2=0;    k2<n2-1; k2++ ) {
  	  	
  	  	  u3      = k3-j3?-.5/h3:.5/h3;
  	  	  
  	  	  v1      = h1*b1[(j1-n1+2)*(n1+1)+k1];
  	  	  v2      = h2*b2[(j2-n2+2)*(n2+1)+k2];
  	  	  v3      = h3*b3[(k3-n3+2)*(n3+1)+j3];
  	  	  
          (*vw)[q]   = v1*v2*u3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[0] + ((k3-n3+1)*(n1-1)+k1)*(n2-1)+k2;
  	  }

  	  /* 0<=k1<=n-2 and k2=n-1,n and 0<=k3<=n-2 */
  	//u1 = 0.;
  	//u3 = 0.;
  	  
  	  for ( k2=n2-1; k2<n2+1; k2++ )
  	  for ( k1=0;    k1<n1-1; k1++ ) 
  	  for ( k3=0;    k3<n3-1; k3++ ) {
  	  	  	  	 
  	  	  u2 = k2-j2?-.5/h2:.5/h2;
  	  	    	 
       	  v1 = h1*b1[(j1-n1+2)*(n1+1)+k1];
       	  v2 = h2*b2[(k2-n2+2)*(n2+1)+j2];
       	  v3 = h3*b3[(j3-n3+2)*(n3+1)+k3];
  	  	  
          (*vw)[q]   = v1*u2*v3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[1] + ((k2-n2+1)*(n1-1)+k1)*(n3-1)+k3;
  	  }
  	  
  	  /* k1=n-1,n and 0<=k2<n-2 and 0<=k3<=n-2 */
  	//u2 = 0.;
  	//u3 = 0.;
  	  
  	  for ( k1=n1-1; k1<n1+1; k1++ ) 
  	  for ( k2=0;    k2<n2-1; k2++ ) 
  	  for ( k3=0;    k3<n3-1; k3++ ) {
  	  	
  	  	  u1 = k1-j1?-.5/h1:.5/h1;
  	  	  	
       	  v1 = h1*b1[(j1-n1+2)*(n1+1)+k1];
       	  v2 = h2*b2[(j2-n2+2)*(n2+1)+k2];
       	  v3 = h3*b3[(j3-n3+2)*(n3+1)+k3];
  	  	  
          (*vw)[q]   = u1*v2*v3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[2] + ((k1-n1+1)*(n2-1)+k2)*(n3-1)+k3;
  	  }
  	  
  	  /* 0<=k1<=n-2 and k2=n-1,n and k3=n-1,n */
  	//u1 = 0.;
  	  
  	  for ( k2=n2-1; k2<n2+1; k2++ ) 
  	  for ( k3=n3-1; k3<n3+1; k3++ ) 
  	  for ( k1=0;    k1<n1-1; k1++ ) {
  	  	
  	      u2 = k2-j2?-.5/h2:.5/h2;
  	      u3 = k3-j3?-.5/h3:.5/h3;
  	      	
       	  v1 = h1*b1[(j1-n1+2)*(n1+1)+k1];
       	  v2 = h2*b2[(j2-n2+2)*(n2+1)+k2];
  	  	  v3 = h3*b3[(j3-n3+2)*(n3+1)+k3];
  	  	  
          (*vw)[q]   = v1*u2*v3 + v1*v2*u3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[3] + ((k2-n2+1)*2+k3-n3+1)*(n1-1)+k1;
  	  }

  	  /* k1=n-1,n and 0<=k2<=n-2 and k3=n-1,n */
  	//u2 = 0.;
  	  
  	  for ( k1=n1-1; k1<n1+1; k1++ )  
  	  for ( k3=n3-1; k3<n3+1; k3++ ) 
  	  for ( k2=0;    k2<n2-1; k2++ ) {
  	  	
  	  	  u1 = k1-j1?-.5/h1:.5/h1;
  	  	  u3 = k3-j3?-.5/h3:.5/h3;
  	  	  	
       	  v1 = h1*b1[(j1-n1+2)*(n1+1)+k1];
  	  	  v2 = h2*b2[(j2-n2+2)*(n2+1)+k2];
  	  	  v3 = h3*b3[(j3-n3+2)*(n3+1)+k3];
  	  	  
          (*vw)[q]   = u1*v2*v3 + v1*v2*u3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[4] + ((k1-n1+1)*2+k3-n3+1)*(n2-1)+k2;
  	  }
  	  
  	  /* k1=n-1,n and k2=n-1,n and 0<=k3<=n-2 */
  	//u3 = 0.;
  	  
  	  for ( k1=n1-1; k1<n1+1; k1++ ) 
  	  for ( k2=n2-1; k2<n2+1; k2++ ) 
  	  for ( k3=0;    k3<n3-1; k3++ ) {
  	  	  
  	  	  u1 = k1-j1?-.5/h1:.5/h1;
  	  	  u2 = k2-j2?-.5/h2:.5/h2;
  	  	  	
       	  v1 = h1*b1[(j1-n1+2)*(n1+1)+k1];
       	  v2 = h2*b2[(j2-n2+2)*(n2+1)+k2];
       	  v3 = h3*b3[(j3-n3+2)*(n3+1)+k3];
  	  	  
          (*vw)[q]   = u1*v2*v3 + v1*u2*v3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[5] + ((k1-n1+1)*2+k2-n2+1)*(n3-1)+k3;
  	  }

  	  /* k1=n-1,n and k2=n-1,n and k3=n-1,n */
  	  
  	  for ( k1=n1-1; k1<n1+1; k1++ ) 
  	  for ( k2=n2-1; k2<n2+1; k2++ ) 
  	  for ( k3=n3-1; k3<n3+1; k3++ ) {
  	  	  
  	  	  u1 = k1-j1?-.5/h1:.5/h1;
  	  	  u2 = k2-j2?-.5/h2:.5/h2;
  	  	  u3 = k3-j3?-.5/h3:.5/h3;
  	  	
       	  v1 = h1*b1[(k1-n1+2)*(n1+1)+j1];
       	  v2 = h2*b2[(k2-n2+2)*(n2+1)+j2];
  	  	  v3 = h3*b3[(k3-n3+2)*(n3+1)+j3];
  	  	  
          (*vw)[q]   = u1*v2*v3 + v1*u2*v3 + v1*v2*u3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jw)[q++] = off[6] + ((k1-n1+1)*2+k2-n2+1)*2+k3-n3+1;
  	  }
  	  
  	  (*iw)[++p] = q;
  }	
  
}

  
  
/*****************************************************************/
void matBcube(int n1, int n2, int n3, R h1, R h2, R h3, R *b1, R *b2, R *b3, R c, int **iz, int **jz, R **vz)
{
 
 int j1,j2,j3, k1,k2,k3, p, q;
 R u1,u2,u3, v1,v2,v3;
 
 
  *iz = (int *) malloc(sizeof(int)*2*(n1*n2+n2*n3+n3*n1+2));
  *jz = (int *) malloc(sizeof(int)*26*(n1-1)*(n2-1)*(n3-1));
  *vz = (R   *) malloc(sizeof(R  )*26*(n1-1)*(n2-1)*(n3-1));
  
  /* number of rows
     = 2*(n1-1)*(n2-1)
     
     DoF 
     =2*(n1-1)*(n2-1)*(n3-1)
     +2*(n1-1)*(n2-1)*(n3-1)
     +2*(n1-1)*(n2-1)*(n3-1)
     +4*(n1-1)*(n2-1)*(n3-1)
     +4*(n1-1)*(n2-1)*(n3-1)
     +4*(n1-1)*(n2-1)*(n3-1)     
     +8*(n1-1)*(n2-1)*(n3-1)
     =26*(n1-1)*(n2-1)*(n3-1)
  */	
	
  (*iz)[0] = 2*(n1*n2+n2*n3+n3*n1+1);
  (*iz)[1] = (n1-1)*(n2-1)*(n3-1);
  
  p = 1;
  (*iz)[2] = q = 0;
  
  for ( j3=n3-1; j3<n3+1; j3++ )
  for ( j1=0;    j1<n1-1; j1++ )
  for ( j2=0;    j2<n2-1; j2++ ) {
  	  
  	  u1 = 1./h1;
  	  u2 = 1./h2;
   /* u3 = 0.; */
  	  
  	  v1 = h1*b1[j1];
  	  v2 = h2*b2[j2];
  	  
  	  for ( k3=0;    k3<n3-1; k3++) { 
  	  	
  	      v3 = h3*b3[(j3-n3+2)*(n3+1)+k3];
  	      
          (*vz)[q]   = u1*v2*v3 + v1*u2*v3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jz)[q++] = (j1*(n2-1)+j2)*(n3-1)+k3;
      }
      
      (*iz)[++p] = q;
  }
  
  for ( j2=n2-1; j2<n2+1; j2++ )
  for ( j1=0;    j1<n1-1; j1++ )
  for ( j3=0;    j3<n3-1; j3++ ) {
  	  
  	  u1 = 1./h1;
   /* u2 = 0.; */
  	  u3 = 1./h3;
  	  
  	  v1 = h1*b1[j1];
  	  v3 = h3*b3[j3];
  	  
  	  for ( k2=0;    k2<n2-1; k2++) { 
  	  	
  	      v2 = h2*b2[(j2-n2+2)*(n2+1)+k2];
  	      
          (*vz)[q]   = u1*v2*v3 + v1*v2*u3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jz)[q++] = (j1*(n2-1)+k2)*(n3-1)+j3;
      }
      
      (*iz)[++p] = q;
  }
  
  for ( j1=n1-1; j1<n1+1; j1++ )
  for ( j2=0;    j2<n2-1; j2++ ) 
  for ( j3=0;    j3<n3-1; j3++ ) {
  	  
   /* u1 = 0.; */
  	  u2 = 1./h2;
  	  u3 = 1./h3;
  	  
  	  v2 = h2*b2[j2];
  	  v3 = h3*b3[j3];
  	  
  	  for ( k1=0;    k1<n1-1; k1++) { 
  	  	
  	      v1 = h1*b1[(j1-n1+2)*(n1+1)+k1];
  	      
          (*vz)[q]   = v1*u2*v3 + v1*v2*u3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jz)[q++] = (k1*(n2-1)+j2)*(n3-1)+j3;
      }
      
      (*iz)[++p] = q;
  }
  
  for ( j2=n2-1; j2<n2+1; j2++ ) 
  for ( j3=n3-1; j3<n3+1; j3++ )  
  for ( j1=0;    j1<n1-1; j1++ ) {
  	  
  	  u1 = 1./h1;
   /* u2 = 0.; */
   /* u3 = 0.; */
  	  
  	  v1 = h1*b1[j1];
  	  
  	  for ( k2=0;    k2<n2-1; k2++)
  	  for ( k3=0;    k3<n3-1; k3++) { 
  	  	
  	  	  v2 = h2*b2[(j2-n2+2)*(n2+1)+k2];
  	      v3 = h3*b3[(j3-n3+2)*(n3+1)+k3];
  	      
          (*vz)[q]   = u1*v2*v3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jz)[q++] = (j1*(n2-1)+k2)*(n3-1)+k3;
      }
      
      (*iz)[++p] = q;
  }
  
  for ( j1=n1-1; j1<n1+1; j1++ ) 
  for ( j3=n3-1; j3<n3+1; j3++ )  
  for ( j2=0;    j2<n2-1; j2++ ) {
  	  
   /* u1 = 0.;*/
  	  u2 = 1./h2;
   /* u3 = 0.; */
  	  
  	  v2 = h2*b2[j2];
  	  
  	  for ( k1=0;    k1<n1-1; k1++)
  	  for ( k3=0;    k3<n3-1; k3++) { 
  	  	
  	  	  v1 = h1*b1[(j1-n1+2)*(n1+1)+k1];
  	      v3 = h3*b3[(j3-n3+2)*(n3+1)+k3];
  	      
          (*vz)[q]   = v1*u2*v3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jz)[q++] = (k1*(n2-1)+j2)*(n3-1)+k3;
      }
      
      (*iz)[++p] = q;
  }

  for ( j1=n1-1; j1<n1+1; j1++ ) 
  for ( j2=n2-1; j2<n2+1; j2++ )  
  for ( j3=0;    j3<n3-1; j3++ ) {
  	  
   /* u1 = 0.;
  	  u2 = 0.;
   */
  	  u3 = 1./h3;
  	  
  	  v3 = h3*b3[j3];
  	  
  	  for ( k1=0;    k1<n1-1; k1++)
  	  for ( k2=0;    k2<n2-1; k2++) { 
  	  	
  	  	  v1 = h1*b1[(j1-n1+2)*(n1+1)+k1];
  	      v2 = h2*b2[(j2-n2+2)*(n2+1)+k2];
  	      
          (*vz)[q]   = v1*v2*u3 + c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jz)[q++] = (k1*(n2-1)+k2)*(n3-1)+j3;
      }
      
      (*iz)[++p] = q;
  }

  for ( j1=n1-1; j1<n1+1; j1++ ) 
  for ( j2=n2-1; j2<n2+1; j2++ )  
  for ( j3=n3-1; j3<n3+1; j3++ ) {
  	  
   /* u1 = 0.;
  	  u2 = 0.;
  	  u3 = 0.;
   */
  	  for ( k1=0;    k1<n1-1; k1++)
  	  for ( k2=0;    k2<n2-1; k2++) 
  	  for ( k3=0;    k3<n3-1; k3++)  { 
  	  	
  	  	  v1 = h1*b1[(j1-n1+2)*(n1+1)+k1];
  	      v2 = h2*b2[(j2-n2+2)*(n2+1)+k2];
  	      v3 = h3*b3[(j3-n3+2)*(n3+1)+k3];
  	      
          (*vz)[q]   = c*v1*v2*v3;  /* dx1, dx2, dx3 */
          (*jz)[q++] = (k1*(n2-1)+k2)*(n3-1)+k3;
      }
      
      (*iz)[++p] = q;
  }

}


void SchurComp()
{
	
  for ( i=0; i<iz[0]; i++ )
  for ( j=0; j<iz[0]; j++ ) {
  	  k = 0;
  	  s = 0.;
      p = iz[i+2];
      q = iz[j+2];
      while ( jz[p]<iz[i+3] && jz[q]<iz[j+3] ) {
      	  
      	  if ( jz[p]<jz[q] ) { p++; continue; }
      	  if ( jz[q]<jz[p] ) { q++; continue; }
      	  	
      	  s += vz[p++]*vz[q++]/a[jz[q]];
      	  k++;
      }
      
      if ( k ) {
      	   
      	
      }
  	
  }
  
  
	
	
	
	
}