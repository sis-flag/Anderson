e#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define F77(a) d ## a ## _
typedef int INT;
typedef double R;

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
   	   (*b)[k]     = t;          /* (L0-L1)/2, psie */
   	   (*b)[k+n+1] = t;          /* (L0+L1)/2, psie */
   }
   
   s = 1./sqrt(90.);
   for ( k=i; k<n-1; k++ ) {     
   	   t           =  s*(*U)[i+(k-i)*(n-1)];
   	   (*b)[k]     =  t;        /* (L0-L1)/2, psio */
   	   (*b)[k+n+1] = -t;        /* (L0+L1)/2, psio */
   }
      
   (*b)[n-1]   =  2./3.;
   (*b)[n]     =  1./3.;
   (*b)[2*n]   =  1./3.;
   (*b)[2*n+1] =  2./3.;

   *b -= (n+1);
   
   return;                                                                     
}


void LegendspMat(INT n, INT **r, INT **c, R **a, R **b, R **U)
{  INT i,j,k;
   R   s,t;
	
   *U = (R *) malloc((n-1)*(n/2)*sizeof(R));
   *r = (INT *) malloc((n+2)*sizeof(INT));
   *c = (INT *) malloc(5*n*sizeof(INT));
   *b = (R   *) malloc(5*n*sizeof(R));
   *a = (R   *) malloc(5*n*sizeof(R));
   
   i   = n/2;
   *b += 3*(n-1);
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
   if ( n>=1 ) F77(steqr)("I", &i, *b,   *b+n,   *U,   &n, *b-3*n, &info);
   if ( n>=2 ) F77(steqr)("I", &j, *b+i, *b+i+n, *U+i, &n, *b-3*n, &info);
   *b -= 3*n;
   
   s = -1./sqrt(6.);
   for ( k=0; k<i; k++ ) {     
   	   (*r)[k]   =  q;
   	   t         =  s*(*U)[k*n];   	
   	   
   	   (*c)[q]   =  k;
   	   (*a)[q]   =  1.;
   	   (*b)[q++] =  (*b)[k+3*n];
   	   
   	   (*c)[q]   =  n;
   	   (*a)[q]   =  0.;
   	   (*b)[q++] =  t;

   	   (*c)[q]   =  n+1;
   	   (*a)[q]   =  0.;
   	   (*b)[q++] =  t;
   }
   
   s = 1./sqrt(90.);
   for ( k=i; k<n; k++ ) {     
   	   (*r)[k]   =  q;
   	   t         =  s*(*U)[i+(k-i)*n];   	
   	   
   	   (*c)[q]   =  k;
   	   (*a)[q]   =  1.;
   	   (*b)[q++] =  (*b)[k+3*n];
   	   
   	   (*c)[q]   =  n;
   	   (*a)[q]   =  0.;
   	   (*b)[q++] =  t;

   	   (*c)[q]   =  n+1;
   	   (*a)[q]   =  0.;
   	   (*b)[q++] = -t;
   }
   
   (*r)[n] = q;
   for ( k=0; k<n; k++ ) {
   	
   	   (*c)[q]   =  k;
   	   (*a)[q]   =  0.;
   	   (*b)[q++] = (*b)[3*k+1];
   }
   
   (*c)[q]   =  n;
   (*a)[q]   =  .5;
   (*b)[q++] =  2./3.;
   
   (*c)[q]   =  n+1;
   (*a)[q]   = -.5;
   (*b)[q++] =  1./3.;
   
   
   (*r)[n+1] = q;
   for ( k=0; k<n; k++ ) {
   	
   	   (*c)[q]   =  k;
   	   (*a)[q]   =  0.;
   	   (*b)[q++] = (*b)[3*k+2];
   }
   
   (*c)[q]   =  n;
   (*a)[q]   = -.5;
   (*b)[q++] =  1./3.;
   
   (*c)[q]   =  n+1;
   (*a)[q]   =  .5;
   (*b)[q++] =  2./3.;
   
   (*r)[n+2] = q;

   return;
	
}

void gdof(int m, int n, int *stif, int *mass)  
{
  *stif = m*(n+2)-(m>1?5:3);
  *mass = 5*m*n-2*m-4*n+1-(m>1?2:0);
}

void 
 off[0] = 0;
 off[1] = off[0] + (n1-1)*(n2-1);
 off[2] = off[1] + (n1-1)*2;
 off[3] = off[2] + (n2-1)*2;
 
 for ( j1=0; j1<n1-1; j1++ )
 for ( j2=0; j2<n2-1; j2++ ) {
 	
 	 /* (k1, k2)  = (j1,j2) */
 	   col[q]   = off[0] + j1*(n2-1)+j2;
     st1[q]   = b2[j2];     /* dx */
     st2[q]   = b1[j1];     /* dy */
     mas[q++] = b1[j1]*b2[j2];
     
     /* (k1, k2) = (j1, {n2-1,n2}) */    
     for ( k2=n2-1; k2<n2+1; k2++ ) {
    	   col[q]   = off[1] + (k2-n2+1)*(n1-1)+j1;
         st1[q]   = b2[j2+(k2-n2+2)*(n2+1)];
         st2[q]   = 0.;
         mas[q++] = b1[j1]*b2[j2+(k2-n2+2)*(n2+1)];
     }
     
     /* (k1,k2) = {n1-1,n1}, j2) */
     for ( k1=n1-1; k1<n1+1; k1++ ) {
    	   col[q]   = off[2] + (k1-n1+1)*(n2-1)+j2;
         st1[q]   = 0.;        
         st2[q]   = b1[j1+(k1-n1+2)*(n1+1)];                 
         mas[q++] = b1[j1+(k1-n1+2)*(n1+1)]*b2[j2];
     } 

     /* (k1,k2) = ({n1-1,n1}, {n2-1,n2}) */ 
     for ( k1=n1-1; k1<n1+1; k1++ )   
     for ( k2=n2-1; k2<n2+1; k2++ ) {
    	   col[q]   = off[3] + (k1-n1+1)*2+(k2-n2+1);
         st1[q]   = 0.;        
         st2[q]   = 0.;                 
         mas[q++] = b1[j1+(k1-n1+2)*(n1+1)]*b2[j2+(k2-n2+2)*(n2+1)];
     }
     row[++p] = q;  
 } 
 
 for ( j2=n2-1; j2<n2+1; j2++ )  
 for ( j1=0;    j1<n1-1; j1++ ) { /* (j1,{n2-1,n2}) */
 	
 	 /* (k1, k2) = (j1, k2) */
 	 for ( k2=0; k2<n2-1; k2++ ) {
    	   col[q]   = off[0] + j1*(n2-1)+k2;
         st1[q]   = b2[k2+(j2-n2+2)*(n2+1)];
         st2[q]   = 0.;
         mas[q++] = b1[j1]*b2[k2+(j2-n2+2)*(n2+1)];
     }
     
     /* (k1, {n2-1,n2}) = (j1, {n2-1,n2}) */
     for ( k2=n2-1; k2<n2+1; k2++ ) {
    	   col[q]   = off[1] + (k2-n2+1)*(n1-1)+j1;
         st1[q]   = b2[k2+(j2-n2+2)*(n2+1)];
         st2[q]   = 0.;
         mas[q++] = b1[j1]*b2[k2+(j2-n2+2)*(n2+1)];
     }
     
     /* ({n1-1,n1},k2) */
     for ( k1=n1-1; k1<n1+1; k1++ )
     for ( k2=0;    k2<n2-1; k2++ ) {
    	   col[q]   = off[2] + (k1-n1+1)*(n2-1)+k2;
         st1[q]   = 0.;
         st2[q]   = 0.;
         mas[q++] = b1[j1+(k1-n1+2)*(n1+1)]*b2[k2+(j2-n2+2)*(n2+1)];
     }

     /* (k1,k2) = ({n1-1,n1}, {n2-1,n2}) */ 
     for ( k1=n1-1; k1<n1+1; k1++ )   
     for ( k2=n2-1; k2<n2+1; k2++ ) {
     	   col[q]   = off[3] + (k1-n1+1)*2+(k2-n2+1);
         st1[q]   = 0.;        
         st2[q]   = b1[j1+(k1-n1+2)*(n1+1)]*(k2==j2?0.5,-0.5);                 
         mas[q++] = b1[j1+(k1-n1+2)*(n1+1)]*b2[k2+(j2-n2+2)*(n2+1)];
     }
     row[++p] = q;  
 }
 
 for ( j1=n1-1; j1<n1+1; j1++ )  
 for ( j2=0;    j2<n2-1; j2++ ) { /* ({n1-1,n1},j2) */
 	
 	 /* (k1, k2) = (k1, j2) */
 	 for ( k1=0; k1<n1-1; k1++ ) {
    	   col[q]   = off[0] + k1*(n2-1)+j2;
         st1[q]   = 0.;
         st2[q]   = b1[k1+(j1-n1+2)*(n1+1)];
         mas[q++] = b1[k1+(j1-n1+2)*(n1+1)]*b2[j2];
     }
     
     /* (k1, {n2-1,n2}) */    
     for ( k2=n2-1; k2<n2+1; k2++ )
     for ( k1=0;    k1<n1-1; k1++ ) {
    	   col[q]   = off[1] + (k2-n2+1)*(n1-1)+k1;
         st1[q]   = 0.;
         st2[q]   = 0.;
         mas[q++] = b1[k1+(j1-n1+2)*(n1+1)]*b2[j2+(k2-n2+2)*(n2+1)];
     }
     
     /* ({n1-1,n1},k2) = ({n1-1,n1},j2) */
     for ( k1=n1-1; k1<n1+1; k1++ ) {
    	   col[q]   = off[2] + (k1-n1+1)*(n2-1)+j2;
         st1[q]   = (j1==k1?.5,-.5)*b2[j2];
         st2[q]   = b1[k1+(j1-n1+2)*(n1+1)];
         mas[q++] = b1[k1+(j1-n1+2)*(n1+1)]*b2[j2];
     }

     /* (k1,k2) = ({n1-1,n1}, {n2-1,n2}) */ 
     for ( k1=n1-1; k1<n1+1; k1++ ) 
     for ( k2=n2-1; k2<n2+1; k1++ ) {
    	   col[q]   = off[3] + (k1-n1+1)*2+(k1-n2+1);
         st1[q]   = (j1==k1?.5,-.5)*b2[j2+(k2-n2+2)*(n2+1)];        
         st2[q]   = 0.;                 
         mas[q++] = b1[k1+(j1-n1+2)*(n1+1)]*b2[j2+(k2-n2+2)*(n2+1)];
     }
     row[++p] = q;  
 } 
 
 
 
void assemble2d()

{
	off[0] = 0;
	off[1] = off[0] + m1*m2(n1-1)*(n2-1);
	off[2] = off[1] + (m2-1)*m1*(n1-1);
	off[3] = off[2] + (m1-1)*m2*(n2-1);
	
	ofe[0] = 0;
	ofe[1] = ofe[0] + (n1-1)*(n2-1);
	ofe[2] = ofe[1] + (n1-1);
	ofe[3] = ofe[2] + (n1-1);
	ofe[4] = ofe[3] + (n2-1);
	ofe[5] = ofe[4] + (n2-1);
	
	/*   ofe[0] + j1*(n2-1)+j2
	 *   ofe[1] + j1
	 *   ofe[2] + j1
	 *   ofe[3] + j2
	 *   ofe[4] + j2
	     
	 */	
	for ( i1=0; i1<m1; i1++ )
	for ( i2=0; i2<m2; i2++ ) {
		  i  = i1*m2+i2;
		  x1 = h2[i2]/h1[i1];
		  x2 = h1[i1]/h2[i2];
		  x3 = h1[i1]*h2[i2]*c[i];
		  
	    
	    for ( j=off[0]; j<off[1]; j++ ) {
	    	
	    	  /* F(i1,i2,j1,j2), */
	    	  for ( l=row[j]; l<row[j+1]; l++ {
	            
	            k     = col[l];
	            va[q] = st1[l]*x1 + st2[l]*x2 + mas[l]*x3;
	            
	    	      if ( k<ofe[1] ) {                                        /* F(i1,i2,j1,j2) */ 
	    	      	
                  ja[q++]	=  off[0] + i*(n1-1)*(n2-1) + k;             /* (i*(n1-1)+j1)*(n2-1)+j2 */  	    
                  
	            } else if ( k<ofe[2] ) {                                 /* E2L(i1,i2,j1,n2-1) */
	            	  if ( i2>0 )
                  ja[q++] = off[1] + ((i2-1)*m1+i1)*(n1-1) + k-ofe[1]; /* ((i2-1)*m1+i1)*(n1-1)+j1 */
                  
	            } else if ( k<ofe[3] ) {                                 /* E2R(i1,i2,j1,n2) */
	            	  if ( i2<m2-1 )
                  ja[q++] = off[1] + (i2*m1+i1)*(n1-1) + k-ofe[2];     /* (i2*m1+i1)*(n1-1)+j1 */
                  
              } else if ( k<ofe[4] ) {                                 /* E1L(i1,i2,n1-1,j2) */
              	  if (i1>0 )
              	  ja[q++] = off[2] + ((i1-1)*m2+i2)*(n2-1) + k-ofe[3]; /* ((i1-1)*m2+i2)*(n2-1)+j2 */
              	  
              } else if ( i1<m1-1 && k<ofe[5] ) {                      /* E1L(i1,i2,n1-1,j2) */
              	  if ( i1<m1-1 )
              	  ja[q++] = off[2] + (i1*m2+i2)*(n2-1) + k-ofe[3];     /* ((i1-1)*m2+i2)*(n2-1)+j2 */

              } else if ( k>=ofe[5] ) {                                /* V(i1,i2,{n1-1,n1},{n2-1,n2}) */
              	
              	                     (i1-1-(k-ofe[5])/2)*(m2-1)+(i2-1+(k-ofe[5])%2);
              	                     
              	  ja[q++] = off[1] + (i1-1)*m2+i2 + k-ofe[5]; /* V(k) */
              }
	            	
	        
	        if ( i2<m2-1 ) for ( ; k<row[j]+2; k++ )  {     /* F(i1,i2,j1,j2), E2L(i1,i2,j1,n2-1) */
              ja[q]	  = ((i2-1)*m1+i1)*(n1-1) + col[k]-ofe[1];         /*  ((i2-1)*m1+i1)*(n1-1)+j1 */; 
	   	        va[q++] = st1[k]*x1 + st2[k]*x2 + mas[k]*x3;
	        }
	        
	        
	        
	        
	        ia[++p] = q;
	    }
  }
	    
	    
 	
 	
 	
}
 
 
 
 
 for ( j1=n1-1; j1<n1+1; j1++ )  
 for ( j2=n2-1; j2<n2+1; j2++ ) { /* ({n1-1,n1},{n2-1,n2}) */
 	
 	 /* (k1, k2) */
 	 for ( k1=0; k1<n1-1; k1++ )
 	 for ( k2=0; k2<n2-1; k2++ ) {
    	   col[q]   = off[0] + k1*(n1-1)+k2;
         st1[q]   = 0.;
         st2[q]   = 0.;
         mas[q++] = b1[k1+(j1-n1+2)*(n1+1)]*b2[k2+(j2-n2+2)*(n2+1)];
     }
     
     /* (k1, {n2-1,n2}) */    
     for ( k2=n2-1; k2<n2+1; k2++ )
     for ( k1=0;    k1<n1-1; k1++ ) {
    	   col[q]   = off[1] + (k2-n2+1)*(n1-1)+k1;
         st1[q]   = 0.;
         st2[q]   = b1[k1+(j1-n1+2)*(n1+1)];
         mas[q++] = b1[k1+(j1-n1+2)*(n1+1)]*b2[j2+(k2-n2+2)*(n2+1)];
     }
     
     /* ({n1-1,n1},k2) */
     for ( k1=n1-1; k1<n1+1; k1++ )
     for ( k2=0;    k2<n2-1; k2++ ) {
    	   col[q]   = off[2] + (k1-n1+1)*(n2-1)+k2;
         st1[q]   = (j1==k1?.5,-.5)*b2[k2+(j2-n2+2)*(n2+1)];
         st2[q]   = 0.;
         mas[q++] = b1[k1+(j1-n1+2)*(n1+1)]*b2[k2+(j2-n2+2)*(n2+1)];
     }

     /* ({n1-1,n1}, {n2-1,n2}) */ 
     for ( k1=n1-1; k1<n1+1; k1++ ) 
     for ( k2=n2-1; k2<n2+1; k2++ ) {
    	   col[q]   = off[2] + (k1-n1+1)*2+(k2+n2+1);
         st1[q]   = (j1==k1?.5,-.5)*b2[k2+(j2-n2+2)*(n2+1)];        
         st2[q]   = b1[k1+(j1-n1+2)*(n1+1)]*(j2==k2?.5,-.5);                 
         mas[q++] = b1[k1+(j1-n1+2)*(n1+1)]*b2[k2+(j2-n2+2)*(n2+1)];
     }
     row[++p] = q;  
 } 
 
 return;
} 





      
{


  g1[0] = 0;
  g1[1] = n1-1;
  g1[2] = n1+1;
  
  g2[0] = 0;
  g2[1] = n2-1;
  g2[2] = n2+1;

  for ( e1=0; e1<2; e1++ ) 
  for ( e2=0; e2<2; e2++ ) 
  for ( j1=g1[e1]; j1<g1[e1+1]; j1++ )
  for ( j2=g2[e2]; j2<g2[e2+1]; j2++ )
  	
      for ( f1=0; f1<2; f1++ ) {
      if ( e1+f1 ) {
          h1[0] = g1[f1];
          h1[1] = g1[f1+1];
      } else {
      	  h1[0] = j1;
      	  h1[1] = j1+1;
      }
      	  	
      for ( f2=0; f2<2; f2++ ) {
      if ( e2+f2 ) {
          h2[0] = g2[f2];
          h2[1] = g2[f2+1];
      } else {
      	  h2[0] = j2;
      	  h2[1] = j2+1;
      }
      	
            		
      for ( k1=h1[0]; k1<h1[1]; k1++ )  
      for ( k2=h2[0]; k2<h2[1]; k2++ ) {
      	
          if ( e1 ) {
          	  i1 = (j1-n1+2)*(n1+1) + k1;
          	   
          	  if ( f1 ) a1 = (j1==k1?.5,-.5);
          	  else      a1 = 0.;	
          else {
          	  if ( f1 ) {
          	   	  i1 = (k1-n1+2)*(n1+1) + j1;
          	   	  a1 = 0.;
          	  } else {
          	   	  i1 = k1;
          	   	  a1 = 1.;
          	  }
          }	
          	
          if ( e2 ) {
          	  i2 = (j2-n2+2)*(n2+1) + k2;
          	   
          	  if ( f2 ) a2 = (j2==k2?.5,-.5);
          	  else      a2 = 0;	   
          else {
          	  if ( f2 ) {
          	   	  i2 = (k2-n2+2)*(n2+1) + j2;
          	   	  a2 = 0.;
          	  } else {
          	      i2 = k2;
          	 	   a2 = 1.;
          	  }
          }
          
          st1[q]   = a1*b2[i2];
      		st2[q]   = b1[i1]*a2;
      		mas[q++] = b1[i1]*b2[i2];
  	  }
  	}}
  }
   	
	return;
	
}


                                                              
/*
 stiff:
   j<n-1,  {(i,j)} 
            1
   j=n,    {(i,n-1), i>0}, {(i,n)}+{(i+1,n-1)}, {(i+1,n), i<m-2}
           -0.5             0.5   +   0.5       -0.5
 mass:
   j<n-1,  {(i,j)}, {(i,n-1), i>0}, {(i,n), i<m-1} 
            b[j]     b[j+n+1]       b[j+2*n+2]
    
            (i,n)              (i+1,n-1)
   j=n,    {(i,k), 0<=k<n-1}, {(i+1,k), 0<=k<n-1,i<m-2}
            b[k+2*n+2]         b[k+n+1]
            
           {(i,n-1), i>0}, {(i,n)}+{(i+1,n)}, {(i+1,n-1), i<m-2}
            b[2*n+1]       b[3*n+2] + b[2*n+1]   b[3*n+1]
            
*/  

void neighbors(INT m, INT i, INT n, INT j, INT *g, R *s, R *t)  
{ 
  /* for 0<=j<n and j=n+1 */
  g[0] = (j>n);

  /*  j=n+1 only */  
  g[1] =  1;
  g[2] = -1;
  if ( i>0 ) {      /* (i,   j,   n) */
  	  g[1]  = 0;
  	  g[2]  = 0;
  	  s[0]  = (g[0]?-.5/h[i]:.0);
  	  t[0]  = b[j+n+2]*h[i];
  }     
  if ( i<m-1 ) {    /* (i,   j,   n+1) */
  	  g[2]  = 1;
  	  s[1]  = (g[0]?.5/h[i]:.0);
  	  t[1]  = b[j+2*n+4]*h[i];
  }     
  if ( g[2]==1 && j==n+1 ) {
  	  s[1] += .5/h[i+1];
      t[1] += b[2*n+2]*h[i+1];  /* (i+1, n, n) */
  	
      if ( i<m-2  ) {           /* (i+1, n, n+1) */
      	  g[2]  = 2;
     	  s[1]  = -.5/h[i+1];
          t[2]  = b[3*n+4]*h[i+1];
      }
  } 
}	

void entries()
{

	/* F, (i1+e1,i2+e2) (j1-e1,j2-e2) 
 	                    (k1,   k2   )  */
	for ( e1=0; e1<=g1[0]; e1++ ) 
	for ( e2=0; e2<=g2[0]; e2++ ) {
	    i  = (i1+e1)*m2+i2+e2;
	    
	    for ( k1=(g1[0]?0:j1); k1<=(g1[0]?n1-1:j1); k1++ ) 
	    for ( k2=(g2[0]?0:j2); k2<=(g2[0]?n2-1:j2); k2++ ) {
	    	
          jv[*q] = off[0] + (i*n1+k1)*n2+k2;
          vv[*q] = a1[g[0]*(2-e1)*(n+2)+k1]*b2[k2] + b1[k1]*a2[k2] + b1[k1]*b2[k2];
          (*q)++;
	    }
	}
	/* k   k+2*(n+2) 
	 * k   k+(n+2)
	 */
	
	for ( e2=g2; e2<=h2; e2++ )
	for ( e1=0;  e1<=f1; e1++ ) {
        i = (i2-1+e2)*m1+i1+e1;
        		
	    for ( k1=r[j1-e1]; k1=ia[j1+1-e1]; k1++ ) {
	    	
	    	jv[*q] = off[1] + i*n1+ja[k1];
	        vv[*q] = a[k1]*t2[e2];
            (*q)++;		    	
	    }
	}	
	
	
	
}

void assemble2d(int m1, int m2, int n1, int n2, R *h1, R *h2, 
               R *a1, R *b1, R *a2, R *b12, 
                int **iu, int **ju, R **vu, 
                int **iv, int **jv, R **vv)
{ int i1,i2,i3, j1,j2,j3, k1,k2,k3, p, q, info;
  
  p   = m*n-1;
  q   = m*(n+2)-(m>1?5:3);
  *iu = (INT *) malloc((p+2)*sizeof(INT));
  *ju = (INT *) malloc(q*sizeof(INT));
  *vu = (R   *) malloc(q*sizeof(R));
  
  (*iu)[0] = p;
  (*iu)++;
  (*iu)[0] = 0;
  p = q = 0;
  
  for ( i1=0; i1<m1;  i1++ ) 
  for ( i2=0; i2<m2;  i2++ ) 
  for ( j1=2; j1<=n1; j1++ )
  for ( j2=2; j2<=n2; j2++ ) {
  	  u1[0] = 0.;             v1[0] = h1[i1]*a1[j1];   
  	  u1[1] = 0.;             v1[1] = h1[i1]*a1[j1+1+n1];
  	  u1[2] = 1./h1[i1];      v1[2] = h1[i1]*a1[j1+2+n1*2];

  	  u2[0] = 0.;             v2[0] = h2[i2]*a2[j2];   
  	  u2[1] = 0.;             v2[1] = h2[i2]*a2[j2+1+n2];
  	  u2[2] = 1./h2[i2];      v2[2] = h2[i2]*a2[j2+2+n2*2];  	    
  	
  	  /* F(i1,i2,j1,j2) */
      (*ju)[q]   = p;
      (*vu)[q++] = u1[2]*v2[2]           /* dx1 */
                 + v1[2]*u2[2]           /* dx2 */
                 + v1[2]*v2[2]*c[i1*m2+i2];

      for ( k2=(i2?0:1); k2<=(m2-1-i2?1:0); k2++ ) { /* F(i1,i2,j1,k2={0,1}) */
      	  	
      	  (*ju)[q]   = off[1] + ( ((i2+k2-1)*m1+i1 )*(n1-1)+j1;
      	  (*vu)[q++] = u1[2]*v2[k2]
      	             /* + v1[2]*u2[k2] */ /* =0 */
      	             + v1[2]*v2[k2]*c[i1*m2+i2];
      }      	             

      for ( k1=(i1?0:1); k1<=(m1-1-i1?1:1); k1++ ) { /* F(i1,i2,k1={0,1},j2) */
      	  	
      	  (*ju)[q]   = off[2] + (((i1+k1-1)*m2+i2)*(n2-1)+j2;
      	  (*vu)[q++] = /* u1[k1]*v2[2]  */ /* =0 */
      	             + v1[k1]*u2[2]
      	             + v1[k1]*v2[2]*c[i1*m2+i2];
      }
      
      for ( k1=(i1?0:1); k1<=(m1-1-i1?1:1); k1++ ) 
      for ( k2=(i2?0:1); k2<=(m2-1-i2?1:0); k2++ ) { 
      	 	
      	  (*ju)[q]   = off[3] + (i1+k1-1)*(m2-1)+(i2+k2-1);
          (*vu)[q++] = /*u1[k1]*v2[k2] */  /*  0 */ 
      	             /*+ v1[k1]*u2[k2] */  /*  0 */
      	             + v1[k1]*v2[k2]*c[i1*m2+i2];
      }
      
      (*iu)[++p] = q;
  }
  
  for ( i2=0; i2<m2-1; i2++ ) /* (i1,i2,j1,1) = (i2,i2+1,j1,0) */
  for ( i1=0; i1<m1;   i1++ ) 
  for ( j1=2; j1<=n1;  j1++ ) { 
  	  u1[0] = 0.;             v1[0] = h1[i1]*a1[j1];   
  	  u1[1] = 0.;             v1[1] = h1[i1]*a1[j1+1+n1];
  	  u1[2] = 1./h1[i1];      v1[2] = h1[i1]*a1[j1+2+n1*2];
  	
  	  for ( e2=0; e2<2;   e2++ )     /*R: (i1,i2+e2,j1,1-e2) */
  	  for ( k2=2; k2<=n2; k2++ ) {   /*   (i1,i2+e2,j1,k2)   */
  	      u2[0] = 0.;         v2[0] = h2[i2+e2]*a2[k2+(1-e2)*(n2+1)];   
  	  	
  	      (*ju)[q]   = off[0] + ((i1*m2+i2+e2)*(n1-1)+j1)*(n2-1)+k2;
  	      (*vu)[q++] = u1[2]*v2[0]
  	      	         /*+ v1[2]*u2[0] */ /* 0 */
                     + v1[2]*v2[0]*c[i1*m2+i2+e2];
      }               

  	  u2[0] = -.5/h2[i2];                 
  	  u2[1] = .5*/h2[i2]+.5/h2[i2+1];  
  	  u2[2] = -.5/h2[i2+1];            
  	  v2[0] = h2[i2]*a2[1];
  	  v2[1] = h2[i2]*a2[n2+2]+h[i2+1]*a2[0];
  	  v2[2] = h2[i2+1]*a2[1];
  	  w2[0] = v2[0]*c[i1*m2+i2];
  	  w2[1] = h2[i2]*a2[n2+2]*c[i1*m2+i2]+h[i2+1]*a2[0]*c[i1*m2+i2+1];
  	  w2[2] = v2[2]*c[i1*m2+i2+1];
  	  
  	  for ( e2=(i2?0:1); e2<(m2-2-i2?2:1); e2++ ) { /*R: (i1,i2+e2,j1,1-e2)    */
  	                                                /*   (i1,i2+e2,j1,k2={0,1} */
  	  	  (*ju)[q]   = off[1] + ((i2-1+e2)*m1+i1)*(n1-1)+j1;
  	  	  (*vu)[q++] = u1[0]*v2[e2]
  	  	             + v1[0]*u2[e2]
                     + v1[0]*w2[e2];
      }
      
      
      for ( e1=(i1?0:1); e1<(m1-1-i1?1:0); e1++ )
  	  for ( e2=(i2?0:1); e2<(m2-2-i2?1:0); e2++ )  /*R: (i1,i2+e2,j1,1-e2)  */
  	  for ( k2=2; k2<=n2; k2++ )                    /*  (i1,i2+e2,{0,1},k2) */
  	  	  (*ju)[q]   = off[2] + ((i1-1+e1)*m2+i2+e2)*(n2-1)+k2;
  	  	  (*vu)[q++] = /* u1[0]*v2[e2] */ /* 0 */
  	  	             /* + v1[0]*u2[e2] */ /* 0 */
                     + h1[i1]*v1[e1]*h2*[i2+e2]*a2[k2+(1-e2)*(n2+1)]*c[i1*m2+i2+e2];
      }
      
      /*R: (i1,i2+e2,j1,1-e2)     */
      /*   (i1,i2+e2,{0,1},{0,1}) */
      
      for ( e1=(i1?0:1); e1<(m1-1-i1?1:0); e1++ )
      for ( e2=(i2?0:1); e2<(m2-2-i2?2:1); e2++ ) { 
  	                                                
  	  	  (*ju)[q]   = off[3] + ((i1-1+e1)*m2+i2-1+e2);
  	  	  (*vu)[q++] = u1[e1]*v2[e2]
  	  	             + v1[e1]*u2[e2]
                     + v1[e1]*w2[e2];
      }
  	   (*iu)[++p] = q;
  }
  	
  	
  for ( i1=0; i1<m1-1; i1++ ) /* (i1,i2,1,j2) = (i1+1,i2,0,j2) */
  for ( i2=0; i2<m2;   i2++ ) 
  for ( j2=2; j2<=n2;  j2++ ) { 
  	  u2[0] = 0.;             v2[0] = h2[i2]*a2[j2];   
  	  u2[1] = 0.;             v2[1] = h2[i2]*a2[j2+1+n2];
  	  u2[2] = 1./h2[i2];      v2[2] = h2[i2]*a2[j2+2+n2*2];
  	
  	  for ( e1=0; e1<2;   e1++ )     /*R: (i1+e1,i2,1-e1,j2) */
  	  for ( k1=2; k1<=n1; k1++ ) {   /*   (i1+e1,i2,k1,j2)   */
  	      u2[0] = 0.;         v2[0] = h2[i2+e2]*a2[k2+(1-e2)*(n2+1)];   
  	  	
  	      (*ju)[q]   = off[0] + ((i1*m2+i2+e2)*(n1-1)+j1)*(n2-1)+k2;
  	      (*vu)[q++] = u1[2]*v2[0]
  	      	         /*+ v1[2]*u2[0] */ /* 0 */
                     + v1[2]*v2[0]*c[i1*m2+i2+e2];
      }               



  
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
	
	R *a, *b, *V;
	
	LegendMat(n, &a, &b, &V);
	
	int *iu, *ju, *iv, *jv;
	R   *vu, *vv;
	
	assemble1d(m, n, h, a, b, &iu, &ju, &vu, &iv, &jv, &vv);
	
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