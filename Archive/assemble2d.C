#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef int INT;
typedef double R;

/*  0<=i1,i2<=m-1, 0<=j1,j2<=n-1  without specifying

    i1*m+i2,                       j1*n+j2
    (i1,  i2 ) (j1, j2 )

    i2*m+i1,                       j1           i2<m-1
    (i1,  i2   ) (j1, n   )
    (i1,  i2+1 ) (j1, n+1 )

    i1*m+i2,                       j2           i1<m-1
    (i1,   i2 ) (n,   j2 )
    (i1+1, i2 ) (n+1, j2 )

    i1*(m-1)+i2                    0            i1,i2<m-1
    (i1,  i2  ) (n,  n  )
    (i1,  i2+1) (n,  n+1)
    (i1+1,i2  ) (n+1,n  )
    (i1+1,i2+1) (n+1,n+1)
*/


INT l2gindex(INT m1, INT m2, INT i1, INT i2, 
             INT n1, INT n2, INT j1, INT j2)
{	     
 INT  off[4], p=-1;

 n1--;
 n2--;

 off[0] = 0;
 off[1] = off[0] + m1*m2*n1*n2;
 off[2] = off[1] + (m2-1)*m1*n1;
 off[3] = off[2] + (m1-1)*m2*m2;

/* local-to-gobal mapping*/
  if ( j1<n1 ) {
  	 if ( j2<n2 ) { 
  	 	 p = off[0] + ((i1*m2+i2)*n1+j1)*n2+j2;    
  	 } else if ( j2==n2   && i2<m2-1 ) { /* F3R */
  	 	 p = off[1] + (i2*m1+i1)*n1+j1;
  	 } else if ( j2==n2+1 && i2>0    ) {  /* F3L */
  	 	 p = off[1] + ((i2-1)*m1+i1)*n1+j1; 	 	 	 
  	 }
  } else if ( j1==n1 && i1<m1-1 ) { 
  	 if ( j2<n2 ) {                       /* F2R */
  	 	 p = off[2] + (i1*m2+i2)*n2+j2;
  	 } else if ( j2==n2   && i2<m2-1 ) {  /* E1RR */
  	 	 p = off[3] + i1*(m2-1)+i2;
  	 } else if ( j2==n2+1 && i2>0    ) {  /* E1RL */
  	 	 p = off[3] + i1*(m2-1)+i2-1;;
  	 }
  } else if ( j1==n1+1 && i1>0 ) {
  	 if ( j2<n2 ) {                       /* F2L */
  	 	 p = off[2] + ((i1-1)*m2+i2)*n2+j2;
  	 } else if ( j2==n2   && i2<m2-1 ) {  /* E1LR */
  	 	 p = off[3] + (i1-1)*(m2-1)+i2;  
  	 } else if ( j2==n2+1 && i2>0    ) {  /* E1LL */
  	 	 p = off[3] + (i1-1)*(m2-1)+i2-1;
  	 }
  }
 
  return p;
}

INT l2gdeg(INT m, INT i, INT n, INT j) 
{	     
 INT  x;

  if ( (j==n-1 && i==m-1) || (j==n && i==0) ) return -1;
  	
  if ( j==n )
  	 x = 7 - 2 * ((n<3)+(n<2)) - (i==1) - (i==m-1);
  else if ( j==n-1 )
  	 x = 7 - 2 * ((n<3)+(n<2)) - (i==0) - (i==m-2);
  else if ( j==n-2 )  
  	 x = 4 - (n<5) - (i==0) - (i==m-1);
  else if ( j==n/2-1 )  
  	 x = 4 - (n<4) - (i==0) - (i==m-1);
 else 
 	 x = 3 - (j==0||j==n/2); 
 
 return x;
}

INT gdof(INT m, INT n) 
{	     
  INT  i,j, x=0;
  
  for ( i=0; i<m; i++ ) 
  for ( j=0; j<(i==m-1?n-1:n); j++ )    
	  x += l2gdeg(m,i,n,j);

 return x;
}



void LegendspMat(INT n, R **va, INT **ia, INT **ja, R **vb, INT **ib, INT **jb)
{ INT i, j, k;

/* basis in order 
    (L[k]-L[k-2])/sqrt(4*k-2)  k=n,  n-2,n-4,..., k>=2,
    (L[k]-L[k-2])/sqrt(4*k-2)  k=n-1,n-3,n-5,..., k>=2,
    (1+x)/2                      
    (1-x)/2                      
*/

  *va = (R *)  malloc((n+3)*sizeof(R));
  *ia = (INT *)malloc((n+3)*sizeof(int));
  *ja = (INT *)malloc((n+3)*sizeof(int));
 
  (*ia)[0] = n+1;
  
  for ( j=k=0; k<n-1; k++ ) {
  	  (*ia)[k+1] = j;
  	  (*ja)[j]   = k;
  	  (*va)[j++]  = 1.;
  }
  
  (*ia)[n]   = j;
  (*ja)[j]   = n-1;
  (*va)[j++] =  .5;
  (*ja)[j]   = n;
  (*va)[j++] = -.5;
  
  (*ia)[n+1] = j;
  (*ja)[j]   = n-1;
  (*va)[j++] = -.5;
  (*ja)[j]   = n;
  (*va)[j++] =  .5;

  (*ia)[n+2] = j;  
  
/****************************************/
  
  *vb = (R *)  malloc((3*n+5)*sizeof(R)*2);
  *ib = (INT *)malloc((n+3)*sizeof(R)*2);
  *jb = (INT *)malloc((3*n+5)*sizeof(R)*2);
  
  (*ib)[0] = n+1;
  (*ib)[1] = j = 0;
  for ( k=n; k>1; k-=2 ) { /* (L[k]-L[k-2])/sqrt(4*k-2) */
  	  i = (n-k)/2;
  	  
  	  if ( k+2<=n ) {
  	  	(*jb)[j]   = i-1;
  	  	(*vb)[j++] = -.25 / ( sqrt( (k-.5) * (k+1.5) ) * (k+0.5) );
  	  }
  	  
  	  (*jb)[j]   = i;
  	  (*vb)[j++] = ( .25/(k+.5) + .25/(k-1.5) ) / (k-0.5);
  	  
  	  if ( k>3 ) {
  	  	(*jb)[j]   = i+1;
  	  	(*vb)[j++] = -.25 / ( sqrt( (k-2.5) * (k-.5) ) * (k-1.5) );
  	  } else if ( k==3 ) {
  	  	(*jb)[j]   = n-1;
  	  	(*vb)[j++] = -1./sqrt(90.);
  	  	(*jb)[j]   = n;
  	  	(*vb)[j++] =  1./sqrt(90.);
  	  } else {
  	  	(*jb)[j]   = n-1;
  	  	(*vb)[j++] = -1./sqrt(6.);
  	  	(*jb)[j]   = n;
  	  	(*vb)[j++] = -1./sqrt(6.);
  	  }
  	  
  	  (*ib)[i+2] = j;
  }
  
  for ( k=n-1; k>1; k-=2 ) { /* (L[k]-L[k-2])/sqrt(4*k-2) */
  	  i = n-k/2-1;
  	  
  	  if ( k+2<=n ) {
  	  	(*jb)[j]   = i-1;
  	  	(*vb)[j++] = -.25 / ( sqrt( (k-.5) * (k+1.5) ) * (k+0.5) );
  	  }
  	  
  	  (*jb)[j]   = i;
  	  (*vb)[j++] = ( .25/(k+.5) + .25/(k-1.5) ) / (k-0.5);
  	  
  	  if ( k>3 ) {
  	  	(*jb)[j]   = i+1;
  	  	(*vb)[j++] = -.25 / ( sqrt( (k-2.5) * (k-.5) ) * (k-1.5) );
  	  } else if ( k==3 ) {
  	  	(*jb)[j]   = n-1;
  	  	(*vb)[j++] = -1./sqrt(90.);
  	  	(*jb)[j]   = n;
  	  	(*vb)[j++] =  1./sqrt(90.);
  	  } else {
  	  	(*jb)[j]   = n-1;
  	  	(*vb)[j++] = -1./sqrt(6.);
  	  	(*jb)[j]   = n;
  	  	(*vb)[j++] = -1./sqrt(6.);
  	  }
  	  
  	  (*ib)[i+2] = j;
  }
  
  /* (L[0]+L[1])/2 */
  if ( n>=2 ) {
     (*jb)[j]   = n/2-1;
     (*vb)[j++] = (k?-1./sqrt(6):-1./sqrt(90.));
  } 
  if ( n>=3 ) {
     (*jb)[j]   = n-2;
     (*vb)[j++] = (k?-1./sqrt(90.):-1./sqrt(6.));
  }
  
  (*jb)[j]   = n-1;
  (*vb)[j++] = 2./3;
  (*jb)[j]   = n;
  (*vb)[j++] = 1./3;
  (*ib)[n+1] = j;
  
  /* (L[0]-L[1])/2 */
  if ( n>=2 ) {
     (*jb)[j]   = n/2-1;
     (*vb)[j++] = (k?-1./sqrt(6):1./sqrt(90.));
  } 
  if ( n>=3 ) {
     (*jb)[j]   = n-2;
     (*vb)[j++] = (k?1./sqrt(90.):-1./sqrt(6.));
  }
  
  (*jb)[j]   = n-1;
  (*vb)[j++] = 1./3;
  (*jb)[j]   = n;
  (*vb)[j++] = 2./3;
  (*ib)[n+2] = j;
 return;
  
}   



int RowEntries(int m1, int m2,  int i1,  int i2, 
               int n1, int n2,  int j1,  int j2,
               int *ia, int *ja, R *a,    int *ib, int *jb, R *b,
               int *q,  int *jv, R *vv)
{
 int k1,k2, f1,f2, g1,g2, h1,h2, e1,e2, i,
     off[4], t1[3],t2[3];
	
 off[0] = 0;
 off[1] = off[0] + m1*m2*n1*n2;
 off[2] = off[1] + (m2-1)*m1*n1;
 off[3] = off[2] + (m1-1)*m2*n2;
 	
   	if  ( j1>n1 || j1<0 ) return 1;
   	else if ( j1<n1 )     f1 = 0;
   	else if ( i1==m1-1 )  return 1;	
   	else                  f1 = 1;
   	
   	if  ( j2>n2 || j2<0 ) return 1;
   	else if ( j2<n2 )     f1 = 0;
   	else if ( i2==m2-1 )  return 1;	
   	else                  f1 = 1;

    g1 =  1;
    h1 = -1;
    if ( i1>0    && (k1=ia[3*j1+2])<ia[3*j1+3] ) {      /* (i,   j,   n+1) */
    	g1     = 0;
    	h1     = 0;
    	t1[0]  = a[k1];
    }     
    if ( i1<m1-1 && (k1=ia[3*j1+1])<ia[3*j1+2] ) {      /* (i,   j,   n  ) */
    	h1     = 1;
        t1[1]  = a[k1];
    }     
    if ( h1==1 && j1==n1 ) {
    	t1[1] += a[ia[3*n1+5]];                         /* (i+1, n+1, n+1) */
    	
        if ( i1<m1-2 && (k1=ia[3*n1+4])<ia[3*n1+5] ) {  /* (i+1, n+1, n  ) */
        	h1     = 2;
            t1[2]  = a[k1];
        }
    }    		
   		
    g2 =  1;
    h2 = -1;
    if ( i2>0    && (k2=ib[3*j2+2])<ib[3*j2+3] ) {      /* (i,   j,   n+1) */
    	g2     = 0;
    	h2     = 0;
    	t2[0]  = b[k2];
    }     
    if ( i2<m2-1 && (k2=ib[3*j2+1])<ib[3*j2+2] ) {      /* (i,   j,   n  ) */
    	h2     = 1;
        t2[1]  = b[k2];
    }     
    if ( h2==1 && j2==n2 ) {
    	t2[1] += b[ib[3*n2+5]];                         /* (i+1, n+1, n+1) */
    	
        if ( i2<m2-2 && (k2=ib[3*n2+4])<ib[3*n2+5] ) {  /* (i+1, n+1, n  ) */
        	h2     = 2;
            t2[2]  = b[k2];
        }
    }    		
   		
	/* F, (i1+e1,i2+e2) (j1+e1,j2+e2) 
 	                    (k1,   k2   )  */
	for ( e1=0; e1<=f1; e1++ ) 
	for ( e2=0; e2<=f2; e2++ ) {
	    i  = (i1+e1)*m2+i2+e2;
	    
	    for ( k1=ia[3*j1+3*e1]; k1=ia[3*j1+3*e1+1]; k1++ ) 
	    for ( k2=ib[3*j2+3*e2]; k2=ib[3*j2+3*e2+1]; k2++ ) {
	    	
          jv[*q] = off[0] + (i*n1+ja[k1])*n2+jb[k2];
          vv[*q] = a[k1]*b[k2];
          (*q)++;
	    }
	}
	
	/* E2,  (i1+e1, i2+e2 ), (j1+e1, j2+e2 )  
           [(i1+e1, i2    ), (j1+e1, j2    ), (k1, n+1 )],  
           [(i1+e1, i2    ), (j1+e1, j2    ), (k1, n   )]
          if j=n and i+1<=m-1
          =[(i1+e1, i2+1,  i3+e3), (j1+e1, j2+1,  j3+e3), (k1, n+1, k3)],
          if j=n and i+1<=m-2
           [(i1+e1, i2+1,  i3+e3), (j1+e1, j2+1,  j3+e3), (k1, n  , k3)]
    */
	for ( e2=g2; e2<=h2; e2++ )
	for ( e1=0;  e1<=f1; e1++ ) {
        i = (i2-1+e2)*m1+i1+e1;
        		
	    for ( k1=ia[3*j1+3*e1]; k1=ia[3*j1+3*e1+1]; k1++ ) {
	    	
	    	jv[*q] = off[1] + i*n1+ja[k1];
	        vv[*q] = a[k1]*t2[e2];
            (*q)++;		    	
	    }
	}	
	
	
	/* F1,  (i1+e1, i2+e2, i3+e3), (j1+e1, j2+e2, j3+e3)  
           [(i1,    i2+e2, i3+e3), (j1,    j2+e2, j3+e3), (n+1, k2, k3)],  
           [(i1,    i2+e2, i3+e3), (j1,    j2+e2, j3+e3), (n  , k2, k3)]
          if j=n and i+1<=m-1                
          =[(i1+1,  i2+e2, i3+e3), (j1+1,  j2+e2, j3+e3), (n+1, k2, k3)],
          if j=n and i+1<=m-2                                   
           [(i1+1,  i2+2e, i3+e3), (j1+1,  j2+e2, j3+e3), (n  , k2, k3)]
    */
	for ( e1=g1; e1<=h1; e1++ )
	for ( e2=0;  e2<=f2; e2++ )  {
        i = (i1-1+e1)*m2+i2+e2;
        		
	    for ( k2=ib[3*j2+3*e2]; k2=ib[3*j2+3*e2+1]; k2++ )  {
	    	
	    	jv[*q] = off[1] + i*n2+ja[k2];
	        vv[*q] = t1[e1]*b[k2];
            (*q)++;		    	
	    }
	}	
	
	for ( e1=g1; e1<=h1; e1++ )
	for ( e2=g2; e1<=h2; e2++ ) {
        i = (i1-1+e1)*(m2-1)+i2-1+e2;
        
	    jv[*q] = off[1] + i;
	    vv[*q] = t1[e1]*t2[e2];
        (*q)++;		    	
	}		
	
}

void csrExtend(int *ia, int *ja, int **ie)
{ int j,k, n=ia[0]-1;
	
  *ie = (INT *) malloc(3*(n+2)*sizeof(INT));
  
  for ( j=0; j<=n; j++ ) {
  	  for ( k=ia[j+1]; k<ia[j+2] && ja[k]<n-1; k++);
  	  (*ie)[3*j  ] = ia[j+1];
  	  (*ie)[3*j+1] = k;
  	  (*ie)[3*j+2] = (k<ia[j+2]?k+1:k);
  }
  (*ie)[3*n+3] = ia[n+2];	    	
}


void assemble2d(int m1,  int m2,   int n1, int n2,
               int *ia,  int *ja,  R *a,   int *ib, int *jb, R *b,
               int **iv, int **jv, R **vv)
{ int i1,i2, j1,j2, k1,k2;
  
  int p   = (m1*n1-1)*(m2*n2-1);
  int q   = gdof(m1, n1) * gdof(m2, n2);
  
  *iv = (INT *) malloc((p+2)*sizeof(INT));
  *jv = (INT *) malloc(q*sizeof(INT));
  *vv = (R   *) malloc(q*sizeof(R));
  
  n1--;
  n2--;
  
  p = q = 0;
  
  /* (i1,i2), (j1,j2) */
  for ( i1=0; i1<m1; i1++ )
  for ( i2=0; i2<m2; i2++ )
  for ( j1=0; j1<n1; j1++ ) 
  for ( j2=0; j2<n2; j2++ ) {
      RowEntries(m1,m2,   i1,i2,   n1,n2,     j1,j2,
                 ia,ja,a, ib,jb,b, &q,*jv,*vv);
      (*iv)[p++] = q;
  }
      
 
  for ( i2=0; i2<m2-1; i2++ )
  for ( i1=0; i1<m1;   i1++ )
  for ( j1=0; j1<n1;   j1++ ) {
      RowEntries(m1,m2,   i1,i2,   n1,n2,     j1,n2,
                 ia,ja,a, ib,jb,b, &q,*jv,*vv);
      (*iv)[p++] = q;
  }

  for ( i1=0; i1<m1-1; i1++ )
  for ( i2=0; i2<m2;   i2++ )
  for ( j2=0; j2<n2;   j2++ ) {
      RowEntries(m1,m2,   i1,i2,   n1,n2,     n1,j2,
                 ia,ja,a, ib,jb,b, &q,*jv,*vv);
      (*iv)[p++] = q;
  }

  for ( i1=0; i1<m1-1; i1++ )
  for ( i2=0; i2<m2-1; i2++ ) {
      RowEntries(m1,m2,   i1,i2,   n1,n2,     n1,n2,
                 ia,ja,a, ib,jb,b, &q,*jv,*vv);
      (*iv)[p++] = q;
  }


	
}               


int main(int argc, char *argv[])
{
	int i1,i2, j1,j2, k1,k2, 
	    p, m1=2,m2=2, n1=10,n2=10;
	
	if (argc>1) m1 = atoi(argv[1]);
	if (argc>2) m2 = atoi(argv[2]);
	if (argc>3) n1 = atoi(argv[3]);
	if (argc>4) n2 = atoi(argv[4]);
	
	int *ia1,*ia2, *ja1,*ja2, *ib1,*ib2, *jb1,*jb2;
	R   *va1,*va2, *vb1,*vb2;
	
	LegendspMat(n1, &va1, &ia1, &ja1, &vb1, &ib1, &jb1);
	LegendspMat(n2, &va2, &ia2, &ja2, &vb2, &ib2, &jb2);
	
	/* checking B */
	printf("local mass:\n");
	for ( j1=0; j1<ib1[0]; j1++ ) {
	    for ( k1=ib1[j1+1]; k1<ib1[j1+2]; k1++ )
	    printf("[%2d %2d], %f\n",j1,jb1[k1],vb1[k1]); 
	}
	
	printf("\n\n");
	
	/* checking A */
	printf("local stiff:\n");
	for ( j1=0; j1<ia1[0]; j1++ ) {
	    for ( k1=ia1[j1+1]; k1<ia1[j1+2]; k1++ )
	    printf("[%2d %2d], %f\n",j1,ja1[k1],va1[k1]); 
	}

	printf("\n\n");
	/* chcking global indices */
	printf("DoF per row:\n");
	for ( i1=0; i1<m1; i1++ ) 
	    for ( j1=0; j1<(i1==m1-1?n1-1:n1); j1++ ) {   
	    	p = l2gdeg(m1,i1,n1,j1);
	    	printf("[%2d %2d], %d\n",i1,j1,p);
	    }
	    
	printf("\n\n");
	
	int *ea1,*ea2, *eb1,*eb2;
    csrExtend(ia1, ja1, &ea1);
    csrExtend(ib1, jb1, &eb1);
    
    /* chcking extend csr */
    printf("extended csr B:\n");
	for ( j1=0; j1<ib1[0]; j1++ ) {
	    printf("[%2d %3d], [%3d %2d], [%3d %2d], [%3d %2d]\n",j1,ib1[j1+1], 
	            eb1[3*j1],  jb1[eb1[3*j1]],
	            eb1[3*j1+1],jb1[eb1[3*j1+1]], eb1[3*j1+2],  jb1[eb1[3*j1+2]]);
	}
	printf("[%2d %3d], %3d --\n\n",j1,ib1[j1+1],eb1[3*j1]);
    
    printf("\nextended csr A:\n");
	for ( j1=0; j1<ia1[0]; j1++ ) {
	    printf("[%2d %3d], [%3d %2d], [%3d %2d], [%3d %2d]\n",j1,ia1[j1+1], 
	            ea1[3*j1],  ja1[ea1[3*j1]],
	            ea1[3*j1+1],ja1[ea1[3*j1+1]], ea1[3*j1+2],  ja1[ea1[3*j1+2]]);
	}
	printf("[%2d %3d], %3d --\n\n",j1,ia1[j1+1],ea1[3*j1]);
	
	    
	
}

