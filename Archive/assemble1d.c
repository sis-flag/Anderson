#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef int INT;
typedef double R;

/*  0<=i<=m-1, 0<=j<=n-1  without specifying

    i,                            j
    i,   j

    i,                            0            i<m-1
    i,   n 
    i+1, n+1
*/


INT l2gindex(INT m, INT i, INT n, INT j)
{	     
  INT  off[2], p=-1;
  
  n--;
  
  off[0] = 0;
  off[1] = off[0] + m*n;

/* local-to-gobal mapping*/
  if ( j<n ) {
  	 p = off[0] + i*n+j;    
  } else if ( j==n && i<m-1 ) {
  	 p = off[1] + i; 
  } else if ( j==n+1 && i>0 ) {
  	 p = off[1] + i-1; 
  }
 
  return p;
}

INT l2gdeg(INT m, INT i, INT n, INT j) 
{	     
  INT  x;

  if ( (j==n-1 && i==m-1) || (j==n && i==0) ) return -1;
  	
  if ( j==n )           /* left  bnd. mode */
  	 x = 7 - 2 * ((n<3)+(n<2)) - (i==1) - (i==m-1); /* DD BCs*/
  else if ( j==n-1 )    /* right bnd. mode */
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
  INT  i,j, x= m*(3*n+4)-13+2*(m==1)-(2*m-4)*((n<3)+(n<2));
  return x;
  
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


int RowEntries(int m, int i, int n, int j, 
               int *ia, int *ja, R *a,
               int *q,  int *jv, R *vv)
{
  int k, f, g, h, e, off[2];
  R   t[3];

  n--;	
  off[0] = 0;
  off[1] = off[0] + m*n;
   
  if  ( j>n || j<0 ) return 1;
  else if ( j<n    ) f = 0;
  else if ( i==m-1 ) return 1;	
  else               f = 1;   

  g =  1;
  h = -1;
  if ( i>0    && (k=ia[3*j+2])<ia[3*j+3] ) {      /* (i,   j,   n+1) */
  	  g     = 0;
  	  h     = 0;
  	  t[0]  = a[k];
  }     
  if ( i<m-1 && (k=ia[3*j+1])<ia[3*j+2] ) {      /* (i,   j,   n  ) */
  	  h     = 1;
      t[1]  = a[k];
  }     
  if ( h==1 && j==n ) {
  	  t[1] += a[ia[3*n+5]];                      /* (i+1, n+1, n+1) */
  	
      if ( i<m-2 && (k=ia[3*n+4])<ia[3*n+5] ) {  /* (i+1, n+1, n  ) */
      	  h     = 2;
          t[2]  = a[k];
      }
  }    		
    	
    //printf("%d %d %d\n",f,g,h);
	for ( e=0; e<=f; e++ )
	for ( k=ia[3*j+3*e]; k<ia[3*j+3*e+1]; k++ )  {
	    
        jv[*q] = off[0] + (i+e)*n+ja[k];
        vv[*q] = a[k];
        (*q)++;
	}
	
	for ( e=g; e<=h; e++ ) {
        		
	    jv[*q] = off[1] + (i-1+e);
	    vv[*q] = t[e];
        (*q)++;	
    }	    	
   

}               
               

/* n: degree of polynomials */



void assemble1d(int m, int n, int *ia,  int *ja,  R *a,
                int **iv, int **jv, R **vv)
{ int i, j, k;
  
  int p   = m*n-1;
  int q   = gdof(m,n);
  
  *iv = (INT *) malloc((p+2)*sizeof(INT));
  *jv = (INT *) malloc(q*sizeof(INT));
  *vv = (R   *) malloc(q*sizeof(R));
  
  (*iv)[0] = p;
  (*iv)[1] = 0;
  p = 1;
  q = 0;
  
  for ( i=0; i<m;   i++ )
  for ( j=0; j<n-1; j++ ) {
      RowEntries(m, i, n, j, ia,ja,a, &q,*jv,*vv);

      (*iv)[++p] = q;
  }
 
  for ( i=0; i<m-1; i++ ) {
      RowEntries(m, i, n, n-1, ia,ja,a, &q,*jv,*vv);
      (*iv)[++p] = q;
  }

  //printf("[m,n]=[%2d,%2d], nnz = %5d %5d\n",m,n,q,q-gdof(m,n));
	
}               


int main(int argc, char *argv[])
{
	int i, j, k, 
	    p, m=2, n=10;
	
	if (argc>1) m = atoi(argv[1]);
	if (argc>2) n = atoi(argv[2]);
	
	int *ia, *ja, *ib, *jb;
	R   *va, *vb;
	
	LegendspMat(n, &va, &ia, &ja, &vb, &ib, &jb);
	
//	/* checking B */
//	printf("local mass:\n");
//	for ( j=0; j<ib[0]; j++ ) {
//	    for ( k=ib[j+1]; k<ib[j+2]; k++ )
//	    printf("[%2d %2d], %f\n",j,jb[k],vb[k]); 
//	}
//	
//	printf("\n\n");
//	
//	/* checking A */
//	printf("local stiff:\n");
//	for ( j=0; j<ia[0]; j++ ) {
//	    for ( k=ia[j+1]; k<ia[j+2]; k++ )
//	    printf("[%2d %2d], %f\n",j,ja[k],va[k]); 
//	}
//
//	printf("\n\n");
//	/* chcking global indices */
//	printf("DoF per row:\n");
//	for ( i=0; i<m; i++ ) 
//	    for ( j=0; j<(i==m-1?n-1:n); j++ ) {   
//	    	p = l2gdeg(m,i,n,j);
//	    	printf("[%2d %2d], %d\n",i,j,p);
//	    }
//	    
//	printf("\n\n");
	
	int *ea, *eb;
    csrExtend(ia, ja, &ea);
    csrExtend(ib, jb, &eb);
    
//    /* chcking extend csr */
//    printf("extended csr B:\n");
//    printf("DoF = %d\n",ib[0]);
//	for ( j=0; j<ib[0]; j++ ) {
//	    printf("[%2d %3d], [%3d %2d], [%3d %2d], [%3d %2d]\n",j,ib[j+1], 
//	            eb[3*j],  jb[eb[3*j]],
//	            eb[3*j+1],jb[eb[3*j+1]], eb[3*j+2],  jb[eb[3*j+2]]);
//	}
//	printf("[%2d %3d], %3d --\n\n",j,ib[j+1],eb[3*j]);
//    
//    printf("\nextended csr A:\n");
//    printf("DoF = %d\n",ia[0]);
//	for ( j=0; j<ia[0]; j++ ) {
//	    printf("[%2d %3d], [%3d %2d], [%3d %2d], [%3d %2d]\n",j,ia[j+1], 
//	            ea[3*j],  ja[ea[3*j]],
//	            ea[3*j+1],ja[ea[3*j+1]], ea[3*j+2],  ja[ea[3*j+2]]);
//	}
//	printf("[%2d %3d], %3d --\n\n",j,ia[j+1],ea[3*j]);
	
	
	int *iv, *jv;
	R *vv;
	//printf("\n\nglobal matrix A:\n");
//	assemble1d(m, n, ea, ja, va, &iv, &jv, &vv);
//	
//	FILE *fp=fopen("stiff","w+");
//	for ( j=0; j<iv[0]; j++ ) 
//	for ( k=iv[j+1]; k<iv[j+2]; k++ )
//	    fprintf(fp,"%4d %4d %.16g\n",j,jv[k],vv[k]);
//	fclose(fp);
//		
//    free(iv);	
//    free(jv);	
//    free(vv);	
	
	//printf("\n\nglobal matrix B:\n");
	assemble1d(m, n, eb, jb, vb, &iv, &jv, &vv);
	
	FILE *fd=fopen("mass","w+");
	for ( j=0; j<iv[0]; j++ ) 
	for ( k=iv[j+1]; k<iv[j+2]; k++ )
	    fprintf(fd,"%4d %4d %.16g\n",j,jv[k],vv[k]);
    fclose(fd);
}

