#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>

typedef double R;
typedef int INT;

void LegendMat(INT N, R **a, INT **ia, INT **ja, R **b, INT **ib, INT **jb);

int main(int argc, char *argv[])
{
INT N = 1, k, n;

INT *ia, *ja, *ib, *jb;
R   *va, *vb;

if (argc>1) N = atoi(argv[1]);
	
		
 LegendMat(N, &va, &ia, &ja, &vb, &ib, &jb);
 
 
 for (n=0; n<ia[0]; n++) {
 	for (k=ia[n+1]; k<ia[n+2]; k++) {
 	    printf("%4d %4d %lf\n",n,ja[k],va[k]);	
 	}
 }
  printf("nnz: %d (%d)\n\n",k,N+3);
	
 for (n=0; n<ib[0]; n++) {
 	for (k=ib[n+1]; k<ib[n+2]; k++) {
 	    printf("%4d %4d %lf\n",n,jb[k],vb[k]);	
 	}
 }
  printf("nnz: %d (%d)\n\n",k,3*N+5);
  
 free(va); free(ia); free(ja); 
 free(vb); free(ib); free(jb); 
	
	return 0;
}


void LegendMat(INT N, R **va, INT **ia, INT **ja, R **vb, INT **ib, INT **jb)
{ INT n, k;

/* basis 
    phi[k] = (L[k]-L[k-2])/sqrt(|4*k-2|)  2<=k<=N
           = (1-x)/2                      k = 0
           = (1+x)/2                      k = 1
*/

  *va = malloc((N+3)*sizeof(R));
  *ia = malloc((N+3)*sizeof(int));
  *ja = malloc((N+3)*sizeof(int));
 
   n = 0;
  (*ia)[0]   = N+1;
  
  (*ia)[1]   = n;
  (*ja)[n]   = 0;
  (*va)[n++] = .5;
  (*ja)[n]   = 1;
  (*va)[n++] = -.5;
  
  (*ia)[2]   = n;
  (*ja)[n]   = 0;
  (*va)[n++] = -.5;
  (*ja)[n]   = 1;
  (*va)[n++] = .5;
  
  for ( k=2; k<=N; k++ ) {
  		(*ia)[k+1] = n;
  		(*ja)[n]   = k;
  		(*va)[n++]  = 1.;
  }
  (*ia)[N+2] = n;
  
  
  *vb = malloc((3*N+5)*sizeof(R)*2);
  *ib = malloc((N+3)*sizeof(R)*2);
  *jb = malloc((3*N+5)*sizeof(R)*2);
  
  n = 0;
  (*ib)[0]   = N+1;
  
  (*ib)[1]   = 0;
  (*jb)[n]   = 0;
  (*vb)[n++] = 2./3.;
  (*jb)[n]   = 1;
  (*vb)[n++] = 1./3.;
  if (2<=N) { 
		(*jb)[n]   = 2;
  	(*vb)[n++] = -1./sqrt(6.);
  }
  if (3<=N) {
  	(*jb)[n]   = 3;
  	(*vb)[n++] = 1./sqrt(90.);
  }
  
  (*ib)[2]   = n;
  (*jb)[n]   = 0;
  (*vb)[n++] = 1./3.;
  (*jb)[n]   = 1;
  (*vb)[n++] = 2./3.;
  if (2<=N) {
  	(*jb)[n]   = 2;
  	(*vb)[n++] = -1./sqrt(6.);
  }
  if (3<=N) {
  	(*jb)[n]  = 3;
  	(*vb)[n++] = -1./sqrt(90.);
  }
  
  if (2<=N) {
  	(*ib)[3]   = n;
  	(*jb)[n]   = 0;
  	(*vb)[n++] = -1./sqrt(6.);
  	(*jb)[n]   = 1;
  	(*vb)[n++] = -1./sqrt(6.);
  	(*jb)[n]   = 2;
  	(*vb)[n++] = 0.4;
  	if (4<=N) {
		(*jb)[n]   = 4;
  		(*vb)[n++] = -sqrt(1./525.);
  	}	
  }
  
  if (3<=N) {
  	(*ib)[4]   = n;
  	(*jb)[n]   = 0;
  	(*vb)[n++] = 1./sqrt(90.);
  	(*jb)[n]   = 1;
  	(*vb)[n++] = -1./sqrt(90.);
  	(*jb)[n]   = 3;
  	(*vb)[n++] = 2./21.;
  	if (5<=N) {
		(*jb)[n]   = 5;
  		(*vb)[n++] = -sqrt(1./2205.);
  	}	
  }  		
  
  for (k=4; k<=N; k++) {
      (*ib)[k+1] = n;
      (*jb)[n]   = k-2;
      (*vb)[n++] = -.25 / ( sqrt( (k-2.5) * (k-.5) ) * (k-1.5) );
      (*jb)[n]   = k;
      (*vb)[n++] = ( .25/(k+.5) + .25/(k-1.5) ) / (k-0.5);
      if ( k+2 <= N ) {
         (*jb)[n]   = k+2;
         (*vb)[n++] = -.25 / ( sqrt( (k-.5) * (k+1.5) ) * (k+0.5) );
      }
  }
  (*ib)[N+2] = n;
  		
 return;
  
} 