#define PI 3.14159265358979323846264338328

typedef double R;
typedef int  INT;

void LegendMat(INT N, R **va, INT **ia, INT **ja, R **vb, INT **ib, INT **jb);


{
INT *ia, *ib, *ja, *jb;
R   *va, *vb;
	
LegendMat(n, &va, &ia, &ja, &vb, &ib, &jb);

for (j=2; j<=n/2; j++) {
  for (k=ib[j]; k<ib[j+1]; k++) {
  	i = (j/2)+(j%2)*(n/2);
  	if (jb[k]==j)   b[i]   = vb[k];
    if (jb[k]==j+2) b[i+n] = vb[k];
  }
}

/* b[n/2+n] = 0.0; */

j = n/2;
i = j*j;
dsteqr_("I", &j, b,   b+n,   z,	  &j, w, &info);
k = (n-1)/2;
dsteqr_("I", &k, b+j, b+n+j, z+i, &k, w, &info);
	
}




DSolve(INT n, INT many, R epsilon, R lambda, R *B, R *h, R *V, R *E, R *DP, R *OP)
{
	
alp = 1.;
bet = 0.;	
no  = n/2;
ne  = n-no;

m   = n*n*n;

eps[0] = epsilon;
lam[0] = lambda;

/* third dimension: Fortran Style */
M = n*n;
for ( j=0; j<many; j++ ) {
	dgemm_("N", "N", &M, ne, ne, &alp, B+M*n*j,    &M,    
	        V, &ne, &bet, C+M*n*j, &M);
	        
	dgemm_("N", "N", &M, no, no, &alp, B+M*n*j+ne, &M,    
	        V+ne*ne, &no, &bet, C+M*n*j+ne, &M);	        
}
memcopy(B,C,sizeof(R)*n*n*n*many);

/* second dimension: Fortran Style */
M = n;
for ( j=0; j<many*n; j++ ) {
	dgemm_("N", "N", &M, ne, ne, &alp, B+M*n*j,    &M,    
	        V, &ne, &bet, C+M*n*j, &M);
	        
	dgemm_("N", "N", &M, no, no, &alp, B+M*n*j+ne, &M,    
	        V+ne*ne, &no, &bet, C+M*n*j+ne, &M);	        
}
memcopy(B,C,sizeof(R)*n*n*n*many);

 	
for ( d=0, M=m/n, N=1; d<dim-1; d++, M/=n, N*=n) {
     
    for ( j=0; j<many; j++) {
    	for ( k=0; k<N; k++ ) {
   	
   			dgemm_("N", "N", &M, ne, ne, &alp, B+M*(n*(j*N+k)),    &M, V,      &ne, &bet, C+M*(n*k),    &M);
   		
   			dgemm_("N", "N", &M, no, no, &alp, B+M*(n*(j*N+k)+ne), &M, V+ne*ne, &no, &bet, C+M*(n*k+ne), &M);
    	}		
    	
   		memcopy(B+m*j,C,sizeof(R)*m);
   	}
   	
   	
    for ( d=0; d<dim-1; d++) {
    	for ( k=0; k<n; k++) {
    	    	
    	}
     
    for ( j=0; j<many; j++) {
    	for ( k=0; k<N; k++ ) {
   	        
   	
   	
	
   	N *= n;
	
}

	
	
if (dim==1) {
   INT q = n/2;
   INT p = n-q;
   INT L = m*n;
   R *E  = (R *) malloc(n*sizeof(R)*2);
   R *F  = E + n;
   
   for (j=0; j<m; ++) {
   	    for (k=0; k<n; k++) {
   	    	E[k] = eps[j]*2./h[0] + lam[j]*h[0]/2*DP[j*n+k];
   	    	F[k] = lam[j]*h[0]/2*OP[j*n+k];
   	    }
   	
   	    dpttrf_( &p, E,   F,   &info );
   	    dpttrf_( &q, E+p, F+p, &info );
   	    
   	    dpttrs_( &p, &howmany, E,   F,   B+j*n,   &L, &info );
   	    dpttrs_( &q, &howmany, E+p, F+p, B+j*n+p, &L, &info );   	
   }
   
   free(E);
   return;
	
} else {	
 	



dgemm_("N", "N", &m, &n, &n, alp, F, &m, V, &n, &bet, C, &m);
	
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