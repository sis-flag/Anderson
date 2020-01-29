

static char help[] = "Legendre spectral elemement method for Laplacian eigenvalues in 2 dimension.\n\n"
  "The command line options are:\n"
  "  -n <n>, where <n> = number of grid subdivisions = matrix dimension.\n\n";

#include <slepceps.h>
#define PI 3.14159265358979323846264338328

typedef PetscReal R;
typedef PetscInt  INT;

void LegendMat(INT N, R **va, INT **ia, INT **ja, R **vb, INT **ib, INT **jb);

#undef __FUNCT__
#define __FUNCT__ "main"
int main(int argc,char **argv)
{
  Mat            A,B;           /* problem matrix */
  EPS            eps;         /* eigenproblem solver context */
  EPSType        type;
  PetscReal      error,tol,re,im;
  PetscScalar    kr,ki;
  Vec            xr,xi;
  PetscInt       n=8,i,j,k,p0,p,nev,maxit,its,nconv;
  PetscErrorCode ierr;
  PetscInt       m=1, np, id, N, *ia, *ja, *ib, *jb, nz, q, col[64];  /* m, number of sub-intervals on each proc*/
  PetscReal      a, *va, *vb, val[64]; /*a: length of interval */
  
     
  SlepcInitialize(&argc,&argv,(char*)0,help);

  ierr = PetscOptionsGetInt(NULL,NULL,"-n",&n,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetInt(NULL,NULL,"-m",&m,NULL);CHKERRQ(ierr);
/*ierr = PetscOptionsGetReal(NULL,NULL,"-a",&a,NULL);CHKERRQ(ierr);
  ierr = PetscOptionsGetIntArray(NULL,NULL,"-n",deg,dmax,NULL);CHKERRQ(ierr);  
  ierr = PetscOptionsGetRealArray(NULL,NULL,"-a",intval,imax,NULL);CHKERRQ(ierr);  
  */

  MPI_Comm_size(PETSC_COMM_WORLD,&np);
  MPI_Comm_rank(PETSC_COMM_WORLD,&id);
  
  np1 = (int) sqrt((R) np); 
  np2 = np1; 
  a1  = 1./(R)(2*np1*m);
  a2  = 1./(R)(2*np2*m);
  
  N1  = np1*m*n-1;
  N2  = np2*m*n-1;
  N   = N1*N2;
  
  id1 = id/np2;
  id2 = id%np2;
  
  M1 = (id1==np1-1?m*n-1:m*n);
  M2 = (id2==np2-1?m*n-1:m*n);
  M  = M1*M2;	

    
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n2-D Laplacian Eigenproblem,\n proc # = %D x %D;\nsubdomains per proc = %D x %D; pol deg = %D x %D\n\n",np1,np2,m,m,n,n);CHKERRQ(ierr);   
   
  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = MatCreate(PETSC_COMM_WORLD,&A);CHKERRQ(ierr);
  ierr = MatSetSizes(A,M,M,PETSC_DETERMINE,PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);
  
  ierr = MatGetSize(A,&j,&k);CHKERRQ(ierr);  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"A: %d x %d\n",j,k);CHKERRQ(ierr);   
  
  ierr = MatCreate(PETSC_COMM_WORLD,&B);CHKERRQ(ierr);
  ierr = MatSetSizes(B,M,M,PETSC_DETERMINE,PETSC_DETERMINE);CHKERRQ(ierr);
  ierr = MatSetFromOptions(B);CHKERRQ(ierr);
  ierr = MatSetUp(B);CHKERRQ(ierr);
  
  ierr = MatGetSize(B,&j,&k);CHKERRQ(ierr);  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"B: %d x %d\n",j,k);CHKERRQ(ierr);    
  
  LegendMat(n, &va, &ia, &ja, &vb, &ib, &jb);

   /* local to global mapping 
   
      local index:
      proc #   local subinterval #   local bais #
      (id,     i,                    j         )
      
      global index:
      g1 = 
      id1*m1*n1 + i1*n1  -  1,         j1 = 0
      id1*m1*n1 + i1*n1  +  n1-1,      j1 = 1
      id1*m1*n1 + i1*n1  +  j1 - 2,    2 <= j1 <= n1 
      
      g2 = 
      id2*m2*n2 + i2*n2  -  1,         j2 = 0
      id2*m2*n2 + i2*n2  +  n2-1,      j2 = 1
      id2*m2*n2 + i2*n2  +  j2 - 2,    2 <= j2 <= n2 
      
      g = g1*N2 + g2
      
   */
   
   
for (i1=0; i1<m1; i1++)  	{
  	p10 = id1*m1*n1 + i1*n1;
  	for (j1=0; j1<=n1; j1++) {
  	   	p1 = p10 + (j1>=2? j1-2 : j1*n1-1); /*local to global mapping */
  	   	if (p1<0 || p1>=N1) continue; /* zero Dirichlet boundary conditions */
  	   	 	
		for (i2=0; i2<m2; i2++)  	{
  	   		p20 = id2*m2*n2 + i2*n2;
  	   		for (j2=0; j2<=n2; j2++) {
  	   	 	p2 = p20 + (j2>=2? j2-2 : j2*n2-1); /*local to global mapping */
  	   	 	if (p2<0 || p2>=N2) continue; /* zero Dirichlet boundary conditions */
  	   	 	 
  	   	 	for (k=ia[j+1],nz=0; k<ia[j+2]; k++) {
  	   	 	  	q = p0 + (ja[k]>=2 ? ja[k]-2 : ja[k]*n-1);
		  	    if (q<0 || q>=N) continue; /* zero Dirichlet boundary conditions */
		  	     	
  	   	 	  	col[nz]   = q;
  	   	 	 	val[nz++] = va[k]/a;
  	   	 	}

	  	    ierr = MatSetValues(A,1, &p, nz, col, val, ADD_VALUES); CHKERRQ(ierr);

  	   	 for (k=ib[j+1],nz=0; k<ib[j+2]; k++) {
  	   	 	  q = p0 + (jb[k]>=2 ? jb[k]-2 : jb[k]*n-1);
		  	     if (q<0 || q>=N) continue; /* zero Dirichlet boundary conditions */
		  	     	
  	   	 	  col[nz]   = q;
  	   	 	  val[nz++] = vb[k]*a;
  	   	 }

	  	    ierr = MatSetValues(B,1, &p, nz, col, val, ADD_VALUES); CHKERRQ(ierr);
  	   	
  	   }
  	   
  }
  
  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyBegin(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(B,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  
/*  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Matrix A\n");CHKERRQ(ierr);
  ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
   
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Matrix B\n");CHKERRQ(ierr);
  ierr = MatView(B,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
*/  

  ierr = MatCreateVecs(A,NULL,&xr);CHKERRQ(ierr);
  ierr = MatCreateVecs(A,NULL,&xi);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Create eigensolver context
  */
  ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);

  /*
     Set operators. In this case, it is a standard eigenvalue problem
  */
  ierr = EPSSetOperators(eps,B,A);CHKERRQ(ierr);
  ierr = EPSSetProblemType(eps,EPS_GHEP); CHKERRQ(ierr);

  /*
     Set solver parameters at runtime
  */
  ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

  ierr = EPSSolve(eps);CHKERRQ(ierr);
  /*
     Optional: Get some information from the solver and display it
  */
  ierr = EPSGetIterationNumber(eps,&its);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of iterations of the method: %D\n",its);CHKERRQ(ierr);
  ierr = EPSGetType(eps,&type);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Solution method: %s\n\n",type);CHKERRQ(ierr);
  ierr = EPSGetDimensions(eps,&nev,NULL,NULL);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of requested eigenvalues: %D\n",nev);CHKERRQ(ierr);
  ierr = EPSGetTolerances(eps,&tol,&maxit);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Stopping condition: tol=%.4g, maxit=%D\n",(double)tol,maxit);CHKERRQ(ierr);

  /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
  /*
     Get number of converged approximate eigenpairs
  */
  ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD," Number of converged eigenpairs: %D\n\n",nconv);CHKERRQ(ierr);

  if (nconv>0) {
    /*
       Display eigenvalues and relative errors
    */
    ierr = PetscPrintf(PETSC_COMM_WORLD,
         "           k          ||Ax-kx||/||kx||\n"
         "   ----------------- ------------------\n");CHKERRQ(ierr);

    for (i=0;i<nconv;i++) {
      /*
        Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and
        ki (imaginary part)
      */
      ierr = EPSGetEigenpair(eps,i,&kr,&ki,xr,xi);CHKERRQ(ierr);
      /*
         Compute the relative error associated to each eigenpair
      */
      ierr = EPSComputeError(eps,i,EPS_ERROR_RELATIVE,&error);CHKERRQ(ierr);

#if defined(PETSC_USE_COMPLEX)
      re = PetscRealPart(kr);
      im = PetscImaginaryPart(kr);
#else
      re = kr;
      im = ki;
#endif
      if (im!=0.0) {
        ierr = PetscPrintf(PETSC_COMM_WORLD," %9f%+9fi %12g\n",(double)re,(double)im,(double)error);CHKERRQ(ierr);
      } else {
        ierr = PetscPrintf(PETSC_COMM_WORLD,"   %12f       %12g\n",1./((double)re*PI*PI),(double)error);CHKERRQ(ierr);
      }
    }
    ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
  }

  /*
     Free work space
  */
  ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = MatDestroy(&B);CHKERRQ(ierr);
  ierr = VecDestroy(&xr);CHKERRQ(ierr);
  ierr = VecDestroy(&xi);CHKERRQ(ierr);
  ierr = SlepcFinalize();
  return ierr;
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



