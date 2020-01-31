static char help[] = 
    "Anderson localization eigenvalue problem with Neumann boundary in 1dimension.\n"
    "Legendre spectral elemement method.\n"
    "The command line options are:\n"
    "    -N <n>, where <n> = degree of polynomials.\n"
    "    -M <m>, where <m> = number of grid subdivisions.\n\n";

// #include <slepceps.h>

#include "lgmat.h"

#define PI 3.14159265358979323846264338328

typedef PetscReal Real;
typedef PetscInt INT;

int main(int argc, char **argv)
{
    // Mat A, B; /* problem matrix */
    // EPS eps;  /* eigenproblem solver context */
    // EPSType type;
    // PetscReal error, tol, re, im;
    // PetscScalar kr, ki;
    // Vec xr, xi;
    // PetscInt n = 8, i, j, k, p0, p, nev, maxit, its, nconv;
    // PetscErrorCode ierr;
    // PetscInt m = 1, np, id, N, *ia, *ja, *ib, *jb, nz, q, col[64]; /* m, number of sub-intervals on each proc*/
    // PetscReal a, *va, *vb, val[64];                                /*a: length of interval */

    SlepcInitialize(&argc, &argv, (char *)0, help);

    Int M, N;
    ierr = PetscOptionsGetInt(NULL, NULL, "-N", &N, NULL); CHKERRQ(ierr);
    ierr = PetscOptionsGetInt(NULL, NULL, "-M", &M, NULL); CHKERRQ(ierr);

    int np, id;
    MPI_Comm_size(PETSC_COMM_WORLD, &np);
    MPI_Comm_rank(PETSC_COMM_WORLD, &id);

    a = 1. / (R)(2 * np * m);
    N = np * m * n - 1;

    ierr = PetscPrintf(PETSC_COMM_WORLD, "\n1-D Laplacian Eigenproblem,\n proc # = %D subint # = %D pol deg = %D\n\n", np, m, n);
    CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     Compute the operator matrix that defines the eigensystem, Ax=kx
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    ierr = MatCreate(PETSC_COMM_WORLD, &A);
    CHKERRQ(ierr);
    ierr = MatSetSizes(A, id == np - 1 ? m * n - 1 : m * n, id == np - 1 ? m * n - 1 : m * n, PETSC_DETERMINE, PETSC_DETERMINE);
    CHKERRQ(ierr);
    ierr = MatSetFromOptions(A);
    CHKERRQ(ierr);
    ierr = MatSetUp(A);
    CHKERRQ(ierr);

    ierr = MatGetSize(A, &j, &k);
    CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "A: %d x %d\n", j, k);
    CHKERRQ(ierr);

    ierr = MatCreate(PETSC_COMM_WORLD, &B);
    CHKERRQ(ierr);
    ierr = MatSetSizes(B, id == np - 1 ? m * n - 1 : m * n, id == np - 1 ? m * n - 1 : m * n, PETSC_DETERMINE, PETSC_DETERMINE);
    CHKERRQ(ierr);
    ierr = MatSetFromOptions(B);
    CHKERRQ(ierr);
    ierr = MatSetUp(B);
    CHKERRQ(ierr);

    ierr = MatGetSize(B, &j, &k);
    CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, "B: %d x %d\n", j, k);
    CHKERRQ(ierr);

    LegendMat(n, &va, &ia, &ja, &vb, &ib, &jb);

    /* local to global mapping 
   
      local index:
      proc #   local subinterval #   local bais #
      (id,     i,                    j         )
      
      global index:
      
      id*m*n + i*n  -  1,      j = 0
      id*m*n + i*n  +  n-1,    j = 1
      id*m*n + i*n  +  j - 2,  2 <= j <= n 
   */
    for (i = 0; i < m; i++)
    {
        p0 = id * m * n + i * n;

        for (j = 0; j <= n; j++)
        {

            p = p0 + (j >= 2 ? j - 2 : j * n - 1); /*local to global mapping */

            if (p < 0 || p >= N)
                continue; /* zero Dirichlet boundary conditions */

            for (k = ia[j + 1], nz = 0; k < ia[j + 2]; k++)
            {
                q = p0 + (ja[k] >= 2 ? ja[k] - 2 : ja[k] * n - 1);
                if (q < 0 || q >= N)
                    continue; /* zero Dirichlet boundary conditions */

                col[nz] = q;
                val[nz++] = va[k] / a;
            }

            ierr = MatSetValues(A, 1, &p, nz, col, val, ADD_VALUES);
            CHKERRQ(ierr);

            for (k = ib[j + 1], nz = 0; k < ib[j + 2]; k++)
            {
                q = p0 + (jb[k] >= 2 ? jb[k] - 2 : jb[k] * n - 1);
                if (q < 0 || q >= N)
                    continue; /* zero Dirichlet boundary conditions */

                col[nz] = q;
                val[nz++] = vb[k] * a;
            }

            ierr = MatSetValues(B, 1, &p, nz, col, val, ADD_VALUES);
            CHKERRQ(ierr);
        }
    }

    ierr = MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);
    ierr = MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);
    CHKERRQ(ierr);

    /*  
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Matrix A\n");CHKERRQ(ierr);
  ierr = MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
   
  ierr = PetscPrintf(PETSC_COMM_WORLD,"Matrix B\n");CHKERRQ(ierr);
  ierr = MatView(B,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  ierr = PetscPrintf(PETSC_COMM_WORLD,"\n");CHKERRQ(ierr);
*/

    ierr = MatCreateVecs(A, NULL, &xr);
    CHKERRQ(ierr);
    ierr = MatCreateVecs(A, NULL, &xi);
    CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                Create the eigensolver and set various options
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /*
     Create eigensolver context
  */
    ierr = EPSCreate(PETSC_COMM_WORLD, &eps);
    CHKERRQ(ierr);

    /*
     Set operators. In this case, it is a standard eigenvalue problem
  */
    ierr = EPSSetOperators(eps, B, A);
    CHKERRQ(ierr);
    ierr = EPSSetProblemType(eps, EPS_GHEP);
    CHKERRQ(ierr);

    /*
     Set solver parameters at runtime
  */
    ierr = EPSSetFromOptions(eps);
    CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                      Solve the eigensystem
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */

    ierr = EPSSolve(eps);
    CHKERRQ(ierr);
    /*
     Optional: Get some information from the solver and display it
  */
    ierr = EPSGetIterationNumber(eps, &its);
    CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, " Number of iterations of the method: %D\n", its);
    CHKERRQ(ierr);
    ierr = EPSGetType(eps, &type);
    CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, " Solution method: %s\n\n", type);
    CHKERRQ(ierr);
    ierr = EPSGetDimensions(eps, &nev, NULL, NULL);
    CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, " Number of requested eigenvalues: %D\n", nev);
    CHKERRQ(ierr);
    ierr = EPSGetTolerances(eps, &tol, &maxit);
    CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, " Stopping condition: tol=%.4g, maxit=%D\n", (double)tol, maxit);
    CHKERRQ(ierr);

    /* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
                    Display solution and clean up
     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - */
    /*
     Get number of converged approximate eigenpairs
  */
    ierr = EPSGetConverged(eps, &nconv);
    CHKERRQ(ierr);
    ierr = PetscPrintf(PETSC_COMM_WORLD, " Number of converged eigenpairs: %D\n\n", nconv);
    CHKERRQ(ierr);

    if (nconv > 0)
    {
        /*
       Display eigenvalues and relative errors
    */
        ierr = PetscPrintf(PETSC_COMM_WORLD,
                           "           k          ||Ax-kx||/||kx||\n"
                           "   ----------------- ------------------\n");
        CHKERRQ(ierr);

        for (i = 0; i < nconv; i++)
        {
            /*
        Get converged eigenpairs: i-th eigenvalue is stored in kr (real part) and
        ki (imaginary part)
      */
            ierr = EPSGetEigenpair(eps, i, &kr, &ki, xr, xi);
            CHKERRQ(ierr);
            /*
         Compute the relative error associated to each eigenpair
      */
            ierr = EPSComputeError(eps, i, EPS_ERROR_RELATIVE, &error);
            CHKERRQ(ierr);

#if defined(PETSC_USE_COMPLEX)
            re = PetscRealPart(kr);
            im = PetscImaginaryPart(kr);
#else
            re = kr;
            im = ki;
#endif
            if (im != 0.0)
            {
                ierr = PetscPrintf(PETSC_COMM_WORLD, " %9f%+9fi %12g\n", (double)re, (double)im, (double)error);
                CHKERRQ(ierr);
            }
            else
            {
                ierr = PetscPrintf(PETSC_COMM_WORLD, "   %12f       %12g\n", 1. / ((double)re * PI * PI), (double)error);
                CHKERRQ(ierr);
            }
        }
        ierr = PetscPrintf(PETSC_COMM_WORLD, "\n");
        CHKERRQ(ierr);
    }

    /*
     Free work space
  */
    ierr = EPSDestroy(&eps);
    CHKERRQ(ierr);
    ierr = MatDestroy(&A);
    CHKERRQ(ierr);
    ierr = MatDestroy(&B);
    CHKERRQ(ierr);
    ierr = VecDestroy(&xr);
    CHKERRQ(ierr);
    ierr = VecDestroy(&xi);
    CHKERRQ(ierr);
    ierr = SlepcFinalize();
    return ierr;
}
