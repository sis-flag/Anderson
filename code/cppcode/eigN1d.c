static char help[] =
    "Anderson localization eigenvalue problem with Neumann boundary in 1dimension.\n"
    "Legendre spectral elemement method.\n"
    "The command line options are:\n"
    "    -N <n>, where <n> = degree of polynomials.\n"
    "    -M <m>, where <m> = number of grid subdivisions in each process.\n\n";

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
    PetscOptionsGetInt(NULL, NULL, "-N", &N, NULL);
    PetscOptionsGetInt(NULL, NULL, "-M", &M, NULL);

    Real h = 1. / (Real)(M); // length of interval

    Int np, id;
    MPI_Comm_size(PETSC_COMM_WORLD, &np);
    MPI_Comm_rank(PETSC_COMM_WORLD, &id);

    PetscPrintf(PETSC_COMM_WORLD, "\n1-D Laplacian Eigenproblem,\n proc # = %D subint # = %D pol deg = %D\n\n", np, M, N);

    Mat A, B;
    MatCreate(PETSC_COMM_WORLD, &A);
    MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, M * N + 1, M * N + 1);
    MatSetFromOptions(A);
    MatSetUp(A);

    MatCreate(PETSC_COMM_WORLD, &B);
    MatSetSizes(B, PETSC_DECIDE, PETSC_DECIDE, M * N + 1, M * N + 1);
    MatSetFromOptions(B);
    MatSetUp(B);

    Int *iAhat, *jAhat, *iBhat, *jBhat, *iFhat;
    Real *vAhat, *vBhat, *vFhat;
    getAhat(N, &vAhat, &iAhat, &jAhat);
    getBhat(N, &vBhat, &iBhat, &jBhat);
    getAhat(N, &vFhat, &iFhat);

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

            MatSetValues(A, 1, &p, nz, col, val, ADD_VALUES);

            for (k = ib[j + 1], nz = 0; k < ib[j + 2]; k++)
            {
                q = p0 + (jb[k] >= 2 ? jb[k] - 2 : jb[k] * n - 1);
                if (q < 0 || q >= N)
                    continue; /* zero Dirichlet boundary conditions */

                col[nz] = q;
                val[nz++] = vb[k] * a;
            }

            MatSetValues(B, 1, &p, nz, col, val, ADD_VALUES);
        }
    }

    MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

    MatAssemblyBegin(B, MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(B, MAT_FINAL_ASSEMBLY);

    PetscPrintf(PETSC_COMM_WORLD, "Matrix A\n");
    MatView(A, PETSC_VIEWER_STDOUT_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "\n");

    PetscPrintf(PETSC_COMM_WORLD, "Matrix B\n");
    MatView(B, PETSC_VIEWER_STDOUT_WORLD);
    PetscPrintf(PETSC_COMM_WORLD, "\n");

    Vec xr, xi;
    MatCreateVecs(A, NULL, &xr);
    MatCreateVecs(A, NULL, &xi);

    EPS eps;
    EPSCreate(PETSC_COMM_WORLD, &eps);
    EPSSetOperators(eps, B, A);
    EPSSetProblemType(eps, EPS_GHEP);
    EPSSetFromOptions(eps);
    EPSSolve(eps);

    // Get some information from the solver and display it
    Int its;
    EPSGetIterationNumber(eps, &its);
    PetscPrintf(PETSC_COMM_WORLD, " Number of iterations of the method: %D\n", its);

    EPSType type;
    EPSGetType(eps, &type);
    PetscPrintf(PETSC_COMM_WORLD, " Solution method: %s\n\n", type);

    Int nev;
    EPSGetDimensions(eps, &nev, NULL, NULL);
    PetscPrintf(PETSC_COMM_WORLD, " Number of requested eigenvalues: %D\n", nev);

    Int maxit;
    EPSGetTolerances(eps, &tol, &maxit);
    PetscPrintf(PETSC_COMM_WORLD, " Stopping condition: tol=%.4g, maxit=%D\n", (double)tol, maxit);

    Int nconv;
    EPSGetConverged(eps, &nconv);
    PetscPrintf(PETSC_COMM_WORLD, " Number of converged eigenpairs: %D\n\n", nconv);

    // Free work space
    EPSDestroy(&eps);
    MatDestroy(&A);
    MatDestroy(&B);
    VecDestroy(&xr);
    VecDestroy(&xi);
    SlepcFinalize();
}