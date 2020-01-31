/*
basic matrices for Legendre spectral elemement method

basis:
    phi[k] = (1-x)/2                        k = 0
           = (L[k+1]-L[k-1])/sqrt(|4*k+2|)  1<=k<=N-1
           = (1+x)/2                        k = N
           
A[i,j] is inner product (phi'[i], phi'[j])
B[i,j] is inner product (phi[i], phi[j])
*/

#include <slepceps.h>

typedef PetscReal Real;
typedef PetscInt Int;

void getAhat(Int N, Real **vA, Int **iA, Int **jA)
{
    *vA = malloc((N + 3) * sizeof(Real));
    *iA = malloc((N + 3) * sizeof(Int));
    *jA = malloc((N + 3) * sizeof(Int));

    Int n = 0;

    (*iA)[n] = 0;
    (*jA)[n] = 0;
    (*vA)[n] = 0.5;
    n++;
    (*iA)[n] = N;
    (*jA)[n] = N;
    (*vA)[n] = 0.5;
    n++;
    (*iA)[n] = 0;
    (*jA)[n] = N;
    (*vA)[n] = -0.5;
    n++;
    (*iA)[n] = N;
    (*jA)[n] = 0;
    (*vA)[n] = -0.5;
    n++;

    for (Int k = 1; k < N; k++)
    {
        (*iA)[n] = k;
        (*jA)[n] = k;
        (*vA)[n] = 1.0;
        n++;
    }
}

void getBhat(Int N, Real **vB, Int **iB, Int **jB)
{
    *vB = malloc((3 * N + 5) * sizeof(Real));
    *iB = malloc((3 * N + 5) * sizeof(Int));
    *jB = malloc((3 * N + 5) * sizeof(Int));

    Int n = 0;

    (*iB)[n] = 0;
    (*jB)[n] = 0;
    (*vB)[n] = 2. / 3.;
    n++;
    (*iB)[n] = N;
    (*jB)[n] = N;
    (*vB)[n] = 2. / 3.;
    n++;
    (*iB)[n] = 0;
    (*jB)[n] = N;
    (*vB)[n] = 1. / 3.;
    n++;
    (*iB)[n] = N;
    (*jB)[n] = 0;
    (*vB)[n] = 1. / 3.;
    n++;
    (*iB)[n] = 0;
    (*jB)[n] = 1;
    (*vB)[n] = -1. / sqrt(6.);
    n++;
    (*iB)[n] = 1;
    (*jB)[n] = 0;
    (*vB)[n] = -1. / sqrt(6.);
    n++;
    (*iB)[n] = 1;
    (*jB)[n] = N;
    (*vB)[n] = -1. / sqrt(6.);
    n++;
    (*iB)[n] = N;
    (*jB)[n] = 1;
    (*vB)[n] = -1. / sqrt(6.);
    n++;
    (*iB)[n] = 0;
    (*jB)[n] = 2;
    (*vB)[n] = 1. / sqrt(90.);
    n++;
    (*iB)[n] = 2;
    (*jB)[n] = 0;
    (*vB)[n] = 1. / sqrt(90.);
    n++;
    (*iB)[n] = 2;
    (*jB)[n] = N;
    (*vB)[n] = -1. / sqrt(90.);
    n++;
    (*iB)[n] = N;
    (*jB)[n] = 2;
    (*vB)[n] = -1. / sqrt(90.);
    n++;

    for (Int k = 1; k < N; k++)
    {
        (*iB)[n] = k;
        (*jB)[n] = k;
        (*vB)[n] = 0.5 / (k - 0.5) / (k + 1.5);
        n++;
    }
    for (Int k = 1; k < N - 2; k++)
    {
        (*iB)[n] = k + 2;
        (*jB)[n] = k;
        (*vB)[n] = -0.25 / (k + 1.5) / sqrt((k + 0.5) * (k + 2.5));
        n++;
        (*iB)[n] = k;
        (*jB)[n] = k + 2;
        (*vB)[n] = -0.25 / (k + 1.5) / sqrt((k + 0.5) * (k + 2.5));
        n++;
    }
}

void getFhat(Int N, Real **vF, Int **iF)
{
    *vF = malloc((N + 3) * sizeof(Real));
    *iF = malloc((N + 3) * sizeof(Int));

    Int n = 0;

    (*iF)[n] = 0;
    (*vF)[n] = 1;
    n++;
    (*iF)[n] = N;
    (*vF)[n] = 1;
    n++;
    (*iF)[n] = 1;
    (*vF)[n] = -2. / sqrt(6.);
    n++;

    return;
}

// H_hat0 = coo_matrix(([1],([0],[0])), shape=(N+1,N+1))
// H_hat1 = coo_matrix(([1],([N],[N])), shape=(N+1,N+1))

// G_hat0 = coo_matrix(([1],([0],[0])), shape=(N+1,1))
// G_hat1 = coo_matrix(([1],([N],[0])), shape=(N+1,1))


// KA = kron(A_hat, B_hat) + kron(B_hat, A_hat)
void getKAhat(Int N, Real **vKA, Int **iKA, Int **jKA)
{
    Int *iA, *jA, *iB, *jB;
    Real *vA, *vB;
    getAhat(N, &vA, &iA, &jA);
    getBhat(N, &vB, &iB, &jB);
    Int nnzA = N + 3;
    Int nnzB = 3 * N + 5;

    *vKA = malloc((2 * nnzA * nnzB) * sizeof(Real));
    *iKA = malloc((2 * nnzA * nnzB) * sizeof(Int));
    *jKA = malloc((2 * nnzA * nnzB) * sizeof(Int));

    Int n = 0;

    // sparse kron
    for (Int k1 = 0; k1 < nnzA; k1++)
    {
        for (Int k2 = 0; k2 < nnzB; k2++)
        {
            (*iKA)[n] = iA[k1] * (N + 1) + iB[k2];
            (*jKA)[n] = jA[k1] * (N + 1) + jB[k2];
            (*vKA)[n] = vA[k1] * vB[k2];
            n++;
        }
    }

    // sparse kron
    for (Int k1 = 0; k1 < nnzA; k1++)
    {
        for (Int k2 = 0; k2 < nnzB; k2++)
        {
            (*iKA)[n] = iA[k1] + iB[k2] * (N + 1);
            (*jKA)[n] = jA[k1] + jB[k2] * (N + 1);
            (*vKA)[n] = vA[k1] * vB[k2];
            n++;
        }
    }
}

// KB = kron(B_hat, B_hat)
void getKBhat(Int N, Real **vKB, Int **iKB, Int **jKB)
{
    Int *iB, *jB;
    Real *vB;
    getBhat(N, &vB, &iB, &jB);
    Int nnzB = 3 * N + 5;

    *vKB = malloc((nnzB * nnzB) * sizeof(Real));
    *iKB = malloc((nnzB * nnzB) * sizeof(Int));
    *jKB = malloc((nnzB * nnzB) * sizeof(Int));

    Int n = 0;

    // sparse kron
    for (Int k1 = 0; k1 < nnzB; k1++)
    {
        for (Int k2 = 0; k2 < nnzB; k2++)
        {
            (*iKB)[n] = iB[k1] * (N + 1) + iB[k2];
            (*jKB)[n] = jB[k1] * (N + 1) + jB[k2];
            (*vKB)[n] = vB[k1] * vB[k2];
            n++;
        }
    }
}

// KF = kron(F_hat, F_hat)
void getKFhat(Int N, Real **vKF, Int **iKF)
{
    Int *iF;
    Real *vF;
    getFhat(N, &vF, &iF);
    Int nnzF = 3;

    *vKF = malloc((nnzF * nnzF) * sizeof(Real));
    *iKF = malloc((nnzF * nnzF) * sizeof(Int));

    Int n = 0;

    // sparse kron
    for (Int k1 = 0; k1 < nnzF; k1++)
    {
        for (Int k2 = 0; k2 < nnzF; k2++)
        {
            (*iKF)[n] = iF[k1] * (N + 1) + iF[k2];
            (*vKF)[n] = vF[k1] * vF[k2];
            n++;
        }
    }
}

// KHB0 = kron(H_hat0, B_hat)
// KHB1 = kron(H_hat1, B_hat)
// KBH0 = kron(B_hat, H_hat0)
// KBH1 = kron(B_hat, H_hat1)

// KGF0 = kron(G_hat0, F_hat)
// KGF1 = kron(G_hat1, F_hat)
// KFG0 = kron(F_hat, G_hat0)
// KFG1 = kron(F_hat, G_hat1)
