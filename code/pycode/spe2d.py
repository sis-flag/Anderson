# !python3
# -*- coding: utf-8 -*-
# author: flag

import numpy as np
import scipy as sp
import scipy.io
import scipy.linalg as spl
import scipy.sparse as sps
import scipy.sparse.linalg as spsl

from scipy.special import legendre

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

#%% parameters
N = 10 # degree of polynomial
Ns = 40 # number of sample points in each region

#%% basis function
def phi(n, x):
    if n == 0:
        return (1-x)/2
    elif n == N:
        return (1+x)/2
    else:
        return (legendre(n+1)(x) - legendre(n-1)(x)) / np.sqrt(4*n+2)

#%% matrix on reference element
rA = [0, N, 0, N]
cA = [0, N, N, 0]
vA = [1/2, 1/2, -1/2, -1/2]
for k in range(1,N):
    rA.append(k)
    cA.append(k)
    vA.append(1)
A_hat = sps.coo_matrix((vA,(rA,cA)), shape=(N+1,N+1))

rB = [0, N, 0, 1, 0, 2, 0, N, 1, N, 2, N]
cB = [0, N, 1, 0, 2, 0, N, 0, N, 1, N, 2]
vB = [2/3, 2/3, -1/np.sqrt(6), -1/np.sqrt(6), 1/np.sqrt(90), 1/np.sqrt(90), \
      1/3, 1/3, -1/np.sqrt(6), -1/np.sqrt(6), -1/np.sqrt(90), -1/np.sqrt(90) ]
for k in range(1,N):
    rB.append(k)
    cB.append(k)
    vB.append(2/(2*k-1)/(2*k+3))
for k in range(1,N-2):
    rB.append(k)
    cB.append(k+2)
    vB.append(-1/(2*k+3)/np.sqrt(2*k+1)/np.sqrt(2*k+5))
    rB.append(k+2)
    cB.append(k)
    vB.append(-1/(2*k+3)/np.sqrt(2*k+1)/np.sqrt(2*k+5))
B_hat = sps.coo_matrix((vB,(rB,cB)), shape=(N+1,N+1))

H_hat0 = sps.coo_matrix(([1],([0],[0])), shape=(N+1,N+1))
H_hat1 = sps.coo_matrix(([1],([N],[N])), shape=(N+1,N+1))

F_hat = sps.coo_matrix(([1,1,-2/np.sqrt(6)],([0,N,1],[0,0,0])), shape=(N+1,1))

G_hat0 = sps.coo_matrix(([1],([0],[0])), shape=(N+1,1))
G_hat1 = sps.coo_matrix(([1],([N],[0])), shape=(N+1,1))

#%% tensor products on 2d
KA = sps.kron(A_hat, B_hat) + sps.kron(B_hat, A_hat)
KB = sps.kron(B_hat, B_hat)
KF = sps.kron(F_hat, F_hat)

KHB0 = sps.kron(H_hat0, B_hat)
KHB1 = sps.kron(H_hat1, B_hat)
KBH0 = sps.kron(B_hat, H_hat0)
KBH1 = sps.kron(B_hat, H_hat1)

KGF0 = sps.kron(G_hat0, F_hat)
KGF1 = sps.kron(G_hat1, F_hat)
KFG0 = sps.kron(F_hat, G_hat0)
KFG1 = sps.kron(F_hat, G_hat1)

#%% samples on reference element
xx = np.linspace(-1, 1, Ns +1)
xx2, xx1 = np.meshgrid(xx, xx) #caution this
yy1, yy2 = [], []
for n in range(N+1):
    yy1.append( phi(n, xx1) )
    yy2.append( phi(n, xx2) )

#%% normalize
def normalize(y):
    if np.max(y) < -np.min(y):
        return y / np.min(y)
    else:
        return y / np.max(y)
    
    
#%% 2d Robin boundary eigen problem
def eigN(V, h, num=5):
    """
    solve the eigen problem
    - laplace u(x) + V(x) u(x) = lam u(x) for x in [0,1] x [0,1]
    u'(x) + h u(x) = 0 for x on boundary
    V(x) is piecewise constant
    """
    M1, M2 = V.shape
    hm1, hm2 = 1/M1, 1/M2
#    assert M1 == M2
    
    # map from elemental to global
    def mmp(m1, m2, n):
        n1, n2 = n%(N+1), n//(N+1)
        u1, u2 = m1*N+n1, m2*N+n2
        return u2*(M1*N+1) + u1
    
    # assemble global matrix
    rA, cA, vA = [], [], []
    rB, cB, vB = [], [], []
    for m1 in range(M1):
        for m2 in range(M2):
            Ae = KA + hm1*hm2/4 * V[m1,m2] * KB
            Be = hm1*hm2/4 * KB
            
            Ae = sps.coo_matrix(Ae)
            Be = sps.coo_matrix(Be)
            
            for k in range(Ae.nnz):
                rA.append( mmp(m1, m2, Ae.row[k]) )
                cA.append( mmp(m1, m2, Ae.col[k]) )
                vA.append( Ae.data[k] )
            
            for k in range(Be.nnz):
                rB.append( mmp(m1, m2, Be.row[k]) )
                cB.append( mmp(m1, m2, Be.col[k]) )
                vB.append( Be.data[k] )

    # left boundary
    m2 = 0
    for m1 in range(M2):
        Ae = sps.coo_matrix(hm2/2 * h * KHB0)
        for k in range(Ae.nnz):
            rA.append( mmp(m1, m2, Ae.row[k]) )
            cA.append( mmp(m1, m2, Ae.col[k]) )
            vA.append( Ae.data[k] )
            
    # right boundary
    m2 = M1-1
    for m1 in range(M2):
        Ae = sps.coo_matrix(hm2/2 * h * KHB1)
        for k in range(Ae.nnz):
            rA.append( mmp(m1, m2, Ae.row[k]) )
            cA.append( mmp(m1, m2, Ae.col[k]) )
            vA.append( Ae.data[k] )
    
    # bottom boundary
    m1 = 0
    for m2 in range(M1):
        Ae = sps.coo_matrix(hm1/2 * h * KBH0)
        for k in range(Ae.nnz):
            rA.append( mmp(m1, m2, Ae.row[k]) )
            cA.append( mmp(m1, m2, Ae.col[k]) )
            vA.append( Ae.data[k] )
            
    # top boundary
    m1 = M2-1
    for m2 in range(M1):
        Ae = sps.coo_matrix(hm1/2 * h * KBH1)
        for k in range(Ae.nnz):
            rA.append( mmp(m1, m2, Ae.row[k]) )
            cA.append( mmp(m1, m2, Ae.col[k]) )
            vA.append( Ae.data[k] )
    
    # solve
    shape = ((M1*N+1)*(M2*N+1), (M1*N+1)*(M2*N+1))
    A = sps.coo_matrix((vA,(rA,cA)), shape=shape).tocsr()
    B = sps.coo_matrix((vB,(rB,cB)), shape=shape).tocsr()
    lam, U = spsl.eigs(A=A, M=B, k=num, sigma=0)
    lam, U = lam.real, U.real

    u = []
    for ind in range(num):    
        # reshape to 2d indicate
        UU = U[:,ind].reshape((M1*N+1, M2*N+1), order="F")
        
        u_t = np.zeros(shape=(Ns*M1+1, Ns*M2+1))
        for m1 in range(M1):
            for m2 in range(M2):
                Uloc = UU[m1*N: (m1+1)*N+1, m2*N: (m2+1)*N+1]
                u_t[m1*Ns: (m1+1)*Ns+1, m2*Ns: (m2+1)*Ns+1] = 0
                for n1 in range(N+1):
                    for n2 in range(N+1):
                        u_t[m1*Ns: (m1+1)*Ns+1, m2*Ns: (m2+1)*Ns+1] += \
                            Uloc[n1,n2] * yy1[n1] * yy2[n2]
                            
        u.append(normalize(u_t))
        
    return lam, u

#%% 2d Nunmman boundary source problem
def solveN(V, g):
    """
    solve the source problem
    - laplace u(x) + V(x) u(x) = 1 for x in [0,1] x [0,1]
    u'(x) = g for x on boundary
    V(x) is piecewise constant
    """
    M1, M2 = V.shape
    hm1, hm2 = 1/M1, 1/M2
#    assert M1 == M2
    
    # map from elemental to global
    def mmp(m1, m2, n):
        n1, n2 = n%(N+1), n//(N+1)
        u1, u2 = m1*N+n1, m2*N+n2
        return u2*(M1*N+1) + u1
    
    # assemble global matrix
    rA, cA, vA = [], [], []
    F = np.zeros(shape=((M1*N+1)*(M2*N+1), 1))
    for m1 in range(M1):
        for m2 in range(M2):
            Ae = KA + hm1*hm2/4 * V[m1,m2] * KB
            Fe = hm1*hm2/4 * KF
            
            Ae = sps.coo_matrix(Ae)
            Fe = sps.coo_matrix(Fe)
            
            for k in range(Ae.nnz):
                rA.append( mmp(m1, m2, Ae.row[k]) )
                cA.append( mmp(m1, m2, Ae.col[k]) )
                vA.append( Ae.data[k] )
            
            for k in range(Fe.nnz):
                i = mmp(m1, m2, Fe.row[k])
                F[i] += Fe.data[k]
                
    # left boundary
    m2 = 0
    for m1 in range(M2):
        Fe = sps.coo_matrix(hm2/2 * g * KGF0)
        for k in range(Fe.nnz):
            i = mmp(m1, m2, Fe.row[k])
            F[i] += Fe.data[k]
            
    # right boundary
    m2 = M1-1
    for m1 in range(M2):
        Fe = sps.coo_matrix(hm2/2 * g * KGF1)
        for k in range(Fe.nnz):
            i = mmp(m1, m2, Fe.row[k])
            F[i] += Fe.data[k]
            
    # bottom boundary
    m1 = 0
    for m2 in range(M1):
        Fe = sps.coo_matrix(hm1/2 * g * KFG0)
        for k in range(Fe.nnz):
            i = mmp(m1, m2, Fe.row[k])
            F[i] += Fe.data[k]
            
    # top boundary
    m1 = M2-1
    for m2 in range(M1):
        Fe = sps.coo_matrix(hm1/2 * g * KFG1)
        for k in range(Fe.nnz):
            i = mmp(m1, m2, Fe.row[k])
            F[i] += Fe.data[k]
                
    # solve
    shape = ((M1*N+1)*(M2*N+1), (M1*N+1)*(M2*N+1))
    A = sps.coo_matrix((vA,(rA,cA)), shape=shape).tocsr()
    U = spsl.spsolve(A, F)
    
    # reshape to 2d indicate
    UU = U.reshape((M1*N+1, M2*N+1), order="F")
    
    w = np.zeros(shape=(Ns*M1+1, Ns*M2+1))
    for m1 in range(M1):
        for m2 in range(M2):
            Uloc = UU[m1*N: (m1+1)*N+1, m2*N: (m2+1)*N+1]
            w[m1*Ns: (m1+1)*Ns+1, m2*Ns: (m2+1)*Ns+1] = 0
            for n1 in range(N+1):
                for n2 in range(N+1):
                    w[m1*Ns: (m1+1)*Ns+1, m2*Ns: (m2+1)*Ns+1] += \
                        Uloc[n1,n2] * yy1[n1] * yy2[n2]
                        
    return w


#%% 2d Mixed boundary source problem
def solveM(V, g, x0=None, y0=None, u0=None):
    """
    solve the source problem
    - laplace u(x) + V(x) u(x) = 1 for x in [0,1] x [0,1]
    u'(x) = g for x on boundary
    enforce u(x0, y0) = u0
    V(x) is piecewise constant
    """
    M1, M2 = V.shape
    hm1, hm2 = 1/M1, 1/M2
#    assert M1 == M2
    
    # map from elemental to global
    def mmp(m1, m2, n):
        n1, n2 = n%(N+1), n//(N+1)
        u1, u2 = m1*N+n1, m2*N+n2
        return u2*(M1*N+1) + u1
    
    # assemble global matrix
    rA, cA, vA = [], [], []
    F = np.zeros(shape=((M1*N+1)*(M2*N+1), 1))
    for m1 in range(M1):
        for m2 in range(M2):
            Ae = KA + hm1*hm2/4 * V[m1,m2] * KB
            Fe = hm1*hm2/4 * KF
            
            Ae = sps.coo_matrix(Ae)
            Fe = sps.coo_matrix(Fe)
            
            for k in range(Ae.nnz):
                rA.append( mmp(m1, m2, Ae.row[k]) )
                cA.append( mmp(m1, m2, Ae.col[k]) )
                vA.append( Ae.data[k] )
            
            for k in range(Fe.nnz):
                i = mmp(m1, m2, Fe.row[k])
                F[i] += Fe.data[k]
                
    # left boundary
    m2 = 0
    for m1 in range(M2):
        Fe = sps.coo_matrix(hm2/2 * g * KGF0)
        for k in range(Fe.nnz):
            i = mmp(m1, m2, Fe.row[k])
            F[i] += Fe.data[k]
            
    # right boundary
    m2 = M1-1
    for m1 in range(M2):
        Fe = sps.coo_matrix(hm2/2 * g * KGF1)
        for k in range(Fe.nnz):
            i = mmp(m1, m2, Fe.row[k])
            F[i] += Fe.data[k]
            
    # bottom boundary
    m1 = 0
    for m2 in range(M1):
        Fe = sps.coo_matrix(hm1/2 * g * KFG0)
        for k in range(Fe.nnz):
            i = mmp(m1, m2, Fe.row[k])
            F[i] += Fe.data[k]
            
    # top boundary
    m1 = M2-1
    for m2 in range(M1):
        Fe = sps.coo_matrix(hm1/2 * g * KFG1)
        for k in range(Fe.nnz):
            i = mmp(m1, m2, Fe.row[k])
            F[i] += Fe.data[k]
            
    if not x0 == None:
        m1, m2 = int(x0 // hm1), int(y0 // hm2)
        mm = mmp(m1, m2, 0)
        # delete a DOF
        i = 0
        while i < len(rA):
            if rA[i] == mm:
                del rA[i], cA[i], vA[i]
            else:
                i += 1
        
        # add enforce in
        rA.append( mm )
        cA.append( mm )
        vA.append( 1 )
        F[mm] = u0
                
    # solve
    shape = ((M1*N+1)*(M2*N+1), (M1*N+1)*(M2*N+1))
    A = sps.coo_matrix((vA,(rA,cA)), shape=shape).tocsr()
    U = spsl.spsolve(A, F)
    
    # reshape to 2d indicate
    UU = U.reshape((M1*N+1, M2*N+1), order="F")
    
    w = np.zeros(shape=(Ns*M1+1, Ns*M2+1))
    for m1 in range(M1):
        for m2 in range(M2):
            Uloc = UU[m1*N: (m1+1)*N+1, m2*N: (m2+1)*N+1]
            w[m1*Ns: (m1+1)*Ns+1, m2*Ns: (m2+1)*Ns+1] = 0
            for n1 in range(N+1):
                for n2 in range(N+1):
                    w[m1*Ns: (m1+1)*Ns+1, m2*Ns: (m2+1)*Ns+1] += \
                        Uloc[n1,n2] * yy1[n1] * yy2[n2]
                        
    return w


#%% 2d Dirichlet boundary eigen problem
def eigD(V, num=5):
    """
    solve the eigen problem
    - laplace u(x) + V(x) u(x) = lam u(x) for x in [0,1] x [0,1]
    u(x) = 0 for x on boundary
    V(x) is piecewise constant
    """    
    M1, M2 = V.shape
    hm1, hm2 = 1/M1, 1/M2
#    assert M1 == M2
    
    # map from elemental to global
    def mmp(m1, m2, n):
        n1, n2 = n%(N+1), n//(N+1)
        u1, u2 = m1*N+n1, m2*N+n2
        if (u1 in [0, M1*N]) or (u2 in [0, M2*N]):
            return -1
        else:
            return (u2-1)*(M1*N-1) + (u1-1)
    
    # assemble global matrix
    rA, cA, vA = [], [], []
    rB, cB, vB = [], [], []
    for m1 in range(M1):
        for m2 in range(M2):
            Ae = KA + hm1*hm2/4 * V[m1,m2] * KB
            Be = hm1*hm2/4 * KB
            
            Ae = sps.coo_matrix(Ae)
            Be = sps.coo_matrix(Be)
            
            for k in range(Ae.nnz):
                if mmp(m1, m2, Ae.row[k]) == -1 or mmp(m1, m2, Ae.col[k]) == -1:
                    continue
                rA.append( mmp(m1, m2, Ae.row[k]) )
                cA.append( mmp(m1, m2, Ae.col[k]) )
                vA.append( Ae.data[k] )
            
            for k in range(Be.nnz):
                if mmp(m1, m2, Be.row[k]) == -1 or mmp(m1, m2, Be.col[k]) == -1:
                    continue
                rB.append( mmp(m1, m2, Be.row[k]) )
                cB.append( mmp(m1, m2, Be.col[k]) )
                vB.append( Be.data[k] )
    
    # solve
    shape = ((M1*N-1)*(M2*N-1), (M1*N-1)*(M2*N-1))
    A = sps.coo_matrix((vA,(rA,cA)), shape=shape).tocsr()
    B = sps.coo_matrix((vB,(rB,cB)), shape=shape).tocsr()
    lam, U = spsl.eigs(A=A, M=B, k=num, sigma=0)
    lam, U = lam.real, U.real
    
    u = []
    for ind in range(num):    
        
        # enforce boundary
        UU = U[:,ind].reshape((M1*N-1, M2*N-1), order="F")
        zero = np.zeros((M1*N-1, 1))
        UU = np.concatenate((zero, UU, zero), axis=1)
        zero = np.zeros((1, M2*N+1))
        UU = np.concatenate((zero, UU, zero), axis=0)
        
        u_t = np.zeros(shape=(Ns*M1+1, Ns*M2+1))
        for m1 in range(M1):
            for m2 in range(M2):
                Uloc = UU[m1*N: (m1+1)*N+1, m2*N: (m2+1)*N+1]
                u_t[m1*Ns: (m1+1)*Ns+1, m2*Ns: (m2+1)*Ns+1] = 0
                for n1 in range(N+1):
                    for n2 in range(N+1):
                        u_t[m1*Ns: (m1+1)*Ns+1, m2*Ns: (m2+1)*Ns+1] += \
                            Uloc[n1,n2] * yy1[n1] * yy2[n2]
                            
        u.append(normalize(u_t))
        
    return lam, u

#%% 2d Dirichlet boundary source problem
def solveD(V):
    """
    solve the source problem
    - laplace u(x) + V(x) u(x) = 1 for x in [0,1] x [0,1]
    u(x) = 0 for x on boundary
    V(x) is piecewise constant
    """
    M1, M2 = V.shape
    hm1, hm2 = 1/M1, 1/M2
#    assert M1 == M2
    
    # map from elemental to global
    def mmp(m1, m2, n):
        n1, n2 = n%(N+1), n//(N+1)
        u1, u2 = m1*N+n1, m2*N+n2
        if (u1 in [0, M1*N]) or (u2 in [0, M2*N]):
            return -1
        else:
            return (u2-1)*(M1*N-1) + (u1-1)
    
    # assemble global matrix
    rA, cA, vA = [], [], []
    F = np.zeros(shape=((M1*N-1)*(M2*N-1), 1))
    for m1 in range(M1):
        for m2 in range(M2):
            Ae = KA + hm1*hm2/4 * V[m1,m2] * KB
            Fe = hm1*hm2/4 * KF
            
            Ae = sps.coo_matrix(Ae)
            Fe = sps.coo_matrix(Fe)
            
            
            for k in range(Ae.nnz):
                if mmp(m1, m2, Ae.row[k]) == -1 or mmp(m1, m2, Ae.col[k]) == -1:
                    continue
                rA.append( mmp(m1, m2, Ae.row[k]) )
                cA.append( mmp(m1, m2, Ae.col[k]) )
                vA.append( Ae.data[k] )
            
            for k in range(Fe.nnz):
                if mmp(m1, m2, Fe.row[k]) == -1:
                    continue
                i = mmp(m1, m2, Fe.row[k])
                F[i] += Fe.data[k]
                
    # solve
    shape = ((M1*N-1)*(M2*N-1), (M1*N-1)*(M2*N-1))
    A = sps.coo_matrix((vA,(rA,cA)), shape=shape).tocsr()
    U = spsl.spsolve(A, F)
    
    # enforce boundary
    UU = U.reshape((M1*N-1, M2*N-1), order="F")
    zero = np.zeros((M1*N-1, 1))
    UU = np.concatenate((zero, UU, zero), axis=1)
    zero = np.zeros((1, M2*N+1))
    UU = np.concatenate((zero, UU, zero), axis=0)
    
    w = np.zeros(shape=(Ns*M1+1, Ns*M2+1))
    for m1 in range(M1):
        for m2 in range(M2):
            Uloc = UU[m1*N: (m1+1)*N+1, m2*N: (m2+1)*N+1]
            w[m1*Ns: (m1+1)*Ns+1, m2*Ns: (m2+1)*Ns+1] = 0
            for n1 in range(N+1):
                for n2 in range(N+1):
                    w[m1*Ns: (m1+1)*Ns+1, m2*Ns: (m2+1)*Ns+1] += \
                        Uloc[n1,n2] * yy1[n1] * yy2[n2]
                        
    return w


#%% test sample
if __name__ == "__main__":
    
    #%% parameters (p=-1 for uniform distribution)
    K = 3000
    p = 0.4
    M = 5
    h = 0.1
    beta = 1
    
    #%% solve
    np.random.seed(0)
    V = np.random.rand(7, 7)
    if p > 0:
        V[V<=p] = 0
        V[V>p] = 1
    KV = K * V
    
    KV = np.ones((7,7))
    
    lam, u = eigN(KV, h)
    wN = solveN(KV, h/beta)
    wM = solveM(KV, h/beta, 0.71, 0.71, 1)
    
#    lam, u = eigD(KV)
#    w = solveD(KV)
    
#    scipy.io.savemat('2.mat',{"K":K, "p":p, "V":V, "h":h, "beta":beta, "lam":lam, "u":u, "w":w})
    
    #%% plot
    plt.rcParams['figure.figsize'] = (4.0, 3.0)
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['figure.dpi'] = 100
    
    x1 = np.linspace(0, 1, KV.shape[0] +1)
    x2 = np.linspace(0, 1, KV.shape[1] +1)
    x2, x1 = np.meshgrid(x2, x1)
    
    plt.figure()
    plt.pcolor(x1, x2, KV, vmin=0, vmax=K, cmap='rainbow')
    plt.colorbar()
    plt.title("potential")
    plt.show()

    x1 = np.linspace(0, 1, wN.shape[0] +1)
    x2 = np.linspace(0, 1, wN.shape[1] +1)
    x2, x1 = np.meshgrid(x2, x1)
    
    plt.figure()
    plt.pcolor(x1, x2, wN, vmin=0, cmap='rainbow')
    plt.colorbar()
    plt.title("landscape1")
    plt.show()
    
    plt.figure()
    plt.pcolor(x1, x2, wM, vmin=0, cmap='rainbow')
    plt.colorbar()
    plt.title("landscape2")
    plt.show()
