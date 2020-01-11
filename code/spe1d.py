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


#%% parameters
N = 6 # degree of polynomial
Ns = 10 # number of sample points in each region

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

#%% samples on reference element
xx = np.linspace(-1, 1, Ns +1)
yy = []
for i in range(N+1):
    yy.append(phi(i, xx))
    
#%% normalize
def normalize(y):
    if np.max(y) < -np.min(y):
        return y / np.min(y)
    else:
        return y / np.max(y)

#%% 1d Robin boundary eigen problem
def eigN(V, h, num=5):
    """
    solve the eigen problem
    - u''(x) + V(x) u(x) = lam u(x) for x in [0,1]
    u'(x) + h u(x) = 0 for x=0 or x=1
    V(x) is piecewise constant
    """
    M = V.shape[0]
    hm = 1 / M
    
    # map from elemental to global
    def mmp(m, n):
        return m*N + n
    
    # assemble global matrix
    rA, cA, vA = [], [], []
    rB, cB, vB = [], [], []
    for m in range(M):
        Ae = 2/hm * A_hat + hm/2 * V[m] * B_hat
        Be = hm/2 * B_hat
        
        Ae = sps.coo_matrix(Ae)
        Be = sps.coo_matrix(Be)
        
        for k in range(Ae.nnz):
            rA.append( mmp(m, Ae.row[k]) )
            cA.append( mmp(m, Ae.col[k]) )
            vA.append( Ae.data[k] )
        
        for k in range(Be.nnz):
            rB.append( mmp(m, Be.row[k]) )
            cB.append( mmp(m, Be.col[k]) )
            vB.append( Be.data[k] )
    
    # boundary
    rA.extend([0, M*N])
    cA.extend([0, M*N])
    vA.extend([h, h])
    
    # solve
    A = sps.coo_matrix((vA,(rA,cA)), shape=(M*N+1, M*N+1)).tocsr()
    B = sps.coo_matrix((vB,(rB,cB)), shape=(M*N+1, M*N+1)).tocsr()
    lam, U = spsl.eigsh(A=A, M=B, k=num, sigma=0)
    
    u = []
    for ind in range(num):
        u_t = np.zeros(shape=(Ns*M+1))
        UU = U[:,ind]
        for m in range(M):
            Uloc = UU[m*N: (m+1)*N+1]
            u_t[m*Ns: (m+1)*Ns+1] = 0
            for n in range(N+1):
                u_t[m*Ns: (m+1)*Ns+1] += Uloc[n] * yy[n]
                
        u.append(normalize(u_t))
        
    return lam, u


#%% 1d Neumman boundary source problem
def solveN(V, g):
    """
    solve the source problem
    - u''(x) + V(x) u(x) = 1 for x in [0,1]
    u'(x) = g for x=0 or x=1
    V(x) is piecewise constant
    """
    M = V.shape[0]
    hm = 1 / M
    
    # map from elemental to global
    def mmp(m, n):
        return m*N + n
    
    # assemble global matrix
    rA, cA, vA = [], [], []
    F = np.zeros(shape=(M*N+1,))
    for m in range(M):
        Ae = 2/hm * A_hat + hm/2 * V[m] * B_hat
        Fe = hm/2 * F_hat
        
        Ae = sps.coo_matrix(Ae)
        Fe = sps.coo_matrix(Fe)
        
        for k in range(Ae.nnz):
            rA.append( mmp(m, Ae.row[k]) )
            cA.append( mmp(m, Ae.col[k]) )
            vA.append( Ae.data[k] )
        
        for k in range(Fe.nnz):
            i = mmp(m, Fe.row[k])
            F[i] += Fe.data[k]
    
    # boundary
    F[0] += g
    F[-1] += g
    
    # solve
    A = sps.coo_matrix((vA,(rA,cA)), shape=(M*N+1, M*N+1)).tocsr()
    U = spsl.spsolve(A, F)
    
    w = np.zeros(shape=(Ns*M+1))
    for m in range(M):
        Uloc = U[m*N: (m+1)*N+1]
        w[m*Ns: (m+1)*Ns+1] = 0
        for n in range(N+1):
            w[m*Ns: (m+1)*Ns+1] += Uloc[n] * yy[n]
        
    return w


#%% 1d Dirichlet boundary eigen problem
def eigD(V, num=5):
    """
    solve the eigen problem
    - u''(x) + V(x) u(x) = lam u(x) for x in [0,1]
    u(x) = 0 for x=0 or x=1
    V(x) is piecewise constant
    """
    M = V.shape[0]
    hm = 1 / M
    
    # map from elemental to global
    def mmp(m, n):
        if m*N+n == 0 or m*N+n == M*N:
            return -1
        else:
            return m*N + n - 1
    
    # assemble global matrix
    rA, cA, vA = [], [], []
    rB, cB, vB = [], [], []
    for m in range(M):
        Ae = 2/hm * A_hat + hm/2 * V[m] * B_hat
        Be = hm/2 * B_hat
        
        Ae = sps.coo_matrix(Ae)
        Be = sps.coo_matrix(Be)
        
        for k in range(Ae.nnz):
            if mmp(m, Ae.row[k]) == -1 or mmp(m, Ae.col[k]) == -1:
                   continue
            rA.append( mmp(m, Ae.row[k]) )
            cA.append( mmp(m, Ae.col[k]) )
            vA.append( Ae.data[k] )
        
        for k in range(Be.nnz):
            if mmp(m, Be.row[k]) == -1 or mmp(m, Be.col[k]) == -1:
                   continue
            rB.append( mmp(m, Be.row[k]) )
            cB.append( mmp(m, Be.col[k]) )
            vB.append( Be.data[k] )
    
    # solve
    A = sps.coo_matrix((vA,(rA,cA)), shape=(M*N-1, M*N-1)).tocsr()
    B = sps.coo_matrix((vB,(rB,cB)), shape=(M*N-1, M*N-1)).tocsr()
    lam, U = spsl.eigsh(A=A, M=B, k=num, sigma=0)
    
    u = []
    for ind in range(num):
        u_t = np.zeros(shape=(Ns*M+1))
        UU = np.concatenate(([0], U[:,ind], [0]), axis=0)
        for m in range(M):
            Uloc = UU[m*N: (m+1)*N+1]
            u_t[m*Ns: (m+1)*Ns+1] = 0
            for n in range(N+1):
                u_t[m*Ns: (m+1)*Ns+1] += Uloc[n] * yy[n]
                
        u.append(normalize(u_t))
        
    return lam, u


#%% 1d Dirichlet boundary source problem
def solveD(V):
    """
    solve the source problem
    - u''(x) + V(x) u(x) = 1 for x in [0,1]
    u(x) = 0 for x=0 or x=1
    V(x) is piecewise constant
    """
    M = V.shape[0]
    hm = 1 / M
    
    # map from elemental to global
    def mmp(m, n):
        if m*N+n == 0 or m*N+n == M*N:
            return -1
        else:
            return m*N + n - 1
    
    # assemble global matrix
    rA, cA, vA = [], [], []
    F = np.zeros(shape=(M*N-1,))
    for m in range(M):
        Ae = 2/hm * A_hat + hm/2 * V[m] * B_hat
        Fe = hm/2 * F_hat
        
        Ae = sps.coo_matrix(Ae)
        Fe = sps.coo_matrix(Fe)
        
        for k in range(Ae.nnz):
            if mmp(m, Ae.row[k]) == -1 or mmp(m, Ae.col[k]) == -1:
                   continue
            rA.append( mmp(m, Ae.row[k]) )
            cA.append( mmp(m, Ae.col[k]) )
            vA.append( Ae.data[k] )
        
        for k in range(Fe.nnz):
            if mmp(m, Fe.row[k]) == -1:
                continue
            i = mmp(m, Fe.row[k])
            F[i] += Fe.data[k]
    
    # solve
    A = sps.coo_matrix((vA,(rA,cA)), shape=(M*N-1, M*N-1)).tocsr()
    U = spsl.spsolve(A, F)
    
    w = np.zeros(shape=(Ns*M+1))
    UU = np.concatenate(([0], U, [0]), axis=0)
    for m in range(M):
        Uloc = UU[m*N: (m+1)*N+1]
        w[m*Ns: (m+1)*Ns+1] = 0
        for n in range(N+1):
            w[m*Ns: (m+1)*Ns+1] += Uloc[n] * yy[n]
        
    return w


#%% test sample
if __name__ == "__main__":
    
    #%% parameters (p=-1 for uniform distribution)
    K = 3000
    p = 0.3
    M = 200
    h = 0
    beta = 0.001
    
    #%% solve
    np.random.seed(0)
    V = np.random.rand(M)
    if p > 0:
        V[V<=p] = 0
        V[V>p] = 1
    KV = K * V
    
    lamN, uN = eigN(KV, h)
    wN = solveN(KV, h/beta)
    
    lamD, uD = eigD(KV)
    wD = solveD(KV)
    
    #%% plot
    plt.rcParams['figure.figsize'] = (4.0, 3.0)
    plt.rcParams['savefig.dpi'] = 300
    plt.rcParams['figure.dpi'] = 100
    
    plt.figure()
    plt.bar(x=np.arange(0,1,1/M)+1/(2*M), height=KV, width=1/M)
    plt.show()
    
    x = np.linspace(0, 1, len(wD))
    
    print("Dirichlet eigen value: ", lamD)
    
    alpha = 1 / np.max(wD)
    
    plt.figure()
    plt.plot(x, wD, 'k-', label='w(x)')
    for ind in range(4):
        plt.plot(x, uD[ind]/lamD[ind], label='u%d' % (ind+1) )
    plt.title("Dirichlet eigen mode")
    plt.legend()
    plt.show()
    
    print("Neumman eigen value: ", lamN)
    
    alpha = 1 / np.max(wN)
    
    plt.figure()
    plt.plot(x, wN, 'k-', label='w(x)')
    for ind in range(4):
        plt.plot(x, uN[ind]/(lamN[ind]+beta+alpha), label='u%d' % (ind+1) )
    plt.title("Neumman eigen mode")
    plt.legend()
    plt.show()
