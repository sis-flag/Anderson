# !python3
# -*- coding: utf-8 -*-
# author: flag

import numpy as np
import scipy as sp
import scipy.linalg as spl
import scipy.sparse as sps
import scipy.sparse.linalg as spsl

from scipy.special import legendre

import matplotlib.pyplot as plt

def phi(n, x):
    if n == 0:
        return (1-x)/2
    elif n == N:
        return (1+x)/2
    else:
        return (legendre(n+1)(x) - legendre(n-1)(x)) / np.sqrt(4*n+2)

N = 4
M = 5
hm = 1/M

Vmax = 1000
a = 0

np.random.seed(0)
#V = Vmax * np.ones(shape=(M,M))
#V = Vmax * np.random.rand(M,M)
V = Vmax * np.random.randint(2, size=(M,M))

fig = plt.figure()
ax = fig.gca(projection='3d')
x = np.arange(0,1,1/M)
X, Y = np.meshgrid(x, x)
ax.bar3d(X.ravel(), Y.ravel(), 0, 1/M, 1/M, V.ravel())
plt.show()

def mmp(m1, m2, n):
    n1, n2 = n//(N+1), n%(N+1)
    return (m2*N+n1)*(M*N+1) + m1*N + n2

A_hat = np.zeros(shape=(N+1, N+1))
A_hat[0,0] = A_hat[N,N] = 1/2
A_hat[N,0] = A_hat[0,N] = -1/2
for k in range(1,N):
    A_hat[k,k] = 1

B_hat = np.zeros(shape=(N+1, N+1))
B_hat[0,0] = B_hat[N,N] = 2/3
B_hat[0,1] = B_hat[1,0] = -1/np.sqrt(6)
B_hat[0,2] = B_hat[2,0] = 1/np.sqrt(90)
B_hat[0,N] = B_hat[N,0] = 1/3
B_hat[1,N] = B_hat[N,1] = -1/np.sqrt(6)
B_hat[2,N] = B_hat[N,2] = -1/np.sqrt(90)
for k in range(1,N):
    B_hat[k,k] = 2/(2*k-1)/(2*k+3)
for k in range(1,N-2):
    B_hat[k+2,k] = -1/(2*k+3)/np.sqrt(2*k+1)/np.sqrt(2*k+5)
    B_hat[k,k+2] = -1/(2*k+3)/np.sqrt(2*k+1)/np.sqrt(2*k+5)
    
A = np.zeros(shape=((M*N+1)**2, (M*N+1)**2))
B = np.zeros(shape=((M*N+1)**2, (M*N+1)**2))

for m1 in range(M):
    for m2 in range(M):
        Ae = 1/4 * (spl.kron(A_hat, B_hat) + spl.kron(B_hat, A_hat)) \
           + hm*hm/4 * V[m1,m2] * spl.kron(B_hat, B_hat)
        Be = hm*hm/4 * spl.kron(B_hat, B_hat)
        
        Ae = sps.coo_matrix(Ae)
        Be = sps.coo_matrix(Be)
        
        for k in range(Ae.nnz):
            r = mmp(m1, m2, Ae.row[k])
            c = mmp(m1, m2, Ae.col[k])
            A[r,c] += Ae.data[k]
        
        for k in range(Be.nnz):
            r = mmp(m1, m2, Be.row[k])
            c = mmp(m1, m2, Be.col[k])
            B[r,c] += Be.data[k]

# 暂时没有第三类边界条件
#H_tilde = np.zeros(shape=(M*N+1, M*N+1))
#H_tilde[0,0] = H_tilde[-1,-1] = a
#H = hm*hm/4 * (spl.kron(H_tilde, B_hat) + spl.kron(B_hat, H_tilde))
#A = A + H

A = sps.csc_matrix(A)
B = sps.csc_matrix(B)
lam, U = spsl.eigs(A=A, M=B, k=4, sigma=0)

print(lam.real)

x = np.linspace(0, 1, 20*M)
x1, x2 = np.meshgrid(x, x)
xx = np.linspace(-1, 1, 20+1)[:-1]
xx1, xx2 = np.meshgrid(xx, xx)
for ind in range(4):    
    y = np.zeros_like(x1)
    UU = U[:,ind].reshape((M*N+1, M*N+1)).real
    for m1 in range(M):
        for m2 in range(M):
            DOF = UU[m1*N: (m1+1)*N+1, m2*N: (m2+1)*N+1]
            for ii in range(N+1):
                for jj in range(N+1):
                    y[m2*20:(m2+1)*20, m1*20:(m1+1)*20] += \
                        DOF[ii,jj] * phi(ii, xx1) * phi(jj, xx2)

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(x1, x2, y, cmap='rainbow')
    plt.show()
    