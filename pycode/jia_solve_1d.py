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


K = 10000
p = 0.5
h = 0
beta = 1

N = 5
M = 100
hm = 1/M

np.random.seed(0)
V = np.random.rand(M)
V[V<=p] = 0
V[V>p] = 1
KV = K * V

plt.bar(x=np.arange(0,1,1/M)+1/(1*M), height=KV, width=1/M)

def mmp(m, n):
    return m*N + n

def phi(n, x):
    if n == 0:
        return (1-x)/2
    elif n == N:
        return (1+x)/2
    else:
        return (legendre(n+1)(x) - legendre(n-1)(x)) / np.sqrt(4*n+2)

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
    
A = np.zeros(shape=(M*N+1, M*N+1))
F = np.zeros(shape=(M*N+1,))

for m in range(M):
    Ae = 2/hm * A_hat + hm/2 * KV[m] * B_hat
    Fe = hm/2 * F_hat
    
    Ae = sps.coo_matrix(Ae)
    Fe = sps.coo_matrix(Fe)
    
    for k in range(Ae.nnz):
        r = mmp(m, Ae.row[k])
        c = mmp(m, Ae.col[k])
        A[r,c] += Ae.data[k]
        
    for k in range(Fe.nnz):
        i = mmp(m, Fe.row[k])
        F[i] += Fe.data[k]

F[0] += h/beta
F[-1] += h/beta

A = sps.csc_matrix(A)
U = spsl.spsolve(A, F)

x = np.linspace(0, 1, 20*M +1)[:-1]
xx = np.linspace(-1, 1, 20 +1)[:-1]
y = np.zeros_like(x)
UU = U
for m in range(M):
    DOF = UU[m*N: (m+1)*N+1]
    for ii in range(N+1):
        y[m*20: (m+1)*20] += DOF[ii] * phi(ii, xx)

plt.figure()
plt.plot(x, y)
plt.show()
    