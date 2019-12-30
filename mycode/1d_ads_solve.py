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

N = 5
M = 100
hm = 1/M

Vmax = 1000
a = 0

np.random.seed(0)
#V = Vmax * np.ones(shape=(M,))
#V = Vmax * np.random.rand(M)
V = Vmax * np.random.randint(2, size=M)

plt.bar(x=np.arange(0,1,1/M)+1/(1*M), height=V, width=1/M)

def mmp(m, n):
    return m*N + n

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

H_hat = np.zeros(shape=(N+1, N+1))
H_hat[0,0] = 1
H_hat[N,N] = 1

F_hat = np.zeros(shape=(N+1,1))
F_hat[0,0] = F_hat[N,0] = 1
F_hat[1,0] = -2/np.sqrt(6)

A = np.zeros(shape=(M*N+1, M*N+1))
F = np.zeros(shape=(M*N+1,))

for m in range(M):
    Ae = 2/hm * A_hat + hm/2 * V[m] * B_hat
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

A[0,0] += a
A[-1,-1] += a

A = sps.csc_matrix(A)
U = spsl.spsolve(A, F)

x = np.linspace(0, 1, 20*M)
xx = np.linspace(-1, 1, 20 +1)[:-1]
y = np.zeros_like(x)
UU = U.real
for m in range(M):
    DOF = UU[m*N: (m+1)*N+1]
    for ii in range(N+1):
        y[m*20: (m+1)*20] += DOF[ii] * phi(ii, xx)

plt.figure()
plt.plot(x, y)
plt.show()
