# !python3
# -*- coding: utf-8 -*-
# author: flag

"""
solve eigen value problem
- laplace u(x) + V u(x) = lam u(x) for x in [-1,1] x [-1,1]
u'(x) + h0 u(x) = 0 for x on boundary
"""

import numpy as np
import scipy as sp
import scipy.linalg as spl
import scipy.sparse as sps
import scipy.sparse.linalg as spsl

from scipy.special import legendre

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

V = 0
h0 = 1
N = 10

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


Ae = sps.kron(A_hat, B_hat) + sps.kron(B_hat, A_hat) \
   + V * sps.kron(B_hat, B_hat) \
   + h0 * (sps.kron(H_hat0, B_hat) + sps.kron(B_hat, H_hat0)) \
   + h0 * (sps.kron(H_hat1, B_hat) + sps.kron(B_hat, H_hat1))
Be = sps.kron(B_hat, B_hat)

Ae = sps.csc_matrix(Ae)
Be = sps.csc_matrix(Be)
lam, U = spsl.eigsh(A=Ae, M=Be, k=6, sigma=0)

print(lam)

x = np.linspace(-1,1,100)
x1, x2 = np.meshgrid(x, x)
for ind in range(6):
    y = 0
    DOF = U[:,ind].reshape((N+1, N+1))
    for ii in range(N+1):
        for jj in range(N+1):
            y += DOF[ii,jj] * phi(ii, x1) * phi(jj, x2)
            
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(x1, x2, y)
    plt.show()
    