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

V = 1
a = 0

N = 10

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


Ae = A_hat + V*B_hat + a*H_hat
Be = B_hat

Ae = sps.csc_matrix(Ae)
Be = sps.csc_matrix(Be)
lam, U = spsl.eigs(A=Ae, M=Be, k=4, sigma=0)

print(lam.real)

x = np.linspace(-1,1,100)
for ind in range(4):
    y = 0
    for ii in range(N+1):
        y += U[ii,ind] * phi(ii, x)
    plt.figure()
    plt.plot(x, y)
    plt.show()