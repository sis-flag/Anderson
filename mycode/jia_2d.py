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
from mpl_toolkits.mplot3d import Axes3D

K = 3000
p = 0.3
h = 0
M = 20
beta = 1

N = 6
hm = 1/M

np.random.seed(0)
V = np.random.rand(M,M)
#V[V<=p] = 0
#V[V>p] = 1
KV = K * V

#fig = plt.figure()
#ax = fig.gca(projection='3d')
#x = np.arange(0,1,1/M)
#X, Y = np.meshgrid(x, x)
#ax.bar3d(X.ravel(), Y.ravel(), 0, 1/M, 1/M, KV.ravel())
#plt.show()

x = np.arange(0,1,1/M)
X, Y = np.meshgrid(x, x)
plt.pcolor(X, Y, V, cmap='rainbow')
plt.show()

def mmp(m1, m2, n):
    n1, n2 = n//(N+1), n%(N+1)
    return (m1*N+n1)*(M*N+1) + m2*N + n2

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


A = np.zeros(shape=((M*N+1)**2, (M*N+1)**2))
B = np.zeros(shape=((M*N+1)**2, (M*N+1)**2))

for m1 in range(M):
    for m2 in range(M):
        Ae = sps.kron(A_hat, B_hat) + sps.kron(B_hat, A_hat) \
           + hm*hm/4 * KV[m1,m2] * sps.kron(B_hat, B_hat)
        Be = hm*hm/4 * sps.kron(B_hat, B_hat)
        
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

# lower boundary
m1 = 0
for m2 in range(M):
    Ae = hm/2 * h * sps.kron(H_hat0, B_hat)
    Ae = sps.coo_matrix(Ae)
    for k in range(Ae.nnz):
        r = mmp(m1, m2, Ae.row[k])
        c = mmp(m1, m2, Ae.col[k])
        A[r,c] += Ae.data[k]
        
# upper boundary
m1 = M-1
for m2 in range(M):
    Ae = hm/2 * h * sps.kron(H_hat1, B_hat)
    Ae = sps.coo_matrix(Ae)
    for k in range(Ae.nnz):
        r = mmp(m1, m2, Ae.row[k])
        c = mmp(m1, m2, Ae.col[k])
        A[r,c] += Ae.data[k]

# left boundary
m2 = 0
for m1 in range(M):
    Ae = hm/2 * h * sps.kron(B_hat, H_hat0)
    Ae = sps.coo_matrix(Ae)
    for k in range(Ae.nnz):
        r = mmp(m1, m2, Ae.row[k])
        c = mmp(m1, m2, Ae.col[k])
        A[r,c] += Ae.data[k]
        
# right boundary
m2 = M-1
for m1 in range(M):
    Ae = hm/2 * h * sps.kron(B_hat, H_hat1)    
    Ae = sps.coo_matrix(Ae)
    for k in range(Ae.nnz):
        r = mmp(m1, m2, Ae.row[k])
        c = mmp(m1, m2, Ae.col[k])
        A[r,c] += Ae.data[k]

A = sps.csc_matrix(A)
B = sps.csc_matrix(B)
lam, U_eig = spsl.eigsh(A=A, M=B, k=6, sigma=0)


A = np.zeros(shape=((M*N+1)**2, (M*N+1)**2))
F = np.zeros(shape=((M*N+1)**2, 1))

for m1 in range(M):
    for m2 in range(M):
        Ae = sps.kron(A_hat, B_hat) + sps.kron(B_hat, A_hat) \
           + hm*hm/4 * KV[m1,m2] * sps.kron(B_hat, B_hat)
        Fe = hm*hm/4 * sps.kron(F_hat, F_hat)
        
        Ae = sps.coo_matrix(Ae)
        Fe = sps.coo_matrix(Fe)
        
        for k in range(Ae.nnz):
            r = mmp(m1, m2, Ae.row[k])
            c = mmp(m1, m2, Ae.col[k])
            A[r,c] += Ae.data[k]
        
        for k in range(Fe.nnz):
            i = mmp(m1, m2, Fe.row[k])
            F[i] += Fe.data[k]
            

# lower boundary
m1 = 0
for m2 in range(M):
    Fe = hm/2 * h/beta * sps.kron(G_hat0, F_hat)
    Fe = sps.coo_matrix(Fe)
    for k in range(Fe.nnz):
        i = mmp(m1, m2, Fe.row[k])
        F[i] += Fe.data[k]
        
# upper boundary
m1 = M-1
for m2 in range(M):
    Fe = hm/2 * h/beta * sps.kron(G_hat1, F_hat)
    Fe = sps.coo_matrix(Fe)
    for k in range(Fe.nnz):
        i = mmp(m1, m2, Fe.row[k])
        F[i] += Fe.data[k]
        
# left boundary
m2 = 0
for m1 in range(M):
    Fe = hm/2 * h/beta * sps.kron(F_hat, G_hat0)
    Fe = sps.coo_matrix(Fe)
    for k in range(Fe.nnz):
        i = mmp(m1, m2, Fe.row[k])
        F[i] += Fe.data[k]
        
# right boundary
m2 = M-1
for m1 in range(M):
    Fe = hm/2 * h/beta * sps.kron(F_hat, G_hat1)
    Fe = sps.coo_matrix(Fe)
    for k in range(Fe.nnz):
        i = mmp(m1, m2, Fe.row[k])
        F[i] += Fe.data[k]


A = sps.csc_matrix(A)
U_solve = spsl.spsolve(A, F)


def normalize(y):
    if np.max(y) < -np.min(y):
        return y / np.min(y)
    else:
        return y / np.max(y)
    

print(lam)

x = np.linspace(0, 1, 7*M +1)[:-1]
x1, x2 = np.meshgrid(x, x)

xx = np.linspace(-1, 1, 7+1)[:-1]
xx1, xx2 = np.meshgrid(xx, xx)
yy1, yy2 = [], []
for ii in range(N+1):
    yy1.append( phi(ii, xx1) )
    yy2.append( phi(ii, xx2) )
    
w = np.zeros_like(x1)
UU = U_solve.reshape((M*N+1, M*N+1)).real
for m1 in range(M):
    for m2 in range(M):
        DOF = UU[m1*N: (m1+1)*N+1, m2*N: (m2+1)*N+1]
        for ii in range(N+1):
            for jj in range(N+1):
                w[m2*7:(m2+1)*7,m1*7:(m1+1)*7] += \
                    DOF[ii,jj] * yy1[ii] * yy2[jj]

u = []
for ind in range(6):    
    u_t = np.zeros_like(x1)
    UU = U_eig[:,ind].reshape((M*N+1, M*N+1))
    for m1 in range(M):
        for m2 in range(M):
            DOF = UU[m1*N: (m1+1)*N+1, m2*N: (m2+1)*N+1]
            for ii in range(N+1):
                for jj in range(N+1):
                    u_t[m2*7:(m2+1)*7,m1*7:(m1+1)*7] += \
                        DOF[ii,jj] * yy1[ii] * yy2[jj]
    u.append(normalize(u_t))
    



#fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.plot_surface(x1, x2, w, cmap='rainbow')
#plt.show()
#
#fig = plt.figure()
#ax = fig.gca()
#cax = ax.pcolor(x1, x2, w, cmap='rainbow')
#fig.colorbar(cax)
#plt.show()
#
#for ind in range(6):
#    fig = plt.figure()
#    ax = fig.gca(projection='3d')
#    ax.plot_surface(x1, x2, normalize(u[ind])/(lam[ind]), cmap='rainbow')
#    plt.show()
#    
#    fig = plt.figure()
#    ax = fig.gca()
#    cax = ax.pcolor(x1, x2, normalize(u[ind])/(lam[ind]), cmap='rainbow')
#    fig.colorbar(cax)
#    plt.show()
    
#%% plot
fig = plt.figure()
ax = fig.gca()
cax = ax.pcolor(x1, x2, w, cmap='rainbow')
fig.colorbar(cax)
plt.show()

alpha = 1 / np.max(w)

for ind in range(6):
#    fig = plt.figure()
#    ax = fig.gca()
#    cax = ax.pcolor(x1, x2, w*lam[ind], vmin=1, vmax=1.001, cmap='rainbow')
#    cax.cmap.set_under('black')
#    cax.cmap.set_over('white')
#    fig.colorbar(cax)
#    plt.show()
    
    fig = plt.figure()
    ax = fig.gca()
    cax = ax.pcolor(x1, x2, u[ind], vmin=-1, vmax=1, cmap='rainbow')
    fig.colorbar(cax)
    plt.show()
