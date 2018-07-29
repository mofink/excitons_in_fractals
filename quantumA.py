#! Outputs Matrix A for Fibb 3-D Potential Well
import math
import cmath
import numpy as np
import pprint as pprint
import matplotlib.pyplot as plt
from IDOS import *
from gen_basis import *
from interaction import *
from matrix_elements import *
from vals_vects import *
from sparse_matrix import *
from scipy.sparse.linalg import eigsh



w_max = 13 # Sequence order
L = 1.0
W_A = 0.01
W_B = 0.0188
n_max = 100
W_Fib = gen_Fibonacci_seq(w_max,W_A,W_B)
M = len(W_Fib)
a = L / M
N = 1000 # number of discrete points
eta = 0.204
h = L/N #size of step
V, w = Fib_potential(W_Fib, L, a, eta, N, h) 
phi = np.array([[cmath.exp(1j*(2*pi*n/L)*k*h) for k in xrange(N)] for n in xrange(n_max)], dtype=np.complex)
vals,vects = get_eigsystem(n_max, V, phi, N, h)

V_0 = 10**(-3)
Np, Ns = 4, 11
vmat = new_find_interaction_matrix(vects,V_0,phi,w,Ns,h)


state_basis = find_states(Np,Ns) #generates state basis

sparse_mtrx = createSparseMatrix(state_basis,vals,vmat,Np)
sparse_mtrx = sparse_mtrx + sparse_mtrx.getH()
neigs = 10
mbvals,mbvects = eigsh(sparse_mtrx,k=neigs,which='SM')

EnS = 1.05**2*100/(2*4.5*1.6*9.1)
aex = 0.8

x = xrange(n_max)
mbx = xrange(neigs)
print mbvals * EnS / (aex*len(W_Fib))**2
f, (ax1, ax2) = plt.subplots(1, 2)
ax1.scatter(x,vals * EnS / (aex*M)**2)
ax1.set_xlim(0, 100)
ax2.scatter(mbx, mbvals * EnS / (aex*M)**2)
plt.show()

#print getMatrixElement(i,j,2,state_basis,vals,vmat,N) #i,j are indices, k is max difference to consider,N is number of particles



