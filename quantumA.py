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
W_A = 0.01
W_B = 0.0188
m_max = n_max = 10
vals,vects,A,phi,w,h = calc_A_matrix(w_max,W_A,W_B,m_max,n_max)
V_0 = 1
vmat = new_find_interaction_matrix(vects,V_0,phi,w,4,h)

N, Ns = 2, 4

state_basis = find_states(N,Ns) #generates state basis

sparse_mtrx = createSparseMatrix(state_basis,vals,vmat,N)
sparse_mtrx = sparse_mtrx + sparse_mtrx.getH()

vals,vects = eigsh(sparse_mtrx,k=10,sigma=0,which='LM')

print vals



#print getMatrixElement(i,j,2,state_basis,vals,vmat,N) #i,j are indices, k is max difference to consider,N is number of particles



