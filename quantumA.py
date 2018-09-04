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
from occupation_numbers import *




w_max = 13 # Sequence order
L = 1.0
W_A = 0.01
W_B = 0.018
n_max = 200
W_Fib = gen_Fibonacci_seq(w_max,W_A,W_B)
M = len(W_Fib)
a = L / M
N = 1000 # number of discrete points
eta = 0.204
h = L/N #size of step

V, w = Fib_potential(W_Fib, L, a, eta, N, h) 
phi = np.array([[cmath.exp(1j*(2*pi*n/L)*k*h) for k in xrange(N)] for n in xrange(n_max)], dtype=np.complex)
vals,vects = get_eigsystem(n_max, V, phi, N, h)



##########
"""

x = xrange(n_max)


EnS = 1.05**2*100/(2*4.5*1.6*9.1)
aex = 0.8

f, (ax1, ax2) = plt.subplots(1, 2)
ax1.scatter(x,vals* EnS / (aex*M)**2)
ax1.set_xlim(-10, 30)
ax1.set_ylim(0, 5000)



W_A = 0.01
W_B = 0.011
W_Fib = gen_Fibonacci_seq(w_max,W_A,W_B)
V, w = Fib_potential(W_Fib, L, a, eta, N, h) 
phi = np.array([[cmath.exp(1j*(2*pi*n/L)*k*h) for k in xrange(N)] for n in xrange(n_max)], dtype=np.complex)
vals1,vects1 = get_eigsystem(n_max, V, phi, N, h)



aix2.scatter(x,vals1* EnS / (aex*M)**2)
ax2.set_xlim(-10, 30)
ax2.set_ylim(0, 5000)
plt.show()
"""
#######


V_0 = 100
#y_vals = np.zeros(250)

print V_0
Np, Ns = 4,6
vmat = new_find_interaction_matrix(vects,V_0,phi,w,Ns,h)

state_basis = find_states(Np,Ns) #generates state basis

sparse_mtrx = createSparseMatrix(state_basis,vals,vmat,Np)
sparse_mtrx = sparse_mtrx + sparse_mtrx.getH()
neigs = 4 #number of states
mbvals,mbvects = eigsh(sparse_mtrx,k=neigs,which='SM',ncv = max(3*neigs + 1, 30))

print mbvals

#y_vals[i] = mbvals[1] - mbvals[0]
#V_0 += 0.5


#with open("report.txt","a") as f:
#	f.write(str(y_vals))



"""
x = xrange(n_max)
mbx = xrange(neigs)
print mbvals * EnS / (aex*len(W_Fib))**2
f, (ax1, ax2) = plt.subplots(1, 2)
ax1.scatter(x,vals * EnS / (aex*M)**2)
ax1.set_xlim(0, 100)
ax2.scatter(mbx, mbvals * EnS / (aex*M)**2)
plt.show()

"""



#print "ENERGY GAP"
#print mbvals[1] - mbvals[0]
"""

print "OCCUPATION"
print occupation(mbvects,state_basis,0,neigs)
print occupation(mbvects,state_basis,1,neigs)
print occupation(mbvects,state_basis,2,neigs)
print occupation(mbvects,state_basis,3,neigs)
print occupation(mbvects,state_basis,4,neigs)
print occupation(mbvects,state_basis,5,neigs)
print occupation(mbvects,state_basis,6,neigs)
print occupation(mbvects,state_basis,7,neigs)
print occupation(mbvects,state_basis,8,neigs)
print occupation(mbvects,state_basis,9,neigs)
print occupation(mbvects,state_basis,10,neigs)
print occupation(mbvects,state_basis,11,neigs)





#print (mbvals[1] - mbvals[0])*EnS / (aex*len(W_Fib))**2

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


######

#For Pouyan

V_0 = 0
L1 = np.zeros(250)
for i in xrange(250):
	print V_0
	Np, Ns = 6, 9
	vmat = new_find_interaction_matrix(vects,V_0,phi,w,Ns,h)

	state_basis = find_states(Np,Ns) #generates state basis

	sparse_mtrx = createSparseMatrix(state_basis,vals,vmat,Np)
	sparse_mtrx = sparse_mtrx + sparse_mtrx.getH()
	neigs = 4 #number of states
	mbvals,mbvects = eigsh(sparse_mtrx,k=neigs,which='SM',ncv = max(3*neigs + 1, 30))
	V_0 += 1
	L1[i] = occupation(mbvects,state_basis,0,neigs)


with open("report.txt","a") as f:
	f.write(str(y_vals))


L1 = occupation(mbvects,state_basis,0,neigs)
L2 = occupation(mbvects,state_basis,1,neigs)
L3 = occupation(mbvects,state_basis,2,neigs)
L4 = occupation(mbvects,state_basis,3,neigs)
L5 = occupation(mbvects,state_basis,4,neigs)




with open("report.txt","a") as f:
	f.write(str(V_0))
	f.write("\n")
	f.write(str(L1))
	f.write("\n")
	f.write(str(L2))
	f.write("\n")
	f.write(str(L3))
	f.write("\n")
	f.write(str(L4))
	f.write("\n")
	f.write(str(L5))
	f.write("\n")
	f.write("\n")
	f.write("\n")
	f.write("\n")
	f.write("\n")
	f.write("\n")

print V_0
print L1
print L2
print L3
print L4
print L5



#print getMatrixElement(i,j,2,state_basis,vals,vmat,N) #i,j are indices, k is max difference to consider,N is number of particles



"""
