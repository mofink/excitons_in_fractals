#! Outputs Matrix A for Fibb 3-D Potential Well

import math
import cmath
import numpy as np
from scipy.sparse.linalg import eigs
from scipy.linalg import eigh


import matplotlib.pyplot as plt


pi = math.pi

# Lets define everything in length of L
w_max = 13 # Sequence order

W_A = 0.01
W_B = 0.0188

X = []
Y = []

#####################




m_max = n_max = 10

#W(x)
W = list()
for i in xrange(w_max): #w_max is sequence order
	W.append(0)
W[0] = 'A' #W_A
W[1] = 'B' #W_B



for i in xrange(w_max):
	if i > 1:
		W[i] = W[i-2] + W[i-1]
	else:
		pass
W = list(W[w_max-1])

#W=25*['A','B']


M = len(W) #Number of letters
L = 1.0
a = L / M

#L = a * M, total length. Each letter is length a


#convert A and B to numeric values
for i in xrange(len(W)):
	if W[i] == 'A':
		W[i] = W_A
	else:
		W[i] = W_B


######## W[i] is a discrete value for each letter 

N = 1000 # number of discrete points



eta = 0.204


h = 1.0/N #size of step

w = list()
for i in xrange(N):
	w.append(0)

	




#Recall: M is number of letters; W[0], ... W[M-1]
#N is number of discrete points; w[0], ... w[N-1]


	
for k in xrange(N): #number of bins
	total = 0
	for i in xrange(M): #number of letters
		total += W[i]*math.exp(-(k*h-(i+0.5)*a)**2/(eta*a)**2)
	
	total += W[0] * math.exp(-(-L + k*h - 0.5 * a)**2/(eta*a)**2) + W[-1] * math.exp(-(k*h + 0.5 * a)**2/(eta*a)**2)  

	w[k] = total

#calculate its derivative
w_ = list()
for i in xrange(N):
	w_.append(0)

for k in xrange(N): #number of bins
	total = 0
	for i in xrange(M): #number of letters

		total += -2.0*(k*h-(i+0.5)*a)/(eta*a)**2 * W[i]*math.exp(-(k*h-(i+0.5)*a)**2/(eta*a)**2)
	
	total += -2.0*(-L + k*h - 0.5 * a) / (eta*a)**2 * W[0] * math.exp(-(-L + k*h - 0.5 * a)**2/(eta*a)**2) - 2.0*(k*h + 0.5*a)/(eta*a)**2 * W[-1] * math.exp(-(k*h + 0.5 * a)**2/(eta*a)**2)
	w_[k] = total

V = list()
for k in xrange(N):

	val = ((pi/(w[k]))**2) + ((pi**2 + 3) / 12)*(w_[k]/w[k])**2

	V.append(val)


A = np.zeros((m_max,n_max),dtype=np.complex)
psi = np.array([[cmath.exp(1j*(2*pi*n/L)*k*h) for k in xrange(N)] for n in xrange(n_max)], dtype=np.complex)

for m in xrange(m_max):
	for n in xrange(n_max):
		
		arg = np.multiply(psi[n],V)
		arg = np.multiply(arg,np.conj(psi[m]))
		A[m, n] = np.trapz(arg, dx = h)
		if m == n:
			A[m, n] = A[m, n] + (2*pi*n)**2





fname = 'eigs.txt'
#	vals, vects = eigs(A, k = 20, which = 'SM')
vals, vects = eigh(A)
np.savetxt(fname, np.real(vals))
EnS = 1.05**2*100/(2*4.5*1.6*9.1)
aex = 0.8

x = xrange(n_max)
"""
f, (ax1, ax2) = plt.subplots(1, 2)
ax1.plot(w)
ax2.scatter(x,vals * EnS / (aex*M)**2)
ax2.set_xlim(0, 100)
#ax2.set_ylim(40.0, 50.0)
plt.show()
"""
def find_IDOS(vals):
	step_size = (vals[len(vals)-1] - vals[0])/1000
	steps = np.arange(np.floor(vals[0]),np.ceil(vals[len(vals)-1]),step=step_size)
	
	eps = list()
	for i in xrange(len(steps)):
		eps.append(0)

	eps[0] = vals[0]
	for i in xrange(len(eps)):
		if i == 0:
			pass
		else:

			eps[i] = eps[i-1] + step_size

	IDOS = list()
	for i in xrange(len(eps)):
		IDOS.append(0)

	index = 0
	for upper in eps:
		total = 0
		
		for energy_level in vals:
			if energy_level < upper:
				total += 1
			else:
				pass
		
		IDOS[index] = float(total)/float(len(vals))
		index += 1



	return eps,IDOS


def find_interaction_matrix(A,V_0,psi,w,num_states):
	w = np.array(w)
	psi = np.array(psi)

	int_mtrx = np.zeros((len(psi),)*4,dtype=complex)


	for i in xrange(len(psi)):
		for j in xrange(len(psi)):
			for k in xrange(len(psi)):
				for l in xrange(len(psi)):
					arg = ((12*pi+9)/(32*pi))*psi[i].conj()*psi[j].conj()*psi[k]*psi[l]*4/w
					int_mtrx[i,j,k,l] = np.trapz(arg, dx = h)

	A = A.transpose()

	mtrx = np.zeros((num_states,)*4,dtype=complex)

	for i in xrange(num_states):
		for j in xrange(num_states):
			for k in xrange(num_states):
				for l in xrange(num_states):
					for m in xrange(len(psi)):
						for n in xrange(len(psi)):
							for o in xrange(len(psi)):
								for p in xrange(len(psi)):
									mtrx[i,j,k,l] = V_0*A[i,m]*A[j,n]*A[k,o]*A[l,p]*int_mtrx[m,n,o,p]
									print i,j,k,l,m,n,o,p

									
									
				

	return mtrx



find_interaction_matrix(vects,1,psi,w,2)




