
import math
import cmath
import numpy as np
from scipy.sparse.linalg import eigs
from scipy.linalg import eigh

pi = math.pi

def calc_A_matrix(w_max,W_A,W_B,m_max,n_max):

	##W(x)
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
	phi = np.array([[cmath.exp(1j*(2*pi*n/L)*k*h) for k in xrange(N)] for n in xrange(n_max)], dtype=np.complex)

	for m in xrange(m_max):
		for n in xrange(n_max):
			
			arg = np.multiply(phi[n],V)
			arg = np.multiply(arg,np.conj(phi[m]))
			A[m, n] = np.trapz(arg, dx = h)
			if m == n:
				A[m, n] = A[m, n] + (2*pi*n)**2





	fname = 'eigs.txt'
	#	vals, vects = eigs(A, k = 20, which = 'SM')
	vals, vects = eigh(A)

	np.savetxt(fname, np.real(vals))

	return vals,vects,A,phi,w,h