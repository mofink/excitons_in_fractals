
import math
import cmath
import numpy as np
from scipy.sparse.linalg import eigs
from scipy.linalg import eigh

pi = math.pi

def gen_Fibonacci_seq(w_max,W_A,W_B):
	##W(x)
	W = [0] * w_max
	W[0] = 'A' #W_A
	W[1] = 'B' #W_B
	for i in xrange(w_max):
		if i > 1:
			W[i] = W[i-2] + W[i-1]

	W = list(W[w_max-1])
	M = len(W) #Number of letters

	#convert A and B to numeric values
	
	
	for i in xrange(len(W)):
		if W[i] == 'A':
			W[i] = W_A
		else:
			W[i] = W_B
	
	return W

def Fib_potential(W_Fib, L, a, eta, N, h):
	M = len(W_Fib)

	w = np.zeros(N, dtype=float)
	w_der = np.zeros(N, dtype = float)


	#Recall: M is number of letters; W[0], ... W[M-1]
	#N is number of discrete points; w[0], ... w[N-1]

		
	for k in xrange(N): #number of bins
		total = 0
		for i in xrange(M): #number of letters
			total += W_Fib[i]*math.exp(-(k*h-(i+0.5)*a)**2/(eta*a)**2)
		
		total += W_Fib[0] * math.exp(-(-L + k*h - 0.5 * a)**2/(eta*a)**2) + W_Fib[-1] * math.exp(-(k*h + 0.5 * a)**2/(eta*a)**2)  

		w[k] = total

	#calculate its derivative
	for k in xrange(N): #number of bins
		total = 0
		for i in xrange(M): #number of letters

			total += -2.0*(k*h-(i+0.5)*a)/(eta*a)**2 * W_Fib[i]*math.exp(-(k*h-(i+0.5)*a)**2/(eta*a)**2)
		
		total += -2.0*(-L + k*h - 0.5 * a) / (eta*a)**2 * W_Fib[0] * math.exp(-(-L + k*h - 0.5 * a)**2/(eta*a)**2) - 2.0*(k*h + 0.5*a)/(eta*a)**2 * W_Fib[-1] * math.exp(-(k*h + 0.5 * a)**2/(eta*a)**2)
		w_der[k] = total

	V =  np.zeros(N, dtype = float)
	for k in xrange(N):
		V[k] = ((pi/(w[k]))**2) + ((pi**2 + 3) / 12)*(w_der[k]/w[k])**2

	return V, w

def get_eigsystem(n_max, V, phi, N, h):
	A = np.zeros((n_max,n_max),dtype=np.complex)

	for m in xrange(n_max):
		for n in xrange(m,n_max):
			arg = np.multiply(V, phi[n])
			arg = np.multiply(arg, np.conj(phi[m]))
			A[m, n] = np.trapz(arg, dx = h)
			if m == n:
				A[m, n] = A[m, n] + (2*pi*n)**2
			else:
				A[n, m] = A[m, n]


	#	fname = 'eigs.txt'
	#	vals, vects = eigs(A, k = 20, which = 'SM')
	vals, vects = eigh(A)

	#np.savetxt(fname, np.real(vals))

	return vals,vects
