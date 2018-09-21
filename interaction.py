from collections import Counter
import math
import cmath
import numpy as np



def new_find_interaction_matrix(A,V_0,phi,w,num_states,h):

	A = A.transpose()
	psi = np.matmul(A[0:num_states,:],phi)

	mtrx = np.zeros(4 * (num_states,),dtype=complex)

	for i in xrange(num_states):
		for j in xrange(i + 1):
			for k in xrange(i + 1):
				for l in xrange(k+1):

					arg = 1.5*psi[i].conj()*psi[j].conj()*psi[k]*psi[l]/w
					result = np.trapz(arg, dx = h)*V_0

					result *= 0.5
					
					mtrx[i,j,k,l] = result
					mtrx[j,i,k,l] = result
					mtrx[i,j,l,k] = result
					mtrx[j,i,l,k] = result
					result = result.conjugate()
					mtrx[k,l,i,j] = result
					mtrx[l,k,i,j] = result
					mtrx[k,l,j,i] = result
					mtrx[l,k,j,i] = result

	
	return mtrx


#use new function instead of this
def find_interaction_matrix(A,V_0,psi,w,num_states,h):

	int_mtrx = np.zeros(4 * (n_max,), dtype=complex)

	for i in xrange(n_max):
		for j in xrange(i + 1):
			for k in xrange(i + 1):
				for l in xrange(k + 1):
					arg = 1.5*psi[i].conj()*psi[j].conj()*psi[k]*psi[l]/w
					result = np.trapz(arg, dx = h)
					int_mtrx[i,j,k,l] = result
					int_mtrx[j,i,k,l] = result
					int_mtrx[i,j,l,k] = result
					int_mtrx[j,i,l,k] = result
					result = result.conjugate()
					int_mtrx[k,l,i,j] = result
					int_mtrx[l,k,i,j] = result
					int_mtrx[k,l,j,i] = result
					int_mtrx[l,k,j,i] = result

	A = A.transpose()

	mtrx = np.zeros(4 * (num_states,),dtype=complex)

	for i in xrange(num_states):
		for j in xrange(i + 1):
			for k in xrange(i + 1):
				for l in xrange(num_states if k < i else j + 1):
					result = 0
					for m in xrange(n_max):
						for n in xrange(n_max):
							for o in xrange(n_max):
								for p in xrange(n_max):
									result += V_0*A[i,m].conjugate()*A[j,n].conjugate()*A[k,o]*A[l,p]*int_mtrx[m,n,o,p]
					result *= 0.5
					mtrx[i,j,k,l] = result
					mtrx[j,i,l,k] = result
					result = result.conjugate()
					mtrx[k,l,i,j] = result
					mtrx[l,k,j,i] = result

	return mtrx

	
def find_V_Args(a,b,N,d,m,n):
	
	if d == 0: #Here a == b
		for i in xrange(len(a)):
			if a[i] != 0:
				for j in xrange(i, len(a)):
					if i == j:
						if a[i] > 1:
							yield i, i
					else:
						if a[j] != 0:
							yield i, j


	elif d == 1:		
		#Recall: m,n are differences between two states by index
		for i in xrange(len(a)): 
			if i == m[0] and a[i] > 1:
				yield i
			elif i == n[0] and b[i] > 1:
				yield i
			elif a[i] != 0 and b[i] != 0:
				yield i
