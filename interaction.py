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
				for l in xrange(num_states if k < i else j + 1):

					arg = 1.5*psi[i].conj()*psi[j].conj()*psi[k]*psi[l]/w
					result = np.trapz(arg, dx = h)*V_0
					
					mtrx[i,j,k,l] = result
					mtrx[j,i,l,k] = result
					result = result.conjugate()
					mtrx[k,l,i,j] = result
					mtrx[l,k,j,i] = result

	
	return mtrx


#use new function instead of this
def find_interaction_matrix(A,V_0,psi,w,num_states):

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
					mtrx[i,j,k,l] = result
					mtrx[j,i,l,k] = result
					result = result.conjugate()
					mtrx[k,l,i,j] = result
					mtrx[l,k,j,i] = result

	return mtrx

	
def find_V_Args(a,b,N,d,m,n):
	
	if d == 0: #Here a == b
		corpus = Counter(a)
		keys = corpus.keys()
		values = corpus.values()

		for i in xrange(len(keys)):
			for j in xrange(i,len(keys)):
				if i == j:
					
					if values[i] < 2:
						pass
					else:
						yield keys[i],keys[j]

				else:
					yield keys[i],keys[j]


	elif d == 1:


		a_dict = Counter(a[0:m[0]] + a[m[0]+1:])
		a_keys = a_dict.keys()
	
		
		#Recall: m,n are differences between two states by index
		for i in a_keys:
			yield i


	else: #d == 2
		
		flag = True
		while flag == True:
			for i in m:
				yield a[i]
			flag = False
		
		for i in n:
			yield b[i]