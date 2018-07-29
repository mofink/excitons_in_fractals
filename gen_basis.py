
import math
import cmath
import numpy as np
from scipy.misc import factorial


def find_states(N,Ns): #N = number of particles, Ns = number of states, outputs basis matrix

	num = int(factorial(N+Ns-1)/(factorial(N)*factorial(Ns-1)))
	mtrx = np.zeros((num, N), dtype=int)
	Nspow = [Ns**(N-1-i) for i in xrange(N)]


	max_idx = 0
	for i in xrange(N):
		if i==N:
			pass
		else:
			max_idx += (Ns-1)*Nspow[i]


	row_idx = 0
	for row in xrange(max_idx+1):	
		state = genState(row,Ns,N,Nspow)
		if sorted(state) == state:
			mtrx[row_idx] = state
			row_idx+=1

	return mtrx

#wrong implementation
def pad(num, N):
	val = str(num).rjust(N)
	return int(val)

#N digit number in base Ns, split by digit into vector. returns vector data type
def genState(num,Ns,N,Nspow):

	state = [0]*N

	i=0
	while(num>0):

		state[i] = num/Nspow[i]
		num = num%Nspow[i]
		i=i+1


	return state

