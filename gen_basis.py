
import math
import cmath
import numpy as np
from scipy.misc import factorial


def find_states(N,Ns): #N = number of particles, Ns = number of states, outputs basis matrix

	num = int(factorial(N+Ns-1)/(factorial(N)*factorial(Ns-1)))
	



	mtrx = [[0 for i in xrange(N)] for j in xrange(num)]


	max_idx = 0
	for i in xrange(N):
		if i==N:
			pass
		else:
			max_idx += (Ns-1)*(Ns)**(N-i-1)



	#max_idx = int("".join(str(digit) for digit in max_idx))
	row_idx = 0
	for row in xrange(max_idx+1):
		
		state = genState(row,Ns,N)
	

		if sorted(state) == state:

			
			mtrx[row_idx] = state
			
			row_idx+=1
			#print "ADDED",state
		else:
			#print "REJECTED",state
			pass

	return mtrx
			
def pad(num, N):
	val = str(num).rjust(N)
	return int(val)

#N digit number in base Ns, split by digit into vector. returns vector data type
def genState(num,Ns,N):

	state = [0]*N

	i=0
	while(num>0):

		state[i] = num/Ns**(N-1-i)	
		num = num%(Ns**(N-1-i))
		i=i+1


	return state

