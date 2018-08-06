from collections import Counter
import numpy as np

def occupation_numbers(state_basis,state_num):
	num_of_times = np.zeros(len(state_basis),dtype=int)
	for i in xrange(len(state_basis)):
		state = state_basis[i]
		state_dict = Counter(state)

		num_of_times[i] = state_dict.get(state_num , 0)

	return num_of_times

def occupation(mbvects,state_basis,state_num,num_states):
	mbvects_t = mbvects.transpose()
	occup = np.zeros(num_states,dtype=complex)
	B = occupation_numbers(state_basis,state_num)
	for i in xrange(num_states):
		occup[i] = np.dot(mbvects_t[i].conj(),np.multiply(B,mbvects_t[i]))

	return np.real(occup)	

