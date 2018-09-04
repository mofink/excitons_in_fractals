import numpy as np
from scipy.sparse import csc_matrix
from matrix_elements import *
from multiprocessing import *

def createSparseMatrix(state_basis,vals,vmat,N):
	dim = 5*len(state_basis)
	basis_size = len(state_basis)

	k = 0
	row_ind = np.zeros(5*dim,dtype=int)
	col_ind = np.zeros(5*dim,dtype=int)
	data = np.zeros(5*dim,dtype=complex)

	for i in xrange(basis_size):
		for j in xrange(i,basis_size):
			matrix_element = getMatrixElement(i,j,2,state_basis,vals,vmat,N)
			if abs(matrix_element) < 10**-7:
				pass
			else:
				if k == len(row_ind):
					row_ind.resize(len(row_ind)*2)
					col_ind.resize(len(col_ind)*2)
					data.resize(len(data)*2)

				row_ind[k] = i
				col_ind[k] = j
				if i == j:
					data[k] = matrix_element/2.0
				else:
					data[k] = matrix_element

				k = k + 1
	
	row_ind.resize(k)
	col_ind.resize(k)
	data.resize(k)

	return csc_matrix((data, (row_ind, col_ind)), shape=(basis_size, basis_size),dtype = complex)

