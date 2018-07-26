import numpy as np
from scipy.sparse import csr_matrix
from matrix_elements import *

def createSparseMatrix(state_basis,vals,vmat,N):
	dim = 5*len(state_basis)

	k = 0
	row_ind = np.zeros(5*dim,dtype=int)
	col_ind = np.zeros(5*dim,dtype=int)
	data = np.zeros(5*dim,dtype=complex)

	for i in xrange(len(state_basis)):
		for j in xrange(i,len(state_basis)):
			matrix_element = getMatrixElement(i,j,2,state_basis,vals,vmat,N)
			if abs(matrix_element) < 10**-7:
				pass
			else:

				if k > len(row_ind):
					row_ind.resize(len(row_ind)*2)
					col_ind.resize(len(col_ind)*2)
					data.resize(len(data)*2,dtype=complex)



				row_ind[k] = i
				col_ind[k] = j
				data[k] = matrix_element
				if i == j:
					data[k] = matrix_element/2

				k = k + 1

	return csr_matrix((data, (row_ind, col_ind)), shape=(k-1, k-1),dtype = complex)

