import numpy as np
from scipy.sparse import csc_matrix
from matrix_elements import *
from multiprocessing import Pool,sharedctypes,RawArray
import itertools
from functools import partial
import warnings


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
"""
def fill_per_window(idxs,state_basis,vals,vmat,N,block_size,shared_array):


	window_x, window_y = idxs
	tmp = np.ctypeslib.as_array(shared_array)

	for idx_x in xrange(window_x, window_x + block_size):
		for idx_y in xrange(window_y, window_y + block_size):
			tmp[idx_x, idx_y] = getMatrixElement(idx_x,idx_y,2,state_basis,vals,vmat,N)



def create_matrix(state_basis,vals,vmat,N):
	dim = len(state_basis)
	block_size = 4
	result = np.ctypeslib.as_ctypes(np.zeros((dim, dim)))
	shared_array = sharedctypes.RawArray(result._type_, result)

	idxs = [(i, j) for i, j in itertools.product(xrange(0, dim,block_size),xrange(0, dim,block_size))]
	func = partial(fill_per_window,state_basis=state_basis,vals=vals,vmat=vmat,N=N,block_size=block_size,shared_array=shared_array)
	p = Pool()
	res = p.map(func, idxs)
	result = np.ctypeslib.as_array(shared_array)

	return result

"""
def fill_per_window(idx,state_basis,vals,vmat,N,block_size,shared_array):

	tmp = shared_array

	dim = len(state_basis)
	

	for idx_x in xrange(idx[0], idx[1]):
		for idx_y in xrange(0, dim):
			tmp[idx_x, idx_y] = getMatrixElement(idx_x,idx_y,2,state_basis,vals,vmat,N)



def create_matrix(state_basis,vals,vmat,N):


	dim = len(state_basis)
	block_size = 4
	shared_array = np.zeros((dim, dim), dtype=complex)
	

	idxs = [(i, min(i + block_size, dim)) for i in xrange(0, dim, block_size)]
	func = partial(fill_per_window,state_basis=state_basis,vals=vals,vmat=vmat,N=N,block_size=block_size,shared_array=shared_array)

	for idx in idxs:
		func(idx)

	#p = Pool()
	#res = p.map(func, idxs)
	result = shared_array

	return result

def create_multiproc_matrix(state_basis,vals,vmat,N):
	global state_basis_sh
	global vals_sh
	global vmat_sh
	global matrix_sh
	with warnings.catch_warnings():
		warnings.simplefilter("ignore")
		state_basis_ct = np.ctypeslib.as_ctypes(state_basis)
		vals_ct = np.ctypeslib.as_ctypes(vals)
		vmat_ct = np.ctypeslib.as_ctypes(vmat.view(dtype='float64'))
		state_basis_sh = RawArray(state_basis_ct._type_, state_basis_ct)
		vals_sh = RawArray(vals_ct._type_, vals_ct)
		vmat_sh = RawArray(vmat_ct._type_, vmat_ct)

		dim = len(state_basis)
		block_size = dim / 20
		matrix = np.zeros((dim, dim), dtype=complex)
		matrix_ct = np.ctypeslib.as_ctypes(matrix.view(dtype='float64'))
		matrix_sh = RawArray(matrix_ct._type_, matrix_ct)

		func = partial(fill_per_window_sh,N=N,block_size=block_size)
		idxs = [(i, min(i + block_size, dim)) for i in xrange(0, dim, block_size)]
	
		p = Pool(processes=10)
		res = p.map(func, idxs)
		np.copyto(matrix, np.ctypeslib.as_array(matrix_sh).view(dtype='complex128'))

	return matrix


def fill_per_window_sh(idx,N,block_size):

	tmp = np.ctypeslib.as_array(matrix_sh).view(dtype='complex128')
	state_basis = np.ctypeslib.as_array(state_basis_sh)
	vals = np.ctypeslib.as_array(vals_sh)
	vmat = np.ctypeslib.as_array(vmat_sh).view(dtype='complex128')

	dim = len(state_basis)

	for idx_x in xrange(idx[0], idx[1]):
		for idx_y in xrange(0, dim):
			tmp[idx_x, idx_y] = getMatrixElement(idx_x,idx_y,2,state_basis,vals,vmat,N)



