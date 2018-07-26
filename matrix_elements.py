from interaction import find_V_Args

def getMatrixElement(i,j,k,state_basis,vals,V,N):
	if i == j: 
		d = 0
		total = 0

		#lone particle terms
		for particle in xrange(N): 
			total += vals[particle]


		#interaction terms
		kets = find_V_Args(state_basis[i],state_basis[j],N,d,0,0) #m = n = 0

		
		for m in kets:
			for n in kets:
				if m == n:
					pass
				else:
					total += V[m,n,m,n] + V[m,n,n,m] + V[n,m,m,n] + V[n,m,n,m]

		return total

	else:
		m,n,d = compare(state_basis[i],state_basis[j],2)
		if m == n == 0:
			pass
		
		else: # d>0

			if d == 1:

				total = 0

				kets = find_V_Args(state_basis[i],state_basis[j],N,d,m,n)
				for x in kets:
					for other in xrange(N):
						if other == x:
							pass
						else:
							total += V[x,other,x,other] + V[other,x,x,other] + V[x,other,other,x] + V[other,x,other,x]
				return total

			elif d == 2:
				total = 0
				kets = find_V_Args(state_basis[i],state_basis[j],N,d,m,n)
				#easier to handle manually than to rework the generator in find_V_args
				matr = [0,0,0,0]
				indx = 0
				for i in kets:
					matr[indx] = i
					indx = indx + 1
				x,y,z,w = matr[0],matr[1],matr[2],matr[3]

				total += V[x,y,z,w] + V[x,y,w,z] + V[y,x,w,z] + V[y,x,z,w]
				
				return total

			else: #d > 2
				return 0



def compare(a,b,k):
	m,n = [-1,-1],[-1,-1] #m,n store differences by index for a,b respectively
	x,y,d,i,j = 0,0,0,0,0

	while i < len(a) and j < len(b):		 
			if a[i] == b[j]:
				i = i + 1
				j = j + 1

			elif a[i] > b[j]:
				if y == k:
					return -1,-1,-1
				n[y] = j
				y = y + 1
				j = j + 1

			else: 
				if x == k:
					return -1,-1,-1
				m[x] = i
				x = x + 1
				i = i + 1

	while i < len(a):
		m[x] = i
		x = x + 1
		i = i + 1

	while j < len(b):
		n[y] = j
		y = y + 1
		j = j + 1


	d = x
	return m,n,d

