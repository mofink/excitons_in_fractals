from interaction import find_V_Args
from collections import Counter
import math

def getMatrixElement(i,j,k,state_basis,vals,V,N):
	if i == j: 
		d = 0
		total = 0

		#lone particle terms
		for state in state_basis[i]: 
			total += vals[state]


		#interaction terms
		kets = find_V_Args(state_basis[i],state_basis[j],N,d,0,0) #m = n = 0

		
		for m,n in kets:
			total += V[m,n,m,n] + V[m,n,n,m] + V[n,m,m,n] + V[n,m,n,m]


	else:
		m,n,d = compare(state_basis[i],state_basis[j],2)
		if d == -1:
			total = 0.0
		
		elif d == 1:
			total = 0
			l1 = state_basis[i,m[0]]
			r1 = state_basis[j,n[0]]

			kets = find_V_Args(state_basis[i],state_basis[j],N,d,m,n)
			for x in kets:
				total += V[l1,x,r1,x] + V[l1,x,x,r1] + V[x,l1,r1,x] + V[x,l1,x,r1]
			

		elif d == 2:
			total = 0
			l1 = state_basis[i,m[0]]
			l2 = state_basis[i,m[1]]
			r1 = state_basis[j,n[0]]
			r2 = state_basis[j,n[1]]
			
			total += V[l1,l2,r1,r2] + V[l1,l2,r2,r1] + V[l2,l1,r1,r2] + V[l2,l1,r2,r1]
		
	
	count1 = Counter(state_basis[i])
	count2 = Counter(state_basis[j])

	for value in count1.values():
		if value == 2:
			total *=math.sqrt(2)
		elif value == 3:
			total *= math.sqrt(3)
		elif value == 4:
			total *=2
		else:
			pass

	for value in count2.values():
		if value == 2:
			total *=math.sqrt(2)
		elif value == 3:
			total *= math.sqrt(3)
		elif value == 4:
			total *=2
		else:
			pass


	return total

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
