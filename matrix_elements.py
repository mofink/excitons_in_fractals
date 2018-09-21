from interaction import find_V_Args
from collections import Counter
import math

def getMatrixElement(i,j,k,state_basis,state_occup,vals,V,N):
	if i == j: 
		d = 0
		total = 0

		#lone particle terms
		for state in state_basis[i]: 
			total += vals[state]


		#interaction terms
		kets = find_V_Args(state_occup[i],state_occup[j],N,d,0,0)

		a = state_occup[i]
		
		for m,n in kets:
			if m != n:
				total += a[m] * a[n] * (V[m,n,m,n] + V[m,n,n,m] + V[n,m,m,n] + V[n,m,n,m])
			else:
				total += a[m] * (a[m] - 1) * V[m,m,m,m] 


	else:
		m,n,d = compare(state_basis[i],state_basis[j],2)
		if d == -1:
			total = 0.0
		
		elif d == 1:
			total = 0

			kets = find_V_Args(state_occup[i],state_occup[j],N,d,m,n)
			a = state_occup[i]
			b = state_occup[j]
			l1 = m[0]
			r1 = n[0]
			for x in kets:
				if x == l1:
					total += math.sqrt(a[x] * (a[x] - 1) * b[x] * b[r1]) * (V[x,x,r1,x] + V[x,x,x,r1])
				elif x == r1:
					total += math.sqrt(a[x] * a[l1] * b[x] * (b[x] - 1)) * (V[l1,x,x,x] + V[x,l1,x,x])
				else:
					total +=math.sqrt(a[x] * a[l1] * b[x] * b[r1]) * (V[l1,x,r1,x] + V[l1,x,x,r1] + V[x,l1,r1,x] + V[x,l1,x,r1])

		elif d == 2:
			total = 0
			a = state_occup[i]
			b = state_occup[j]
			l1 = m[0]
			l2 = m[1]
			r1 = n[0]
			r2 = n[1]

			if l1 == l2 and r1 == r2:
				total += math.sqrt(a[l1] * (a[l1] - 1) * b[r1] * (b[r1] - 1)) * V[l1,l1,r1,r1]
			elif l1 == l2:
				total += math.sqrt(a[l1] * (a[l1] - 1) * b[r1] * b[r2]) * (V[l1,l1,r1,r2] + V[l1,l1,r2,r1])
			elif r1 == r2:
				total += math.sqrt(a[l1] * a[l2] * b[r1] * (b[r1] - 1)) * (V[l1,l2,r1,r1] + V[l2,l1,r1,r1])
			else:
				total += math.sqrt(a[l1] * a[l2] * b[r1] * b[r2]) * V[l1,l2,r1,r2] + V[l1,l2,r2,r1] + V[l2,l1,r1,r2] + V[l2,l1,r2,r1]
		
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
				n[y] = b[j]
				y = y + 1
				j = j + 1

			else: 
				if x == k:
					return -1,-1,-1
				m[x] = a[i]
				x = x + 1
				i = i + 1

	while i < len(a):
		m[x] = a[i]
		x = x + 1
		i = i + 1

	while j < len(b):
		n[y] = b[j]
		y = y + 1
		j = j + 1


	d = x
	return m,n,d
