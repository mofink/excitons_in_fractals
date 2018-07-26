def find_IDOS(vals):
	step_size = (vals[len(vals)-1] - vals[0])/1000
	steps = np.arange(np.floor(vals[0]),np.ceil(vals[len(vals)-1]),step=step_size)
	
	eps = list()
	for i in xrange(len(steps)):
		eps.append(0)

	eps[0] = vals[0]
	for i in xrange(len(eps)):
		if i == 0:
			pass
		else:

			eps[i] = eps[i-1] + step_size

	IDOS = list()
	for i in xrange(len(eps)):
		IDOS.append(0)

	index = 0
	for upper in eps:
		total = 0
		
		for energy_level in vals:
			if energy_level < upper:
				total += 1
			else:
				pass
		
		IDOS[index] = float(total)/float(len(vals))
		index += 1



	return eps,IDOS
