import matplotlib.pyplot as plt
from vals_vects import *
import numpy as np

with open('report.txt') as f:
	words = f.read().split()
	for point in words:
		y_vals = [float(x) for x in words]



w_max = 13 # Sequence order
L = 1.0
W_A = 0.01
W_B = 0.018
n_max = 200
W_Fib = gen_Fibonacci_seq(w_max,W_A,W_B)
M = len(W_Fib)
a = L / M
N = 1000 # number of discrete points
eta = 0.204
h = L/N #size of step


EnS = 1.05**2*100/(2*4.5*1.6*9.1)
aex = 0.8

for i in xrange(len(y_vals)):
	y_vals[i] = y_vals[i] * EnS / (aex*M)**2

inc = 10 #increment when generating data

x = np.linspace(0,len(y_vals)/2,50)
plt.scatter(x, y_vals)
plt.show()
