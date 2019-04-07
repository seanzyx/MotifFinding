import numpy as np
from scipy.optimize import fsolve

def entropy(x):
	if np.abs(x)<0.000001:
		x = 0
	if x == 0:
		return 0
	if x < 0 or x > 1:
		# print('error', x)
		return 'error'
	return x*np.log2(x/0.25)

def calculateICPC(a, b, c):
	# print(a,b,c)
	return entropy(a)+entropy(b)+entropy(c)+entropy(1-a-b-c)

def equation(c, a, b, ICPC):
	return (calculateICPC(a, b, c) - ICPC)

def randomAB():
	a = np.random.rand()
	b = np.random.rand()*(1-a)
	return a, b

def solveC(a, b, ICPC):
	criterion = ICPC - 2 - a*np.log2(a) - b * np.log2(b)
	if criterion < (1-a-b) * np.log2((1-a-b)/2) or criterion > (1-a-b)*(np.log2(1-a-b)):
		# print('a b will not fit ICPC', (1-a-b) * np.log2((1-a-b)/2), criterion, (1-a-b)*(np.log2(1-a-b)), a, b, ICPC)
		return -1
	else:
		guess = (1-a-b)/2
	try:
		c = fsolve(equation, guess, args=(a, b, ICPC), factor = 1,xtol=1.49012e-20)
	except:
		# print('no sol')
		return -1
	else:
		# print(a, b, c, 1-a-b-c)
		# print(calculateICPC(a,b,c))
		return c

def getColumn(ICPC):
	if ICPC == 2:
		return np.array([0,0,0,1])
	if ICPC == 0:
		return np.array([0.25]*4)
	while True:
		a, b = randomAB()
		result = solveC(a, b, ICPC)
		if result == -1:
			continue
		return np.array([a,b,result,1-a-b-result])

def getPWM(ICPC, ML):
	pwm = np.zeros([ML, 4])
	for i in range(ML):
		column = getColumn(ICPC)
		column = np.random.permutation(column)
		pwm[i, :] = column
	return pwm