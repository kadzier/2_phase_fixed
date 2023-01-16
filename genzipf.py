N = 100000
s = 3
l = 8861572
qi = []

def genZipf(N, a):
	cursum = 0	
	for i in range(N):
		qi.append(1/((i+1)**a))
		cursum += qi[i]

	for i in range(N):
		qi[i] /= cursum

def eta(N, l):
	sum = 0
	for i in range(N):
		sum += qi[i] * (1 - qi[i])**(l-1)
	return sum


genZipf(N,s)
for i in range(10):
	ind = i + int(N*i/10)
	print("i =",ind,":",qi[ind])
print(eta(N,l))