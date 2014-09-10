import numpy as np

def NMI(z1, z2):
	N = len(z1)
	assert N == len(z2)

	p1 = np.bincount(z1)/N
	p1[p1 == 0] = 1
	H1 = (-p1*np.log(p1)).sum()

	p2 = np.bincount(z2)/N
	p2[p2 == 0] = 1
	H2 = (-p2*np.log(p2)).sum()

	joint = np.histogram2d(z1,z2,[range(0,z1.max()+2), range(0,z2.max()+2)],normed=True)
	joint_p = joint[0]
	pdiv = joint_p/np.outer(p1,p2)
	pdiv[joint_p == 0] = 1
	MI = (joint_p*np.log(pdiv)).sum()

	if MI == 0:
		NMI = 0
	else:
		NMI = MI/np.sqrt(H1*H2)

	return NMI

if __name__ == '__main__':
	print('NMI test: ')
	print('   [0 1 2 3], [0 1 2 3]: ' + str(NMI(np.array([0,1,2,3]),np.array([0,1,2,3]))))
	print('   [0 1 2 3], [1 1 2 5]: ' + str(NMI(np.array([0,1,2,3]),np.array([1,1,2,5]))))