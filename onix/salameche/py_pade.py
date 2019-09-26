from scipy import linalg

def pade(At, N):

	expAt = linalg.expm(At)

	return expAt.dot(N)