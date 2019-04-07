from sklearn.decomposition import PCA
import numpy as np
import numpy.linalg as alg
# define a matrix
def pca_eigen(A):
	M = np.mean(A.T, axis=1)
	C = A - M
	V = np.cov(C.T)
	values, vectors = alg.eig(V)
	P = vectors.T.dot(C.T)
	return P.T

def pca_svd(A,numcomponents):
	A = A - np.mean(A,axis=0)
	[U,D,V] = alg.svd(A)
	V = V.T
	V = V[:,:numcomponents]
	return np.dot(A,V)

def pca(A,numcomponents):
	pca = PCA(numcomponents)
	pca.fit(A)
	B = pca.transform(A)
	return B
def main():
	A = [[-11, 2], [9, 4], [-1, 8]]

	print("PCA Eigen :")
	PCA_Eigen = pca_eigen(np.array(A))	
	print(PCA_Eigen)

	print("PCA SVD :")
	PCA_SVD   = pca_svd(np.array(A),numcomponents=2)
	print(PCA_SVD)

	"""print("PCA In-Built")
	PCA = pca(np.array(A),numcomponents=2)
	print(PCA)
	"""
if __name__ == '__main__':
	main()