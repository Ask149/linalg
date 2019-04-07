def cofactor(A,temp,p,q,n):
	i=0
	j=0
	for row in range(0,n):
		for col in range(0,n):
			if row!=0 and col!=q:
				temp[i][j] = A[row][col]
				j+=1
				if j== (n-1):
					j=0
					i+=1
	return temp

def determinant(A,n):
	D=0
	if n==1:
		return A[0][0]
	temp = [[0 for col in range(n)] for row in range(n)]
	sign = 1
	for i in range(n):
		temp = cofactor(A,temp,0,i,n)
		D += sign*A[0][i]*determinant(temp,n-1)
		sign=-sign
	return D

def adjoint(A):
	n=len(A)
	adj = [[0 for col in range(n)] for row in range(n)]
	if n==1:
		adj[0][0]=1
		return
	else:
		sign = 1
		temp = [[0 for col in range(n)] for row in range(n)]
		for i in range(0,n):
			for j in range(0,n):
				temp = cofactor(A,temp,i,j,n)
				sign = -1
				if (i+j)%2==0:
					sign = 1
				adj[j][i] = (sign*(determinant(temp,n-1)))
	return adj

def inverse(A):
	n=len(A)
	inverse = [[0 for col in range(n)] for row in range(n)]
	det=determinant(A,n)
	if det==0:
		return False
	
	adj = adjoint(A)
	for i in range(0,n):
		for j in range(0,n):
			inverse[i][j] = adj[i][j]/float(det)
	return inverse

def main():
	A = [[1,2,3,4],[-1,-2,-2,-1],[1,2,2,1],[4,3,2,1]]
	#A = [[5,-2,2,7],[1,0,0,3],[-3,1,5,0],[3,-1,-9,4]]
	print(inverse(A))

if __name__ == '__main__':
	main()