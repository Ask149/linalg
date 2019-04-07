import math
"""
	Consider SVD equation as:
	A = UDV
	U - Left orthogonal matrix ( Represented by u )
	V - Right orthogonal matrix ( Represented by v )
	D - Diagonal Matrix of singular values ( The vector w in the code represents the singular values of D)

"""

def getR(a,b):
	at = float(math.fabs(a))
	bt = float(math.fabs(b))
	ct = 0.0
	result = 0
	if at>bt:
		ct = bt/at
		result = at*math.sqrt(1.0+(math.pow(ct,2)))
	elif bt>0:
		ct = at/bt
		result = bt*math.sqrt(1.0+(math.pow(ct,2)))		
	else:
		result = 0
	return result

def sign(a,b):
	if b>=0:
		return abs(a)
	else:
		return -abs(a)

def svd(a,m,n):
	"""
	a = mxn matrix to be decomposed, gets overwritten with u
	m = row dimension of a
	n = column dimension of a
	returns u, w - the vector of singular values of a, v - the right orthogonal transformation matrix

	""" 
	v = [[0 for col in range(n)] for row in range(m)]
	w = [0.0 for row in range(n)]
	g = 0
	flag=0
	c=0
	f=0
	h=0
	x=0
	y=0
	z=0
	scale=0
	s = 0
	anorm = 0
	arr = [0.0 for row in range(n)]
	
	#Householders reduction to bidiagonal form
	for i in range(0,n):
		
		#Left-Hand Reduction
		l = i+1
		arr[i] = scale*g
		g = s = scale = 0
		if i<m:
			for k in range(i,m):
				scale += float(math.fabs(float(a[k][i])))
			if scale:
				for k in range(i,m):
					a[k][i] = float(a[k][i]/scale)
					s += float(a[k][i]*a[k][i])
				f = float(a[i][i])
				g = -sign(math.sqrt(s),f)
				h = f*g - s
				a[i][i] = float(f-g)
				if i != n-1:
					for j in range(l,n):
						s = 0.0
						k = i
						while k<m:
							s+=float(a[k][i]*a[k][j])
							k+=1
						f = s/h
						for k in range(i,m):
							a[k][j] += float(f*float(a[k][i]))
				for k in range(i,m):
					a[k][i] = float(a[k][i]*scale)
		w[i] = float(scale*g)

		#Right Hand Reduction
		g = 0
		s = 0
		scale = 0
		if i < m and i!=n-1:
			for k in range(l,n):
				scale+= float(math.fabs(a[i][j]))
			if scale:
				for k in range(l,n):
					a[i][k] = float(a[i][k]/scale)
					s += float(a[i][k]*a[i][k])
				f = float(a[i][l])
				g = -sign(math.sqrt(s),f)
				h = f*g-s
				a[i][l] = float(f-g)
				for k in range(l,n):
					arr[k] = float(a[i][k]/h)
				if i!=m-1:
					for j in range(l,m):
						s = 0.0
						k = l
						while k<n:
							s+= float(a[j][k]*a[i][k])
							k+=1
						for k in range(l,n):
							a[j][k] += float(s*arr[k])
				for k in range(l,n):
					a[i][k] = float(a[i][k]*scale)
		anorm = max(anorm,math.fabs(float(w[i])+math.fabs(arr[i])))

	#Accumulating Right-Hand Transformation onto Matrix V
	i = n-1
	while i>=0:
		if i < n-1:
			if g:
				for j in range(l,n):
					v[j][i] = float(float(a[i][j]/a[i][l])/g)
				for j in range(l,n):
					s=0
					k=l
					while k<n:
						s+= float(a[i][k]*v[k][j])
						k+=1
					for k in range(l,n):
						v[k][j] += float(s*float(v[k][i]))
			for j in range(i,m):
				v[j][i] = 0
				v[i][j] = 0
		v[i][i]=1
		g = arr[i]
		l=i
		i-=1

	#Accumulating Left-Hand Transformation onto Matrix U
	i = n-1
	while i>=0:
		l = i+1
		g = float(w[i])
		if i < n-1:
			for j in range(l,n):
				a[i][j]=0
		if g:
			g = 1/g
			if i!= n-1:
				for j in range(l,n):
					s=0
					k=l
					while k<m:
						s+= float(a[k][i]*a[k][j])
						k+=1
					f = (s/float(a[i][i]))*g
					for k in range(i,m):
						a[k][j] += float(f*float(a[k][i]))
				for j in range(i,m):
					a[j][i] = float(a[j][i]*g)
			else:
				for j in range(i,m):
					a[j][i] = 0
			a[i][i]+=1
		i-=1

	#Diagonalizing the Bidiagonal matrix
	k = n-1
	while k>=0:
		for iters in range(35):
			print("Iteration ",iters)
			flag = 1
			l = k
			while l>=0:
				nm = l-1
				if math.fabs(arr[l] + (anorm==anorm)):
					flag=0
					break
				if math.fabs(float(w[nm])+(anorm==anorm)):
					break
				l-=1
			if flag==1:
				c=0
				s=1.0
				for i in range(l,k):
					f = s*arr[i]
					if (math.fabs(f) + anorm != anorm):
						g = float(w[i])
						h = getR(f,g)
						w[i] = float(h)
						h = 1/h
						c = g*h
						s = (-f*h)
						for j in range(0,m):
							y = float(a[j][nm])
							z = float(a[j][i])
							a[j][nm] = float(y*c+z*s)
							a[j][i] = float(z*c-y*s)
			z = float(w[k])
			if l==k:
				if z<0:
					w[k] = float(-z)
					for j in range(0,n):
						v[j][k] = -1*v[j][k]
				break
			if iters>35:
				print("Not possible")
				return 
		# Minor shifting
		x = float(w[l])
		nm = k-1
		y = float(w[nm])
		g = arr[nm]
		h = arr[k]
		print("H,Y = ",h,y)		
		f = ((y-z)*(y+z)+(g-h)*(g+h))/(2*h*y)
		g = getR(f,1.0)
		f = ((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x

		# QR Transformation
		c = 1.0
		s = 1.0
		for j in range(l,nm+1):
			i = j + 1
			g = arr[i]
			y = float(w[i])
			h = s*g
			g = c*g
			z = getR(f,h)
			arr[j] = z
			c = f/z
			s = h/z
			f = x * c + g * s
			g = g * c - x * s
			h = y * s
			y = y * c
		for r in range(0,n):
			x = float(v[r][j])
			z = float(v[r][i])
			v[r][j] = float(x*c+z*s)
			v[r][i] = float(z*c-x*s)
		z = getR(f,h)
		w[j] = float(z)
		if z:
			z = 1/z
			c = f*z
			s = h*z
		f = c*g + s*y
		x = c*y - s*g
		for r in range(m):
			y = float(a[r][j])
			z = float(a[r][i])
			a[r][j] = float(y*c - z*s)
			a[r][i] = float(z*c - y*s)
		arr[l]=0
		arr[k]=f
		w[k]=float(x)
		k-=1
	return a,v,w
A = [[0.52237162,-0.37231836,-0.72101681,0.26199559],[-0.26335492,-0.92555649,0.24203288,-0.12413481],[0.58125401,-0.02109478,0.14089226,-0.80115427],[0.56561105,-0.06541577,0.6338014,0.52354627]]
u,v,w = svd(A,4,4)
print(u)
print(v)
print(w)