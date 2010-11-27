
#!/usr/bin/python
import numpy
import matplotlib.pyplot as plot
import scipy.optimize

def read_file(f):
	file = open(f, 'r')

	length = -1
	file.readline()
	while 1:
		line = file.readline()
		if line == '':
			break
		split_line = line.split()
		if(float(split_line[0]) < 5000):
			length += 1

	#length = 650

	freq = numpy.zeros(length)
	com = numpy.zeros(length, complex)

	file.seek(0)
	file.readline()
	for i in xrange(length):
		line = file.readline()
		split_line = line.split()

		freq[i] = split_line[0]
		com[i] = -1.0/complex(float(split_line[1]), float(split_line[2]))
	return (freq, com)

def transfer_func(s, n, m):
	length = len(s[0])
	
	r_n = 10.0
	
	out = numpy.zeros(length, complex)
	for x in range(length):
		f = 1.0/(1.0/s[1][x] - 1.0/n[1][x])
		out[x] = r_n * f* (1.0/s[1][x] - 1.0/m[1][x])
	return out

def plot_complex(freq, data):
	a = plot.subplot(311)
	a.plot(numpy.log10(freq), numpy.real(data), 'o')
	a = plot.subplot(312)
	a.plot(numpy.log10(freq), numpy.imag(data), 'o')
	a = plot.subplot(313)
	a.plot(numpy.real(data), numpy.imag(data), 'o')
	#plot.axes().set_aspect('equal')

	plot.show()

def predict_func_matrix(coeff, spacing):
	for i in coeff:
		if i < 0:
			return 100000;
	

	r_0 = coeff[0]
	b_i = coeff[1]
	l_i = coeff[2]
	l = 1.0 #coeff[3]
	
	g = coeff[3]
	c = coeff[4]
	
	i_0 =  1.0 #coeff[6]
	r_l = 1.0 #coeff[7]
	
	g1 = coeff[5]
	c1 = coeff[6]

	t_i = c/(g*(1-l_i))
	i = (0+1j)

	predict = numpy.zeros(len(spacing), complex)
	loop = 0
	for w in spacing:
		#matrix = numpy.array([
		#	[(r_l + r_0*(1 + b_i))/l + i*w, l_i*g/i_0/l],
		#	[-i_0*r_0*(2+b_i)/c, 1/t_i + i*w]
		#])
		#matrix = numpy.array([
		#	[(r_l + r_0*(1 + b_i))/l + i*w, l_i*g/i_0/l, 0],
		#	[-i_0*r_0*(2+b_i)/c, 1/t_i + i*w, g1/c1],
		#	[0, -g1/c1, i*w]
		#])
		matrix = numpy.array([
			[(r_l + r_0*(1 + b_i))/l + i*w, l_i*g/i_0/l, 0],
			[-i_0*r_0*(2+b_i)/c, g/c -g*l_i/c + i*w, -g/c],
			[0, -g/c1, g/c1 + g1/c1 + i*w]
		])
		inv = numpy.linalg.inv(matrix)
		predict[loop] = l/inv[0][0] - r_l - i*w*l

		loop += 1
	
	return predict

def predict_func(coeff, spacing):
	r0 = coeff[0]
	b_i = coeff[1]
	l_i = coeff[2]
	t_i = coeff[3]

	i = (0+1j)
	predict = r0*(1+b_i) + (r0*l_i)/(1-l_i)*(2+b_i)/(1+i*t_i*spacing)
	return predict

def error_func(coeff, spacing, source):
	predict = predict_func_matrix(coeff, spacing)
	error = numpy.sum(numpy.power(numpy.abs(predict - source), 2))
	print(error)
	return error

def fit_nist(s, out):
	args = (s[0], out)
	coeff = scipy.optimize.fmin_powell(error_func, numpy.array([10.0, 1.0, 15.0, 1000.0, 41.0, 1.0, 1.0]), args)
	#coeff = numpy.array([391.0, 0.9, 0.004, 878.0, 285.0])
	#coeff = numpy.array([10.0, 0.9, 15.0, 3857.0, 45.0])
	print "Coeffecients: ", coeff
	print "Error: ", error_func(coeff, s[0], out)
	predict = predict_func_matrix(coeff, s[0])

	a = plot.subplot(311)
	a.plot(numpy.log10(s[0]), numpy.real(out), 'o')
	a.plot(numpy.log10(s[0]), numpy.real(predict), 'o')
	a = plot.subplot(312)
	a.plot(numpy.log10(s[0]), numpy.imag(out), 'o')
	a.plot(numpy.log10(s[0]), numpy.imag(predict), 'o')
	a = plot.subplot(313)
	a.plot(numpy.real(out), numpy.imag(out), 'o')
	a.plot(numpy.real(predict), numpy.imag(predict), 'o')
	plot.show()



base_dir = 'ABS_Pod23_241010/rs.04/'

s = read_file(base_dir + 'const_complex_na_0_0.38_0.250/complex.dat')
n = read_file(base_dir + 'const_complex_na_0_0.60_0.250/complex.dat')
m = read_file(base_dir + 'const_complex_na_2_0.38_0.240/complex.dat')

out = transfer_func(s, n, m)

#plot_complex(s[0], out)
fit_nist(s, out)

#Coeffecients:  [  5.98340389e+00   6.23565491e-01   1.61887268e-01   1.52429108e+04
  # 5.61840916e+00   1.75209612e+01   2.07340763e+01]

