
#!/usr/bin/python
import numpy
import matplotlib.pyplot as plot
import scipy.optimize

def read_file(f):
	file = open(f, 'r')

	length = -1
	while file.readline():
		length += 1

	#length = 531

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


base_dir = 'ABS_Pod23_241010/rs.04/'

s = read_file(base_dir + 'const_complex_na_0_0.38_0.250/complex.dat')
n = read_file(base_dir + 'const_complex_na_0_0.60_0.250/complex.dat')
m = read_file(base_dir + 'const_complex_na_2_0.38_0.240/complex.dat')

out = transfer_func(s, n, m)

#plot_complex(s[0], out)

def predict_func(coeff, spacing):
	r0 = coeff[0]
	b_i = coeff[1]
	l_i = coeff[2]
	t_i = coeff[3]

	i = (0+1j)
	predict = r0*(1+b_i) + (r0*l_i)/(1-l_i)*(2+b_i)/(1+i*t_i*spacing)
	return predict

def error_func(coeff, spacing, source):
	predict = predict_func(coeff, spacing)
	error = numpy.sum(numpy.power(numpy.abs(predict - source), 2))
	return error

args = (s[0], out)
coeff = scipy.optimize.fmin_powell(error_func, numpy.array([10.0,0.0,0.0,0.0]), args)
print "Coeffecients: ", coeff
print "Error: ", error_func(coeff, s[0], out)
predict = predict_func(coeff, s[0])

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



