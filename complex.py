
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

spacing = numpy.linspace(0,10, 100)
source = 5*numpy.sin(spacing)
plot.plot(spacing, source)

def error_func(coeff, spacing, source):
	predict = coeff[0] * numpy.sin(coeff[1]*spacing)
	error = numpy.sum(numpy.power(predict - source, 2))
	print error
	return error

args = (spacing, source)
output = scipy.optimize.fmin(error_func, numpy.array([1.0,1.0]), args)
print output

plot.show()



