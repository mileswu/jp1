
#!/usr/bin/python
import numpy
import matplotlib.pyplot as plot

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


base_dir = 'ABS_Pod23_241010/rs.04/'

s = read_file(base_dir + 'const_complex_na_0_0.38_0.250/complex.dat')
n = read_file(base_dir + 'const_complex_na_0_0.60_0.250/complex.dat')
m = read_file(base_dir + 'const_complex_na_2_0.38_0.240/complex.dat')
length = len(s[0])

r_n = 10.0

out = numpy.zeros(length, complex)
out2 = numpy.zeros(length, complex)
for x in range(length):
	f = 1.0/(1.0/s[1][x] - 1.0/n[1][x])
	out[x] = r_n * f* (1.0/s[1][x] - 1.0/m[1][x])
	


a = plot.subplot(311)
a.plot(numpy.log10(s[0]), numpy.real(out), 'o')
a = plot.subplot(312)
a.plot(numpy.log10(s[0]), numpy.imag(out), 'o')
a = plot.subplot(313)
a.plot(numpy.real(out), numpy.imag(out), 'o')
#plot.axes().set_aspect('equal')

plot.show()
