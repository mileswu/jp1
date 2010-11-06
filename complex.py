
#!/usr/bin/python
import numpy
import matplotlib.pyplot as plot

file = open('complex.dat', 'r')

length = -1
while file.readline():
	length += 1

freq = numpy.zeros(length)
com = numpy.zeros(length, complex)

file.seek(0)
file.readline()
for i in xrange(length):
	line = file.readline()
	split_line = line.split()

	freq[i] = split_line[0]
	com[i] = 1.0 / complex(float(split_line[1]), float(split_line[2]))

plot.plot(numpy.real(com), numpy.imag(com))

plot.show()
