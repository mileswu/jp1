
#!/usr/bin/python
import numpy
import matplotlib.pyplot as plot
import scipy.optimize
import os
import re

def read_file(f):
	file = open(f + '/complex.dat', 'r')

	length = -1
	file.readline()
	while 1:
		line = file.readline()
		if line == '':
			break
		split_line = line.split()
		if(float(split_line[0]) < 9000):
			length += 1

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

def plot_complexn(freq, data):
	titles = ["10% Rn","50% Rn","90% Rn"]

	plot.subplots_adjust(wspace=0.4, hspace=0.5)
	a = plot.subplot(121)
	#plot.title("Z_TES")
	for i in data:
		a.semilogx(freq, numpy.real(i), ',')
	plot.xlabel("Frequency /Hz")
	plot.ylabel("Re(Z) /mOhms")
	plot.xlim(xmin=1)

	a = plot.subplot(122)
	for i in data:
		a.semilogx(freq, numpy.imag(i), ',')
	plot.xlabel("Frequency /Hz")
	plot.ylabel("Im(Z) /mOhms")
	plot.xlim(xmin=1)
	plot.show()
	
	for i in data:
		plot.plot(numpy.real(i), numpy.imag(i), ',')
	plot.axes().set_aspect('equal')
	plot.xlabel("Re(Z) /mOhms")
	plot.ylabel("Im(Z) /mOhms")
	
	plot.legend(titles, numpoints=100, loc=2, ncol=3, bbox_to_anchor=(-0.0, -0.2))
	plot.show()

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
	predict = []
	coeffs = []
	errors = []
	
	for i in out:
		args = (s[0], i)
		#coeff = scipy.optimize.fmin_powell(error_func, numpy.array([10.0, 1.0, 15.0, 1000.0, 41.0]), args)
		coeff = scipy.optimize.fmin_powell(error_func, numpy.array([10.0, 1.0, 15.0, 1000.0, 41.0, 1.0, 1.0]), args)
		predict.append( predict_func_matrix(coeff, s[0]) )
		coeffs.append(coeff)
		errors.append(error_func(coeff, s[0], i))
	
	return predict, coeffs, errors

def plot_fit(original, predictions):
	titles = ["10% Rn","50% Rn","90% Rn"]
	for i in original:
		plot.plot(numpy.real(i), numpy.imag(i), ',')
	for i in predictions:
		plot.plot(numpy.real(i), numpy.imag(i), 'k-')
	plot.axes().set_aspect('equal')
	plot.xlabel("Re(Z) /mOhms")
	plot.ylabel("Im(Z) /mOhms")
	
	plot.legend(titles, numpoints=100, loc=2, ncol=3, bbox_to_anchor=(-0.0, -0.2))
	plot.show()



base_dir = 'ABS_Pod15_081110/rs.01/'

files = []
temps = []
for i in os.listdir(base_dir):
	match = re.match("const_complex_na_(.*?)_(.*?)_", i)
	if match:
		files.append(i)
		temps.append(match.group(2))

temps = list(set(temps)) #uniq
temps.sort()
files.sort()

#have temp choosing function
temp = temps[0]

n_f = [x for x in files if re.match("const_complex_na_0_" + temps[-1], x)][0]
s_f = [x for x in files if re.match("const_complex_na_0_" + temp, x)][0]
m_f = [x for x in files if re.match("const_complex_na_[^0]_" + temp, x)]
print m_f

if len(files) == 0:
	print "No data found"
	sys.exit()




s = read_file(base_dir + s_f)
n = read_file(base_dir + n_f)
m = []
for i in m_f:
	temp = read_file(base_dir + i)
	out = transfer_func(s, n, temp)
	m.append(out)

#plot_complexn(s[0], m)
predict, coeffs, errors = fit_nist(s, m)
plot_fit(m, predict)

print coeffs
print errors

#Coeffecients:  [  5.98340389e+00   6.23565491e-01   1.61887268e-01   1.52429108e+04
  # 5.61840916e+00   1.75209612e+01   2.07340763e+01]

