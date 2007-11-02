import struct
import numpy as N

def getFilteredSlice(imageName='junk.raw', bytes=2, rows=512, columns=512):
	image = open(imageName, 'rb')
	slice = image.read(rows*columns*bytes)
	if bytes == 2:
		values = struct.unpack('h'*rows*columns, slice)
	else:
		values = struct.unpack('i'*rows*columns, slice)

	ImageSlice = N.array(values, dtype=float).reshape(rows, columns)
	return (ImageSlice)

def getSlice(mySlice, rows=512, columns=512, bytes=2):
	image = open('C:\PythonStuff\CardiacCT.vol', 'rb')
	image.seek(mySlice*rows*columns*bytes)
	slice = image.read(rows*columns*bytes)
	values = struct.unpack('h'*rows*columns, slice)
	ImageSlice = N.array(values, dtype=float).reshape(rows, columns)
	return (ImageSlice+2048)

def saveSlice(mySlice, filename='junk.raw', rows=512, columns=512, bytes=2):
	# just save the slice to a fixed file
	slice = mySlice.astype(int)
	image = open(filename, 'wb')
	image.write(slice)
	image.close()
	return


def buildDGaussKernel(gWidth, sigma):

	kernel  = zeros((1+2*(gWidth-1)), dtype=float)
	indices = range(1, gWidth)  

	i = 0
	kernel[gWidth-1]  = math.exp(((-i*i)/(2.0 * sigma * sigma)))
	kernel[gWidth-1] *= -(i / (sigma * sigma))
	for i in indices:
		kernel[gWidth-1+i]  = math.exp(((-i*i)/(2.0 * sigma * sigma)))
		kernel[gWidth-1+i] *= -(i / (sigma * sigma))
		kernel[gWidth-1-i]  = -kernel[gWidth-1+i]	# since using correlated1d that doesn't
								# flip the kernel like in convolution

	return kernel


def build2DKernel(aperature, hiFilterCutoff):
	rad = math.pi / 180.0
	HalfFilterTaps = (aperature-1) / 2
	kernel = zeros((aperature), dtype=float)

	LC = 0.0
	HC = hiFilterCutoff * rad 
	t2 = 2.0 * math.pi
	t1 = 2.0 * HalfFilterTaps + 1.0
	indices = range(-HalfFilterTaps, HalfFilterTaps+1, 1)  
	# indices   = arange(-HalfFilterTaps, HalfFilterTaps+1, 1)  
	# indicesNZ = indices[indices.nonzero()]
	# indicesP  = indices[indices>0]

	j = 0
	for i in indices:
	    if i == 0:
		tLOW  = LC
	        tHIGH = HC
	    else:
		tLOW  = math.sin(i*LC)/i
	        tHIGH = math.sin(i*HC)/i
	    t3 = 0.54 + 0.46*(math.cos(i*t2/t1))
	    t4 = t3*(tHIGH-tLOW)
	    kernel[j] = t4
	    j += 1

	# normalize the kernel
	sum = kernel.sum()
	kernel /= sum

	return kernel





