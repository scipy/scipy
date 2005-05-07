
# Various utilities and data for color processing

# RGB is the NTSC receiver primary color coordinate system
# CIERGB is the CIE, spectral RGB primary system 

CIE_P1 = 700
CIE_P2 = 546.1
CIE_P3 = 435.8

XYX_from_CIERGB = [[0.490, 0.310, 0.200],
                   [0.177, 0.813, 0.011],
                   [0.000, 0.010, 0.990]]

RGB_from_XYZ = [[1.910, -0.533, -0.288],
                [-0.985, 2.000, -0.028],
                [0.058, -0.118, 0.896]]

YIQ_from_RGB = [[0.299, 0.587, 0.114],
                [0.596, -0.274, -0.322],
                [0.211, -0.523, 0.312]]

UVW_from_XYZ = [[2.0/3.0, 0, 0],
                [0,1,0],
                [-0.5,1.5,0.5]]

# Read XYZ spectral matching curves from file

varnames = ['XYZ_CIE31','XYZ_CIE64','XYZ_CIEJV']
k=-1
thisdict = globals()
for name in ['ciexyz31_1.txt','ciexyz64_1.txt','ciexyzjv.txt']:
    k = k + 1
    afile = open(name)
    lines = afile.readlines()
    afile.close()
    wlen = []
    xl = []
    yl = []
    zl = []
    for line in lines:
        this = line.split(',')
        if this[0][0] not in '0123456789':
            break
        wlen.append(int(this[0].strip()))
        xl.append(float(this[1].strip()))
        yl.append(float(this[2].strip()))
        zl.append(float(this[3].strip()))

    thisdict[varnames[k]] = (wlen,xl,yl,zl)

del thisdict, wlen, xl, yl, zl, afile, lines, this





