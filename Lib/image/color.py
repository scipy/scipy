
import scipy_base as sb
import scipy

# Various utilities and data for color processing
#  See http://www.cvrl.org/ for much of the data

# RGB  the linear sRGB color space using D65 as a white-point
#      (IEC 61966-2-1).  Represents standard monitor (w/o gamma correction).
# sRGB is the non-linear color-space (w/ gamma correction)
# RGBntsc is the NTSC receiver primary color coordinate system
# RGBcie is the CIE, monochromatic RGB primary system
# RGBsb is the Stiles & Burch (1955) 2-deg color coordinate system

# Primaries for the coordinate systems
cie_primaries = [700, 546.1, 435.8]

sb_primaries = [1./155 * 1e5, 1./190 * 1e5, 1./225 * 1e5]


# Matrices from Jain 

xyz_from_rgbcie = [[0.490, 0.310, 0.200],
                   [0.177, 0.813, 0.011],
                   [0.000, 0.010, 0.990]]

rgbcie_from_xyz = scipy.linalg.inv(XYZ_from_RGBcie)

rgbntsc_from_xyz = [[1.910, -0.533, -0.288],
                    [-0.985, 2.000, -0.028],
                    [0.058, -0.118, 0.896]]

yiq_from_rgbntsc = [[0.299, 0.587, 0.114],
                    [0.596, -0.274, -0.322],
                    [0.211, -0.523, 0.312]]

uvw_from_xyz = [[2.0/3.0, 0, 0],
                [0,1,0],
                [-0.5,1.5,0.5]]

# From sRGB specification

xyz_from_rgb =  [[0.412453, 0.357580, 0.180423],
                 [0.212671, 0.715160, 0.072169],
                 [0.019334, 0.119193, 0.950227]]

rgb_from_xyz = scipy.linalg.inv(XYZ_from_RGB)

# LMS color space spectral matching curves provide the
#  spectral response curves of three types of cones.
# 
# 
# Vos, Estevez, and Walraven (1990)
# with alteration in S-cone sensitivity from
#  Stockman and Sharpe (2000)
# scaled so that sum(LMS) has a peak of 1
#  just like LMS_from_XYZ

lms_from_rgbsb = [[0.14266235473644004, 0.49009667755566039,
                   0.028959576047175539],
                  [0.013676614570405768, 0.35465861798651171,
                   0.029062883056895625],
                  [0.0, 0.00029864360424843419, 0.01837806004659253]]

# use Judd, Vos CIE color matching curves XYZJV
#  with Stockman and Sharpe(2000) S-cone alteration.
# scaled so that sum(LMS,axis=0) has a peak of 1
#  based on CIE standard observer

lms_from_xyz = [[0.15513920309034629, 0.54298741130344153,
                 -0.037010041369525896],
                [-0.15513920309034629, 0.45684891207177714,
                 0.029689739651154123],
                [0.0, 6.3686624249879016e-05, 0.0073203016383768691]]
        
# Read spectral matching curves from file
# XYZJV and RGBsb55 are most modern curves to use
# LMScvrl are the cone response curves from www.cvrl.org
#   (normalized to peak at 1 for all three cones)

varnames = ['xyz31','xyz64','xyzjv','rgbsb55','lmscvrl']
k=-1
thisdict = globals()
for name in ['ciexyz31_1.txt','ciexyz64_1.txt','ciexyzjv.txt',
             'sbrgb2.txt','linss2_10e_1.txt']:
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
        if this[0].strip()[0] not in '0123456789':
            break
        wlen.append(int(this[0].strip()))
        xl.append(float(this[1].strip()))
        yl.append(float(this[2].strip()))
        try:
            zl.append(float(this[3].strip()))
        except ValueError, inst:
            msg = inst.args[0]
            if msg.startswith("empty string"):
                zl.append(0.0)
            else:
                raise inst

    thisdict[varnames[k]] = (wlen,xl,yl,zl)

del thisdict, wlen, xl, yl, zl, afile, lines, this

# XYZ white-point coordinates
#  from http://www.aim-dtp.net/aim/technology/cie_xyz/cie_xyz.htm

whitepoints = {'CIE A': ['Normal incandescent', 0.4476, 0.4074],
               'CIE B': ['Direct sunlight', 0.3457, 0.3585],
               'CIE C': ['Average sunlight', 0.3101, 0.3162],
               'CIE E': ['Normalized reference', 1.0/3, 1.0/3],
               'D50' : ['Bright tungsten', 0.3457, 0.3585],
               'D55' : ['Cloudy daylight', 0.3324, 0.3474],
               'D65' : ['Daylight', 0.312713, 0.329016],
               'D75' : ['?', 0.299, 0.3149],
               'D93' : ['low-quality old CRT', 0.2848, 0.2932]
               }
# convert to X,Y,Z white-point

def triwhite(chrwhite):
    x,y = chrwhite
    X = x / y
    Y = 1.0
    Z = (1-x-y)/y
    return X,Y,Z

for key in whitepoints.keys():
    whitepoints[key].append(triwhite(whitepoints[key][1:]))
del key


def tri2chr(tri,axis=None):
    """Convert tristimulus values to chromoticity values"""
    tri = asarray(tri)
    n = len(tri.shape)
    if axis is None:
        axis = coloraxis(tri.shape)
    slices = []
    for k in range(n):
        slices.append(slice(None))
    slices[axis] = NewAxis
    norm = scipy.sum(tri,axis=axis)[slices]
    slices[axis] = slice(None,2)
    out = tri[slices]/norm
    return out


# find the lowest dimension of size 3
def coloraxis(shape):
    for k, val in enumerate(shape):
        if val == 3:
            return k
    raise ValueError, "No Color axis found."

def convert(matrix,TTT,axis=None):
    TTT = asarray(TTT)
    if axis is None:
        axis = coloraxis(TTT.shape)
    if (axis != 0):
        TTT = sb.swapaxes(TTT,0,axis)
    oldshape = TTT.shape        
    TTT = reshape(TTT,(3,-1))
    OUT = dot(matrix, TTT)
    OUT.shape = oldshape
    if (axis != 0):
        OUT = sb.swapaxes(OUT,axis,0)
    return OUT

def xyz2rgbcie(xyz,axis=None):
    return convert(rgbcie_from_xyz, xyz, axis)

def xyz2rgb(xyz,axis=None):
    return convert(rgb_from_xyz, xyz, axis)

def rgb2xyz(rgb, axis=None):
    return convert(xyz_from_rgb, rgb, axis)

def makeslices(n):
    slices = []
    for k in range(n):
        slices.append(slice(None))
    return slices

def separate_colors(xyz,axis=None):
    if axis is None:
        axis = coloraxis(xyz.shape)
    n = len(xyz.shape)
    slices = makeslices(n)
    slices[axis] = 0
    x = xyz[slices]
    slices[axis] = 1
    y = xyz[slices]
    slices[axis] = 2
    z = xyz[slices]
    return x, y, z, axis
    
def join_colors(c1,c2,c3,axis):
    c1,c2,c3 = asarray(c1),asarray(c2),asarray(c3)
    newshape = c1.shape[:axis] + (1,) + c1.shape[axis:]
    c1.shape = newshape
    c2.shape = newshape
    c3.shape = newshape
    return sb.concatenate((c1,c2,c3),axis=axis)

def xyz2lab(xyz, axis=None, wp=whitepoints['D65'][-1], doclip=True):
    x,y,z,axis = separate_colors(xyz,axis)
    xn,yn,zn = x/wp[0], y/wp[1], z/wp[2]
    def f(t):
        eps = 216/24389.
        kap = 24389/27.        
        return sb.where(t > eps,
                        sb.power(t, 1.0/3),
                        (kap*t + 16.0)/116)
    fx,fy,fz = f(xn), f(yn), f(zn)
    L = 116*fy - 16
    a = 500*(fx - fy)
    b = 200*(fy - fz)
    if doclip:
        L = sb.clip(L, 0.0, 100.0)
        a = sb.clip(a, -500.0, 500.0)
        b = sb.clip(b, -200.0, 200.0)
    return join_colors(L,a,b,axis)

def lab2xyz(lab, axis=None, wp=whitepoints['D65'][-1]):
    lab = asarray(lab)
    L,a,b,axis = separate_colors(lab,axis)
    fy = (L+16)/116.0
    fz = fy - b / 200.
    fx = a/500.0 + fy
    def finv(y):
        eps3 = (216/24389.)**3
        kap = 24389/27.        
        return sb.where(y > eps3,
                        sb.power(y,3),
                        (116*y-16)/kap)
    xr, yr, zr = finv(fx), finv(fy), finv(fz)
    return join_colors(xr*wp[0],yr*wp[1],zr*wp[2],axis)
    
def rgb2lab(rgb):
    return xyz2lab(rgb2xyz(rgb))

def lab2rgb(lab):
    return xyz2rgb(lab2xyz(lab))

# sRGB <-> sR'G'B'  equations from
#   http://www.w3.org/Graphics/Color/sRGB
#   http://www.srgb.com/basicsofsrgb.htm

# rgb2rgbp gives the nonlinear (gamma corrected) sR'G'B' from
#    linear sRGB values
#   approximately the same as rgb**(1.0/2.2)
def rgb2rgbp(rgb):
    rgb = asarray(rgb)
    eps = 0.0031308
    return where(rgb < eps, 12.92*rgb,
                 1.055*rgb**(1.0/2.4) - 0.055)

# rgbp2rgb gives linear sRGB values from nonlinear sR'G'B' values
#
#  approximately the same as rgbp**2.2
def rgbp2rgb(rgbp,axis=None):
    rgbp = asarray(rgbp)
    eps = 0.04045
    return where(rgbp <= eps, rgbp / 12.92,
                 power((rgbp + 0.055)/1.055,2.4))


    
    







