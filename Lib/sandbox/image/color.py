
import numpy as np
from scipy import linalg
import os

# Various utilities and data for color processing
#  See http://www.cvrl.org/ for text file of the data

# rgb  the linear sRGB color space using D65 as a white-point
#      (IEC 61966-2-1).  Represents standard monitor (w/o gamma correction).
# rgbp is the non-linear color-space (w/ gamma correction)
# rgbntsc is the NTSC receiver primary color coordinate system
# rgbcie is the CIE, monochromatic RGB primary system
# rgbsb is the Stiles & Burch (1955) 2-deg color coordinate system

# Primaries for the coordinate systems
cie_primaries = [700, 546.1, 435.8]

sb_primaries = [1./155 * 1e5, 1./190 * 1e5, 1./225 * 1e5]


# Matrices from Jain

xyz_from_rgbcie = [[0.490, 0.310, 0.200],
                   [0.177, 0.813, 0.011],
                   [0.000, 0.010, 0.990]]

rgbcie_from_xyz = linalg.inv(xyz_from_rgbcie)

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

rgb_from_xyz = linalg.inv(xyz_from_rgb)

# From http://www.mir.com/DMG/ycbcr.html

ycbcr_from_rgbp = [[0.299, 0.587, 0.114],
                   [-0.168736, -0.331264, 0.5],
                   [0.5, -0.418688, -0.081312]]

rgbp_from_ycbcr = linalg.inv(ycbcr_from_rgbp)


# LMS color space spectral matching curves provide the
#  spectral response curves of three types of cones.
#
#
# Vos, Estevez, and Walraven (1990)
# with alteration in S-cone sensitivity from
#  Stockman and Sharpe (2000)
# scaled so that sum(LMS,axis=0) has a peak of 1
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
    name = os.path.join(os.path.dirname(__file__),name)
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

del thisdict, wlen, xl, yl, zl, afile, lines, this, line, k, msg, name
del varnames, inst

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
    tri = np.asarray(tri)
    n = len(tri.shape)
    if axis is None:
        axis = coloraxis(tri.shape)
    slices = []
    for k in range(n):
        slices.append(slice(None))
    slices[axis] = np.newaxis
    norm = np.sum(tri,axis=axis)[slices]
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
    TTT = np.asarray(TTT)
    if axis is None:
        axis = coloraxis(TTT.shape)
    if (axis != 0):
        TTT = np.swapaxes(TTT,0,axis)
    oldshape = TTT.shape
    TTT = np.reshape(TTT,(3,-1))
    OUT = np.dot(matrix, TTT)
    OUT.shape = oldshape
    if (axis != 0):
        OUT = np.swapaxes(OUT,axis,0)
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
    c1,c2,c3 = np.asarray(c1),np.asarray(c2),np.asarray(c3)
    newshape = c1.shape[:axis] + (1,) + c1.shape[axis:]
    c1.shape = newshape
    c2.shape = newshape
    c3.shape = newshape
    return np.concatenate((c1,c2,c3),axis=axis)

def xyz2lab(xyz, axis=None, wp=whitepoints['D65'][-1], doclip=1):
    x,y,z,axis = separate_colors(xyz,axis)
    xn,yn,zn = x/wp[0], y/wp[1], z/wp[2]
    def f(t):
        eps = 216/24389.
        kap = 24389/27.
        return np.where(t > eps,
                        np.power(t, 1.0/3),
                        (kap*t + 16.0)/116)
    fx,fy,fz = f(xn), f(yn), f(zn)
    L = 116*fy - 16
    a = 500*(fx - fy)
    b = 200*(fy - fz)
    if doclip:
        L = np.clip(L, 0.0, 100.0)
        a = np.clip(a, -500.0, 500.0)
        b = np.clip(b, -200.0, 200.0)
    return join_colors(L,a,b,axis)

def lab2xyz(lab, axis=None, wp=whitepoints['D65'][-1]):
    lab = np.asarray(lab)
    L,a,b,axis = separate_colors(lab,axis)
    fy = (L+16)/116.0
    fz = fy - b / 200.
    fx = a/500.0 + fy
    def finv(y):
        eps3 = (216/24389.)**3
        kap = 24389/27.
        return np.where(y > eps3,
                        np.power(y,3),
                        (116*y-16)/kap)
    xr, yr, zr = finv(fx), finv(fy), finv(fz)
    return join_colors(xr*wp[0],yr*wp[1],zr*wp[2],axis)

def rgb2lab(rgb):
    return xyz2lab(rgb2xyz(rgb))

def lab2rgb(lab):
    return xyz2rgb(lab2xyz(lab))

def _uv(x, y, z):
    """ The u, v formulae for CIE 1976 L*u*v* computations.
    """
    denominator = (x + 15*y + 3*z)
    zeros = (denominator == 0.0)
    denominator[zeros] = 1.0
    u_numerator = 4 * x
    u_numerator[zeros] = 4.0
    v_numerator = 9 * y
    v_numerator[zeros] = 9.0 / 15.0

    return u_numerator/denominator, v_numerator/denominator

def xyz2luv(xyz, axis=None, wp=whitepoints['D65'][-1]):
    x, y, z, axis = separate_colors(xyz, axis)
    xn, yn, zn = x/wp[0], y/wp[1], z/wp[2]
    Ls = 116.0 * np.power(yn, 1./3) - 16.0
    small_mask = (y <= 0.008856*wp[1])
    Ls[small_mask] = 903.0 * y[small_mask] / wp[1]
    unp, vnp = _uv(*wp)
    up, vp = _uv(x, y, z)
    us = 13 * Ls * (up - unp)
    vs = 13 * Ls * (vp - vnp)

    return join_colors(Ls, us, vs, axis)

def luv2xyz(luv, axis=None, wp=whitepoints['D65'][-1]):
    Ls, us, vs, axis = separate_colors(luv, axis)
    unp, vnp = _uv(*wp)
    small_mask = (Ls <= 903.3 * 0.008856)
    y = wp[1] * ((Ls + 16.0) / 116.0) ** 3
    y[small_mask] = Ls * wp[1] / 903.0
    up = us / (13*Ls) + us
    vp = vs / (13*Ls) + vs
    x = 9.0 * y * up / (4.0 * vp)
    z = -x / 3.0 - 5.0 * y + 3.0 * y/vp

    return join_colors(x, y, z, axis)

def rgb2luv(rgb):
    return xyz2luv(rgb2xyz(rgb))

def luv2rgb(luv):
    return xyz2rgb(luv2xyz(luv))

#  RGB values that will be displayed on a screen are always
#  R'G'B' values.  To get the XYZ value of the color that will be
#  displayed you need a calibrated monitor with a profile
#  -- someday we should support reading and writing such profiles and
#     doing color conversion with them.
#  But, for quick-and-dirty calculation you can often assume the sR'G'B'
#   coordinate system for your computer, and so the rgbp2rgb will
#   put you in the linear coordinate system (assuming normalized to [0,1]
#   sR'G'B' coordiates)
#

# sRGB <-> sR'G'B'  equations from
#   http://www.w3.org/Graphics/Color/sRGB
#   http://www.srgb.com/basicsofsrgb.htm


# Macintosh displays are usually gamma = 1.8

# These transformations are done with normalized [0,1.0] coordinates

# when gamma is None:
#    rgb2rgbp gives the nonlinear (gamma corrected) sR'G'B' from
#       linear sRGB values
#    approximately the same as rgb**(1.0/2.2)
# otherwise do a simple gamma calculation
#    rgbp = rgb**(1.0/gamma)

def rgb2rgbp(rgb,gamma=None):
    rgb = np.asarray(rgb)
    if gamma is None:
        eps = 0.0031308
        return np.where(rgb < eps, 12.92*rgb,
                     1.055*rgb**(1.0/2.4) - 0.055)
    else:
        return rgb**(1.0/gamma)

# when gamma is None:
#     rgbp2rgb gives linear sRGB values from nonlinear sR'G'B' values
#     approximately the same as rgbp**2.2
# otherwise do a simple gamma coorection
#     rgb = rgbp**gamma
#
def rgbp2rgb(rgbp,gamma=None):
    rgbp = np.asarray(rgbp)
    if gamma is None:
        eps = 0.04045
        return np.where(rgbp <= eps, rgbp / 12.92,
                     np.power((rgbp + 0.055)/1.055,2.4))
    else:
        return rgbp**gamma

# The Y'CbCr coordinate system is useful because
#
# Y'CbCr information from here
#  http://www.mir.com/DMG/ycbcr.html
# This transforms from rgbp coordinates to normalized
#  y' cb cr coordinates  y' in [0,1], cb and cr in [-0.5,0.5]
#
# To convert to 8-bit use (according to the web-page cited)
#  Y' = y'*219 + 16   => [16,235]
#  Cb = cb*224 + 128  => [16,240]
#  Cr = cr*224 + 128  => [16,240]

def rgbp2ycbcr(rgbp,axis=None):
    return convert(ycbcr_from_rgbp, rgbp, axis)

def ycbcr2rgbp(ycbcr,axis=None):
    return convert(rgbp_from_ycbcr, ycbcr, axis)

def rgb2ycbcr(rgb,gamma=None,axis=None):
    return rgbp2ycbcr(rgb2rgbp(rgb,gamma),axis)

def ycbcr2rgb(ycbcr,gamma=None,axis=None):
    return rgbp2rgb(ycbcr2rgbp(ycbcr,axis),gamma)

def ycbcr_8bit(ycbcr,axis=None):
    y,cb,cr,axis = separate_colors(ycbcr,axis)
    Y = np.asarray((y*219 + 16),np.uint8)
    Cb = np.asarray((cb*224 + 128),np.uint8)
    Cr = np.asarray((cr*224 + 128),np.uint8)
    return join_colors(Y,Cb,Cr,axis)

def ycbcr_norm(YCbCr,axis=None):
    Y,Cb,Cr,axis = separate_colors(YCbCr,axis)
    y = (Y-16.)/219
    cb = (Cb-128.)/224
    cr = (Cr-128.)/224
    return join_colors(y,cb,cr,axis)
