import scipy_base as sb
import scipy
pi = scipy.pi

# fast discrete cosine transforms of real sequences (using the fft)
#  These implement the DCT-II and inverse DCT-II (DCT-III)
#  described at http://en.wikipedia.org/wiki/Discrete_cosine_transform

def dct(x,axis=-1):
    n = len(x.shape)
    N = x.shape[axis]
    if (N%2 != 0):
        raise ValueError, "Length of sequence must be even."
    xtilde = 0.0*x
    slices = [None]*4
    for k in range(4):
        slices[k] = []
        for j in range(n):
            slices[k].append(slice(None))
    slices[0][axis] = slice(None,N/2)
    slices[1][axis] = slice(None,None,2)
    slices[2][axis] = slice(N/2,None)
    slices[3][axis] = slice(N,None,-2)
    for k in range(4):
        slices[k] = tuple(slices[k])
    xtilde[slices[0]] = x[slices[1]]
    xtilde[slices[2]] = x[slices[3]]
    Xt = scipy.fft(xtilde,axis=axis)
    pk = exp(-1j*pi*arange(N)/(2*N))
    newshape = ones(n)
    newshape[axis] = N
    pk.shape = newshape
    return squeeze(real(Xt*pk))
    
def idct(v,axis=-1):
    n = len(v.shape)
    N = v.shape[axis]    
    if (N%2 != 0):
        raise ValueError, "Length of sequence must be even."
    k = arange(N)
    ak = sb.r_[1.0,[2]*(N-1)]*exp(1j*pi*k/(2*N))
    newshape = ones(n)
    newshape[axis] = N
    ak.shape = newshape
    xhat = real(scipy.ifft(v*ak,axis=axis))
    x = 0.0*v
    slices = [None]*4
    for k in range(4):
        slices[k] = []
        for j in range(n):
            slices[k].append(slice(None))
    slices[0][axis] = slice(None,None,2)
    slices[1][axis] = slice(None,N/2)
    slices[2][axis] = slice(N,None,-2)
    slices[3][axis] = slice(N/2,None) 
    for k in range(4):
        slices[k] = tuple(slices[k])
    x[slices[0]] = xhat[slices[1]]
    x[slices[2]] = xhat[slices[3]]
    return x

def dct2(x,axes=(-1,-2)):
    return dct(dct(x,axis=axes[0]),axis=axes[1])

def idct2(v,axes=(-1,-2)):
    return idct(idct(v,axis=axes[0]),axis=axes[1])

def dctn(x,axes=None):
    if axes is None:
        axes = arange(len(x.shape))
    res = x
    for k in axes:
        res = dct(res,axis=k)
    return res

def idctn(v,axes=None):
    if axes is None:
        axes = arange(len(v.shape))
    res = v
    for k in axes:
        res = idct(res,axis=k)
    return res
    
def makeC(N):
    n,l = ogrid[:N,:N]
    C = cos(pi*(2*n+1)*l/(2*N))
    return C

def dct2raw(x):
    M,N = x.shape
    CM = makeC(M)
    CN = makeC(N)
    return dot(transpose(CM),dot(x,CN))

def idct2raw(v):
    M,N = v.shape
    iCM = linalg.inv(makeC(M))
    iCN = linalg.inv(makeC(N))
    return dot(transpose(iCM),dot(v,iCN))
