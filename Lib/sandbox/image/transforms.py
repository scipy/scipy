import scipy.base as sb
import scipy
pi = scipy.pi

# fast discrete cosine transforms of real sequences (using the fft)
#  These implement the DCT-II and inverse DCT-II (DCT-III)
#  described at http://en.wikipedia.org/wiki/Discrete_cosine_transform

def dct(x,axis=-1):
    """Discrete cosine transform based on the FFT.

    For even-length signals it uses an N-point FFT
    For odd-length signals it uses a 2N-point FFT.
    """
    n = len(x.shape)
    N = x.shape[axis]
    even = (N%2 == 0)
    slices = [None]*4
    for k in range(4):
        slices[k] = []
        for j in range(n):
            slices[k].append(slice(None))    
    if even:
        xtilde = 0.0*x
        slices[0][axis] = slice(None,N/2)
        slices[1][axis] = slice(None,None,2)
        slices[2][axis] = slice(N/2,None)
        slices[3][axis] = slice(N,None,-2)
    else:
        newshape = list(x.shape)
        newshape[axis] = 2*N
        xtilde = sb.empty(newshape,sb.Float)
        slices[0][axis] = slice(None,N)
        slices[2][axis] = slice(N,None)
        slices[3][axis] = slice(None,None,-1)

    for k in range(4):
        slices[k] = tuple(slices[k])
    xtilde[slices[0]] = x[slices[1]]
    xtilde[slices[2]] = x[slices[3]]
    Xt = scipy.fft(xtilde,axis=axis)
    pk = sb.exp(-1j*pi*sb.arange(N)/(2*N))
    newshape = sb.ones(n)
    newshape[axis] = N
    pk.shape = newshape

    if not even:
        pk /= 2;
        Xt = Xt[slices[0]]

    return sb.real(Xt*pk)

    
def idct(v,axis=-1):
    n = len(v.shape)
    N = v.shape[axis]
    even = (N%2 == 0)
    slices = [None]*4
    for k in range(4):
        slices[k] = []
        for j in range(n):
            slices[k].append(slice(None))    
    k = arange(N)
    if even:
        ak = sb.r_[1.0,[2]*(N-1)]*exp(1j*pi*k/(2*N))
        newshape = ones(n)
        newshape[axis] = N
        ak.shape = newshape
        xhat = real(scipy.ifft(v*ak,axis=axis))
        x = 0.0*v
        slices[0][axis] = slice(None,None,2)
        slices[1][axis] = slice(None,N/2)
        slices[2][axis] = slice(N,None,-2)
        slices[3][axis] = slice(N/2,None) 
        for k in range(4):
            slices[k] = tuple(slices[k])
        x[slices[0]] = xhat[slices[1]]
        x[slices[2]] = xhat[slices[3]]
        return x
    else:
        ak = 2*sb.exp(1j*pi*k/(2*N))
        newshape = ones(n)
        newshape[axis] = N
        ak.shape = newshape
        newshape = list(v.shape)
        newshape[axis] = 2*N
        Y = zeros(newshape,sb.Complex)
        #Y[:N] = ak*v
        #Y[(N+1):] = conj(Y[N:0:-1])
        slices[0][axis] = slice(None,N)
        slices[1][axis] = slice(None,None)
        slices[2][axis] = slice(N+1,None)
        slices[3][axis] = slice((N-1),0,-1)
        Y[slices[0]] = ak*v
        Y[slices[2]] = conj(Y[slices[3]])
        x = real(scipy.ifft(Y,axis=axis))[slices[0]]
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

def makeS(N):
    n,k = ogrid[:N,:N]
    C = sin(pi*(k+1)*(n+1)/(N+1))
    return C

# DST-I 
def dst(x,axis=-1):
    """Discrete Sine Transform (DST-I)

    Implemented using 2(N+1)-point FFT
    xsym = r_[0,x,0,-x[::-1]]
    DST = (-imag(fft(xsym))/2)[1:(N+1)]

    adjusted to work over an arbitrary axis for entire n-dim array
    """
    n = len(x.shape)
    N = x.shape[axis]
    slices = [None]*3
    for k in range(3):
        slices[k] = []
        for j in range(n):
            slices[k].append(slice(None))    
    newshape = list(x.shape)
    newshape[axis] = 2*(N+1)
    xtilde = sb.zeros(newshape,sb.Float)
    slices[0][axis] = slice(1,N+1)
    slices[1][axis] = slice(N+2,None)
    slices[2][axis] = slice(None,None,-1)
    for k in range(3):
        slices[k] = tuple(slices[k])
    xtilde[slices[0]] = x
    xtilde[slices[1]] = -x[slices[2]]
    Xt = scipy.fft(xtilde,axis=axis)
    return (-sb.imag(Xt)/2)[slices[0]]

def idst(v,axis=-1):
    n = len(v.shape)
    N = v.shape[axis]
    slices = [None]*3
    for k in range(3):
        slices[k] = []
        for j in range(n):
            slices[k].append(slice(None))    
    newshape = list(v.shape)
    newshape[axis] = 2*(N+1)
    Xt = sb.zeros(newshape,sb.Complex)
    slices[0][axis] = slice(1,N+1)
    slices[1][axis] = slice(N+2,None)
    slices[2][axis] = slice(None,None,-1)
    val = 2j*v
    for k in range(3):
        slices[k] = tuple(slices[k])
    Xt[slices[0]] = -val
    Xt[slices[1]] = val[slices[2]]
    xhat = real(scipy.ifft(Xt,axis=axis))
    return xhat[slices[0]]

def dst2(x,axes=(-1,-2)):
    return dst(dst(x,axis=axes[0]),axis=axes[1])

def idst2(v,axes=(-1,-2)):
    return idst(idst(v,axis=axes[0]),axis=axes[1])

def dstn(x,axes=None):
    if axes is None:
        axes = arange(len(x.shape))
    res = x
    for k in axes:
        res = dst(res,axis=k)
    return res

def idstn(v,axes=None):
    if axes is None:
        axes = arange(len(v.shape))
    res = v
    for k in axes:
        res = idst(res,axis=k)
    return res

def digitrevorder(x,base):
    x = asarray(x)
    rem = N = len(x)
    L = 0
    while 1:
        if rem < base:
            break
        intd = rem // base
        if base*intd != rem:
            raise ValueError, "Length of data must be power of base."
        rem = intd
        L += 1
    vec = r_[[base**n for n in range(L)]]
    newx = x[NewAxis,:]*vec[:,NewAxis]
    # compute digits 
    for k in range(L-1,-1,-1):
        newx[k] = x // vec[k]
        x = x - newx[k]*vec[k]
    # reverse digits
    newx = newx[::-1,:]
    x = 0*x
    # construct new result from reversed digts
    for k in range(L):
        x += newx[k]*vec[k]
    return x
               

def bitrevorder(x):
    return digitrevorder(x,2)
    

# needs to be fixed
def wht(data):
    """Walsh-Hadamaard Transform (sequency ordered)

    adapted from MATLAB algorithm published on the web by
    Author: Gylson Thomas
    e-mail: gylson_thomas@yahoo.com
    Asst. Professor, Electrical and Electronics Engineering Dept.
    MES College of Engineering Kuttippuram,
    Kerala, India, February 2005.
    copyright 2005.
    Reference: N.Ahmed, K.R. Rao, "Orthogonal Transformations for 
    Digital Signal Processing" Spring Verlag, New York 1975. page-111.
    """
    N = len(data)
    L=log2(N);
    if ((L-floor(L)) > 0.0):
        raise ValueError, "Length must be power of 2"
    x=bitrevorder(data);

    k1=N; k2=1; k3=N/2;
    for i1 in range(1,L+1):  #Iteration stage 
        L1=1;
        for i2 in range(1,k2+1):
            for i3 in range(1,k3+1):
                i=i3+L1-1; j=i+k3;
                temp1= x[i-1]; temp2 = x[j-1]; 
                if (i2 % 2) == 0:
                  x[i-1] = temp1 - temp2;
                  x[j-1] = temp1 + temp2;
                  x[i-1] = temp1 + temp2;
                  x[j-1] = temp1 - temp2;
                L1=L1+k1;
            k1 = k1/2;  k2 = k2*2;  k3 = k3/2;
    x = x*1.0/N; # Delete this line for inverse wht
    return x
    
