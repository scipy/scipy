
from sigtools import *
import MLab

modedict = {'valid':0, 'same':1, 'full':2}


def convolveND(volume,kernel,mode='full'):
    """out = convolveND(in,kernel{,mode}) returns the true convolution of an N-D
    input (in) with the N-D kernel (with index reversal) and zero-padded edges.
    The size of the return depends on the third argument, mode:

    'valid' (0):  The output consists only of those elements that do not rely
                  on the zero-padding.
    'same'  (1):  The output is the same size as the input centered with respect
                  to the 'full' output.
    'full'  (2):  The output is the full discrete linear convolution with the
                  input.  (Default)
    """
    volume = MLab.asarray(volume)
    kernel = MLab.asarray(kernel)
    if (MLab.product(kernel.shape) > MLab.product(volume.shape)):
        temp = kernel
        kernel = volume
        volume = temp
        del temp

    slice_obj = [slice(None,None,-1)]*len(kernel.shape)
    try:
        val = modedict[mode]
    except KeyError:
        val = mode
        
    return correlateND(volume,kernel[slice_obj],val)



def medfiltND(volume,kernel_size=None):
    """out = medfiltND(input{,kernel_size}) returns a median filtered version of
    input where the median is taken over a window of size kernel_size whose
    elements should be odd (kernel_size defaults to 3 along each axis).
    """
    volume = MLab.asarray(volume)
    if kernel_size == None:
        kernel_size = [3] * len(volume.shape)
    kernel_size = MLab.asarray(kernel_size)

    domain = MLab.ones(kernel_size)

    numels = MLab.product(kernel_size)
    order = numels/2
    if not (numels % 2):    # Even number in window
        return (order_filterND(volume,domain,order-1) + order_filterND(volume,domain,order))
    else:                   # Odd number in window 
        return order_filterND(volume,domain,order)


def wienerND(im,mysize=None,noise=None):
    """out = wienerND(im{,kernel_size,noise_power}) returns a wiener filtered
    version of im with optional kernel size (default is 3 along each axis)
    and noise power given.
    """
    im = asarray(im)
    if mysize == None:
        mysize = [3] * len(im.shape)
    mysize = MLab.asarray(mysize);

    # Estimate the local mean
    lMean = correlateND(im,ones(mysize),1) / MLab.prod(mysize)

    # Estimate the local variance
    lVar = correlateND(im**2,ones(mysize),1) / MLab.prod(mysize) - lMean**2

    # Estimate the noise power if needed.
    if noise==None:
        noise = MLab.mean(ravel(lVar))

    # Compute result
    # f = lMean + (maximum(0, lVar - noise) ./
    #               maximum(lVar, noise)) * (im - lMean) 
    #
    out = im - lMean
    im = lVar - noise
    im = MLab.maximum(im,0)
    lVar = MLab.maximum(lVar,noise)
    out = out / lVar
    out = out * im
    out = out + lMean

    return out



def test():
    a = [3,4,5,6,5,4]
    b = [1,2,3]
    c = convolveND(a,b)
    if (MLab.product(equal(c,[3,10,22,28,32,32,23,12]))==0):
        print "Error in convolveND."

    f = [[3,4,5],[2,3,4],[1,2,5]]
    d = medfiltND(f)
    if (MLab.product(ravel(equal(d,[[0,3,0],[2,3,3],[0,2,0]])))==0):
        print "Error in medfiltND."

    g = MLab.array([[5,6,4,3],[3,5,6,2],[2,3,5,6],[1,6,9,7]],'d')
    correct = MLab.array([[2.16374269,3.2222222222, 2.8888888889, 1.6666666667],[2.666666667, 4.33333333333, 4.44444444444, 2.8888888888],[2.222222222, 4.4444444444, 5.4444444444, 4.801066874837],[1.33333333333, 3.92735042735, 6.0712560386, 5.0404040404]])
    h = wienerND(g)
    if (MLab.abs(MLab.product(MLab.ravel(h-correct)))> 1e-7):
        print "Error in wienerND."

    return

if __name__ == "__main__":
    test()


        
       

    




