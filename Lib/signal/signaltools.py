
import sigtools
import MLab

modedict = {'valid':0, 'same':1, 'full':2}


def correlateND(volume, kernel, mode='full'):
    """ correlateND(in1, in2, mode='full')  Cross-correlation of in1 with in2.

  Description:

     Cross-correlate in1 and in2 with the output size determined by mode.

  Inputs:

    in1 -- an N-dimensional array.
    in2 -- an array with the same number of dimensions as in1.
    mode -- a flag indicating the size of the output
            'valid'  (0): The output consists only of those elements that
                            do not rely on the zero-padding.
            'same'   (1): The output is the same size as the input centered
                            with respect to the 'full' output.
            'full'   (2): The output is the full discrete linear
                            cross-correlation of the inputs. (Default)

  Outputs:  (out,)

    out -- an N-dimensional array containing a subset of the discrete linear
           cross-correlation of in1 with in2.
 
    """
    # Code is faster if kernel is smallest array.
    volume = MLab.asarray(volume)
    kernel = MLab.asarray(kernel)
    if (MLab.product(kernel.shape) > MLab.product(volume.shape)):
        temp = kernel
        kernel = volume
        volume = temp
        del temp

    try:
        val = modedict[mode]
    except KeyError:
        if val not in [0,1,2]:
            raise ValueError, "Acceptable mode flags are 'valid' (0), 'same' (1), or 'full' (2)."
        val = mode

    return sigtools._correlateND(volume, kernel, val)

def convolveND(volume,kernel,mode='full'):
    """ convolveND(in1, in2, mode='full')  Convolution of in1 with in2.

  Description:

     Convolve in1 and in2 with output size determined by mode.

  Inputs:

    in1 -- an N-dimensional array.
    in2 -- an array with the same number of dimensions as in1.
    mode -- a flag indicating the size of the output
            'valid'  (0): The output consists only of those elements that
                            do not rely on the zero-padding.
            'same'   (1): The output is the same size as the input centered
                            with respect to the 'full' output.
            'full'   (2): The output is the full discrete linear convolution
                            of the inputs. (Default)

  Outputs:  (out,)

    out -- an N-dimensional array containing a subset of the discrete linear
           convolution of in1 with in2.

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
        if val not in [0,1,2]:
            raise ValueError, "Acceptable mode flags are 'valid' (0), 'same' (1), or 'full' (2)."
        val = mode
        
    return sigtools._correlateND(volume,kernel[slice_obj],val)

def order_filterND(a, domain, order):
    """
 order_filterND(in, domain, rank)  Perform an order filter on in.

  Description:

    Perform an order filter on the array in.  The domain argument acts as a
    mask centered over each pixel.  The non-zero elements of domain are
    used to select elements surrounding each input pixel which are placed
    in a list.   The list is sorted, and the output for that pixel is the
    element corresponding to rank in the sorted list.
    
  Inputs:

    in -- an N-dimensional input array.
    domain -- a mask array with the same number of dimensions as in.  Each
              dimension should have an odd number of elements.
    rank -- an non-negative integer which selects the element from the sorted
            list (0 corresponds to the largest element, 1 is the next largest
            element, etc.)

  Output: (out,)

    out -- the results of the order filter in an array with the same
           shape as in.
          
    """
    domain = MLab.asarray(domain)
    size = domain.shape
    for k in range(len(size)):
        if (size[k] % 2) != 1:
            raise ValueError, "Each dimension of domain argument should have an odd number of elements."
    return sigtools._orderfilterND(a, domain, rank)
   

def medfiltND(volume,kernel_size=None):
    """
 medfiltND(in, kernel_size=3)  Perform a median filter on input array.

  Description:

    Apply a median filter to the input array using a local window-size
    given by kernel_size.

  Inputs:

    in -- An N-dimensional input array.
    kernel_size -- A scalar or an N-length list giving the size of the
                   median filter window in each dimension.  Elements of
                   kernel_size should be odd.  If kernel_size is a scalar,
                   then this scalar is used as the size in each dimension.

  Outputs: (out,)

    out -- An array the same size as input containing the median filtered
           result.
  
    """
    volume = MLab.asarray(volume)
    if kernel_size == None:
        kernel_size = [3] * len(volume.shape)
    kernel_size = MLab.asarray(kernel_size)
    if len(kernel_size.shape) == 0:
        kernel_size = [kernel_size.toscalar()] * len(volume.shape)
    kernel_size = MLab.asarray(kernel_size)

    for k in range(len(volume.shape)):
        if (kernel_size[k] % 2) != 1:
            raise ValueError, "Each element of kernel_size should be odd." 

    domain = MLab.ones(kernel_size)

    numels = MLab.product(kernel_size)
    order = numels/2
    return sigtools._order_filterND(volume,domain,order)


def wienerND(im,mysize=None,noise=None):
    """
 wienerND(in, kernel_size=3, noise_power=None)  Perform a wiener filter.

  Description:

    Apply a wiener filter to the N-dimensional array in.

  Inputs:

    in -- an N-dimensional array.
    kernel_size -- A scalar or an N-length list giving the size of the
                   median filter window in each dimension.  Elements of
                   kernel_size should be odd.  If kernel_size is a scalar,
                   then this scalar is used as the size in each dimension.
    noise -- The noise-power to use.  If None, then noise is estimated as
             the average of the local variance of the input.

  Outputs: (out,)

    out -- Wiener filtered result with the same shape as in.

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

def convolve2d():
    pass

def correlate2d():
    pass

def medfilt2d():
    pass
    

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


        
       

    




