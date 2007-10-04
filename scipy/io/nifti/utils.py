### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#    Utility function for PyNifti
#
#    Copyright (C) 2007 by
#    Michael Hanke <michael.hanke@gmail.com>
#
#    This is free software; you can redistribute it and/or
#    modify it under the terms of the MIT License.
#
#    This package is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the COPYING
#    file that comes with this package for more details.
#
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##

import numpy
import scipy.io.nifti

def time2vol( t, tr, lag=0.0, decimals=0 ):
    """ Translates a time 't' into a volume number. By default function returns
    the volume number that is closest in time. Volumes are assumed to be
    recorded exactly (and completely) after tr/2, e.g. if 'tr' is 2 secs the
    first volume is recorded at exactly one second.

    't' might be a single value, a sequence or an array.

    The repetition 'tr' might be specified directly, but can also be a 
    NiftiImage object. In the latter case the value of 'tr' is determined from
    the 'rtime' property of the NiftiImage object.

    't' and 'tr' can be given in an arbitrary unit (but both have to be in the
    same unit).

    The 'lag' argument can be used to shift the times by constant offset.

    Please note that numpy.round() is used to round to interger value (rounds
    to even numbers). The 'decimals' argument will be passed to numpy.round().
    """
    # transform to numpy array for easy handling
    tmp = numpy.array(t)
    
    # determine tr if NiftiImage object
    if isinstance( tr, nifti.NiftiImage ):
        tr = tr.rtime

    vol = numpy.round( ( tmp + lag + tr/2 ) / tr, decimals )

    return vol


def applyFxToVolumes( ts, vols, fx, **kwargs ):
    """ Apply a function on selected volumes of a timeseries.

    'ts' is a 4d timeseries. It can be a NiftiImage or a numpy array.
    In case of a numpy array one has to make sure that the time is on the
    first axis. 'ts' can actually be of any dimensionality, but datasets aka
    volumes are assumed to be along the first axis.

    'vols' is either a sequence of sequences or a 2d array indicating which 
    volumes fx should be applied to. Each row defines a set of volumes.

    'fx' is a callable function to get an array of the selected volumes as
    argument. Additonal arguments may be specified as keyword arguments and
    are passed to 'fx'.

    The output will be a 4d array with one computed volume per row in the 'vols'
    array.
    """
    # get data array from nifti image or assume data array is
    # already present
    if isinstance( ts, nifti.NiftiImage ):
        data = ts.data
    else:
        data = ts

    out = []

    for vol in vols:
        out.append( fx( data[ numpy.array( vol ) ], **kwargs ) )

    return numpy.array( out )


def cropImage( nimg, bbox ):
    """ Crop an image.

    'bbox' has to be a sequency of (min,max) tuples (one for each image
    dimension).

    The function returns the cropped image. The data is not shared with the
    original image, but is copied.
    """

    # build crop command
    cmd = 'nimg.data.squeeze()['
    cmd += ','.join( [ ':'.join( [ str(i) for i in dim ] ) for dim in bbox ] )
    cmd += ']'

    # crop the image data array
    cropped = eval(cmd).copy()

    # return the cropped image with preserved header data
    return nifti.NiftiImage(cropped, nimg.header)


def getPeristimulusTimeseries( ts, onsetvols, nvols = 10, fx = numpy.mean ):
    """ Returns 4d array with peristimulus timeseries.

    Parameters:
        ts        - source 4d timeseries
        onsetvols - sequence of onsetvolumes to be averaged over
        nvols     - length of the peristimulus timeseries in volumes
                    (starting from onsetvol)
        fx        - function to be applied to the list of corresponding
                    volumes. Typically this will be mean(), so it is default,
                    but it could also be var() or something different. The
                    supplied function is to be able to handle an 'axis=0'
                    argument similiar to NumPy's mean(), var(), ...
    """
    selected = [ [ o + offset for o in onsetvols ] \
                    for offset in range( nvols ) ]

    if fx == tuple:
        return applyFxToVolumes( ts, selected, fx )
    else:
        return applyFxToVolumes( ts, selected, fx, axis=0 )

