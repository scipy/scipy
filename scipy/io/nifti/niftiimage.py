### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ##
#
#    Python interface to the NIfTI file format
#
#    Copyright (C) 2006-2007 by
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

# the swig wrapper if the NIfTI C library
import numpy
from scipy.io.nifti import nifticlib

class NiftiImage(object):
    """Wrapper class for convenient access to NIfTI data.

    The class can either load an image from file or convert a NumPy
    array into a NIfTI file structure. Either way is automatically determined
    by the type of the 'source' argument (string == filename, array == Numpy).

    One can optionally specify whether the image data should be loaded into
    memory when opening NIfTI data from files ('load'). When converting a NumPy
    array one can optionally specify a dictionary with NIfTI header data as
    available via the 'header' attribute.
    """

    filetypes = [ 'ANALYZE', 'NIFTI', 'NIFTI_PAIR', 'ANALYZE_GZ', 'NIFTI_GZ',
                  'NIFTI_PAIR_GZ' ]

    # map NumPy datatypes to NIfTI datatypes
    numpy2nifti_dtype_map = { numpy.uint8: nifticlib.NIFTI_TYPE_UINT8,
                              numpy.int8 : nifticlib.NIFTI_TYPE_INT8,
                              numpy.uint16: nifticlib.NIFTI_TYPE_UINT16,
                              numpy.int16 : nifticlib.NIFTI_TYPE_INT16,
                              numpy.uint32: nifticlib.NIFTI_TYPE_UINT32,
                              numpy.int32 : nifticlib.NIFTI_TYPE_INT32,
                              numpy.uint64: nifticlib.NIFTI_TYPE_UINT64,
                              numpy.int64 : nifticlib.NIFTI_TYPE_INT64,
                              numpy.float32: nifticlib.NIFTI_TYPE_FLOAT32,
                              numpy.float64: nifticlib.NIFTI_TYPE_FLOAT64,
                              numpy.complex128: nifticlib.NIFTI_TYPE_COMPLEX128
                            }


    def numpydtype2niftidtype(array):
        """ Return the NIfTI datatype id for a corrsponding numpy array
        datatype.
        """
        # get the real datatype from numpy type dictionary
        dtype = numpy.typeDict[str(array.dtype)]

        if not NiftiImage.numpy2nifti_dtype_map.has_key(dtype):
            raise ValueError, "Unsupported datatype '%s'" % str(array.dtype)

        return NiftiImage.numpy2nifti_dtype_map[dtype]
    numpydtype2niftidtype = staticmethod(numpydtype2niftidtype)

    def splitFilename(filename):
        """ Split a NIfTI filename and returns a tuple of basename and
        extension. If no valid NIfTI filename extension is found, the whole
        string is returned as basename and the extension string will be empty.
        """

        parts = filename.split('.')

        if parts[-1] == 'gz':
            if not parts[-2] in [ 'nii', 'hdr', 'img' ]:
                return filename, ''
            else:
                return '.'.join(parts[:-2]), '.'.join(parts[-2:])
        else:
            if not parts[-1] in [ 'nii', 'hdr', 'img' ]:
                return filename, ''
            else:
                return '.'.join(parts[:-1]), parts[-1]
    splitFilename = staticmethod(splitFilename)

    def nhdr2dict(nhdr):
        """ Convert a NIfTI header struct into a python dictionary.

        While most elements of the header struct will be translated
        1:1 some (e.g. sform matrix) are converted into more convenient
        datatypes (i.e. 4x4 matrix instead of 16 separate values).
        """
        h = {}

        # the following header elements are converted in a simple loop
        # as they do not need special handling
        auto_convert = [ 'session_error', 'extents', 'sizeof_hdr',
                         'slice_duration', 'slice_start', 'xyzt_units',
                         'cal_max', 'intent_p1', 'intent_p2', 'intent_p3',
                         'intent_code', 'sform_code', 'cal_min', 'scl_slope',
                         'slice_code', 'bitpix', 'descrip', 'glmin', 'dim_info',
                         'glmax', 'data_type', 'aux_file', 'intent_name',
                         'vox_offset', 'db_name', 'scl_inter', 'magic',
                         'datatype', 'regular', 'slice_end', 'qform_code',
                         'toffset' ]


        # now just dump all attributes into a dict
        for attr in auto_convert:
            h[attr] = eval('nhdr.' + attr)

        # handle a few special cases
        # handle 'pixdim': carray -> list
        pixdim = nifticlib.floatArray_frompointer(nhdr.pixdim)
        h['pixdim'] = [ pixdim[i] for i in range(8) ]

        # handle dim: carray -> list
        dim = nifticlib.shortArray_frompointer(nhdr.dim)
        h['dim'] = [ dim[i] for i in range(8) ]

        # handle sform: carrays -> (4x4) numpy array
        srow_x = nifticlib.floatArray_frompointer( nhdr.srow_x )
        srow_y = nifticlib.floatArray_frompointer( nhdr.srow_y )
        srow_z = nifticlib.floatArray_frompointer( nhdr.srow_z )

        h['sform'] = numpy.array( [ [ srow_x[i] for i in range(4) ],
                                    [ srow_y[i] for i in range(4) ],
                                    [ srow_z[i] for i in range(4) ],
                                    [ 0.0, 0.0, 0.0, 1.0 ] ] )

        # handle qform stuff: 3 numbers -> list
        h['quatern'] = [ nhdr.quatern_b, nhdr.quatern_c, nhdr.quatern_d ]
        h['qoffset'] = [ nhdr.qoffset_x, nhdr.qoffset_y, nhdr.qoffset_z ]

        return h
    nhdr2dict = staticmethod(nhdr2dict)

    def updateNiftiHeaderFromDict(nhdr, hdrdict):
        """ Update a NIfTI header struct with data from a dictionary.

        The supplied dictionary might contain additonal data elements
        that do not match any nifti header element. These are silently ignored.

        Several checks are performed to ensure validity of the resulting
        nifti header struct. If any check fails a ValueError exception will be
        thrown. However, some tests are still missing.
        """
        # this function is still incomplete. add more checks

        if hdrdict.has_key('data_type'):
            if len(hdrdict['data_type']) > 9:
                raise ValueError, \
                      "Nifti header property 'data_type' must not be longer " \
                      + "than 9 characters."
            nhdr.data_type = hdrdict['data_type']
        if hdrdict.has_key('db_name'):
            if len(hdrdict['db_name']) > 79:
                raise ValueError, "Nifti header property 'db_name' must " \
                                  + "not be longer than 17 characters."
            nhdr.db_name = hdrdict['db_name']

        if hdrdict.has_key('extents'):
            nhdr.extents = hdrdict['extents']
        if hdrdict.has_key('session_error'):
            nhdr.session_error = hdrdict['session_error']

        if hdrdict.has_key('regular'):
            if len(hdrdict['regular']) > 1:
                raise ValueError, \
                      "Nifti header property 'regular' has to be a single " \
                      + "character."
            nhdr.regular = hdrdict['regular']
        if hdrdict.has_key('dim_info'):
            if len(hdrdict['dim_info']) > 1:
                raise ValueError, \
                      "Nifti header property 'dim_info' has to be a " \
                      + "single character."
            nhdr.dim_info = hdrdict['dim_info']

        if hdrdict.has_key('dim'):
            dim = nifticlib.shortArray_frompointer(nhdr.dim)
            for i in range(8): dim[i] = hdrdict['dim'][i]
        if hdrdict.has_key('intent_p1'):
            nhdr.intent_p1 = hdrdict['intent_p1']
        if hdrdict.has_key('intent_p2'):
            nhdr.intent_p2 = hdrdict['intent_p2']
        if hdrdict.has_key('intent_p3'):
            nhdr.intent_p3 = hdrdict['intent_p3']
        if hdrdict.has_key('intent_code'):
            nhdr.intent_code = hdrdict['intent_code']
        if hdrdict.has_key('datatype'):
            nhdr.datatype = hdrdict['datatype']
        if hdrdict.has_key('bitpix'):
            nhdr.bitpix = hdrdict['bitpix']
        if hdrdict.has_key('slice_start'):
            nhdr.slice_start = hdrdict['slice_start']
        if hdrdict.has_key('pixdim'):
            pixdim = nifticlib.floatArray_frompointer(nhdr.pixdim)
            for i in range(8): pixdim[i] = hdrdict['pixdim'][i]
        if hdrdict.has_key('vox_offset'):
            nhdr.vox_offset = hdrdict['vox_offset']
        if hdrdict.has_key('scl_slope'):
            nhdr.scl_slope = hdrdict['scl_slope']
        if hdrdict.has_key('scl_inter'):
            nhdr.scl_inter = hdrdict['scl_inter']
        if hdrdict.has_key('slice_end'):
            nhdr.slice_end = hdrdict['slice_end']
        if hdrdict.has_key('slice_code'):
            nhdr.slice_code = hdrdict['slice_code']
        if hdrdict.has_key('xyzt_units'):
            nhdr.xyzt_units = hdrdict['xyzt_units']
        if hdrdict.has_key('cal_max'):
            nhdr.cal_max = hdrdict['cal_max']
        if hdrdict.has_key('cal_min'):
            nhdr.cal_min = hdrdict['cal_min']
        if hdrdict.has_key('slice_duration'):
            nhdr.slice_duration = hdrdict['slice_duration']
        if hdrdict.has_key('toffset'):
            nhdr.toffset = hdrdict['toffset']
        if hdrdict.has_key('glmax'):
            nhdr.glmax = hdrdict['glmax']
        if hdrdict.has_key('glmin'):
            nhdr.glmin = hdrdict['glmin']

        if hdrdict.has_key('descrip'):
            if len(hdrdict['descrip']) > 79:
                raise ValueError, \
                      "Nifti header property 'descrip' must not be longer " \
                      + "than 79 characters."
            nhdr.descrip = hdrdict['descrip']
        if hdrdict.has_key('aux_file'):
            if len(hdrdict['aux_file']) > 23:
                raise ValueError, \
                      "Nifti header property 'aux_file' must not be longer " \
                      + "than 23 characters."
            nhdr.aux_file = hdrdict['aux_file']

        if hdrdict.has_key('qform_code'):
            nhdr.qform_code = hdrdict['qform_code']

        if hdrdict.has_key('sform_code'):
            nhdr.sform_code = hdrdict['sform_code']

        if hdrdict.has_key('quatern'):
            if not len(hdrdict['quatern']) == 3:
                raise ValueError, \
                      "Nifti header property 'quatern' must be float 3-tuple."

            nhdr.quatern_b = hdrdict['quatern'][0]
            nhdr.quatern_c = hdrdict['quatern'][1]
            nhdr.quatern_d = hdrdict['quatern'][2]

        if hdrdict.has_key('qoffset'):
            if not len(hdrdict['qoffset']) == 3:
                raise ValueError, \
                      "Nifti header property 'qoffset' must be float 3-tuple."

            nhdr.qoffset_x = hdrdict['qoffset'][0]
            nhdr.qoffset_y = hdrdict['qoffset'][1]
            nhdr.qoffset_z = hdrdict['qoffset'][2]

        if hdrdict.has_key('sform'):
            if not hdrdict['sform'].shape == (4,4):
                raise ValueError, \
                      "Nifti header property 'sform' must be 4x4 matrix."

            srow_x = nifticlib.floatArray_frompointer(nhdr.srow_x)
            for i in range(4): srow_x[i] = hdrdict['sform'][0][i]
            srow_y = nifticlib.floatArray_frompointer(nhdr.srow_y)
            for i in range(4): srow_y[i] = hdrdict['sform'][1][i]
            srow_z = nifticlib.floatArray_frompointer(nhdr.srow_z)
            for i in range(4): srow_z[i] = hdrdict['sform'][2][i]

        if hdrdict.has_key('intent_name'):
            if len(hdrdict['intent_name']) > 15:
                raise ValueError, \
                      "Nifti header property 'intent_name' must not be " \
                      + "longer than 15 characters."
            nhdr.intent_name = hdrdict['intent_name']

        if hdrdict.has_key('magic'):
            if hdrdict['magic'] != 'ni1' and hdrdict['magic'] != 'n+1':
                raise ValueError, \
                      "Nifti header property 'magic' must be 'ni1' or 'n+1'."
            nhdr.magic = hdrdict['magic']
    updateNiftiHeaderFromDict = staticmethod(updateNiftiHeaderFromDict)

    def __init__(self, source, header = {}, load=False ):
        """ Create a Niftifile object.

        This method decides whether to load a nifti image from file or create
        one from array data, depending on the datatype of 'source'. If source
        is a string, it is assumed to be a filename and an attempt will be made
        to open the corresponding NIfTI file. If 'load' is set to True the image
        data will be loaded into memory.

        If 'source' is a numpy array the array data will be used for the to be
        created nifti image and a matching nifti header is generated. Additonal
        header data might be supplied in a dictionary. However, dimensionality
        and datatype are determined from the numpy array and not taken from
        a header dictionary.

        If an object of a different type is supplied as 'source' a ValueError
        exception will be thrown.
        """

        self.__nimg = None

        if type( source ) == numpy.ndarray:
            self.__newFromArray( source, header )
        elif type ( source ) == str:
            self.__newFromFile( source, load )
        else:
            raise ValueError, \
                  "Unsupported source type. Only NumPy arrays and filename " \
                  + "string are supported."


    def __del__(self):
        """ Do all necessary cleanups by calling __close().
        """
        self.__close()


    def __close(self):
        """Close the file and free all unnecessary memory.
        """
        if self.__nimg:
            nifticlib.nifti_image_free(self.__nimg)
            self.__nimg = None


    def __newFromArray(self, data, hdr = {}):
        """ Create a nifti image struct from a numpy array and optional header
        data.
        """

        # check array
        if len(data.shape) > 7:
            raise ValueError, \
                  "NIfTI does not support data with more than 7 dimensions."

        # create template nifti header struct
        niptr = nifticlib.nifti_simple_init_nim()
        nhdr = nifticlib.nifti_convert_nim2nhdr(niptr)

        # intermediate cleanup
        nifticlib.nifti_image_free(niptr)

        # convert virgin nifti header to dict to merge properties
        # with supplied information and array properties
        hdic = NiftiImage.nhdr2dict(nhdr)

        # copy data from supplied header dict
        for k, v in hdr.iteritems():
            hdic[k] = v

        # finally set header data that is determined by the data array
        # convert numpy to nifti datatype
        hdic['datatype'] = self.numpydtype2niftidtype(data)

        # make sure there are no zeros in the dim vector
        # especially not in #4 as FSLView doesn't like that
        hdic['dim'] = [ 1 for i in hdic['dim'] ]

        # set number of dims
        hdic['dim'][0] = len(data.shape)

        # set size of each dim (and reverse the order to match nifti format
        # requirements)
        for i, s in enumerate(data.shape):
            hdic['dim'][len(data.shape)-i] = s

        # set magic field to mark as nifti file
        hdic['magic'] = 'n+1'

        # update nifti header with information from dict
        NiftiImage.updateNiftiHeaderFromDict(nhdr, hdic)

        # make clean table
        self.__close()

        # convert nifti header to nifti image struct
        self.__nimg = nifticlib.nifti_convert_nhdr2nim(nhdr, 'pynifti_none')

        if not self.__nimg:
            raise RuntimeError, "Could not create nifti image structure."

        # kill filename for nifti images from arrays
        self.__nimg.fname = ''
        self.__nimg.iname = ''

        # allocate memory for image data
        if not nifticlib.allocateImageMemory(self.__nimg):
            raise RuntimeError, "Could not allocate memory for image data."

        # assign data
        self.data[:] = data[:]


    def __newFromFile(self, filename, load=False):
        """Open a NIfTI file.

        If there is already an open file it is closed first. If 'load' is True
        the image data is loaded into memory.
        """
        self.__close()
        self.__nimg = nifticlib.nifti_image_read( filename, int(load) )

        if not self.__nimg:
            raise RuntimeError, "Error while opening nifti header."

        if load:
            self.load()


    def save(self, filename=None, filetype = 'NIFTI'):
        """Save the image.

        If the image was created using array data (not loaded from a file) one
        has to specify a filename.

        Setting the filename also determines the filetype (NIfTI/ANALYZE).
        Please see the documentation of the setFilename() method for some 
        details on the 'filename' and 'filetype' argument.

        Calling save() without a specified filename on a NiftiImage loaded
        from a file, will overwrite the original file.

        If not yet done already, the image data will be loaded into memory
        before saving the file.

        Warning: There will be no exception if writing fails for any reason,
        as the underlying function nifti_write_hdr_img() from libniftiio does
        not provide any feedback. Suggestions for improvements are appreciated.
        """

        # If image data is not yet loaded, do it now.
        # It is important to do it already here, because nifti_image_load
        # depends on the correct filename set in the nifti_image struct
        # and this will be modified in this function!
        if not self.__haveImageData():
            self.load()

        # set a default description if there is none
        if not self.description:
            self.description = 'Created with PyNIfTI'

        # update header information
        self.updateCalMinMax()

        # saving for the first time?
        if not self.filename or filename:
            if not filename:
                raise ValueError, \
                      "When saving an image for the first time a filename " \
                      + "has to be specified."

            self.setFilename(filename, filetype)

        # now save it
        nifticlib.nifti_image_write_hdr_img(self.__nimg, 1, 'wb')
        # yoh comment: unfortunately return value of nifti_image_write_hdr_img
        # can't be used to track the successful completion of save
        # raise IOError, 'An error occured while attempting to save the image
        # file.'


    def __haveImageData(self):
        """Returns true if the image data was loaded into memory.
        or False if not.

        See: load(), unload()
        """
        self.__ensureNiftiImage()

        if self.__nimg.data:
            return True
        else:
            return False


    def load(self):
        """Load the image data into memory.

        It is save to call this method several times.
        """
        self.__ensureNiftiImage()

        if nifticlib.nifti_image_load( self.__nimg ) < 0:
            raise RuntimeError, "Unable to load image data."


    def unload(self):
        """Unload image data and free allocated memory.
        """
        # if no filename is se, the data will be lost and cannot be recovered
        if not self.filename:
            raise RuntimeError, "No filename is set, unloading the data would " \
                              + "loose it completely without a chance of recovery. "
        self.__ensureNiftiImage()

        nifticlib.nifti_image_unload(self.__nimg)


    def getDataArray(self):
        """ Calls asarray(False) to return the NIfTI image data wrapped into
        a NumPy array.

        Attention: The array shares the data with the NiftiImage object. Any
        resize operation or datatype conversion will most likely result in a
        fatal error. If you need to perform such things, get a copy
        of the image data by using asarray(copy=True).

        The 'data' property is an alternative way to access this function.
        """
        return self.asarray(False)


    def asarray(self, copy = True):
        """Convert the image data into a multidimensional array.

        Attention: If copy == False (the default) the array only wraps
        the image data. Any modification done to the array is also done
        to the image data.

        If copy is true the array contains a copy of the image data.

        Changing the shape, size or data of a wrapping array is not supported
        and will most likely result in a fatal error. If you want to data
        anything else to the data but reading or simple value assignment
        use a copy of the data by setting the copy flag. Later you can convert
        the modified data array into a NIfTi file again.
        """
        self.__ensureNiftiImage()

        if not self.__haveImageData():
            self.load()

        a = nifticlib.wrapImageDataWithArray(self.__nimg)

        if copy:
            return a.copy()
        else:
            return a


    def getScaledData(self):
        """ Returns a copy of the data array scaled by multiplying with the
        slope and adding the intercept that is stored in the NIfTI header.
        """
        data = self.asarray(copy = True)

        return data * self.slope + self.intercept


    def __ensureNiftiImage(self):
        """Check whether a NIfTI image is present.

        Returns True if there is a nifti image file structure or False
        otherwise. One can create a file structure by calling open().
        """
        if not self.__nimg:
            raise RuntimeError, "There is no NIfTI image file structure."


    def updateCalMinMax(self):
        """ Update the image data maximum and minimum value in the
        nifti header.
        """
        self.__nimg.cal_max = float(self.data.max())
        self.__nimg.cal_min = float(self.data.min())


    def getVoxDims(self):
        """ Returns a 3-tuple a voxel dimensions/size in (x,y,z).

        The 'voxdim' property is an alternative way to access this function.
        """
        return ( self.__nimg.dx, self.__nimg.dy, self.__nimg.dz )


    def setVoxDims(self, value):
        """ Set voxel dimensions/size.

        This method takes a 3-tuple of floats as argument. The qform matrix
        and its inverse will be recalculated automatically.

        Besides reading it is also possible to set the voxel dimensions by
        assigning to the 'voxdim' property.
        """
        if len(value) != 3:
            raise ValueError, 'Requires 3-tuple.'

        self.__nimg.dx = float(value[0])
        self.__nimg.dy = float(value[1])
        self.__nimg.dz = float(value[2])

        self.updateQFormFromQuaternion()


    def setPixDims(self, value):
        """ Set the pixel dimensions.

        The methods takes a sequence of up to 7 values (max. number of
        dimensions supported by the NIfTI format.

        The supplied sequence can be shorter than seven elements. In this case
        only present values are assigned starting with the first dimension
        (spatial: x). Calling setPixDims() with a length-3 sequence equals
        calling setVoxDims().
        """
        if len(value) > 7:
            raise ValueError, \
                  'The Nifti format does not support more than 7 dimensions.'

        pixdim = nifticlib.floatArray_frompointer( self.__nimg.pixdim )

        for i in value:
            pixdim[i+1] = float(i)


    def getPixDims(self):
        """ Returns the pixel dimensions on all 7 dimensions.

        The function is similar to getVoxDims(), but instead of the 3d spatial
        dimensions of a voxel it returns the dimensions of an image pixel on
        all 7 dimensions supported by the NIfTI dataformat.
        """
        return tuple( 
                    [ nifticlib.floatArray_frompointer(self.__nimg.pixdim)[i] 
                      for i in range(1,8) ]
                )


    def getExtent(self):
        """ Returns a tuple describing the shape (size in voxel/timepoints)
        of the dataimage.

        The order of dimensions is (x,y,z,t,u,v,w). If the image has less 
        dimensions than 7 the return tuple will be shortened accordingly.

        Please note that the order of dimensions is different from the tuple
        returned by calling NiftiImage.data.shape!

        See also getVolumeExtent() and getTimepoints().

        The 'extent' property is an alternative way to access this function.
        """
        # wrap dim array in nifti image struct
        dims_array = nifticlib.intArray_frompointer(self.__nimg.dim)
        dims = [ dims_array[i] for i in range(8) ]

        return tuple( dims[1:dims[0]+1] )


    def getVolumeExtent(self):
        """ Returns the size/shape of the volume(s) in the image as a tuple.

        This method returns either a 3-tuple or 2-tuple or 1-tuple depending 
        on the available dimensions in the image.

        The order of dimensions in the tuple is (x [, y [, z ] ] ).

        The 'volextent' property is an alternative way to access this function.
        """

        # it is save to do this even if self.extent is shorter than 4 items
        return self.extent[:3]


    def getTimepoints(self):
        """ Returns the number of timepoints in the image.

        In case of a 3d (or less dimension) image this method returns 1.

        The 'timepoints' property is an alternative way to access this
        function.
        """

        if len(self.extent) < 4:
            return 1
        else:
            return self.extent[3]


    def getRepetitionTime(self):
        """ Returns the temporal distance between the volumes in a timeseries.

        The 'rtime' property is an alternative way to access this function.
        """
        return self.__nimg.dt


    def setRepetitionTime(self, value):
        """ Set the repetition time of a nifti image (dt).
        """
        self.__nimg.dt = float(value)


    def getHeader(self):
        """ Returns the header data of the nifti image in a dictionary.

        Note, that modifications done to this dictionary do not cause any 
        modifications in the NIfTI image. Please use the updateHeader() method
        to apply changes to the image.

        The 'header' property is an alternative way to access this function. 
        But please note that the 'header' property cannot be used like this:

            nimg.header['something'] = 'new value'

        Instead one has to get the header dictionary, modify and later reassign
        it:

            h = nimg.header
            h['something'] = 'new value'
            nimg.header = h
        """
        h = {}

        # Convert nifti_image struct into nifti1 header struct.
        # This get us all data that will actually make it into a
        # NIfTI file.
        nhdr = nifticlib.nifti_convert_nim2nhdr(self.__nimg)

        return NiftiImage.nhdr2dict(nhdr)


    def updateHeader(self, hdrdict):
        """ Update NIfTI header information.

        Updated header data is read from the supplied dictionary. One cannot
        modify dimensionality and datatype of the image data. If such
        information is present in the header dictionary it is removed before
        the update. If resizing or datatype casting are required one has to 
        convert the image data into a separate array 
        ( NiftiImage.assarray(copy=True) ) and perform resize and data 
        manipulations on this array. When finished, the array can be converted
        into a nifti file by calling the NiftiImage constructor with the
        modified array as 'source' and the nifti header of the original
        NiftiImage object as 'header'.

        It is save to call this method with and without loaded image data.

        The actual update is done by NiftiImage.updateNiftiHeaderFromDict().

        Besides reading it is also possible to set the header data by assigning
        to the 'header' property. Please see the documentation of the 
        getHeader() method for important information about the special usage
        of the 'header' property.
        """
        # rebuild nifti header from current image struct
        nhdr = nifticlib.nifti_convert_nim2nhdr(self.__nimg)

        # remove settings from the hdrdict that are determined by
        # the data set and must not be modified to preserve data integrity
        if hdrdict.has_key('datatype'):
            del hdrdict['datatype']
        if hdrdict.has_key('dim'):
            del hdrdict['dim']

        # update the nifti header
        NiftiImage.updateNiftiHeaderFromDict(nhdr, hdrdict)

        # if no filename was set already (e.g. image from array) set a temp
        # name now, as otherwise nifti_convert_nhdr2nim will fail
        have_temp_filename = False
        if not self.filename:
            self.filename = 'pynifti_updateheader_temp_name'
            have_temp_filename = True

        # recreate nifti image struct
        new_nimg = nifticlib.nifti_convert_nhdr2nim(nhdr, self.filename)
        if not new_nimg:
            raise RuntimeError, \
                  "Could not recreate NIfTI image struct from updated header."

        # replace old image struct by new one
        # be careful with memory leak (still not checked whether successful)

        # rescue data ptr
        new_nimg.data = self.__nimg.data

        # and remove it from old image struct
        self.__nimg.data = None

        # to be able to call the cleanup function without lossing the data
        self.__close()

        # assign the new image struct
        self.__nimg = new_nimg

        # reset filename if temp name was set
        if have_temp_filename:
            self.filename = ''


    def setSlope(self, value):
        """ Set the slope attribute in the NIfTI header.

        Besides reading it is also possible to set the slope by assigning
        to the 'slope' property.
        """
        self.__nimg.scl_slope = float(value)


    def setIntercept(self, value):
        """ Set the intercept attribute in the NIfTI header.

        Besides reading it is also possible to set the intercept by assigning
        to the 'intercept' property.
        """
        self.__nimg.scl_inter = float(value)


    def setDescription(self, value):
        """ Set the description element in the NIfTI header.

        Descriptions must not be longer than 79 characters.

        Besides reading it is also possible to set the description by assigning
        to the 'description' property.
        """
        if len(value) > 79:
            raise ValueError, \
                  "The NIfTI format only supports descriptions shorter than " \
                  + "80 chars."

        self.__nimg.descrip = value


    def getSForm(self):
        """ Returns the sform matrix.

        The 'sform' property is an alternative way to access this function.

        Please note, that the returned SForm matrix is not bound to the
        NiftiImage object. Therefore it cannot be successfully modified
        in-place. Modifications to the SForm matrix can only be done by setting
        a new SForm matrix either by calling setSForm() or by assigning it to
        the sform attribute.
        """
        return nifticlib.mat442array(self.__nimg.sto_xyz)


    def setSForm(self, m):
        """ Sets the sform matrix.
        The supplied value has to be a 4x4 matrix. The matrix elements will be
        converted to floats. By definition the last row of the sform matrix has
        to be (0,0,0,1). However, different values can be assigned, but will
        not be stored when the niftifile is saved.

        The inverse sform matrix will be automatically recalculated.

        Besides reading it is also possible to set the sform matrix by
        assigning to the 'sform' property.
        """
        if m.shape != (4,4):
            raise ValueError, "SForm matrix has to be of size 4x4."

        # make sure it is float
        m = m.astype('float')

        nifticlib.set_mat44( self.__nimg.sto_xyz,
                         m[0,0], m[0,1], m[0,2], m[0,3],
                         m[1,0], m[1,1], m[1,2], m[1,3],
                         m[2,0], m[2,1], m[2,2], m[2,3],
                         m[3,0], m[3,1], m[3,2], m[3,3] )

        # recalculate inverse
        self.__nimg.sto_ijk = \
            nifticlib.nifti_mat44_inverse( self.__nimg.sto_xyz )


    def getInverseSForm(self):
        """ Returns the inverse sform matrix.

        The 'sform_inv' property is an alternative way to access this function.

        Please note, that the inverse SForm matrix cannot be modified in-place.
        One needs to set a new SForm matrix instead. The corresponding inverse
        matrix is then re-calculated automatically.
        """
        return nifticlib.mat442array(self.__nimg.sto_ijk)


    def getQForm(self):
        """ Returns the qform matrix.

        The 'qform' property is an alternative way to access this function.

        Please note, that the returned QForm matrix is not bound to the
        NiftiImage object. Therefore it cannot be successfully modified
        in-place. Modifications to the QForm matrix can only be done by setting
        a new QForm matrix either by calling setSForm() or by assigning it to
        the sform attribute.
        """
        return nifticlib.mat442array(self.__nimg.qto_xyz)


    def getInverseQForm(self):
        """ Returns the inverse qform matrix.

        The 'qform_inv' property is an alternative way to access this function.

        Please note, that the inverse QForm matrix cannot be modified in-place.
        One needs to set a new QForm matrix instead. The corresponding inverse
        matrix is then re-calculated automatically.
        """
        return nifticlib.mat442array(self.__nimg.qto_ijk)


    def setQForm(self, m):
        """ Sets the qform matrix.
        The supplied value has to be a 4x4 matrix. The matrix will be converted
        to float.

        The inverse qform matrix and the quaternion representation will be
        automatically recalculated.

        Besides reading it is also possible to set the qform matrix by
        assigning to the 'qform' property.
        """
        if m.shape != (4,4):
            raise ValueError, "QForm matrix has to be of size 4x4."

        # make sure it is float
        m = m.astype('float')

        nifticlib.set_mat44( self.__nimg.qto_xyz,
                         m[0,0], m[0,1], m[0,2], m[0,3],
                         m[1,0], m[1,1], m[1,2], m[1,3],
                         m[2,0], m[2,1], m[2,2], m[2,3],
                         m[3,0], m[3,1], m[3,2], m[3,3] )

        # recalculate inverse
        self.__nimg.qto_ijk = \
            nifticlib.nifti_mat44_inverse( self.__nimg.qto_xyz )

        # update quaternions
        ( self.__nimg.quatern_b, self.__nimg.quatern_c, self.__nimg.quatern_d,
          self.__nimg.qoffset_x, self.__nimg.qoffset_y, self.__nimg.qoffset_z,
          self.__nimg.dx, self.__nimg.dy, self.__nimg.dz,
          self.__nimg.qfac ) = \
            nifticlib.nifti_mat44_to_quatern( self.__nimg.qto_xyz )


    def updateQFormFromQuaternion(self):
        """ Recalculates the qform matrix (and the inverse) from the quaternion
        representation.
        """
        # recalculate qform
        self.__nimg.qto_xyz = nifticlib.nifti_quatern_to_mat44 (
          self.__nimg.quatern_b, self.__nimg.quatern_c, self.__nimg.quatern_d,
          self.__nimg.qoffset_x, self.__nimg.qoffset_y, self.__nimg.qoffset_z,
          self.__nimg.dx, self.__nimg.dy, self.__nimg.dz,
          self.__nimg.qfac )


        # recalculate inverse
        self.__nimg.qto_ijk = \
            nifticlib.nifti_mat44_inverse( self.__nimg.qto_xyz )


    def setQuaternion(self, value):
        """ Set Quaternion from 3-tuple (qb, qc, qd).

        The qform matrix and its inverse are re-computed automatically.

        Besides reading it is also possible to set the quaternion by assigning
        to the 'quatern' property.
        """
        if len(value) != 3:
            raise ValueError, 'Requires 3-tuple.'

        self.__nimg.quatern_b = float(value[0])
        self.__nimg.quatern_c = float(value[1])
        self.__nimg.quatern_d = float(value[2])

        self.updateQFormFromQuaternion()


    def getQuaternion(self):
        """ Returns a 3-tuple containing (qb, qc, qd).

        The 'quatern' property is an alternative way to access this function.
        """
        return( ( self.__nimg.quatern_b, 
                  self.__nimg.quatern_c, 
                  self.__nimg.quatern_d ) )


    def setQOffset(self, value):
        """ Set QOffset from 3-tuple (qx, qy, qz).

        The qform matrix and its inverse are re-computed automatically.

        Besides reading it is also possible to set the qoffset by assigning
        to the 'qoffset' property.
        """
        if len(value) != 3:
            raise ValueError, 'Requires 3-tuple.'

        self.__nimg.qoffset_x = float(value[0])
        self.__nimg.qoffset_y = float(value[1])
        self.__nimg.qoffset_z = float(value[2])

        self.updateQFormFromQuaternion()


    def getQOffset(self):
        """ Returns a 3-tuple containing (qx, qy, qz).

        The 'qoffset' property is an alternative way to access this function.
        """
        return( ( self.__nimg.qoffset_x,
                  self.__nimg.qoffset_y,
                  self.__nimg.qoffset_z ) )


    def setQFac(self, value):
        """ Set qfac.

        The qform matrix and its inverse are re-computed automatically.

        Besides reading it is also possible to set the qfac by assigning
        to the 'qfac' property.
        """
        self.__nimg.qfac = float(value)
        self.updateQFormFromQuaternion()


    def getQOrientation(self, as_string = False):
        """ Returns to orientation of the i,j and k axis as stored in the
        qform matrix.

        By default NIfTI orientation codes are returned, but if 'as_string' is
        set to true a string representation ala 'Left-to-right' is returned
        instead.
        """
        codes = nifticlib.nifti_mat44_to_orientation(self.__nimg.qto_xyz)
        if as_string:
            return [ nifticlib.nifti_orientation_string(i) for i in codes ]
        else:
            return codes


    def getSOrientation(self, as_string = False):
        """ Returns to orientation of the i,j and k axis as stored in the
        sform matrix.

        By default NIfTI orientation codes are returned, but if 'as_string' is
        set to true a string representation ala 'Left-to-right' is returned
        instead.
        """
        codes = nifticlib.nifti_mat44_to_orientation(self.__nimg.sto_xyz)
        if as_string:
            return [ nifticlib.nifti_orientation_string(i) for i in codes ]
        else:
            return codes


    def getBoundingBox(self):
        """ Get the bounding box of the image.

        This functions returns a tuple of (min, max) tuples. It contains as
        many tuples as image dimensions. The order of dimensions is identical
        to that in the data array.

        The 'bbox' property is an alternative way to access this function.
        """
        nz = self.data.squeeze().nonzero()

        bbox = []

        for dim in nz:
            bbox.append( ( dim.min(), dim.max() ) )

        return tuple(bbox)


    def setFilename(self, filename, filetype = 'NIFTI'):
        """ Set the filename for the NIfTI image.

        Setting the filename also determines the filetype. If the filename
        ends with '.nii' the type will be set to NIfTI single file. A '.hdr'
        extension can be used for NIfTI file pairs. If the desired filetype
        is ANALYZE the extension should be '.img'. However, one can use the
        '.hdr' extension and force the filetype to ANALYZE by setting the
        filetype argument to ANALYZE. Setting filetype if the filename
        extension is '.nii' has no effect, the file will always be in NIFTI
        format.

        If the filename carries an additional '.gz' the resulting file(s) will
        be compressed.

        Uncompressed NIfTI single files are the default filetype that will be
        used if the filename has no valid extension. The '.nii' extension is
        appended automatically. The 'filetype' argument can be used to force a
        certain filetype when no extension can be used to determine it. 
        'filetype' can be one of the nifticlibs filtetypes or any of 'NIFTI',
        'NIFTI_GZ', 'NIFTI_PAIR', 'NIFTI_PAIR_GZ', 'ANALYZE', 'ANALYZE_GZ'.

        Setting the filename will cause the image data to be loaded into memory
        if not yet done already. This has to be done, because without the
        filename of the original image file there would be no access to the
        image data anymore. As a side-effect a simple operation like setting a
        filename may take a significant amount of time (e.g. for a large 4d
        dataset).

        By passing an empty string or none as filename one can reset the
        filename and detach the NiftiImage object from any file on disk.

        Examples:

          Filename          Output of save()
          ----------------------------------
          exmpl.nii         exmpl.nii (NIfTI)
          exmpl.hdr         exmpl.hdr, exmpl.img (NIfTI)
          exmpl.img         exmpl.hdr, exmpl.img (ANALYZE)
          exmpl             exmpl.nii (NIfTI)
          exmpl.hdr.gz      exmpl.hdr.gz, exmpl.img.gz (NIfTI)

        ! exmpl.gz          exmpl.gz.nii (uncompressed NIfTI)

        Setting the filename is also possible by assigning to the 'filename'
        property.
        """
        # If image data is not yet loaded, do it now.
        # It is important to do it already here, because nifti_image_load
        # depends on the correct filename set in the nifti_image struct
        # and this will be modified in this function!
        if not self.__haveImageData():
            self.load()

        # if no filename is given simply reset it to nothing
        if not filename:
            self.__nimg.fname = ''
            self.__nimg.iname = ''
            return

        # separate basename and extension
        base, ext = NiftiImage.splitFilename(filename)

        # if no extension default to nifti single files
        if ext == '': 
            if filetype == 'NIFTI' \
               or filetype == nifticlib.NIFTI_FTYPE_NIFTI1_1:
                ext = 'nii'
            elif filetype == 'NIFTI_PAIR' \
                 or filetype == nifticlib.NIFTI_FTYPE_NIFTI1_2:
                ext = 'hdr'
            elif filetype == 'ANALYZE' \
                 or filetype == nifticlib.NIFTI_FTYPE_ANALYZE:
                ext = 'img'
            elif filetype == 'NIFTI_GZ':
                ext = 'nii.gz'
            elif filetype == 'NIFTI_PAIR_GZ':
                ext = 'hdr.gz'
            elif filetype == 'ANALYZE_GZ':
                ext = 'img.gz'
            else:
                raise RuntimeError, "Unhandled filetype."

        # Determine the filetype and set header and image filename
        # appropriately.

        # nifti single files are easy
        if ext == 'nii.gz' or ext == 'nii':
            self.__nimg.fname = base + '.' + ext
            self.__nimg.iname = base + '.' + ext
            self.__nimg.nifti_type = nifticlib.NIFTI_FTYPE_NIFTI1_1
        # uncompressed nifti file pairs
        elif ext in [ 'hdr', 'img' ]:
            self.__nimg.fname = base + '.hdr'
            self.__nimg.iname = base + '.img'
            if ext == 'hdr' and not filetype.startswith('ANALYZE'):
                self.__nimg.nifti_type = nifticlib.NIFTI_FTYPE_NIFTI1_2
            else:
                self.__nimg.nifti_type = nifticlib.NIFTI_FTYPE_ANALYZE
        # compressed file pairs
        elif ext in [ 'hdr.gz', 'img.gz' ]:
            self.__nimg.fname = base + '.hdr.gz'
            self.__nimg.iname = base + '.img.gz'
            if ext == 'hdr.gz' and not filetype.startswith('ANALYZE'):
                self.__nimg.nifti_type = nifticlib.NIFTI_FTYPE_NIFTI1_2
            else:
                self.__nimg.nifti_type = nifticlib.NIFTI_FTYPE_ANALYZE
        else:
            raise RuntimeError, "Unhandled filetype."


    def getFilename(self):
        """ Returns the filename.

        To be consistent with setFilename() the image filename is returned
        for ANALYZE images while the header filename is returned for NIfTI
        files.

        The 'filename' property is an alternative way to access this function.
        """
        if self.__nimg.nifti_type == nifticlib.NIFTI_FTYPE_ANALYZE:
            return self.__nimg.iname
        else:
            return self.__nimg.fname

    # class properties
    # read only
    nvox =          property(fget=lambda self: self.__nimg.nvox)
    max =           property(fget=lambda self: self.__nimg.cal_max)
    min =           property(fget=lambda self: self.__nimg.cal_min)
    data =          property(fget=getDataArray)
    sform_inv =     property(fget=getInverseSForm)
    qform_inv =     property(fget=getInverseQForm)
    extent =        property(fget=getExtent)
    volextent =     property(fget=getVolumeExtent)
    timepoints =    property(fget=getTimepoints)
    bbox =          property(fget=getBoundingBox)

    # read and write
    filename =      property(fget=getFilename, fset=setFilename)
    slope =         property(fget=lambda self: self.__nimg.scl_slope,
                             fset=setSlope)
    intercept =     property(fget=lambda self: self.__nimg.scl_inter,
                             fset=setIntercept)
    voxdim =        property(fget=getVoxDims, fset=setVoxDims)
    pixdim =        property(fget=getPixDims, fset=setPixDims)
    description =   property(fget=lambda self: self.__nimg.descrip,
                             fset=setDescription)
    header =        property(fget=getHeader, fset=updateHeader)
    sform =         property(fget=getSForm, fset=setSForm)
    qform =         property(fget=getQForm, fset=setQForm)
    quatern =       property(fget=getQuaternion, fset=setQuaternion)
    qoffset =       property(fget=getQOffset, fset=setQOffset)
    qfac =          property(fget=lambda self: self.__nimg.qfac, fset=setQFac)
    rtime =         property(fget=getRepetitionTime, fset=setRepetitionTime)

