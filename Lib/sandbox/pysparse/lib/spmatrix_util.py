import math, random
import spmatrix
import pysparse_version

def bytesToString(n):
    if n < 1024:
        return '%d Bytes' % n
    n /= 1024.0
    if n < 1024:
        return '%.1f Kbytes' % n
    n /= 1024.0
    if n < 1024:
        return '%.1f Mbytes' % n
    n /= 1024.0
    return '%.1f Gbytes' % n


def printInfo(mat, name):
    
    if type(mat) == spmatrix.LLMatType:
        if mat.issym:
            typeName = 'LL symmetric'
        else:
            typeName = 'LL general'
        storage = bytesToString(4*mat.shape[0] + 16*mat.nnz)
    elif type(mat) == spmatrix.SSSMatType:
        typeName = 'SSS'
        storage = bytesToString(12*mat.shape[0] + 12*(mat.nnz - mat.shape[0]))
    elif type(mat) == spmatrix.CSRMatType:
        typeName = 'CSR'
        storage = bytesToString(4*mat.shape[0] + 12*mat.nnz)
    else:
        typeName = 'Unknown'
        storage = 'Unknown'


    print 'Matrix statistics:'
    print '%-20s: %s' % ('name', name)
    print '%-20s: %s' % ('type', typeName)
    print '%-20s: %dx%d' % ('dimensions', mat.shape[0], mat.shape[1])
    print '%-20s: %d' % ('#non-zeros', mat.nnz)
    print '%-20s: %s' % ('storage', storage)
    print

def ll_mat_rand(n, m, density):
    """return a ll_mat object representing a general n-by-m sparse matrix filled with random non-zero values

    The number of non-zero entries is less or equal than
    n*m*density. The values of the non-zero entries are in the range
    [0.0, 1.0)."""
    nnz = int(density*n*m)
    A = spmatrix.ll_mat(n, m, nnz)
    for k in xrange(nnz):
        i = random.randrange(n)
        j = random.randrange(m)
        A[i, j] = random.random()
    return A

def exportVtk(mat, fileName):
    "export matrix to a VTK data file"

    print 'Write VTK file...'
    # write VTK file
    f = open(fileName, 'w')
    f.write('# vtk DataFile Version 3.0\n')
    comment = 'generated using pysparse-%s\n' % (pysparse_version.version, )
    f.write(comment[:256])
    f.write('ASCII\n\n')
    f.write('DATASET STRUCTURED_POINTS\n')
    f.write('DIMENSIONS %d %d 1\n' % (mat.shape[0], mat.shape[1]))
    f.write('ORIGIN 1 1 1\n')
    f.write('SPACING 1 1 1\n\n')
    f.write('\nPOINT_DATA %d\n' % (mat.shape[0]*mat.shape[1]))
    f.write('SCALARS entries float\n')
    f.write('LOOKUP_TABLE default\n\n')
    for i in xrange(mat.shape[0]):
        for j in xrange(mat.shape[1]):
            v = mat[i,j]
            if v <> 0.0:
                v = math.log(math.fabs(v))
            f.write('%lf\n' % v)
    f.close()

def viewVtk(mat):
    import vtk

    imageData = vtk.vtkImageData()
    imageData.SetDimensions(mat.shape[0], mat.shape[1], 1)
    imageData.SetScalarTypeToUnsignedShort()
    imageData.SetNumberOfScalarComponents(1)
    imageData.AllocateScalars()
    
    viewer = vtk.vtkImageViewer()
    viewer.SetInput(imageData)
    viewer.SetZSlice(0)
    
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(viewer.GetRenderWindow())

    iren.Initialize()
    iren.Start()
    
    imageData = vtk.vtkImageSinusoidSource()
    imageData.SetWholeExtent(0, 300, 0, 300, 0, 10)
    imageData.SetAmplitude(63)
    imageData.SetDirection(1, 1, 0)
    imageData.SetPeriod(25)

    viewer = vtk.vtkImageViewer()
    viewer.SetInput(imageData.GetOutput())
    viewer.SetColorWindow(126)
    viewer.SetColorLevel(0)
    viewer.SetZSlice(0)

    def hamschti(obj, event):
        print 'Haam:'
        print obj
        print event

    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(viewer.GetRenderWindow())
    iren.Initialize()

    interactor = vtk.vtkInteractorStyleImage()
    ##interactor.AddObserver('LeftButtonPressEvent', hamschti)
    iren.SetInteractorStyle(interactor)

    iren.Start()
