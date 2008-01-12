
""" Tests for BuildGrid.c module.
To view the grids I used QuikGrid from
http://www.PerspectiveEdge.com 
"""

import sys,math,random,time
from scipy.sandbox import buildgrid

random.seed(7)

global ncol,nrow,step,xmin,ymin,simple
xmin,ymin,step = 0.0,0.0,1.0 

# Surface model:
simple = True   #  simple: (a1+b1*x)*(a2+b2*y)
                # This is almost a plain.

simple = False  # !simple: cos(a*x*x)*cos(b*y*y)
                # In this case only left bottom part of surface
                # should be believable.

def gridsizes(nodes,ratio=1.0):
    # ratio = ncol/nrow
    nodes = max(100,nodes)
    ratio = max(0.1,min(ratio,10.0))
    nrow = int(math.sqrt(nodes/ratio))
    ncol = nodes/nrow+1
    return ncol,nrow

ncol,nrow = gridsizes(10**5) # I used 10**3 - 10**7 nodes

ax,ay = math.pi/(4*ncol+1),math.pi/(4*nrow+1)

percent = 10.2 # how many data points w.r.t. grid nodes:
datapoints = int(0.01 * percent * ncol * nrow) or 5
print "Testing",("complex","simple")[simple],
print "surface. Grid size %d x %d (%d nodes)" % (ncol,nrow,ncol*nrow)
print "About %d datapoints (%.2f%% of nodes)" % (datapoints,percent)

# Test flag
build_only = True  #  no output to files
build_only = False #  output to files

# Trimming distance
trimdist = 5*step  # make trimming
trimdist = 0.0     # do not make trimming
if trimdist < step:
    print "No trimming"
else:
    print "With trimming",trimdist

def surface(x,y):
    if simple:
        return 100+100*(1+float(x)/ncol)*(1+float(y)/nrow)
    return 100+50*(math.cos(ax*x*x)*math.cos(ay*y*y)+1)
    
def surfaceGen():
    while True:
        if trimdist < step:
            x = random.uniform(-1,ncol) # some points may be 
            y = random.uniform(-1,nrow) # out of grid
        else:
            trdist = int(trimdist/step+0.5)
            x = random.uniform(-1+trdist,ncol-trdist) # all the points 
            y = random.uniform(-1+trdist,nrow-trdist) # are inside the grid
        z = surface(x,y)
        yield x,y,z

def makeInputXYZ(outfile,nvals):
    surgen = surfaceGen()
    x,y,z = [],[],[]
    fo = open(outfile,"wt") # for fromfile()
    while(nvals > 0):
        xv,yv,zv = surgen.next()
        fo.write("%.3f %.3f %.3f\n" % (xv,yv,zv))
        x.append(xv) # for fromxyz()
        y.append(yv)
        z.append(zv)
        nvals -= 1
    fo.close()
    return x,y,z

# create input data file with
xyzfile = "inpdata.xyz" # input for fromfile()
xp,yp,zp = makeInputXYZ(xyzfile, datapoints)



# get file statistics
nvals,xmin,xmax,ymin,ymax,zmin,zmax,zavrg,zstnd = \
    buildgrid.filestatistics(xyzfile)
print xyzfile,'statistics.  Number of values',nvals
print 'X min %9.1f   max %9.1f' % (xmin,xmax)
print 'Y min %9.1f   max %9.1f' % (ymin,ymax)
print 'S min %9.1f   max %9.1f' % (zmin,zmax)
print 'S avrg %7.2f   stnd %7.3f' % (zavrg,zstnd)

# build grid 'fromfile'
t0 = time.clock()
grid = buildgrid.fromfile(xyzfile=xyzfile,
    nx=ncol, ny=nrow, step=step, xmin=xmin, ymin=ymin, 
    method='Good',      # or 'Best' - not implemented
    trimming=trimdist,  # if no trimming - full grid
    unvalue=123.321,    # will be used if trimming > 0
    abserror=0.001,     # may be == 0
    relerror=0.0001)    # ignores if abserror > 0
print 'fromfile():',xyzfile,len(grid)-ncol*nrow,"(%.2f sec, %.0f nodes/sec)" %\
    (time.clock()-t0,ncol*nrow/(time.clock()-t0))


# put grid data to xyz file
if not build_only:
    outfile = "outdata1.xyz"
    t0 = time.clock()
    rv = buildgrid.tofile(filename=outfile, griddata=grid, gridtype='xyz',
        nx=ncol, ny=nrow, step=step, xmin=xmin, ymin=ymin, unvalue=1234.4321)
    print 'tofile():',outfile,rv,"(%.2f sec)" % (time.clock()-t0)


# build grid from xyz lists
t0 = time.clock()
grid2 = buildgrid.fromxyz(xdata=xp,ydata=yp,zdata=zp,
    nx=ncol, ny=nrow, step=step, xmin=xmin, ymin=ymin, 
    method='Good',      # or 'Best' (not implemented)
    trimming=trimdist,  # if no trimming - full grid
    unvalue=123.321,    # will be used if trimming > 0
    abserror=0.001,     # may be == 0
    relerror=0.0001)    # ignores if abserror > 0
print 'fromxyz():',len(grid)-ncol*nrow,"(%.2f sec, %.0f nodes/sec)" %\
    (time.clock()-t0,ncol*nrow/(time.clock()-t0))

# put grid to file
if not build_only:
    outfile = "outdata2.xyz"
    t0 = time.clock()
    rv = buildgrid.tofile(filename=outfile, griddata=grid2, gridtype='xyz',
        nx=ncol, ny=nrow, step=step, xmin=xmin, ymin=ymin, unvalue=1234.4321)
    print 'tofile():',outfile,rv,"(%.2f sec)" % (time.clock()-t0)


# write full exact grid
if not build_only:
    exactfile = "exactdata.xyz"
    t0 = time.clock()
    fo = open(exactfile,"wt")
    for nr in xrange(nrow):
        for nc in xrange(ncol):
            fo.write("%.2f %.2f %.3f\n" % (nc,nr,surface(nc,nr)))
    fo.close()
    print 'exact file:',exactfile,"(%.2f sec)" % (time.clock()-t0)


if build_only:
    sys.exit()

# Compare Input & Output:
layer = max(5,min(ncol,nrow)/10)
def avst(left,rite,low,top):
    av = st = 0.0
    nv = 0
    for nr in xrange(low,top):
        for nc in xrange(left,rite):
            d = grid[nr*ncol+nc]-surface(nc,nr)
            av += d
            st += d*d
            nv += 1
    av /= nv
    st = st/nv-av*av
    if st > 0.0: st = math.sqrt(st)
    return av,st

print "Comparing %d x %d squares" % (layer,layer),"of grid lines"
av,st = avst(0,layer,0,layer)
print "Bottom Left:   av = %8.3f  st = %8.3f" % (av,st)
av,st = avst(ncol-layer,ncol,0,layer)
print "Bottom Right:  av = %8.3f  st = %8.3f" % (av,st)
av,st = avst(0,layer,nrow-layer,nrow)
print "Top Left:      av = %8.3f  st = %8.3f" % (av,st)
av,st = avst(ncol-layer,ncol,nrow-layer,nrow)
print "Top Right:     av = %8.3f  st = %8.3f" % (av,st)
av,st = avst(ncol/2-layer/2,ncol/2+layer/2,nrow/2-layer/2,nrow/2+layer/2)
print "Middle:        av = %8.3f  st = %8.3f" % (av,st)
    
