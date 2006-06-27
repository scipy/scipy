/* 
  Build regular grids from scattered 2D data.

  Author: Eugene Druker, 2006 (eugene.druker@gmail.com)

Summary
========  
Input:
    a set of arbitrary displaced points in 2D plane,
    grid geometry: x/y min/max and cell size, 
    a few (optional) control parameters 

Output:
    to a python list or to file

Method: 
    Put every input data point into nearest grid node (take average if more
    than one), this is main source of errors. Apply empirical method of 
    building: make sequence of averagings with full interpolation. Use 
    distances from data points to remove far field values.
    The method is fair for large data volumes (and grid sizes). 
    The method is irrelevant if (x,y) must be treated as exact values.


Details
========
Input data and parameters:
   - (x,y,z) in file or in lists. Text XYZ-file contains lines of x, y, z 
     values, space divided. Lists are separate for x, y, z data. No special 
     order of (x,y,z) triples is supposed. Sizes, distances, coordinates 
     are in same units.
   - grid boundaries (xmin,...,ymax) and (square) cell size
   - grid building method: 'Good' (default) or 'Best' (not implemented)
   - acceptable error level - absolute or relative. Default is 0.0 (but it
     is not an 'interpolation')
   - trimming distance: throw out (i.e. replace by non-values) all the values 
     with distances from sources larger than trimdistance, starting from grid 
     boundaries (i.e. internal nodes may be untouched)
   - unvalue, use for output, to mark nodes without values (default: 12345678.9)
   - a few others ?

Registered functions (to call from python):

   filestatistics(xyzfile=filename)
		learn file statistics, on x,y,z.
		Return tuple (nvals,xmin,xmax,ymin,ymax,zmin,zmax,zavrg,zstnd) of
        statistics: nvals = number of values, *min = minimum, *max = maximum,
        *avrg = average, *stnd = standard, where * is x/y/z.
		Raise IOError -- could not open file.

   fromfile(
        xyzfile=filename,       // input file of lines (x,y,z)
		nx=xnodes, ny=ynodes,   // grid sizes in nodes
		step=cellsize,          // cell sizes
		xmin=minX, ymin=minY,   // left bottom, x - east, y - north
		method='Good',          // or 'Best' (not implemented yet)
		trimming=trimdistance,  // distance to trim, ignored if trimming < step
		unvalue=non_value,      // caller unvalue, ignored if no trimming
        abserror=AbsError,      // acceptable error if > 0
        relerror=RelError,      // the same, as part of standard (if abserror==0)
                                // relerror ignored if abserror > 0.0
		)
        Return grid as list, rows from bottom, columns from left
		Raise ValueError/IOError for errors with file/memory

   fromxyz(
        xdata=x, ydata=y, zdata=z, // lists of input data (of same length)
		nx=xnodes, ny=ynodes,   // grid sizes in nodes
		step=cellsize,          // cell sizes
		xmin=minX, ymin=Ymin,   // left bottom, x - east, y - north
        # Optional parameters:
		method='Good',          // or 'Best' (not implemented)
		trimming=trimdistance,  // distance to trim, ignored if trimming < step
		unvalue=non_value,      // caller unvalue, ignored if no trimming
        abserror=AbsError,      // acceptable error if > 0
        relerror=RelError,      // the same, as part of standard (if abserror==0)
                                // relerror ignored if abserror > 0.0
        )
		Return grid as list, rows from bottom, columns from left
		Raise ValueError/IOError for errors with values/memory
        
   tofile(
        outfile=filename,       // output file name
		griddata=xyz,           // list as returned by fromfile() or fromxyz()
		nx=xnodex, ny=ynodes,   // grid sizes in nodes
        filetype=grid_type,     // 'xyz' or 'gxf' or 'grd' etc (not implemented) 
        # Optional parameters:
		step=cellsize,          // cell sizes, default: 1.0
		xmin=minX, ymin=Ymin,   // left bottom, x - east, y - north
                                // default: xmin=ymin=0.0
		unvalue=non_value,      // caller unvalue, default: 12345678.9
                                // grid nodes with unvalue will not be put to 'xyz' file
        )
		Return 0 for Ok
		Raise ValueError/IOError for errors with file/memory

Method - sizes/memory/speed:
    minimal grid size is 9 nodes (in X and Y), maximal size is about 10000 
    (in this case it will take about 300 MB of memory, for data only).
    Output data take 8*ncol*nrow bytes - it is minimum (for doubles in nodes).
    Maximum is about 30*ncol*nrow bytes - during the building.
    Speed is about 100000 nodes/sec (1.6 GHz, 1 GB RAM), there is also a weak
    dependence on input data volume, precision, etc. The most important 
    parameter for speed is trimming (currently).

Tests:
    1. this text was successfully compiled as C (not C++) with MS VC++ 6.0,
        MinGW Dev.Studio 2.05 and gcc-4 (ubuntu).
    2. python test program is BuildGrid.py
        run it to see numerical/graphical results for two examples
*/

#include "Python.h"

#include <stdio.h>   /* FILE, fgets(), printf() */
#include <math.h>    /* fabs() */
#include <string.h>  /* memset() */

#define NoValue                 32767     /* unreachable distance for short */
#define Margin                      7     /* not used (so far) */
#define MinGridSize                 9     /* minimal grid size in nodes */
#define MaxGridSize   (NoValue/3-Margin)  /* maximal grid size in nodes */
#define MinStep                 0.001     /* minimal grid cell size */

/* unvalue should be: unvalue != 1.01 * unvalue, for floats and doubles */
#define UnValue        12345678.9         /* double and float */
#define isUnvalue(v) (fabs((v)-UnValue) < 0.01)

#define TopLimitValue  (1.2E34)
#define LowLimitValue  (-1.2E34)

#define TextLen     256 /* text lines */

/* Grid sizes are limited by about 9*3**ZIO nodes in X & Y */
#define ZIO  12       /* Number of levels */

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

//=====================================================
//=== input XYZ file statistics =======================
//=====================================================
typedef struct _File_Statistics_ 
{
	int nvals;
	double xmin,xmax,ymin,ymax,zmin,zmax,zavrg,zstnd;
} FileStatistics;

static void
zerostatistics(FileStatistics * fs)
{
	fs->nvals = 0;
	fs->zavrg = fs->zstnd = 0.0;
	fs->xmin = fs->ymin = fs->zmin = TopLimitValue;
	fs->xmax = fs->ymax = fs->zmax = LowLimitValue;
}

static void
standard(int nvals, double * avrg, double * stnd)
{
	if(nvals > 0) {
		*avrg /= nvals;
		*stnd = *stnd/nvals - (*avrg) * (*avrg);
		if(*stnd > 0.0) *stnd = sqrt(*stnd);
        else *stnd = 0.0;
	} //else *stnd = 0.0;
}
    

static int
getfilestatistics(char * xyzfile, FileStatistics * fs)
{
	char text[TextLen];
	double x,y,z;
	int n;
	FILE *fi = fopen(xyzfile,"rt");
	if(!fi) return -1;

	while(fgets(text,TextLen,fi)) { // until EOF
		n = sscanf(text,"%lf %lf %lf",&x,&y,&z);
		if(n != 3) continue;
		++fs->nvals;
		fs->zavrg += z;
		fs->zstnd += z*z;
		fs->zmin = MIN(fs->zmin,z);
		fs->zmax = MAX(fs->zmax,z);
		fs->ymin = MIN(fs->ymin,y);
		fs->ymax = MAX(fs->ymax,y);
		fs->xmin = MIN(fs->xmin,x);
		fs->xmax = MAX(fs->xmax,x);
	}
	fclose(fi);
    standard(fs->nvals, &fs->zavrg, &fs->zstnd);
	return 0;
}

//--- Registered function ---------------------------
//--- get input XYZ file statistics -----------------
//---------------------------------------------------
static PyObject *
filestatistics(PyObject * self, PyObject * args, PyObject *kwds)
{
	int rv;
	FileStatistics fs;
	char *xyzfile;

	if(!PyArg_ParseTuple(args, "s", &xyzfile)) {
		PyErr_SetString(PyExc_IOError,"Filename is a must");
		return NULL;
	}
	zerostatistics(&fs);
	rv = getfilestatistics(xyzfile, &fs);
	if(rv < 0) {
		PyErr_SetString(PyExc_IOError,"No such file");
		return NULL;
	}
	return Py_BuildValue("idddddddd",fs.nvals,fs.xmin,fs.xmax,
		fs.ymin,fs.ymax,fs.zmin,fs.zmax,fs.zavrg,fs.zstnd);
}


//=====================================================
//=== build grid ======================================
//=====================================================
typedef struct _Grid_Field_ {
	int nx, ny;    // grid sizes in nodes
	double step,   // grid mesh is square
		xmin,ymin; // coordinates in X/Y directions

	double *field; // result
	short *index;  // input data index, then distance

	char method;     // 'G'ood or 'B'est (not implemented)
	double trimdist; // distance to trim the field
	double unvalue;  // caller unvalue - use for output only
	double maydiff;  // user-set tolerance

	double *grids[ZIO];
	short *indxs[ZIO];
	int nxs[ZIO], nys[ZIO];

	int indat; // input data points
	int indup; // duplicates
	int inout; // out of grid
	int inunv; // unvals
	int ingrd; // number of initial values in grid

} GridField;

/* The only global variable to keep data between calls.
   See initialization in initBuild() below */
static GridField gf;

// debug function
static void
debug_write(char *file, int version, double *field, int nx, int ny)
{
    //return;
    FILE *fo;
    int kx,ky;
    char filename[80];

    sprintf(filename,"outdata_%s_%d.fld",file,version);
    fo = fopen(filename,"wt");
    if(!fo) { printf("fail write to debug %s\n",filename); return; }
    for(ky=ny-1; ky >= 0; --ky) {
        for(kx=0; kx < nx; ++kx) 
            fprintf(fo,"%d %d %f\n",kx,ky,field[nx*ky+kx]);
    }
    fclose(fo);
}
    

static void
initGridField() // for initialization and after freeing memory
{
    int n;
	gf.method = 'G'; // Good
	gf.nx = gf.ny = 0;
	gf.step = 0.0;
	gf.unvalue = UnValue;
	gf.trimdist = 0.0;  // no trimming
    gf.maydiff = 0.0;   // exactly
	gf.xmin = gf.ymin = TopLimitValue;
	gf.indat = gf.inunv = gf.indup = gf.inout = gf.ingrd = 0;
	for(n=0; n < ZIO; ++n) {
		gf.grids[n] = NULL;
		gf.indxs[n] = NULL;
	}
	gf.index = NULL;
	gf.field = NULL;
}

static void
freeGridFieldMemory(int freefield)
{
	if(freefield & 1) {
		if(gf.field) { free(gf.field); gf.field = NULL; }
	}
	if(freefield & 2) {
        if(gf.index) { free(gf.index); gf.index = NULL; }
    }
	if(freefield & 4) {
        int n;
        for(n=1; n < ZIO; ++n) {
            if(gf.grids[n]) { free(gf.grids[n]); gf.grids[n] = NULL; }
            if(gf.indxs[n]) { free(gf.indxs[n]); gf.indxs[n] = NULL; }
        }
	}
}

static int
initGridFieldMemory()
{
    int k,kx,ky;
	freeGridFieldMemory(7);

	gf.field = (double *)malloc(gf.nx*gf.ny*sizeof(double));
	if(!gf.field) return -1;
	gf.index = (short *)malloc(gf.nx*gf.ny*sizeof(short));
	if(!gf.index) { freeGridFieldMemory(1); return -1; }
    memset((void *)gf.field,0,gf.nx*gf.ny*sizeof(double));
    memset((void *)gf.index,0,gf.nx*gf.ny*sizeof(short));
	gf.grids[0] = gf.field;
	gf.indxs[0] = gf.index;
	gf.nxs[0] = gf.nx;
	gf.nys[0] = gf.ny;
	for(k=1; k < ZIO; ++k) {
		kx = gf.nxs[k] = (gf.nxs[k-1]+2)/3;
		ky = gf.nys[k] = (gf.nys[k-1]+2)/3;
		if(kx < MinGridSize || ky < MinGridSize) break;
		gf.grids[k] = (double *)malloc(kx*ky*sizeof(double));
		gf.indxs[k] = (short *)malloc(kx*ky*sizeof(short));
		if(!gf.grids[k] || !gf.indxs[k]) return -1;
        memset((void *)gf.grids[k],0,kx*ky*sizeof(double));
        memset((void *)gf.indxs[k],0,kx*ky*sizeof(short));
	}
	return 0;
}

static int
getInputFileData(char *xyzfile, double * stnd) 
{
	char text[TextLen];
	int n,nn,indxy, xn,yn;
	double x,y,z, av,st;

	FILE *fi = fopen(xyzfile,"rt");
	if(!fi) return -1;
    av = st = 0.0;

	while(fgets(text,TextLen,fi)) {
		n = sscanf(text,"%lf %lf %lf",&x,&y,&z);
		if(n < 3) continue; // ignore other lines/errors
		++gf.indat; // number of input data points
		xn = (int)((x-gf.xmin)/gf.step+0.5);
		yn = (int)((y-gf.ymin)/gf.step+0.5);
		if(xn < 0 || xn >= gf.nx || yn < 0 || yn >= gf.ny) {
			++gf.inout; // out of grid
		} else {
            ++gf.ingrd; // in grid
            av += z;
            st += z*z;
			nn = xn+gf.nx*yn;
			indxy = gf.index[nn];
			if(indxy > 0) {
				z = (z+gf.field[nn]*indxy)/(1+indxy);
				++gf.indup; // repeated node
			}
			gf.field[nn] = z;  // field value
			gf.index[nn] += 1; // source node
			if(gf.index[nn] > NoValue/2) gf.index[nn] /= 2; // just for safety
		}
	}
	fclose(fi);
    standard(gf.ingrd,&av,&st);
    *stnd = st;
	return 0;
}

static void 
cleanIndex() 
{
    int k;
	gf.ingrd = 0;
	for(k=gf.nx*gf.ny-1; k >= 0; --k) {
		if(gf.index[k]) { ++gf.ingrd; gf.index[k] = 3; } // data
		//else gf.index[k] = 0; // no data
	}
}

static void 
zoom_out(double *field1, short *index1, int nx1, int ny1,	// src large   
		 double *field2, short *index2, int nx2, int ny2) { // => dst small
	int rw[] = { 4,2,1 };
	int jx,jy, ix,iy, mx,my,mm, px,py, pp,qq;
	double t,w,s;

    for(jy=0; jy < ny2; ++jy) { // dst
        my = 3*jy+1;            // src
        if(my >= ny1) continue;
        for(jx=0; jx < nx2; ++jx) { // dst
            mx = 3*jx+1;            // src
            if(mx >= nx1) continue;
			mm = my*nx1+mx;
			if(index1[mm] == 3) {
                qq = nx2*jy+jx;
                field2[qq] = field1[mm];
                index2[qq] = 3; 
                continue;
            }
            s = w = 0.0; // average of input data
            for(ix=-1; ix <= 1; ++ix) {
                px = mx+ix; // src
                if(px < 0 || px >= nx1) continue;
                for(iy=-1; iy <= 1; ++iy) {
                    py = my+iy; // src
                    if(py < 0 || py >= ny1) continue;
					pp = nx1*py+px;
                    if(index1[pp] > 0) {
                        t = rw[abs(ix)+abs(iy)];
                        w += t;
                        s += t * field1[pp];
					}
				}
			}
            if(w > 0.0) {
				qq = nx2*jy+jx;
                field2[qq] = s/w;
                index2[qq] = 1; // added data
			}
		}
	}
}

static void 
zoom_in(double *grid1, short *indx1, int nx1, int ny1,	// src small
		double *grid2, short *indx2, int nx2, int ny2) { // => dst large
	int jx,jy,mx,my, pp;

    for(jx=0; jx < nx1; ++jx) { // from
        mx = 3*jx+1; // to
        if(mx >= nx2) continue;
        for(jy=0; jy < ny1; ++jy) { // from
            my = 3*jy+1; // to
            if(my >= ny2) continue;
			pp = nx2*my+mx;
            if(indx2[pp] > 0) continue;
            indx2[pp] = 1;
            grid2[pp] = grid1[nx1*jy+jx];
		}
	}
}

static int 
fillup1(short *indx, int nx, int ny) 
{ 
    int k, n = 0;
    for(k=nx*ny-1; k >= 0; --k) if(indx[k] > 0) ++n;
    return n;
}

static int 
fillup2(double *fgrid, short *indx, int nx, int ny, int inx, int iny, double wmin) { 
	int rw[] = { 4,2,1 };
	int jx,jy,pp,qq,ix,iy,px,py, kk = 0;
	double t,w,s;
    for(jx=inx; jx < nx; jx+=2) {
        for(jy=iny; jy < ny; jy+=2) {
			pp = nx*jy+jx;
			if(indx[pp] > 0) continue;
			s = w = 0.0;
			for(ix=-1; ix <= 1; ++ix) {
				px = jx+ix;
				if(px < 0 || px >= nx) continue;
				for(iy=-1; iy <= 1; ++iy) {
					py = jy+iy;
					if(py < 0 || py >= ny) continue;
					qq = nx*py+px;
					if(indx[qq] > 0) {
						t = rw[abs(ix)+abs(iy)];
						w += t;
						s += t * fgrid[qq];
					}
				}
            }
            if(w > wmin) {
                fgrid[pp] = s/w;
                indx[pp] = 1; // added node
                ++kk;
			}
		}
	}
	return kk;
}

static void 
fillup3(double *fgrid, short *indx, int nx, int ny) { 
	int kk = 0, kwas;
    double wmin = 3.5;
    while(kk < nx*ny) {
        kwas = kk;  // to decrease direction dependence:
        kk = fillup1(indx, nx, ny);
		kk += fillup2(fgrid, indx, nx, ny, 0, 0, wmin);
		kk += fillup2(fgrid, indx, nx, ny, 1, 1, wmin);
		kk += fillup2(fgrid, indx, nx, ny, 0, 1, wmin);
		kk += fillup2(fgrid, indx, nx, ny, 1, 0, wmin);
        if(kk <= kwas) {
			if(wmin > 1.0) wmin -= 1.0;
			else { // debug: my error
                printf("! kk <= kwas && wmin < 1.0\n");
                exit(1); // ?
            }
		}
	}
}

static void
buildFast() 
{
	int k;
	for(k=1; k < ZIO; ++k) {
		if(gf.nxs[k] < MinGridSize || gf.nys[k] < MinGridSize) break;
        memset((void *)gf.indxs[k],0,gf.nxs[k]*gf.nys[k]*sizeof(short));
		zoom_out(gf.grids[k-1],gf.indxs[k-1],gf.nxs[k-1],gf.nys[k-1],
				 gf.grids[k],  gf.indxs[k],  gf.nxs[k],  gf.nys[k]);
	}
	--k;
	fillup3(gf.grids[k],gf.indxs[k],gf.nxs[k],gf.nys[k]);
	for( ; k > 0; --k) {
		zoom_in(gf.grids[k],  gf.indxs[k],  gf.nxs[k],  gf.nys[k],
				gf.grids[k-1],gf.indxs[k-1],gf.nxs[k-1],gf.nys[k-1]);
		fillup3(gf.grids[k-1],gf.indxs[k-1],gf.nxs[k-1],gf.nys[k-1]);
	}
}

static void 
cosmoothxy(double *fgrid, int nx, int ny, double a, double b) {
	int jx,jy,jj;
    for(jy=0; jy < ny; ++jy) {
        for(jx=1; jx < nx; ++jx) {
			jj = nx*jy+jx;
            fgrid[jj] = a*fgrid[jj]+b*fgrid[jj-1];
            //fgrid[jj] -= b*(fgrid[jj]-fgrid[jj-1]); -- slower
		}
        for(jx=nx-2; jx >= 0; --jx) {
			jj = nx*jy+jx;
            fgrid[jj] = a*fgrid[jj]+b*fgrid[jj+1];
		}
	}
    for(jx=0; jx < nx; ++jx) {
        for(jy=1; jy < ny; ++jy) {
			jj = nx*jy+jx;
            fgrid[jj] = a*fgrid[jj]+b*fgrid[jj-nx];
		}
        for(jy=ny-2; jy >= 0; --jy) {
			jj = nx*jy+jx;
            fgrid[jj] = a*fgrid[jj]+b*fgrid[jj+nx];
		}
	}
}

static double
factor()
{
    double b10 = 0.85, b1000 = 0.96;
    double c = (b1000-b10)/(1/10.-1/1000.);
    double a = b10+c/10;
    return a-c/MIN(gf.nx,gf.ny); // or some average
}

static int
buildGood() {
	int k, nn = gf.nx*gf.ny;
    double vmin = TopLimitValue, vmax = LowLimitValue;
    double bcos, initdiff, curdiff, prevdiff;
	double *ogrid, *ofield; // 'o' for output
	short *ondex;

	ogrid = (double *)malloc(nn*sizeof(double));
	if(ogrid) {
		ofield = (double *)malloc(nn*sizeof(double));
		if(ofield) {
			ondex = (short *)malloc(nn*sizeof(short));
			if(!ondex) { free(ofield); free(ogrid); return -1; }
		} else { free(ogrid); return -1; }
	} else return -1;

	for(k=0; k < nn; ++k) {
		ofield[k] = 0.0;
		ondex[k] = gf.index[k];
		if(gf.field[k] < vmin) vmin = gf.field[k];
		if(gf.field[k] > vmax) vmax = gf.field[k];
	}
    //debug_write("good_init",0,gf.field,gf.nx,gf.ny);
	bcos = factor();
    prevdiff = initdiff = vmax-vmin;
	while (bcos > 0.00001) {
		buildFast(); // gf.field => gf.grids => gf.field
		for(k=0; k < nn; ++k) ogrid[k] = gf.field[k];
		cosmoothxy(ogrid,gf.nx,gf.ny,1-bcos,bcos); // ogrid => ogrid
		vmin = TopLimitValue;
		vmax = LowLimitValue;
		for(k=0; k < nn; ++k) {
			ofield[k] += ogrid[k];
			gf.field[k] -= ogrid[k];
			if(gf.field[k] < vmin) vmin = gf.field[k];
			if(gf.field[k] > vmax) vmax = gf.field[k];
			gf.index[k] = ondex[k];
		}
		curdiff = vmax-vmin;
        if(curdiff <= gf.maydiff) break;
		if(curdiff > 0.999*prevdiff && bcos < 0.001) break; // hard conditions ?
		prevdiff = curdiff;
        bcos *= sqrt(bcos);
	}
	for(k=0; k < nn; ++k)  // copy to output
		gf.field[k] = ofield[k];
	free(ofield);
	free(ogrid);
	free(ondex);
	return 0;
}

static int
buildBest() {
    // In comparison with buildGood() above this function will temporary
    // use enlarged grid (see margin) with smaller cells (lesser step).
    // Also, will apply more accurate using error limits (try to reach 
    // in each datanode).
    return 0;
}

static void
fillDistRadius() 
{
	int nx = gf.nx, ny = gf.ny;
	short *dist = gf.index;
    int kx,ky,kk,k,d;

    for(ky=0; ky < ny; ++ky) {
        for(kx=0; kx < nx; ++kx) {
            kk = nx*ky+kx;
            if(dist[kk]) continue; // only from source nodes
            for(d=2,k=1; kx+k < nx; ++k,d+=2) { // E
				dist[kk+k] = MIN(dist[kk+k],d); }
			for(d=2,k=-1; kx+k >= 0; --k,d+=2) { // W
				dist[kk+k] = MIN(dist[kk+k],d);	}
			for(d=2,k=1; ky+k < ny; ++k,d+=2) { // N
				dist[kk+k*nx] = MIN(dist[kk+k*nx],d); }
			for(d=2,k=-1; ky+k >= 0; --k,d+=2) { // S
				dist[kk+k*nx] = MIN(dist[kk+k*nx],d); }
			for(d=3,k=1; kk+k < nx && ky+k < ny; ++k,d+=3) { // NE
				dist[kk+k*(nx+1)] = MIN(dist[kk+k*(nx+1)],d); }
			for(d=3,k=1; kx-k >= 0 && ky+k < ny; ++k,d+=3) { // NW
				dist[kk+k*(nx-1)] = MIN(dist[kk+k*(nx-1)],d); }
			for(d=3,k=1; kx+k < nx && ky-k >= 0; ++k,d+=3) { // SE
				dist[kk+k*(-nx+1)] = MIN(dist[kk+k*(-nx+1)],d); }
			for(d=3,k=1; kx-k >= 0 && ky-k >= 0; ++k,d+=3) { // SW
				dist[kk+k*(-nx-1)] = MIN(dist[kk+k*(-nx-1)],d); }
		}
	}
}

static void
fillDistNodes() 
{
	int nx = gf.nx, ny = gf.ny;
	short *dist = gf.index;
    int kx,ky,kk, cdist, cdist2, cdist3, ndist;

    for(ndist=2; ndist < NoValue; ++ndist) {
        int news = 0;
        for(ky=0; ky < ny; ++ky) {
			for(kx=0; kx < nx; ++kx) {
                kk = kx+nx*ky;
                if(dist[kk] < ndist) continue;
                cdist = dist[kk];
                cdist2 = cdist+2; cdist3 = cdist+3;
                if(0 < kx) {
                    if(dist[kk-1] > cdist2) { dist[kk-1] = cdist2; ++news; }
                    if(ky > 0 && dist[kk-1-nx] > cdist3) {
                        dist[kk-1-nx] = cdist3; ++news; }
                    if(ky < ny-1 && dist[kk-1+nx] > cdist3) {
                        dist[kk-1+nx] = cdist3; ++news; }
                }
                if(kx+1 < nx) {
                    if(dist[kk+1] > cdist2) { dist[kk+1] = cdist2; ++news; }
                    if(ky > 0 && dist[kk+1-nx] > cdist3) {
                        dist[kk+1-nx] = cdist3; ++news; }
                    if(ky < ny-1 && dist[kk+1+nx] > cdist3) {
                        dist[kk+1+nx] = cdist3; ++news; }
                }
                if(ky > 0 && dist[kk-nx] > cdist2) {
                    dist[kk-nx] = cdist2; ++news; }
                if(ky < ny-1 && dist[kk+nx] > cdist2) {
                    dist[kk+nx] = cdist2; ++news; }
            }
        }
        if(!news) break;
    }
}

static void
trimDistances() 
{
	int nx = gf.nx, ny = gf.ny;
	short *dist = gf.index;
	int kx,ky,kk, news=0;
	int xdist = (int)(2.001*gf.trimdist/gf.step+0.5);
	int nearn[8] = { 1, 1+nx, nx, nx-1, -1, -1-nx, -nx, 1-nx };
    int n, nears = 8; // == sizeof(nearn)/sizeof(int)

	for(kx=0; kx < nx; ++kx) {	// most exterior layer
		kk = kx;
		if(dist[kk] >= xdist) { dist[kk] = NoValue; ++news; }
		kk = nx*(ny-1)+kx;
		if(dist[kk] >= xdist) { dist[kk] = NoValue; ++news; }
	}
	for(ky=0; ky < ny; ++ky) {
		kk = nx*ky;
		if(dist[kk] >= xdist) { dist[kk] = NoValue; ++news; }
		kk = nx*ky+(nx-1);
		if(dist[kk] >= xdist) { dist[kk] = NoValue; ++news; }
	}
    
	while(news) {
		news = 0;
		for(ky=1; ky < ny-1; ++ky) {
			for(kx=1; kx < nx-1; ++kx) {
				kk = nx*ky+kx;
				if(dist[kk] < xdist || dist[kk] == NoValue) continue;
                for(n=0; n < nears; ++n) {
                    if(dist[kk+nearn[n]] == NoValue) {
                        dist[kk] = NoValue; ++news; break;
                    }
                }
			}
		}
	}
}

static void 
useDistances() 
{
    int n;
    short *dist = gf.index;
    for(n=gf.nx*gf.ny-1; n >= 0; --n) {
        if(gf.index[n] > 0) gf.index[n] = 0;
        else gf.index[n] = NoValue;
    }
	fillDistRadius();
	fillDistNodes();
	trimDistances();
	for(n = gf.nx*gf.ny-1; n >= 0; --n) {
		if(dist[n]==NoValue) {
            gf.field[n] = UnValue;
        }
	}
}

static int
goodgeometry()
{
	if(gf.nx < MinGridSize || gf.ny < MinGridSize) {
		PyErr_SetString(PyExc_ValueError,"Grid size in nodes too small");
		return 0;
	}
	if(gf.nx > MaxGridSize || gf.ny > MaxGridSize) {
		PyErr_SetString(PyExc_ValueError,"Grid size in nodes too large");
		return 0;
	}
	if(gf.step <= MinStep) {
		PyErr_SetString(PyExc_ValueError,"Cell size too small");
		return 0;
	}
    return 1;
}


static PyObject *
gridresult()
{
    int r,c, rv;
    PyObject *done = PyList_New(gf.nx*gf.ny);
	if(!done) {
		PyErr_SetString(PyExc_ValueError,"Failure create output list");
		return done;
	}
	for(r=0; r < gf.ny; ++r) {
		for(c=0; c < gf.nx; ++c) {
			double v = gf.field[r*gf.nx+c];
	        if(isUnvalue(v)) v = gf.unvalue; // caller unvalue
	        {
                PyObject * value = PyFloat_FromDouble(v);
                PyList_SetItem(done,r*gf.nx+c,value);
                //debug:
                //PyObject * value = PyFloat_FromDouble(v);
                //if(!value) printf("PyFloat_FromDouble %d %d %f\n",r,c,v);
                //rv = PyList_SetItem(done,r*gf.nx+c,value);
                //if(rv) { printf("list setitem %d\n",r); }
            }
        }
	}
    freeGridFieldMemory(7); // free gf.field
	return done;
}
	
//--- (Registered function.) ---------------------------
static PyObject *
fromfile(PyObject * self, PyObject * args, PyObject *kwds)
{
	static char *kwlist[] = {
		"xyzfile",		// input file, x,y,z always
		"nx","ny",      // grid sizes in nodes
		"step",         // cell sizes
		"xmin","ymin",  // low bottom
		"unvalue",      // caller unvalue
		"method",       // good/best 
		"trimming",     // distance to trim grid
        "abserror",     // acceptable error, given
        "relerror",     // the same, part of standard (if abserror==0)
        NULL
    };
    int rv;
    double st; // for tolerance
	char *xyzfile, *method="";
    double abserr=0.0, relerr=0.0;
    PyObject *grid;

	freeGridFieldMemory(7); // for the case of a new call
	initGridField();       // set default values
	if(!PyArg_ParseTupleAndKeywords(args, kwds,
			"siiddd|dsddd", kwlist,
			&xyzfile, &gf.nx, &gf.ny, &gf.step, &gf.xmin, &gf.ymin,
			&gf.unvalue, &method, &gf.trimdist, &abserr, &relerr)) {
        return NULL;
	}
	// Check main input parameters:
    if(!goodgeometry()) return NULL;
	if(method[0] == 'B') gf.method = 'B'; // Best/Good
    else gf.method = 'G';

	rv = initGridFieldMemory(); // also initialize
	if(rv) { PyErr_NoMemory(); goto badexit; }

	// read input XYZ data, also put to dist grid
	rv = getInputFileData(xyzfile,&st);
	if(rv < 0) {
		PyErr_SetString(PyExc_IOError,"Could not open input file");
		goto badexit;
	}
	if(gf.ingrd < 1) {
		PyErr_SetString(PyExc_IOError,"No input data");
		goto badexit;
	}
    if(abserr > 0.0) gf.maydiff = MIN(abserr,st/2);
    else             gf.maydiff = MIN(relerr*st,st/2);

	cleanIndex();  // mark nodes 0/3; set gf.innodes
	rv = buildGood();
	if(rv) { PyErr_NoMemory(); goto badexit; }
	if(gf.method=='B')
		buildBest(); // not implemented
	freeGridFieldMemory(4); // free gf.aux grids

	// distances
	if(gf.trimdist >= gf.step)
		useDistances();
	freeGridFieldMemory(2); // free gf.index

    // Create resulting list
    grid = gridresult();
    if(grid) 
        return Py_BuildValue("O",grid);

badexit:
	freeGridFieldMemory(7);
	return NULL;
}

static void
getInputXYZdata(PyObject * xfast, PyObject * yfast, PyObject * zfast, double * stnd) 
{
	int k, nn, indxy, xn,yn, nvals = PyObject_Length(xfast);
	double x,y,z, av,st;
    av = st = 0.0;

	for(k=0; k < nvals; ++k) {
        PyObject *xf = PyNumber_Float(PySequence_GetItem(xfast,k));
        PyObject *yf = PyNumber_Float(PySequence_GetItem(yfast,k));
        PyObject *zf = PyNumber_Float(PySequence_GetItem(zfast,k));
        if(!xf || !yf || !zf) {
            //PyErr_SetString(PyExc_ValueError,"Not a Float Value in xdata/ydata/zdata");
            //return;
            continue;
        }
        x = PyFloat_AsDouble(xf);
        y = PyFloat_AsDouble(yf);
        z = PyFloat_AsDouble(zf);
		++gf.indat; // number of input data points
		xn = (int)((x-gf.xmin)/gf.step+0.5);
		yn = (int)((y-gf.ymin)/gf.step+0.5);
		if(xn < 0 || xn >= gf.nx || yn < 0 || yn >= gf.ny) {
			++gf.inout; // out of grid
		} else {
            ++gf.ingrd;
            av += z;
            st += z*z;
			nn = xn+gf.nx*yn;
			indxy = gf.index[nn];
			if(indxy > 0) {
				z = (z+gf.field[nn]*indxy)/(1+indxy);
				++gf.indup; // repeated node
			}
			gf.field[nn] = z;  // field value
			gf.index[nn] += 1; // source node
			if(gf.index[nn] > 16300) gf.index[nn] /= 2; // just for safety
		}
	}
    standard(gf.ingrd,&av,&st);
    *stnd = st;
}


//--- (Registered function.) ---------------------------
static PyObject *
fromxyz(PyObject * self, PyObject * args, PyObject *kwds)
{
	static char *kwlist[] = {
		"xdata","ydata","zdata", // input lists/tuples/... x,y,z
        // other parameters are the same as in fromfile():
		"nx","ny",      // grid sizes in nodes
		"step",         // cell sizes
		"xmin","ymin",  // low bottom
		"unvalue",      // caller unvalue
		"method",       // good/best 
		"trimming",     // distance to trim grid
        "abserror",     // acceptable error, given
        "relerror",     // the same, part of standard (if abserror==0)
        NULL
    };
    int rv;
    double st; // for tolerance
	char *method="";
    PyObject *xdata, *ydata, *zdata, *xfast, *yfast, *zfast;
    PyObject *grid;
    double abserr=0.0, relerr=0.0;

	freeGridFieldMemory(7); // for the case of a new call
	initGridField();       // set default values
	if(!PyArg_ParseTupleAndKeywords(args, kwds,
			"OOOiiddd|dsddd", kwlist,
			&xdata, &ydata, &zdata, &gf.nx, &gf.ny, &gf.step, &gf.xmin, &gf.ymin,
			&gf.unvalue, &method, &gf.trimdist, &abserr, &relerr)) {
        return NULL;
	}
	// Check main input parameters:
    if(!goodgeometry()) return NULL;
	if(method[0] == 'B') gf.method = 'B'; // Best/Good
    else gf.method = 'G';

   	xfast = PySequence_Fast(xdata,"xdata must be sequence"); // TypeError
   	yfast = PySequence_Fast(ydata,"ydata must be sequence"); // TypeError
   	zfast = PySequence_Fast(zdata,"zdata must be sequence"); // TypeError
	if(!xfast || !yfast || !zfast) return NULL;

    if(PyObject_Length(xfast) != PyObject_Length(yfast) ||
        PyObject_Length(xfast) != PyObject_Length(zfast)) {
        PyErr_SetString(PyExc_ValueError,"x/y/z data of different sizes");
		return NULL;
	}

	rv = initGridFieldMemory(); // also initialize
	if(rv) { PyErr_NoMemory(); goto badexit; }

	// read input XYZ data, also put to dist grid
	getInputXYZdata(xfast,yfast,zfast,&st);
	if(gf.ingrd < 1) {
		PyErr_SetString(PyExc_IOError,"No input data");
		goto badexit;
	}
    if(abserr > 0.0) gf.maydiff = MIN(abserr,st/2);
    else             gf.maydiff = MIN(relerr*st,st/2);

	cleanIndex();  // mark nodes 0/3; set gf.innodes
	rv = buildGood();
	if(rv) { PyErr_NoMemory(); goto badexit; }
	if(gf.method=='B')
		buildBest(); // not implemented
	freeGridFieldMemory(4); // free gf.aux grids

	// distances
	if(gf.trimdist >= gf.step)
		useDistances();
	freeGridFieldMemory(2); // free gf.index

    // Create resulting list
    grid = gridresult();
    if(grid) 
        return Py_BuildValue("O",grid);

badexit:
	freeGridFieldMemory(7);
	return NULL;
}

static int
makeGXF(char *filename, PyObject *fast, int nx, int ny, 
		double xmin, double ymin, double step, double unvalue)
{
    return 0;
}

static int
makeGRD(char *filename, PyObject *fast, int nx, int ny, 
		double xmin, double ymin, double step, double unvalue)
{
    return 0;
}

static int
makeXYZ(char *filename, PyObject *fast, int nx, int ny, 
		double xmin, double ymin, double step, double unvalue,
		char *header, char *format)
{
	int rv, kx,ky;
	double uneps = 0.01*fabs(unvalue); // hope: unvalue != unvalue+uneps
	char xyzformat[50];
	FILE *fo = fopen(filename,"wt");
	if(!fo) {
		PyErr_SetString(PyExc_ValueError,"failure open output file");
		return -1;
	}
	if(strlen(header) > 0) {
		rv = fprintf(fo,"%s\n",header);
		if(rv < 1) { fclose(fo); return -2; }
	}

	if(strlen(format) > 8) strcpy(xyzformat,format); // at least "%f %f %f\n"
	else strcpy(xyzformat,"%.2f %.2f %.3f\n");
	for(ky=0; ky < ny; ++ky) {
		for(kx=0; kx < nx; ++kx) {
			PyObject *value = PySequence_Fast_GET_ITEM(fast,ky*nx+kx);
			double v = PyFloat_AsDouble(value);
			if(fabs(v-unvalue) <= uneps) continue;
			rv = fprintf(fo,xyzformat,xmin+kx*step,ymin+ky*step,v);
			if(rv < 1) { fclose(fo); return -2; }
		}
	}
	fclose(fo);
	return 0;
}
	
//--- (Registered function.) ---------------------------
static PyObject *
tofile(PyObject * self, PyObject * args, PyObject *kwds)
{
	static char *kwlist[] = {
		"filename",		// output file name
		"gridtype",		// output grid type: "xyz"/"gxf"/...
		"griddata",     // list or similar of values in nodes
		"nx","ny",      // grid sizes in nodes
		"step",         // cell sizes
		"xmin","ymin",  // low bottom
		"unvalue",      // caller unvalue to skip in output
		"header",       // text string
		"format",       // output line format
        NULL
    };
	char *filename, *gridtype, *header="", *format="";
	PyObject *grid, *fast;
	int nx=0, ny=0, rv;
	double step=1.0, xmin=0.0, ymin=0.0, unvalue=UnValue;

	if(!PyArg_ParseTupleAndKeywords(args, kwds,
			"ssOii|ddddss", kwlist,
			&filename, &gridtype, &grid, &nx, &ny,
			&step, &xmin, &ymin, &unvalue, &header, &format)) {
		return NULL;
	}
	// Check main input parameters:
	if(nx < 1 || ny < 1) {
		PyErr_SetString(PyExc_ValueError,"No grid size");
		return NULL;
	}
	fast = PySequence_Fast(grid,"griddata must be sequence"); // TypeError
	if(!fast) return NULL;

	if(!strncmp(gridtype,"gxf",3)) {
        // I mean:
        // http://www.geosoft.com/resources/technotes/pdfs/gxfr3d9_1.pdf
        // though there are others, e.g.
        // http://www.euronav.co.uk/Downloads/GXF%20spec/GXF_Spec_123.pdf
		rv = makeGXF(filename,fast,nx,ny,xmin,ymin,step,unvalue);
	} else
	if(!strncmp(gridtype,"grd",3)) {
        // There are a lot of them. I mean one or two of Surfer's:
        // http://www.goldensoftware.com/
		rv = makeGRD(filename,fast,nx,ny,xmin,ymin,step,unvalue);
	} else
	if(!strncmp(gridtype,"xyz",3)) {
        // Just text lines with triples:  x y z
		rv = makeXYZ(filename,fast,nx,ny,xmin,ymin,step,unvalue,header,format);
	} else {
        // Another formats ?
		PyErr_SetString(PyExc_ValueError,"unsupported output grid type");
		return NULL;
	}
	return Py_BuildValue("i",rv);
}

//--- registered methods -------------------------------------
static PyMethodDef buildmethods[] = {
	{ "fromfile", (PyCFunction)fromfile,
				METH_VARARGS | METH_KEYWORDS,	" build field from xyz file" },
	{ "fromxyz", (PyCFunction)fromxyz,
				METH_VARARGS | METH_KEYWORDS,	" build field from xyz data" },
	//{ "frominfo", (PyCFunction)fromafter,
	//			METH_VARARGS | METH_KEYWORDS,	" info after building" },
	{ "tofile", (PyCFunction)tofile,
				METH_VARARGS | METH_KEYWORDS,	" write field to xyz/gxf/... file" },
	{ "filestatistics", (PyCFunction)filestatistics,
				METH_VARARGS | METH_KEYWORDS,	" find file data limits" },
	{ NULL, NULL, 0, NULL }
};

//--- registration --------------------------------------------
PyMODINIT_FUNC
initbuild_grid(void) {
	Py_InitModule("build_grid", buildmethods);
	initGridField();	// initialize global 'gf'
}

