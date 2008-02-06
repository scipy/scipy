#include<stdio.h>
#include<stdlib.h>

float trilinear_A(unsigned char *pVolume, int x, int y, int z, float dx, float dy, float dz, int dims[]){

	// Vxyz for [0,1] values of x, y, z
	int V000;
	int V100;
	int V010;
	int V001;
	int V011;
	int V101;
	int V110;
	int V111;

	int ptr_x0;
	int ptr_y0;
	int ptr_z0;

	int ptr_x1;
	int ptr_y1;
	int ptr_z1;

	float valueXYZ;

	ptr_x0 = x;
	ptr_y0 = y * dims[0];
	ptr_z0 = z * dims[1];

	ptr_x1 = ptr_x0 + 1;
	ptr_y1 = ptr_y0 + dims[0];
	ptr_z1 = ptr_z0 + dims[1];

	V000 = pVolume[ptr_x0+ptr_y0+ptr_z0];
	V100 = pVolume[ptr_x1+ptr_y0+ptr_z0];
	V010 = pVolume[ptr_x0+ptr_y1+ptr_z0];
	V001 = pVolume[ptr_x0+ptr_y0+ptr_z1];
	V011 = pVolume[ptr_x0+ptr_y1+ptr_z1];
	V101 = pVolume[ptr_x1+ptr_y0+ptr_z1];
	V110 = pVolume[ptr_x1+ptr_y1+ptr_z0];
	V111 = pVolume[ptr_x1+ptr_y1+ptr_z1];

	// dx, dy, dz are increments in x, y, z
	// dx = 0 is x = 1 as x, y and z are [0, 1] in range

	valueXYZ = 
		V000 * (1.0-dx) * (1.0 - dy) * (1.0 - dz) +
		V100 * (dx)     * (1.0 - dy) * (1.0 - dz) +
		V010 * (1.0-dx) * (dy)       * (1.0 - dz) +
		V001 * (1.0-dx) * (1.0 - dy) * (dz)       +
		V101 * (dx)     * (1.0 - dy) * (dz)       +
		V011 * (1.0-dx) * (dy)       * (dz)       +
		V110 * (dx)     * (dy)       * (1.0 - dz) +
		V111 * (dx)     * (dy)       * (dz);


	return(valueXYZ);

}

float trilinear_B(unsigned char *pVolume, float dx, float dy, float dz, int corners[]){

	// Vxyz for [0,1] values of x, y, z
	int V000;
	int V100;
	int V010;
	int V001;
	int V011;
	int V101;
	int V110;
	int V111;

	int ptr_x0 = corners[0];
	int ptr_y0 = corners[1];
	int ptr_z0 = corners[2];

	int ptr_x1 = corners[3];
	int ptr_y1 = corners[4];
	int ptr_z1 = corners[5];

	float valueXYZ;

	V000 = pVolume[ptr_x0+ptr_y0+ptr_z0];
	V100 = pVolume[ptr_x1+ptr_y0+ptr_z0];
	V010 = pVolume[ptr_x0+ptr_y1+ptr_z0];
	V001 = pVolume[ptr_x0+ptr_y0+ptr_z1];
	V011 = pVolume[ptr_x0+ptr_y1+ptr_z1];
	V101 = pVolume[ptr_x1+ptr_y0+ptr_z1];
	V110 = pVolume[ptr_x1+ptr_y1+ptr_z0];
	V111 = pVolume[ptr_x1+ptr_y1+ptr_z1];

	// dx, dy, dz are increments in x, y, z
	// dx = 0 is x = 1 as x, y and z are [0, 1] in range

	valueXYZ = 
		V000 * (1.0-dx) * (1.0 - dy) * (1.0 - dz) +
		V100 * (dx)     * (1.0 - dy) * (1.0 - dz) +
		V010 * (1.0-dx) * (dy)       * (1.0 - dz) +
		V001 * (1.0-dx) * (1.0 - dy) * (dz)       +
		V101 * (dx)     * (1.0 - dy) * (dz)       +
		V011 * (1.0-dx) * (dy)       * (dz)       +
		V110 * (dx)     * (dy)       * (1.0 - dz) +
		V111 * (dx)     * (dy)       * (dz);


	return(valueXYZ);

}

int NI_Histogram2D(int layersF, int rowsF, int colsF, int layersG, int rowsG, int colsG, 
		   int *dimSteps, double *M, unsigned char *imageG, unsigned char *imageF, double *H)
{

	int status;
	int seed;
	int dimsF[3];
	int dimsG[3];
	int dims_F[2];
	int dims_G[2];
	int ivf, ivg;
	float ran_x, ran_y, ran_z;
	float vf, delta;
	float x, y, z;
	float dx, dy, dz;
	float xp, yp, zp;
	float rx, ry, rz;

	dimsF[0] = colsF;
	dimsF[1] = rowsF;
	dimsF[2] = layersF;
	dimsG[0] = colsG;
	dimsG[1] = rowsG;
	dimsG[2] = layersG;

	dims_G[0] = dimsG[0];
	dims_G[1] = dimsG[0]*dimsG[1];
	dims_F[0] = dimsF[0];
	dims_F[1] = dimsF[0]*dimsF[1];

	seed = 1000;
	srand(seed);

	/* because of stochastic sampling, subtract 1 from upper bounds */
	for(z = 0.0; z < layersG-dimSteps[2]-1; z += dimSteps[2]){
	    for(y = 0.0; y < rowsG-dimSteps[1]-1; y += dimSteps[1]){
	        for(x = 0.0; x < colsG-dimSteps[0]-1; x += dimSteps[0]){
		    /* positive jitter the x, y, z values */
		    ran_x = 1.0 * rand()/((float)RAND_MAX);
		    ran_y = 1.0 * rand()/((float)RAND_MAX);
		    ran_z = 1.0 * rand()/((float)RAND_MAX);
		    dx = x + ran_x*dimSteps[0];
		    dy = y + ran_y*dimSteps[1];
		    dz = z + ran_z*dimSteps[2];

		    /* get the 'from' coordinates */
		    xp = M[0]*dx + M[1]*dy + M[2]*dz  + M[3];
		    yp = M[4]*dx + M[5]*dy + M[6]*dz  + M[7];
		    zp = M[8]*dx + M[9]*dy + M[10]*dz + M[11];
		    /* clip the resample window */
		    if((zp >= 0.0 && zp < layersF-dimSteps[2]) && 
		       (yp >= 0.0 && yp < rowsF-dimSteps[1])   && 
		       (xp >= 0.0 && xp < colsF-dimSteps[0])){
		        /* resample the coordinates using a trilinear interpolation */
			/* resample imageF using the rotated-jittered xyz coordinates */
			rx = xp - (int)xp; 
			ry = yp - (int)yp; 
			rz = zp - (int)zp; 
			//vf = trilinear_A(imageF, (int)dx, (int)dy, (int)dz, rx, ry, rz, dims_F);
			vf = trilinear_A(imageF, (int)xp, (int)yp, (int)zp, rx, ry, rz, dims_F);
			/* floor */
			ivf = (int)vf;
			delta = vf - ivf;
			/* resample imageG using the jittered xyz coordinates */
			rx = dx - (int)dx; 
			ry = dy - (int)dy; 
			rz = dz - (int)dz; 
			ivg = (int)trilinear_A(imageG, (int)dx, (int)dy, (int)dz, rx, ry, rz, dims_G);
			//ivg = (int)trilinear_A(imageG, (int)xp, (int)yp, (int)zp, rx, ry, rz, dims_G);
			/* ivf will be < 255 as 8 bit data and trilinear doesn't ring */
			H[ivf+256*ivg] += 1.0 - delta;
			if(ivf < 255){
			    H[ivf+1+256*ivg] += delta;
			}
		    }
	        }
	    }
	}

	status = 1;

	return status;

}


int NI_Histogram2DLite(int layersF, int rowsF, int colsF, int layersG, int rowsG, int colsG,
		       int *dimSteps, double *M, unsigned char *imageG, unsigned char *imageF, double *H)
{

	int i;
	int status;
	int sliceG;
	int rowG;
	int sliceSizeG;
	int dimsF[3];
	int dimsG[3];
	int dims[2];
	int ivf, ivg;
	float vf, delta;
	float x, y, z;
	float xp, yp, zp;
	float dx, dy, dz;

	int ptr_x0;
	int ptr_y0;
	int ptr_z0;
	int ptr_x1;
	int ptr_y1;
	int ptr_z1;
	//
	// Vxyz for [0,1] values of x, y, z
	//
	int V000;
	int V100;
	int V010;
	int V001;
	int V011;
	int V101;
	int V110;
	int V111;
	float valueXYZ;

	//
	// G is fixed; F is rotated
	//
	sliceSizeG = rowsG * colsG;
	dimsF[0] = colsF;
	dimsF[1] = rowsF;
	dimsF[2] = layersF;
	dimsG[0] = colsG;
	dimsG[1] = rowsG;
	dimsG[2] = layersG;

	dims[0] = dimsF[0];
	dims[1] = dimsF[0]*dimsF[1];

	for(z = 0.0; z < layersG-dimSteps[2]-1; z += dimSteps[2]){
	    sliceG = (int)z * sliceSizeG;
	    for(y = 0.0; y < rowsG-dimSteps[1]-1; y += dimSteps[1]){
		rowG = (int)y * colsG;
	        for(x = 0.0; x < colsG-dimSteps[0]-1; x += dimSteps[0]){
		    // get the 'from' coordinates 
		    xp = M[0]*x + M[1]*y + M[2]*z  + M[3];
		    yp = M[4]*x + M[5]*y + M[6]*z  + M[7];
		    zp = M[8]*x + M[9]*y + M[10]*z + M[11];
		    // clip the resample window 
		    if((zp >= 0.0 && zp < layersF-dimSteps[2]) && 
		       (yp >= 0.0 && yp < rowsF-dimSteps[1]) && 
		       (xp >= 0.0 && xp < colsF-dimSteps[0])){

			// corners of the 3D unit volume cube
	    		ptr_z0 = (int)zp * dims[1];
	    		ptr_z1 = ptr_z0 + dims[1];
			ptr_y0 = (int)yp * dims[0];
			ptr_y1 = ptr_y0 + dims[0];
		    	ptr_x0 = (int)xp;
		    	ptr_x1 = ptr_x0 + 1;
			dx = xp - (int)xp; 
			dy = yp - (int)yp; 
			dz = zp - (int)zp; 

			// imageG is not rotated. sample the given x,y,z
			ivg = imageG[sliceG+rowG+(int)x];
			// imageF IS rotated. sample the rotated xp,yp,zp
			V000 = imageF[ptr_x0+ptr_y0+ptr_z0];
			V100 = imageF[ptr_x1+ptr_y0+ptr_z0];
			V010 = imageF[ptr_x0+ptr_y1+ptr_z0];
			V001 = imageF[ptr_x0+ptr_y0+ptr_z1];
			V011 = imageF[ptr_x0+ptr_y1+ptr_z1];
			V101 = imageF[ptr_x1+ptr_y0+ptr_z1];
			V110 = imageF[ptr_x1+ptr_y1+ptr_z0];
			V111 = imageF[ptr_x1+ptr_y1+ptr_z1];
			
		        vf = V000 * (1.0-dx) * (1.0 - dy) * (1.0 - dz) +
			     V100 * (dx)     * (1.0 - dy) * (1.0 - dz) +
			     V010 * (1.0-dx) * (dy)       * (1.0 - dz) +
			     V001 * (1.0-dx) * (1.0 - dy) * (dz)       +
			     V101 * (dx)     * (1.0 - dy) * (dz)       +
			     V011 * (1.0-dx) * (dy)       * (dz)       +
			     V110 * (dx)     * (dy)       * (1.0 - dz) +
			     V111 * (dx)     * (dy)       * (dz);

			ivf = (int)(vf);
			H[ivf+256*ivg] += 1.0;
		    }
	        }
	    }
	}

	status = 1;

	return status;

}


int NI_LinearResample(int layersF, int rowsF, int colsF, int layersG, int rowsG, int colsG,
		       int *dimSteps, double *M, unsigned char *imageG, unsigned char *imageF)
{

	int i;
	int status;
	int sliceG;
	int rowG;
	int sliceSizeG;
	int dimsF[3];
	int dimsG[3];
	int dims[2];
	int ivf, ivg;
	float vf, delta;
	float x, y, z;
	float xp, yp, zp;
	float dx, dy, dz;

	int ptr_x0;
	int ptr_y0;
	int ptr_z0;
	int ptr_x1;
	int ptr_y1;
	int ptr_z1;
	//
	// Vxyz for [0,1] values of x, y, z
	//
	int V000;
	int V100;
	int V010;
	int V001;
	int V011;
	int V101;
	int V110;
	int V111;
	float valueXYZ;

	//
	// G is fixed; F is rotated
	//
	sliceSizeG = rowsG * colsG;
	dimsF[0] = colsF;
	dimsF[1] = rowsF;
	dimsF[2] = layersF;
	dimsG[0] = colsG;
	dimsG[1] = rowsG;
	dimsG[2] = layersG;

	dims[0] = dimsF[0];
	dims[1] = dimsF[0]*dimsF[1];

	for(z = 0.0; z < layersG-dimSteps[2]-1; z += dimSteps[2]){
	    sliceG = (int)z * sliceSizeG;
	    for(y = 0.0; y < rowsG-dimSteps[1]-1; y += dimSteps[1]){
		rowG = (int)y * colsG;
	        for(x = 0.0; x < colsG-dimSteps[0]-1; x += dimSteps[0]){
		    // get the 'from' coordinates 
		    xp = M[0]*x + M[1]*y + M[2]*z  + M[3];
		    yp = M[4]*x + M[5]*y + M[6]*z  + M[7];
		    zp = M[8]*x + M[9]*y + M[10]*z + M[11];
		    // clip the resample window 
		    if((zp >= 0.0 && zp < layersF-dimSteps[2]) && 
		       (yp >= 0.0 && yp < rowsF-dimSteps[1]) && 
		       (xp >= 0.0 && xp < colsF-dimSteps[0])){

			// corners of the 3D unit volume cube
	    		ptr_z0 = (int)zp * dims[1];
	    		ptr_z1 = ptr_z0 + dims[1];
			ptr_y0 = (int)yp * dims[0];
			ptr_y1 = ptr_y0 + dims[0];
		    	ptr_x0 = (int)xp;
		    	ptr_x1 = ptr_x0 + 1;
			dx = xp - (int)xp; 
			dy = yp - (int)yp; 
			dz = zp - (int)zp; 

			// imageF IS rotated. sample the rotated xp,yp,zp
			// and stored in imageG
			V000 = imageF[ptr_x0+ptr_y0+ptr_z0];
			V100 = imageF[ptr_x1+ptr_y0+ptr_z0];
			V010 = imageF[ptr_x0+ptr_y1+ptr_z0];
			V001 = imageF[ptr_x0+ptr_y0+ptr_z1];
			V011 = imageF[ptr_x0+ptr_y1+ptr_z1];
			V101 = imageF[ptr_x1+ptr_y0+ptr_z1];
			V110 = imageF[ptr_x1+ptr_y1+ptr_z0];
			V111 = imageF[ptr_x1+ptr_y1+ptr_z1];
			
		        vf = V000 * (1.0-dx) * (1.0 - dy) * (1.0 - dz) +
			     V100 * (dx)     * (1.0 - dy) * (1.0 - dz) +
			     V010 * (1.0-dx) * (dy)       * (1.0 - dz) +
			     V001 * (1.0-dx) * (1.0 - dy) * (dz)       +
			     V101 * (dx)     * (1.0 - dy) * (dz)       +
			     V011 * (1.0-dx) * (dy)       * (dz)       +
			     V110 * (dx)     * (dy)       * (1.0 - dz) +
			     V111 * (dx)     * (dy)       * (dz);

			imageG[sliceG+rowG+(int)x] = (int)vf;

		    }
	        }
	    }
	}

	status = 1;

	return status;

}


