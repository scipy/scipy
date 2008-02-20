#include<stdio.h>
#include<stdlib.h>

float tri_cubic_convolve(unsigned char *pVolume, int x, int y, int z, float xp, float yp,
	                 float zp, int colsG, int rowsG, int layersG, int sliceSizeG){

	int i, j, k;
	int layerOffsets[4];
	int rowOffsets[4];
	float ps1, ps2, ps3;
	float Y[4], NewRow[4], NewLayer[4];
	float R, C, L, D, T;
	float valueXYZ = 0.0;
	float dataCube[4][4][4];
	/*            [cols][rows][layers] */

	rowOffsets[0]   = (y-1)*colsG;
	rowOffsets[1]   = (y  )*colsG;
	rowOffsets[2]   = (y+1)*colsG;
	rowOffsets[3]   = (y+2)*colsG;

	layerOffsets[0] = (z-1)*sliceSizeG;
	layerOffsets[1] = (z  )*sliceSizeG;
	layerOffsets[2] = (z+1)*sliceSizeG;
	layerOffsets[3] = (z+2)*sliceSizeG;

	/* get numerator for interpolation */
	C = xp - (float)x;
	R = yp - (float)y;
	L = zp - (float)z;
	D = (float)0.002;

	/* get 4x4 window over all 4 layers */
	for(i = 0; i < 4; ++i){
	    for(j = 0; j < 4; ++j){
		dataCube[0][j][i] = (float)pVolume[layerOffsets[i]+rowOffsets[j]+x-1];
		dataCube[1][j][i] = (float)pVolume[layerOffsets[i]+rowOffsets[j]+x];
		dataCube[2][j][i] = (float)pVolume[layerOffsets[i]+rowOffsets[j]+x+1];
		dataCube[3][j][i] = (float)pVolume[layerOffsets[i]+rowOffsets[j]+x+2];
	    }
	}

	for(i = 0; i < 4; ++i){
	    /* interpolate 4 rows in all 4 layers */
	    for(j = 0; j < 4; ++j){
		if(C > D){
		    Y[0] = dataCube[0][j][i];
		    Y[1] = dataCube[1][j][i];
		    Y[2] = dataCube[2][j][i];
		    Y[3] = dataCube[3][j][i];
		    ps1       = Y[2] - Y[0];
		    ps2       = (float)2.0*(Y[0] - Y[1]) + Y[2] - Y[3];
		    ps3       = -Y[0] + Y[1] - Y[2] + Y[3];
		    NewRow[j] = Y[1]+C*(ps1+C*(ps2+C*ps3));
		}
		else{
		    NewRow[j] = dataCube[1][j][i];
		}
	    }
	    /* interpolate across 4 columns */
	    if(R > D){
		Y[0] = NewRow[0];
		Y[1] = NewRow[1];
		Y[2] = NewRow[2];
		Y[3] = NewRow[3];
		ps1  = Y[2] - Y[0];
		ps2  = (float)2.0*(Y[0] - Y[1]) + Y[2] - Y[3];
		ps3  = -Y[0] + Y[1] - Y[2] + Y[3];
		T    = (Y[1]+R*(ps1+R*(ps2+R*ps3)));
		NewLayer[i] = T;
	    }
	    else{
		T = NewRow[1];
		NewLayer[i] = T;
	    } 
	}
	/* interpolate across 4 layers */
	if(L > D){
	    Y[0] = NewLayer[0];
	    Y[1] = NewLayer[1];
	    Y[2] = NewLayer[2];
	    Y[3] = NewLayer[3];
	    ps1  = Y[2] - Y[0];
	    ps2  = (float)2.0*(Y[0] - Y[1]) + Y[2] - Y[3];
	    ps3  = -Y[0] + Y[1] - Y[2] + Y[3];
	    T    = (Y[1]+L*(ps1+L*(ps2+L*ps3)));
	    valueXYZ = T;
	}
	else{
	    T = NewLayer[1];
	    valueXYZ = T;
	} 

	return(valueXYZ);

}

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



int NI_CubicResample(int layersF, int rowsF, int colsF, int layersG, int rowsG, int colsG,
	             int *dimSteps, double *M, unsigned char *imageG, unsigned char *imageF)
{

	int i;
	int status;
	int sliceG;
	int rowG;
	int sliceSizeG;
	int ivf;
	float vf;
	float x, y, z;
	float xp, yp, zp;

	sliceSizeG = rowsG * colsG;
	for(z = 1.0; z < layersG-dimSteps[2]-2; z += dimSteps[2]){
	    sliceG = (int)z * sliceSizeG;
	    for(y = 1.0; y < rowsG-dimSteps[1]-2; y += dimSteps[1]){
		rowG = (int)y * colsG;
	        for(x = 1.0; x < colsG-dimSteps[0]-2; x += dimSteps[0]){
		    // get the 'from' coordinates 
		    xp = M[0]*x + M[1]*y + M[2]*z  + M[3];
		    yp = M[4]*x + M[5]*y + M[6]*z  + M[7];
		    zp = M[8]*x + M[9]*y + M[10]*z + M[11];
		    // clip the resample window 
		    if((zp >= 1.0 && zp < layersF-dimSteps[2]-2) && 
		       (yp >= 1.0 && yp < rowsF-dimSteps[1]-2) && 
		       (xp >= 1.0 && xp < colsF-dimSteps[0]-2)){
			vf = tri_cubic_convolve(imageF, (int)xp, (int)yp, (int)zp, xp, yp,
				          	zp, colsG, rowsG, layersG, sliceSizeG);
			/* clip at hard edges */
			if(vf < 0.0) vf = 0.0;
			imageG[sliceG+rowG+(int)x] = (int)vf;
		    }
	        }
	    }
	}

	status = 1;

	return status;

}


