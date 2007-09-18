/*********************************************************************
 *
 * Very simple code snippets to read/write nifti1 files
 * This code is placed in the public domain.
 *
 * If you are the type who doesn't want to use a file format unless
 * you can write your own i/o code in less than 30minutes, this
 * example is for you.
 *
 * This code does not deal with wrong-endian data, compressed data,
 * the new qform/sform orientation codes, parsing filenames, volume-
 * wise or timecourse-wise data access or any of a million other very useful
 * things that are in the niftilib i/o reference libraries.
 * We encourage people to use the niftilib reference library and send
 * feedback/suggestions, see http://niftilib.sourceforge.net/
 * But, if that is too much to tackle and you just want to jump in, this
 * code is a starting point.
 * This code was written for maximum readability, not for the greatest
 * coding style.
 *
 *
 * If you are already a little familiar with reading/writing Analyze
 * files of some flavor, and maybe even have some of your own code, here
 * are the most important things to be aware of in transitioning to nifti1:
 *
 * 1. nii vs .hdr/.img
 *      nifti1 datasets can be stored either in .hdr/.img pairs of files
 *      or in 1 .nii file.  In a .nii file the data will start at the byte
 *      specified by the vox_offset field, which will be 352 if no extensions
 *      have been added.  And, nifti1 really does like that magic field set
 *      to "n+1" for .nii and "ni1" for .img/.hdr
 *
 * 2. scaling
 *      nifti1 datasets can contain a scaling factor.  You need to check the
 *      scl_slope field and if that isn't 0, scale your data by 
 *      Y * scl_slope  + scl_inter
 *
 * 3. extensions
 *      nifti1 datasets can have some "extension data" stuffed after the 
 *      regular header.  You can just ignore it, but, be aware that a
 *      .hdr file may be longer than 348 bytes, and, in a .nii file
 *      you can't just jump to byte 352, you need to use the vox_offset
 *      field to get the start of the image data.
 *
 * 4. new datatypes
 *      nifti1 added a few new datatypes that were not in the Analyze 7.5
 *      format from which nifti1 is derived.  If you're just working with
 *      your own data this is not an issue but if you get a foreign nifti1
 *      file, be aware of exotic datatypes like DT_COMPLEX256 and mundane
 *      things like DT_UINT16.
 *
 * 5. other stuff
 *     nifti1 really does like the dim[0] field set to the number of
 *     dimensions of the dataset.  Other Analyze flavors might not
 *     have been so scrupulous about that.
 *     nifti1 has a bunch of other new fields such as intent codes,
 *     qform/sform, etc. but, if you just want to get your hands on
 *     the data blob you can ignore these.  Example use of these fields
 *     is in the niftilib reference libraries.
 *
 *
 *
 * To compile:
 * You need to put a copy of the nifti1.h header file in this directory.
 * It can be obtained from the NIFTI homepage  http://nifti.nimh.nih.gov/
 * or from the niftilib SourceForge site http://niftilib.sourceforge.net/
 * 
 * cc -o nifti1_read_write nifti1_read_write.c
 * 
 * 
 * To run:
 * nifti1_read_write -w abc.nii abc.nii
 * nifti1_read_write -r abc.nii abc.nii
 * 
 * 
 * The read method is hardcoded to read float32 data.  To change
 * to your datatype, just change the line:
 * typedef float MY_DATATYPE;
 *
 * The write method is hardcoded to write float32 data.  To change
 * to your datatype, change the line:
 * typedef float MY_DATATYPE;
 * and change the lines:
 * hdr.datatype = NIFTI_TYPE_FLOAT32;
 * hdr.bitpix = 32;
 *
 *
 * Written by Kate Fissell, University of Pittsburgh, May 2005.
 *
 *********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "nifti1.h"

typedef float MY_DATATYPE;

#define MIN_HEADER_SIZE 348
#define NII_HEADER_SIZE 352


main(argc,argv) 
int argc;
char *argv[];
{

char *hdr_file, *data_file;
short do_read=0;
short do_write=0;


/********** process commandline parameters */
if (argc != 4) {
        fprintf(stderr, "\nUsage: %s <-r|-w> <header file> <data file>\n",argv[0]);
        exit(1);
}

if (!strncmp(argv[1],"-r",2))
        do_read=1;
else if (!strncmp(argv[1],"-w",2))
        do_write=1;
else {
        fprintf(stderr, "\nUsage: %s <-r|-w> <header file> <data file>\n",argv[0]);
        exit(1);
}

hdr_file = argv[2];
data_file = argv[3];



/********** do the simple read or write */
if (do_read)
        read_nifti_file(hdr_file, data_file);
else if (do_write)
        write_nifti_file(hdr_file, data_file);


exit(0);
}

/**********************************************************************
 *
 * read_nifti_file
 *
 **********************************************************************/
int read_nifti_file(hdr_file, data_file)
char *hdr_file, *data_file;
{
nifti_1_header hdr;
FILE *fp;
int ret,i;
double total;
MY_DATATYPE *data=NULL;


/********** open and read header */
fp = fopen(hdr_file,"r");
if (fp == NULL) {
        fprintf(stderr, "\nError opening header file %s\n",hdr_file);
        exit(1);
}
ret = fread(&hdr, MIN_HEADER_SIZE, 1, fp);
if (ret != 1) {
        fprintf(stderr, "\nError reading header file %s\n",hdr_file);
        exit(1);
}
fclose(fp);


/********** print a little header information */
fprintf(stderr, "\n%s header information:",hdr_file);
fprintf(stderr, "\nXYZT dimensions: %d %d %d %d",hdr.dim[1],hdr.dim[2],hdr.dim[3],hdr.dim[4]);
fprintf(stderr, "\nDatatype code and bits/pixel: %d %d",hdr.datatype,hdr.bitpix);
fprintf(stderr, "\nScaling slope and intercept: %.6f %.6f",hdr.scl_slope,hdr.scl_inter);
fprintf(stderr, "\nByte offset to data in datafile: %ld",(long)(hdr.vox_offset));
fprintf(stderr, "\n");


/********** open the datafile, jump to data offset */
fp = fopen(data_file,"r");
if (fp == NULL) {
        fprintf(stderr, "\nError opening data file %s\n",data_file);
        exit(1);
}

ret = fseek(fp, (long)(hdr.vox_offset), SEEK_SET);
if (ret != 0) {
        fprintf(stderr, "\nError doing fseek() to %ld in data file %s\n",(long)(hdr.vox_offset), data_file);
        exit(1);
}


/********** allocate buffer and read first 3D volume from data file */
data = (MY_DATATYPE *) malloc(sizeof(MY_DATATYPE) * hdr.dim[1]*hdr.dim[2]*hdr.dim[3]);
if (data == NULL) {
        fprintf(stderr, "\nError allocating data buffer for %s\n",data_file);
        exit(1);
}
ret = fread(data, sizeof(MY_DATATYPE), hdr.dim[1]*hdr.dim[2]*hdr.dim[3], fp);
if (ret != hdr.dim[1]*hdr.dim[2]*hdr.dim[3]) {
        fprintf(stderr, "\nError reading volume 1 from %s (%d)\n",data_file,ret);
        exit(1);
}
fclose(fp);


/********** scale the data buffer  */
if (hdr.scl_slope != 0) {
        for (i=0; i<hdr.dim[1]*hdr.dim[2]*hdr.dim[3]; i++)
                data[i] = (data[i] * hdr.scl_slope) + hdr.scl_inter;
}


/********** print mean of data */
total = 0;
for (i=0; i<hdr.dim[1]*hdr.dim[2]*hdr.dim[3]; i++)
        total += data[i];
total /= (hdr.dim[1]*hdr.dim[2]*hdr.dim[3]);
fprintf(stderr, "\nMean of volume 1 in %s is %.3f\n",data_file,total);


return(0);
}


/**********************************************************************
 *
 * write_nifti_file
 * 
 * write a sample nifti1 (.nii) data file
 * datatype is float32
 * XYZT size is 64x64x16x10
 * XYZ voxel size is 1mm
 * TR is 1500ms
 *
 **********************************************************************/
int write_nifti_file(hdr_file, data_file)
char *hdr_file, *data_file;
{
nifti_1_header hdr;
nifti1_extender pad={0,0,0,0};
FILE *fp;
int ret,i;
MY_DATATYPE *data=NULL;
short do_nii;


/********** make sure user specified .hdr/.img or .nii/.nii */
if ( (strlen(hdr_file) < 4) || (strlen(data_file) < 4) ) {
        fprintf(stderr, "\nError: write files must end with .hdr/.img or .nii/.nii extension\n");
        exit(1);
}

if ( (!strncmp(hdr_file+(strlen(hdr_file)-4), ".hdr",4)) &&
     (!strncmp(data_file+(strlen(data_file)-4), ".img",4)) ) {
        do_nii = 0;
}
else if ( (!strncmp(hdr_file+(strlen(hdr_file)-4), ".nii",4)) &&
     (!strncmp(data_file+(strlen(data_file)-4), ".nii",4)) ) {
        do_nii = 1;
}
else {
        fprintf(stderr, "\nError: file(s) to be written must end with .hdr/.img or .nii/.nii extension\n");
        exit(1);
}
        

/********** fill in the minimal default header fields */
bzero((void *)&hdr, sizeof(hdr));
hdr.sizeof_hdr = MIN_HEADER_SIZE;
hdr.dim[0] = 4;
hdr.dim[1] = 64;
hdr.dim[2] = 64;
hdr.dim[3] = 16;
hdr.dim[4] = 10;
hdr.datatype = NIFTI_TYPE_FLOAT32;
hdr.bitpix = 32;
hdr.pixdim[1] = 1.0;
hdr.pixdim[2] = 1.0;
hdr.pixdim[3] = 1.0;
hdr.pixdim[4] = 1.5;
if (do_nii)
        hdr.vox_offset = (float) NII_HEADER_SIZE;
else
        hdr.vox_offset = (float)0;
hdr.scl_slope = 100.0;
hdr.xyzt_units = NIFTI_UNITS_MM | NIFTI_UNITS_SEC;
if (do_nii)
        strncpy(hdr.magic, "n+1\0", 4);
else
        strncpy(hdr.magic, "ni1\0", 4);


/********** allocate buffer and fill with dummy data  */
data = (MY_DATATYPE *) malloc(sizeof(MY_DATATYPE) * hdr.dim[1]*hdr.dim[2]*hdr.dim[3]*hdr.dim[4]);
if (data == NULL) {
        fprintf(stderr, "\nError allocating data buffer for %s\n",data_file);
        exit(1);
}

for (i=0; i<hdr.dim[1]*hdr.dim[2]*hdr.dim[3]*hdr.dim[4]; i++)
        data[i] = i / hdr.scl_slope;


/********** write first 348 bytes of header   */
fp = fopen(hdr_file,"w");
if (fp == NULL) {
        fprintf(stderr, "\nError opening header file %s for write\n",hdr_file);
        exit(1);
}
ret = fwrite(&hdr, MIN_HEADER_SIZE, 1, fp);
if (ret != 1) {
        fprintf(stderr, "\nError writing header file %s\n",hdr_file);
        exit(1);
}


/********** if nii, write extender pad and image data   */
if (do_nii == 1) {

        ret = fwrite(&pad, 4, 1, fp);
        if (ret != 1) {
                fprintf(stderr, "\nError writing header file extension pad %s\n",hdr_file);
                exit(1);
        }

        ret = fwrite(data, (size_t)(hdr.bitpix/8), hdr.dim[1]*hdr.dim[2]*hdr.dim[3]*hdr.dim[4], fp);
        if (ret != hdr.dim[1]*hdr.dim[2]*hdr.dim[3]*hdr.dim[4]) {
                fprintf(stderr, "\nError writing data to %s\n",hdr_file);
                exit(1);
        }

        fclose(fp);
}


/********** if hdr/img, close .hdr and write image data to .img */
else {

        fclose(fp);     /* close .hdr file */

        fp = fopen(data_file,"w");
        if (fp == NULL) {
                fprintf(stderr, "\nError opening data file %s for write\n",data_file);
                exit(1);
        }
        ret = fwrite(data, (size_t)(hdr.bitpix/8), hdr.dim[1]*hdr.dim[2]*hdr.dim[3]*hdr.dim[4], fp);
        if (ret != hdr.dim[1]*hdr.dim[2]*hdr.dim[3]*hdr.dim[4]) {
                fprintf(stderr, "\nError writing data to %s\n",data_file);
                exit(1);
        }

        fclose(fp);
}



return(0);
}
