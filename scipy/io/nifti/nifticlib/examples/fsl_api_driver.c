/*******************************************************************
 *
 * fsl_api_driver.c
 *
 * Test fslio API 
 *
 * Usage: fsl_api_driver <command> <dataset(s)> [params]
 *
 * 
 *      print <dataset>
 *              print dataset header fields
 *      peek <dataset> X Y Z T
 *              print value at voxel location X Y Z T (0-based index)
 *      timecourse <dataset> X Y Z
 *              print timecourse at voxel location X Y Z (0-based index)
 *
 * Kate Fissell, University of Pittsburgh
 * 12/04
 *******************************************************************/


#include <stdio.h>
#include <strings.h>

#include <nifti1.h>
#include <fslio.h>

void pusage(char *cmd);

int main(int argc, char * argv[])
{

FSLIO *fslio;
void *buffer;
char *f1name;
double ***vol;
int x,y,z,t;


/*** process commandline parameters */
if (argc < 2) {
        pusage(argv[0]);
        exit(1);
}


/************************* PRINT ***************************/
if (!strncmp(argv[1],"print",5)) {
        if (argc != 3) {
                fprintf(stderr, "\nError, print command takes one parameter: print <dataset>\n");
                exit(1);
        }

        f1name = argv[2];
        /** open nifti dataset */
        fslio = FslInit();
        buffer = FslReadAllVolumes(fslio,f1name);
        if (buffer == NULL) {
                fprintf(stderr, "\nError opening and reading %s.\n",f1name);
                exit(1);
        }

        nifti_image_infodump(fslio->niftiptr);

        exit(0);
}



/************************* PEEK ***************************/
if (!strncmp(argv[1],"peek",4)) {
        if (argc != 7) {
                fprintf(stderr, "\nError, peek command takes five parameters: peek <dataset> X Y Z T\n");
                exit(1);
        }

        /**** get inputs */
        f1name = argv[2];
        x = atoi(argv[3]);
        y = atoi(argv[4]);
        z = atoi(argv[5]);
        t = atoi(argv[6]);

        /** open nifti dataset header */
        fslio = FslReadHeader(f1name);
        if (fslio == NULL) {
                fprintf(stderr, "\nError, could not read header info for %s.\n",f1name);
                exit(1);
        }

        /**** check inputs */
        if ( (x<0) || (x>=fslio->niftiptr->nx) ) {
                fprintf(stderr, "\nError: x index (%d) out of range [0..%d]\n",x,fslio->niftiptr->nx-1);
                exit(1);
        }
        if ( (y<0) || (y>=fslio->niftiptr->ny) ) {
                fprintf(stderr, "\nError: y index (%d) out of range [0..%d]\n",y,fslio->niftiptr->ny-1);
                exit(1);
        }
        if ( (z<0) || (z>=fslio->niftiptr->nz) ) {
                fprintf(stderr, "\nError: z index (%d) out of range [0..%d]\n",z,fslio->niftiptr->nz-1);
                exit(1);
        }
        if ( (t<0) || (t>=fslio->niftiptr->nt) ) {
                fprintf(stderr, "\nError: t index (%d) out of range [0..%d]\n",t,fslio->niftiptr->nt-1);
                exit(1);
        }

        /*** get volume data as scaled doubles */

        vol = FslGetVolumeAsScaledDouble(fslio,t);
        if (vol == NULL) {
                fprintf(stderr, "\nError accessing %s\n",f1name);
                exit(1);
        }
        else {
                fprintf(stderr, "\nLocation %d %d %d %d: %.4f\n",x,y,z,t,vol[z][y][x]);
                exit(0);
        }

}



fprintf(stderr, "\nError, unrecognized command %s\n",argv[1]);
pusage(argv[0]);


exit(1);
}


void pusage(char *cmd) {
        fprintf(stderr, "\n%s is a small driver program to test out the fslio API.\n",cmd);
        fprintf(stderr, "\nUsage: %s <command> <dataset(s)> <parameters>",cmd);
        fprintf(stderr, "\n\n\tCommands:");
        fprintf(stderr, "\n\tprint <dataset>\t\t\tprint dataset header");
        fprintf(stderr, "\n\tpeek <dataset> X Y Z T\t\tprint dataset value at location (0-based) (x,y,z,t).");

        fprintf(stderr, "\n");
return;
}
