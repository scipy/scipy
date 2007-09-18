/*
 * nifti_test2 -- this exercises a particular bug, whereby nifti_copy_nim_info was
 * not working correctly, which left the copy of an image pointing to the same
 * extension list, which could cause all sorts of grief.
 */
#include <nifti1_io.h>
int main (int argc, char *argv[])
{
  /*
   * create a 'dummy' image
   */
  nifti_image *i1 = nifti_simple_init_nim();
  nifti_image *i2;
  
  /*
   * add an extension to the dummy
   */
  static char ext[] = "THIS IS A TEST";
  nifti_add_extension(i1,ext,sizeof(ext),NIFTI_ECODE_COMMENT);
  /*
   * make a new nim from the dummy
   */
  i2 = nifti_copy_nim_info(i1);
  /*
   * if the bug isn't fixed in niftilib, the second nifti_image_free
   * will fail because both nims point to the same extensions. With gcc
   * this will abort inside the standard library
   */
  nifti_image_free(i1);
  nifti_image_free(i2);
  return 0;
}
