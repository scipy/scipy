
/*--------------------------------------------------------------------*/

/*
 *  This Quickselect routine is based on the algorithm described in
 *  "Numerical recipies in C", Second Edition,
 *  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 */

void f_medfilt2(float*,float*,int*,int*);
void d_medfilt2(double*,double*,int*,int*);
void b_medfilt2(unsigned char*,unsigned char*,int*,int*);
extern char *check_malloc (int);

#define ELEM_SWAP(a,b) { register float t=(a);(a)=(b);(b)=t; }

float f_quick_select(float arr[], int n) 
{
    int low, high ;
    int median;
    int middle, ll, hh;

    low = 0 ; high = n-1 ; median = (low + high) / 2;
    for (;;) {
        if (high <= low) /* One element only */
            return arr[median] ;

        if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                ELEM_SWAP(arr[low], arr[high]) ;
            return arr[median] ;
        }

    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high])    ELEM_SWAP(arr[middle], arr[high]) ;
    if (arr[low] > arr[high])       ELEM_SWAP(arr[low], arr[high]) ;
    if (arr[middle] > arr[low])     ELEM_SWAP(arr[middle], arr[low]) ;

    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAP(arr[middle], arr[low+1]) ;

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
        do ll++; while (arr[low] > arr[ll]) ;
        do hh--; while (arr[hh]  > arr[low]) ;

        if (hh < ll)
        break;

        ELEM_SWAP(arr[ll], arr[hh]) ;
    }

    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP(arr[low], arr[hh]) ;

    /* Re-set active partition */
    if (hh <= median)
        low = ll;
        if (hh >= median)
        high = hh - 1;
    }
}

#undef ELEM_SWAP


#define ELEM_SWAP(a,b) { register double t=(a);(a)=(b);(b)=t; }

double d_quick_select(double arr[], int n) 
{
    int low, high ;
    int median;
    int middle, ll, hh;

    low = 0 ; high = n-1 ; median = (low + high) / 2;
    for (;;) {
        if (high <= low) /* One element only */
            return arr[median] ;

        if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                ELEM_SWAP(arr[low], arr[high]) ;
            return arr[median] ;
        }

    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high])    ELEM_SWAP(arr[middle], arr[high]) ;
    if (arr[low] > arr[high])       ELEM_SWAP(arr[low], arr[high]) ;
    if (arr[middle] > arr[low])     ELEM_SWAP(arr[middle], arr[low]) ;

    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAP(arr[middle], arr[low+1]) ;

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
        do ll++; while (arr[low] > arr[ll]) ;
        do hh--; while (arr[hh]  > arr[low]) ;

        if (hh < ll)
        break;

        ELEM_SWAP(arr[ll], arr[hh]) ;
    }

    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP(arr[low], arr[hh]) ;

    /* Re-set active partition */
    if (hh <= median)
        low = ll;
        if (hh >= median)
        high = hh - 1;
    }
}

#undef ELEM_SWAP

#define ELEM_SWAP(a,b) { register unsigned char t=(a);(a)=(b);(b)=t; }

unsigned char b_quick_select(unsigned char arr[], int n) 
{
    int low, high ;
    int median;
    int middle, ll, hh;

    low = 0 ; high = n-1 ; median = (low + high) / 2;
    for (;;) {
        if (high <= low) /* One element only */
            return arr[median] ;

        if (high == low + 1) {  /* Two elements only */
            if (arr[low] > arr[high])
                ELEM_SWAP(arr[low], arr[high]) ;
            return arr[median] ;
        }

    /* Find median of low, middle and high items; swap into position low */
    middle = (low + high) / 2;
    if (arr[middle] > arr[high])    ELEM_SWAP(arr[middle], arr[high]) ;
    if (arr[low] > arr[high])       ELEM_SWAP(arr[low], arr[high]) ;
    if (arr[middle] > arr[low])     ELEM_SWAP(arr[middle], arr[low]) ;

    /* Swap low item (now in position middle) into position (low+1) */
    ELEM_SWAP(arr[middle], arr[low+1]) ;

    /* Nibble from each end towards middle, swapping items when stuck */
    ll = low + 1;
    hh = high;
    for (;;) {
        do ll++; while (arr[low] > arr[ll]) ;
        do hh--; while (arr[hh]  > arr[low]) ;

        if (hh < ll)
        break;

        ELEM_SWAP(arr[ll], arr[hh]) ;
    }

    /* Swap middle item (in position low) back into correct position */
    ELEM_SWAP(arr[low], arr[hh]) ;

    /* Re-set active partition */
    if (hh <= median)
        low = ll;
        if (hh >= median)
        high = hh - 1;
    }
}

#undef ELEM_SWAP

/* 2-D median filter with zero-padding on edges. */
void d_medfilt2(in, out, Nwin, Ns)
double *in, *out;
int *Nwin, *Ns;
{ 
  int nx, ny, hN[2];
  int pre_x, pre_y, pos_x, pos_y;
  int subx, suby, k, totN;
  double *myvals, *fptr1, *fptr2, *ptr1, *ptr2;

  totN = Nwin[0] * Nwin[1];
  myvals = (double *) check_malloc( totN * sizeof(double));

  hN[0] = Nwin[0] >> 1;
  hN[1] = Nwin[1] >> 1;
  ptr1 = in;
  fptr1 = out;
  for (ny = 0; ny < Ns[0]; ny++)
    for (nx = 0; nx < Ns[1]; nx++) {
      pre_x = hN[1];
      pre_y = hN[0];
      pos_x = hN[1];
      pos_y = hN[0];
      if (nx < hN[1]) pre_x = nx;
      if (nx >= Ns[1] - hN[1]) pos_x = Ns[1] - nx - 1;
      if (ny < hN[0]) pre_y = ny;
      if (ny >= Ns[0] - hN[0]) pos_y = Ns[0] - ny - 1;
      fptr2 = myvals;
      ptr2 = ptr1 - pre_x - pre_y*Ns[1];
      for (suby = -pre_y; suby <= pos_y; suby++) {
	for (subx = -pre_x; subx <= pos_x; subx++) 	  
	  *fptr2++ = *ptr2++;
	ptr2 += Ns[1] - (pre_x + pos_x + 1);
      }
      ptr1++;

      /* Zero pad */
      for (k = (pre_x + pos_x + 1)*(pre_y + pos_y + 1); k < totN; k++)
	*fptr2++ = 0.0;

      /*      *fptr1++ = median(myvals,totN); */
      *fptr1++ = d_quick_select(myvals,totN);
    }
}


/* 2-D median filter with zero-padding on edges. */
void f_medfilt2(in, out, Nwin, Ns)
float *in, *out;
int *Nwin, *Ns;
{ 
  int nx, ny, hN[2];
  int pre_x, pre_y, pos_x, pos_y;
  int subx, suby, k, totN;
  float *myvals, *fptr1, *fptr2, *ptr1, *ptr2;

  totN = Nwin[0] * Nwin[1];
  myvals = (float *) check_malloc( totN * sizeof(float));

  hN[0] = Nwin[0] >> 1;
  hN[1] = Nwin[1] >> 1;
  ptr1 = in;
  fptr1 = out;
  for (ny = 0; ny < Ns[0]; ny++)
    for (nx = 0; nx < Ns[1]; nx++) {
      pre_x = hN[1];
      pre_y = hN[0];
      pos_x = hN[1];
      pos_y = hN[0];
      if (nx < hN[1]) pre_x = nx;
      if (nx >= Ns[1] - hN[1]) pos_x = Ns[1] - nx - 1;
      if (ny < hN[0]) pre_y = ny;
      if (ny >= Ns[0] - hN[0]) pos_y = Ns[0] - ny - 1;
      fptr2 = myvals;
      ptr2 = ptr1 - pre_x - pre_y*Ns[1];
      for (suby = -pre_y; suby <= pos_y; suby++) {
	for (subx = -pre_x; subx <= pos_x; subx++) 	  
	  *fptr2++ = *ptr2++;
	ptr2 += Ns[1] - (pre_x + pos_x + 1);
      }
      ptr1++;

      /* Zero pad */
      for (k = (pre_x + pos_x + 1)*(pre_y + pos_y + 1); k < totN; k++)
	*fptr2++ = 0.0;

      /*      *fptr1++ = median(myvals,totN); */
      *fptr1++ = f_quick_select(myvals,totN);
    }
}


/* 2-D median filter with zero-padding on edges. */
void b_medfilt2(in, out, Nwin, Ns)
unsigned char *in, *out;
int *Nwin, *Ns;
{ 
  int nx, ny, hN[2];
  int pre_x, pre_y, pos_x, pos_y;
  int subx, suby, k, totN;
  unsigned char *myvals, *fptr1, *fptr2, *ptr1, *ptr2;

  totN = Nwin[0] * Nwin[1];
  myvals = (unsigned char *) check_malloc( totN * sizeof(unsigned char));

  hN[0] = Nwin[0] >> 1;
  hN[1] = Nwin[1] >> 1;
  ptr1 = in;
  fptr1 = out;
  for (ny = 0; ny < Ns[0]; ny++)
    for (nx = 0; nx < Ns[1]; nx++) {
      pre_x = hN[1];
      pre_y = hN[0];
      pos_x = hN[1];
      pos_y = hN[0];
      if (nx < hN[1]) pre_x = nx;
      if (nx >= Ns[1] - hN[1]) pos_x = Ns[1] - nx - 1;
      if (ny < hN[0]) pre_y = ny;
      if (ny >= Ns[0] - hN[0]) pos_y = Ns[0] - ny - 1;
      fptr2 = myvals;
      ptr2 = ptr1 - pre_x - pre_y*Ns[1];
      for (suby = -pre_y; suby <= pos_y; suby++) {
	for (subx = -pre_x; subx <= pos_x; subx++) 	  
	  *fptr2++ = *ptr2++;
	ptr2 += Ns[1] - (pre_x + pos_x + 1);
      }
      ptr1++;

      /* Zero pad */
      for (k = (pre_x + pos_x + 1)*(pre_y + pos_y + 1); k < totN; k++)
	*fptr2++ = 0.0;

      /*      *fptr1++ = median(myvals,totN); */
      *fptr1++ = b_quick_select(myvals,totN);
    }
}
