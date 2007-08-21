/*
 * numfmt.c -- $Id$
 * compute size, alignment, byte order, and floating point format
 *
 * Copyright (c) 1998.  See accompanying LEGAL file for details.
 */

/*
 *   The ANSI C standard requires:
 *    sizeof(char)  == 1
 *    sizeof(short) >= 2
 *    sizeof(int)   >= 2
 *    sizeof(long)  >= 4
 *    float exponent >= 8-bits    (1.0e-37 to 1.0e37 representable)
 *    float mantissa >= 18-bits   (1.0+1.0e-5 distinct from 1.0)
 *       --> sizeof(float) >= 4
 *    double exponent >= 8-bits   (1.0e-37 to 1.0e37 representable)
 *    double mantissa >= 31-bits  (1.0+1.0e-9 distinct from 1.0)
 *       --> sizeof(double) >= 5
 *
 *   Assume that:
 *   (1) char is an 8-bit byte
 *   (2) The order of the bytes can be described in terms of a
 *       word size (possibly more than one byte) in which the
 *       bytes are laid out either from most-to-least or least-to-most
 *       significant, and the words are also either most-to-least or
 *       least-to-most ordered.  (The VAX floating point formats are
 *       the only known examples of opposite ordering between words
 *       within the object and bytes within the words.)
 *   (3) The floating point format consists of an overall sign bit,
 *       a group of 31 or fewer contiguous exponent bits from which
 *       a fixed bias value can be subtracted to get the binary exponent,
 *       and contiguous mantissa bits.  The leading 1 bit on mantissas
 *       may be explicit or implicit.  Floating point zero must be
 *       represented as all bytes 0.
 *   The only known machines that this excludes are those which use
 *   hex exponents instead of binary exponents.
 */

/* data type size, alignment, and byte order, all in bytes */
static int type_size[7], type_alignment[8], type_order[6];

struct fp_layout {
  /* Addresses are in bits with 0 128-bit of the most significant byte
   * of the object, 8 the 128-bit of the 2nd most significant byte,
   * and so on until 8*size-1 is the 1-bit of the least significant
   * byte.  The actual storage order of the bytes is given by the
   * order member of the DataLayout.
   */
  int sgn_addr;              /* bit address of overall sign */
  int exp_addr, exp_size;    /* bit address and size of exponent */
  int man_addr, man_size;    /* bit address and size of mantissa */
  int man_norm; /* if non-zero, address of leading 1 in mantissa */
  long exp_bias;             /* exponent bias */
};
static struct fp_layout type_layout[2];

extern int p_num_formats(void);

static void swap_bytes(char *a, int nbytes, int ntot);
static int bit_addr(char *a, char *b, int nbytes);
static int integer_order(char *a, int nbytes);
static int fltpt_order(char *data, int nbytes, struct fp_layout *layout);

#ifndef STAND_ALONE
#include <stdio.h>

int
main(int argc, char *argv[])
{
  FILE *f = 0;
  if (argc==2) {
    if (argv[1][0]=='-' && !argv[1][1]) f = stdout;
    else if (argv[1][0]!='-') f = fopen(argv[1], "w");
  } else if (argc==1) {
    f = fopen("numfmt.h", "w");
  }
  if (!f) {
    puts("numfmt: failed to open output file");
    puts("numfmt usage:   numfmt [outfilename or '-' for stdout]");
    return 1;
  }

  p_num_formats();

  fprintf(f,"#define P_CHAR_ALIGN %d\n\n", type_alignment[0]);

  fprintf(f,"#define P_SHORT_SIZE %dL\n", type_size[1]);
  fprintf(f,"#define P_SHORT_ALIGN %d\n", type_alignment[1]);
  fprintf(f,"#define P_SHORT_ORDER %d\n\n", type_order[1]);

  fprintf(f,"#define P_INT_SIZE %dL\n", type_size[2]);
  fprintf(f,"#define P_INT_ALIGN %d\n", type_alignment[2]);
  fprintf(f,"#define P_INT_ORDER %d\n\n", type_order[2]);

  fprintf(f,"#define P_LONG_SIZE %dL\n", type_size[3]);
  fprintf(f,"#define P_LONG_ALIGN %d\n", type_alignment[3]);
  fprintf(f,"#define P_LONG_ORDER %d\n\n", type_order[3]);

  fprintf(f,"#define P_FLOAT_SIZE %dL\n", type_size[4]);
  fprintf(f,"#define P_FLOAT_ALIGN %d\n", type_alignment[4]);
  fprintf(f,"#define P_FLOAT_ORDER %d\n", type_order[4]);
  fprintf(f,"#define P_FLOAT_LAYOUT %d, %d, %d, %d, %d, %d, %ldL\n\n",
          type_layout[0].sgn_addr,
          type_layout[0].exp_addr, type_layout[0].exp_size,
          type_layout[0].man_addr, type_layout[0].man_size,
          type_layout[0].man_norm, type_layout[0].exp_bias);

  fprintf(f,"#define P_DOUBLE_SIZE %dL\n", type_size[5]);
  fprintf(f,"#define P_DOUBLE_ALIGN %d\n", type_alignment[5]);
  fprintf(f,"#define P_DOUBLE_ORDER %d\n", type_order[5]);
  fprintf(f,"#define P_DOUBLE_LAYOUT %d, %d, %d, %d, %d, %d, %ldL\n\n",
          type_layout[1].sgn_addr,
          type_layout[1].exp_addr, type_layout[1].exp_size,
          type_layout[1].man_addr, type_layout[1].man_size,
          type_layout[1].man_norm, type_layout[1].exp_bias);

  fprintf(f,"#define P_POINTER_SIZE %dL\n", type_size[6]);
  fprintf(f,"#define P_POINTER_ALIGN %d\n\n", type_alignment[6]);

  fprintf(f,"#define P_STRUCT_ALIGN %d\n", type_alignment[7]);

  return 0;
}
#endif

int
p_num_formats(void)
{
  if (!type_size[0]) {
    /* following data structures exhibit the definition of "alignment" */
    struct { char x; char y[1]; }               align_c;
    struct { char x; short y[1]; }              align_s;
    struct { char x; int y[1]; }                align_i;
    struct { char x; long y[1]; }               align_l;
    struct { char x; float y[1]; }              align_f;
    struct { char x; double y[1]; }             align_d;
    struct { char x; char *y[1]; }              align_p;
    struct { char x; struct { char x; } y[1]; } align_st;

    struct ppc_tester { double d; char c; };
    struct ppc_atester { double d[1]; char c; };
    struct ppc_stester { char c0; struct ppc_tester x; };
    unsigned int ppc_check = sizeof(struct ppc_tester);

    /* data block holds numbers used to decode floating point format */
    union {
      char c[600];
      float f[75];
      double d[75];  /* dimension must be at least 45+sizeof(double) */
    } data;

    type_size[0] = 1; /* ==sizeof(char) by definition */
    type_size[1] = sizeof(short);
    type_size[2] = sizeof(int);
    type_size[3] = sizeof(long);
    type_size[4] = sizeof(float);
    type_size[5] = sizeof(double);
    type_size[6] = sizeof(char *);

    /* fill in correct data alignments */
    type_alignment[0] = (int)((char *)(align_c.y) - (char *)(&align_c));
    type_alignment[1] = (int)((char *)(align_s.y) - (char *)(&align_s));
    type_alignment[2] = (int)((char *)(align_i.y) - (char *)(&align_i));
    type_alignment[3] = (int)((char *)(align_l.y) - (char *)(&align_l));
    type_alignment[4] = (int)((char *)(align_f.y) - (char *)(&align_f));
    type_alignment[5] = (int)((char *)(align_d.y) - (char *)(&align_d));
    type_alignment[6] = (int)((char *)(align_p.y) - (char *)(&align_p));
    type_alignment[7] = (int)((char *)(align_st.y) - (char *)(&align_st));

    if (ppc_check!=sizeof(double)+type_alignment[5] &&
        ppc_check==2*sizeof(double)) { /* cope with PowerPC idiocy */
      ppc_check = sizeof(struct ppc_stester);
      if (ppc_check == 2*sizeof(double)+type_alignment[5])
        type_alignment[7] = -1;  /* ibm misinterprets their own idiocy */
      else if ((ppc_check == 3*sizeof(double)) &&
               (sizeof(struct ppc_atester) ==
                sizeof(double)+type_alignment[5]))
        type_alignment[7] = -2;  /* gcc interpretation of ibm idiocy */
    }

    /* byte ordering for integers is directly testable -- just need to
     * compare set values in successive bytes by the shift left operator,
     * then re-interpret them as characters to find the byte order.  */
    {
      union {
        char c[8];
        short s[1];
        int i[1];
        long l[1];
      } u;
      int i;
      type_order[0] = 0;  /* char order not meaningful */
      for (i=0,u.l[0]=0 ; i<sizeof(short) ; i++) {u.s[0] |= (i+1)<<(8*i);}
      type_order[1] = integer_order(u.c, type_size[1]);
      for (i=0,u.l[0]=0 ; i<sizeof(int) ; i++)   {u.i[0] |= (i+1)<<(8*i);}
      type_order[2] = integer_order(u.c, type_size[2]);
      for (i=0,u.l[0]=0 ; i<sizeof(long) ; i++)  {u.l[0] |= (i+1UL)<<(8*i);}
      type_order[3] = integer_order(u.c, type_size[3]);
    }

    {
      float x;
      int n, i = 0;
      data.f[i++] = -1.0f;
      data.f[i++] = 2.0f;
      data.f[i++] = 65536.f*65536.f*65536.f*65536.f*65536.f*65536.f; /* 2^96 */
      data.f[i++] = 1.0f/(65536.f*65536.f*65536.f*65536.f*65536.f*65536.f);
      data.f[i++] = 1.0f;
      for (n=0,x=0.5f ; n<sizeof(double) ; n++,x*=(1.0f/256.0f))
        data.f[i+n] = data.f[i+n-1]+x;
      for (n=0,x=0.5f ; n<sizeof(double) ; n++,x*=(1.0f/256.0f))
        if (bit_addr((char*)&data.f[i+n],
                     (char*)&data.f[i+n+1],type_size[4])<0) break;
      i += sizeof(double);
      for (n=0 ; n<8 ; n++,x*=0.5f) data.f[i++] = 1.0f+x;
      for (n=0,x=1.5f ; n<16 ; n++,x*=2.0f) {
        data.f[i++] = x;
        data.f[i++] = -1.75f*x;
      }
    }
    type_order[4] = fltpt_order(data.c, type_size[4], &type_layout[0]);

    {
      double x;
      int n, i = 0;
      data.d[i++] = -1.0;
      data.d[i++] = 2.0;
      data.d[i++] = 65536.0*65536.0*65536.0*65536.0*65536.0*65536.0; /* 2^96 */
      data.d[i++] = 1.0/(65536.0*65536.0*65536.0*65536.0*65536.0*65536.0);
      data.d[i++] = 1.0;
      for (n=0,x=0.5 ; n<sizeof(double) ; n++,x*=(1.0/256.0))
        data.d[i+n] = data.d[i+n-1]+x;
      for (n=0,x=0.5 ; n<sizeof(double) ; n++,x*=(1.0/256.0))
        if (bit_addr((char*)&data.d[i+n],
                     (char*)&data.d[i+n+1],type_size[5])<0) break;
      i += sizeof(double);
      for (n=0 ; n<8 ; n++,x*=0.5) data.d[i++] = 1.0+x;
      for (n=0,x=1.5 ; n<16 ; n++,x*=2.0) {
        data.d[i++] = x;
        data.d[i++] = -1.75*x;
      }
    }
    type_order[5] = fltpt_order(data.c, type_size[5], &type_layout[1]);
  }

  return type_order[1]==0 || type_order[2]==0 ||
    type_order[3]==0 || type_order[4]==0 || type_order[5]==0;
}

static void
swap_bytes(char *a, int nbytes, int ntot)
{
  char tmp;
  int i;
  for (ntot-=nbytes; ntot>0 ; ntot-=nbytes, a+=nbytes)
    for (i=0 ; i<nbytes-1-i ; i++) {
      tmp = a[nbytes-1-i];
      a[nbytes-1-i] = a[i];
      a[i] = tmp;
    }
}

static int
bit_addr(char *a, char *b, int nbytes)
{
  int addr;
  for (addr=0 ; addr<nbytes ; addr++) if (a[addr]!=b[addr]) break;
  if (addr<nbytes) {
    unsigned int mask, diff = (a[addr] ^ b[addr]) & 0xff;
    addr<<= 3;
    for (mask=0x80 ; !(mask&diff) ; mask>>=1) addr++;
  } else {
    addr = -1;
  }
  return addr;
}

/* a[i] contains bytes of the number 0x0807060504030201
 * caution: Cray short uses 24-bit register but stores in 64-bit word
 *  - there may be other similar cases in which some bytes show as 0 */
static int
integer_order(char *a, int nbytes)
{
  int i, order=0;

  /* scan for 1 */
  for (i=0 ; i<nbytes ; i++) if (a[i]==1) break;

  if (i<nbytes) {
    int wdsz;
    if (i+i < nbytes) { /* least significant word first */
      wdsz = i+1;
      order = wdsz==1? -1 : wdsz;
    } else {        /* most significant word first */
      wdsz = nbytes-i;
      order = wdsz==1? 1 : -wdsz;
    }
    if (wdsz!=1) {
      /* no known integer format gets here
       * - branch is on analogy with VAX middle-endian floating point
       *   organized as words in one order consisting of bytes in other */
      if (nbytes%wdsz) order = 0;
      else swap_bytes(a, wdsz, nbytes);
    }
    if (order) {
      /* final check to be sure format is exactly as expected */
      if (i+i < nbytes) {
        for (i=0 ; i<nbytes ; i++) if (a[i]!=i+1) break;
        if (i<nbytes) for (; i<nbytes ; i++) if (a[i]) break;
      } else {
        for (i=0 ; i<nbytes ; i++) if (a[nbytes-1-i]!=i+1) break;
        if (i<nbytes) for (; i<nbytes ; i++) if (a[nbytes-1-i]) break;
      }
      if (i<nbytes) order = 0;
    }
  }

  return order;
}

static int
fltpt_order(char *data, int nbytes, struct fp_layout *layout)
{
  char *a, *b;
  int i, addr0, addr1, step, step1, step2, order;
  int nswap = (45 + (int)sizeof(double))*nbytes;
  long bias;

  /* mantissa order determined by watching bits during successive adds
   * of (1.0/256.0), which is one bit moving one byte per pass */
  for (order=1 ; ; ) {
    a = data + 4*nbytes;  /* a[i] = 1.0 + sum(0.5/256.0^i) */
    addr1 = bit_addr(&a[0],&a[nbytes],nbytes);
    a += nbytes;
    addr0 = bit_addr(&a[0],&a[nbytes],nbytes);
    a += nbytes;
    step1 = addr0 - addr1;
    step2 = step = 0;
    for (i=2 ; i<sizeof(double) ; i++, a+=nbytes) {
      addr1 = bit_addr(&a[0],&a[nbytes],nbytes);
      if (addr1 < 0) break;
      step = addr1 - addr0;
      if (step!=step1) {
        if (!step2) {
          step2 = step1;
          step1 = step;
        } else if (step!=step2) {
          return 0;
        }
      }
      addr0 = addr1;
    }

    if (step2) { /* middle-endian (VAX) */
      if (order!=1 || (step1<0)==(step2<0) ||
          ((step1<0? -step1 : step1)!=8 &&
           (step2<0? -step2 : step2)!=8))
        return 0;

      step = step2 - step1;
      if (step < 0) step = -step;
      order = step/16; /* word size in bytes (step=2*bit size) */
      if ((step&15) || (nbytes%order))
        return 0;

      /* swap bytes within putative words and try again */
      swap_bytes(data, order, nswap);

    } else {
      if (i>=sizeof(double)) return 0;
      if (step==-8)
        order = -order;
      else if (step!=8)
        return 0;
      break;
    }
  }

  /* convert entire data series to big-endian (msb first) order */
  if (order<0) swap_bytes(data, nbytes, nswap);

  /* mantissa size, look at scan to find where 1.0+x!=1.0
   * a[0] = 1.0 + sum(0.5/256.0^i), but next number is same */
  step = 8*i-8;  /* x = 0.5^step still has 1.0+x!=1.0 */
  a = data + 4*nbytes;  /* 1.0 */
  b = data + (5+(int)sizeof(double))*nbytes;
  for (i=0 ; i<8 ; i++, b+=nbytes)
    if (bit_addr(a,b,nbytes)<0) break;
  layout->man_size = step+i;

  /* mantissa address, compare 1.0 and 1.5 */
  a = data + 4*nbytes;  /* 1.0 */
  b = data + 5*nbytes;  /* 1.5 */
  layout->man_addr = bit_addr(a, b, nbytes);

  /* mantissa normalization, test if bit before man_addr always set */
  layout->man_norm = 0;
  if (layout->man_addr > 0) {
    int ba = (layout->man_addr-1) / 8;
    int mask = 0x80U >> ((layout->man_addr-1) % 8);
    a = data + (5+(int)sizeof(double)+8)*nbytes;
    for (i=0 ; i<16 ; i++,a+=2*nbytes)
      if ((a[ba]&mask)==0 || (a[nbytes+ba]&mask)==0) break;
    if (i==16) {
      layout->man_norm = 1;
      layout->man_addr--;
      layout->man_size++;
    }
  }

  /* sign bit address, compare -1.0 and 1.0 */
  a = data;             /* -1.0 */
  b = data + 4*nbytes;  /*  1.0 */
  layout->sgn_addr = bit_addr(a, b, nbytes);

  /* exponent address, compare 2.0^96, 0.5^96 */
  a = data + 2*nbytes;  /* 2.0^96 */
  b = data + 3*nbytes;  /* 0.5^96 */
  layout->exp_addr = bit_addr(a, b, nbytes);

  /* exponent size, compare 1.0 and 2.0 in reverse order */
  a = data + 4*nbytes;  /* 1.0 */
  b = data + nbytes;    /* 2.0 */
  for (i=nbytes-1 ; i>=0 ; i--)
    if (a[i]!=b[i]) {
      step1 = a[i]^b[i];
      for (addr1=0,step=1 ; addr1<8 ; addr1++,step<<=1)
        if (step1 & step) break;
      addr1 = 8*i + 7-addr1;
      break;
    }
  layout->exp_size = addr1 - layout->exp_addr + 1;

  /* exponent bias, look at exponent for 1.0 */
  a = data + 4*nbytes;  /* 1.0 */
  addr0 = layout->exp_addr / 8;
  bias = a[addr0] & ((1 << (8-(layout->exp_addr%8))) - 1);
  step = addr1 % 8;
  addr1 = addr1 / 8;
  for (addr0++ ; addr0<addr1 ; addr0++)
    bias = (bias<<8) | a[addr0];
  bias = (bias<<(step+1)) | (((unsigned char)a[addr0]) >> (7-step));
  layout->exp_bias = bias;
  if (layout->man_norm) layout->exp_bias--;

  /* check overlaps against things like hex exponent machines */

  if ((layout->sgn_addr>=layout->exp_addr &&
       layout->sgn_addr<layout->exp_addr+layout->exp_size) ||
      (layout->sgn_addr>=layout->man_addr &&
       layout->sgn_addr<layout->man_addr+layout->man_size)) return 0;

  if (layout->exp_addr==layout->man_addr) return 0;
  if (layout->exp_addr<layout->man_addr) {
    if (layout->exp_addr+layout->exp_size>layout->man_addr) return 0;
  } else {
    if (layout->man_addr+layout->man_size>layout->exp_addr) return 0;
  }

  return order;
}
