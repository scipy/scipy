extern double MACHEP, PIO2, PI;
double ellie(), ellpe(), floor(), fabs();
double ellie2();

main()
{
double y, m, phi, e, E, phipi, y1;
int i, j, npi;

/* dprec();  */
m = 0.9;
E = ellpe(0.1);
for( j=-10; j<=10; j++ )
  {
    printf( "%d * PIO2\n", j );
    for( i=-2; i<=2; i++ )
      {
	phi = PIO2 * j + 50 * MACHEP * i;
	npi = floor(phi/PIO2);
	if( npi & 1 )
		npi += 1;
	phipi = phi - npi * PIO2;
	npi = floor(phi/PIO2);
	if( npi & 1 )
		npi += 1;
	phipi = phi - npi * PIO2;
	printf( "phi %.9e npi %d ", phi, npi );
	y1 = E * npi + ellie(phipi,m);
	y = ellie2( phi, m );
	printf( "y %.9e ", y );
	e = fabs(y - y1);
	if( y1 != 0.0 )
	  e /= y1;
	printf( "e %.4e\n", e );
      }
  }
}
