/* ---------------------------------------------------------------------- */

void filter

     (float    **ur,        /* real part of Fourier coeffs, changed */
      float    **ui,        /* imag. part of Fourier coeffs, changed */
      float    **hr,        /* real part of Fourier kernel, unchanged */
      float    **hi,        /* imag. part of Fourier kernel, unchanged */
      float    param,       /* filter parameter */
      long     nx,          /* pixel number in x direction */
      long     ny)          /* pixel number in y direction */

/*
 allows to filter the Fourier coefficients
*/

{
  long   i, j;              /* loop variables */
  float  N;                 /* denominator */
  float  vr,vi;             /* auxiliary variables for cplx arithm. */

/* ---- compute filtered coefficients ---- */

/** ---- INSERT CODE HERE ---- **/
for (i=0; i<nx; i++)
{
  for (j=0; j<ny; j++)
  {
    vr = ur[i][j];
    vi = ui[i][j];

    N = hr[i][j] * hr[i][j] + hi[i][j] * hi[i][j] + param;
    ur[i][j] = ( hr[i][j]*vr + hi[i][j]*vi ) / N;
    ui[i][j] = ( hr[i][j]*vi - hi[i][j]*vr ) / N;
  }
}

return;
}
