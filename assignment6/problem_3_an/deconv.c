#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "deconv-filter.c"


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                  DECONVOLUTION VIA THE FOURIER DOMAIN                    */
/*                                                                          */
/*         (Copyright Joachim Weickert and Martin Welk, 12/2004)            */
/*                                                                          */
/*--------------------------------------------------------------------------*/


/* 
 features:
 - DFT or FFT, whatever is possible 
 - reflecting extension at the boundaries 
*/


/*--------------------------------------------------------------------------*/

void alloc_vector

     (float **vector,   /* vector */
      long  n)          /* size */

     /* allocates storage for a vector of size n */


{
*vector = (float *) malloc (n * sizeof(float));
if (*vector == NULL)
   {
   printf("alloc_vector: not enough storage available\n");
   exit(1);
   }
return;
}

/*--------------------------------------------------------------------------*/

void alloc_matrix

     (float ***matrix,  /* matrix */
      long  nx,         /* size in x direction */
      long  ny)         /* size in y direction */

     /* allocates storage for matrix of size nx * ny */


{
long i;

*matrix = (float **) malloc (nx * sizeof(float *));
if (*matrix == NULL)
   {
   printf("alloc_matrix: not enough storage available\n");
   exit(1);
   }
for (i=0; i<nx; i++)
    {
    (*matrix)[i] = (float *) malloc (ny * sizeof(float));
    if ((*matrix)[i] == NULL)
       {
       printf("alloc_matrix: not enough storage available\n");
       exit(1);
       }
    }
return;
}

/*--------------------------------------------------------------------------*/

void disalloc_vector

     (float *vector,    /* vector */
      long  n)          /* size */

     /* disallocates storage for a vector of size n */

{
free(vector);
return;
}

/*--------------------------------------------------------------------------*/

void disalloc_matrix

     (float **matrix,   /* matrix */
      long  nx,         /* size in x direction */
      long  ny)         /* size in y direction */

     /* disallocates storage for matrix of size nx * ny */

{
long i;
for (i=0; i<nx; i++)
    free(matrix[i]);
free(matrix);
return;
}

/*--------------------------------------------------------------------------*/

long mylog2 

     (long n)               /* should be positive */

/*
 returns ld(n) if n is a power of 2; 
 returns -1    in all other cases. 
*/

{
long  ld;     /* ld(n) */
long  m;      /* auxiliary variable */

if (n <= 0)
   ld = -1;
else if (n == 1)
   ld = 0;
else
   {
   ld = 0;
   m  = 1; 
   do {
      m = 2 * m; 
      ld = ld + 1;
      }
   while (m < n);
   if (m > n) ld = -1;
   }

return (ld);
}

/*--------------------------------------------------------------------------*/

void FFT 

     (float    *vr,         /* real part of signal / Fourier coeff. */
      float    *vi,         /* imaginary part of signal / Fourier coeff. */
      long     n)           /* signal length, has to be power of 2 */ 

/*
 Fast Fourier Transform of a (complex) 1-D signal. 
 Based on the description in the Bronstein book.
 The signal length has to be a power of 2.
*/

{
const    float fa = sqrt(0.5);    /* frequently used renormalization factor */
long     p, q, r, j, k;           /* variables for indices and sizes */
long     nh, qq, qh;              /* variables for indices and sizes */
long     jq, jqh, jv, jn;         /* variables for indices and sizes */
long     logn;                    /* ld(n) */
long     m;                       /* auxiliary variable */
float    rh, ih, ch, chih;        /* for intermediate results */
float    *scrr, *scri, *exh;      /* auxiliary vectors */
float*   srr;                     /* point at source arrays, real part */ 
float*   sri;                     /* point at source arrays, imag. part */ 
float*   der;                     /* point at dest. arrays, real part */ 
float*   dei;                     /* point at dest. arrays, imag. part */
float*   swpp;                    /* used for pointer swapping */


/* ---- memory allocations ---- */

alloc_vector (&scrr, n);
alloc_vector (&scri, n);
alloc_vector (&exh,  n);


/* ---- initializations ----*/

/* init pointers */
srr = vr; 
sri = vi; 
der = scrr; 
dei = scri; 

/* logn := ld(n) */
m = n;
logn = -1;
while (m >= 1)
      {
      m >>= 1;              /* m = m / 2 */
      logn++;               
      }

/* init sizes */
nh = n>>1;                  /* n / 2 */ 
qh = nh>>1;                 /* n / 4 */ 

/* trigonometric values */
exh[0]  = 1.0;   exh[nh]    = 0.0;    /* exp (2 pi i 0.0 ) */
exh[qh] = 0.0;   exh[nh+qh] = 1.0;    /* exp (2 pi i 0.25) */
ch = -1.0;                            /* cos pi */
/* further trigonometric values will be computed by recursion */

/* other initializations */
qq = n; 
q  = nh; 
qh = q>>1; 
p  = 1;


/* ---- loop through all levels ----*/

for (r=logn; r>=1; r--) 
    {
    /* iterate through levels */
    if (r < logn) 
       ch = sqrt (0.5 * (1.0 + ch));    /* recursion for cosines */
    for (j=0; j<p; j++) 
        {         
        /* iterate through columns */
        jq = j * qq; 
        jqh = jq >> 1;
        if ((j&1==1) && (r<logn-1)) 
           {       
           /* recursion for exp(i*angle) */
           chih = 0.5 / ch;                    /* cosine inverse half */
           jv = jqh - q;                       
           jn = jqh + q; 
           if (jn == nh) 
              {                   
              /* use half-angle formula for exp */
              exh[jqh]    = (exh[jv] - 1.0) * chih;
              exh[jqh+nh] = exh[jv+nh] * chih;
              }
           else 
              {
              exh[jqh]    = (exh[jv]    + exh[jn]   ) * chih;
              exh[jqh+nh] = (exh[jv+nh] + exh[jn+nh]) * chih;
              }
           } /* if */
        for (k=0; k<q; k++) 
            {               
            /* iterate through rows */
            rh =  exh[jqh]    * srr[jq+k+q] + exh[jqh+nh] * sri[jq+k+q];
            ih = -exh[jqh+nh] * srr[jq+k+q] + exh[jqh]    * sri[jq+k+q];
            der[jqh+k]    = fa * (srr[jq+k] + rh);
            dei[jqh+k]    = fa * (sri[jq+k] + ih);
            der[jqh+nh+k] = fa * (srr[jq+k] - rh);
            dei[jqh+nh+k] = fa * (sri[jq+k] - ih);
            }
        } /* for j */
        
        /* swap array pointers */
        swpp = srr;   srr = der;   der = swpp;        
        swpp = sri;   sri = dei;   dei = swpp;

        /* adjust sizes */ 
        qq = q; 
        q = qh; 
        qh >>= 1; 
        p <<= 1;          
    }  /* for r */

/* copy f^ back, if ld(n) is odd */
if (logn&1 == 1)                             
   for (j=0; j<n; j++) 
       {                 
       der[j] = srr[j]; 
       dei[j] = sri[j]; 
       }


/* ---- disallocate memory ----*/

disalloc_vector (scrr, n);
disalloc_vector (scri, n);
disalloc_vector (exh,  n);

return;

} /* FFT */
  
/*--------------------------------------------------------------------------*/

void DFT  

     (float    *vr,         /* real part of signal / Fourier coeff. */
      float    *vi,         /* imaginary part of signal / Fourier coeff. */
      long     n)           /* signal length (>0) */ 


/* 
 Discrete Fourier transform of a (complex) 1-D signal. Quadratic complexity.
 Does not require powers of 2 as signal length.
*/

{
long    i, j;              /* loop variables */
float   help1, help2;      /* time savers */
float   help3, c, s;       /* time savers */
float   *fr, *fi;          /* auxiliary vectors (real / imaginary part) */
     
 
/* ---- allocate storage ---- */

alloc_vector (&fr, n);
alloc_vector (&fi, n);


/* ---- copy (vr,vi) into (fr,fi) ---- */

for (i=0; i<=n-1; i++)
    {
    fr[i] = vr[i];
    fi[i] = vi[i];
    }


/* ---- time savers ---- */

help1 = - 2.0 * 3.1415927 / (float)n;
help2 = 1.0 / sqrt ((float)n);

 
/* ---- perform DFT ---- */

for (i=0; i<=n-1; i++)
    {
    vr[i] = 0.0;
    vi[i] = 0.0;
    for (j=0; j<=n-1; j++)
        {
        help3 = help1 * i * j;
        c     = cos (help3);
        s     = sin (help3);
        vr[i] = vr[i] + fr[j] * c - fi[j] * s;
        vi[i] = vi[i] + fi[j] * c + fr[j] * s;
        }
    vr[i] = vr[i] * help2;
    vi[i] = vi[i] * help2;
    }


/* ---- disallocate storage ---- */

disalloc_vector (fr, n);
disalloc_vector (fi, n);
return;

} /* DFT */

/* ---------------------------------------------------------------------- */

void FT2D  

     (float    **ur,        /* real part of image / Fourier coeff. */
      float    **ui,        /* imaginary part of image / Fourier coeff. */
      long     nx,          /* pixel number in x direction */ 
      long     ny)          /* pixel number in y direction */ 


/* 
 Two-dimensional discrete Fourier transform of a (complex) image.
 This algorithm exploits the separability of the Fourier transform. 
 Uses FFT when the pixel numbers are powers of 2, DFT otherwise.
*/


{
long   i, j;              /* loop variables */
long   n;                 /* max (nx, ny) */
long   logn;              /* ld(n) */
float  *vr, *vi;          /* real / imaginary signal or Fourier data */


/* ---- allocate auxiliary vectors vr, vi ---- */

if (nx > ny) 
   n = nx; 
else 
   n = ny;

alloc_vector (&vr, n);
alloc_vector (&vi, n);


/* ---- transform along x direction ---- */

logn = mylog2 (nx);
for (j=0; j<=ny-1; j++)
    {
    /* write in 1-D vector */
    for (i=0; i<=nx-1; i++)
        {
        vr[i] = ur[i][j];
        vi[i] = ui[i][j];
        }

    /* apply Fourier transform */
    if (logn == -1)
       /* nx is not a power of 2; use DFT */
       DFT (vr, vi, nx); 
    else
       /* nx is a power of 2; use FFT */
       FFT (vr, vi, nx); 

    /* write back in 2-D image */
    for (i=0; i<=nx-1; i++)
        {
        ur[i][j] = vr[i];
        ui[i][j] = vi[i];
        }
    }


/* ---- transform along y direction ---- */

logn = mylog2 (ny);
for (i=0; i<=nx-1; i++)
    {
    /* write in 1-D vector */
    for (j=0; j<=ny-1; j++)
        {
        vr[j] = ur[i][j];
        vi[j] = ui[i][j];
        }

    /* apply Fourier transform */
    if (logn == -1) 
       /* ny is not a power of 2; use DFT */
       DFT (vr, vi, ny);
    else
       /* ny is a power of 2; use FFT */
       FFT (vr, vi, ny); 

    /* write back in 2-D image */
    for (j=0; j<=ny-1; j++)
        {
        ur[i][j] = vr[j];
        ui[i][j] = vi[j];
        }
    }


/* ---- disallocate storage ---- */

disalloc_vector (vr, n);
disalloc_vector (vi, n);

return;
}

/* ---------------------------------------------------------------------- */

void periodic_shift  

     (float    **u,         /* image, changed */
      long     nx,          /* pixel number in x direction */ 
      long     ny,          /* pixel number in y direction */
      long     xshift,      /* shift in x direction */ 
      long     yshift)      /* shift in y direction */ 

/*
 shifts an image u by the translation vector (xshift,yshift) 
 with 0 <= xshift <= nx-1 and 0 <= yshift <= ny-1.
*/

{
long    i, j;         /* loop variables */
float   **f;          /* auxiliary image */

/* allocate storage */
alloc_matrix (&f, nx, ny);

/* shift in x direction */
for (i=0; i<=nx-1; i++)
 for (j=0; j<=ny-1; j++)
     if (i-xshift >= 0)
        f[i][j] = u[i-xshift][j];
     else 
        f[i][j] = u[i+nx-xshift][j];

/* shift in y direction */
for (i=0; i<=nx-1; i++)
 for (j=0; j<=ny-1; j++)
     if (j-yshift >= 0)
        u[i][j] = f[i][j-yshift];
     else 
        u[i][j] = f[i][j+ny-yshift];

/* disallocate storage */
disalloc_matrix (f, nx, ny);

return;
}

/* ---------------------------------------------------------------------- */

void FT_deconv  

     (float    **ur,        /* input image, changed */
      float    **hr,        /* kernel image, changed */
      float    param,       /* parameter for filtering */
      long     nx,          /* pixel number in x direction, unchanged */ 
      long     ny)          /* pixel number in y direction, unchanged */ 
      

/*
 computes the Gaussian derivative of order (xord, yord) and writes
 it back into **ur
*/


{
long   i, j;           /* loop variables */
float  **ui;           /* imaginary image or Fourier data (image)  */
float  **hi;           /* imaginary image or Fourier data (kernel) */
float  hweight;        /* weight of kernel, for normalisation */


/* ---- allocate storage ---- */

alloc_matrix (&ui, nx, ny);
alloc_matrix (&hi, nx, ny);


/* ---- initialise imaginary part with 0 ---- */

for (j=0; j<ny; j++)
 for (i=0; i<nx; i++) {
     ui[i][j] = 0.0;
     hi[i][j] = 0.0;
 }


/* ---- compute discrete Fourier transformation ---- */

printf ("computing Fourier transformation of image\n");
FT2D (ur, ui, nx, ny);

printf ("computing Fourier transformation of kernel\n");
FT2D (hr, hi, nx, ny);

/* normalise kernel */

hweight=hr[0][0];
for (j=0; j<ny; j++)
 for (i=0; i<nx; i++) {
     hr[i][j] = hr[i][j] / hweight;
     hi[i][j] = hi[i][j] / hweight;
 }


/* ---- shift lowest frequency in the centre ----*/

/* periodic_shift (ur, nx, ny, nx/2, ny/2); */
/* periodic_shift (ui, nx, ny, nx/2, ny/2); */


/* ---- filter the Fourier coefficients ---- */

printf ("filtering the Fourier coefficients\n");
filter (ur, ui, hr, hi, param, nx, ny);

/* ---- shift lowest frequency back to the corners ----*/

/* periodic_shift (ur, nx, ny, nx-nx/2, ny-ny/2); */
/* periodic_shift (ui, nx, ny, nx-nx/2, ny-ny/2); */


/* ---- compute discrete Fourier backtransformation ---- */

printf ("computing Fourier backtransformation\n");

/* backtransformation = DFT of complex conjugated Fourier coefficients */
for (i=0; i<=nx-1; i++)
 for (j=0; j<=ny-1; j++)
     ui[i][j] = - ui[i][j];
FT2D (ur, ui, nx, ny);


/* ---- free storage ---- */

disalloc_matrix (ui, nx, ny);

return;

} /* FT_deconv */

/* ---------------------------------------------------------------------- */

void analyse

     (float   **u,         /* image, unchanged */
      long    nx,          /* pixel number in x direction */
      long    ny,          /* pixel number in x direction */
      float   *min,        /* minimum, output */
      float   *max,        /* maximum, output */
      float   *mean,       /* mean, output */ 
      float   *stdv)       /* standard deviation, output */

/*
 calculates minimum, maximum, mean and variance of an image u
*/

{
long    i, j;       /* loop variables */
float   vari;       /* variance */
float   help;       /* auxiliary variable */
double  help2;      /* auxiliary variable */

*min  = u[0][0];
*max  = u[0][0];
help2 = 0.0;
for (i=0; i<=nx-1; i++)
 for (j=0; j<=ny-1; j++)
     {
     if (u[i][j] < *min) *min = u[i][j];
     if (u[i][j] > *max) *max = u[i][j];
     help2 = help2 + (double)u[i][j];
     }
*mean = (float)help2 / (nx * ny);

vari = 0.0;
for (i=0; i<=nx-1; i++)
 for (j=0; j<=ny-1; j++)
     {
     help  = u[i][j] - *mean;
     vari = vari + help * help;
     }
vari = vari / (nx * ny);
*stdv = sqrt (vari);

return;

} /* analyse */

/* ---------------------------------------------------------------------- */

void rescale
 
     (float   **u,         /* image, changed */
      long    nx,          /* pixel number in x direction */
      long    ny,          /* pixel number in x direction */
      float   u1,          /* grey value of u */
      float   u2,          /* grey value of u */
      float   v1,          /* corresponding grey value of v */
      float   v2)          /* corresponding grey value of v */

/*
 rescales image u such that the grey value u1 is mapped into v1,
 and u2 into v2.
*/

{
long    i, j;       /* loop variables */
float   m, b;       /* parameters of the transformation */

m = (v2 - v1) / (u2 - u1);
b = v1 - m * u1;

for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     u[i][j] = m * u[i][j] + b;

return;

} /* rescale */

/* ---------------------------------------------------------------------- */

int main ()

{
char   row[300];             /* for reading data */
char   in[300];              /* input file name */
char   out[300];             /* output file name */
float  **ur;                 /* image data */
float  **hr;                 /* kernel data */
long   i, j;                 /* loop variables */
long   nx, ny;               /* image size in x, y direction */
long   hnx, hny;             /* kernel image size in x, y direction */
float  param;                /* filter parameter */
float  max, min;             /* largest, smallest grey value */
float  mean;                 /* average grey value */
float  stdv;                 /* standard deviation */
FILE   *inimage, *outimage;  /* input file, output file */
unsigned char byte;          /* for data conversion */


printf("\n");
printf("DECONVOLUTION IN THE FOURIER DOMAIN\n\n");
printf("**********************************************************\n\n");
printf("    Copyright 2004 by Joachim Weickert and Martin Welk\n");
printf("    Faculty of Mathematics and Computer Science\n");
printf("    Saarland University, Germany                      \n\n");
printf("    All rights reserved. Unauthorized usage, copying, \n");
printf("    hiring, and selling prohibited.                   \n\n");
printf("    Send bug reports to bruhn@vis.uni-stuttgart.de \n\n");
printf("**********************************************************\n\n");


/* ---- read input image (pgm format P5) ---- */

/* read image name */
printf("input image:                      ");
fgets (in, 300, stdin);
if (in[strlen(in)-1] == '\n') 
   in[strlen(in)-1] = 0;

/* open pgm file and read header */
inimage = fopen(in,"rb");
fgets (row, 300, inimage);
fgets (row, 300, inimage);
while (row[0]=='#') fgets(row, 300, inimage);
sscanf (row, "%ld %ld", &nx, &ny);
fgets (row, 300, inimage);

/* allocate storage */
alloc_matrix (&ur, nx, ny);

/* read image data */
for (j=0; j<ny; j++)
 for (i=0; i<nx; i++)
     ur[i][j] = (float) getc (inimage);
fclose(inimage);


/* ---- read kernel image (pgm format P5) ---- */

/* read image name */
printf("input kernel:                     ");
fgets (in, 300, stdin);
if (in[strlen(in)-1] == '\n') 
   in[strlen(in)-1] = 0;

/* open pgm file and read header */
inimage = fopen(in,"rb");
fgets (row, 300, inimage);
fgets (row, 300, inimage);
while (row[0]=='#') fgets(row, 300, inimage);
sscanf (row, "%ld %ld", &hnx, &hny);
fgets (row, 300, inimage);

/* allocate storage */
alloc_matrix (&hr, nx, ny);
for (j=0; j<ny; j++)
 for (i=0; i<nx; i++) 
  hr[i][j]=0.0;
   
/* read kernel data */
for (j=0; j<hny; j++)
 for (i=0; i<hnx; i++) {
   if (i<nx && j<ny) {
     hr[i][j] = (float) getc (inimage);
   }
 }
fclose(inimage);


/* ---- shift kernel to the corner ----*/

periodic_shift (hr, nx, ny, nx-hnx/2, ny-hny/2);


/* ---- read other parameters ---- */

printf("filter parameter:                 ");
fgets (row, 300, stdin);
if (row[strlen(row)-1] == '\n') 
   row[strlen(row)-1] = 0;
sscanf (row, "%f", &param);

printf("output image:                     ");
fgets(out, 300, stdin);
if (out[strlen(out)-1] == '\n') 
   out[strlen(out)-1] = 0;

printf("\n");


/* ---- compute filter in the Fourier domain ---- */

FT_deconv (ur, hr, param, nx, ny);


/* ---- analyse the outcome ---- */

analyse (ur, nx, ny, &min, &max, &mean, &stdv);
printf("\n");
printf("minimum:       %8.2f \n", min);
printf("maximum:       %8.2f \n", max);
printf("mean:          %8.2f \n", mean);
printf("st. dev.:      %8.2f \n\n", stdv);


/* ---- write output image (pgm format P5) ---- */

/* open file and write header (incl. filter parameters) */
outimage = fopen (out, "wb");
fprintf (outimage, "P5 \n");
fprintf (outimage, "# deconvolution via the Fourier domain\n");
fprintf (outimage, "# initial image:  %s\n", in);
fprintf (outimage, "# param:          %8.4f\n", param);
fprintf (outimage, "# min:            %8.2f\n", min);
fprintf (outimage, "# max:            %8.2f\n", max);
fprintf (outimage, "# mean:           %8.2f\n", mean);
fprintf (outimage, "# st. dev.:       %8.2f\n", stdv);
fprintf (outimage, "%ld %ld \n255\n", nx, ny);

/* write image data and close file */
for (j=0; j<=ny-1; j++)
 for (i=0; i<=nx-1; i++)
     {
     if (ur[i][j] < 0.0)
        byte = (unsigned char)(0.0);
     else if (ur[i][j] > 255.0)
        byte = (unsigned char)(255.0);
     else
        byte = (unsigned char)(ur[i][j]);
     fwrite (&byte, sizeof(unsigned char), 1, outimage);
     }
fclose(outimage);
printf("output image %s successfully written\n\n", out);


/* ---- disallocate storage ---- */

disalloc_matrix (hr, nx, ny);
disalloc_matrix (ur, nx, ny);
return(0);
}
