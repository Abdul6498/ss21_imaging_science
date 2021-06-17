#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                         POINT TRANSFORMATIONS                            */
/*                                                                          */
/*                 (Copyright by Joachim Weickert, 12/1999)                 */
/*                 (alterations by Kai Hagenburg, 11/2010)                  */
/*                                                                          */
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/

void alloc_vector

     (double **vector,   /* vector */
      long   n)          /* size */

     /* allocates storage for a vector of size n */


{
*vector = (double *) malloc (n * sizeof(double));
if (*vector == NULL)
   {
   printf("alloc_vector: not enough storage available\n");
   exit(1);
   }
return;
}

/*--------------------------------------------------------------------------*/

void alloc_matrix

     (double ***matrix,  /* matrix */
      long   nx,         /* size in x direction */
      long   ny)         /* size in y direction */

     /* allocates storage for matrix of size nx * ny */


{
long i;

*matrix = (double **) malloc (nx * sizeof(double *));
if (*matrix == NULL)
   {
   printf("alloc_matrix: not enough storage available\n");
   exit(1);
   }
for (i=0; i<nx; i++)
    {
    (*matrix)[i] = (double *) malloc (ny * sizeof(double));
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

     (double *vector,    /* vector */
      long   n)          /* size */

     /* disallocates storage for a vector of size n */

{
free(vector);
return;
}

/*--------------------------------------------------------------------------*/

void disalloc_matrix

     (double **matrix,   /* matrix */
      long   nx,         /* size in x direction */
      long   ny)         /* size in y direction */

     /* disallocates storage for matrix of size nx * ny */

{
long i;
for (i=0; i<nx; i++)
    free(matrix[i]);
free(matrix);
return;
}

/*--------------------------------------------------------------------------*/

void rescale 

     (double  **u,        /* input image, range [0,255] */
      long    nx,         /* size in x direction */
      long    ny,         /* size in y direction */
      double  a,          /* smallest transformed grey level */
      double  b,          /* largest transformed grey level */
      double  *g)         /* transformed grey levels */

/* 
 affine rescaling of the grey values of u such that 
 min(u) -> a, and max(u) -> b. 
*/

{
long    i, j, k;    /* loop variables */
double  min, max;   /* extrema of u */
double  factor;     /* time saver */

/* determine extrema of u */

/*
 INSERT CODE HERE
*/
min= 0;
max=255;
for(i=min;i<=max; i++)
{
  g[i]= a*i+b;
}

/* rescale */

/*
 INSERT CODE HERE
*/
for(j=0;j<nx;j++)
{
  for(k=0;k<ny;k++)
  {
        u[i][j] = g[(long)(u[i][j])];
  }

}
return;
}

/*--------------------------------------------------------------------------*/
void gamma_correct

     (double  gamma,      /* gamma correction factor */
      double  *g)         /* transformed grey levels */

/* 
 applies gamma correction to the 256 grey levels that may appear
 in byte wise coded images
*/

{
long  k;   /* loop variable */

/*
 INSERT CODE HERE
*/
for(k=0;k<256;k++)
{
	g[k]= 255 * pow(k/255.0, 1/gamma);
   printf("%f",g[k]);
}
return;
}


/*--------------------------------------------------------------------------*/

int main ()

{
char   row[80];              /* for reading data */
char   in[80];               /* for reading data */
char   out[80];              /* for reading data */
double **u;                  /* image */
double *g;                   /* grey level mapping */
long   i, j;                 /* loop variables */ 
long   nx, ny;               /* image size in x, y direction */ 
long   transform;            /* type of point transformation */
double a, b;                 /* rescaling bounds */
double gamma;                /* gamma correction factor */
FILE   *inimage, *outimage;  /* input file, output file */
unsigned char byte;          /* for data conversion */

printf("\n");
printf("POINT TRANSFORMATIONS\n\n");
printf("*************************************************\n\n");
printf("    Copyright 2002 by Joachim Weickert       \n");
printf("    Faculty of Mathematics and Computer Science\n");
printf("    Saarland University, Germany              \n\n");
printf("    All rights reserved. Unauthorized usage, \n");
printf("    copying, hiring, and selling prohibited. \n\n");
printf("    Send bug reports to                      \n");
printf("    bruhn@vis.uni-stuttgart.de             \n\n");
printf("*************************************************\n\n");


/* ---- read input image (pgm format P5) ---- */

/* read image name */
printf("input image:                      ");
gets (in);

/* open pgm file and read header */
inimage = fopen(in,"rb");
fgets (row, 80, inimage);
fgets (row, 80, inimage);
while (row[0]=='#') fgets (row, 80, inimage);
sscanf (row, "%ld %ld", &nx, &ny);
fgets (row, 80, inimage);

/* allocate storage */
alloc_matrix (&u, nx+2, ny+2);

/* read image data */
for (j=1; j<=ny; j++)
 for (i=1; i<=nx; i++)
     u[i][j] = (double) getc (inimage);
fclose(inimage);


/* ---- read other parameters ---- */

printf("available point transformations:\n");
printf("(0) affine rescaling\n");
printf("(1) gamma correction\n");
printf("your choice:                      ");
gets(row);  sscanf(row, "%ld", &transform);
if (transform == 0)
   {
   printf("smallest grey value:              ");
   gets(row);  sscanf(row, "%lf", &a);
   printf("largest  grey value:              ");
   gets(row);  sscanf(row, "%lf", &b);
   }
if (transform == 1)
   {
   printf("gamma correction factor:          ");
   gets(row);  sscanf(row, "%lf", &gamma);
   }
printf("output image:                     ");
gets(out);
printf("\n");


/* ---- greyscale transformation ---- */

/* allocate storage for greyscale transformation vector */
alloc_vector (&g, 256);

/* calculate greyscale transformation vector */
if (transform == 0) 
   rescale (u, nx, ny, a, b, g);
if (transform == 1) 
   gamma_correct (gamma, g);

/* apply greyscale transformation to the image */
for (i=1; i<=nx; i++)
 for (j=1; j<=ny; j++)
     u[i][j] = g[(long)(u[i][j])];
 

/* ---- write output image (pgm format P5) ---- */

/* open file and write header (incl. filter parameters) */
outimage = fopen (out, "wb");
fprintf (outimage, "P5 \n");
fprintf (outimage, "# input image:  %s\n", in);
if (transform == 0)
   {
   fprintf (outimage, "# affine rescaling\n");
   fprintf (outimage, "# a: %6.2lf\n", a);
   fprintf (outimage, "# b: %6.2lf\n", b);
   }
if (transform == 1)
   {
   fprintf (outimage, "# gamma correction\n");
   fprintf (outimage, "# gamma: %4.2lf\n", gamma);
   }
fprintf (outimage, "%ld %ld \n255\n", nx, ny);

/* write image data and close file */
for (j=1; j<=ny; j++)
 for (i=1; i<=nx; i++)
     {
     byte = (unsigned char)(u[i][j]);
     fwrite (&byte, sizeof(unsigned char), 1, outimage);
     }
fclose(outimage);
printf("output image %s successfully written\n\n", out);


/* ---- disallocate storage ---- */

disalloc_vector (g, 256);
disalloc_matrix (u, nx+2, ny+2);
return(0);
}
