#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#define PI 3.141592654

/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                              Quantisation                                */
/*                                                                          */
/*         (Copyright Kai Hagenburg, Andres Bruhn,                          */
/*                    Markus Mainberger and Joachim Weickert, 10/2010)      */
/*                                                                          */
/*--------------------------------------------------------------------------*/

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
/*                             Uniform noise                                */

float uniform_noise(float a, float b)
{

  double n_1;
  // INSERT CODE HERE
  n_1=a + (float)rand() / RAND_MAX * (b-a);
  return n_1;

}


/*                            Gaussian noise                                */
/*--------------------------------------------------------------------------*/

float gaussian_noise(float sigma, float mu)
{
  float U,V;         /* Random Variables U,V */
  float n_1, n_2;


  /* compute random variables with normal distribution */
  // INSERT CODE HERE
  U=(float)rand() / RAND_MAX;
  V=(float)rand() / RAND_MAX;
  n_1=sqrt(-2*log(U))*cos(2*PI*V);
  n_2=sqrt(-2*log(U))*sin(2*PI*V);

  /* compute random variables with normal distribution sigma and mean mu */
  // INSERT CODE HERE
  
  //mean mu=0
  n_1= n_1*sigma;

  return n_1;
}

/*			Quantization Mean Value	     			    */
/*--------------------------------------------------------------------------*/
float Mean_Value(int intervalStart, int intervalEnd, float d)
{
float mean=0;
    while(intervalStart<=intervalEnd)
    {
	mean+=intervalStart;
	intervalStart++;
    }
mean=mean/d;
int roundedValue=mean;
mean=roundedValue;
return mean;
}

/*--------------------------------------------------------------------------*/
int main ()

{
char   row[80];              /* for reading data */
char   in[80];               /* for reading data */
char   out[80];              /* for reading data */
float  **g;                  /* original image */
float  **f;                  /* filtered image */
long   i, j;                 /* loop variables */ 
long   nx, ny;               /* image size in x, y direction */ 
FILE   *inimage, *outimage;  /* input file, output file */
int    q;                    /* number of bits used to represent a value
                                in the output image */
int    noise;                /* Noise model */
float  n_1;                  /* Noise variable */
float  sigma;                /* standard deviation */
float  mu  = 0;              /* mean */
float  mse;                  /* mean squared error (MSE) */
float  psnr;                 /* peek-signal-to-noise ratio (PSNR) */
unsigned char byte;          /* for data conversion */

printf("\n");
printf("QUANTISATION\n\n");
printf("****************************************************************\n");
printf("\n");
printf("    Copyright 2010 by Kai Hagenburg, Andres Bruhn               \n");
printf("                   and Markus Mainberger and Joachim Weickert   \n");
printf("    Faculty of Mathematics and Computer Science\n");
printf("    Saarland University, Germany\n");
printf("\n");
printf("    All rights reserved. Unauthorized usage,\n");
printf("    copying, hiring, and selling prohibited.\n");
printf("\n");
printf("    Send bug reports to\n");
printf("    bruhn@vis.uni-stuttgart.de\n");
printf("\n");
printf("****************************************************************\n\n");

/* ---- read input image (pgm format P5) ---- */

/* read image name */
printf("input image:                                            ");
gets (in);

/* open pgm file and read header */
inimage = fopen(in,"rb");
fgets (row, 80, inimage);
fgets (row, 80, inimage);
while (row[0]=='#') fgets(row, 80, inimage);
sscanf (row, "%ld %ld", &nx, &ny);
fgets (row, 80, inimage);

/* allocate storage */
alloc_matrix (&f, nx+2, ny+2);
alloc_matrix (&g, nx+2, ny+2);

/* read image data */
for (j=1; j<=ny; j++)
 for (i=1; i<=nx; i++)
     g[i][j] = (float) getc (inimage);
fclose(inimage);

/* ---- read other parameters ---- */

printf("number of bits used to represent a quantised value (q): ");
gets(row);  sscanf(row, "%i", &q);
printf("noise model, uniform (0), Gaussian (1), else no noise : ");
gets(row);  sscanf(row, "%i", &noise);
printf("noise parameter sigma:                                  ");
gets(row);  sscanf(row, "%f", &sigma);
printf("output image:                                           ");
gets(out);
printf("\n");


/* ---- quantisation ---- */
float d;


n_1 = 0;

/* compute interval size d */
//INSERT CODE HERE
d=pow(2,8-q);
int intervalStart=0;
int intervalEnd=d;
for(intervalEnd=d;intervalEnd<=256;intervalEnd+=d)
{
    for (j=1; j<=ny; j++)
	for (i=1; i<=nx; i++)
         {

        /* Switch between noise levels */

         if (noise==0)
             n_1 = uniform_noise(-sigma,sigma);
         else if (noise==1) 
             n_1 = gaussian_noise(sigma,mu);
         else 
             n_1 = 0.0;

    /* compute quantisation */
//INSERT CODE HERE
         if(g[i][j]>=intervalStart && g[i][j]<intervalEnd)

		f[i][j]=Mean_Value(intervalStart, intervalEnd-1, d)+n_1;
 	 }
intervalStart+=d;
}



/* ---- error measures ---- */

mse  = 0.0;
psnr = 0.0;

/* compute MSE und PSNR */
// INSERT CODE HERE
float substraction=0;
for (j=1; j<=ny; j++)
    for (i=1; i<=nx; i++)
    {
	substraction=f[i][j]-g[i][j];
	mse+=pow(substraction,2);
    }
mse= mse/(ny*nx);
psnr=10*log10(255*255/mse);


/* ---- error measures ---- */

printf("MSE: %f, PSNR: %f\n\n",mse, psnr );

/* ---- write output image (pgm format P5) ---- */

/* open file and write header (incl. filter parameters) */
outimage = fopen (out, "wb");
fprintf (outimage, "P5 \n");
fprintf (outimage, "# Quantisation\n");
fprintf (outimage, "# initial image:  %s\n", in);
fprintf (outimage, "# q:              %i\n", q);
fprintf (outimage, "%ld %ld \n255\n", nx, ny);

/* write image data and close file */
for (j=1; j<=ny; j++)
 for (i=1; i<=nx; i++)
     {
     if (f[i][j] < 0.0)
        byte = (unsigned char)(0.0);
     else if (f[i][j] > 255.0)
        byte = (unsigned char)(255.0);
     else
        byte = (unsigned char)(f[i][j]);
     fwrite (&byte, sizeof(unsigned char), 1, outimage);
     }
fclose(outimage);
printf("output image %s successfully written\n\n", out);


/* ---- disallocate storage ---- */

disalloc_matrix (f, nx+2, ny+2);
disalloc_matrix (g, nx+2, ny+2);

return(0);
}