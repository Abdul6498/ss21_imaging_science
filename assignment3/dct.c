#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>


/*--------------------------------------------------------------------------*/
/*                                                                          */
/*                       DISCRETE COSINE TRANSFORM                          */
/*                                                                          */
/*         (Copyright Joachim Weickert and Andres Bruhn, 11/2007)           */
/*                                                                          */
/*--------------------------------------------------------------------------*/


/*--------------------------------------------------------------------------*/

void alloc_vector

     (float **vector,   /* vector */
      long  n)          /* size */     

     /* allocates storage for vector of size n */


{
*vector = (float *) malloc (n* sizeof(float));
if (*vector == NULL)
   {
   printf("alloc_vector: not enough storage available\n");
   exit(1);
   }

return;
}

/*--------------------------------------------------------------------------*/

void disalloc_vector

     (float *vector,    /* vector */
      long  n)          /* size */  

     /* disallocates storage for vector of size n */

{
free(vector);
return;
}

/*--------------------------------------------------------------------------*/

void alloc_matrix

     (float ***matrix,  /* matrix */
      long  nx,         /* size in x-direction */
      long  ny)         /* size in y-direction */

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
      long  nx,         /* size in x-direction */
      long  ny)         /* size in y-direction */

     /* disallocates storage for matrix of size nx * ny */

{
long i;
for (i=0; i<nx; i++)
    free(matrix[i]);
free(matrix);
return;
}

/*--------------------------------------------------------------------------*/

void analyse

     (float   **u,         /* image, unchanged */
      long    nx,          /* pixel number in x-direction */
      long    ny,          /* pixel number in y-direction */
      float   *min,        /* minimum, output */
      float   *max,        /* maximum, output */
      float   *mean,       /* mean, output */
      float   *vari)       /* variance, output */

/*
 calculates minimum, maximum, mean and variance of an image u
*/

{
long    i, j;       /* loop variables */
float   help;       /* auxiliary variable */
double  help2;      /* auxiliary variable */

*min  = u[1][1];
*max  = u[1][1];
help2 = 0.0;
for (i=0; i<nx; i++)
 for (j=0; j<ny; j++)
     {
     if (u[i][j] < *min) *min = u[i][j];
     if (u[i][j] > *max) *max = u[i][j];
     help2 = help2 + (double)u[i][j];
     }
*mean = (float)help2 / (nx * ny);

*vari = 0.0;
for (i=0; i<nx; i++)
 for (j=0; j<ny; j++)
     {
     help  = u[i][j] - *mean;
     *vari = *vari + help * help;
     }
*vari = *vari / (nx * ny);

return;

} /* analyse */


/*--------------------------------------------------------------------------*/

void jpeg_multiply_block

     (float   **c_block)         /* coefficients of the DCT */                

/*
  weights 8x8 coefficient block with JPEG weighting matrix
*/

{
long    i, j;        /* loop variables */

c_block[0][0] *= 10;
c_block[0][1] *= 15;
c_block[0][2] *= 25;
c_block[0][3] *= 37;
c_block[0][4] *= 51;
c_block[0][5] *= 66;
c_block[0][6] *= 82;
c_block[0][7] *= 100;

c_block[1][0] *= 15;
c_block[1][1] *= 19;
c_block[1][2] *= 28;
c_block[1][3] *= 39;
c_block[1][4] *= 52;
c_block[1][5] *= 67;
c_block[1][6] *= 83;
c_block[1][7] *= 101;

c_block[2][0] *= 25;
c_block[2][1] *= 28;
c_block[2][2] *= 35;
c_block[2][3] *= 45;
c_block[2][4] *= 58;
c_block[2][5] *= 72;
c_block[2][6] *= 88;
c_block[2][7] *= 105;

c_block[3][0] *= 37;
c_block[3][1] *= 39;
c_block[3][2] *= 45;
c_block[3][3] *= 54;
c_block[3][4] *= 66;
c_block[3][5] *= 79;
c_block[3][6] *= 94;
c_block[3][7] *= 111;

c_block[4][0] *= 51;
c_block[4][1] *= 52;
c_block[4][2] *= 58;
c_block[4][3] *= 66;
c_block[4][4] *= 76;
c_block[4][5] *= 89;
c_block[4][6] *= 103;
c_block[4][7] *= 119;

c_block[5][0] *= 66;
c_block[5][1] *= 67;
c_block[5][2] *= 72;
c_block[5][3] *= 79;
c_block[5][4] *= 89;
c_block[5][5] *= 101;
c_block[5][6] *= 114;
c_block[5][7] *= 130;

c_block[6][0] *= 82;
c_block[6][1] *= 83;
c_block[6][2] *= 88;
c_block[6][3] *= 94;
c_block[6][4] *= 103;
c_block[6][5] *= 114;
c_block[6][6] *= 127;
c_block[6][7] *= 142;

c_block[7][0] *= 100;
c_block[7][1] *= 101;
c_block[7][2] *= 105;
c_block[7][3] *= 111;
c_block[7][4] *= 119;
c_block[7][5] *= 130;
c_block[7][6] *= 142;
c_block[7][7] *= 156;

return;
}


/*--------------------------------------------------------------------------*/

void equal_multiply_block

     (float   **c_block,      /* coefficients of the DCT */     
      long    factor)                   

/*
  multiplies an 8x8 coefficient block with a given factor 
*/

{
long    i, j;        /* loop variables */

for(i=0;i<=7;i++)
  for(j=0;j<=7;j++)
    c_block[i][j] *= factor;

return;
}


/*--------------------------------------------------------------------------*/

void jpeg_divide_block

     (float   **c_block)         /* coefficients of the DCT */             

/*
  weights 8x8 coefficient block with inverse JPEG weighting matrix
*/

{
long    i, j;        /* loop variables */

c_block[0][0] /= 10;
c_block[0][1] /= 15;
c_block[0][2] /= 25;
c_block[0][3] /= 37;
c_block[0][4] /= 51;
c_block[0][5] /= 66;
c_block[0][6] /= 82;
c_block[0][7] /= 100;

c_block[1][0] /= 15;
c_block[1][1] /= 19;
c_block[1][2] /= 28;
c_block[1][3] /= 39;
c_block[1][4] /= 52;
c_block[1][5] /= 67;
c_block[1][6] /= 83;
c_block[1][7] /= 101;

c_block[2][0] /= 25;
c_block[2][1] /= 28;
c_block[2][2] /= 35;
c_block[2][3] /= 45;
c_block[2][4] /= 58;
c_block[2][5] /= 72;
c_block[2][6] /= 88;
c_block[2][7] /= 105;

c_block[3][0] /= 37;
c_block[3][1] /= 39;
c_block[3][2] /= 45;
c_block[3][3] /= 54;
c_block[3][4] /= 66;
c_block[3][5] /= 79;
c_block[3][6] /= 94;
c_block[3][7] /= 111;

c_block[4][0] /= 51;
c_block[4][1] /= 52;
c_block[4][2] /= 58;
c_block[4][3] /= 66;
c_block[4][4] /= 76;
c_block[4][5] /= 89;
c_block[4][6] /= 103;
c_block[4][7] /= 119;

c_block[5][0] /= 66;
c_block[5][1] /= 67;
c_block[5][2] /= 72;
c_block[5][3] /= 79;
c_block[5][4] /= 89;
c_block[5][5] /= 101;
c_block[5][6] /= 114;
c_block[5][7] /= 130;

c_block[6][0] /= 82;
c_block[6][1] /= 83;
c_block[6][2] /= 88;
c_block[6][3] /= 94;
c_block[6][4] /= 103;
c_block[6][5] /= 114;
c_block[6][6] /= 127;
c_block[6][7] /= 142;

c_block[7][0] /= 100;
c_block[7][1] /= 101;
c_block[7][2] /= 105;
c_block[7][3] /= 111;
c_block[7][4] /= 119;
c_block[7][5] /= 130;
c_block[7][6] /= 142;
c_block[7][7] /= 156;

return;
}


/*--------------------------------------------------------------------------*/

void round_block_coeff

     (float   **c_block)         /* coefficients of the DCT */    

/*
  round entries of a 8x8 coefficient block
*/

{
long    i, j;        /* loop variables */

for(i=0;i<=7;i++)
  for (j=0;j<=7;j++)
    {    
      c_block[i][j] = rintf(c_block[i][j]);          
    }

return;
}


/*--------------------------------------------------------------------------*/

void equal_divide_block

     (float   **c_block,      /* coefficients of the DCT */     
      long    factor)                   

/*
  divides an 8x8 coefficient block by a given factor 
*/

{
long    i, j;        /* loop variables */

for(i=0;i<=7;i++)
  for(j=0;j<=7;j++)
    c_block[i][j] /= factor;

return;
}

/*--------------------------------------------------------------------------*/

void DCT_2d

     (float   **u,         /* image, unchanged */
      float   **c,         /* coefficients of the DCT */        
      long    nx,          /* pixel number in x-direction */
      long    ny)          /* pixel number in y-direction */      

/*
  calulates DCT of input image
*/

{
long    i, j, m, p;  /* loop variables */
float   nx_1;        /* time saver */
float   ny_1;        /* time saver */
float   pi;          /* variable pi */
float   **tmp;       /* temporary image */
float   *cx,*cy;     /* arrays for coefficients */


/* ---- derive pi ---- */

pi = 2.0 * asinf(1.0);


/* ---- set time savers ---- */

nx_1 = pi / (2.0 * nx);
ny_1 = pi / (2.0 * ny);


/* ---- allocate memory ---- */

alloc_matrix (&tmp, nx, ny);
alloc_vector (&cx, nx);
alloc_vector (&cy, ny);


/* ---- set up coefficients ---- */

cx[0] = sqrt(1.0/nx);
for (p=1;p<nx;p++)
  cx[p] = sqrt(2.0/nx);

cy[0] = sqrt(1.0/ny);
for (p=1;p<ny;p++)
  cy[p] = sqrt(2.0/ny);


/* ---- DCT in y-direction ---- */ 

for (i=0; i<nx; i++)
  for (p=0; p<ny; p++)
    {
      /*
	SUPPLEMENT CODE HERE
      */
    float tmp_c = 0;
    for (m=0; m<ny; m++)
      tmp_c += u[m][i] * cy[p] * cos(ny_1 * (2*m+1) * p);
    tmp[p][i] = tmp_c;
    }


/* ---- DCT in x-direction ---- */ 

for (p=0; p<nx; p++)
  for (j=0; j<ny; j++)
    {
      /*
	SUPPLEMENT CODE HERE
      */
    float tmp_c = 0;
    for (m=0; m<nx; m++)
      tmp_c += tmp[j][m] * cx[p] * cos(nx_1 * (2*m+1) * p);
    c[j][p] = tmp_c;            
    }


/* ---- free memory ---- */

disalloc_matrix (tmp, nx, ny);
disalloc_vector (cx, nx);
disalloc_vector (cy, ny);

return;
}


/*--------------------------------------------------------------------------*/

void IDCT_2d

     (float   **u,         /* image, unchanged */
      float   **c,         /* coefficients of the DCT */        
      long    nx,          /* pixel number in x-direction */
      long    ny)          /* pixel number in y-direction */      

/*
  calulates inverse DCT of input image
*/

{
long    i, j, m, p;  /* loop variables */
float   nx_1;        /* time saver */
float   ny_1;        /* time saver */
float   pi;          /* variable pi */
float   **tmp;       /* temporary image */
float   *cx,*cy;     /* arrays for coefficients */


/* ---- derive pi ---- */

pi = 2.0 * asinf(1.0);


/* ---- set time savers ---- */

nx_1 = pi / (2.0 * nx);
ny_1 = pi / (2.0 * ny);


/* ---- allocate memory ---- */

alloc_matrix (&tmp, nx, ny);
alloc_vector (&cx, nx);
alloc_vector (&cy, ny);


/* ---- set up coefficients ---- */

cx[0] = sqrt(1.0/nx);
for (m=1; m<nx; m++)
    cx[m] = sqrt (2.0 / nx);

cy[0] = sqrt(1.0/ny);
for (m=1; m<ny; m++)
    cy[m] = sqrt (2.0 / ny);


/* ---- DCT in y-direction ---- */ 

for (i=0; i<nx; i++)
  for (m=0; m<ny; m++)
    {
      /*
	SUPPLEMENT CODE HERE
      */
    float tmp_f = 0;
    for(p=0; p<ny; p++)
      tmp_f += c[p][i] * cy[p] * cos(ny_1 * (2*m+1) * p);
    tmp[m][i] = tmp_f;      
    }
 

/* ---- DCT in x-direction ---- */ 

for (m=0; m<nx; m++)
  for (j=0; j<ny; j++)
    {
      /*
	SUPPLEMENT CODE HERE
      */
    float tmp_f = 0;
    for(p=0; p<nx; p++)
      tmp_f += tmp[j][p] * cx[p] * cos(nx_1 * (2*m+1) * p);
    u[j][m] = tmp_f;
    }


/* ---- free memory ---- */

disalloc_matrix (tmp, nx, ny);
disalloc_vector (cx, nx);
disalloc_vector (cy, ny);

return;
}


/*--------------------------------------------------------------------------*/

void remove_freq_2d

     (float   **c,         /* in and out: coefficients of the DCT */        
      long    nx,          /* pixel number in x-direction */
      long    ny)          /* pixel number in y-direction */      

/*
 remove frequencies
*/

{
long    i, j;        /* loop variables */

/* sets frequencies to zero */ 
for (i=0;i<nx;i++)
  for (j=0;j<ny;j++)
    if ((i >= nx/sqrt(10)) || (j >= ny/sqrt(10)))
      c[i][j]=0;
 
return;
}


/*--------------------------------------------------------------------------*/

void blockwise_DCT_2d

     (float   **u,         /* in: image */
      float   **c,         /* out: coefficients of the DCT */        
      long    nx,          /* pixel number in x-direction */
      long    ny)          /* pixel number in y-direction */      

/*
 calculates DCT in image blocks of size 8x8
*/

{
long    i, j;       /* loop variables */
long    k, l;    
float   **u_block;    /* 8x8 blocks */
float   **c_block;

/* allocate memory for 8x8 block */
alloc_matrix(&u_block, 8, 8);
alloc_matrix(&c_block, 8, 8);

for (i=0;i<nx;i+=8)
 for (j=0;j<ny;j+=8)
     {
     /* copy 8x8 block */
     for (k=0;k<=7;k++)
      for (l=0;l<=7;l++)
	  u_block[k][l] = u[i+k][j+l];
      
     /* DCT of 8x8 block */
     DCT_2d(u_block, c_block, 8, 8);
      
     /* copy back coefficients */
     for (k=0;k<=7;k++)
      for (l=0;l<=7;l++)
	  c[i+k][j+l] = c_block[k][l];           
     }

/* free memory for 8x8 block */
disalloc_matrix (u_block, 8, 8);
disalloc_matrix (c_block, 8, 8);

return;
}

/*--------------------------------------------------------------------------*/

void blockwise_IDCT_2d

     (float   **u,         /* out: image */
      float   **c,         /* in: coefficients of the DCT */        
      long    nx,          /* pixel number in x-direction */
      long    ny)          /* pixel number in y-direction */      
/*
 calculates inverse DCT in image blocks of size 8x8
*/

{
long    i, j;       /* loop variables */
long    k, l;    
float   **u_block;    /* 8x8 blocks */
float   **c_block;

/* allocate memory for 8x8 block */
alloc_matrix(&u_block, 8, 8);
alloc_matrix(&c_block, 8, 8);

for (i=0;i<nx;i+=8)
 for (j=0;j<ny;j+=8)
     {
     /* copy 8x8 block */
     for (k=0;k<=7;k++)
      for (l=0;l<=7;l++)
	  c_block[k][l] = c[i+k][j+l];
      
     /* inverse DCT of 8x8 block */
     IDCT_2d(u_block, c_block, 8, 8);
      
     /* copy back coefficients */
     for (k=0;k<=7;k++)
      for (l=0;l<=7;l++)
	  u[i+k][j+l] = u_block[k][l];           
     }

/* free memory for 8x8 block */
disalloc_matrix (u_block, 8, 8);
disalloc_matrix (c_block, 8, 8);

return;
}


/*--------------------------------------------------------------------------*/

void blockwise_remove_freq_2d

     (float   **c,         /* in and out: coefficients of the DCT */        
      long    nx,          /* pixel number in x-direction */
      long    ny)          /* pixel number in y-direction */      
/*
 removes frequencies within 8x8 block
*/

{
long    i, j;         /* loop variables */
long    k, l;         /* loop variables */ 
float   **c_block;    /* 8x8 block */


/* allocate memory for 8x8 block */
alloc_matrix(&c_block, 8, 8);

/* scale coefficients blockwise */
for (i=0;i<nx;i+=8)
 for (j=0;j<ny;j+=8)
     {
       /* copy 8x8 block */
       for (k=0;k<=7;k++)
	 for (l=0;l<=7;l++)
           c_block[k][l] = c[i+k][j+l];
 
       /* set frequencies to zero */
       remove_freq_2d(c_block,8,8);      

       /* copy back coefficients */
       for (k=0;k<=7;k++)
	 for (l=0;l<=7;l++)
	   c[i+k][j+l] = c_block[k][l];       
     }

/* free memory for 8x8 block */
disalloc_matrix (c_block, 8, 8);

return;
}


/*--------------------------------------------------------------------------*/

void blockwise_quantisation_jpeg_2d

     (float   **c,         /* in and out: coefficients of the DCT */        
      long    nx,          /* pixel number in x-direction */
      long    ny)          /* pixel number in y-direction */      
/*
 quantises coefficients of 8x8 block
*/

{
long    i, j;         /* loop variables */
long    k, l;         /* loop variables */ 
float   **c_block;    /* 8x8 block */


/* allocate memory for 8x8 block */
alloc_matrix(&c_block, 8, 8);

/* scale coefficients blockwise */
for (i=0;i<nx;i+=8)
 for (j=0;j<ny;j+=8)
     {
       /* copy 8x8 block */
       for (k=0;k<=7;k++)
	 for (l=0;l<=7;l++)
           c_block[k][l] = c[i+k][j+l];
       
       /* scale coefficients of 8x8 block */
       jpeg_divide_block(c_block);       
 
       /* round coefficients */
       round_block_coeff(c_block);
      
       /* rescale coefficients of 8x8 block */
       jpeg_multiply_block(c_block);    

       /* copy back coefficients */
       for (k=0;k<=7;k++)
	 for (l=0;l<=7;l++)
	   c[i+k][j+l] = c_block[k][l];       
     }

/* free memory for 8x8 block */
disalloc_matrix (c_block, 8, 8);

return;
}


/*--------------------------------------------------------------------------*/

void blockwise_quantisation_equal_2d

     (float   **c,         /* in and out: coefficients of the DCT */        
      long    nx,          /* pixel number in x-direction */
      long    ny)          /* pixel number in y-direction */      
/*
 quantises coefficients of 8x8 block
*/

{
long    i, j;         /* loop variables */
long    k, l;         /* loop variables */ 
long    count;        /* counter */
float   **c_block;    /* 8x8 block */


/* allocate memory for 8x8 block */
alloc_matrix(&c_block, 8, 8);

/* scale coefficients blockwise */
for (i=0;i<nx;i+=8)
 for (j=0;j<ny;j+=8)
     {
       /* copy 8x8 block */
       for (k=0;k<=7;k++)
	 for (l=0;l<=7;l++)
           c_block[k][l] = c[i+k][j+l];
       
       /* scale coefficients of 8x8 block */     
       equal_divide_block(c_block,40);
 
       /* round coefficients */
       round_block_coeff(c_block);
      
       /* rescale coefficients of 8x8 block */     
       equal_multiply_block(c_block,40);

       /* copy back coefficients */
       for (k=0;k<=7;k++)
	 for (l=0;l<=7;l++)
	   c[i+k][j+l] = c_block[k][l];       
     }

/* free memory for 8x8 block */
disalloc_matrix (c_block, 8, 8);

return;
}


/*--------------------------------------------------------------------------*/

int main ()

{
char   row[80];              /* for reading data */
char   in[80];               /* for reading data */
char   out1[80];             /* for writing data */
char   out2[80];             /* for writing data */
float  **u;                  /* image */
float  **c;                  /* DCT coefficients */
long   i, j;                 /* loop variables */ 
long   nx, ny;               /* image size in x, y direction */ 
long   flag;                 /* processing flag */
long   counter;              /* count variable */
FILE   *inimage, *outimage;  /* input file, output file */
float  sigma;                /* standard deviation */
float  max, min;             /* largest, smallest grey value */
float  mean;                 /* average grey value */
float  vari;                 /* variance */
unsigned char byte;          /* for data conversion */

printf("\n");
printf("DISCRETE COSINE TRANSFORM\n\n");
printf("*************************************************\n");
printf("\n");
printf("    Copyright 2007 by\n");
printf("    Andres Bruhn and Joachim Weickert\n");
printf("    Faculty of Mathematics and Computer Science\n");
printf("    Saarland University, Germany\n");
printf("\n");
printf("    All rights reserved. Unauthorized usage,\n");
printf("    copying, hiring, and selling prohibited.\n");
printf("\n");
printf("    Send bug reports to\n");
printf("    bruhn@vis.uni-stuttgart.de\n");
printf("\n");
printf("*************************************************\n\n");

/* ---- read input image (pgm format P5) ---- */

/* read image name */
printf("input image:                          ");
gets (in);

/* open pgm file and read header */
inimage = fopen(in,"rb");
fgets (row, 80, inimage);
fgets (row, 80, inimage);
while (row[0]=='#') fgets(row, 80, inimage);
sscanf (row, "%ld %ld", &nx, &ny);
fgets (row, 80, inimage);


/* check if image can be devided in blocks of size 8x8 */
if ((nx%8!=0) || (ny%8!=0))
   {
   printf("\n\n Image size does not allow decomposition! \n\n");
   return(0);
   }

/* allocate storage */
alloc_matrix (&u, nx+1, ny+1);
alloc_matrix (&c, nx+1, ny+1);

/* read image data */
for (j=0; j<ny; j++)
 for (i=0; i<nx; i++)
     u[i][j] = (float) getc (inimage);
fclose(inimage);


/* ---- analyse image ---- */

analyse (u, nx, ny, &min, &max, &mean, &vari);
printf("\n");
printf("minimum:       %8.2f \n", min);
printf("maximum:       %8.2f \n", max);
printf("mean:          %8.2f \n", mean);
printf("std. dev.:     %8.2f \n\n", sqrt(vari));


/* ---- read other parameters ---- */

printf("spectrum image:                       ");
gets(out1);
printf("output image:                         ");
gets(out2);
printf("\n");
printf("\n");
printf("\nyou have the following options:"); 
printf("\n");
printf("\n   (1) DCT/IDCT of the whole image   ");
printf("\n   (2) DCT/IDCT of 8x8 blocks        ");
printf("\n   (3) DCT/IDCT of the whole image   ");
printf("\n       with removal of frequencies");
printf("\n   (4) DCT/IDCT of 8x8 blocks        ");
printf("\n       with removal of frequencies");
printf("\n   (5) DCT/IDCT of 8x8 blocks        ");
printf("\n       with equal quantisation       ");
printf("\n   (6) DCT/IDCT of 8x8 blocks        ");
printf("\n       with JPEG quantisation        ");
printf("\n");
 printf("\nchose the processing mode (1-6):     ");
gets(row);  sscanf(row, "%ld", &flag);
printf("\n\n");


/* ---- process image ---- */

switch(flag) {
 case 1 :
   /* perform DCT and IDCT for the whole image */
   DCT_2d (u, c, nx, ny);   
   IDCT_2d (u, c, nx, ny);   
   break;
 case 2 :      
   /* perform DCT and IDCT in 8x8 blocks */
   blockwise_DCT_2d (u, c, nx, ny);      
   blockwise_IDCT_2d (u, c, nx, ny);
   break;
 case 3 :
   /* perform DCT and IDCT for the whole image */
   /* remove frequencies */   
   DCT_2d (u, c, nx, ny);   
   remove_freq_2d(c, nx, ny);
   IDCT_2d (u, c, nx, ny);   
   break;
 case 4 :      
   /* perform DCT and IDCT in 8x8 blocks */
   /* remove frequencies */
   blockwise_DCT_2d (u, c, nx, ny);      
   blockwise_remove_freq_2d(c, nx, ny);
   blockwise_IDCT_2d (u, c, nx, ny);
   break;
 case 5 :
   /* perform DCT and IDCT in 8x8 blocks */
   /* and use equal quantisation */
   blockwise_DCT_2d (u, c, nx, ny);    
   blockwise_quantisation_equal_2d(c, nx, ny);
   blockwise_IDCT_2d (u, c, nx, ny);
   break;
 case 6 :
   /* perform DCT and IDCT in 8x8 blocks */
   /* and use JPEG quantisation */
   blockwise_DCT_2d (u, c, nx, ny);     
   blockwise_quantisation_jpeg_2d(c, nx, ny);
   blockwise_IDCT_2d (u, c, nx, ny);
   break;
 default : 
   printf("option (%ld) not available! \n\n\n",flag);
   return(0);
}


/* count zero entries in the DCT coefficients */
counter = 0;
 
for (j=0; j<ny; j++)
  for (i=0; i<nx; i++)
    if (c[i][j]==0) 
      counter++;

printf("percentage of zero coefficients %f\n\n",
       (float)counter/(nx*ny));

   
/* ---- compute logarithmised spectrum of c ---- */

for (j=0; j<ny; j++)
  for (i=0; i<nx; i++)
    c[i][j]=logf(1.0f + fabsf(c[i][j]));


/* ---- normalize spectrum of c ---- */

analyse (c, nx, ny, &min, &max, &mean, &vari);


if (max!=0)
  {
    for (j=0; j<ny; j++)
      for (i=0; i<nx; i++)
	c[i][j]=c[i][j]*255.0/max;
  }


/* ---- analyse processed image ---- */

analyse (u, nx, ny, &min, &max, &mean, &vari);
printf("minimum:       %8.2f \n", min);
printf("maximum:       %8.2f \n", max);
printf("mean:          %8.2f \n", mean);
printf("std. dev.:     %8.2f \n\n", sqrt(vari));


/* ---- write output image 1 (pgm format P5) ---- */

/* open file and write header (incl. filter parameters) */
outimage = fopen (out1, "wb");
fprintf (outimage, "P5 \n");
fprintf (outimage, "# Discrete Cosine Transform (spectrum)\n");
fprintf (outimage, "# initial image:  %s\n", in);
fprintf (outimage, "# menu option:  %ld\n", flag);
fprintf (outimage, "%ld %ld \n255\n", nx, ny);

/* write image data and close file */
for (j=0; j<ny; j++)
  for (i=0; i<nx; i++)
     {
     if (c[i][j] < 0.0)
        byte = (unsigned char)(0.0);
     else if (c[i][j] > 255.0)
        byte = (unsigned char)(255.0);
     else
        byte = (unsigned char)(c[i][j]);
     fwrite (&byte, sizeof(unsigned char), 1, outimage);
     }
fclose(outimage);
printf("output image %s successfully written\n", out1);



/* ---- write output image 2 (pgm format P5) ---- */

/* open file and write header (incl. filter parameters) */
outimage = fopen (out2, "wb");
fprintf (outimage, "P5 \n");
fprintf (outimage, "# Discrete Cosine Transform (image)\n");
fprintf (outimage, "# initial image:  %s\n", in);
fprintf (outimage, "# menu option:  %ld\n", flag);
fprintf (outimage, "%ld %ld \n255\n", nx, ny);

/* write image data and close file */
for (j=0; j<ny; j++)
  for (i=0; i<nx; i++)
     {
     if (u[i][j] < 0.0)
        byte = (unsigned char)(0.0);
     else if (u[i][j] > 255.0)
        byte = (unsigned char)(255.0);
     else
        byte = (unsigned char)(u[i][j]);
     fwrite (&byte, sizeof(unsigned char), 1, outimage);
     }
fclose(outimage);
printf("output image %s successfully written\n\n", out2);



/* ---- disallocate storage ---- */

disalloc_matrix (u, nx+1, ny+1);
disalloc_matrix (c, nx+1, ny+1);

return(0);
}
