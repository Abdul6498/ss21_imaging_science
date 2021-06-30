/******************************************************************************/
/*                                                                            */
/*       Copyright 11/2008 by Markus Mainberger and Joachim Weickert          */
/*     Faculty of Mathematics and Computer Science, Saarland University,      */
/*                           Saarbruecken, Germany.                           */
/*                                                                            */
/*                Canny edge detector for grey value images                   */
/*                                                                            */
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define EPSILON 10E-7


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

void presmooth 

     (float    **f,       /* input: original image */
      float    **u,       /* output: smoothed */
      long     nx,        /* image dimension in x direction */ 
      long     ny,        /* image dimension in y direction */ 
      float    sigma)     /* standard deviation of Gaussian */


      
/* 
 Gaussian convolution. Copyright by Joachim Weickert 5/2000
*/


{
long    i, j, p;              /* loop variables */
long    length;               /* convolution vector: 0..length */
float   sum;                  /* for summing up */
float   *conv;                /* convolution vector */
float   *help;                /* row or column with dummy boundaries */
      
     
/* ------------------------ diffusion in x direction -------------------- */

/* calculate length of convolution vector */
length = (long)(3.0 * sigma) + 1;
if (length > nx)
   {
   printf("gauss_conv: sigma too large \n"); 
   exit(0);
   }

/* allocate storage for convolution vector */
alloc_vector (&conv, length+1);

/* calculate entries of convolution vector */
for (i=0; i<=length; i++)
    conv[i] = 1 / (sigma * sqrt(2.0 * 3.1415926)) 
              * exp (- (i * i) / (2.0 * sigma * sigma));

/* normalization */
sum = conv[0];
for (i=1; i<=length; i++)
    sum = sum + 2.0 * conv[i];
for (i=0; i<=length; i++)
    conv[i] = conv[i] / sum;

/* allocate storage for a row */
alloc_vector (&help, nx+length+length);

for (j=1; j<=ny; j++)
    {
    /* copy in row vector */
    for (i=1; i<=nx; i++)
        help[i+length-1] = f[i][j];

    /* assign reflecting boundary conditions */       
    for (p=1; p<=length; p++)
      {
	help[length-p]      = help[length+p-1];
	help[nx+length-1+p] = help[nx+length-p];
      }
    
    /* convolution step */
    for (i=length; i<=nx+length-1; i++)
        {
        /* calculate convolution */
        sum = conv[0] * help[i];
        for (p=1; p<=length; p++)
            sum = sum + conv[p] * (help[i+p] + help[i-p]);
        /* write back */
        f[i-length+1][j] = sum;
        }
    } /* for j */

/* disallocate storage for a row */
disalloc_vector (help, nx+length+length);

/* disallocate convolution vector */
disalloc_vector (conv, length + 1);


/* ------------------------ diffusion in y direction -------------------- */

/* calculate length of convolution vector */
length = (long)(3.0 * sigma) + 1;
if (length > ny)
   {
   printf("gauss_conv: sigma too large \n"); 
   exit(0);
   }

/* allocate storage for convolution vector */
alloc_vector (&conv, length + 1);

/* calculate entries of convolution vector */
for (j=0; j<=length; j++)
    conv[j] = 1 / (sigma * sqrt(2.0 * 3.1415927)) 
              * exp (- (j * j) / (2.0 * sigma * sigma));

/* normalization */
sum = conv[0];
for (j=1; j<=length; j++)
    sum = sum + 2.0 * conv[j];
for (j=0; j<=length; j++)
    conv[j] = conv[j] / sum;

/* allocate storage for a row */
alloc_vector (&help, ny+length+length);

for (i=1; i<=nx; i++)
    {
    /* copy in column vector */
    for (j=1; j<=ny; j++)
        help[j+length-1] = f[i][j];

    /* assign reflecting boundary conditions */   
    for (p=1; p<=length; p++)
      {
	help[length-p]      = help[length+p-1];
	help[ny+length-1+p] = help[ny+length-p];
      }
    
    /* convolution step */
    for (j=length; j<=ny+length-1; j++)
        {
        /* calculate convolution */
        sum = conv[0] * help[j];
        for (p=1; p<=length; p++)
            sum = sum + conv[p] * (help[j+p] + help[j-p]);
        /* write back */
        u[i][j-length+1] = sum;
        }
    } /* for i */

/* disallocate storage for a row */
disalloc_vector (help, ny+length+length);

/* disallocate convolution vector */
disalloc_vector (conv, length+1);


return;

} /* gauss_conv */

/*----------------------------------------------------------------------------*/
void apply_nms
(
    float x1,       // in:  left neighbour
    float x2,       // in:  central value
    float *x3,      // out: new central value
    float x4        // in:  right neighbour
)

// sets a point to 0 when one of its neighbours (x1,x4) is larger

{
    if (x1 > x2 || x4 > x2)
        *x3 = 0.0f;
}


/*----------------------------------------------------------------------------*/
float get_direction
(
    float x,        // in: first component of vector
    float y         // in: second component of vector
)

// computes the edge direction in radian for a vector (x,y)^T
// used for grey value images

{
    float direction;
    if (fabsf(x) < EPSILON)
        if (fabsf(y) < EPSILON)
            direction = 0;
        else
            direction = 1.5708;
    else
        direction = atanf(y/x);

    return direction;
}


/*----------------------------------------------------------------------------*/
void nonmaxima_suppression
(
    float **u,      // in/out;  image
    float **dx,     // in:      x-derivatives
    float **dy,     // in:      y-derivatives
    int nx,         // in:      x-dimension of u
    int ny          // in:      y-dimension of u
)

// applies nonmaxima suppression

{
    int i,j;

    float **tmpu;
    alloc_matrix(&tmpu,nx+2,ny+2);

    // copy image
    for (i=0; i<nx+2; ++i)
        for (j=0; j<ny+2; ++j)
            tmpu[i][j] = u[i][j];

    for (i=1; i<nx+1; ++i)
        for (j=1; j<ny+1; ++j)
            if (tmpu[i][j] !=0.0f)
            {
                float direction;

                // get edge direction in pixel(i,j)
                direction = get_direction(dx[i][j],dy[i][j]);

                // apply nonmaxima suppression considering neighbours
                // in edge direction
                if (direction >= -0.3927 && direction <= 0.3927)
                    apply_nms(tmpu[i-1][j],tmpu[i][j],
                                                     &u[i][j],tmpu[i+1][j]);
                else if (direction > 0.3927 && direction <= 1.1781)
                    apply_nms(tmpu[i-1][j-1],tmpu[i][j],
                                                   &u[i][j],tmpu[i+1][j+1]);
                else if ((direction > 1.1781 && direction <= 1.5708) ||
                         (direction < -1.1781 && direction >= -1.5708))
                    apply_nms(tmpu[i][j-1],tmpu[i][j],
                                                     &u[i][j],tmpu[i][j+1]);
                else if (direction < -0.3927 && direction >= -1.1781)
                    apply_nms(tmpu[i+1][j-1],tmpu[i][j],
                                                   &u[i][j],tmpu[i-1][j+1]);
                else
                {
                    printf("Error in nonmaxima_suppression(): ");
                    printf("\nUndefined direction\n");
                    exit(-1);
                }
            }

    disalloc_matrix(tmpu,nx+2,ny+2);
}


/*----------------------------------------------------------------------------*/
void traceedge
(
    int i,          // in:      x component of current pixel
    int j,          // in:      y component of current pixel
    float **u,      // in/out:  image
    float T1,       // in:      lower threshold
    float T2        // in:      higher threshold
)

// Consider the neighbours of pixel (i,j). If a neighbour is under the higher 
// threshold T2 and over the lower threshold T1 add this pixel to the 
// edge pixels and repeat all steps for that pixel.

{
    u[i][j] = 255;

    if (u[i+1][j+1] <= T2 && u[i+1][j+1] > T1)
        traceedge(i+1,j+1,u,T1,T2);

    if (u[i+1][j] <= T2 && u[i+1][j] > T1)
        traceedge(i+1,j,u,T1,T2);

    if (u[i+1][j-1] <= T2 && u[i+1][j-1] > T1)
        traceedge(i+1,j-1,u,T1,T2);

    if (u[i-1][j+1] <= T2 && u[i-1][j+1] > T1)
        traceedge(i-1,j+1,u,T1,T2);

    if (u[i-1][j] <= T2 && u[i-1][j] > T1)
        traceedge(i-1,j,u,T1,T2);

    if (u[i-1][j-1] <= T2 && u[i-1][j-1] > T1)
        traceedge(i-1,j-1,u,T1,T2);

    if (u[i][j+1] <= T2 && u[i][j+1] > T1)
        traceedge(i,j+1,u,T1,T2);

    if (u[i][j-1] <= T2 && u[i][j-1] > T1)
        traceedge(i,j-1,u,T1,T2);
}


/*----------------------------------------------------------------------------*/
void getderivatives
(
    float **u,          // in:  image
    float **dx,         // out: derivatives in x-direction
    float **dy,         // out: derivatives in y-direction
    int nx,             // in:  x-dimension of u
    int ny              // in:  y-dimension of u
)

// compute the derivatives in x and y direction;
// the grid size h is assumed to be 1

{
    int i,j,p;
    // compute the derivatives of u using the sobel operator
    /*
      INSERT MISSING CODE HERE
    */
    float deri[3] = {-1.0/2, 0.0, 1.0/2}; // mirrored derivate kernel
    float bino[3] = {1.0/4, 2.0/4 , 1.0/4}; // (mirrored) binomial kernel
    int length = 3;
    float *help_x;
    float *help_y_1;
    float *help_y_2;
    alloc_vector (&help_x, nx+2);

    // derivative convolution in the x direction for dx
    // binomial convolution in the x direction for dy
    for (j=1; j<=ny; j++)
    {
        // copy into row vector
        for(p=0; p<=nx+1; p++)
            help_x[p] = u[p][j];

        // convolution along row (x-direction)
        for (p=1; p<=nx; p++)
        {
            float sum_dx = 0;
            float sum_by = 0;
            for(int k=0; k<length; k++)
            {
                sum_dx += help_x[p+k-1]*deri[k];
                sum_by += help_x[p+k-1]*bino[k];
            }
            dx[p][j] = sum_dx;
            dy[p][j] = sum_by;            
        }
    }
    disalloc_vector (help_x, nx+2);

    alloc_vector (&help_y_1, ny+2);
    alloc_vector (&help_y_2, ny+2);
    for (i=1; i<=nx; i++)
    {
        // copy column vector
        for (p=0; p<=ny+1; p++)
        {
            help_y_1[p] = dx[i][p];
            help_y_2[p] = dy[i][p];
        }

        // convolution along column (y-direction)
        for (p=1; p<=ny; p++)
        {
            float sum_bx = 0;
            float sum_dy = 0;            
            for(int k=0; k<length; k++)
            {
                sum_bx += help_y_1[p+k-1]*bino[k];
                sum_dy += help_y_2[p+k-1]*deri[k];
            }
            dx[i][p] = sum_bx;
            dy[i][p] = sum_dy;
        }
    }
    disalloc_vector (help_y_1, ny+2);
    disalloc_vector (help_y_2, ny+2);
}


/*----------------------------------------------------------------------------*/
void hysteresis_thresholding
(
    float **u,      // in/out:  image
    float T1,       // in:      lower threshold
    float T2,       // in:      higher threshold
    int nx,         // in:  x-dimension of u
    int ny          // in:  y-dimension of u
)

// applies hysteresis thresholding and creates the binary edge image

{
    int i,j;

    // use pixels that are larger than T2 as seed points to find
    // additional edge pixels that are larger than T1
    for (i=1; i<nx+1; ++i)
        for (j=1; j<ny+1; ++j)
        {
            if(u[i][j]>=T2) 
                traceedge(i,j,u,T1,T2);
        }

    // create a binary image
    for (i=1; i<nx+1; ++i)
        for (j=1; j<ny+1; ++j)
        {
            if(u[i][j]<255) 
                u[i][j]=0.0f;
        }
}
/*----------------------------------------------------------------------------*/



/*----------------------------------------------------------------------------*/
int main(int argc, char* argv[])
{
    // input parameters
    float sigma;    // Gaussian standard deviation
    float T1;       // lower threshold
    float T2;       // higher threshold

    char input[80];              // name of input image
    char output[80];             // name of edge image
    char row[80];                // for reading data
    unsigned char byte;          // for data conversion
    FILE   *inimage, *outimage;  // input file, output file


    int nx;             // image size in x direction
    int ny;             // image size in y direction

    float **f;       // original input image
    float **u;       // binary edge image

    float**  dx;        // derivative in x direction
    float**  dy;        // derivative in y direction


    int i,j;            // loop variables


    /*------------------------------------------------------------------------*/
    /* read input arguments, i.e. sigma, T1 and T2                            */

    if (argc == 6)
    {
       strcpy(input,argv[1]);
       strcpy(output,argv[2]);
       sigma = atof(argv[3]);
       T1 = atof(argv[4]);
       T2 = atof(argv[5]);
    }
    else
    {
        printf("\nSyntax:   %s <input.pgm/.ppm> <output.pgm> sigma T1 T2\n",argv[0]);
	    printf("Note: T2 should be larger than T1! \n\n");
        exit(0);
    }

    /*------------------------------------------------------------------------*/
    /* read original image                                                    */

    /* open pgm file and read header */
    inimage = fopen(input,"rb");
    fgets (row, 80, inimage);
    fgets (row, 80, inimage);

    while (row[0]=='#') fgets(row, 80, inimage);
    sscanf (row, "%i %i", &nx, &ny);
    fgets (row, 80, inimage);

    /* allocate storage */
    alloc_matrix(&f,nx+2,ny+2);

    /* read image data */
    for (j=1; j<=ny; j++)
        for (i=1; i<=nx; i++)
            f[i][j] = (float) getc (inimage);
    fclose(inimage);

    alloc_matrix(&u,nx+2,ny+2);

    // initialise
    for (i=0; i<nx+2; i++)
        for (j=0; j<ny+2; j++)
            u[i][j] = 0.0f;

    /*------------------------------------------------------------------------*/
    /* apply canny edge detection                                             */

    alloc_matrix(&dx,nx+2,ny+2);
    alloc_matrix(&dy,nx+2,ny+2);

    // smooth the image with a Gaussian of standard deviation sigma
    presmooth(f,u,nx,ny,sigma);

    // mirror_boundaries
    for (i=1; i<=nx; i++)
    {
        u[i][0]=u[i][1];
        u[i][ny+1]=u[i][ny];
    }
    for (j=0; j<=ny+1; j++)
    {
        u[0][j]=u[1][j];
        u[nx+1][j]=u[nx][j];
    }

    // get the image derivatives
    getderivatives(u,dx,dy,nx,ny);

    // set boundaries to zero
    for (i=1; i<=nx; i++)
    {
        u[i][0]=0.0f;
        u[i][ny+1]=0.0f;
    }
    for (j=0; j<=ny+1; j++)
    {
        u[0][j]=0.0f;
        u[nx+1][j]=0.0f;
    }

    // compute the gradient magnitude and threshold
    // with T1 to get edge candidates
    for (i=1; i<nx+1; ++i)
        for (j=1; j<ny+1; ++j)
        {
            u[i][j] = 0.0f;
            u[i][j] += dx[i][j]*dx[i][j] + dy[i][j]*dy[i][j];
            u[i][j] = sqrtf(u[i][j]);
            u[i][j] = (u[i][j] >= T1)? u[i][j]:0.0f;
        }

    // apply nonmaxima suppression
    nonmaxima_suppression(u,dx,dy,nx,ny);

    // apply hysteresis thresholding
    hysteresis_thresholding(u,T1,T2,nx,ny);

    /*------------------------------------------------------------------------*/
    /* write binary edge image and free memory                                */

    /* open file and write header (incl. filter parameters) */
    outimage = fopen (output, "wb");
    fprintf (outimage, "P5 \n");
    fprintf (outimage, "%i %i \n255\n", nx, ny);

    /* write image data and close file */
    for (j=1; j<ny+1; j++)
        for (i=1; i<nx+1; i++)
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
    printf("output image %s successfully written\n\n", output);

    disalloc_matrix(dx,nx+2,ny+2);
    disalloc_matrix(dy,nx+2,ny+2);

    disalloc_matrix(f,nx+2,ny+2);
    disalloc_matrix(u,nx+2,ny+2);

    return 0;
}
