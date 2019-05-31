#include <stdlib.h>
#include <stdio.h>
#include <string.h>  /**for get_tmp_uname_filename**/
#include <math.h>

#include "coco.h"    /**for DOUBLE**/

/***************************** floatprint ************************************/
/* Prints to 1) fp  value of 2) val  only to precision of 1 or 2 to omit 0's */
void floatprint (FILE *fp, double val) {

    if  (fabs (val * 10.0 - (double)((int)(val * 10.0))) < 0.0000001)
        fprintf (fp, "%.1f", val);
    else {
        if  (fabs (val * 100.0 - (double)((int)(val * 100.0))) < 0.0000001)
            fprintf (fp, "%.2f", val);
        else {
            if  (fabs (val * 1000.0 - (double)((int)(val * 1000.0))) < 0.0000001)
                fprintf (fp, "%.3f", val);
            else {
                fprintf (fp, "%g", val);
            }
        }
    }
}

/******************************* d_vector ************************************/
/* Allocate a DOUBLE vector with subscript range v[0..len].                  */
DOUBLE *d_vector (int len) {
    DOUBLE *v;

    if  ((v = (DOUBLE *)malloc (len * sizeof (DOUBLE))) == NULL)
        fprintf (stderr, "\nallocation failure in d_vector()\n");

    return (v);
}

/***************************** f_vector **************************************/
/* Allocate a float vector with subscript range v[0..len].                    */
float *f_vector (int len) {
    float *v;

    if  ((v = (float *)malloc (len * sizeof (float))) == NULL)
        fprintf (stderr, "\nallocation failure in f_vector()\n");

    return (v);
}

/******************************* d_matrix ************************************/
/* Allocate a DOUBLE matrix with subscript range m[0..rows][0..cols].        */
DOUBLE **d_matrix (long rows, long cols) {
    DOUBLE **m;

    size_t r = rows;
    size_t c = cols;

//fprintf (stderr, "\nbegin d_matrix; r = %d, c = %d\n", r, c);

    if  ((m = (DOUBLE **)malloc (r * sizeof (DOUBLE *))) == NULL)
        fprintf (stderr, "\nallocation failure 1 in d_matrix()\n");

    for (unsigned int i = 0; i < r; i++) {
        if  ((m[i] = (DOUBLE *)calloc (c, sizeof (DOUBLE))) == NULL) /**<-- malloc replaced by calloc to init activations with zero**/
            fprintf (stderr, "\nallocation failure 2 in d_matrix()\n");
    }
//fprintf (stderr, "...end d_matrix\n");

    return (m);
}

/******************************* f_matrix ************************************/
/* Allocate a float matrix with subscript range m[0..rows][0..cols].        */
float **f_matrix (long rows, long cols) {
    int i;
    float **m;
    int err = 0;

    if  ((m = (float **)malloc (rows * sizeof (float *))) == NULL)
        err = 1;

    for (i = 0; i < rows; ++i)
        if  ((m[i] = (float *)malloc (cols * sizeof (float))) == NULL)
            err = 2;

    if  (err == 1)
        fprintf (stderr, "\nallocation failure 1 in f_matrix()\n");
    if  (err == 2)
        fprintf (stderr, "\nallocation failure 2 in f_matrix()\n");
    if  (err)
        fprintf (stderr, "cannot allocate %ld*%ld f_elements\n", rows, cols);

    return (m);
}

/******************************* uc_matrix ***********************************/
/* Allocate an unsigned char matrix with subscript range m[0..rows][0..cols].*/
unsigned char **uc_matrix (long rows, long cols) {
    int i;
    unsigned char **m;
    int err = 0;

    if  ((m = (unsigned char **)malloc (rows * sizeof (unsigned char *)))
         == NULL)
        err = 1;

    for (i = 0; i < rows; ++i)
        if  ((m[i] = (unsigned char *)malloc (cols * sizeof (unsigned char)))
             == NULL)
            err = 2;

    if  (err == 1)
        fprintf (stderr, "\nallocation failure 1 in uc_matrix()\n");
    if  (err == 2)
        fprintf (stderr, "\nallocation failure 2 in uc_matrix()\n");
    if  (err)
        fprintf (stderr, "cannot allocate %ld*%ld elements\n", rows, cols);

    return (m);
}

/***************************** i_vector **************************************/
/* Allocate an int vector with subscript range v[0..len].                    */

int *i_vector (long len) {
    int *v;

    if  ((v = (int *)malloc (len * sizeof (int))) == NULL)
        fprintf (stderr, "\nallocation failure in i_vector()\n");

    return (v);
}

/******************************* i_matrix ************************************/
/* Allocate an int matrix with subscript range m[0..rows][0..cols].          */

int **i_matrix (long rows, long cols) {
    int i;
    int **m;
    int err = 0;

    if  ((m = (int **)malloc (rows * sizeof (int *))) == NULL)
        err = 1;

    for (i = 0; i < rows; ++i)
        if  ((m[i] = i_vector(cols)) == NULL)
            err = 2;

    if  (err == 1)
        fprintf (stderr, "\nallocation failure 1 in i_matrix()\n");
    if  (err == 2)
        fprintf (stderr, "\nallocation failure 2 in i_matrix()\n");
    if  (err)
        fprintf (stderr, "cannot allocate %ld*%ld elements\n", rows, cols);

    return (m);
}

/******************************** i_tensor ***********************************/
/* Allocate int tensor, subscript ranges m[0..depth][0..rows][0..cols].      */

int ***i_tensor (long depth, long rows, long cols) {
    int k;
    int ***t;

    if  ((t = (int ***)malloc (depth * sizeof (int **))) == NULL)
        fprintf (stderr, "\nallocation failure 1 in i_tensor()\n");

    for (k = 0; k < depth; ++k)
        if  ((t[k] = i_matrix (rows, cols)) == NULL)
            fprintf (stderr, "\nallocation failure 2 in i_tensor()\n");

    return (t);
}

/******************************** i_4tensor **********************************/
/* Allocate int 4D-tensor, ranges m[0..time][0..depth][0..rows][0..cols].    */

int ****i_4tensor (long time, long depth, long rows, long cols) {
    int ****t;

    if  ((t = (int ****)malloc (time * sizeof (int ***))) == NULL)
        fprintf (stderr, "\nallocation failure 1 in i_4tensor()\n");

    for (int k = 0; k < time; ++k)
        if  ((t[k] = i_tensor (depth, rows, cols)) == NULL)
            fprintf (stderr, "\nallocation failure 2 in i_4tensor()\n");

    return (t);
}

/***************************** free_i_matrix *********************************/
void free_i_matrix (int **m, long rows) {

    for (int i = 0; i < rows; ++i) {
        free (m[i]);
    }
    free (m);
}

/***************************** free_i_tensor *********************************/
void free_i_tensor (int ***m, long depth, long rows) {

    for (int i = 0; i < depth; ++i)
        free_i_matrix (m[i], rows);
    free (m);
}

/***************************** free_i_4tensor *********************************/
void free_i_4tensor (int ****m, long time, long depth, long rows) {

    for (int k = 0; k < time; ++k)
        free_i_tensor (m[k], depth, rows);
    free (m);
}



/******************************** d_tensor ***********************************/
/* Allocate DOUBLE tensor, subscript ranges m[0..depth][0..rows][0..cols].   */
DOUBLE ***d_tensor (long depth, long rows, long cols) {
    int k;
    DOUBLE ***t;

//    fprintf (stderr, "\nd_tensor: depth = %ld, rows = %ld, cols = %ld\n", depth, rows, cols);

    if  ((t = (DOUBLE ***)malloc (depth * sizeof (DOUBLE **))) == NULL)
        fprintf (stderr, "\nallocation failure in tensor()\n");

    for (k = 0; k < depth; ++k)
        t[k] = d_matrix (rows, cols);

    return (t);
}


/***************************** free_d_matrix *********************************/
/* free a DOUBLE 1) matrix with 2) rows prevoiusly allocated by matrix ().   */

void free_d_matrix (DOUBLE **m, long rows) {

    for (int i = 0; i < rows; ++i)
        free (m[i]);
    free (m);
}


/***************************** free_d_tensor *********************************/
/* free a DOUBLE 1) tensor prevoiusly allocated by d_tensor.                 */

void free_d_tensor (DOUBLE ***m, long depth, long rows) {

    for (int i = 0; i < depth; ++i) {
        for (int j = 0; j < rows; ++j)
            free (m[i][j]);
        free (m[i]);
    }
    free (m);
}


/***************************** init_matrix ***********************************/
/* Uses drand48 () !                                                         */
int init_matrix (DOUBLE **S, int rows, int cols, DOUBLE low, DOUBLE high) {
    int i, j;

    for (i = 0; i < rows; ++i)
        for (j = 0; j < cols; ++j)
            S[i][j] = low + drand48 () * (high - low);

    return (0);
}


/***************************** init_matrix_zero ******************************/
int init_matrix_zero (DOUBLE **S, int rows, int cols) {
    int i, j;

    for (i = 0; i < rows; ++i)
        for (j = 0; j < cols; ++j)
            S[i][j] = 0.0;

    return (0);
}


/*************************** symmetrize_matrix *******************************/
int symmetrize_matrix (DOUBLE **S, int rows, int cols) {
    int i, j;

    if  (rows != cols)
        fprintf (stderr, "\nonly square matrices to by symmetrized\n");

    for (i = 0; i < rows; ++i)                 /**average values in one half**/
        for (j = 0; j < i; ++j)
            S[i][j] = 0.5 * (S[i][j] + S[j][i]);

    for (i = 0; i < rows; ++i)                  /**same values in other half**/
        for (j = rows - 1; j > i; --j)
            S[i][j] = S[j][i];

    return (0);
}



/*************************** quad_length *************************************/
DOUBLE quad_length (DOUBLE *vec, int dim) {
  int i;
  DOUBLE out = 0.0;

  for (i = 0; i < dim; ++i)
      out += vec[i] * vec[i];

  return (out);
}


/**************************** normalize **************************************/
/* Normalizes 1) vector of size 2) dim to 3) length_soll. Uses quad_length.  */

void normalize (DOUBLE *vec, int dim, DOUBLE length_soll) {
    DOUBLE length_ist;
    int k;

    static int firsttime = 1;

    if  (length_soll == 0.0) {
        for (k = 0; k < dim; ++k)
            vec[k] = 0.0;
    } else {
        length_ist = sqrt (quad_length (vec, dim));
        if  (length_ist > 0.00001) {
            for (k = 0; k < dim; ++k)
                vec[k] *= length_soll / length_ist;
        } else {
            if  (firsttime == 1) {
                fprintf (stderr, "--no normalize possible--");
                fprintf (stderr, " next printing only \"-\" to stdout--\n");
                firsttime = 0;
            } else {
                 printf ("-");
            }
        }
    }
}


void normalize_two (DOUBLE *vec1, int dim1, DOUBLE lsoll1,
                           DOUBLE *vec2, int dim2, DOUBLE lsoll2) {
  DOUBLE lenof1 = 0.0, lenof2 = 0.0, lenall;
  int k;

  if  (lsoll1 != 0.0)
      lenof1 = quad_length (vec1, dim1);
  if  (lsoll2 != 0.0)
      lenof2 = quad_length (vec2, dim2);
  lenall = lenof1 + lenof2;
  lenof1 = sqrt (lenof1);
  lenof2 = sqrt (lenof2);
  lenall = sqrt (lenall);

  if  (lenall == 0.0) {
      fprintf (stderr, "  length of W-M too small   ");
      return;
  }

  if  (lsoll1 == lsoll2) {
      for (k = 0; k < dim1; ++k)
          vec1[k] *= lsoll1 / lenall;
      for (k = 0; k < dim2; ++k)
          vec2[k] *= lsoll2 / lenall;
  } else {
      if  (lsoll1 == 0.0) {
          for (k = 0; k < dim2; ++k)
              vec2[k] *= lsoll2 / lenall;
          for (k = 0; k < dim1; ++k)
              vec1[k] = 0.0;

      } else {
          if  (lsoll2 == 0.0) {
              for (k = 0; k < dim1; ++k)
                  vec1[k] *= lsoll1 / lenall;
              for (k = 0; k < dim2; ++k)
                  vec2[k] = 0.0;
          } else {
              fprintf (stderr, "\nfunny parameters if W-M normalized! ");
          }
      }
  }
}

/*************************** sub_mean_vector *********************************/
int sub_mean_vector (DOUBLE *vec, int dim) {
    int i;
    DOUBLE sp;

    sp = 0.0;                            /**Mittelwert abziehen**/
    for (i = 0; i < dim; ++i)
        sp += vec[i];
    sp /= (DOUBLE)(dim);
    for (i = 0; i < dim; ++i)
        vec[i] -= sp;

    return (0);
}


/*************************** spherize_vector *********************************/
/* Sets variance of componetents to one. Always use sub_mean_of_vector first!*/
/* Returns value by which vector elements are divided.                       */

DOUBLE spherize_vector (DOUBLE *vec, int dim) {
    int i;
    DOUBLE sig = 0.0;

    static int firsttime = 1;

    if  (quad_length (vec, dim) < 0.00001) {
        if  (firsttime == 1) {
            fprintf (stderr, "  --small vector in spherize_vector -");
            fprintf (stderr, " next printing only \"-\" to stdout--\n");
            firsttime = 0;
        } else {
            printf ("-");
        }
        return sig;
    }

    /**cheap version of spherize**/
    for (i = 0; i < dim; ++i)
        sig += vec[i] * vec[i];
    sig /= (DOUBLE)(dim);
    sig = sqrt (sig);

    if  (sig > 0.00001)
        for (i = 0; i < dim; ++i)           /**pixels (coordinates)**/
            vec[i] /= sig;                  /**will have 1 variance**/
    else
        fprintf (stderr, "--flat vector--");

    return sig;
}



/**************************** discretize *************************************/
/* 1) val becomes next smaller value (as abslolute) in units of 1/steps 2).  */
DOUBLE discretize (DOUBLE val, DOUBLE stepsize) {

 if  (stepsize == 0.0)
     return (val);
 else                          /**stepsize = discretization value:**/
     return ((DOUBLE)(floor(val / stepsize - 0.5) + 1.0) * stepsize);

/**  steps = number of steps within an interval of 1:
     return ((DOUBLE)(floor(val * steps - 0.5) + 1.0) / steps);
**/
}



/************************** vec_init_xxxx ************************************/
/* Functions used to initialize a vector 1) vec of length 2) dim.            */
/* double arguments 3), 4) may be dummys.                                    */
/* Return values not interesting.                                            */

/* Uses drand48 () !                                                         */
int vec_init_homo (double *vec, int dim, double low,   double hig) {
    int i;

    for (i = 0; i < dim; ++i)
        vec[i] = low + drand48 () * (hig - low);

    return (0);
}

int vec_init_gaus (double *vec, int dim, double sigma,double dumm) {
    fprintf (stderr, "\nvec_init_gaus not implemented\n");

    return (0);
}

/* vec gets values -1, 0, +1 such that distribution has kurtosis kurt        */
int vec_init_gibb (double *vec, int dim, double kurt, double dumm) {
    int i;

    if  (dumm > 3.0) {  /**only positive values for init**/
        for (i = 0; i < dim; ++i)
            vec[i] = ((drand48 () * (kurt + 1.0)) > 1.0) ? 0.0 : 1.0;
    } else {
        for (i = 0; i < dim; ++i)
            vec[i] = ((drand48 () * (kurt + 2.0)) > 2.0)
                   ? 0.0
                   : (drand48 () > 0.5 ? 1.0 : -1.0);
    }

    return (0);
}


/* Differential of Gaussian; x-value a, integral v, sigma s, mu m.           */
double dgauss (double a, double v, double s, double m) {

    return (-v * sqrt(2.0) * (a-m) * exp (- (a-m)*(a-m) / (s*s))
                                   / (2.0 * sqrt(M_PI) * s*s));
}


/**************************** f_kurtosis *************************************/
/* Returns non-normalized Fisher kurtosis on elements of 1) vec  with 2) dim.*/
/* Assumes that mean is zero.                                                */

double f_kurtosis (double *vec, int dim) {
    int i;
    double merk, zaehler = 0.0, nenner = 0.0;

    if  (dim == 0)  fprintf (stderr, "\nf_kurtosis: dim = 0 ! ");

    for (i = 0; i < dim; ++i) {
        merk = vec[i] * vec[i];
        zaehler += merk * merk;
        nenner  += merk;
    }

    zaehler /= (double)(dim);
    nenner  /= (double)(dim);
    nenner  *= nenner;

    if  (nenner == 0.0) {
        /* fprintf (stderr, ** "\nf_kurtosis: nenner = 0.0 ! " ** "n"); !!*/

        return (0.0);
    }

    return (zaehler / nenner - 3.0);
}



/************************* get_pwd *******************************************/
/* Alloc mem and return as string the current directory using system(pwd).   */

char *get_pwd () {
  FILE *fp;
  char myname[256];
  char *returnstring;

  system ("pwd > /tmp/l4vPz3Nrb7b");
  fp = fopen ("/tmp/l4vPz3Nrb7b", "r");
  if  (1 != fscanf (fp, "%s", myname)) {
      fprintf (stderr, "\nsth wrong reading %s from /tmp/l4vPz3Nrb7b", myname);
      exit (1);
  }
  fclose (fp);
  system ("rm /tmp/l4vPz3Nrb7b");

  if  (NULL == (returnstring = (char *)malloc ((strlen (myname) + 1) * sizeof (char)))) {
      fprintf (stderr, "\nno malloc possible in get_pwd");
      exit (1);
  }
  strcpy (returnstring, myname);

  fprintf (stderr, "\nget_pwd returns %s  ", returnstring);

  return (returnstring);
}
