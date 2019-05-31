#ifndef _utils_H
#define _utils_H

void   floatprint (FILE *fp, double val);
DOUBLE *d_vector (int len);
int    *i_vector (long len);
float  *f_vector (int len);
int    **i_matrix (long rows, long cols);
DOUBLE **d_matrix (long rows, long cols);
float  **f_matrix (long rows, long cols);
unsigned char   **uc_matrix (long rows, long cols);
int    ***i_tensor (long depth, long rows, long cols);
int    ****i_4tensor (long time, long depth, long rows, long cols);
DOUBLE ***d_tensor (long depth, long rows, long cols);
void   free_i_matrix (int **m, long rows);
void   free_i_tensor (int ***m, long depth, long rows);
void   free_i_4tensor (int ****m, long time, long depth, long rows);
void   free_d_matrix (DOUBLE **m, long rows);
void   free_d_tensor (DOUBLE ***m, long depth, long rows);
int    init_matrix (DOUBLE **S, int rows, int cols, DOUBLE low, DOUBLE high);

int    init_matrix_zero (DOUBLE **S, int rows, int cols);
int    symmetrize_matrix (DOUBLE **S, int rows, int cols);
DOUBLE quad_length (DOUBLE *vec, int dim);
void   normalize (DOUBLE *vec, int dim, DOUBLE length_soll);
void   normalize_two (DOUBLE *vec1, int dim1, DOUBLE lsoll1,
                           DOUBLE *vec2, int dim2, DOUBLE lsoll2);
int    sub_mean_vector (DOUBLE *vec, int dim);
DOUBLE spherize_vector (DOUBLE *vec, int dim);
DOUBLE discretize (double val, double stepsize);
int    vec_init_homo (double *vec, int dim, double low,   double hig);
int    vec_init_gaus (double *vec, int dim, double sigma,double dumm);
int    vec_init_gibb (double *vec, int dim, double kurt, double dumm);
double dgauss (double a, double v, double s, double m);
double f_kurtosis (double *vec, int dim);

char   *get_pwd ();

#endif
