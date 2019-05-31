#ifndef _iter_H
#define _iter_H

#define DOUBLE float


typedef struct GLOBPTR {

  char **words;
  int  num_words;
  void *data;
  int  int_val;
  float float_val; /**used only once in total_spherize -- might be dropped again if no further use**/

} GLOBPTR;


typedef struct PARAMS {                              /**global parameters  g**/

#undef CONVERT_GLO
#undef CONVERT_SCA
#undef CONVERT_VEC
#undef CONVERT_ARR
#define CONVERT_GLO(a,b)  a b;
#define CONVERT_SCA(a,b)  /*nothing*/
#include "parameters.h"

 int     num_pointers_at_global;             /**good to know for printparams**/
 int     num_pointers;
 char    **ch_pointers;                                /**each can be a word**/
 GLOBPTR **pointers;       /**each which is used here will show to a GLOBPTR**/

} PARAMS;


typedef struct AREA {                                     /**all states  A[]**/

#undef CONVERT_GLO
#undef CONVERT_SCA
#define CONVERT_GLO(a,b)  /*nothing*/
#define CONVERT_SCA(a,b)  a b;
#include "parameters.h"

  int    d_n;       /**d_n = d_a * d_b**/
  int    *shuffle;  /**neuron list allows shuffled random permutation update**/

  DOUBLE **A;
  DOUBLE **B;
  DOUBLE **C;
  DOUBLE **D;
  DOUBLE **E;
  DOUBLE **F;
  DOUBLE **G;
  DOUBLE **H;
  DOUBLE **I;
  DOUBLE **J;
  DOUBLE **K;
  DOUBLE **L;
  DOUBLE **M;
  DOUBLE **N;
  DOUBLE **O;
  DOUBLE **P;
  DOUBLE **Q;
  DOUBLE **R;
  DOUBLE **S;
  DOUBLE **T;
  DOUBLE **U;
  DOUBLE **V;
  DOUBLE **W;
  DOUBLE **X;
  DOUBLE **Y;
  DOUBLE **Z;

} AREA;



/**obsolete -- use pointers:**
typedef struct STATE {                                     **observables  o**
  double *kurtosis;
  double *variance;
  int    bins;
  double range;
  double **distr;
} STATE;
**/


#define REPORT_ERR 0
#define ERR(x)  if  (REPORT_ERR)  fprintf (stderr, "%s ", x);

#endif
