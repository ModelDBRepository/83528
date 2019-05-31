#include <stdlib.h>
#include <stdio.h>


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

/***************************** free_i_4tensor *********************************/
/* free a int 1) 4tensor prevoiusly allocated by i_4tensor.                   */

void free_i_4tensor (int ****m, long time, long depth, long rows) {

    fprintf (stderr, "\nfree_i_4tensor: time=%ld depth=%ld rows=%ld\n", time, depth, rows);

    for (int k = 0; k < time; ++k) {
        for (int j = 0; j < depth; ++j) {
            for (int i = 0; i < rows; ++i) {
                fprintf (stderr, " free m[%d][%d][%d] ", k, j, i);
                free (m[k][j][i]);
                fprintf (stderr, " . ");
            }
        }
    }

    for (int k = 0; k < time; ++k) {
        for (int j = 0; j < depth; ++j) {
            fprintf (stderr, " free m[%d][%d] ", k, j);
            free (m[k][j]);
            fprintf (stderr, " . ");
        }
    }


    for (int k = 0; k < time; ++k) {
        fprintf (stderr, " free m[%d] ", k);
        free (m[k]);
        fprintf (stderr, " . ");
    }

    fprintf (stderr, " free m ");
    free(m);
    fprintf (stderr, " . ");
}


int main () {

 fprintf (stderr, "\nbla 1     ");
 int ****connected = i_4tensor(4, 1, 4, 1);
 fprintf (stderr, "\nbla 2     ");
 free_i_4tensor (connected, 4, 1, 4);
 fprintf (stderr, "\nbla 3     ");

 return 0;
}


/*
[cs0cwe@marvin kernel]$ g++ test.c
[cs0cwe@marvin kernel]$ ./a.out

bla 1
bla 2
free_i_4tensor: time=4 depth=1 rows=4
 free m[0][0][0]  .  free m[0][0][1]  .  free m[0][0][2]  .  free m[0][0][3]  .
 free m[1][0][0]  .  free m[1][0][1]  .  free m[1][0][2]  .  free m[1][0][3]  .
 free m[2][0][0]  .  free m[2][0][1]  .  free m[2][0][2]  .  free m[2][0][3]  .
 free m[3][0][0]  .  free m[3][0][1]  .  free m[3][0][2]  .  free m[3][0][3]  .
 free m[0][0]  .  free m[1][0]  .  free m[2][0]  .  free m[3][0]  .  free m[0] Segmentation fault

*/
