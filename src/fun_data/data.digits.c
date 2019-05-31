#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "iter.h"
#include "series.h"
#include "utils.h"
#include "data.h"


/******************************** torus_shift ********************************/
/* Writes to 1) imout a randomly shifted version of 2) im_in. Periodic bound.*/

void torus_shift (double *image, int d_a, int d_b) {
     int i, j;
     int shi = (int)(drand48() * d_a);
     int shj = (int)(drand48() * d_b);

     double tmp[64][64];

     if  ((d_a > 64) || (d_b > 64))
         fprintf (stderr, "\n\nplease enlarge torus_shift for this area!\n\n");

     for (i = 0; i < d_a; ++i)
         for (j = 0; j < d_b; ++j)
             tmp[i][j] = image[((i+shi)%d_a)*d_b + (j+shj)%d_b];

     for (i = 0; i < d_a; ++i)
         for (j = 0; j < d_b; ++j)
             image[i*d_b + j] = tmp[i][j];
}

/******************************** init_digit *********************************/
/* cmd->S_target is initialized with a digit of size z->d_a * z->d_b.        */
/* Usage like init_image but additional teacher value.                       */
/* Uses sub_mean_vector for each image part (but NOT spherize_vector).       */
/* cmd->offset == 1/-1 -> only then get new image.                           */
/* quantum[0][0] == 1: get new image.                                        */
/* quantum[1][0] scales the image pixels (for bounded transfkt of backprop). */
/* quantum[2][0] == 1: sub_mean_vector of each datapoint.                    */
/* quantum[3][0] == 1: torus_shift.                                          */
/* quantum[4][0] == 1: teacher.                                              */
/* Suited for resolutions: 16x15, 8x8.                                       */

void init_digit (AGENT *z, COMMAND *cmd, DATA *d, int begin, int end) {

    static int firsttime = 1;
    static double **Alldigits;
    static int select;
    double *image;
    int ct_t, i, j, d_r;
    int teacher = (int)cmd->quantum[4][0];

    if  (teacher) {
        if  (z->d_b != 10)
            fprintf (stderr, "\n\ninit_digit teacher wants input size 10\n\n");
    } else {
        if  (! (((z->d_a == 16) && (z->d_b == 15))
             || ((z->d_a == 8)  && (z->d_b == 8 ))))
            fprintf (stderr, "\n\ninit_digit input size: 15x16 or 8x8\n\n");
    }
    if  (cmd->anz_quantums != 5)
        fprintf (stderr, "\n5 parameters for init_digit, please!\n");
    if  (cmd->quantum[1][0] == 0.0)
        fprintf (stderr, "\n\ninit_digit wants quantum[1][0] set!\n\n");

    if  (firsttime) {
        FILE *fp;
        char dummy;

        Alldigits = d_matrix (2000, 16 * 15);

        if  ((fp=fopen ("/home/ni/cweber/p/co99/d/mfeat/mfeat-pix","r")) == 0)
            fprintf (stderr, "\n\ninit_digit cannot find file mfeat-pix\n\n");

        for (i = 0; i < 2000; ++i)
            for (j = 0; j < 16 * 15; ++j) {
                fscanf (fp, "%lf", &Alldigits[i][j]);
                Alldigits[i][j] /= 6.0;           /**because data range 0..6**/
            }

        fscanf (fp, "%c", &dummy);  /**read the last return**/
        fscanf (fp, "%c", &dummy);  /**read the eof char**/
        if  (!feof (fp))
            fprintf (stderr, "\n\ninit_digit: no correct file end found!\n\n");

        fclose (fp);
        fprintf (stderr, "  -- read the multiple-feature digit file\n");
        firsttime = 0;
    }

    image = cmd->S_target[0];

    /**get new digit**/
    if  (cmd->quantum[0][0] == 1)
        select = (int)(drand48() * 2000);

    d_r = z->d_a * z->d_b;

    if  (! teacher) {
        if  (z->d_a == 16)       /**for size 16x15**/
            for (j = 0; j < 16 * 15; ++j)
                image[j] = Alldigits[select][j];
        if  (z->d_a == 8) {      /**for size 8x8**/
            int i;
            for (i = 0; i < 8; ++i)
            for (j = 0; j < 8; ++j)
                if  (j*2+1 != 15)
                    image[i*8+j] = ( Alldigits[select][2*i*15+j*2]
                                   + Alldigits[select][2*i*15+j*2+1]
                                   + Alldigits[select][(2*i+1)*15+j*2]
                                   + Alldigits[select][(2*i+1)*15+j*2+1])*0.25;
                else
                    image[i*8+j] = ( Alldigits[select][2*i*15+j*2]
                                   + Alldigits[select][(2*i+1)*15+j*2])*0.375;
        }
    } else {
        for (j = 0; j < z->d_b; ++j)
            for (i = 0; i < z->d_a; ++i)  /**fill redundantly along ho_a**/
                if  (j == select / 200)
                    image[z->d_b * i + j] = cmd->quantum[1][0]; /**scale**/
                else
                    image[z->d_b * i + j] = 0.0;
    }

                                              /**from init_image. there:**/
    if  (cmd->quantum[2][0])                  /**must NOT be set for ICA**/
        sub_mean_vector (image, d_r);         /**but should for BM / old**/

    /* spherize_vector (image, d_r); */

    if  ((int)cmd->quantum[3][0])
        torus_shift (image, z->d_a, z->d_b);

    /**copy digit to target at all times**/
    for (ct_t = begin; ct_t < end; ++ct_t)
        for (j = 0; j < d_r; ++j)
            cmd->S_target[ct_t][j] = image[j] * cmd->quantum[1][0]; /**scale**/
}




/******************************** init_mnist *********************************/
/* Like init_digit but reads mnist data set.                                 */
/* Suited for resolutions: 28x28, 14x14, 7x7; not yet 8x8 because 7x7 bad.   */

void init_mnist (AGENT *z, COMMAND *cmd, DATA *d, int begin, int end) {

    static int firsttime = 1;
    static unsigned char **Alldigits;
    static int *allteach;
    static int select;
    double *image;
    int ct_t, i, j, d_r;
    const int dummy = (cmd->anz_quantums == 5) ? 0 :
                  fprintf (stderr, "\n5 parameters for init_mnist, please!\n");
    int teacher = (int)cmd->quantum[4][0];

    if  (teacher) {
        if  (z->d_b != 10)
            fprintf (stderr, "\n\ninit_digit teacher wants input size 10\n\n");
    } else {
        if  (! (((z->d_a == 28) && (z->d_b == 28))
             || ((z->d_a == 14) && (z->d_b == 14))
             || ((z->d_a == 7)  && (z->d_b == 7 ))))
            fprintf (stderr,
                   "\n\ninit_mnist input size: 28x28 / 14x14 / 7x7 / 8x8\n\n");
    }
    if  (cmd->quantum[1][0] == 0.0)
        fprintf (stderr, "\ninit_mnist wants quantum[1][0] set!%d\n\n", dummy);

    if  (firsttime) {
        FILE *fp;

        /**read the images**/
        Alldigits = uc_matrix (5000, 28 * 28);   /**only first & better half**/

        if  ((fp=fopen("d/mnist/t10k-images-idx3-ubyte", "r")) == 0)
            fprintf (stderr, "\n\ninit_mnist cannot find mnist images\n\n");

        /**read the header**/
        fread (&i, sizeof(int), 1, fp);
        fprintf (stderr, "\nmagicnumber=%d ", i);
        fread (&i, sizeof(int), 1, fp);
        fprintf (stderr, " images:%d ", i);
        fread (&i, sizeof(int), 1, fp);
        fprintf (stderr, " rows:%d ", i);
        fread (&i, sizeof(int), 1, fp);
        fprintf (stderr, " columns:%d  ", i);

        for (i = 0; i < 5000; ++i)
            for (j = 0; j < 28 * 28; ++j)
                fscanf (fp, "%c", Alldigits[i]+j); /**data range 0..256 !**/

        fclose (fp);
        fprintf (stderr, " -- read the mnist digits");

        /**read the labels**/
        allteach = i_vector (5000);

        if  ((fp=fopen("d/mnist/t10k-labels-idx1-ubyte",
                       "r")) == 0)
            fprintf (stderr, "\n\ninit_mnist cannot find mnist labels\n\n");

        /**read the header**/
        fread (&i, sizeof(int), 1, fp);
        fprintf (stderr, "\nmagicnumber=%d ", i);
        fread (&i, sizeof(int), 1, fp);
        fprintf (stderr, " items:%d ", i);

        for (i = 0; i < 5000; ++i) {
            char c;
            fscanf (fp, "%c", &c);
            allteach[i] = (int)c;
        }

        fclose (fp);
        fprintf (stderr, " -- read the mnist labels\n");

        firsttime = 0;
    }

    image = cmd->S_target[0];

    /**get new digit**/
    if  (cmd->quantum[0][0] == 1)
        select = (int)(drand48() * 5000);

    d_r = z->d_a * z->d_b;

    if  (! teacher) {
        if  (z->d_a == 28)        /**for size 28x28**/
            for (j = 0; j < 28 * 28; ++j)
                image[j] = (double)Alldigits[select][j] / 256.0;
        if  (z->d_a == 14) {      /**for size 14x14**/
            int i;
            double fac = 1.0 / 4.0 / 256.0;
            for (i = 0; i < 14; ++i)
            for (j = 0; j < 14; ++j)
                image[i*14+j] = ((double)Alldigits[select][2*i*28+j*2]
                                + Alldigits[select][2*i*28+j*2+1]
                                + Alldigits[select][(2*i+1)*28+j*2]
                                + Alldigits[select][(2*i+1)*28+j*2+1])*fac;
        }
        if  (z->d_a == 7) {      /**for size 7x7 <-- digits don't look good!**/
            int i;
            double fac = 1.0 / 16.0 / 256.0;
            for (i = 0; i < 7; ++i)
            for (j = 0; j < 7; ++j)
                image[i*7+j] = ((double)Alldigits[select][4*i*28+j*4]
                                + Alldigits[select][4*i*28+j*4+1]
                                + Alldigits[select][4*i*28+j*4+2]
                                + Alldigits[select][4*i*28+j*4+3]
                                + Alldigits[select][(4*i+1)*28+j*4]
                                + Alldigits[select][(4*i+1)*28+j*4+1]
                                + Alldigits[select][(4*i+1)*28+j*4+2]
                                + Alldigits[select][(4*i+1)*28+j*4+3]
                                + Alldigits[select][(4*i+2)*28+j*4]
                                + Alldigits[select][(4*i+2)*28+j*4+1]
                                + Alldigits[select][(4*i+2)*28+j*4+2]
                                + Alldigits[select][(4*i+2)*28+j*4+3]
                                + Alldigits[select][(4*i+3)*28+j*4]
                                + Alldigits[select][(4*i+3)*28+j*4+1]
                                + Alldigits[select][(4*i+3)*28+j*4+2]
                                + Alldigits[select][(4*i+3)*28+j*4+3])*fac;
        }
    } else {
        for (j = 0; j < z->d_b; ++j)
            for (i = 0; i < z->d_a; ++i)  /**fill redundantly along ho_a**/
                if  (j == allteach[select])
                    image[z->d_b * i + j] = cmd->quantum[1][0]; /**scale**/
                else
                    image[z->d_b * i + j] = 0.0;
    }

                                              /**from init_image. there:**/
    if  (cmd->quantum[2][0])                  /**must NOT be set for ICA**/
        sub_mean_vector (image, d_r);         /**but should for BM / old**/

    /* spherize_vector (image, d_r); */

    if  ((int)cmd->quantum[3][0])
        torus_shift (image, z->d_a, z->d_b);

    /**copy digit to target at all times**/
    for (ct_t = begin; ct_t < end; ++ct_t)
        for (j = 0; j < d_r; ++j)
            cmd->S_target[ct_t][j] = image[j]
                                   * cmd->quantum[1][0]; /**scale**/
}
