#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "iter.h"
#include "series.h"
#include "utils.h"      /**for init_matrix_zero, d_tensor**/
#include "vehicle.h"    /**for print_command to debug**/
#include "local.h"      /**for local_mean_01 in weight_rec_mean_01**/

/************** functions pointed to by cmd->weightfunc **********************/

/*************************** weight_full_hebb ********************************/
/* Updates dW_target of neuron 5) ct_n  from cmd->area==n_from1[0].          */
/* Pre: cmd->S_from2; Post (here): cmd->S_from1[0]. (NICHT kommutativ!)      */
/* Stepsize weighted by cmd->moment.                                         */

void weight_full_hebb (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int inarea, ct_l, ct_in;
    double *dW, *Pre;

    const int area = cmd->area;

    const double post = cmd->S_from1[0][ct_t][ct_n];    /**assumed as "here"**/

    if  ((cmd->anz_from1 != 1) || (cmd->n_from1[0] != cmd->area))
        fprintf (stderr, "wrong use of weight_hebb");

    for (ct_l = 0; ct_l < cmd->anz_from2; ++ct_l) {

        inarea = cmd->n_from2[ct_l];

        dW = cmd->dW_target[inarea][ct_n];
        /* W_delta = A[area].W_delta[inarea][ct_n]; */

        Pre = cmd->S_from2[ct_l][ct_t];

        for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)

            dW[ct_in] += cmd->moment * Pre[ct_in] * post;

        if  (z->ch_Theta)

            A[area].Theta_delta[ct_n] -= cmd->moment * post;
    }
}



/*************************** weight_full_kohonen *****************************/
/* Like weight_hebb, but Pre[] replaced by (Pre[]-W[]).                      */
/* Expects neighborhood function as post.                                    */

void weight_full_kohonen (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int inarea, ct_l, ct_in;
    double *dW, *W, *Pre;

    const int area = cmd->area;

    const double post = cmd->S_from1[0][ct_t][ct_n];    /**assumed as "here"**/

    if  ((cmd->anz_from1 != 1) || (cmd->n_from1[0] != area))
        fprintf (stderr, "wrong use of weight_kohonen");

    for (ct_l = 0; ct_l < cmd->anz_from2; ++ct_l) {

        inarea = cmd->n_from2[ct_l];

        dW = cmd->dW_target[inarea][ct_n];
        W  = cmd->W_target[inarea][ct_n];

        Pre = cmd->S_from2[ct_l][ct_t];

        for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in) {

            dW[ct_in] += cmd->moment * (Pre[ct_in] - W[ct_in]) * post;

            if  (W[ct_in] > 256.0) {
                fprintf (stderr, "\nweight_kohonen exits because weights are presumably too large .. too large data?\n");
                exit (1);
            }
        }

        if  (z->ch_Theta)

            fprintf (stderr, "no thresholds for weight_kohonen, please");
    }
}





/****************************** weight_full_self_zero ************************/
/* Sets diagonal of self-connecting weight matrix (autapses) to zero.        */

void weight_full_self_zero (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    if  ((cmd->anz_from1 != 1) || (cmd->anz_from2 != 0))
        fprintf (stderr, "\nweight_self_zero: too many areas\n");
    if  (cmd->n_from1[0] != cmd->area)
        fprintf (stderr, "\nweight_self_zero: only self-connections!\n");

    cmd->W_target[cmd->area][ct_n][ct_n] = 0.0;
}






/****************************** weight_init **********************************/
/* Allocates memory for W, dW and initializes. No need for scaff_W anymore.  */
/* To be written later, if possible. Problems:                               */
/* - a->W[in_area] etc must be allocated (have addresses)                    */
/*   BEFORE cmd->W_target is assigned to A[cmd->area].W .                    */
/* - areas is unknown --> hand over from x to A                           */
/*
void weight_init (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

in weight.c:
    W  = cmd->W_target[inarea][ct_n];

in vehicle.c  alloc_a:
    for (ct_ar = 0; ct_ar < areas; ++ct_ar) {
        a = A + ct_ar;
        z = Z + ct_ar;
        for (in_ar = 0; in_ar < x->areas; ++in_ar) {
            **weights to areas spezified in z->scaff_W/_V**
            if  (z->scaff_W[in_ar] != 0)
                a->W[in_ar]  = d_matrix (a->d_r, A[in_ar].d_r);
            else
                a->W[in_ar]  = NULL;

in vehicle.c  choose_pointers:
       if  (cmd->ch_target == 'W')
           cmd->W_target  = A[cmd->area].W;
******************************************************************************/


    A = (AREA *)malloc (x->areas * sizeof (AREA));

    /**essential within-one area values (must be ready for next loop)**/
    for (ct_ar = 0; ct_ar < x->areas; ++ct_ar) {

        /**weight pointers**/
        a->W  = (double ***)malloc (x->areas * sizeof (double **));
        a->dW = (double ***)malloc (x->areas * sizeof (double **));
    }

    /**between-area values**/
    for (ct_ar = 0; ct_ar < x->areas; ++ct_ar) {

        a = A + ct_ar;
        z = Z + ct_ar;

        for (in_ar = 0; in_ar < x->areas; ++in_ar) {

            /**weights to areas spezified in z->scaff_W/_V**/
            if  (z->scaff_W[in_ar] != 0) {
                a->W[in_ar]  = d_matrix (a->d_r, A[in_ar].d_r);
                a->dW[in_ar] = d_matrix (a->d_r, A[in_ar].d_r);
            } else {
                a->W[in_ar]  = NULL;
                a->dW[in_ar] = NULL;
            }
        }




/******************************** init_a *************************************/
/* Initializes and normalizes weights, thresholds or sets to zero.           */
/* previously in vehicle.c                                                   */

void init_a (PARAMS *x, AGENT *Z, AREA *A, char dir[256]) {
    AREA  *a;
    AGENT *z;
    int ct_ar, in_ar, ct_n;
    char fullname[256];

    /**all weights W**/
    for (in_ar = 0; in_ar < x->areas; ++in_ar) {

        for (ct_ar = 0; ct_ar < x->areas; ++ct_ar) {

            a = A + ct_ar;
            z = Z + ct_ar;

            if  (z->scaff_W[in_ar] == 1) {

                /**initialize with random values**/
                init_matrix (a->W[in_ar], a->d_r, A[in_ar].d_r, -0.01, 0.01);
                                                               /***!!!***/

                /**normalize weights**
                for (ct_n = 0; ct_n < a->d_r; ++ct_n)
                    normalize (a->W[in_ar][ct_n], A[in_ar].d_r, 1.0);
                 !!! **/
            }

            if  (z->scaff_W[in_ar] == 2) {

                sprintf (fullname, "%s/obs_W_%d_%d.ext", dir, ct_ar, in_ar);

                /*if file.ext does not exist then take file.pnm*/
                if  (fopen (fullname, "r") == NULL)
                    sprintf (fullname, "%s/obs_W_%d_%d.pnm", dir, ct_ar, in_ar);

                importP36_matrix (a->W[in_ar],
                               z->d_a, z->d_b, (Z[in_ar].d_a), (Z[in_ar].d_b),
                               fullname);
            }

            if  (z->scaff_W[in_ar] == 3) {

                init_weights_arbor (a->W[in_ar],
                               z->d_a, z->d_b, (Z[in_ar].d_a), (Z[in_ar].d_b));
            }

            if  (z->scaff_W[in_ar] == 4) {

                init_matrix_zero (a->W[in_ar], a->d_r, A[in_ar].d_r);
            }

            if  (z->scaff_W[in_ar] == 5) {

                init_weights_arbor_color (a->W[in_ar],
                               z->d_a, z->d_b, (Z[in_ar].d_a), (Z[in_ar].d_b));
            }
        }
    }
}



/********************** previously not for pointer ***************************/



/************************* weight_initdelta **********************************/
void weight_initdelta (PARAMS *x, AGENT *Z, AREA *A) {

    int ct_ar, in_ar, ct_n;
    AREA *a;
    AGENT *z;

    for (ct_ar = 0; ct_ar < x->areas; ++ct_ar) {

        a = A + ct_ar;
        z = Z + ct_ar;

        for (in_ar = 0; in_ar < x->areas; ++in_ar) {

            if  (z->scaff_W[in_ar] != 0)
                init_matrix_zero (a->dW[in_ar], a->d_r, A[in_ar].d_r);

            if  (z->scaff_V[in_ar] != 0)
                init_matrix_zero (a->dV[in_ar], a->d_r, A[in_ar].d_r);

            if  (z->scaff_X[in_ar] != 0)
                init_matrix_zero (a->dX[in_ar], a->d_r, A[in_ar].d_r);

            if  (z->scaff_Y[in_ar] != 0)
                init_matrix_zero (a->dY[in_ar], a->d_r, A[in_ar].d_r);
        }

        /**no weights U, U_delta yet!!**/

        if  (z->ch_Theta)
            for (ct_n = 0; ct_n < a->d_r; ++ct_n)
                a->Theta_delta[ct_n] = 0.0;
    }
}


/************************** weight_update ************************************/
/* Add W_delta to the weight using parameter  z->eps_W.                      */

void weight_update (PARAMS *x, AGENT *Z, AREA *A) {

    int ct_ar, in_ar, ct_n, ct_in;
    AREA *a;
    AGENT *z;

    for (ct_ar = 0; ct_ar < x->areas; ++ct_ar) {

        a = A + ct_ar;
        z = Z + ct_ar;

        for (in_ar = 0; in_ar < x->areas; ++in_ar) {

            if  (z->scaff_W[in_ar] != 0)
                for (ct_n = 0; ct_n < a->d_r; ++ct_n)
                    for (ct_in = 0; ct_in < A[in_ar].d_r; ++ct_in)
                        a->W[in_ar][ct_n][ct_in]
                        += z->eps_W[in_ar] * a->dW[in_ar][ct_n][ct_in];

            if  (z->scaff_V[in_ar] != 0)
                for (ct_n = 0; ct_n < a->d_r; ++ct_n)
                    for (ct_in = 0; ct_in < A[in_ar].d_r; ++ct_in)
                        a->V[in_ar][ct_n][ct_in]
                        += z->eps_V[in_ar] * a->dV[in_ar][ct_n][ct_in];

            if  (z->scaff_X[in_ar] != 0)
                for (ct_n = 0; ct_n < a->d_r; ++ct_n)
                    for (ct_in = 0; ct_in < A[in_ar].d_r; ++ct_in)
                        a->X[in_ar][ct_n][ct_in]
                        += z->eps_X[in_ar] * a->dX[in_ar][ct_n][ct_in];

            if  (z->scaff_Y[in_ar] != 0)
                for (ct_n = 0; ct_n < a->d_r; ++ct_n)
                    for (ct_in = 0; ct_in < A[in_ar].d_r; ++ct_in)
                        a->Y[in_ar][ct_n][ct_in]
                        += z->eps_Y[in_ar] * a->dY[in_ar][ct_n][ct_in];
        }

        /**no weights U, U_delta yet!!**/

        if  (z->ch_Theta)
            for (ct_n = 0; ct_n < a->d_r; ++ct_n)
                a->Theta[ct_n] += z->eps_Theta * a->Theta_delta[ct_n];
    }
}

