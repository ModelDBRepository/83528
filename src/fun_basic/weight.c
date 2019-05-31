#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "iter.h"
#include "series.h"
#include "utils.h"      /**for init_matrix_zero, d_tensor**/
#include "vehicle.h"    /**for print_command to debug**/
#include "local.h"      /**for local_mean_01 in weight_rec_mean_01**/

/************** functions pointed to by cmd->weightfunc **********************/

/****************************** weight_hebb **********************************/
/* Updates dW_target of neuron 5) ct_n  from cmd->area==n_from1[0].          */
/* Pre: cmd->S_from2; Post (here): cmd->S_from1[0]. (NICHT kommutativ!)      */
/* Stepsize weighted by cmd->moment.                                         */

void weight_hebb (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

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


/****************************** weight_antihebb ******************************/
/* As above but different sign.                                              */

void weight_antihebb (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

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

            dW[ct_in] -= cmd->moment * Pre[ct_in] * post;

        if  (z->ch_Theta)

            A[area].Theta_delta[ct_n] += cmd->moment * post;
    }
}

/****************************** weight_kohonen *******************************/
/* Like weight_hebb, but Pre[] replaced by (Pre[]-W[]).                      */
/* Expects neighborhood function as post.                                    */

void weight_kohonen (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

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



/****************************** weight_heuristic *****************************/

void weight_heuristic (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int inarea, ct_l, ct_in;
    double *dW, *Pre;

    const double post = cmd->S_from1[0][ct_t][ct_n];    /**assumed as "here"**/

    if  ((cmd->anz_from1 != 1) || (cmd->n_from1[0] != cmd->area))
        fprintf (stderr, "wrong use of weight_heuristic");

    for (ct_l = 0; ct_l < cmd->anz_from2; ++ct_l) {

        inarea = cmd->n_from2[ct_l];

        dW = cmd->dW_target[inarea][ct_n];
        /* W_delta = A[area].W_delta[inarea][ct_n]; */

        Pre = cmd->S_from2[ct_l][ct_t];

        for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)

            dW[ct_in] += cmd->moment * (cmd->quantum[0][0] * cmd->quantum[0][0] - Pre[ct_in] * post)
                                * fabs (cmd->quantum[0][0] * cmd->quantum[0][0] - Pre[ct_in] * post);
    }
}



/****************************** weight_log ***********************************/
/* From weight_hebb, but rule is dW = log ...                                */
/* Relatively stupid idea without long-lasting trace ... to stop divergence  */
/* through large growth below, small decay above the correct value.          */
/* Trivial (i.e. no effect) if only two discrete activation values are used. */

/* ** dist = q[1][0]: negative distance of the log-singularity to the origin *
   **                 in multiples of pp = q[0][1]^2.                        *
   ** dist=0: normal log;  dist>>1: linear function like original Foldiak    *
   gnuplot
   dist = 0.1
   p   = 0.03
   pp  = p*p
   set xrange [0:2*pp]
   plot 0, pp-x, -log((x+pp*dist)/(pp*(1+dist)))*pp*(1.0+dist)
*/

void weight_log (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int inarea, ct_l, ct_in;
    double *dW, *Pre;

    const double post = cmd->S_from1[0][ct_t][ct_n];    /**assumed as "here"**/

    if  ((cmd->anz_from1 != 1) || (cmd->n_from1[0] != cmd->area))
        fprintf (stderr, " wrong use of weight_log ");
    if  (cmd->anz_quantums != 2)
        fprintf (stderr, " weight_log wants 2 parameters! ");
    if  (cmd->quantum[0][0] == 0.0)
        fprintf (stderr, " weight_log wants the sign given ");

    for (ct_l = 0; ct_l < cmd->anz_from2; ++ct_l) {

        inarea = cmd->n_from2[ct_l];

        dW = cmd->dW_target[inarea][ct_n];
        /* W_delta = A[area].W_delta[inarea][ct_n]; */

        Pre = cmd->S_from2[ct_l][ct_t];

        for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in) {

            double psquare = cmd->quantum[0][0] * cmd->quantum[0][0];

            dW[ct_in] -= cmd->moment * log( (Pre[ct_in] * post + psquare * cmd->quantum[1][0])    /**<-- shift of the argument**/
                                          / (psquare * (1.0 + cmd->quantum[1][0])))               /**<-- move back the zero-point to psquare**/
                                     * psquare * (1.0 + cmd->quantum[1][0]);                      /**<-- scale to have slope 1 at zero-point**/

/*          **for comparison, the original Foldiak rule**
            dW[ct_in] += cmd->moment * Pre[ct_in] * post - cmd->quantum[0][0] * cmd->quantum[0][0];
*/
/*          **like Foldiak, but parameter q[1][0] to increase positive growth (a bend in the linear function)**
            if  ((cmd->quantum[0][0] * cmd->quantum[0][0] - Pre[ct_in] * post) < 0.0)
                dW[ct_in] += cmd->moment * (cmd->quantum[0][0] * cmd->quantum[0][0] - Pre[ct_in] * post);     **negative growth**
            else
                dW[ct_in] += cmd->moment * (cmd->quantum[0][0] * cmd->quantum[0][0] - Pre[ct_in] * post) * cmd->quantum[1][0];  **(more) positive growth**
*/
        }

        if  (z->ch_Theta)

            fprintf (stderr, " no theta rule for weight_log ");
    }
}


/****************************** weight_ica ***********************************/
/* Updates dW_target of neuron 5) ct_n  from cmd->area==n_from1[0]=n_from1[1]*/
/* Here: Post: cmd->S_from1[0] and all inner activation cmd->S_from1[1].     */
/* No use of cmd->S_from2 (but cmd->n_from2 as source area).                 */
/* Stepsize weighted by cmd->moment.                                         */
/* Attention0: Works after use of transfer function "local_mean_01" only!    */
/* Attention1: Do not use within the same sweep as functions for activations!*/
/* Attention2: Do not subtract CM of single image data point (turns vector!)!*/

void weight_ica (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    static int firsttime = 1;
    static int d_here;
    static int d_in;
    static double *C;
    static double *D;
    double **Wmat, *dW;
    int i, j;

    const int inarea  = cmd->n_from2[0];
    const int area    = cmd->area;
    const double post = cmd->S_from1[0][ct_t][ct_n];
    const double *U   = cmd->S_from1[1][ct_t];

    if  (firsttime) {
        d_here    = z->d_a * z->d_b;
        d_in      = z->d_a * z->d_b;
        C         = d_vector (d_here);
        D         = d_vector (d_here);
        firsttime = 0;
    }

    if  (d_here != z->d_a * z->d_b)
        fprintf (stderr, "\ncan't use weight_ica twice\n");
    if  (d_here != A[inarea].d_r)
        fprintf (stderr, "\narea mismatch in weight_ica\n");
    if  (cmd->anz_from1 != 2)
        fprintf (stderr, "\nwrong use of weight_ica\n");
    if  (cmd->anz_from2 != 1)
        fprintf (stderr, "\nwrong use of weight_ica\n");
    if  ((cmd->n_from1[0] != area) || (cmd->n_from1[1] != area))
        fprintf (stderr, "\nwrong areas in weight_ica\n");

    Wmat = cmd->W_target[inarea];
    dW   = cmd->dW_target[inarea][ct_n];

    for (j = 0; j < d_in; ++j) {
        D[j] = 0.0;
        for (i = 0; i < d_here; ++i)
            D[j] += U[i] * Wmat[i][j]; /*rueck-summierte Ausgabe von Quelle j*/
        D[j] *= 1.0 - 2.0 * post;     /*multipl. mit Erschoepfung des Ziels i*/
        dW[j] += cmd->moment * (Wmat[ct_n][j] + D[j]);
    }

    if  (z->ch_Theta)
        A[area].Theta_delta[ct_n] -= cmd->moment * (1.0 - 2.0 * post);
                               /* * A[area].Theta[ct_n] * A[area].Theta[ct_n]*/
}


/****************************** weight_rec_mean_01 ***************************/
/* For recurrent lateral weights; recursive rule like recurrent backprop.    */
/* 1:o {W; weight_rec_mean_01;, ; (U^mu-U^l), U(t-1); $onemineta, 8+1+0+0+1} */
/* inputs: weights (as target), difference, u(t), (1-eta), transfkt params   */
/* sideff: weights changed, static r[i][k][l] := d/dw_kl u^l[i](t) changed   */
/* Uses fixed transfer function: local_mean_01, and its derivative!          */
/* ATTENTION: Use me only in one context!                                    */

/* k <-> ct_n  ;  l <-> ct_in */

void weight_rec_mean_01 (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int i, j, ct_in;
    int inarea = cmd->n_from1[0];
    int d_r = A[cmd->area].d_r;
    double *dW = cmd->dW_target[inarea][ct_n];
    double  *W = cmd->W_target[inarea][ct_n];
    double (*func)(double *, double, double)      = local_mean_01;
    double (*func_diff)(double *, double, double) = local_mean_01_diff;
    double par[5];
    static double *w_ftick_oneminuseta;
    static double ***r;
    static double ***r_new;
    static int firsttime = 1;
    static int OK = 0;

    if  (firsttime) {
        r     = d_tensor (d_r, d_r, d_r);
        r_new = d_tensor (d_r, d_r, d_r);
        w_ftick_oneminuseta = d_vector (d_r);
        firsttime = 0;
    }

    par[0] = cmd->quantum[1][0];
    par[1] = cmd->quantum[1][1];
    par[2] = cmd->quantum[1][2];
    par[3] = cmd->quantum[1][3];
    par[4] = cmd->quantum[1][4];

    if  (inarea != cmd->area)
        fprintf (stderr, "\narea error in weight_rec_mean_01!\n");

    if  (ct_t == cmd->quantum[2][0]) {                                                              /**!!!**/
        for (i = 0; i < d_r; ++i)
            for (ct_in = 0; ct_in < d_r; ++ct_in)
                r[i][ct_n][ct_in] = /* (i == ct_n)
                                  ? func (par, cmd->S_from2[0][ct_t][ct_in], 0.0)
                                  : */ 0.0;
        OK = 1;

        if  (ct_n == 0)
            fprintf (stderr, "\nweight_rec_mean_01: 1-eta=%.2f, par=%.2f %.2f %.2f %.2f %.2f, rset=%d",
                         cmd->quantum[0][0], par[0], par[1], par[2], par[3], par[4], (int)cmd->quantum[2][0]);
    }

    if  (!OK)
        fprintf (stderr, "\ncould not init r in weight_rec_mean_01 !\n");

    for (j = 0; j < d_r; ++j)
        w_ftick_oneminuseta[j] = W[j]
                               * func_diff (par, cmd->S_from2[0][ct_t - 1][j], 0.0)
                               * cmd->quantum[0][0];

    /**loop over target neurons (k) is outside of this function, in a stay**/
    /**all ingoing weights (l)**/
    for (ct_in = 0; ct_in < d_r; ++ct_in) {

        /**"output" neurons where later will be the error (i)**/
        for (i = 0; i < d_r; ++i) {

            r_new[i][ct_n][ct_in] = 0.0;
            for (j = 0; j < d_r; ++j)
                r_new[i][ct_n][ct_in] += w_ftick_oneminuseta[j]
                                       * r[j][ct_n][ct_in];

            if  (i == ct_n)
                r_new[i][ct_n][ct_in] += func (par, cmd->S_from2[0][ct_t - 1][ct_in], 0.0);

            if  (ct_t > cmd->quantum[2][0])                                                               /**!!!**/
                dW[ct_in] += cmd->S_from1[0][ct_t][i] * r_new[i][ct_n][ct_in];
        }
    }

    for (ct_in = 0; ct_in < d_r; ++ct_in)
        for (i = 0; i < d_r; ++i)
            r[i][ct_n][ct_in] = r_new[i][ct_n][ct_in];
}






/****************************** weight_copy_topo *****************************/
/* Copies activation vector to the lateral weights receptive fields.         */

void weight_copy_topo (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int a, b, a_eff, b_eff, ct_a, ct_b;
    double *W, *Act;

    if  ((cmd->anz_from1 != 1) || (cmd->anz_from2 != 0))
        fprintf (stderr, "\nweight_copy_topo: too many areas\n");
    if  (cmd->n_from1[0] != cmd->area)
        fprintf (stderr, "\nweight_copy_topo: only self-connections!\n");

    W  = cmd->W_target[0][ct_n];
    Act = cmd->S_from1[0][ct_t];

    ct_a = ct_n / z->d_b;
    ct_b = ct_n % z->d_b;

    for (a = 0; a < z->d_a; ++a) {
        a_eff = (a - ct_a + z->d_a) % z->d_a;
        for (b = 0; b < z->d_b; ++b) {
            b_eff = (b - ct_b + z->d_b) % z->d_b;
            W[a * z->d_b + b] = Act[a_eff * z->d_b + b_eff];
        }
    }
}


/****************************** weight_add_diag ******************************/
/* Adds activation vector to the diagonal of the (square) weight matrix.     */

void weight_add_diag (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    double *W, *Act;
    static int firstcount = 0;
    int inarea;

    if  ((cmd->anz_from1 != 1) || (cmd->anz_from2 != 0))
        fprintf (stderr, "\nweight_add_diag: too many areas\n");
    if  (cmd->n_from1[0] != cmd->area)
        fprintf (stderr, "\nweight_add_diag now no self-connections!\n");

    inarea  = cmd->n_from1[0];

    W  = cmd->W_target[inarea][ct_n];    /**inarea was previously 0**/
    Act = cmd->S_from1[0][ct_t];         /**first index like ct_l**/

    if  (cmd->quantum[0][0] == 1)
        if  (firstcount < z->d_a * z->d_b) {
            W[ct_n] += Act[ct_n];
            firstcount++;
        }

    if  (cmd->quantum[0][0] == 0)
        W[ct_n] += Act[ct_n];
}


/****************************** dweight_self_set *****************************/
/* Set self-connections of dW to q[0][0] (init learn value of autapses).     */

void dweight_self_set (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    if  ((cmd->anz_from1 != 1) || (cmd->anz_from2 != 1))
        fprintf (stderr, "\ndweight_self_set: too many areas\n");
    if  ((cmd->n_from1[0] != cmd->area) || (cmd->n_from2[0] != cmd->area))
        fprintf (stderr, "\ndweight_self_set: only self-connections!\n");

    cmd->dW_target[cmd->area][ct_n][ct_n] = cmd->quantum[0][0];
}


/****************************** weight_self_zero *****************************/
/* Sets diagonal of self-connecting weight matrix (autapses) to zero.        */

void weight_self_zero (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    if  ((cmd->anz_from1 != 1) || (cmd->anz_from2 != 0))
        fprintf (stderr, "\nweight_self_zero: too many areas\n");
    if  (cmd->n_from1[0] != cmd->area)
        fprintf (stderr, "\nweight_self_zero: only self-connections!\n");

    cmd->W_target[cmd->area][ct_n][ct_n] = 0.0;
}


/****************************** weight_seung *********************************/
/* Updates dW_target of neuron 5) ct_n  from cmd->area==n_from1[0].          */
/* Pre: cmd->S_from2; Post (here): cmd->S_from1[0]. (NICHT kommutativ!)      */
/* Stepsize weighted by cmd->moment.                                         */

void weight_seung (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int inarea, ct_l, ct_in;
    double *dW, *W, *Pre;

    const double post = cmd->S_from1[0][ct_t][ct_n];    /**assumed as "here"**/

    if  ((cmd->anz_from1 != 1) || (cmd->n_from1[0] != cmd->area))
        fprintf (stderr, "wrong use of weight_hebb");

    for (ct_l = 0; ct_l < cmd->anz_from2; ++ct_l) {

        inarea = cmd->n_from2[ct_l];

        dW = cmd->dW_target[inarea][ct_n];
        W  = cmd->W_target[inarea][ct_n];

        Pre = cmd->S_from2[ct_l][ct_t];

        for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)

            dW[ct_in] += W[ct_in] * (cmd->moment * Pre[ct_in] * post - 1.0); /**originally batch: w->Pre*post*w now additive: w->w + (Pre*post - 1)*w**/

        if  (z->ch_Theta)

            fprintf (stderr, "\nweight_seung: cannot deal with Theta's\n");
    }
}





/****************************** weight_respcorr ******************************/
/* Linear response correction for mean-field BM (additive to Hebb).          */

void weight_respcorr (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int inarea, ct_l, ct_in;
    double *dW, *W, *Pre;

    const double post = cmd->S_from1[0][ct_t][ct_n];    /**assumed as "here"**/

    if  ((cmd->anz_from1 != 1) || (cmd->n_from1[0] != cmd->area))
        fprintf (stderr, "wrong use of weight_respcorr");

    for (ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        inarea = cmd->n_from1[ct_l];

        dW = cmd->dW_target[inarea][ct_n];
        W  = cmd->W_target[inarea][ct_n];

        Pre = cmd->S_from1[ct_l][ct_t];

        for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)

            dW[ct_in] += W[ct_in]    * (1.0 - Pre[ct_in] * Pre[ct_in])
                       * cmd->moment * (1.0 - post * post);
    }
}

/****************************** weight_antirespcorr **************************/
/* As above but different sign.                                              */

void weight_antirespcorr (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int inarea, ct_l, ct_in;
    double *dW, *W, *Pre;

    const double post = cmd->S_from1[0][ct_t][ct_n];    /**assumed as "here"**/

    if  ((cmd->anz_from1 != 1) || (cmd->n_from1[0] != cmd->area))
        fprintf (stderr, "wrong use of weight_antirespcorr");

    for (ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        inarea = cmd->n_from1[ct_l];

        dW = cmd->dW_target[inarea][ct_n];
        W  = cmd->W_target[inarea][ct_n];

        Pre = cmd->S_from1[ct_l][ct_t];

        for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)

            dW[ct_in] -= W[ct_in]    * (1.0 - Pre[ct_in] * Pre[ct_in])
                       * cmd->moment * (1.0 - post * post);
    }
}



/****************************** weight_decay *********************************/
/* Decrease W_delta.                                                         */
/* q[0][.] decay parameters (individual for each input area n_from2[..]).    */
/* q[1][0] =0 pos&neg together; =1 only pos; else(=-1) only neg weights.     */
/* q[2][0] =0 global; =1 local.                                              */
/* q[3][0] =0 no act; =1 mult by S_from1[0][ct_t][ct_n].                     */
/* q[4][0] =0 normal equation; =1 quad; =2 decay is just -w.                 */

void weight_decay (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int inarea, ct_l, ct_in;
    double *dW, *W, dec_W;
    double tot_square = 0.0;

    const int sign_only = cmd->quantum[1][0] == 0.0 ? 0 :
                         (cmd->quantum[1][0] == 1.0 ? 1 : -1);
    const int nonglob   = cmd->quantum[2][0] == 0.0 ? 0 : 1;
    const double post   = cmd->quantum[3][0] == 0.0 ? 1.0
                        : cmd->S_from1[0][ct_t][ct_n];     /**assumed as "here"**/
    const int select    = (int)(cmd->quantum[4][0]);

    if  (cmd->anz_quantums != 5)
        fprintf (stderr, "\n5 parameters for weight_decay, please!\n");
    if  (cmd->anz_quant[0] != cmd->anz_from2)
        fprintf (stderr, "\none weight_decay strength for each input area!\n");
    if  ((cmd->anz_from1 != 1) || (cmd->n_from1[0] != cmd->area))
        fprintf (stderr, "wrong use of weight_decay");

    if  (! nonglob)
        for (ct_l = 0; ct_l < cmd->anz_from2; ++ct_l) {

            inarea = cmd->n_from2[ct_l];
            W  = cmd->W_target[inarea][ct_n];

            if  (sign_only != 0) {
                for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)
                    if  (sign_only * W[ct_in] > 0)
                        tot_square += W[ct_in] * W[ct_in];
            } else {
                for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)
                    tot_square += W[ct_in] * W[ct_in];
            }
        }

    for (ct_l = 0; ct_l < cmd->anz_from2; ++ct_l) {

        inarea = cmd->n_from2[ct_l];
        dW = cmd->dW_target[inarea][ct_n];
        W  = cmd->W_target[inarea][ct_n];
        /* W_delta = A[cmd->area].W_delta[inarea][ct_n]; */
        /* W       = A[cmd->area].W      [inarea][ct_n]; */

        if  (nonglob) {
            tot_square = 0.0;
            if  (sign_only != 0) {
                for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)
                    if  (sign_only * W[ct_in] > 0)
                        tot_square += W[ct_in] * W[ct_in];
            } else {
                for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)
                    tot_square += W[ct_in] * W[ct_in];
            }
        }

        dec_W = cmd->quantum[0][ct_l];
                              /**here stood inarea instead of ct_l**/

        /**the normal term**/
        if  (select == 0) {
            if  (sign_only != 0) {
                for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)
                    if  (sign_only * W[ct_in] > 0)
                        dW[ct_in] -= cmd->moment * W[ct_in] * tot_square * post
                                   * dec_W;
            } else {
                for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)
                                      /** ^ must not be ct_l**/
                        dW[ct_in] -= cmd->moment * W[ct_in] * tot_square * post
                                   * dec_W;
            }
        }

        /**the quadratic term**/
        if  (select == 1) {
            if  (sign_only != 0) {
                for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)
                    if  (sign_only * W[ct_in] > 0)
                        dW[ct_in] -= cmd->moment * W[ct_in] * tot_square * post
                                   * dec_W /** fabs(W[ct_in]) * fabs(W[ct_in]) !!!*/;
            } else {
                for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)
                    dW[ct_in] -= cmd->moment * W[ct_in] * tot_square * post
                               * dec_W /** fabs(W[ct_in]) * fabs(W[ct_in]) !!!*/;
            }
        }
        /**decay is just: -w ; can also be activity dependent (*post)**/
        if  (select == 2) {
            for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)
                dW[ct_in] -= cmd->moment * W[ct_in] * dec_W * post;
        }

        /**decay is a constant ; (note the sign!; may be act dependent)**/
        if  (select == 3) {
            for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)
                dW[ct_in] -= cmd->moment * dec_W * post;
        }
    }
}



/****************************** weight_decay_topo ****************************/
/* Decrease W_delta propto topographic act-scaffold S_from1 on _this_ area.  */
/* Scaffold function must be 0-centered and includes all parameters.         */
/* q[0][0] =0 pos&neg together; =1 only pos; else(=-1) only neg weights.     */

void weight_decay_topo (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int inarea, ct_in, diff;
    double *dW, *W;
    double *scaffold;
    const int sign_only = cmd->quantum[0][0] == 0.0 ? 0 :
                         (cmd->quantum[0][0] == 1.0 ? 1 : -1);

    if  (cmd->anz_quantums != 1)
        fprintf (stderr, "\none parameter for weight_decay_topo, please!\n");
    if  ((cmd->anz_from1 != 1) || (cmd->n_from1[0] != cmd->area))
        fprintf (stderr, "wrong use of weight_decay_topo");

    scaffold = cmd->S_from1[0][0]; /**assumed: here, now**/

    inarea = cmd->area;
    dW = cmd->dW_target[inarea][ct_n];
    W  = cmd->W_target[inarea][ct_n];

    for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in) {
        diff = abs (ct_n - ct_in);

        if  (sign_only == 0)
            dW[ct_in] -= cmd->moment * W[ct_in] * scaffold[diff];

        if  (sign_only == 1)
            if  (W[ct_in] > 0.0)
                dW[ct_in] -= cmd->moment * W[ct_in] * scaffold[diff];

        if  (sign_only == -1)
            if  (W[ct_in] < 0.0)
                dW[ct_in] -= cmd->moment * W[ct_in] * scaffold[diff];
    }
}



/****************************** weight_incr_glob *****************************/
/* Anti weight decay.                                                        */

void weight_incr_glob (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int inarea, ct_l, ct_in;
    double tot_square = 0.0;
    double *dW, *W, dec_W;

    for (ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        inarea = cmd->n_from1[ct_l];

        tot_square += quad_length (cmd->W_target[inarea][ct_n], A[inarea].d_r);
    }

    if  (tot_square != 0.0)
        for (ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

            inarea = cmd->n_from1[ct_l];

            dW = cmd->dW_target[inarea][ct_n];
            W  = cmd->W_target[inarea][ct_n];

            dec_W = cmd->quantum[0][inarea];

            for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)

                dW[ct_in] += dec_W * W[ct_in] / tot_square;
        }
    else
        fprintf (stderr, " not yet weight_incr ");
/***
    dW = 0.0002; dV = 0.1
    f(x,y) = dW * x * y * y * abs(x); g(x,y) = dV * x / (y * y)
    set xrange [ 0.0:1.0]; set yrange [ 1.0:10.0]; set isosamples 25, 25
    set cntrparam levels discrete 0.0;; set contour;; splot -f(x,y) + g(x,y)
***/
}


/****************************** weight_incr_pos_glob *************************/
/* Anti weight decay for positive weights only:                              */
/* neglect negative weights, also in calculating length.                     */

void weight_incr_pos_glob (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int inarea, ct_l, ct_in;
    double tot_square = 0.0;
    double *dW, *W, dec_W;

    for (ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        inarea = cmd->n_from1[ct_l];
        W  = cmd->W_target[inarea][ct_n];

        for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)
            if  (W[ct_in] > 0.0)
                tot_square += W[ct_in] * W[ct_in];
    }

    if  (tot_square != 0.0)
        for (ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

            inarea = cmd->n_from1[ct_l];
            dW = cmd->dW_target[inarea][ct_n];
            W  = cmd->W_target[inarea][ct_n];

            dec_W = cmd->quantum[0][inarea];

            for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)
                if  (W[ct_in] >= 0.0)
                    dW[ct_in] += dec_W * W[ct_in] / tot_square;
        }
    else
        fprintf (stderr, " not yet weight_incr_pos ");
}



/****************************** weight_act_decay *****************************/
/* Decrease all dW dependent on the activation of a cell.                    */
/* Parameters z->dec_W can, however, be chosen individually for each inarea. */
/* Stepsize weighted by cmd->moment.                                         */

void weight_act_decay (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int inarea, ct_l, ct_in;
    double *dW, *W, dec_W;
    const double post = cmd->S_from1[0][ct_t][ct_n];    /**assumed as "here"**/

    if  ((cmd->anz_from1 != 1) || (cmd->n_from1[0] != cmd->area))
        fprintf (stderr, "wrong use of weight_act_decay");

    for (ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        inarea = cmd->n_from1[ct_l];

        dW = cmd->dW_target[inarea][ct_n];
        W  = cmd->W_target[inarea][ct_n];

        dec_W = cmd->quantum[0][inarea];

        if  (post != 0.0)

            for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)

                dW[ct_in] -= dec_W * W[ct_in] /* W[ct_in]*W[ct_in]*/
                                   * post * cmd->moment;

        if  (z->ch_Theta)
            fprintf (stderr, "\nno decay for A[area].Theta_delta[ct_n]!\n\n");
    }
}



/****************************** weight_norm_glob *****************************/
/* Normalize total (all areas) input of weight vector A->W of neuron 5) ct_n */
/* to length z->dec_W. All dec_W for this area must better be equal!         */

void weight_norm_glob (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int ct_l, inarea, ct_in;
    double *W, dec_W, tot_square = 0.0;

    for (ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        inarea = cmd->n_from1[ct_l];

        tot_square += quad_length (cmd->W_target[inarea][ct_n], A[inarea].d_r);
    }

    if  (tot_square != 0.0) {

        tot_square = sqrt (tot_square);

        for (ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

            inarea = cmd->n_from1[ct_l];

            dec_W = cmd->quantum[0][ct_l];    /**was cmd->quantum[0][inarea which is wrong!**/

            W  = cmd->W_target[inarea][ct_n];

            for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)

                W[ct_in] *= dec_W / tot_square;
        }
    } else {
        fprintf (stderr, " tot_square=0 in weight_norm_glob ");
    }
}


/****************************** weight_norm_each *****************************/
/* Normalize each area input of weight vector A->W of neuron 5) ct_n         */
/* to its length z->dec_W[inarea]. Thus, different inputs are "independent". */

void weight_norm_each (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int ct_l, inarea;
    double dec_W;

    for (ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        inarea = cmd->n_from1[ct_l];

        dec_W = cmd->quantum[0][inarea];

        normalize (cmd->W_target[inarea][ct_n], A[inarea].d_r, dec_W);
               /**nasty bug here ^ (was ct_l) and here ^    **/
    }
}


/****************************** weight_rectify *******************************/
/* If parameters 0,1: set all negative weights A->W/V/X/Y of neuron 5) to 0. */

void weight_rectify (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int ct_l, inarea, in_n;

    if  (cmd->anz_quantums != 2)
        fprintf (stderr, "\nuse weight_rectify e.g. like 0, 1\n");
    if  (cmd->quantum[1][0] == 0.0)
        fprintf (stderr, "\nweight rectify wants 1 or -1\n");
    if  (cmd->anz_arguments != 2)
      fprintf (stderr, "\nweight rectify wants area and inarea (like local)");

    for (ct_l = 0; ct_l < cmd->anz_from2; ++ct_l) {  /**<-- was n_from1[ct_l] error!!!**/

        inarea = cmd->n_from2[ct_l];  /**<-- was n_from1[ct_l] error!!!**/

        for (in_n = 0; in_n < A[inarea].d_r; ++in_n)

            if  (cmd->quantum[1][0] == 1.0) {
                if  (cmd->W_target[inarea][ct_n][in_n] < cmd->quantum[0][0])
                    cmd->W_target[inarea][ct_n][in_n] = cmd->quantum[0][0];
            } else { /**-1.0**/
                if  (cmd->W_target[inarea][ct_n][in_n] > cmd->quantum[0][0])
                    cmd->W_target[inarea][ct_n][in_n] = cmd->quantum[0][0];
            }
    }
}


/****************************** weight_rect_eps ******************************/
/* Set all negative weights A->W/V/X/Y of neuron 5) to small random values.  */

void weight_rect_eps (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int ct_l, inarea, in_n;

    for (ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        inarea = cmd->n_from1[ct_l];

        for (in_n = 0; in_n < A[inarea].d_r; ++in_n)

            if  (cmd->W_target[inarea][ct_n][in_n] < 0.0)

                cmd->W_target[inarea][ct_n][in_n] = 0.000001 * drand48();
    }
}


/****************************** weight_add_rand ******************************/
/* Add random value to weights A->W/V/X/Y of neuron 5).                      */

void weight_add_rand (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int ct_l, inarea, in_n;
    const double fac = cmd->quantum[0][0];

    if  (cmd->anz_quantums != 1)
        fprintf (stderr, "\n1 parameter for weight_add_rand, please!\n");

    for (ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        inarea = cmd->n_from1[ct_l];

        for (in_n = 0; in_n < A[inarea].d_r; ++in_n)

            cmd->W_target[inarea][ct_n][in_n] += (drand48() - 0.5) * fac;
    }
}


/****************************** weight_bound *********************************/
/* Set bound for each weight of A->W/V/X/Y of neuron 5).                     */

void weight_bound (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

  fprintf (stderr, "\nuse weight_rectify instead of weight_bound!\n");
/*
    int ct_l, inarea, in_n;
    double dec_W;

    for (ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        inarea = cmd->n_from1[ct_l];

        dec_W = cmd->quantum[0][ct_l];
                               *here stood inarea instead of ct_l!*

        for (in_n = 0; in_n < A[inarea].d_r; ++in_n)

            if  (cmd->W_target[inarea][ct_n][in_n] > dec_W)
                cmd->W_target[inarea][ct_n][in_n] = dec_W;
            else
                if  (cmd->W_target[inarea][ct_n][in_n] < -dec_W)
                    cmd->W_target[inarea][ct_n][in_n] = -dec_W;
    }
*/
}


/****************************** weight_copy **********************************/
/* Copies A->W going to inarea(s) to here inputs. For symmetric weights.     */

void weight_copy (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int inarea, ct_l, ct_in;
    double *W_area, **W_inarea;

    for (ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        inarea = cmd->n_from1[ct_l];

        /**point to inputs from inarea to here: ct_n**/
        W_area = cmd->W_target[inarea][ct_n];
            /* = A[area].W[inarea][ct_n]; */

        /**point to toinarea-fromthisarea-weightmatrix**/
        W_inarea = cmd->W_from1[ct_l];
              /* = A[inarea].W[area]; */

        for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)

            W_area[ct_in] = W_inarea[ct_in][ct_n];
    }
}



/****************************** weight_hack_init *****************************/
/* Multiplies W_target of neuron 5) ct_n  from cmd->area==n_from1[0]         */
/* by a factor given by q[0][0] | q[0][1]; all weighted by cmd->moment.      */
/* Attention: Can only be called once!!                                      */

void weight_hack_init (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int inarea, ct_l, ct_in;
    double *W;

    const int area = cmd->area;

    const double factor = (ct_n < z->d_a * z->d_b / 3)                /**!!!**/
                        ? cmd->quantum[0][0] : cmd->quantum[0][1];

    static int firsttime = 1;
    static int *first_this;

    if  (firsttime) {
        int i;
        first_this = i_vector (A[area].d_r);
        for (i = 0; i < A[area].d_r; ++i)
            first_this[i] = 1;
        firsttime  = 0;
    }

    if  (first_this[ct_n] == 0)
        return;
    first_this[ct_n] = 0;

    fprintf (stderr, " %d ", ct_n);

    if  ((cmd->anz_from1 != 1) || (cmd->n_from1[0] != cmd->area))
        fprintf (stderr, "wrong use of weight_hack_slow");

    if  (cmd->moment != 1.0) {
        fprintf (stderr, "weight_hack_init: cmd->moment must be 1!");
        fprintf (stderr, " But ct_t does not play a role");
    }

    for (ct_l = 0; ct_l < cmd->anz_from2; ++ct_l) {

        inarea = cmd->n_from2[ct_l];

        W = cmd->W_target[inarea][ct_n];

        for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)

            W[ct_in] *= factor;

        if  (z->ch_Theta)

            A[area].Theta[ct_n] *= factor;
    }
}

/****************************** weight_hack_slow *****************************/
/* Multiplies dW_target of neuron 5) ct_n  from cmd->area==n_from1[0]        */
/* by a factor given by q[0][0] | q[0][1];  All devided by cmd->moment.      */

void weight_hack_slow (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int inarea, ct_l, ct_in;
    double *dW;

    const int area = cmd->area;

    const double factor = (ct_n < z->d_a * z->d_b / 3)                /**!!!**/
                        ? cmd->quantum[0][0] : cmd->quantum[0][1];

    if  ((cmd->anz_from1 != 1) || (cmd->n_from1[0] != cmd->area))
        fprintf (stderr, "wrong use of weight_hack_slow");

    if  (cmd->moment != 1.0) {
        fprintf (stderr, "weight_hack_slow: cmd->moment must be 1!");
        fprintf (stderr, " But ct_t does not play a role");
    }

    for (ct_l = 0; ct_l < cmd->anz_from2; ++ct_l) {

        inarea = cmd->n_from2[ct_l];

        dW = cmd->dW_target[inarea][ct_n];

        for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)

            dW[ct_in] *= factor;

        if  (z->ch_Theta)

            A[area].Theta_delta[ct_n] *= factor;
    }
}


/****************************** weight_hack_slow_inv *************************/
/* Analogeous to weight_hack_slow but for Weights of opposite direction.     */
/* I.e. hacked parameter changes among inarea, not here.                     */

void weight_hack_slow_inv (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int inarea, ct_l, ct_in;
    double *dW;

    if  ((cmd->anz_from1 != 1) || (cmd->n_from1[0] != cmd->area))
        fprintf (stderr, "wrong use of weight_hack_slow_inv");

    if  (cmd->moment != 1.0) {
        fprintf (stderr, "weight_hack_slow_inv: cmd->moment must be 1!");
        fprintf (stderr, " But ct_t does not play a role");
    }

    for (ct_l = 0; ct_l < cmd->anz_from2; ++ct_l) {

        inarea = cmd->n_from2[ct_l];

        dW = cmd->dW_target[inarea][ct_n];

        for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in) {

            const double factor = (ct_in < A[inarea].d_r / 3)         /**!!!**/
                        ? cmd->quantum[0][0] : cmd->quantum[0][1];

            dW[ct_in] *= factor;
        }
    }
}



/****************************** init_theta_const *****************************/
/* Should be use in a sweep with ilen==elen==1. For initialization only.     */

void init_theta_const (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    if  (! z->ch_Theta)
      fprintf (stderr, "\n\ninit_theta_const should have ch_Theta != 0\n\n");

    A[cmd->area].Theta[ct_n] = cmd->quantum[0][0];
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
       if  (cmd->b_target == 'W')
           cmd->W_target  = A[cmd->area].W;
******************************************************************************/




/************** other functions for weights not for pointer ******************/



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

