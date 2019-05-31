#include <stdio.h>
#include <stdlib.h>    /**for total_four**/
#include <math.h>
#include <string.h>    /**for phoneme2array**/
#include <time.h>      /**for total_pause**/
#include <unistd.h>    /**for total_pause; usleep**/
#include "../kernel/coco.h"
#include "../kernel/series.h"
#include "../kernel/utils.h"
#include "local.h"     /**for local_atan_to_2pi used in total_retina_to_SC**/



/* This is to be used for all neurons together in a stay with mode total!*/



/****************************** total_cut_peaks ******************************/
/* Acts at a local max/min are cut to their next extreme neighbor's value.   */
/* Inarea size must = target area size -- doesn't check!                     */
/* S_TARGET SHOULD BE DIFFERENT FROM S_FROM !                                */

DOUBLE total_cut_peaks (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int area = cmd->area;
    int d_b  = A[area].d_b;

    for (int X = 0; X < A[area].d_a; X++)
        for (int Y = 0; Y < d_b; Y++) {

            DOUBLE min_neigh = cmd->S_from1[0][ct_t][X * d_b + Y] + 99.9;
            DOUBLE max_neigh = cmd->S_from1[0][ct_t][X * d_b + Y] - 99.9;

            for (int XX = X-1; XX <= X+1; ++XX)
            for (int YY = Y-1; YY <= Y+1; ++YY)
            if  ((XX >= 0) && (YY >= 0) && (XX < A[area].d_a) && (YY < d_b) && (XX != X) && (YY != Y)) {

                min_neigh = (cmd->S_from1[0][ct_t][XX * d_b + YY] < min_neigh) ? cmd->S_from1[0][ct_t][XX * d_b + YY] : min_neigh;
                max_neigh = (cmd->S_from1[0][ct_t][XX * d_b + YY] > max_neigh) ? cmd->S_from1[0][ct_t][XX * d_b + YY] : max_neigh;
            }

            cmd->S_target[ct_t][X * d_b + Y] = cmd->S_from1[0][ct_t][X * d_b + Y];    /**default: just copy**/

            if  (cmd->S_from1[0][ct_t][X * d_b + Y] > max_neigh)
                cmd->S_target[ct_t][X * d_b + Y] = max_neigh;

            if  (cmd->S_from1[0][ct_t][X * d_b + Y] < min_neigh)
                cmd->S_target[ct_t][X * d_b + Y] = min_neigh;
        }
 
    return (DOUBLE)(0);
}



/****************************** total_rectify_save ***************************/
/* Remembers values 1st time, then cuts values if higher (lower).            */
/* q[0][0] = -1 don't allow higher; = 1 don't allow smaller; = 0 do nothing  */
/* Could instead be done with local_biggerof, but found this more convenient */

DOUBLE total_rectify_save (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    static int firsttime = 1;
    static DOUBLE *save_vals = NULL;
    int area = cmd->area;
    int dim = A[area].d_n;

    if  (firsttime) {
        save_vals = d_vector (dim);

        for (int i = 0; i < dim; ++i)
            save_vals[i] = cmd->S_from1[0][ct_t][i];

        firsttime = 0;

    } else {

        if  (cmd->quantum[0][0] == -1)
            for (int i = 0; i < dim; ++i)
                if  (cmd->S_from1[0][ct_t][i] > save_vals[i])
                    cmd->S_target[ct_t][i] = save_vals[i];
                else
                    cmd->S_target[ct_t][i] = cmd->S_from1[0][ct_t][i];

        if  (cmd->quantum[0][0] == 1)
            for (int i = 0; i < dim; ++i)
                if  (cmd->S_from1[0][ct_t][i] < save_vals[i])
                    cmd->S_target[ct_t][i] = save_vals[i];
                else
                    cmd->S_target[ct_t][i] = cmd->S_from1[0][ct_t][i];
    }

    return (DOUBLE)(0);
}



/****************************** total_sub_mean *******************************/
/* S_target will equal S_from1 - mean (of all activations).                  */

DOUBLE total_sub_mean (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int area = cmd->area;

    for (int i = 0; i < A[area].d_n; ++i)
        cmd->S_target[ct_t][i] = cmd->S_from1[0][ct_t][i];

    sub_mean_vector (cmd->S_target[ct_t], A[area].d_n);

    return (DOUBLE)(0);
}


/****************************** total_sub_mean_col ***************************/
/* Partitions area into 3 subsections. For color retina.                     */
/* S_target will equal S_from1 - mean (of all activations).                  */

DOUBLE total_sub_mean_col (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int area = cmd->area;
    int retina_size = A[area].d_n / 3;

    for (int i = 0; i < A[area].d_n; ++i)
        cmd->S_target[ct_t][i] = cmd->S_from1[0][ct_t][i];

    sub_mean_vector (cmd->S_target[ct_t], retina_size);
    sub_mean_vector (cmd->S_target[ct_t] + retina_size, retina_size);
    sub_mean_vector (cmd->S_target[ct_t] + 2 * retina_size, retina_size);

    return (DOUBLE)(0);
}



/****************************** total_average ****************************/
/* All activations on S_target are assigned the mean of acts on S_from1. */
/* May be different areas.                                               */

DOUBLE total_average (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int inarea = cmd->n_from1[0];
    int area = cmd->area;
    DOUBLE sp = 0.0;

    for (int i = 0; i < A[inarea].d_n; ++i)
        sp += cmd->S_from1[0][ct_t][i];
    sp /= (DOUBLE)(A[inarea].d_n);
 
    for (int i = 0; i < A[area].d_n; ++i)
        cmd->S_target[ct_t][i] = sp;

    return (DOUBLE)(0);
}



/****************************** total_sub_min ********************************/
/* S_target will equal S_from1 - min (of all activations).                   */

DOUBLE total_sub_min (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int area = cmd->area;
    DOUBLE min = 99999.9;

    for (int i = 0; i < A[area].d_n; ++i)
        min = cmd->S_from1[0][ct_t][i] < min ? cmd->S_from1[0][ct_t][i] : min;

    for (int i = 0; i < A[area].d_n; ++i)
        cmd->S_target[ct_t][i] = cmd->S_from1[0][ct_t][i] - min;

    return (DOUBLE)(0);
}



/****************************** total_set_mean *******************************/
/* S_target will all be the mean of all neurons activations (fixed time).    */

DOUBLE total_set_mean (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int i;
    double mean = 0.0;
    int area = cmd->area;
    int inarea = cmd->n_from1[0];

    for (i = 0; i < A[inarea].d_n; ++i)
        mean += cmd->S_from1[0][ct_t][i];
    mean /= (double)(A[inarea].d_n);

    for (i = 0; i < A[area].d_n; ++i)
        cmd->S_target[ct_t][i] = mean;

    return (DOUBLE)(0);
}


/****************************** total_time_mean ******************************/
/* S_target will be the local mean of last q[0][0] activations of the sweep. */

DOUBLE total_time_mean (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int i, t;
    const int t_time   = (int)(cmd->quantum[0][0]);
    int area = cmd->area;

    if  (ct_t + 1 < t_time)
        fprintf (stderr, "\ntotal_time_mean: no start before mean-time\n");

    if  (cmd->S_target == cmd->S_from1[0])
        fprintf (stderr, "\ntotal_time_mean: target must differ from source!");

    for (i = 0; i < A[area].d_n; ++i)
        cmd->S_target[ct_t][i] = 0.0;

    for (t = ct_t + 1 - t_time; t < ct_t + 1; ++t)
        for (i = 0; i < A[area].d_n; ++i)
            cmd->S_target[ct_t][i] += cmd->S_from1[0][t][i];

    for (i = 0; i < A[area].d_n; ++i)
        cmd->S_target[ct_t][i] /= (double)t_time;

    return (DOUBLE)(0);
}


/******************************** all_time_mean ******************************/
/* S_target at all times will be the local mean of all act's of the sweep.   */
/* (-> non causality by considering values of the future!)                   */
/* !alltime sweep only!                                                      */

DOUBLE all_time_mean (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {

    int i, t;
    int area = cmd->area;

    for (i = 0; i < A[area].d_n; ++i)
        cmd->S_target[begin][i] = cmd->S_from1[0][begin][i];

    for (t = begin + 1; t < end; ++t)
        for (i = 0; i < A[area].d_n; ++i)
            cmd->S_target[begin][i] += cmd->S_from1[0][t][i];

    for (i = 0; i < A[area].d_n; ++i)
        cmd->S_target[begin][i] /= (double)(end - begin);

    for (t = begin + 1; t < end; ++t)
        for (i = 0; i < A[area].d_n; ++i)
            cmd->S_target[t][i] = cmd->S_target[begin][i];

    return (DOUBLE)(0);
}


/****************************** total_set_attime *****************************/
/* Set S_target to q[1][0] after q[0][0] invocations of this function.       */
/* Use only once per script because of counter!                              */
/* MIND (OR EVEN USE) THE FACT THAT RELAXATIONS INVOKE IT SEVERAL TIMES !    */

DOUBLE total_set_attime (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

   static int count = 0;
    int area = cmd->area;

   if  (count == (int)(cmd->quantum[0][0])) {

       fprintf (stderr, "\ntotal_set_time setting values to %f.    ",
                cmd->quantum[1][0]);

       for (int i = 0; i < A[area].d_n; ++i)
           cmd->S_target[ct_t][i] = cmd->quantum[1][0];
   }

   count += 1;

    return (DOUBLE)(0);
}


/****************************** total_epoche_mean ****************************/
/* S_target will be the running mean over items of epoche; "decay" q[0][0].  */
/* Local on neurons and same time within each relaxation.                    */
/* ATTENTION: Don't use S_target anywhere else!                              */

DOUBLE total_epoche_mean (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int i;
    double lambda = cmd->quantum[0][0];
    int area = cmd->area;

    for (i = 0; i < A[area].d_n; ++i)
        cmd->S_target[ct_t][i] = (1.0 - lambda) * cmd->S_target[ct_t][i]
                                      + lambda  * cmd->S_from1[0][ct_t][i];

    return (DOUBLE)(0);
}


/****************************** total_epoche_x_square ************************/
/* S_target will be the running mean squared activation over items of epoche.*/
/* Local on neurons and same time within each relaxation.                    */
/* ATTENTION: Don't use S_target anywhere else!                              */
/* Used to compute the kurtosis.                                             */

DOUBLE total_epoche_x_square (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int i;
    double lambda = cmd->quantum[0][0];
    int area = cmd->area;

    for (i = 0; i < A[area].d_n; ++i)
        cmd->S_target[ct_t][i] = (1.0 - lambda) * cmd->S_target[ct_t][i]
                                      + lambda  * cmd->S_from1[0][ct_t][i] * cmd->S_from1[0][ct_t][i];

    return (DOUBLE)(0);
}

/****************************** total_epoche_x_fourth ************************/
/* S_target will be the running mean x^4 activation over items of epoche.    */
/* Local on neurons and same time within each relaxation.                    */
/* ATTENTION: Don't use S_target anywhere else!                              */
/* Used to compute the kurtosis.                                             */

DOUBLE total_epoche_x_fourth (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int i;
    double lambda = cmd->quantum[0][0];
    int area = cmd->area;

    for (i = 0; i < A[area].d_n; ++i)
        cmd->S_target[ct_t][i] = (1.0 - lambda) * cmd->S_target[ct_t][i]
                                      + lambda  * cmd->S_from1[0][ct_t][i] * cmd->S_from1[0][ct_t][i] * cmd->S_from1[0][ct_t][i] * cmd->S_from1[0][ct_t][i];

    return (DOUBLE)(0);
}

/****************************** total_epoche_kurt ****************************/
/* Get kurtosis conveniently from averages: S_from1=fourth, S_from2=square.  */
/* Assumes that mean is zero.                                                */

DOUBLE total_epoche_kurt (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int i;
    int area = cmd->area;

    for (i = 0; i < A[area].d_n; ++i)
        cmd->S_target[ct_t][i] = cmd->S_from1[0][ct_t][i] / (cmd->S_from2[0][ct_t][i] * cmd->S_from2[0][ct_t][i]);

    return (DOUBLE)(0);
}



/****************************** total_import *********************************/
/* ct_t has to match the imported file size! Thus use at e.g. ct_t = rlen.   */
/* q[0][0]=2: reads from pipe and writes feedback-pipe.                      */
/*
DOUBLE total_import (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {
    char fullname[256];
    static int firsttime = 1;
    static char *tmp_uname;
    if  (firsttime) {
        tmp_uname = get_tmp_uname ();
        firsttime = 0;
    }
    if  ((cmd->quantum[0][0] == 2) || (cmd->quantum[0][0] == 3)) {
        FILE *pp;
        **pipe**
        sprintf (fullname, "%s/obs_%c_%d.pipe", tmp_uname, cmd->b_target, cmd->area);
        importP36_matrix (cmd->S_target, 1, ct_t, A[area].d_a, A[area].d_b, fullname);
        **feedback-pipe**
        if  (cmd->quantum[0][0] == 2) {
            sprintf (fullname, "%s/obs_%c_%d.back", tmp_uname, cmd->b_target, cmd->area);
            pp = fopen (fullname, "w");
            fprintf (pp, "9");
            fclose (pp);
        }
    } else {
        sprintf (fullname, "%s/obs_%c_%d.pnm", tmp_uname, cmd->b_target, cmd->area);
        importP36_matrix (cmd->S_target, 1, ct_t, A[area].d_a, A[area].d_b, fullname);  **see vehicle.c for import of Theta's**
    }

    return (DOUBLE)(0);
}
*/


/****************************** total_back_pipe ******************************/
/* q[0][0]=3: writes feedback-pipe. (When total_import does not do it.)      */
/* q[0][0]=2: writes forward-pipe. (Don't use together with total_import!)   */
/*
DOUBLE total_back_pipe (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {
    char fullname[256], halfname[256];
    FILE *pp;
    static int firsttime = 1;
    static char *tmp_uname;
    if  (firsttime) {
        tmp_uname = get_tmp_uname ();
        firsttime = 0;
    }
    if  (cmd->anz_quantums != 2) {
        fprintf (stderr, "\ntotal_back_pipe now wants 2 quantums\n");
        exit (1);
    }
    if  (cmd->quantum[1][0] == 0)
        sprintf (halfname, "%s/obs_%c_%d", tmp_uname, cmd->b_target, cmd->area);
    else
        sprintf (halfname, "%s/tcl", tmp_uname);
    **feedback-pipe**
    if  (cmd->quantum[0][0] == 3) {
        sprintf (fullname, "%s.back", halfname);
        pp = fopen (fullname, "w");
        fprintf (pp, "9");
        fclose (pp);
    }
    **forward-pipe**
    if  (cmd->quantum[0][0] == 2) {
        sprintf (fullname, "%s.pipe", halfname);
        pp = fopen (fullname, "r");
        fscanf (pp, "9");
        fclose (pp);
    }
    return (DOUBLE)(0);
}
*/


/****************************** total_copy_back ******************************/
/* Set S_target of all recent times to now-value of S_from1.                 */

DOUBLE total_copy_back (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int i, t;
    int area = cmd->area;

    for (t = 0; t <= ct_t; ++t)
        for (i = 0; i < A[area].d_n; ++i)
            cmd->S_target[t][i] = cmd->S_from1[0][ct_t][i];

    return (DOUBLE)(0);
}

/****************************** total_embed **********************************/
/* Embeds one area's acts into the middle of another. "total_copy, _cut_at"  */
/* If target area is larger, then zero padding at edge.                      */
/*          "       smaller,  "   cutoff of outer rim of S_from area.        */
/* Size differences must be even, so that edge is same on either sides.      */
/* (see also total_cut_at, below)                                            */

DOUBLE total_embed (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int area   = cmd->area;
    int inarea = cmd->n_from1[0];
    int edge_a = A[area].d_a - A[inarea].d_a;
    int edge_b = A[area].d_b - A[inarea].d_b;

    if  ((edge_a % 2 != 0) || (edge_b % 2 != 0))
        fprintf (stderr, "\ntotal_embed wants even size differences!\n");
    if  (edge_a * edge_b < 0)
        fprintf (stderr, "\ntotal_embed wants same-sign size differences!\n");

    edge_a /= 2;
    edge_b /= 2;

    /**works for both, target area smaller or larger than input area**/
    for (int X_area = 0; X_area < A[area].d_a; ++X_area)
        for (int Y_area = 0; Y_area < A[area].d_b; ++Y_area) {
            int X_inarea = X_area - edge_a;
            int Y_inarea = Y_area - edge_b;
            if  ((X_inarea >= 0) && (X_inarea < A[inarea].d_a)  && (Y_inarea >= 0) && (Y_inarea < A[inarea].d_b))
                cmd->S_target[ct_t][X_area * A[area].d_b + Y_area] = cmd->S_from1[0][ct_t][X_inarea * A[inarea].d_b + Y_inarea];
            else
                cmd->S_target[ct_t][X_area * A[area].d_b + Y_area] = 0.0;
        }

    return (DOUBLE)(0);
}

/****************************** total_embed_tiles ****************************/
/* Embeds one area's acts at (0,0) coord of another.                         */
/* If target area is larger, then repeats 2-dim act pattern (hence, tiles).  */
/*          "       smaller   "   cutoff.                                    */
/* (see also: _copy, _cut_at; total_set_all is special case, if unit 0 used) */

DOUBLE total_embed_tiles (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int area   = cmd->area;
    int inarea = cmd->n_from1[0];

    for (int X_area = 0; X_area < A[area].d_a; ++X_area)
        for (int Y_area = 0; Y_area < A[area].d_b; ++Y_area) {
            int X_inarea = X_area % A[inarea].d_a;
            int Y_inarea = Y_area % A[inarea].d_b;

            cmd->S_target[ct_t][X_area * A[area].d_b + Y_area] = cmd->S_from1[0][ct_t][X_inarea * A[inarea].d_b + Y_inarea];
        }

    return (DOUBLE)(0);
}

/****************************** total_spread_time ****************************/
/* Copys act's from a small area over t-times to a t-fold larger area at t=0.*/

DOUBLE total_spread_time (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int t, i, k;
    const int inarea   = cmd->n_from1[0];
    const int t_time   = (int)(cmd->quantum[0][0]);
    int area = cmd->area;

    if  (A[area].d_n != A[inarea].d_n * t_time)
        fprintf (stderr, "\ntotal_spread_time: wrong area sizes!\n");

    if  (ct_t < t_time)
        fprintf (stderr, "\ntotal_spread_time: no start before integr-time\n");

    k = 0;
    for (t = ct_t - t_time; t < ct_t; ++t) {
        for (i = 0; i < A[inarea].d_n; ++i) {
            cmd->S_target[ct_t][k] = cmd->S_from1[0][t][i];
            k += 1;
        }
    }
 /* for (t = 0; t < t_time; ++t)                         */
 /*     ...                                              */
 /*         cmd->S_target[0][k] = cmd->S_from1[0][t][i]; */

    return (DOUBLE)(0);
}


/****************************** total_spread_xy ******************************/
/* Spread act's from a 1-dim input area to the 2nd dim along target area.    */
/* q[0][0] == 1: spread along x (d_a); then must d_b=A[inarea].d_n           */
/* q[0][0] == 2: spread along y (d_b); then must d_a=A[inarea].d_n           */

DOUBLE total_spread_xy (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    const int inarea   = cmd->n_from1[0];
    int i, j;
    int area = cmd->area;

    if  (cmd->quantum[0][0] == 0)
        fprintf (stderr, "\ntotal_spread_xy wants your decision!\n");
    if  (cmd->quantum[0][0] == 1)
        if  (A[inarea].d_n != A[area].d_b)
            fprintf (stderr, "\ntotal_spread_xy: input dim must match d_b!\n");
    if  (cmd->quantum[0][0] == 2)
        if  (A[inarea].d_n != A[area].d_a)
            fprintf (stderr, "\ntotal_spread_xy: input dim must match d_a!\n");

    for (i = 0; i < A[area].d_a; ++i) {
        for (j = 0; j < A[area].d_b; ++j) {
            int k = (cmd->quantum[0][0] == 1) ? j : i;
            cmd->S_target[ct_t][i*A[area].d_b + j] = cmd->S_from1[0][ct_t][k];
        }
    }

    return (DOUBLE)(0);
}


/****************************** total_set_all ********************************/
/* Set all activations to that of unit q[0][0] of S_from1.                   */

DOUBLE total_set_all (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    const int inarea = cmd->n_from1[0];
    int area = cmd->area;

    int in_neuron = (int)(cmd->quantum[0][0]);

    if  (in_neuron >= A[inarea].d_n)
        fprintf (stderr, "\n\ntotal_set_all: argument out of range\n\n");

    for (int i = 0; i < A[area].d_n; ++i)
        cmd->S_target[ct_t][i] = cmd->S_from1[0][ct_t][in_neuron];

    return (DOUBLE)(0);
}


/****************************** total_collapse_xy ****************************/
/* Sum act's over 2nd dim of input area write to a 1-dim target area.        */
/* q[0][0] == 1: collapse x (d_a); then must d_b=A[inarea].d_n               */
/* q[0][0] == 2: collapse y (d_b); then must d_a=A[inarea].d_n               */
/* q[1][0] == d_a(in_area) which is elsewise not available. d_b not needed.  */

DOUBLE total_collapse_xy (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    const int inarea   = cmd->n_from1[0];
    int i = 0, j = 0, in_d_a, in_d_b;
    int area = cmd->area;

    if  (cmd->anz_quantums != 2) {
        fprintf (stderr, "\ntotal_collapse_xy wants 2 quantums\n");
        exit (1);
    }
    if  (cmd->anz_quant[1] != 2) {
        fprintf (stderr, "\ntotal_collapse_xy wants 2 quant[1]\n");
        exit (1);
    }

    in_d_a = (int)(cmd->quantum[1][0]);
    in_d_b = (int)(cmd->quantum[1][1]);

    if  (cmd->quantum[0][0] == 0)
        fprintf (stderr, "\ntotal_collapse_xy wants your decision!\n");
    if  (cmd->quantum[0][0] == 1)
        if  (in_d_b != A[area].d_b)
            fprintf (stderr, "\ntotal_collapse_xy: d_b's must match\n");
    if  (cmd->quantum[0][0] == 2)
        if  (in_d_a != A[area].d_a)
            fprintf (stderr, "\ntotal_collapse_xy: d_a's must match\n");
    if  (in_d_a * in_d_b != A[inarea].d_n)
        fprintf (stderr, "\ntotal_collapse_xy got inconsistent dimensions\n");

   if  (cmd->quantum[0][0] == 1) /**sum over d_a**/
       for (j = 0; j < A[area].d_b; ++j) {
           cmd->S_target[ct_t][j] = 0.0;
           for (i = 0; i < A[area].d_a; ++i)
               cmd->S_target[ct_t][j] += cmd->S_from1[0][ct_t][i * in_d_b + j];
       }

   if  (cmd->quantum[0][0] == 2) /**sum over d_b**/
       for (i = 0; i < A[area].d_a; ++i) {
           cmd->S_target[ct_t][i] = 0.0;
           for (j = 0; j < in_d_b; ++j)
               cmd->S_target[ct_t][i] += cmd->S_from1[0][ct_t][i * in_d_b + j];
       }

    return (DOUBLE)(0);
}




/****************************** total_mean_to_one ****************************/
/* S_target will scaled so that mean=1.                                      */
/* This is a useful equivalent of total_spherize if all acts are >= 0        */
/* From total_set_mean.                                                      */

DOUBLE total_mean_to_one (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int i;
    double mean = 0.0;
    int area = cmd->area;

    for (i = 0; i < A[area].d_n; ++i)
        mean += cmd->S_from1[0][ct_t][i];
    mean /= (double)(A[area].d_n);

    if  (mean == 0.0)
        mean = 0.0000001;                  /** !! UNBELIEVABLY DIRTY !! **/

    for (i = 0; i < A[area].d_n; ++i)
        cmd->S_target[ct_t][i] /= mean;

    if  (cmd->anz_pointers == 1)
        cmd->pointers[0]->float_val = mean;  /**same function as "scale" in total_spherize**/

    return (DOUBLE)(0);
}


/****************************** total_spherize *******************************/
/* Spherizes (make element-variance 1). Does NOT subtract mean. Do it first! */
/* cmd-pointer's float_val will be value by which elements have been divided.*/
/* q[0][0] threshold under which NOT to spherize.                            */

DOUBLE total_spherize (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int i;
    int area = cmd->area;

    for (i = 0; i < A[area].d_n; ++i)
        cmd->S_target[ct_t][i] = cmd->S_from1[0][ct_t][i];

    DOUBLE scale = spherize_vector (cmd->S_target[ct_t], A[area].d_n);

    /**don't spherize (reverse the effects) if |S_from1| < q[0][0]**/
    if  (quad_length (cmd->S_from1[0][ct_t], A[area].d_n) < cmd->quantum[0][0])
        for (i = 0; i < A[area].d_n; ++i)
            cmd->S_target[ct_t][i] = cmd->S_from1[0][ct_t][i];

    if  (cmd->anz_pointers == 1)
        cmd->pointers[0]->float_val = scale;

    return (DOUBLE)(0);
}


/****************************** total_normalize *******************************/
/* Normalizes to length quantum[0][0]. Uses quad_length. Subtract mean first! */

DOUBLE total_normalize (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int i;
    int area = cmd->area;

    for (i = 0; i < A[area].d_n; ++i)
        cmd->S_target[ct_t][i] = cmd->S_from1[0][ct_t][i];

    normalize (cmd->S_target[ct_t], A[area].d_n, cmd->quantum[0][0]);

    return (DOUBLE)(0);
}


/****************************** total_scale_to_max ****************************/
/* Multiplies by a factor so that the maximum is q[0][0].                     */

DOUBLE total_scale_to_max (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int i;
    DOUBLE scale_factor, max = 0.0;
    int area = cmd->area;

    for (i = 0; i < A[area].d_n; ++i)
        max = cmd->S_from1[0][ct_t][i] > max ? cmd->S_from1[0][ct_t][i] : max;

    scale_factor = cmd->quantum[0][0] / max;

    for (i = 0; i < A[area].d_n; ++i)
        cmd->S_target[ct_t][i] = cmd->S_from1[0][ct_t][i] * scale_factor;

    return (DOUBLE)(0);
}


/****************************** total_threshold *******************************/
/* Sets values <= (q[0][0] * max) to zero.                                    */

DOUBLE total_threshold (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    DOUBLE max = 0.0;
    int area = cmd->area;

    for (int i = 0; i < A[area].d_n; ++i)
        max = cmd->S_from1[0][ct_t][i] > max ? cmd->S_from1[0][ct_t][i] : max;

    const DOUBLE threshold = cmd->quantum[0][0] * max;

    for (int i = 0; i < A[area].d_n; ++i)
        cmd->S_target[ct_t][i] = (cmd->S_from1[0][ct_t][i] > threshold)
                               ?  cmd->S_from1[0][ct_t][i] : 0.0;

    return (DOUBLE)(0);
}


/****************************** total_true_softmax ****************************/
/* P(i=ON) = exp(beta a_i) / (sum_j exp(beta a_j))                            */
/* Returns continuous values. (maybe use something such as local_as_rand)     */
/* q[0][0] = beta                                                             */

DOUBLE total_true_softmax (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int i;
    DOUBLE denominator = 0.0;
    int area = cmd->area;

    double beta = cmd->quantum[0][0];

    for (i = 0; i < A[area].d_n; ++i)
        denominator += exp (beta * cmd->S_from1[0][ct_t][i]);

    for (i = 0; i < A[area].d_n; ++i)
        cmd->S_target[ct_t][i] = exp (beta * cmd->S_from1[0][ct_t][i]) / denominator;

    return (DOUBLE)(0);
}

/****************************** total_true_softmax_row ************************/
/* q[0][0] = beta                                                             */

DOUBLE total_true_softmax_row (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int area = cmd->area;

    double beta = cmd->quantum[0][0];

    for (int X = 0; X < A[area].d_a; X++) {
        DOUBLE denominator = 0.0;
        for (int Y = 0; Y < A[area].d_b; Y++)
            denominator += exp (beta * cmd->S_from1[0][ct_t][X * A[area].d_b + Y]);

        for (int Y = 0; Y < A[area].d_b; Y++)
            cmd->S_target[ct_t][X * A[area].d_b + Y] = exp (beta * cmd->S_from1[0][ct_t][X * A[area].d_b + Y]) / denominator;
    }

    return (DOUBLE)(0);
}



/****************************** total_four ***********************************/
/* Fourier-transforms source, overwrites target of this AND next time step!  */
/* Attention: target must exist for rlen=2!                                  */

DOUBLE total_four (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

     int i, j;
     int OK = -1;
     static int d_a_save;
     static int d_b_save;
     static int firsttime = 1;
     static float *c_fl_field;
     int direction;
    int area = cmd->area;

     if  (ct_t != 0)
         fprintf (stderr, "\ntotal_four: rlen=0! Target must exist at 0,1!\n");

     if  (  (cmd->anz_arguments != 1)
         || (cmd->anz_from1 != 1)
         || (cmd->area != cmd->n_from1[0])
         || (cmd->anz_quantums != 1)
         || (cmd->anz_quant[0] != 1))
         fprintf (stderr, "\nwrong arguments in total_four\n");

     for (i = 1; i <= 2048; i *= 2) {
         if  (A[area].d_a == i)
             OK += 1;
         if  (A[area].d_b == i)
             OK += 1;
     }
     if  (OK != 1)
         fprintf (stderr, "\n\ntotal_four: check dimensions to be 2^n !\n\n");


     if  (firsttime) {
         c_fl_field  = (float *)malloc (2 * A[area].d_n * sizeof (float));
         d_a_save = A[area].d_a;
         d_b_save = A[area].d_b;
         firsttime = 0;
         fprintf (stderr, "\nReminder: total_four target must exist at rlen 0 and 1!\n");
     }

     if  ((A[area].d_a != d_a_save) || (A[area].d_b != d_b_save))
         fprintf (stderr, "\ntotal_four used with one dim for ever only\n");

     if  ((cmd->quantum[0][0] != -1) && (cmd->quantum[0][0] != 1))
         fprintf (stderr, "\ntotal_four wants useful direction given\n");

     direction = (int)(cmd->quantum[0][0]);

     /**assume only real; fill complex with zero's**/
     if  (direction == 1)
//         make_complex_do_fl (cmd->S_from1[0][ct_t], c_fl_field, A[area].d_a, A[area].d_b);
;

     /**assume real and complex in target[ct_t] and target[ct_t + 1]**/
     if  (direction == -1)
         for (i = 0, j = 0; i < A[area].d_n; i+=1, j+=2) {

             /**real values for the transform**/
             c_fl_field[j]     = (float)(cmd->S_from1[0][ct_t][i]);

             /**real values for the transform**/
             c_fl_field[j + 1] = (float)(cmd->S_from1[0][ct_t + 1][i]);
         }

//     do_fft (c_fl_field, A[area].d_a, A[area].d_b, direction);

     /**write real and complex into target[ct_t] and target[ct_t + 1]**/
     for (i = 0, j = 0; i < A[area].d_n; i+=1, j+=2) {

         /**real values of the transform**/
         cmd->S_target[ct_t][i]     = (double)(c_fl_field[j]);

         /**complex values of the transform**/
         cmd->S_target[ct_t + 1][i] = (double)(c_fl_field[j + 1]);
     }

    return (DOUBLE)(0);
}



/****************************** total_mult_compl *****************************/
/* Complex multiplication.                                                   */
/* act[ct_t=0]=real, act[ct_t+1]=imaginary in both target and sources.       */

DOUBLE total_mult_compl (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int area = cmd->area;
     int i;
     int      dim = A[area].d_n;
     DOUBLE **one = cmd->S_from1[0];
     DOUBLE **two = cmd->S_from2[0];

     if  (  (cmd->anz_arguments != 2)
         || (cmd->anz_from1 != 1)
         || (cmd->anz_from2 != 1)
         || (cmd->area != cmd->n_from1[0])
         || (ct_t != 0))
         fprintf (stderr, "\nwrong arguments in total_mult_compl\n");

     /**real result**/
     for (i = 0; i < dim; ++i)
         cmd->S_target[0][i] = one[0][i] * two[0][i] - one[1][i] * two[1][i];

     /**complex result**/
     for (i = 0; i < dim; ++i)
         cmd->S_target[1][i] = one[1][i] * two[0][i] + one[0][i] * two[1][i];

    return (DOUBLE)(0);
}



/****************************** total_absq_compl *****************************/
/* Absolute square value of complex. Target is real, act[ct=0] only.         */
/* act[ct_t=0]=real, act[ct_t+1]=imaginary in source.                        */

DOUBLE total_absq_compl (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int area = cmd->area;
     int i;
     int      dim = A[area].d_n;
     DOUBLE **one = cmd->S_from1[0];

     if  (  (cmd->anz_arguments != 1)
         || (cmd->anz_from1 != 1)
         || (cmd->area != cmd->n_from1[0])
         || (ct_t != 0))
         fprintf (stderr, "\nwrong arguments in total_absq_compl\n");

     /**result (real)**/
     for (i = 0; i < dim; ++i)
         cmd->S_target[0][i] = one[0][i] * one[0][i] + one[1][i] * one[1][i];

    return (DOUBLE)(0);
}


/*************************** total_geom_vector_sum ***************************/
/* Computes, e.g. weighted mean orientation in area. Adds all vectors in 2D. */
/* S_from1 = ori (0..2PI) (scale from 0..PI to 0..2PI first!)                */
/* S_from2 = strength                                                        */
/* q[0][0] = 0: returns orientation of vector sum                            */
/* q[0][0] = 1/2: returns length/average of vector sum                       */

DOUBLE total_geom_vector_sum (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int area = cmd->area;
    int inarea = cmd->n_from1[0];

    if  (  (cmd->anz_from2 != 1)
        || (cmd->n_from2[0] != inarea))
        fprintf (stderr, "\nwrong arguments in total_geom_vector_sum\n");

    DOUBLE sumX = 0.0;
    DOUBLE sumY = 0.0;
    for (int i = 0; i < A[inarea].d_n; ++i) {
        DOUBLE ori = cmd->S_from1[0][ct_t][i];
        DOUBLE len = cmd->S_from2[0][ct_t][i];
        DOUBLE X = len * cos(ori);
        DOUBLE Y = len * sin(ori);
        sumX += X;
        sumY += Y;
    }

    DOUBLE dum;
    if  ((int)(cmd->quantum[0][0]) == 0)
        for (int i = 0; i < A[area].d_n; ++i)
            cmd->S_target[ct_t][i] = local_atan_to_2pi (&dum, sumY, sumX);

    if  ((int)(cmd->quantum[0][0]) == 1)
        for (int i = 0; i < A[area].d_n; ++i)
            cmd->S_target[ct_t][i] = sqrt (sumX*sumX + sumY*sumY);

    if  ((int)(cmd->quantum[0][0]) == 2)
        for (int i = 0; i < A[area].d_n; ++i)
            cmd->S_target[ct_t][i] = sqrt (sumX*sumX + sumY*sumY) / (DOUBLE)(A[inarea].d_n);

    return (DOUBLE)(0);
}


/****************************** total_zhang_m_fix ****************************/
/* Solves x = w * f(x) where w is determined as the mean of cmd->W_from1.    */
/* Parameters are the same as for local function f, here local_zhang.        */
/*
DOUBLE total_zhang_m_fix (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {
    double xl = 0.0, xh = 0.0, xm, yl, yh, ym;
    double sum_w = 0.0, dumm = 0.0;
    double (*func)(double *, double, double);
    int i;
    func = local_zhang;
    if  (cmd->anz_quantums == 2)
        if  (cmd->quantum[1][0] == 2)
            func = local_zhang2;
    **get mean(W); use only first neurons, as all fields are equal (neglecting the shift)**
    for (i = 0; i < A[area].d_n; ++i)
        sum_w += cmd->W_from1[0][0][i];
    fprintf (stderr, "\ntotal_zhang_m_fix: sum_w=%f meanW=%f", sum_w, sum_w / (double)(A[area].d_n));
    **get two starting points, f(xl) < 0 and f(xh) > 0**
    do  {
        xl -= 1.0;
        xh += 1.0;
        yl = sum_w * func (cmd->quantum[0], xl, dumm) - xl;
        yh = sum_w * func (cmd->quantum[0], xh, dumm) - xh;
    } while (yl * yh > 0.0);
    **confine the interval**
    do  {
        xm = 0.5 * (xl + xh);
        ym = sum_w * func (cmd->quantum[0], xm, dumm) - xm;
        if  (yl * ym > 0.0) {
            xl = xm;
            yl = sum_w * func (cmd->quantum[0], xl, dumm) - xl;
        }
        if  (yh * ym > 0.0) {
            xh = xm;
            yh = sum_w * func (cmd->quantum[0], xh, dumm) - xh;
        }
    } while (xh - xl > 0.0000001);
    fprintf (stderr, "  mean act: %f -> mean out: %f", 0.5 * (xl + xh), func (cmd->quantum[0], 0.5 * (xl + xh), dumm));
    if  ((fabs (yl) > 0.000001) || (fabs (yh) > 0.000001))
        fprintf (stderr, "\nWarning: yl=%f yh=%f are too large\n\n", yl, yh);
    **write the x-value into all S_target entries**
    for (i = 0; i < A[area].d_n; ++i)
         cmd->S_target[0][i] = 0.5 * (xl + xh);

    return (DOUBLE)(0);
}
*/





/******************************* total_orange ********************************/
/* From init_orange, but data-func didn't know inarea, move was impossible.  */
/* cmd->S_target (t=0 .. ct_t) initialized with orange, pasted into S_from1. */
/* Color area assumed if orange, i.e.: d_a 3x long and range: 0 .. 255.      */
/* q[0][0] == 1: new rand locat. == 2: location from last according to input.*/
/* q[1][0] == radius of orange / sigma of Gaussian.                          */
/* q[2][0] == 0: Gaussian (greyscale); == 1: orange (RGB format).            */

DOUBLE total_orange (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int a, b, count, d_r, height;
    static double ho_center;
    static double br_center;
    double radius = cmd->quantum[1][0];
    static double red;
    static double green;
    static double blue;
    int area = cmd->area;

    double C[39][3] = {/*R    G    B         orange-likeliness              from file   **/
                       {253, 133,  80},   /**brighter                       fruit_02.jpg**/
                       {247, 122,  30},   /**brighter**/
                       {244, 126,  39},   /**brighter**/
                       {241, 109,  24},   /**prototype**/
                       {238, 144, 106},   /**bright shiny spot**/
                       {230,  88,  16},   /**nice prototype**/
                       {227,  86,  14},   /**nice medium**/
                       {221,  92,  26},   /**nice prototype**/
                       {206,  69,  14},   /**darker**/
                       {203,  59,  25},   /**darker, reddish**/
                       {175,  66,  33},   /**darker**/
                       {161,  39,  34},   /**dark red-brown**/
                       {254, 182,  72},   /**pale                           orangegrapefruitgiftbox.jpg**/
                       {242, 147,  55},   /**brighter, prototype region**/
                       {240, 175,   0},   /**grapefruit?**/
                       {234, 146,  36},   /**prototype**/
                       {199,  86,   8},   /**darker region**/
                       {195, 124,  20},   /**darker pale**/

                       { 89,  48,  37},   /**shadow      snap1.ppm**/
                       {136,  57,  27},   /**  to  **/
                       {203, 102,  36},   /**  ..  **/
                       {238, 130,  41},   /**bright**/
                       {115,  54,  34},   /**shadow      snap2.ppm**/
                       {129,  57,  12},   /**  to  **/
                       {180,  85,  35},   /**  ..  **/
                       {231, 124,  53},   /**bright**/
                       {110,  47,  30},   /**shadow      snap3.ppm**/
                       {144,  66,  28},   /**  to  **/
                       {158,  77,  40},   /**  ..  **/
                       {202,  91,  46},   /**bright (not brightest)**/
                       { 85,  40,  19},   /**shadow      snap6.ppm**/
                       {153, 105,  39},   /**  to  **/
                       {198, 151,  60},   /**bright**/
                       { 77,  35,  15},   /**shadow      snap7.ppm**/
                       {125,  64,  35},   /**  to  **/
                       {194, 126,  69},   /**bright**/
                       {106,  56,  28},   /**shadow, other orange**/
                       {170, 119,  70},   /**  to  **/
                       {222, 177, 128}    /**bright**/
                      };

    if  (cmd->anz_quantums != 3)
        fprintf (stderr, "\n3 parameters for init_orange, please!\n");

    d_r = A[area].d_n; /**true also for color**/

    if  (cmd->quantum[2][0] == 1)
        height = A[area].d_a / 3;   /**because color**/
    else
        height = A[area].d_a;       /**will only be used for coose location, next line**/

    /**choose new random orange location**/
    if  (cmd->quantum[0][0] == 1.0) {
        ho_center = radius + drand48 () * ((double)(height) - radius - radius);    /**d_a || hoehe:  langsamer index**/
        br_center = radius + drand48 () * ((double)(A[area].d_b) - radius - radius);    /**d_b || breite: schneller index**/
    }

    /**choose orange location according to input area**/ /**<-- cannot be done with a data-function, that's why it is here**/
    if  (cmd->quantum[0][0] == 2.0) {

        int inarea;

        if  (cmd->anz_from2 != 1) {
            fprintf (stderr, "Dirty function total_orange wants here anz_from2 = 1!\n");
            exit (1);
        }
        inarea = cmd->n_from2[0];
        if  (A[inarea].d_n != 2) {
            fprintf (stderr, "\ntotal_orange with q[0][0]=2 wants two units in inarea, area %d!\n", inarea);
            exit (1);
        }
        ho_center -= cmd->S_from2[0][ct_t][0];
        br_center -= cmd->S_from2[0][ct_t][1];

fprintf (stderr, "\ntotal_orange: ho_center moved minus %f,  br_center moved minus %f ",
                                    cmd->S_from2[0][ct_t][0], cmd->S_from2[0][ct_t][1]);

        if  (ho_center < radius) {
            ho_center = radius;
            fprintf (stderr, " upper out ");
        }
        if  (br_center < radius) {
            br_center = radius;
            fprintf (stderr, " left out ");
        }
        if  (ho_center >= height - radius) {
            ho_center = height - 1 - radius;
            fprintf (stderr, " lower out ");
        }
        if  (br_center >= A[area].d_b - radius) {
            br_center = A[area].d_b - 1 - radius;
            fprintf (stderr, " right out ");
        }
    }


    /**draw orange**/
    if  (cmd->quantum[2][0] == 1.0) {

        /**choose color only if new _random_ orange location**/
        if  (cmd->quantum[0][0] == 1.0) {

            int proto1 = (int)(drand48() * 39);
            int proto2 = (int)(drand48() * 39);
            double colmix = drand48();
            red    = colmix * C[proto1][0] + (1.0 - colmix) * C[proto2][0];
            green  = colmix * C[proto1][1] + (1.0 - colmix) * C[proto2][1];
            blue   = colmix * C[proto1][2] + (1.0 - colmix) * C[proto2][2];
        }

        /**draw orange disc (for all times) and maintain other background**/
        for (count = 0; count <= ct_t; ++count) {
            for (a = 0; a < height; ++a)
            for (b = 0; b < A[area].d_b; ++b)
            if  (((double)(a)-ho_center)*((double)(a)-ho_center) + ((double)(b)-br_center)*((double)(b)-br_center) < (radius*radius))
                cmd->S_target[count][a*A[area].d_b + b] = red;
            else
                cmd->S_target[count][a*A[area].d_b + b] = cmd->S_from1[0][count][a*A[area].d_b + b];

            for (a = height; a < 2 * height; ++a)
            for (b = 0; b < A[area].d_b; ++b)
            if  (((double)(a-height)-ho_center)*((double)(a-height)-ho_center) + ((double)(b)-br_center)*((double)(b)-br_center)
                < (radius*radius))
                cmd->S_target[count][a*A[area].d_b + b] = green;
            else
                cmd->S_target[count][a*A[area].d_b + b] = cmd->S_from1[0][count][a*A[area].d_b + b];

            for (a = 2 * height; a < 3 * height; ++a)
            for (b = 0; b < A[area].d_b; ++b)
            if  (((double)(a-2*height)-ho_center)*((double)(a-2*height)-ho_center) + ((double)(b)-br_center)*((double)(b)-br_center)
                < (radius*radius))
                cmd->S_target[count][a*A[area].d_b + b] = blue;
            else
                cmd->S_target[count][a*A[area].d_b + b] = cmd->S_from1[0][count][a*A[area].d_b + b];
        }

    } else {  /**draw greyscale(!) Gaussian**/

        for (count = 0; count <= ct_t; ++count) {
            for (a = 0; a < A[area].d_a; ++a)
            for (b = 0; b < A[area].d_b; ++b)
                    cmd->S_target[count][a*A[area].d_b + b]
                = 1.0 * exp (- (((double)a-ho_center)*((double)a-ho_center)+((double)b-br_center)*((double)b-br_center))
                            /  (2.0 * radius * radius));
        }
    }

    return (DOUBLE)(0);
}



/****************************** total_cosinus ********************************/
/* Sine waves, one centered (_in), one surround (_out). Later: move in time. */
/* $cmx+$cmy, 100      1.0+0,    $freq+0   $phase+0   $ori+0 0.0+0.0  1      */
/* center     diameter amplitude frequency phase(t=0) angle  velocity sphere */
/* See init_image_cos for a wrapper function that randomises some parameters */

DOUBLE total_cosinus (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    double cm_x, amplitude_in,  frequency_in,  phase_in,  angle_in,  velocity_in,
           cm_y, amplitude_out, frequency_out, phase_out, angle_out, velocity_out;
    double sq_diameter, sphere, sq_dist, directed_dist, diffA, diffB,
           e_A_in, e_B_in, e_A_out, e_B_out;
    int X, Y;
    int area = cmd->area;

    if  (cmd->anz_quantums != 8)
        fprintf (stderr, "\nwrong no of quantums in total_cosinus");

    cm_x          = cmd->quantum[0][0];
    cm_y          = cmd->quantum[0][1];
    sq_diameter   = cmd->quantum[1][0] * cmd->quantum[1][0];
    amplitude_in  = cmd->quantum[2][0];
    amplitude_out = cmd->quantum[2][1];
    frequency_in  = cmd->quantum[3][0];
    frequency_out = cmd->quantum[3][1];
    phase_in      = cmd->quantum[4][0];
    phase_out     = cmd->quantum[4][1];
    angle_in      = cmd->quantum[5][0];
    angle_out     = cmd->quantum[5][1];
    velocity_in   = cmd->quantum[6][0];
    velocity_out  = cmd->quantum[6][1];
    sphere        = (int)(cmd->quantum[7][0]);  /**not implemented**/

    /**Einheitsvektor in Richtung der Welle**/
    e_A_in  = cos (angle_in);
    e_B_in  = sin (angle_in);
    e_A_out = cos (angle_out);
    e_B_out = sin (angle_out);

    for (X = 0; X < A[area].d_a; X++)
        for (Y = 0; Y < A[area].d_b; Y++) {

            /**for periodic boundary**/
            diffA = (double)(X) - cm_x;
            diffB = (double)(Y) - cm_y;

            sq_dist = (double)(diffA * diffA + diffB * diffB);

            if  (sq_dist < sq_diameter) {

                directed_dist = (double)diffA * e_A_in + (double)diffB * e_B_in;

                directed_dist += velocity_in * (double)ct_t;   /**for inside**/

                cmd->S_target[ct_t][X * A[area].d_b + Y] 
                = amplitude_in
                * cos (directed_dist * 2.0 * M_PI * frequency_in + phase_in);    /**"+ phase_in"  changed recently from "- phase_in"**/

            } else {

                cmd->S_target[ct_t][X * A[area].d_b + Y] = 0.0;
            }
        }

    return (DOUBLE)(0);
}


#define ANZ_PAR 4

/******************************* gabor ***************************************/
/* Returns Gabor fkt at pos (2,3) x,y. Parameters in (1) par                 */
/* Exactly taken from gabor_analyze.c                                        */
/* r[0], r[1] = cm_x, cm_y                                                   */
/* r[2]       = height                                                       */
/* r[3]       = sigma of gaussian                                            */
/* r[4]       = frequency                                                    */
/* r[5]       = ori (between -pi/2 and pi/2)                                 */
/* r[6]       = phase                                                        */
/* XXX       = offset here assumed zero    OUT-COMMENTED!                    */

DOUBLE gabor (DOUBLE *par, double x, double y) {

  return (par[2] * exp (-0.5 * ((x-par[0])*(x-par[0]) + (y-par[1])*(y-par[1]))
                             / (par[3]*par[3])
                       )
              /* * cos (2.0 * M_PI * par[4] * ( (x-par[0]) * cos (par[5])
                                 + (y-par[1]) * sin (par[5]))
                       + par[6]
                       )*/
        /*+ par[7]*/
         );

    return (DOUBLE)(0);
}



/****************************** total_gabor **********************************/

DOUBLE total_gabor (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int X, Y;
    int area = cmd->area;
    int inarea = cmd->n_from1[0];

    if  ((cmd->anz_quant[0] != ANZ_PAR) && (A[inarea].d_n != ANZ_PAR))
        fprintf (stderr, "total_gabor: rethink how to communicate parameters!");

    for (X = 0; X < A[area].d_a; X++)
        for (Y = 0; Y < A[area].d_b; Y++) {

            /**non-periodic boundary**/
            double diffA = (double)(X) - cmd->quantum[0][0];
            double diffB = (double)(Y) - cmd->quantum[0][1];

            if  (cmd->anz_quant[0] == ANZ_PAR)
                cmd->S_target[ct_t][X * A[area].d_b + Y] = gabor (cmd->quantum[0], diffA, diffB);

            if  (A[inarea].d_n == ANZ_PAR)
                cmd->S_target[ct_t][X * A[area].d_b + Y] = gabor (cmd->S_from1[0][ct_t], (DOUBLE)X, (DOUBLE)Y);
        }

    return (DOUBLE)(0);
}
#undef ANZ_PAR

#define ANZ_PAR 6
/******************************* gauss_elliptic ******************************/
/* Returns Gauss fkt at pos (2,3) x,y. Parameters in (1) par                 */
/* p[0], p[1] = cm_x, cm_y                                                   */
/* p[2]       = height                                                       */
/* p[3], p[4] = sigma along one, and other half-axis; "a", "b"               */
/* p[5]       = orientation of half-axis "a" (not necessarily longer axis)   */

DOUBLE gauss_elliptic (DOUBLE *par, double x, double y) {

  double phi_xy = (x-par[0] == 0.0) ? 0.5*M_PI : atan ((y-par[1])/(x-par[0]));
  double phi_a  = par[5];
  double dist   = sqrt ((x-par[0])*(x-par[0]) + (y-par[1])*(y-par[1]));

  DOUBLE diff_a = dist * cos (phi_xy - phi_a);
  DOUBLE diff_b = dist * sin (phi_xy - phi_a);

  return (par[2] * exp (-0.5 * (diff_a*diff_a) / (par[3]*par[3]))
                 * exp (-0.5 * (diff_b*diff_b) / (par[4]*par[4])));

  return (DOUBLE)(0);
}

/****************************** total_gauss_elliptic *************************/
/* Takes six parameters either from q[0] or S_from1. Calls gauss_elliptic.   */

DOUBLE total_gauss_elliptic (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int X, Y;
    int area = cmd->area;
    int inarea = cmd->n_from1[0];

    if  ((cmd->anz_quant[0] != ANZ_PAR) && (A[inarea].d_n != ANZ_PAR))
        fprintf (stderr, "total_gauss_elliptic: rethink how to communicate parameters!");

    for (X = 0; X < A[area].d_a; X++)
        for (Y = 0; Y < A[area].d_b; Y++) {

            /**non-periodic boundary  double diffA = (double)(X) - cmd->quantum[0][0];
                                      double diffB = (double)(Y) - cmd->quantum[0][1]; <-- was previously used in first version for some reason ... !!!**/

            if  (cmd->anz_quant[0] == ANZ_PAR)
                cmd->S_target[ct_t][X * A[area].d_b + Y] = gauss_elliptic (cmd->quantum[0], (DOUBLE)X, (DOUBLE)Y);

            if  (A[inarea].d_n == ANZ_PAR)
                cmd->S_target[ct_t][X * A[area].d_b + Y] = gauss_elliptic (cmd->S_from1[0][ct_t], (DOUBLE)X, (DOUBLE)Y);
        }

    return (DOUBLE)(0);
}
#undef ANZ_PAR

#define ANZ_PAR 3
/******************************* my_cos **************************************/
/* p[0] = orientation (angle to the line that crosses the waves)             */
/* p[1] = frequency                                                          */
/* p[2] = phase                                                              */
/* "center" is the origin (X=0,Y=0) of the area.                             */

DOUBLE my_cos (DOUBLE *par, int X, int Y) {

    double angle  = par[0];
    double freq   = par[1];
    double phase  = par[2];

    /**Einheitsvektor in Richtung der Welle**/
    double e_A_in = cos (angle);
    double e_B_in = sin (angle);

    /**Projektion auf Einheitsvektor**/
    double directed_dist = (double)X * e_A_in + (double)Y * e_B_in;

    return (cos (directed_dist * 2.0 * M_PI * freq + phase));
}

/****************************** total_cos ************************************/
/* Takes parameters from q[0] or S_from1 (must be at least ANZ_PAR units).   */

DOUBLE total_cos (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int X, Y;
    int area = cmd->area;
    int inarea = cmd->n_from1[0];

    if  ((cmd->anz_quant[0] != ANZ_PAR) && (A[inarea].d_n < ANZ_PAR))
        fprintf (stderr, "total_cos: rethink how to communicate parameters!");

    for (X = 0; X < A[area].d_a; X++)
        for (Y = 0; Y < A[area].d_b; Y++) {

            if  (cmd->anz_quant[0] == ANZ_PAR)
                cmd->S_target[ct_t][X * A[area].d_b + Y] = my_cos (cmd->quantum[0], X, Y);

            else
                cmd->S_target[ct_t][X * A[area].d_b + Y] = my_cos (cmd->S_from1[0][ct_t], X, Y);
        }

    return (DOUBLE)(0);
}
#undef ANZ_PAR



/******************************** total_gauss_color **************************/
/* NOT YET TESTED! DESIGNED TO TRAIN A HUGE V1 WITH REPETITIVE PATTERN.      */

DOUBLE total_gauss_color (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {
     int X, Y;
     int area = cmd->area;
     double sigma = A[area].d_b / 3.0;

     /**red**/
     for (X = 0; X < A[area].d_a / 3; ++X)
         for (Y = 0; Y < A[area].d_b; ++Y) {
             double diffA = X - (double)(A[area].d_a / 3) / 2;
             double diffB = X - (double)A[area].d_b / 2;

             cmd->S_target[ct_t][X * A[area].d_b + Y] = exp (-0.5 * (diffA*diffA+diffB*diffB)/(sigma*sigma))
                                                      - exp (-0.5 * (A[area].d_b*A[area].d_b)/(sigma*sigma));
             if  (cmd->S_target[ct_t][X * A[area].d_b + Y] < 0.0)
                 cmd->S_target[ct_t][X * A[area].d_b + Y] = 0.0;
         }
     /*gnuplot  range=3.0; set xrange [-range:range]; plot exp(-0.5*x**2) - exp(-0.5*range**2)*/

     /**green -- just copy**/
     for (X = 0; X < A[area].d_a / 3; ++X)
         for (Y = 0; Y < A[area].d_b; ++Y)
             cmd->S_target[ct_t][(X + A[area].d_a / 3) * A[area].d_b + Y] = cmd->S_target[ct_t][X * A[area].d_b + Y];
     /**blue -- just copy**/       /*^^^^^^^^^^*/
     for (X = 0; X < A[area].d_a / 3; ++X)
         for (Y = 0; Y < A[area].d_b; ++Y)
             cmd->S_target[ct_t][(X + 2 * A[area].d_a / 3) * A[area].d_b + Y] = cmd->S_target[ct_t][X * A[area].d_b + Y];

    return (DOUBLE)(0);
}


/******************************** total_torus_shift_color ********************/
/* Adapted from torus_shift in data.digits.c.                                */
/* S_target must be different from S_from1!                                  */
/* NOT YET TESTED! DESIGNED TO TRAIN A HUGE V1 WITH REPETITIVE PATTERN.      */

DOUBLE total_torus_shift_color (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {
     int X, Y;
    int area = cmd->area;
     int sh_a = (int)(drand48() * A[area].d_a);
     int sh_b = (int)(drand48() * A[area].d_b);

     /**red**/
     for (X = 0; X < A[area].d_a / 3; ++X)
         for (Y = 0; Y < A[area].d_b; ++Y)
             cmd->S_target[ct_t][X * A[area].d_b + Y] = cmd->S_from1[0][ct_t][((X+sh_a)%(A[area].d_a / 3)) * A[area].d_b + (Y+sh_b)%A[area].d_b];

     /**green -- just copy**/
     for (X = 0; X < A[area].d_a / 3; ++X)
         for (Y = 0; Y < A[area].d_b; ++Y)
             cmd->S_target[ct_t][(X + A[area].d_a / 3) * A[area].d_b + Y] = cmd->S_target[ct_t][X * A[area].d_b + Y];
     /**blue -- just copy**/       /*^^^^^^^^^^*/
     for (X = 0; X < A[area].d_a / 3; ++X)
         for (Y = 0; Y < A[area].d_b; ++Y)
             cmd->S_target[ct_t][(X + 2 * A[area].d_a / 3) * A[area].d_b + Y] = cmd->S_target[ct_t][X * A[area].d_b + Y];

    return (DOUBLE)(0);
}


/******************************** total_cut_at *******************************/
/* Cut out a patch from S_from1 of size as S_target at position S_from2[0/1] */
/* S_target must be different from S_from1!                                  */
/* q[0][0] = 1: cut out w.r.t. middle of areas; else w.r.t. (0,0)-corner.    */
/* (see also total_embed, above, and total_shift_torus, below)               */

DOUBLE total_cut_at (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {
     int X, Y;
     int area    = cmd->area;
     int inarea  = cmd->n_from1[0];
     int inarea2 = cmd->n_from2[0];
     int at_a;
     int at_b;

     if  (cmd->anz_from2 != 1)
         fprintf (stderr, "\ntotal_cut_at needs a S_from2!\n");
     if  (A[inarea2].d_n < 2)
         fprintf (stderr, "\ntotal_cut_at: 2nd inarea must have 2 units for useful function!");
     if  (A[area].d_n > A[inarea].d_n)
         fprintf (stderr, "\ntotal_cut_at: target area should be smaller (or maybe equal) than input area!");

     if  (cmd->quantum[0][0] == 1) { /**take w.r.t. middles of areas**/
         at_a = (int)(cmd->S_from2[0][ct_t][0] + A[inarea].d_a / 2.0 - A[area].d_a / 2.0);
         at_b = (int)(cmd->S_from2[0][ct_t][1] + A[inarea].d_b / 2.0 - A[area].d_b / 2.0);
     } else {                        /**take w.r.t. origin (0,0) of areas**/
         at_a = (int)cmd->S_from2[0][ct_t][0];
         at_b = (int)cmd->S_from2[0][ct_t][1];
     }

     /**check whether out of boundaries**/
     if  (at_a < 0) {
         fprintf (stderr, "\ntotal_cut_at: at_a=%d, set to zero  ", at_a);
         at_a = 0;
     }
     if  (at_b < 0) {
         fprintf (stderr, "\ntotal_cut_at: at_b=%d, set to zero  ", at_b);
         at_b = 0;
     }
     if  (at_a + A[area].d_a >= A[inarea].d_a) {
         fprintf (stderr, "\ntotal_cut_at: at_a=%d, set to %d ", at_a, A[inarea].d_a - A[area].d_a - 1);
         at_a = A[inarea].d_a - A[area].d_a - 1;
     }
     if  (at_b + A[area].d_b >= A[inarea].d_b) {
         fprintf (stderr, "\ntotal_cut_at: at_b=%d, set to %d ", at_b, A[inarea].d_b - A[area].d_b - 1);
         at_b = A[inarea].d_b - A[area].d_b - 1;
     }

     for (X = 0; X < A[area].d_a; ++X)
         for (Y = 0; Y < A[area].d_b; ++Y)
             cmd->S_target[ct_t][X * A[area].d_b + Y] = cmd->S_from1[0][ct_t][(X+at_a) * A[inarea].d_b + (Y+at_b)];

    return (DOUBLE)(0);
}


/******************************** total_shift_torus **************************/
/* Shift the whole S_from1 act vector *randomly* along the 2-dim grid.       */
/* S_target must have same size as S_from1.                                  */
/* ! Target act must be different from source act (may be same area) !       */
/* Used in order to present toroidally-shifted images.                       */
/* (see also total_cut_at and total_embed, above)                            */

DOUBLE total_shift_torus (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {
     int X, Y;
     int area    = cmd->area;
     int inarea  = cmd->n_from1[0];
     int d_a     = A[area].d_a;
     int d_b     = A[area].d_b;

     if  ((d_a != A[inarea].d_a) || (d_b != A[inarea].d_b))
         fprintf (stderr, "\ntotal_shift_torus: target area must match input area!");

     if  (cmd->ch_target == cmd->ch_from1[0])
         fprintf (stderr, "\ntotal_shift_torus: target act must NOT be inarea act!");

     /**random shift**/
     int shift_a = (int)(drand48() * d_a);
     int shift_b = (int)(drand48() * d_b);

     for (X = 0; X < d_a; ++X)
     for (Y = 0; Y < d_b; ++Y) {
         int X_from = (X + shift_a) % d_a;
         int Y_from = (Y + shift_b) % d_b;
         cmd->S_target[ct_t][X * d_b + Y] = cmd->S_from1[0][ct_t][X_from * d_b + Y_from];
     }

    return (DOUBLE)(0);
}


/****************************** total_winner *********************************/
/* Simplified from total_neigh_winner. Normal mode: activate only the winner.*/
/* q[0][0] = 1: find winner; = -1: find looser                               */
/* q[1][0] = 2: Write x,y-coord of winning unit to two target area units.    */
/* Extra service: returns cmd->pointers[0] = 1 if maximum at edge of area.   */

DOUBLE total_winner (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int X, Y, X_winner = 0, Y_winner = 0;
    double max_act = 0.0;
    double min_act = 9999999.9;
    int area = cmd->area;
    int inarea = cmd->n_from1[0];

    if  ((cmd->quantum[0][0] != -1.0) && (cmd->quantum[0][0] != 1.0))
        fprintf (stderr, "\nfirst quant in total_winner wrong!\n");

    /**find winner: highest activation**/
    if  (cmd->quantum[0][0] == 1.0)
        for (X = 0; X < A[inarea].d_a; X++)
            for (Y = 0; Y < A[inarea].d_b; Y++)
                if  (cmd->S_from1[0][ct_t][Y + A[inarea].d_b * X] > max_act) {
                    max_act = cmd->S_from1[0][ct_t][Y + A[inarea].d_b * X];
                    X_winner = X;
                    Y_winner = Y;
                }

    /**find winner: lowest activation**/
    if  (cmd->quantum[0][0] == -1.0)
        for (X = 0; X < A[inarea].d_a; X++)
            for (Y = 0; Y < A[inarea].d_b; Y++)
                if  (cmd->S_from1[0][ct_t][Y + A[inarea].d_b * X] < min_act) {
                    min_act = cmd->S_from1[0][ct_t][Y + A[inarea].d_b * X];
                    X_winner = X;
                    Y_winner = Y;
                }


    if  ((cmd->anz_quantums == 2) && (cmd->quantum[1][0] == 2)) {

        /**write winner coordinates**/
        cmd->S_target[ct_t][0] = (DOUBLE)X_winner;
        cmd->S_target[ct_t][1] = (DOUBLE)Y_winner;

    } else {

        /**here activate only the winner**/
        for (X = 0; X < A[area].d_a; X++)
            for (Y = 0; Y < A[area].d_b; Y++)
                if  ((X == X_winner) && (Y == Y_winner))
                    cmd->S_target[ct_t][X * A[area].d_b + Y] = 1;
                else
                    cmd->S_target[ct_t][X * A[area].d_b + Y] = 0;
    }


    /**if winner at border, return pointer's int_val as 1**/
    if  (cmd->anz_pointers == 1) {
        if  ((X_winner == 0) || (Y_winner == 0) || (X_winner == A[inarea].d_a-1) || (Y_winner == A[inarea].d_b-1))
            cmd->pointers[0]->int_val = 1;
        else
            cmd->pointers[0]->int_val = 0;            
    }

    /**ABUSE OF FUNCTION, for square area: if winner is outside of a circle of diameter of area, then also, return pointer's int_val as 1**/
    if  (cmd->anz_pointers == 1)
        if  (A[inarea].d_a == A[inarea].d_b) {
            float half = (A[inarea].d_a - 1.0) * 0.5;
            float dist = sqrt ((X_winner-half)*(X_winner-half)+(Y_winner-half)*(Y_winner-half));
            if  (dist > half)
                cmd->pointers[0]->int_val = 1;
        }

    return (DOUBLE)(0);
}



/****************************** total_winner_per_row *************************/
/* Now with Gaussian along corresponding row around each winner per row.     */
/* q[0][0] = sigma of Gaussian                                               */

DOUBLE total_winner_per_row (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int X, Y;
    int area = cmd->area;
    int inarea = cmd->n_from1[0];
    static int firsttime = 1;
    static int *winners = NULL;

    if  (firsttime) {
        winners = (int *) malloc (A[inarea].d_a * sizeof (int));
        firsttime = 0;
    }

    /**find winner: highest activation**/
    for (X = 0; X < A[inarea].d_a; X++) {

        double max_act = 0.0;
        for (Y = 0; Y < A[inarea].d_b; Y++)
            if  (cmd->S_from1[0][ct_t][Y + A[inarea].d_b * X] > max_act) {
                max_act = cmd->S_from1[0][ct_t][Y + A[inarea].d_b * X];
                winners[X] = Y;
             }
    }

    /**activate only the winner in each row**/
    for (X = 0; X < A[area].d_a; X++)
        for (Y = 0; Y < A[area].d_b; Y++) {

            /*
            if  (Y == winners[X])
                cmd->S_target[ct_t][X * A[area].d_b + Y] = 1;
            else
                cmd->S_target[ct_t][X * A[area].d_b + Y] = 0;
            */

            double diff = (double)(Y - winners[X]);
            double sigma = cmd->quantum[0][0];
            cmd->S_target[ct_t][X * A[area].d_b + Y] = exp (-0.5 * (diff*diff) / (sigma*sigma));
        }

    return (DOUBLE)(0);
}


/****************************** total_neigh_winner ***************************/
/* Finds winner and sets Gaussian (non-periodic boundary) around. (Kohonen). */
/* q[0][0] = 1: find winner; = -1: find looser (good for min euclid distance)*/
/* q[1][0] = sigma at beginning                                              */
/* q[1][1] = sigma at first lap                                              */
/* q[1][2] = sigma from second lap on (then const until end)                 */
/* q[2][0] = time of first lap                                               */
/* q[2][1] =     "   second "                                                */
/* q[3][0]: scale of time, if function used 2ce !! TIME COUNTS INTERNALLY !! */

DOUBLE total_neigh_winner (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int X, Y, X_winner = 0, Y_winner = 0;
    double max_act = 0.0;
    double min_act = 9999999.9;
    double sigma = 0.0, timefactor = 0.0;
    static double time = 0.0;
    int area = cmd->area;

    if  (cmd->anz_quantums != 4)
        fprintf (stderr, "\n3 quants in total_neigh_winner!");
    if  ((cmd->quantum[0][0] != -1.0) && (cmd->quantum[0][0] != 1.0))
        fprintf (stderr, "\nfirst quant in total_neigh_winner wrong!\n");

    double sigma0 = cmd->quantum[1][0];
    double sigma1 = cmd->quantum[1][1];
    double sigma2 = cmd->quantum[1][2];
    double time1  = cmd->quantum[2][0];
    double time2  = cmd->quantum[2][1];

    time += cmd->quantum[3][0];

    if  (time <= time1) {
        timefactor = time / time1;
        sigma = sigma0 + timefactor * (sigma1 - sigma0);
    }
    if  ((time > time1) && (time <= time2)) {
        timefactor = (time - time1) / (time2 - time1);
        sigma = sigma1 + timefactor * (sigma2 - sigma1);
    }
    if  (time > time2)
        sigma = sigma2;

    if  ((int)time % 100 == 0)
        fprintf (stderr, "\ntotal_neigh_winner: time=%.1f sigma=%.2f  ", time, sigma);

    /**find winner: highest activation**/
    if  (cmd->quantum[0][0] == 1.0)
        for (X = 0; X < A[area].d_a; X++)
            for (Y = 0; Y < A[area].d_b; Y++)
                if  (cmd->S_from1[0][ct_t][Y + A[area].d_b * X] > max_act) {              /***was wrong (A[area].d_a instead A[area].d_b) !!!***/
                    max_act = cmd->S_from1[0][ct_t][Y + A[area].d_b * X];                 /***was wrong (A[area].d_a instead A[area].d_b) !!!***/
                    X_winner = X;
                    Y_winner = Y;
                }

    /**find winner: lowest activation**/
    if  (cmd->quantum[0][0] == -1.0)
        for (X = 0; X < A[area].d_a; X++)
            for (Y = 0; Y < A[area].d_b; Y++)
                if  (cmd->S_from1[0][ct_t][Y + A[area].d_b * X] < min_act) {
                    min_act = cmd->S_from1[0][ct_t][Y + A[area].d_b * X];
                    X_winner = X;
                    Y_winner = Y;
                }

    /**Gaussian around the winner**/
    for (X = 0; X < A[area].d_a; X++)
        for (Y = 0; Y < A[area].d_b; Y++) {

            /**non-periodic boundary**/
            double diffA = (double)(X - X_winner);
            double diffB = (double)(Y - Y_winner);

            cmd->S_target[ct_t][X * A[area].d_b + Y]
            = exp (-0.5 * (diffA*diffA + diffB*diffB) / (sigma*sigma));
        }

    return (DOUBLE)(0);
}

/****************************** total_neigh_winner_3D ************************/
/* Finds winner and sets 3D-Gaussian (non-periodic boundary) around.(Kohonen)*/
/* q[0][0] = 1: find winner; = -1: find looser (good for min euclid distance)*/
/* q[1][0] = sigma at beginning                                              */
/* q[1][1] = sigma at first lap                                              */
/* q[1][2] = sigma from second lap on (then const until end)                 */
/* q[2][0] = time of first lap                                               */
/* q[2][1] =     "   second "                                                */
/* q[3][0]: scale of time, if function used 2ce !! TIME COUNTS INTERNALLY !! */
/* q[4][0/1/2]: dimensions of area <-- only difference to total_neigh_winner */

DOUBLE total_neigh_winner_3D (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int X, Y, X_winner = 0, Y_winner = 0;
    double max_act = 0.0;
    double min_act = 9999999.9;
    double sigma = 0.0, timefactor = 0.0;
    static double time = 0.0;
    int area = cmd->area;

    if  (cmd->anz_quantums != 5)
        fprintf (stderr, "\nneed 5 quants in total_neigh_winner_3D!");
    if  ((cmd->quantum[0][0] != -1.0) && (cmd->quantum[0][0] != 1.0))
        fprintf (stderr, "\nfirst quant in total_neigh_winner wrong!\n");

    double sigma0 = cmd->quantum[1][0];
    double sigma1 = cmd->quantum[1][1];
    double sigma2 = cmd->quantum[1][2];
    double time1  = cmd->quantum[2][0];
    double time2  = cmd->quantum[2][1];

    time += cmd->quantum[3][0];

    if  (time <= time1) {
        timefactor = time / time1;
        sigma = sigma0 + timefactor * (sigma1 - sigma0);
    }
    if  ((time > time1) && (time <= time2)) {
        timefactor = (time - time1) / (time2 - time1);
        sigma = sigma1 + timefactor * (sigma2 - sigma1);
    }
    if  (time > time2)
        sigma = sigma2;

    if  ((int)time % 100 == 0)
        fprintf (stderr, "\ntotal_neigh_winner: time=%.1f sigma=%.2f  ", time, sigma);


    /**find winner: highest activation**/
    if  (cmd->quantum[0][0] == 1.0)
        for (X = 0; X < A[area].d_a; X++)
            for (Y = 0; Y < A[area].d_b; Y++)
                if  (cmd->S_from1[0][ct_t][Y + A[area].d_b * X] > max_act) {
                    max_act = cmd->S_from1[0][ct_t][Y + A[area].d_b * X];
                    X_winner = X;
                    Y_winner = Y;
                }

    /**find winner: lowest activation**/
    if  (cmd->quantum[0][0] == -1.0)
        for (X = 0; X < A[area].d_a; X++)
            for (Y = 0; Y < A[area].d_b; Y++)
                if  (cmd->S_from1[0][ct_t][Y + A[area].d_b * X] < min_act) {
                    min_act = cmd->S_from1[0][ct_t][Y + A[area].d_b * X];
                    X_winner = X;
                    Y_winner = Y;
                }

    int xdim = (int)cmd->quantum[4][0];
    int ydim = (int)cmd->quantum[4][1]; /**should be the width A[area].d_b of the area**/
    int zdim = (int)cmd->quantum[4][2]; /**should be: xdim*zdim=A[area].d_a**/

    if  (ydim != A[area].d_b) {
        fprintf (stderr, "\ntotal_neigh_winner_3D: y-dimension doesn't fit to d_b\n");
        exit (1);
    }
    if  (xdim*zdim != A[area].d_a) {
        fprintf (stderr, "\ntotal_neigh_winner_3D: x,z-dimensions don't fit to d_a\n");
        exit (1);
    }


    /**3D-Gaussian around the winner**/
    for (X = 0; X < A[area].d_a; X++)
        for (Y = 0; Y < A[area].d_b; Y++) {

            /**non-periodic boundary**/
            double diffx = (double)(X % (ydim) - X_winner % (ydim));
            double diffy = (double)(Y - Y_winner);
            double diffz = (double)((X*ydim + Y) / (xdim*ydim) - (X_winner*ydim + Y_winner) / (xdim*ydim));

            cmd->S_target[ct_t][X * A[area].d_b + Y]
            = exp (-0.5 * (diffx*diffx + diffy*diffy + diffz*diffz) / (sigma*sigma));
        }

    return (DOUBLE)(0);
}


/****************************** total_fit_gauss ******************************/
/* Finds softmax and sets Gaussian (both(!) non-periodic boundary) around it.*/
/* q[0][0]/[1]/[2] = softmax sigma at beginning / first lap / then till end  */
/* q[1][0]/[1]/[2] = replace sigma at beginning / first lap / then till end  */
/* q[2][0]/[1] = end of first / second lap, with linear progression between  */
/* q[3][0]: scale of time, if function used 2ce !! TIME COUNTS INTERNALLY !! */
/* q[4][0] < 0.0: normalize to (pos) volume, > 0.0: "sat" <- replace-gaussian*/
/* q[5][0] cut softmax-gauss beyond this range (in x and y, thus on square)  */

DOUBLE total_fit_gauss (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int i, j, X, Y, X_winner = 0, Y_winner = 0;
    double max_act = 0.0;
    double timefactor, normfactor;
    static double time = 0.0;
    static int firsttime = 1;
    static int cut = 0;
    static DOUBLE **G;
    static int ct_t_old = ct_t;
    int area = cmd->area;

    if  (cmd->anz_quantums != 6)
        fprintf (stderr, "\nwrong no of quants in total_fit_gauss");
    if  ((cmd->anz_quant[2] != 2) || (cmd->anz_quant[0] != 3) || (cmd->anz_quant[1] != 3))
        fprintf (stderr, "\n\ntotal_fit_gauss: please update parameters to new version!\n\n");

    if  (firsttime) {
        cut = (int)(cmd->quantum[5][0]);
        G = d_matrix (cut+1, cut+1);
        firsttime = 0;
    }

    /**timefactor increases from 0 to 1 at the beginning of each new relaxation**/
    if  (ct_t_old > ct_t)
        time += 1.0;
    timefactor = time / (cmd->quantum[2][0] * cmd->quantum[3][0]); /**need not be multiplied by rlen because line above**/

    double sigmasoft = 0, sigmarepl = 0;
    double sigmasoft0 = cmd->quantum[0][0];
    double sigmasoft1 = cmd->quantum[0][1];
    double sigmasoft2 = cmd->quantum[0][2];
    double sigmarepl0 = cmd->quantum[1][0];
    double sigmarepl1 = cmd->quantum[1][1];
    double sigmarepl2 = cmd->quantum[1][2];
    double time1  = cmd->quantum[2][0];
    double time2  = cmd->quantum[2][1];

    time += cmd->quantum[3][0];

    if  (time <= time1) {
        timefactor = time / time1;
        sigmasoft = sigmasoft0 + timefactor * (sigmasoft1 - sigmasoft0);
        sigmarepl = sigmarepl0 + timefactor * (sigmarepl1 - sigmarepl0);
    }
    if  ((time > time1) && (time <= time2)) {
        timefactor = (time - time1) / (time2 - time1);
        sigmasoft = sigmasoft1 + timefactor * (sigmasoft2 - sigmasoft1);
        sigmarepl = sigmarepl1 + timefactor * (sigmarepl2 - sigmarepl1);
    }
    if  (time > time2) {
        sigmasoft = sigmasoft2;
        sigmarepl = sigmarepl2;
    }

    if  ((int)time % 100 == 0)
        fprintf (stderr, "\ntotal_fit_gauss: time=%.1f sigmasoft=%.2f sigmarepl=%.2f  ", time, sigmasoft, sigmarepl);

    /**set up Gaussian for softmax (normfactor wouldn't play a role for softmax-finding!)**/
    for (i = 0; i <= cut; ++i)
        for (j = 0; j <= cut; ++j)
            G[i][j] = exp (-0.5 * (double)(i*i+j*j) / (sigmasoft * sigmasoft));

    /**find winner**/
    max_act = 0.0;
    for (X = 0; X < A[area].d_a; ++X)
        for (Y = 0; Y < A[area].d_b; ++Y) {

            double sum = 0.0;
            for (i = X-cut; i <= X+cut; ++i)
                if  ((i >= 0) && (i < A[area].d_a))
                    for (j = Y-cut; j <= Y+cut; ++j)
                        if  ((j >= 0) && (j < A[area].d_b))
                            sum += G[abs(i-X)][abs(j-Y)] * cmd->S_from1[0][ct_t][j + A[area].d_b * i];

            if  (sum > max_act) {
                max_act = sum;
                X_winner = X;
                Y_winner = Y;
            }
        }

    /**Gaussian around the winner**/
    if  (cmd->quantum[4][0] >= 0)
        normfactor = cmd->quantum[4][0];
    else
        normfactor = -cmd->quantum[4][0] / (sqrt (2.0 * M_PI) * sigmarepl);

    for (X = 0; X < A[area].d_a; X++)
        for (Y = 0; Y < A[area].d_b; Y++) {

            /**non-periodic boundary**/
            double diffA = (double)(X - X_winner);
            double diffB = (double)(Y - Y_winner);

            cmd->S_target[ct_t][X * A[area].d_b + Y]
            = normfactor * exp (-0.5 * (diffA*diffA + diffB*diffB) / (sigmarepl*sigmarepl));
        }

    return (DOUBLE)(0);
}


/****************************** total_gauss_at *******************************/
/* Sets Gaussian (non-periodic boundary) on a 1-, 2- or "3"-dim area.        */
/* Place where to set is given in S_from1. Dimension given in anz_quant[0].  */
/* q[0][0,1,2] = sigmas                                                      */
/* q[1][0,1,2] = resolutions; where q11*q12 must be d_b if "3"-dim           */
/* q[2][0,1,2] = lower bounds of inputs                                      */
/* q[3][0,1,2] = upper bounds of inputs                                      */

DOUBLE total_gauss_at (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int i, X, Y;
    int area = cmd->area;
    double sigma[3];
    double res[3];
    double cm[3];
    double d[3];

    if  (cmd->anz_quantums != 4)
        fprintf (stderr, "\nwrong no of quantums in total_gauss_at");

    int dim = cmd->anz_quant[0];

    if  ((dim < 1) || (dim > 3))
        fprintf (stderr, "\nwrong dim in total_gauss_at");

    for (i = 0; i < dim; ++i) {
        sigma[i]     = cmd->quantum[0][i];
        res[i]       = cmd->quantum[1][i];
        double lower = cmd->quantum[2][i];
        double upper = cmd->quantum[3][i];
        cm[i]  = (cmd->S_from1[0][ct_t][i] - lower) / (upper - lower) * res[i];
    }

    if  (dim == 1) {
        for (Y = 0; Y < A[area].d_a * A[area].d_b; Y++) {

            d[0] = (double)(Y - cm[0]);
            cmd->S_target[ct_t][Y] = exp (-0.5 * (d[0]*d[0]) / (sigma[0]*sigma[0]));
        }
    }

    if  (dim == 2) {
        for (X = 0; X < A[area].d_a; X++)
        for (Y = 0; Y < A[area].d_b; Y++) {

            d[0] = (double)(X - cm[0]);
            d[1] = (double)(Y - cm[1]);

            cmd->S_target[ct_t][X * A[area].d_b + Y]
            = exp (-0.5 * (d[0]*d[0]) / (sigma[0]*sigma[0]))
            * exp (-0.5 * (d[1]*d[1]) / (sigma[1]*sigma[1]));
        }
    }

    if  (dim == 3) {

        int xdim = (int)res[0];
        int ydim = (int)res[1]; /**should be the width A[area].d_b of the area**/
        int zdim = (int)res[2]; /**should be: xdim*zdim=A[area].d_a**/

        if  (ydim != A[area].d_b) {
            fprintf (stderr, "\ntotal_gauss_at: y-dimension doesn't fit to d_b\n");
            exit (1);
        }
        if  (xdim*zdim != A[area].d_a) {
            fprintf (stderr, "\ntotal_gauss_at: x,z-dimensions don't fit to d_a\n");
            exit (1);
        }

        /**3D-Gaussian**/
        for (X = 0; X < A[area].d_a; X++)
        for (Y = 0; Y < A[area].d_b; Y++) {

            d[0] = (double)(X % (ydim) - cm[0]);
            d[1] = (double)(Y - cm[1]);
            d[2] = (double)((X*ydim + Y) / (xdim*ydim) - cm[2]);

            cmd->S_target[ct_t][X * A[area].d_b + Y]
            = exp (-0.5 * (d[0]*d[0]) / (sigma[0]*sigma[0]))
            * exp (-0.5 * (d[1]*d[1]) / (sigma[1]*sigma[1]))
            * exp (-0.5 * (d[2]*d[2]) / (sigma[2]*sigma[2]));
        }

    }

    return (DOUBLE)(0);
}


/****************************** total_softmax ********************************/
/* Finds softmax and writes x,y position to 1st two units of another area.   */
/* q[0][0] = sigma of Gaussian                                               */
/* q[1][0] cut softmax-gauss beyond this range (in x and y, thus on square)  */

DOUBLE total_softmax (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int i, j, X, Y, X_winner = 0, Y_winner = 0;
    double max_act = 0.0;
    double sigma;
    static int firsttime = 1;
    static int cut = 0;
    static DOUBLE **G;

    int inarea = cmd->n_from1[0];

    if  (cmd->anz_quantums != 2)
        fprintf (stderr, "\nwrong no of quants in total_softmax");

    if  (firsttime) {
        cut = (int)(cmd->quantum[1][0]);
        G = d_matrix (cut+1, cut+1);
        firsttime = 0;
    }

    sigma = cmd->quantum[0][0];

    /**set up Gaussian for softmax (normfactor wouldn't play a role for softmax-finding!)**/
    for (i = 0; i <= cut; ++i)
        for (j = 0; j <= cut; ++j)
            G[i][j] = exp (-0.5 * (double)(i*i+j*j) / (sigma*sigma));

    /**find winner**/
    max_act = 0.0;
    for (X = 0; X < A[inarea].d_a; ++X)
        for (Y = 0; Y < A[inarea].d_b; ++Y) {

            double sum = 0.0;
            for (i = X-cut; i <= X+cut; ++i)
                if  ((i >= 0) && (i < A[inarea].d_a))
                    for (j = Y-cut; j <= Y+cut; ++j)
                        if  ((j >= 0) && (j < A[inarea].d_b)) {
                            int dist_i = (i-X >= 0) ? i-X : X-i;
                            int dist_j = (j-Y >= 0) ? j-Y : Y-j;
                            sum += G[dist_i][dist_j] * cmd->S_from1[0][ct_t][j + A[inarea].d_b * i];
                        }

            if  (sum > max_act) {
                max_act = sum;
                X_winner = X;
                Y_winner = Y;
            }
        }

    cmd->S_target[ct_t][0] = X_winner;
    cmd->S_target[ct_t][1] = Y_winner;

    return (DOUBLE)(0);
}

/****************************** total_error_at *******************************/
/* At the x,y-position given by the values of S_from2[0][t][0,1],            */
/* write the value S_from2[0][t][x/y] - S_from2[1][t][x/y]. Where x^=0 y^=1. */ /*ATTENTION! NO DIFFERENCE TAKEN; INCOMPATIBLE WITH OLDER PROG!*/
/* At all other positions write 0.                                           */
/* Used in conjunction with total_softmax in order to show the deviation of  */
/* two softmax hills w.r.t. two different activity patterns.                 */
/* i.e.: where one softmax is placed, how far away is the other from there.  */
/* q[0][0] = 0: write x-error; = 1: write y-error; =2: counter to normalise. */

DOUBLE total_error_at (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int d_a = A[cmd->area].d_a;
    int d_b = A[cmd->area].d_b;

    for (int X = 0; X < d_a; X++)
        for (int Y = 0; Y < d_b; Y++)
            cmd->S_target[ct_t][X * d_b + Y] = 0.0;

    int target_pos_x = (int)(cmd->S_from2[0][ct_t][0]); /**both S_from2 input areas have two units denoting x,y pos**/
    int target_pos_y = (int)(cmd->S_from2[0][ct_t][1]);
    int attrac_pos_x = (int)(cmd->S_from2[1][ct_t][0]);
    int attrac_pos_y = (int)(cmd->S_from2[1][ct_t][1]);

    /**error in x-direction**/
    if  (cmd->quantum[0][0] == 0)
        cmd->S_target[ct_t][target_pos_x * d_b + target_pos_y] = /*** !!! target_pos_x - !!! ***/   attrac_pos_x;

    /**error in y-direction**/
    if  (cmd->quantum[0][0] == 1)
        cmd->S_target[ct_t][target_pos_x * d_b + target_pos_y] = /*** !!! target_pos_y - !!! ***/   attrac_pos_y;

    /**count if this point is a target**/
    if  (cmd->quantum[0][0] == 2)
        cmd->S_target[ct_t][target_pos_x * d_b + target_pos_y] = 1;

    return (DOUBLE)0;
}

/****************************** gnuplot_error_at *****************************/
/* To be used after total_softmax and total_error_at. Inputs (on each unit): */
/* S_from1[0]/[1]: x/y-deviations, S_from1[2]: count for normalisation.      */
/* Prints out commands to draw a gnuplot flowfield via cut&paste.            */
/* Can be used more generally, with just x/y inputs (without normalisation). */
/* Can be used in "total" or "alltime" sweep, where it then uses only "begin"*/
/* q[0][0] = scaling factor for the difference                               */

DOUBLE gnuplot_error_at (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    if  (cmd->anz_from1 < 2) {
        fprintf (stderr, "\n\nwrong use of gnuplot_error_at\n\n");
        exit (1);
    }

    int d_a = A[cmd->area].d_a;
    int d_b = A[cmd->area].d_b;

    float scale = cmd->quantum[0][0];

    char fullname[512];
    FILE *fp;

    sprintf (fullname, "%s/gnuplot.dat", cmd->pointers[0]->words[0]);                             /**first pointer gives the directory!**/
    if  ((fp = fopen (fullname, "w")) == 0)
        fprintf (stderr, "\n\n\nError opening file %s/gnuplot.dat", cmd->pointers[0]->words[0]);

    fprintf (fp, "set xrange [-0.5:%d.5]\n", d_a-1);
    fprintf (fp, "set yrange [-0.5:%d.5]\n", d_b-1);

    for (int X = 0; X < d_a; X++)
        for (int Y = 0; Y < d_b; Y++) {

            float X_dev = cmd->S_from1[0][ct_t][X * d_b + Y];
            float Y_dev = cmd->S_from1[1][ct_t][X * d_b + Y];
            float norm = 1.0;

            if  (cmd->anz_from1 == 3)
                norm = cmd->S_from1[2][ct_t][X * d_b + Y];
            if  (norm == 0.0)
                norm = 1.0;

            fprintf (fp, "set arrow from %d,%d to %f,%f\n", X, Y, X + scale * X_dev / norm, Y + scale * Y_dev / norm);
        }

    fprintf (fp, "set size 0.721,1.0\n");
    fprintf (fp, "set nokey\n");
    fprintf (fp, "plot -1\n");

    fclose (fp);

    fprintf (stderr, "\ngnuplot\n");
    fprintf (stderr, "load \"%s/gnuplot.dat\"\n", cmd->pointers[0]->words[0]);

    return (DOUBLE)0;
}




/****************************** total_dist_xy_middle *************************/
/* Returns the x,y-deviation (with sign) of the max act neuron to the middle.*/
/* q[0][0]/[1] scaling factor, vertical/horizontal along d_a/d_b             */

DOUBLE total_dist_xy_middle (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int X, Y, X_winner = 0, Y_winner = 0;
    double max_act = 0.0;
    int d_a, d_b;
    int inarea = cmd->n_from1[0];

    if  (cmd->anz_quantums != 1)
        fprintf (stderr, "\ntotal_dist_middle needs only q[0][0/1] - are you sure the input area sizes need to be given in addition? ");

    if  (cmd->anz_quant[0] != 2)
        fprintf (stderr, "\nwrong no of quant in total_dist_middle");

    if  (A[cmd->area].d_n != 2) {
        fprintf (stderr, "\ntotal_dist_middle must give values to 2 neurons!");
        exit (1);
    }

    d_a = A[inarea].d_a;
    d_b = A[inarea].d_b;

    /**find winner**/
    for (X = 0; X < d_a; X++)
        for (Y = 0; Y < d_b; Y++)
            if  (cmd->S_from1[0][ct_t][Y + d_b * X] > max_act) {
                max_act = cmd->S_from1[0][ct_t][Y + d_b * X];
                X_winner = X;
                Y_winner = Y;
            }

/*fprintf (stderr, "\ndist %c: X_winner=%d dist=%.2f   Y_winner=%d dist=%.2f  ", cmd->b_target, X_winner, d_a*0.5 - X_winner, Y_winner, d_b*0.5 - Y_winner);*/

    cmd->S_target[ct_t][0] = ((double)(d_a)*0.5 - (double)(X_winner))
                           * cmd->quantum[0][0];
    cmd->S_target[ct_t][1] = ((double)(d_b)*0.5 - (double)(Y_winner))
                           * cmd->quantum[0][1];

    return (DOUBLE)(0);
}




/****************************** total_rand ***********************************/
/* Like local rand but all neurons get same activation.                      */

DOUBLE total_rand (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int X, Y;
    DOUBLE val, *par;
    int area = cmd->area;

    if  (cmd->anz_quant[0] != 2)
        fprintf (stderr, "\nwrong no of quant[0] in total_rand");

    par = cmd->quantum[0];

    val = par[0] + (par[1] - par[0]) * drand48();

    for (X = 0; X < A[area].d_a; X++)
        for (Y = 0; Y < A[area].d_b; Y++)
            cmd->S_target[ct_t][X * A[area].d_b + Y] = val;

    return (DOUBLE)(0);
}


/******************************* total_as_rand *******************************/
/* Like local_as_rand, but continues until at least one neuron is ON.        */
/* par[0] scales input                                                       */
/* par[1] scales output (like sat)                                           */
/* q[0][2]=2: per row, i.e. make sure that one in each row is ON             */
/* Requires different output vector than input!                              */

DOUBLE total_as_rand (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int area = cmd->area;
    int counter = 0;
    int is_set = 0;

    if  (cmd->ch_target == cmd->ch_from1[0])
        fprintf (stderr, "\n\n total_as_rand needs other target vector than source vector!\n");

    if  ((cmd->anz_quant[0] > 2) && (cmd->quantum[0][2]) == 2) {

        /**make sure at least one unit _per_row_ is ON**/

        for (int X = 0; X < A[area].d_a; X++) {          /**for each row**/
            is_set = 0;
            counter = 0;
            do  {
                for (int Y = 0; Y < A[area].d_b; Y++)
                    if  (drand48() < (cmd->quantum[0][0] * cmd->S_from1[0][ct_t][X * A[area].d_b + Y])) {
                        cmd->S_target[ct_t][X * A[area].d_b + Y] = cmd->quantum[0][1];
                        is_set += 1;
                    } else {
                        cmd->S_target[ct_t][X * A[area].d_b + Y] = 0.0;
                    }
                counter += 1;
            } while ((is_set == 0) && (counter < 1000));    /**until a unit is ON**/

           /**swich random ones OFF _per_row_ until only one unit is ON**/
           int on;
           for (on = is_set; on > 1; ) {
               int off = X * A[area].d_b + (int)(drand48() * A[area].d_b);
               if  (cmd->S_target[ct_t][off] == cmd->quantum[0][1]) {
                   cmd->S_target[ct_t][off] = 0.0;
                   on -= 1;
               }
           }
           
        }

    } else {

        /**make sure at least one unit is ON**/

        is_set = 0;
        counter = 0;
        do  {
            for (int X = 0; X < A[area].d_a; X++)
                for (int Y = 0; Y < A[area].d_b; Y++)
                    if  (drand48() < (cmd->quantum[0][0] * cmd->S_from1[0][ct_t][X * A[area].d_b + Y])) {
                        cmd->S_target[ct_t][X * A[area].d_b + Y] = cmd->quantum[0][1];
                        is_set += 1;
                    } else {
                        cmd->S_target[ct_t][X * A[area].d_b + Y] = 0.0;
                    }
            counter += 1;
        } while ((is_set == 0) && (counter < 1000));
    }

    /**swich random ones OFF until only one unit is ON
    int on;
    for (on = is_set; on > 1; ) {
        int off = (int)(drand48() * A[area].d_n);
        if  (cmd->S_target[ct_t][off] == cmd->quantum[0][1]) {
            cmd->S_target[ct_t][off] = 0.0;
            on -= 1;
        }
    }
    **/

    if  (counter > 500)
        fprintf (stderr, "\n total_as_rand tried %d times and got %d result\n", counter, is_set);

    return (DOUBLE)(0);
}


/****************************** total_any_nonzero ****************************/
/* Returns a cmd->pointer where int_val is 1 if any input unit is active.    */

DOUBLE total_any_nonzero (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int nonzero = 0;
    int inarea = cmd->n_from1[0];

    for (int i = 0; i < A[inarea].d_n; ++i)
        if  (cmd->S_from1[0][ct_t][i] != 0.0)
            nonzero = 1;

    if  (nonzero)
        cmd->pointers[0]->int_val = 1;
    else
        cmd->pointers[0]->int_val = 0;

    return (DOUBLE)(0);
}




/****************************** total_pause **********************************/
/* wait / sleep / halt / continue after time cmd->quantum[0][0] milliseconds.*/

DOUBLE total_pause (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

  //time_t now = time (NULL);
  //do  { *** nothing *** }  while (time (NULL) < now + cmd->quantum[0][0]);

  usleep ((long)(1000 * cmd->quantum[0][0]));

    return (DOUBLE)(0);
}





/******************************************************************************
MirrorBot      Phonemes                         Phonemes
Word           Sphinx/Festival                  CEDEX
======= =======================                 ===========

BACKWARD B AE K W AXR D                          b & k w @ d
BALL     B AO L                                  b O: l
BLACK    B L AE K                                b l & k
BLUE     B L UW                                  b l u:
BODY     B AA DX IY                              b O d I
BOT      B AA T                                  b O t
BROWN    B R AW N                                b r a U n
CAT      K AE T                                  k & t
CUP      K AH P                                  k V p
DESK     D EH S K                                d E s k
DOG      D AO G                                  d O g
DOWN     D AW N                                  d a U n
DROP     D R AA P                                d r O p
DROP(2)  D R AO P                                d r O p
FORWARD  F AO R W AXR D                          f O: w @ d
GO       G OW                                    g @ U
HEAD     HH EH D                                 h E d
LEFT     L EH F T                                l E f t
LIFT     L IH F T                                l I f t
MOVE     M UW V                                  m u: v
NUT      N AH T                                  n V t
PICK     P IH K                                  p I k
PLUM     P L AH M                                p l V m
PUT      P UH T                                  p U t
RIGHT    R AY T                                  r a I t
SAM      S AE M                                  s & m
SHOW     SH OW                                   S @ U
STOP     S T AA P                                s t O p
TOUCH    T AH CH                                 t V tS
TURN     T ER N                                  t 3: n
UP       AH P                                    V p
WALL     W AO L                                  w O: l
WHITE    HH W AY T                               h w a I t
WHITE(2) W AY T                                  w a I t

******************************************************************************/


#define word_length 4       /**4 or 6   undefined after total_behaviour !!! **/


/****************************** phoneme2array ********************************/
/* not a simulator function                                                  */
/* used in total behaviour, for language input as Andreas Knoblauch does     */

DOUBLE phoneme2array (char *first, char *second, char *third,
                    char *fourth, char *fifth, char *sixth, DOUBLE *target) {

#define num_phonemes 39

    char *phonemes[num_phonemes] = {"AA", "AE", "AH", "AO", "AW", "AY",
    "B", "CH", "D", "DH", "EH", "ER", "EY", "F", "G", "HH", "IH", "IY",
    "JH", "K", "L", "M", "N", "NG", "OW", "OY", "P", "R", "S", "SH",
    "T", "TH", "UH", "UW", "V", "W", "Y", "Z", "ZH"};

    int j;

    for (j = 0; j < num_phonemes; ++j) {

        if  (! strcmp (first, phonemes[j]))
            target[j] = 1.0;
        else
            target[j] = 0.0;

        if  (! strcmp (second, phonemes[j]))
            target[j + num_phonemes] = 1.0;
        else
            target[j + num_phonemes] = 0.0;

        if  (! strcmp (third, phonemes[j]))
            target[j + 2 * num_phonemes] = 1.0;
        else
            target[j + 2 * num_phonemes] = 0.0;

        if  (! strcmp (fourth, phonemes[j]))
            target[j + 3 * num_phonemes] = 1.0;
        else
            target[j + 3 * num_phonemes] = 0.0;

      if  (word_length > 4) {

        if  (! strcmp (fifth, phonemes[j]))
            target[j + 4 * num_phonemes] = 1.0;
        else
            target[j + 4 * num_phonemes] = 0.0;

        if  (! strcmp (sixth, phonemes[j]))
            target[j + 5 * num_phonemes] = 1.0;
        else
            target[j + 5 * num_phonemes] = 0.0;
      }
    }

#undef num_phonemes

    return (DOUBLE)(0);
}

/******************************* CEDEX2array *********************************/
/* not a simulator function                                                  */
/* used in total behaviour, for language input as Fermin / Andreas does      */

DOUBLE CEDEX2array (char *first, char *second, char *third,
                    char *fourth, char *fifth, char *sixth, DOUBLE *target) {

#define num_phonemes 46
#define phoneme_dim 20

    char *phonemes[num_phonemes] = {"p", "b", "t", "d", "k", "g", "N", "m",
    "n", "l", "r", "f", "v", "T", "D", "s", "z", "S", "Z", "j", "x", "h", "w",
    "tS", "dZ", "N,", "m,", "n,", "l,", "r*", "I", "E", "&", "V", "O", "U",
    "@", "i", "A", "u", "3", "e", "a", "&~", "O~", "A~"};

    int ivector[num_phonemes][phoneme_dim] = {
        {1,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0},
        {1,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0},
        {1,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0},
        {1,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0},
        {1,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0},
        {1,0,0,1,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0},
        {1,1,0,1,0,0,0,1,0,0,0,0,0,1,0,1,0,0,0,0},
        {1,1,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0},
        {1,1,0,1,0,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0},
        {1,1,1,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0},
        {1,1,1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0},
        {1,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {1,0,0,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
        {1,0,0,0,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0},
        {1,0,0,1,0,0,1,0,0,0,1,1,0,0,0,0,0,0,0,0},
        {1,0,0,0,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,0},
        {1,0,0,1,0,0,1,0,0,0,1,0,1,0,0,0,0,0,0,0},
        {1,0,0,0,0,0,1,0,0,0,0,1,1,1,0,0,0,0,0,0},
        {1,0,0,1,0,0,1,0,0,0,0,1,1,1,0,0,0,0,0,0},
        {0,1,1,1,0,0,1,0,0,0,1,1,0,1,0,0,0,0,0,0},
        {1,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,0,0,0},
        {0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,1,0,0,0,0},
        {0,1,1,1,0,0,1,0,0,1,0,0,0,1,0,0,0,0,0,0},
        {1,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0},
        {1,0,0,0,0,0,0,0,0,0,1,1,1,1,0,0,0,0,0,0},
        {1,1,0,1,0,0,0,1,0,0,0,0,0,1,0,0,0,0,1,0},
        {1,1,0,1,0,0,0,1,0,0,0,1,0,0,0,0,0,0,1,0},
        {1,1,0,1,0,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0},
        {1,1,1,1,0,0,0,0,1,0,1,0,0,0,0,0,0,0,1,0},
        {1,1,1,1,0,0,1,0,0,0,1,0,0,0,0,0,0,0,0,0},
        {0,1,1,1,0,0,1,0,0,0,0,0,0,1,0,0,0,0,1,1},
        {0,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1},
        {0,1,1,1,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,1},
        {0,1,1,1,0,0,1,0,0,0,0,0,0,0,1,0,0,0,1,0},
        {0,1,1,1,0,0,1,0,0,1,0,0,0,0,1,1,0,0,1,0},
        {0,1,1,1,0,0,1,0,0,1,0,0,0,1,0,1,0,0,1,0},
        {0,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,0},
        {0,1,1,1,0,0,1,0,0,0,0,0,0,1,0,0,1,0,1,1},
        {0,1,1,1,0,0,1,0,0,0,0,0,0,0,1,1,1,0,1,0},
        {0,1,1,1,0,0,1,0,0,1,0,0,0,1,0,1,1,0,1,0},
        {0,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,1,0,1,0},
        {0,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,1,0,1,1},
        {0,1,1,1,0,0,1,0,0,0,0,0,0,0,1,0,1,0,1,0},
        {0,1,1,1,0,0,1,1,0,0,0,0,0,0,1,0,0,0,1,1},
        {0,1,1,1,0,0,1,1,0,1,0,0,0,0,1,1,0,0,1,0},
        {0,1,1,1,0,0,1,1,0,0,0,0,0,0,1,1,1,0,1,0}};

    int i, j, k;

    /**init with zero**/
    for (i = 0; i < word_length * phoneme_dim; ++i)
        target[i] = 0.0;

    for (j = 0; j < num_phonemes; ++j) {

        if  (! strcmp (first, phonemes[j]))
            for (k = 0; k < phoneme_dim; ++k)
                target[0 * phoneme_dim + k] = ivector[j][k];

        if  (! strcmp (second, phonemes[j]))
            for (k = 0; k < phoneme_dim; ++k)
                target[1 * phoneme_dim + k] = ivector[j][k];

        if  (! strcmp (third, phonemes[j]))
            for (k = 0; k < phoneme_dim; ++k)
                target[2 * phoneme_dim + k] = ivector[j][k];

        if  (! strcmp (fourth, phonemes[j]))
            for (k = 0; k < phoneme_dim; ++k)
                target[3 * phoneme_dim + k] = ivector[j][k];

      if  (word_length > 4) {

        if  (! strcmp (fifth, phonemes[j]))
            for (k = 0; k < phoneme_dim; ++k)
                target[4 * phoneme_dim + k] = ivector[j][k];

        if  (! strcmp (sixth, phonemes[j]))
            for (k = 0; k < phoneme_dim; ++k)
                target[5 * phoneme_dim + k] = ivector[j][k];
      }
    }

#undef num_phonemes
#undef phoneme_dim

    return (DOUBLE)(0);
}



/****************************** settargetangle *******************************/
/* not a simulator function                                                  */
/* used in total behaviour, for wander/carry to bounce off walls             */

double settargetangle (double a, double b) {
    double randval = drand48();
    double prob_to_turn_randomly = 0.0; /**see also below!**/

    if  (randval < prob_to_turn_randomly)
        return (a + drand48() * (b - a));
    else
        return b;

    return (DOUBLE)(0);
}



/****************************** total_behaviour ******************************/

DOUBLE total_behaviour (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {
    int area = cmd->area;

   double x          = cmd->S_from1[0][ct_t][0];
   double y          = cmd->S_from1[0][ct_t][1];
   double phi        = cmd->S_from1[0][ct_t][2];
   double gripopen   = cmd->S_from1[0][ct_t][3];
   double gripheight = cmd->S_from1[0][ct_t][4];

   double proximity;

   const double pi = 3.1415926535;
   double newphi = (180/pi) * phi;
   const double epsangle = 0.11 * 180 / pi;

   if  (fabs(newphi) > 180)
       fprintf (stderr, "\nwarning: phi has unexpected values!\n");

   double d     = sqrt (x*x + y*y);   /**distance to target (also == a*a + b*b)**/
   double theta = atan (y/x) - phi;   /**seen angle to target in visual field**/

   double a = (double)(16) - d * cos (theta);
   double b = (double)(24) / 2.0 + d * sin (theta);

   int motor     = cmd->quantum[0][0] == 1 ? 1 : 0;
   int vision    = cmd->quantum[0][0] == 2 ? 1 : 0;
   int gripper   = cmd->quantum[0][0] == 3 ? 1 : 0;
   int language  = cmd->quantum[0][0] == 4 ? 1 : 0;
   int reinforce = cmd->quantum[0][0] == 5 ? 1 : 0;

   static int counter = 0;

   enum Behaviour {wander, docking, grasp, carry, put, release};

   static Behaviour behaviour           = wander;
   static Behaviour behavioursuggestion = wander;

/* non-static used for performance test in 32.EpiRob04_testperformance.enr:
   Behaviour behaviour           = docking;
   Behaviour behavioursuggestion = docking;
*/
   const double avoidancedistance = 5.0;
   const double ywall = 28 - 1; /*16;*/  /**if values change then**/
   const double xwall = 56 + 10; /*48;*/  /**check also above and below!**/

   /**new: second argument restricts the selection of allowed behaviours**/
   int num_allowed_behaviours = (cmd->anz_quantums == 1) ? 6 : cmd->anz_quant[1];
   int allowed_behaviours[6];
   for (int i = 0; i < num_allowed_behaviours; ++i) {
       if  (num_allowed_behaviours < 6)
           allowed_behaviours[i] = (int)(cmd->quantum[1][i]);
       else
           allowed_behaviours[i] = i;
   }


   /**set output to zero as default**/
   for (int ct_n = 0; ct_n < A[area].d_n; ++ct_n)
       cmd->S_target[ct_t][ct_n] = 0.0;


   /**tests**/
   if  (gripper)
       if  ((A[area].d_n != 5) && (A[area].d_n != 0))   /**correct target area size?**/
           fprintf (stderr, "\ntotal_behaviour: gripper area has incorrect size!\n");

   /**tests**/
   if  (motor) {
       /**two areas as input? (second is motor for copying)**/
       if  (cmd->anz_from1 != 2)
           fprintf (stderr, "\ntotal_behaviour: need 2 anz_from1!\n");
       /**second input area == target area?**/
       if  (cmd->n_from1[1] != cmd->area)
           fprintf (stderr, "\ntotal_behaviour: different areas! (wrong here)\n");
       /**size of current (target) area == size of (second) input area?**/
       if  ((A[area].d_n) != (A[cmd->n_from1[1]].d_n))
           fprintf (stderr, "\ntotal_behaviour: incompatible area sizes!\n");
   }


   if  (behaviour == docking) {    /**coordinate origin visible AND angle in interval -45...45 deg <-- consider somehow !!!**/

       int timeout = 200;
       static int total_time = 0;
       total_time += 1;

       if  (motor) {
           if  (total_time < timeout)
               for (int ct_n = 0; ct_n < A[area].d_n; ++ct_n)
                   cmd->S_target[ct_t][ct_n] = cmd->S_from1[1][ct_t][ct_n];   /**overtake docking wheel motor command**/
           else
               cmd->S_target[ct_t][3] = 1.0;   /**turn for a while if timeout**/
       }

       if  (gripper) {
           if  (gripopen == 0.0)
               cmd->S_target[ct_t][1] = 1.0;   /**open gripper if closed**/

           if (gripheight > 0.5)
               cmd->S_target[ct_t][3] = 1.0;   /**lower gripper to table height**/
       }

       if  (language) {
           if  (A[area].d_n == word_length * 39)
               phoneme2array ("P", "IH", "K", "-", "-", "-", cmd->S_target[ct_t]);
           else
           if  (A[area].d_n == word_length * 20)
               CEDEX2array ("p", "I", "k", "-", "-", "-", cmd->S_target[ct_t]);
           else
               for (int i = 0; i < A[area].d_b; ++i)
                   cmd->S_target[ct_t][1 * A[area].d_b + i] = 1.0;     /**behaviour report to language**/
       }                                                          /**low level action reports are same for all behaviours -- see below**/

       /**transitions**/
       if  ((x < 3) && (fabs(y) < 4)) {     /**at fruit, analogeously to break-beam**/

           if  (gripper)
               cmd->S_target[ct_t][5] = 1.0;       /**break-beam sensor, just at last time step**/

           behavioursuggestion = grasp;
           total_time = 0;
       }
       if  (x < 1) {                          /**avoid robot catastrophe**/
           behavioursuggestion = grasp;
           total_time = 0;
       }
       if  (a < 0 || b < 0 || a >= 16 || b >= 24) {   /**orange is not visible**/
           behavioursuggestion = wander;
           total_time = 0;
       }
       if  (total_time > timeout + 20) {                     /**if timeout**/
           behavioursuggestion = wander;
           total_time = 0;

           fprintf (stderr, "\ntimeout in docking; x=%f, y=%f, phi=%f   ", x, y, phi);
       }
   }

   else if  (behaviour == grasp) {

       if  (gripper) {

           cmd->S_target[ct_t][0] = 1.0;       /**close gripper**/

           if  (gripheight < 1.0)
               cmd->S_target[ct_t][2] = 1.0;   /**lift gripper**/

           cmd->S_target[ct_t][5] = 1.0;       /**break-beam sensor**/
       }

       if  (motor) {
           static int dotime = 0;
           static int dotimemax = 20 + (int)(drand48() * 10.0); /**changed to 10 from 20 !!!**/

           double prob_to_turn_randomly = 0.0;   /**see also above!**/

           //go backward 30cm (10 pixels)
        /* if  (dotime < 16) { */

           if  (x < 14) {                  /**new: start turning at a certain distance, not after certain time !!!**/

               cmd->S_target[ct_t][1] = 1.0;                   /**go backward**/
           } else {
               if  (newphi > 0)
                   cmd->S_target[ct_t][2] = 1.0;   /**turn right**/
               else
                   cmd->S_target[ct_t][3] = 1.0;   /**turn left**/
           }

           dotime += 1;

           if  (dotime >= dotimemax) {
               behavioursuggestion = carry;          /**transition to carrybehaviour**/

               if  (num_allowed_behaviours == 3)
                   behavioursuggestion = wander;         /*** !!! THIS OVERRIDES THE SIX BEHAVIOURS SO THAT ONLY THREE ARE USED !!! ***/

               dotime = 0;
               if  (drand48() < prob_to_turn_randomly)
                   dotimemax = 20 + (int)(drand48() * 20.0);
               else
                   dotimemax = 40;
           }
       }

       if  (language) {
           if  (A[area].d_n == word_length * 39)
               phoneme2array ("L", "IH", "F", "T", "-", "-", cmd->S_target[ct_t]);
           else
           if  (A[area].d_n == word_length * 20)
               CEDEX2array ("l", "I", "f", "t", "-", "-", cmd->S_target[ct_t]);
           else
               for (int i = 0; i < A[area].d_b; ++i)
                   cmd->S_target[ct_t][2 * A[area].d_b + i] = 1.0;     /**behaviour report to language**/
       }
   }

   else if  (behaviour == put) {
       //(coordinate origin visible AND angle in zone)

       static int total_time = 0;
       total_time += 1;

       if  (motor) {
           if  (total_time < 1000)
               for (int ct_n = 0; ct_n < A[area].d_n; ++ct_n)
                   cmd->S_target[ct_t][ct_n] = cmd->S_from1[1][ct_t][ct_n];   /**overtake docking wheel motor command**/
           else
               cmd->S_target[ct_t][3] = 1.0;   /**turn for a while if timeout**/
       }

       if  (gripper) {

           if  (gripheight > 0.5)
               cmd->S_target[ct_t][3] = 1.0;   /**lower gripper to table height**/

           cmd->S_target[ct_t][5] = 1.0;       /**break-beam sensor**/
       }

       if  (language)
           for (int i = 0; i < A[area].d_b; ++i)
               cmd->S_target[ct_t][4 * A[area].d_b + i] = 1.0;      /**behaviour report to language**/


       /**grasping at coordinate origin**/
       if  ((x < 3) && (fabs(y) < 4)) {    /*at fruit, analogeously to break-beam*/
           behavioursuggestion = release;
       }
       if  (x < 1) {                        /**avoid robot catastrophe**/
           behavioursuggestion = release;
       }
       if  (a < 0 || b < 0 || a >= 16 || b >= 24) {   /**target is not visible**/
           behavioursuggestion = carry;
       }
       if  (total_time > 1020) {                     /**if timeout**/
           behavioursuggestion = carry;
           total_time = 0;

           fprintf (stderr, "\ntimeout in put; x=%f, y=%f, phi=%f   ", x, y, phi);
       }
    }


    else if (behaviour == release) {

       if (gripper) {

           cmd->S_target[ct_t][1] = 1.0;   /**open gripper**/

           if  (gripheight < 1.0)
               cmd->S_target[ct_t][2] = 1.0;   /**lift gripper**/
       }

       if  (motor) {
	   static int dotime = 0;
	   static int dotimemax = 20 + (int)(drand48() * 40.0);

           //go backward 30cm (10 pixels)
           if  (dotime < 16)
               cmd->S_target[ct_t][1] = 1.0;   /**go backward**/
	   else
               cmd->S_target[ct_t][2] = 1.0;

           dotime += 1;

	   if  (dotime >= dotimemax) {
               behavioursuggestion = wander;          /**transition to carrybehaviour**/
               dotime = 0;
               dotimemax = 20 + (int)(drand48() * 40.0);
           }
       }

       if  (language)
           for (int i = 0; i < A[area].d_b; ++i)
               cmd->S_target[ct_t][5 * A[area].d_b + i] = 1.0;      /**behaviour report to language**/

   }


   else if  (behaviour == wander || behaviour == carry) {

      if  (motor) {

          int flag = 0;

	  if ((x < avoidancedistance) && ((ywall - fabs(y)) < avoidancedistance) &&  (y < 0.0)) {   /**robot near front wall and near right wall**/

                 /* fprintf (stderr, "\n near front wall and near right wall"); */

		 static double targetangle = -135;
                 static int turning = 0;

		 if  (newphi > 0 && newphi < 180)
                     turning = 1;
                 else if (newphi > -90 && newphi <= 0)
                     turning = -1;

                 if  (fabs (newphi - targetangle) < epsangle) {
                     cmd->S_target[ct_t][0] = 1.0; /**go forward **/
		     turning = 0;
		 }

		 if  (turning == 1) {
                     cmd->S_target[ct_t][2] = 1.0; /**turn right**/
                     flag = 1;
		}
		 if  (turning == -1) {
                     cmd->S_target[ct_t][3] = 1.0;  /**turn left**/
                     flag = 1;
		}

                /* fprintf (stderr, "\nturning=%d, flag=%d, targetangle=%.2f x=%f, y=%f, phi=%f", turning, flag, targetangle, x, y, phi); */
	  }

	  else if ((x < avoidancedistance) && ((ywall - fabs (y)) < avoidancedistance) && (y > 0.0)) {   /**robot near front wall and left wall**/

                 /* fprintf (stderr, "\n near front wall and left wall"); */

                 static double targetangle = 135; /** Turn until reach 135 degrees **/
                 static int turning = 0;

		 if  (newphi > 0 && newphi < 90)
                     turning = 1;
                 else if (newphi > -180 && newphi <= 0)
                     turning = -1;
                 if  (fabs (newphi - targetangle) < epsangle) {
                     cmd->S_target[ct_t][0] = 1.0; /**go forward **/
		     turning = 0;
		 }

		 if  (turning == 1) {
                     cmd->S_target[ct_t][2] = 1.0; /**turn right**/
                     flag = 1;
		}
		 if  (turning == -1) {
                     cmd->S_target[ct_t][3] = 1.0;  /**turn left**/
                     flag = 1;
		}

                /* fprintf (stderr, "\nturning=%d, flag=%d, targetangle=%.2f x=%f, y=%f, phi=%f", turning, flag, targetangle, x, y, phi); */
	  }

	  else if (((xwall - x) < avoidancedistance) && ((ywall - fabs(y)) < avoidancedistance) &&  (y < 0.0)) {   /**robot near end wall and right wall**/

                 /* fprintf (stderr, "\n near end wall and right wall"); */

		 static double targetangle = 45;
                 static int turning = 0;

		 if  ((newphi <= -90.0) && (newphi >= -180.0))
                     turning = 1;
                 else if ((newphi >= 0) && (newphi <= 180.0))
                     turning = -1;

		 if  (fabs (newphi - targetangle) < epsangle) {
                    cmd->S_target[ct_t][0] = 1.0;  /** go   forward **/
                    turning = 0;
		 }

		 if  (turning == 1) {
                     cmd->S_target[ct_t][2] = 1.0; /**turn right once**/
                     flag = 1;
		 }
		 if  (turning == -1) {
                     cmd->S_target[ct_t][3] = 1.0;  /**turn left**/
                     flag = 1;
		 }

          }

          else if (((xwall - x) < avoidancedistance) && ((ywall - fabs (y)) < avoidancedistance) && (y > 0.0)){   /**robot near end wall and left wall**/

                 /* fprintf (stderr, "\n near end wall and left wall"); */

		 static double targetangle = -45;
                 static int turning = 0;

		 if  ((newphi <= 0) && (newphi >= -180.0))
                     turning = 1;
                 else if ((newphi >= 90.0) && (newphi <= 180.0))
                     turning = -1;

		 if  (fabs (newphi - targetangle) < epsangle) {
                    cmd->S_target[ct_t][0] = 1.0;  /** go   forward **/
                    turning = 0;
		 }

		 if  (turning == 1) {
                     cmd->S_target[ct_t][2] = 1.0; /**turn right once**/
                     flag = 1;
		 }
		 if  (turning == -1) {
                     cmd->S_target[ct_t][3] = 1.0;  /**turn left**/
                     flag = 1;
		 }
	  }

          else if (x < avoidancedistance) {   /**robot near front wall**/

                 /* fprintf (stderr, "\n near front wall "); */

                 static double targetangle; /* = (drand48() > 0.5) ? 95.0 + drand48() * 85.0 : -95.0 - drand48() * 85.0; */
                 static int turning = 0;

		 if  (newphi > 0 && newphi < 90) {
                     turning = 1;
                     targetangle = settargetangle (90, 120);
                 }

                 else if (newphi > -90 && newphi <= 0) {
                     turning = -1;
                     targetangle = settargetangle (-90, -120);
                 }

                 if  (fabs (newphi - targetangle) < epsangle) {
                     cmd->S_target[ct_t][0] = 1.0; /**go forward **/
		     turning = 0;
		 }

		 if  (turning == 1) {
                     cmd->S_target[ct_t][2] = 1.0; /**turn right**/
                     flag = 1;
		}
		 if  (turning == -1) {
                     cmd->S_target[ct_t][3] = 1.0;  /**turn left**/
                     flag = 1;
		}

                /* fprintf (stderr, "\nturning=%d, flag=%d, targetangle=%.2f x=%f, y=%f, phi=%f", turning, flag, targetangle, x, y, phi); */
	  }

	  else if ((xwall - x) < avoidancedistance ){   /**robot near end wall**/

                 /* fprintf (stderr, "\n near end wall "); */

		 static double targetangle; /* = -85.0 + 170.0 * drand48(); */
                 static int turning = 0;

		 if  ((newphi <= -90.0) && (newphi >= -180.0)) {
                     turning = 1;
                     targetangle = settargetangle (-90, -60);
                 }
                 else if ((newphi >= 90.0) && (newphi <= 180.0)) {
                     turning = -1;
                     targetangle = settargetangle (90, 60);
                 }

		 if  (fabs (newphi - targetangle) < epsangle) {
                    cmd->S_target[ct_t][0] = 1.0;  /** go   forward **/
                    turning = 0;
		 }

		 if  (turning == 1) {
                     cmd->S_target[ct_t][2] = 1.0; /**turn right once**/
                     flag = 1;
		 }
		 if  (turning == -1) {
                     cmd->S_target[ct_t][3] = 1.0;  /**turn left**/
                     flag = 1;
		 }
	  }

	  else if  (((ywall - fabs(y)) < avoidancedistance) &&  (y < 0.0)) {   /**robot near right wall**/

                 /* fprintf (stderr, "\n near right wall "); */

		 static double targetangle; /* = -5.0 - 170.0 * drand48(); */
                 static int turning = 0;

                 if ((newphi >= 0.0) && (newphi < 90.0)) {
                     turning = -1;
                     targetangle = settargetangle (0, -30);
                 }
                 else if ((newphi >= 90.0) && (newphi <= 180.0)) {
                     turning = 1;
                     targetangle = settargetangle (-180, -150);
                 }

		  if (fabs (newphi - targetangle) < epsangle) {
		      cmd->S_target[ct_t][0] = 1.0; /** go   forward **/    
                      turning = 0;
		  }

		  if  (turning == 1) {
		      cmd->S_target[ct_t][2] = 1.0; /**turn right**/
		      flag = 1;
		  }
		  if  (turning == -1) {
                      cmd->S_target[ct_t][3] = 1.0;  /**turn left**/
                      flag = 1;
		  }
	  }

	  else if  (((ywall - fabs (y)) < avoidancedistance) && (y > 0.0)) {   /**robot near left wall**/

                  /* fprintf (stderr, "\n near left wall "); */

                  static double targetangle; /* = 5.0 + 170.0 * drand48(); */
                  static int turning = 0;

                  if  ((newphi < -90.0) && (newphi >= -180.0)) {
                      turning = -1;
                      targetangle = settargetangle (180, 150);
                  }
                  if  ((newphi >= -90.0) && (newphi <= 0.0)) {
                      turning = 1;
                      targetangle = settargetangle (0, 30);
                  }

		  if (fabs (newphi - targetangle) < epsangle) {
		      cmd->S_target[ct_t][0] = 1.0; /** go   forward **/    
                      turning = 0;
		  }

		  if  (turning == 1) {
		      cmd->S_target[ct_t][2] = 1.0; /**turn right**/
		      flag = 1;
		  }
		  if  (turning == -1) {
                      cmd->S_target[ct_t][3] = 1.0;  /**turn left**/
                      flag = 1;
		  }
	  }

	  if (flag == 0)
	     cmd->S_target[ct_t][0] = 1.0; /** go   forward **/

          /**don't allow forward AND turn -- then only turn**/
          if  (cmd->S_target[ct_t][0] == 1.0 && (cmd->S_target[ct_t][2] == 1.0 || cmd->S_target[ct_t][3] == 1.0))
              cmd->S_target[ct_t][0] = 0.0;
       }

       if  (gripper) {

           if  (behaviour == carry){

               cmd->S_target[ct_t][4] = 1.0;   /**set break-beam sensor**/
           }
       }

       if  (behaviour == wander) {

      /**  if  (fabs (newphi) < 45 && a >= 0 && a < 16 && b >= 0 && b < 24)  **/   /**for normal mode, i.e. goto docking if orange visible**/
                                                                                   /** !!! CHANGE THIS BACK !!! **/
      /**/ if  ((x > ywall - y) || (x > y + ywall) || (x > xwall / 2))      /**/   /**for "jump" mode, i.e. goto docking if gone away from wall (no physical coherence)**/
           /**front wall is NOT closest anymore**/

           {

               behavioursuggestion = docking;
           }

           if  (language) {
               if  (A[area].d_n == word_length * 39)
                   phoneme2array ("G", "OW", "-", "-", "-", "-", cmd->S_target[ct_t]);
               else
               if  (A[area].d_n == word_length * 20)
                   CEDEX2array ("g", "@", "U", "-", "-", "-", cmd->S_target[ct_t]);
               else
                   for (int i = 0; i < A[area].d_b; ++i)
                       cmd->S_target[ct_t][0 * A[area].d_b + i] = 1.0;     /**behaviour report to language**/
           }
       }



       if  (behaviour == carry) {
           if  (fabs (newphi) < 45 && a >= 0 && a < 16 && b >= 0 && b < 24)
               behavioursuggestion = put;

           if  (language) {
               if  (A[area].d_n == word_length * 39)
                   phoneme2array ("G", "OW", "-", "-", "-", "-", cmd->S_target[ct_t]);
               else
               if  (A[area].d_n == word_length * 20)
                   CEDEX2array ("g", "@", "U", "-", "-", "-", cmd->S_target[ct_t]);
               else {
                   for (int i = 0; i < A[area].d_b; ++i)
                       cmd->S_target[ct_t][3 * A[area].d_b + i] = 1.0;     /**behaviour report to language**/
                   for (int i = 0; i < A[area].d_b; ++i)
                       cmd->S_target[ct_t][0 * A[area].d_b + i] = 1.0;     /** !!! !!! TEMPORARY !!! CARRY == WANDER + CARRY !!! !!! **/
               }
           }
       }
   }



   /**vision is outside of any behaviour**/
   if  (vision) {

       /**orange-proximity, approach-speed, depart-speed, x, y_left, y_right**/

       static double old_distance = d; 

       proximity=(-1.0/16.0)*d+1.0;  /** proximity will be between 0 and 1**/

       double speed = d - old_distance;

       old_distance = d;
     
       if (proximity<0.0)
          proximity = 0.0;           /** if orange cannot be seen proximity will be set again to 0 **/

       cmd->S_target[ct_t][0] = proximity; /** proximity **/

       cmd->S_target[ct_t][1] = (d / 16.0) > 1.0 ? 1.0 : (d / 16.0);  /** scaled and cut distance **/

       if (speed > 0)
          cmd->S_target[ct_t][2] = speed;
       else if (speed < 0)      
          cmd->S_target[ct_t][3] = -speed;
     
       cmd->S_target[ct_t][4] = x / xwall;            /** scale x between 0 and 1 **/
 
       if  (y > 0) 
           cmd->S_target[ct_t][5] =  y / ywall;       /** scale y between 0 and 1 **/
       else if (y < 0)
           cmd->S_target[ct_t][6] = -y / ywall;       /** scale y between 0 and 1 **/
   }


   /**low level action reports -- outside of any behaviour**/
   if  (language) {

       if  ((A[area].d_a > 6) && (A[area].d_n != word_length * 39) && (A[area].d_n != word_length * 20)) {      /**don't write low-level action reports if phoneme coding!!!**/

           for (int i = 0; i < 4; ++i)
               cmd->S_target[ct_t][6 + i] = cmd->S_from1[1][ct_t][i];    /**the 4 wheel motor actions**/
           for (int i = 0; i < 4; ++i)
               cmd->S_target[ct_t][10 + i] = cmd->S_from1[2][ct_t][i];   /**the 4 gripper actions**/
       }
   }



   /**for the "reinforcement signal", value area; notify that a transition is made by giving a "reward"**/
   if  (reinforce) {
       if  (behavioursuggestion != behaviour) {
           cmd->S_target[ct_t][0] = 1.0;
           fprintf (stderr, "setting reinforcement");
       } else {
           cmd->S_target[ct_t][0] = 0.0;
       }
   }

   counter++;

   if  (counter == 4) {

       /**exit if escape from arena**/
       if  ((fabs (x) > xwall) || (fabs (y) > ywall))
           fprintf (stderr, "\nRobot escaped! behaviour=%d behavioursuggestion=%d, x=%f, y=%f, phi=%f", behaviour, behavioursuggestion, x, y, phi);
       if  ((fabs (x) > xwall + avoidancedistance) || (fabs (y) > ywall + avoidancedistance))
           exit (1);


       /**for restricted list of behaviours -- here only allow transitions between wander and docking**/
       {
           int behaviour_is_allowed = 0;
           for (int i = 0; i < num_allowed_behaviours; ++i)
               if  (behavioursuggestion == allowed_behaviours[i])
                   behaviour_is_allowed = 1;

           /**THIS ONLY WORKS FOR WANDER- AND DOCKING-BEHAVIOURS**/
           if  (num_allowed_behaviours == 2)
               if  (behaviour_is_allowed == 0) {
                   if  (behaviour == docking)
                       behavioursuggestion = wander;
                   if  (behaviour == wander)
                       behavioursuggestion = docking;
               }               
           if  (num_allowed_behaviours == 1)
               if  (allowed_behaviours[0] == 0)
                   behavioursuggestion = wander;
               if  (allowed_behaviours[0] == 1)
                   behavioursuggestion = docking;
       }


       if  (behavioursuggestion != behaviour)
           fprintf (stderr, "\ntransition from behav %d to behav %d", behaviour, behavioursuggestion);

       behaviour = behavioursuggestion;

       counter = 0;
   }

    return (DOUBLE)(0);
}

#undef word_length


/****************************** total_real_world *****************************/

DOUBLE total_real_world (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    double x   = cmd->S_from1[0][ct_t][0];
    double y   = cmd->S_from1[0][ct_t][1];
    double phi = cmd->S_from1[0][ct_t][2];

    double gripperheight = cmd->S_from1[0][ct_t][4];

    double advance = 0.0;
    double phi_advance = 0.0;
    const double pi = 3.1415926535;

    /**change the position randomly to avoid deadlocks during docking**/
    x += -0.05 * drand48() * 0.1;
    y += -0.05 * drand48() * 0.1;

/*
    fprintf (stderr, "\nt_r_w: x=%.2f y=%.2f phi=%.2f ", x, y, phi);
    fprintf (stderr, "f=%.2f b=%.2f l=%.2f r=%.2f ", cmd->S_from1[1][ct_t][0], cmd->S_from1[1][ct_t][1], cmd->S_from1[1][ct_t][2], cmd->S_from1[1][ct_t][3]);
*/

    if  (cmd->S_from1[1][ct_t][0] == 1.0)
        advance = cmd->quantum[0][0];
    if  (cmd->S_from1[1][ct_t][1] == 1.0)
        advance = -1.0 * cmd->quantum[0][0];
    if  (cmd->S_from1[1][ct_t][2] == 1.0)
        phi_advance = cmd->quantum[0][1];
    if  (cmd->S_from1[1][ct_t][3] == 1.0)
        phi_advance = -1.0 * cmd->quantum[0][1];

    cmd->S_target[ct_t][0] = x - cos (phi) * advance;   /**distance to table**/
    cmd->S_target[ct_t][1] = y - sin (phi) * advance;   /**lateral offset   **/
    cmd->S_target[ct_t][2] = phi + phi_advance;         /**angle            **/
    if  (cmd->S_target[ct_t][2] > pi)
        cmd->S_target[ct_t][2] -= 2.0 * pi;
    if  (cmd->S_target[ct_t][2] < -pi)
        cmd->S_target[ct_t][2] += 2.0 * pi;

    /**gripperopen**/
    cmd->S_target[ct_t][3] = cmd->S_from1[0][ct_t][3];
    if  (cmd->S_from1[2][ct_t][0] == 1.0)               /**closing gripper action**/
        cmd->S_target[ct_t][3] = 0.0;                   /**gripperopen = 0**/
    if  (cmd->S_from1[2][ct_t][1] == 1.0)               /**opening gripper action**/
        cmd->S_target[ct_t][3] = 1.0;                   /**gripperopen = 1**/

    /**gripperheight**/
    if  (cmd->S_from1[2][ct_t][2] == 1.0)               /**lifting gripper action**/
        gripperheight += 0.1;                           /**lifting speed**/
    if  (cmd->S_from1[2][ct_t][3] == 1.0)               /**lowering gripper action**/
        gripperheight -= 0.1;                           /**lowering speed (relevant only whether 0.5)**/

    /**just in case ...; should never be executed**/
    if  (gripperheight > 1.1) {
        fprintf (stderr, "\ntotal_real_world finds too high gripperheight!\n");
        gripperheight = 1.0;
    }

    cmd->S_target[ct_t][4] = gripperheight;

    return (DOUBLE)(0);
}



/****************************** total_init_behaviour *************************/
/* initialize in the real world                                              */

DOUBLE total_init_behaviour (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    cmd->S_target[ct_t][0] = 15.0;                       /**x**/
    cmd->S_target[ct_t][1] = -5.0 + drand48() * 10.0;    /**y**/
    cmd->S_target[ct_t][2] = 0.0;                        /**phi**/
    cmd->S_target[ct_t][3] = 1.0;                        /**gripperopen**/
    cmd->S_target[ct_t][4] = 0.75;                       /**0.75**/

    return (DOUBLE)(0);
}



/****************************** total_high_vision ****************************/
/* From total_imag_coord; Compute high-level vision input given (x,y,phi).   */

DOUBLE total_high_vision (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    double x   = cmd->S_from1[0][ct_t][0];
    double y   = cmd->S_from1[0][ct_t][1];

    double d     = sqrt (x*x + y*y);
    double d_max = sqrt (12*12 + 16*16);

    double phi   = cmd->S_from1[0][ct_t][2];
    double theta = atan (y/x) - phi;

    int area = cmd->area;

    if  (A[area].d_b == 1) {

        if  ((A[area].d_n != 7) && (A[area].d_n != 10) && (A[area].d_n != 12))
            fprintf (stderr, "\ntotal_high_vision: wrong dimensions!\n\n\n");

        cmd->S_target[ct_t][0] = 1.0 - d / d_max;                   /**orange-proximity**/
        cmd->S_target[ct_t][1] = d / d_max;                         /**orange-distance**/
        cmd->S_target[ct_t][2] = 0.0;                               /**approach-speed**/
        cmd->S_target[ct_t][3] = 0.0;                               /**leaving-speed**/
        cmd->S_target[ct_t][4] = x / 16.0;                          /**x/xwall**/
        cmd->S_target[ct_t][5] = (y > 0) ? y / 12.0 : 0.0;          /**y/ywall**/
        cmd->S_target[ct_t][6] = (y < 0) ? -y / 12.0 : 0.0;         /**-y/ywall**/

        if  (A[area].d_n == 10) {
            cmd->S_target[ct_t][7] = 1.0 - cmd->S_target[ct_t][4];  /**next also first: 1-x/xwall**/
            cmd->S_target[ct_t][8] = 1.0 - cmd->S_target[ct_t][5];  /**1-y/ywall**/
            cmd->S_target[ct_t][9] = 1.0 - cmd->S_target[ct_t][6];  /**1--y/ywall**/
        }

        if  (A[area].d_n == 12) {
            cmd->S_target[ct_t][10] = (theta < 0.0) ? -theta / 1.5707963 : 0.0;  /**orange-to-left**/
            cmd->S_target[ct_t][11] = (theta > 0.0) ?  theta / 1.5707963 : 0.0;  /**orange-to-right**/
        }
    }


    if  ((A[area].d_b == 2) || (A[area].d_b == 3) || (A[area].d_b == 4)) {
        double sigma = 3.0;   /** !!! was 1.0 !!! **/

        for (int a = 0; a < A[area].d_a; a++) {
            double diffA = 999;
            double diffB = 999;
            double x_max     = 16;
            double y_max     = 24;
            double theta_max = 1.5707963;

            /**all non-periodic boundary**/

            /**absolute x-coordinate**/
            if  (cmd->quantum[0][0] == 2) {

                diffA = (double)(a) - (x/x_max * (double)(A[area].d_a));                             /**using x_max designed after retina for docking**/
                diffB = (double)(a) - ((y+y_max)/(2.0*y_max) * (double)(A[area].d_a));

            } else {

                /**relative x-coordinate to nearest wall**/
                if  (cmd->quantum[0][0] != 3)
                    fprintf (stderr, "\ntotal_high_vision wants proper paramter\n");

                const double ywall = 28 - 1; /*16;*/  /**if values change then**/
                const double xwall = 56 + 10; /*48;*/  /**check also 2 occurrences above!**/

                /**front wall is closest**/
                if  ((x < ywall - y) && (x < y + ywall) && (x < xwall / 2)) {
                    diffA = (double)(a) - (x/xwall * 2.0) * (double)(A[area].d_a);                  /**using xwall designed after arena for wander; x < xwall/2 in this part of the arena**/
                    diffB = (double)(a) - ((y/ywall)/2.0+0.5) * (double)(A[area].d_a);
                }

                /**end wall is closest**/
                if  ((xwall - x < ywall - y) && (xwall - x < y + ywall) && (x > xwall / 2)) {
                    diffA = (double)(a) - (1.0 - (x/xwall - 0.5) * 2.0) * (double)(A[area].d_a);     /**x > xwall/2 in end half of arena**/
                    diffB = (double)(a) - (1.0 - ((y/ywall)/2.0+0.5)) * (double)(A[area].d_a);
                }

                /**left wall is closest**/
                if  ((x > ywall - y) && (xwall - x > ywall - y) && (y > 0.0)) {
                    diffA = (double)(a) - (ywall - y)/ywall * (double)(A[area].d_a);                 /**y/Xwall makes sense only if xwall>2ywall! Use y/Ywall to use all neurons**/
                    diffB = (double)(a) - (x/xwall) * (double)(A[area].d_a);
                }

                /**right wall is closest**/
                if  ((x > y + ywall) && (xwall - x > y + ywall) && (y < 0.0)) {
                    diffA = (double)(a) - (ywall + y)/ywall * (double)(A[area].d_a);                 /**y < 0 in this part of the arena**/
                    diffB = (double)(a) - (1.0 - x/xwall) * (double)(A[area].d_a);
                }
            }

            if  (diffA == 999)
                fprintf (stderr, "\ntotal_high_vision warning: diffA not set!");

            /**x**/
            cmd->S_target[ct_t][a * A[area].d_b + 0] = exp (-0.5 * (diffA*diffA) / (sigma*sigma));

            /**y**/
            cmd->S_target[ct_t][a * A[area].d_b + 1] = exp (-0.5 * (diffB*diffB) / (sigma*sigma));

            /**distance**/
            if  (A[area].d_b >= 3) {
                diffA = (double)(a) - (d/d_max * (double)(A[area].d_a));
                cmd->S_target[ct_t][a * A[area].d_b + 2] = exp (-0.5 * (diffA*diffA) / (sigma*sigma));
            }

            /**theta (visually perceived angle)**/
            if  (A[area].d_b == 4) {
                diffA = (double)(a) - ((theta+theta_max)/(2.0*theta_max) * (double)(A[area].d_a));
                cmd->S_target[ct_t][a * A[area].d_b + 3] = exp (-0.5 * (diffA*diffA) / (sigma*sigma));
            }
        }
  
    }

    return (DOUBLE)(0);
}


/****************************** total_mult_const_if **************************/
/* Multiply by q[0][0] if any activation in input area is non-zero.          */

DOUBLE total_mult_const_if (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

   int scale_condition = 0;
   int inarea = cmd->n_from1[0];
   int comparea = cmd->n_from2[0];

   if  (inarea != cmd->area)
       fprintf (stderr, "\nwrong use of total_mult_const_if!\n");

   for (int i = 0; i < A[comparea].d_n; ++i)
       if  (cmd->S_from2[0][ct_t][i] != 0.0)
           scale_condition = 1;

   if  (scale_condition)
       fprintf (stderr, "s");

   if  (scale_condition)
       for (int i = 0; i < A[inarea].d_n; ++i)
           cmd->S_target[ct_t][i] = cmd->S_from1[0][ct_t][i] * cmd->quantum[0][0];
   else
       for (int i = 0; i < A[inarea].d_n; ++i)
           cmd->S_target[ct_t][i] = cmd->S_from1[0][ct_t][i];               /**NOT YET TESTED!!!**/

    return (DOUBLE)(0);
}



/****************************** total_set_nth ********************************/
/* Set the q[0][j] units to q[1][0] or q[1][j].                              */
/* E.g. quantums like: "3+5, 1.0" sets units indexed 3 and 5 to 1.0.         */
/*                 or: "3+5, 1.0+0.5" sets unit 3 to 1.0 and unit 5 to 0.5.  */

DOUBLE total_set_nth (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {
    int area = cmd->area;

   /**check**/
   for (int j = 0; j < cmd->anz_quant[0]; ++j)
       if  ((int)(cmd->quantum[0][j]) >= A[area].d_n)
           fprintf (stderr, "\ntotal_set_nth: n larger than area size!\n");

   if  ((cmd->anz_quant[1] != 1) && (cmd->anz_quant[1] != cmd->anz_quant[0]))
       fprintf (stderr, "\ntotal_set_nth: anz_quant[1,0] inconsistent!\n");

   /**does NOT anymore init (the others) to zero**/

   /**set values**/
   if  (cmd->anz_quant[1] == 1) /**all get same value**/
       for (int j = 0; j < cmd->anz_quant[0]; ++j)
           cmd->S_target[ct_t][(int)(cmd->quantum[0][j])] = cmd->quantum[1][0];
   else                         /**each gets indiviual value**/
       for (int j = 0; j < cmd->anz_quant[0]; ++j)
           cmd->S_target[ct_t][(int)(cmd->quantum[0][j])] = cmd->quantum[1][j];

   return (DOUBLE)(0);
}


/****************************** total_swap_nth *******************************/
/* Swap the indeces from S_from1 and write result to S_target.               */
/* Act's with index q[j][0] and q[j][1] are swapped.                         */
/* E.g. quantums: "0+3, 1+2" swaps the act's of unit 0 with 3 and 1 with 2.  */

DOUBLE total_swap_nth (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {
    int area = cmd->area;

   /**check**/
   for (int j = 0; j < cmd->anz_quantums; ++j) {
       if  ((int)(cmd->anz_quant[j]) != 2)
           fprintf (stderr, "\ntotal_swap_nth: wrong usage!\n");
       if  (((int)(cmd->quantum[j][0]) >= A[area].d_n) || ((int)(cmd->quantum[j][1]) >= A[area].d_n))
           fprintf (stderr, "\ntotal_swap_nth: n larger than area size!\n");
   }

   /**init**/
   for (int i = 0; i < A[area].d_n; ++i)
       cmd->S_target[ct_t][i] = cmd->S_from1[0][ct_t][i];

   /**swap values -- use temporary values so that same source and target is possible**/
   for (int j = 0; j < cmd->anz_quantums; ++j) {
       int    left_idx  = (int)(cmd->quantum[j][0]);
       int    right_idx = (int)(cmd->quantum[j][1]);
       DOUBLE left_val  = cmd->S_from1[0][ct_t][left_idx];
       DOUBLE right_val = cmd->S_from1[0][ct_t][right_idx];
       cmd->S_target[ct_t][left_idx]  = right_val;
       cmd->S_target[ct_t][right_idx] = left_val;
   }

   return (DOUBLE)(0);
}


/****************************** total_copy_nth *******************************/
/* Copies acts with indeces given from S_from1 to S_target.                  */
/* E.g. quantums "0+3" copys act's of units 0,3 from S_from1 to S_target.    */
/* In special mode, e.g. quantums "0+3, 1+2" copies acts 0->1 and 3->2.      */

DOUBLE total_copy_nth (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {
   int area = cmd->area;

   /**normal mode: unit indexes are the same for source and target**/
   if  (cmd->anz_quantums == 1) {

       /**check**/
       for (int j = 0; j < cmd->anz_quant[0]; ++j) {
           if  ((int)(cmd->quantum[0][j]) >= A[area].d_n)
              fprintf (stderr, "\ntotal_copy_nth: n larger than area size!\n");
       }

       /**copy values**/
       for (int j = 0; j < cmd->anz_quant[0]; ++j)
           cmd->S_target[ct_t][(int)(cmd->quantum[0][j])] = cmd->S_from1[0][ct_t][(int)(cmd->quantum[0][j])];
   }

   /**special mode: individual indexes for source and target**/
   if  (cmd->anz_quantums == 2) {

       /**no checking at all!**/

       /**copy values**/
       for (int j = 0; j < cmd->anz_quant[0]; ++j)
           cmd->S_target[ct_t][(int)(cmd->quantum[1][j])] = cmd->S_from1[0][ct_t][(int)(cmd->quantum[0][j])];
   }


   return (DOUBLE)(0);
}


/****************************** total_rand_nth *******************************/
/* Set units q[0][j] to rand between q[1][j] .. q[2][j], all others to zero. */
/* E.g. quantums like: "0+2, 0+-3.1415, 1+3.1415"                            */
/* sets unit 0 between 0 .. 1 and unit 2 between -PI .. PI (and unit 1 to 0).*/

DOUBLE total_rand_nth (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {
    int area = cmd->area;

   /**check**/
   for (int j = 0; j < cmd->anz_quant[0]; ++j)
       if  ((int)(cmd->quantum[0][j]) >= A[area].d_n)
           fprintf (stderr, "\ntotal_rand_nth: n larger than area size!\n");

   if  ((cmd->anz_quant[1] != cmd->anz_quant[0]) || (cmd->anz_quant[2] != cmd->anz_quant[0]))
       fprintf (stderr, "\ntotal_rand_nth: anz_quant[1,0] or anz_quant[2,0] inconsistent!\n");

   /**init to zero**/
   for (int i = 0; i < A[area].d_n; ++i)
       cmd->S_target[ct_t][i] = 0.0;

   /**set values**/
   for (int j = 0; j < cmd->anz_quant[0]; ++j)
       cmd->S_target[ct_t][(int)(cmd->quantum[0][j])] = cmd->quantum[1][j]
                                                      + drand48() * (cmd->quantum[2][j] - cmd->quantum[1][j]);

    return (DOUBLE)(0);
}


/****************************** total_compare_with ***************************/
/* Tests whether S_from1[0] and S_from1[1] match the same of S_from1[1..n].  */
/* (Used to check whether behaviours are recognised correctly by language.)  */
/* Writes to stderr every then and when. Dummy target.                       */

DOUBLE total_compare_with (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

  if  (cmd->anz_from1 != 2)
      fprintf (stderr, "\nwrong use of total_compare_with\n");

    int area = cmd->area;

  double dist_min;
  int best_match_0 = 999, best_match_1 = 998;
  static int matches[10]    = {0,0,0,0,0,0,0,0,0,0};
  static int mismatches[10][10] = {{0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},
                                   {0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},
                                   {0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},
                                   {0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0},
                                   {0,0,0,0,0,0,0,0,0,0},{0,0,0,0,0,0,0,0,0,0}};
  static int count = 0;

  dist_min = 999999999.0;
  for (int ct_l = 0; ct_l < cmd->anz_from2; ++ct_l) {

      double dist = 0.0;
      for (int ct_n = 0; ct_n < A[area].d_n; ++ct_n) {
          dist += (cmd->S_from1[0][ct_t][ct_n] - cmd->S_from2[ct_l][ct_t][ct_n])
                * (cmd->S_from1[0][ct_t][ct_n] - cmd->S_from2[ct_l][ct_t][ct_n]);
      }

      if  (dist < dist_min) {
          dist_min = dist;
          best_match_0 = ct_l;
      }
  }

  dist_min = 999999999.0;
  for (int ct_l = 0; ct_l < cmd->anz_from2; ++ct_l) {

      double dist = 0.0;
      for (int ct_n = 0; ct_n < A[area].d_n; ++ct_n) {
          dist += (cmd->S_from1[1][ct_t][ct_n] - cmd->S_from2[ct_l][ct_t][ct_n])
                * (cmd->S_from1[1][ct_t][ct_n] - cmd->S_from2[ct_l][ct_t][ct_n]);
      }

      if  (dist < dist_min) {
          dist_min = dist;
          best_match_1 = ct_l;
      }
  }


  for (int ct_l = 0; ct_l < cmd->anz_from2; ++ct_l) {

      if  (best_match_0 == ct_l) {
          if  (best_match_1 == ct_l)
              matches[ct_l] += 1;
          else
              mismatches[ct_l][best_match_1] += 1;
      }
  }

  if  (count % 100 == 0)
      for (int ct_l = 0; ct_l < cmd->anz_from2; ++ct_l) {
          fprintf (stderr, "\nbehav %d: matches=%d, mismatches", ct_l, matches[ct_l]);
          for (int i = 0; i < cmd->anz_from2; ++i)
              fprintf (stderr, " to %d = %d  ", i, mismatches[ct_l][i]);
      }

  count += 1;

    return (DOUBLE)(0);
}


/****************************** total_compare_count **************************/
/* Tests whether S_from1[0] and S_from2[0] are the same.                     */
/* (Used to check whether behaviours are performed correctly by learner.)    */
/* Writes to stderr every then and when. Dummy target.                       */

DOUBLE total_compare_count (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {
    int area = cmd->area;

  static int count = 0;
  static int total_count = 0;
  double max_dist_threshold = 0.0001;

  double dist = 0.0;
  for (int ct_n = 0; ct_n < A[area].d_n; ++ct_n)
      dist += (cmd->S_from1[0][ct_t][ct_n] - cmd->S_from2[0][ct_t][ct_n])
            * (cmd->S_from1[0][ct_t][ct_n] - cmd->S_from2[0][ct_t][ct_n]);
  dist = sqrt (dist);

  if  (dist > max_dist_threshold)
      count += 1;

  total_count += 1;

  if  (total_count % 100 == 0)
      fprintf (stderr, "\ntotal_count = %d (examples); count = %d (errors) ", total_count, count);



  static int pluserrors[4] = {0,0,0,0};
  static int minuserrors[4] = {0,0,0,0};
  static int onleft[4] = {0,0,0,0};
  static int onright[4] = {0,0,0,0};

  for (int ct_n = 0; ct_n < 4; ++ct_n) {
      if  (cmd->S_from1[0][ct_t][ct_n] - cmd->S_from2[0][ct_t][ct_n] > max_dist_threshold)
          pluserrors[ct_n] += 1;
      if  (cmd->S_from2[0][ct_t][ct_n] - cmd->S_from1[0][ct_t][ct_n] > max_dist_threshold)
          minuserrors[ct_n] += 1;
  }

  for (int ct_n = 0; ct_n < 4; ++ct_n) {
      if  (cmd->S_from1[0][ct_t][ct_n] > max_dist_threshold)
          onleft[ct_n] += 1;
      if  (cmd->S_from2[0][ct_t][ct_n] > max_dist_threshold)
          onright[ct_n] += 1;
  }

  if  (total_count % 100 == 0)
      fprintf (stderr, "\npluserrors = %d %d %d %d,  minuserrors = %d %d %d %d  ",
                    pluserrors[0], pluserrors[1], pluserrors[2], pluserrors[3],
                    minuserrors[0], minuserrors[1], minuserrors[2], minuserrors[3]);
  if  (total_count % 100 == 0)
      fprintf (stderr, "\nonleft = %d %d %d %d,  onright = %d %d %d %d  ",
                    onleft[0], onleft[1], onleft[2], onleft[3],
                    onright[0], onright[1], onright[2], onright[3]);

    return (DOUBLE)(0);
}




/****************************** total_retina_to_SC ***************************/
/* Projects retina image log-polar-like to SC.                               */
/* Averages between 4 pixels -- BAD INTERPOLATION -- needs more care:        */
/* Difficult: e.g. one SC cell might receive input from just one foveal cell,*/
/* while another might receive input from dozens of peripheral cells.        */
/* Rudimentary weight-matrix needed. Do that in weight_list.c?               */

DOUBLE total_retina_to_SC (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    float s = 0.25;
    float q_rc = 5.0 / log(s*90+1);  /*=1.5837819*/
    float q_ml = 5.0 / 180.0;

    int area = cmd->area;
    int inarea = cmd->n_from1[0];
    static int *apixels_a = NULL, *apixels_b = NULL;      /**anchor pixels of all SC neurons on the retina**/
    static int *aunits_a = NULL, *aunits_b = NULL;        /**anchor units of all retina neurons in the SC**/
    static int firsttime = 1;
    static int *scaffold_SC;
    static DOUBLE *push_act_SC;

     if  (firsttime) {

         /**Method 1: for each SC cell, define an anchor pixel in the retina**/
         apixels_a = (int *)malloc (A[area].d_a * A[area].d_b * sizeof (int));
         apixels_b = (int *)malloc (A[area].d_a * A[area].d_b * sizeof (int));

         for (int a = 0; a < A[area].d_a; ++a)
         for (int b = 0; b < A[area].d_b; ++b) {

             /**get point position on SC**/
             float SC_rc =  5.0 * (float)a / (float)(A[area].d_a-0.000001);            /**range 0..5 [mm]**/
             float SC_ml = -5.0 + 10.0 * (float)b / (float)(A[area].d_b-0.000001);     /**range -5..5 [mm]**/

             /**map from SC to retina**/
             float ret_r = 1 / s * (exp(SC_rc / q_rc) - 1.0);               /**range 0..90 [deg]**/
             float ret_t = SC_ml / q_ml;                                    /**range -180..180 [deg]**/

             /**convert to cartesian coordinates**/
             float ret_x = ret_r * cos (ret_t / 180.0 * M_PI);
             float ret_y = ret_r * sin (ret_t / 180.0 * M_PI);

             /**determine anchor points**/
             float point_a = (ret_x + 90.0) / 180.0 * (float)(A[inarea].d_a-0.000001);
             float point_b = (ret_y + 90.0) / 180.0 * (float)(A[inarea].d_b-0.000001);

             /**determine anchor pixels**/
             apixels_a[a * A[area].d_b + b] = (int)point_a;
             apixels_b[a * A[area].d_b + b] = (int)point_b;
         }

         /** ... vice versa ... **/

         /**Method 2: for each retina cell, determine a target unit on SC**/
         aunits_a = (int *)malloc (A[inarea].d_a * A[inarea].d_b * sizeof (int));
         aunits_b = (int *)malloc (A[inarea].d_a * A[inarea].d_b * sizeof (int));
         scaffold_SC = (int *)calloc (A[area].d_n, sizeof (int));

         for (int a = 0; a < A[inarea].d_a; ++a)
         for (int b = 0; b < A[inarea].d_b; ++b) {

             /**get point position on retina**/
             float ret_x = -90.0 + 180.0 * (float)a / (float)(A[inarea].d_a-1);     /**range -90..90 [deg]**/
             float ret_y = -90.0 + 180.0 * (float)b / (float)(A[inarea].d_b-1);     /**range -90..90 [deg]**/

             /**convert to polar coordinates**/
             float ret_r = sqrt (ret_x*ret_x + ret_y*ret_y);                        /**range 0..90 [deg]**/
           /*float ret_t = atan (ret_y/ret_x) / M_PI * 2.0 * 180.0;*/               /**range -180..180 [deg]**/
             DOUBLE *dummy=NULL;
             float ret_t = local_atan_to_2pi (dummy, ret_y, ret_x) / M_PI * 180.0;
             if  (ret_t >= 180.0)
                 ret_t -= 360.0;

             /**determine target position**/
             float sc_rc = q_rc * log (s * ret_r + 1.0);                            /**range 0..5 [mm]**/
             float sc_ml = q_ml * ret_t;                                            /**range -5..5 [mm]**/

             /**scale**/
             sc_rc = sc_rc / 5.0 * (float)(A[area].d_a-1);                    /**range 0..A[area].d_a ("-1" despite cast at next step?!)**/    /** !!! hampered with !!! **/
             sc_ml = (5.0 + sc_ml) / 10.0 * (float)(A[area].d_b-0.000001);           /**range 0..A[area].d_b                    "              **/

             int ct_in = a * A[inarea].d_b + b;

             /**determine anchor pixels**/
             if  (ret_r <= 90.0)
                 aunits_a[ct_in] = (int)sc_rc;
             else
                 aunits_a[ct_in] = -1;                              /**round retina; at corners set "-1" as flag**/

             aunits_b[ct_in] = (int)sc_ml;

             if  (aunits_a[ct_in] != -1)
                 scaffold_SC[  aunits_a[ct_in] * A[area].d_b + aunits_b[ct_in]  ] += 1;

             /**check**/
             if  (((int)sc_rc >= A[area].d_a) && (aunits_a[a * A[inarea].d_b + b] != -1))
                 fprintf (stderr, "\n(int)sc_rc=%d out of bounds  d_a=%d\n", (int)sc_rc, A[area].d_a);
             if  (((int)sc_ml >= A[area].d_b) || ((int)sc_ml < 0))
                 fprintf (stderr, "\n(int)sc_ml=%d out of bounds  d_b=%d\n", (int)sc_ml, A[area].d_b);
         }

         push_act_SC = (DOUBLE *)malloc (A[area].d_n * sizeof (DOUBLE));

         firsttime = 0;
     }


     if  (cmd->quantum[0][0] == 1) {

         /**Method 1: for all SC units, get one pixel from retina <-- not all retinal pixels are used**/
         for (int ct_n = 0; ct_n < A[area].d_n; ++ct_n)
             cmd->S_target[ct_t][ct_n] = cmd->S_from1[0][ct_t][  apixels_a[ct_n] * A[inarea].d_b + apixels_b[ct_n]  ];

     } else {

         /**Method 2: for all retinal pixels, push activation to one SC unit <-- not all SC units are used**/
         for (int ct_n = 0; ct_n < A[area].d_n; ++ct_n)
             push_act_SC[ct_n] = 0.0;

         /**push retinal activations to matching SC cells**/
         for (int ct_in = 0; ct_in < A[inarea].d_n; ++ct_in)
             if  (aunits_a[ct_in] != -1)
                 push_act_SC[  aunits_a[ct_in] * A[area].d_b + aunits_b[ct_in]  ] += cmd->S_from1[0][ct_t][ct_in];

         for (int ct_n = 0; ct_n < A[area].d_n; ++ct_n) {
             if  (scaffold_SC[ct_n] != 0)
                 cmd->S_target[ct_t][ct_n] = push_act_SC[ct_n] / (double)(scaffold_SC[ct_n]) ;
             else
                 if  (push_act_SC[ct_n] != 0.0)
                     fprintf (stderr, "\n inconsistency with scaffold_SC in total_retina_to_SC");
         }

         /**Method 3: average between Methods 1 and 2**/
         if  (cmd->quantum[0][0] == 0) {
             for (int ct_n = 0; ct_n < A[area].d_n; ++ct_n) {

                 if  (scaffold_SC[ct_n] == 0) { /**use only Method 1 if Method 2 doesn't fill the unit**/
                     cmd->S_target[ct_t][ct_n] = cmd->S_from1[0][ct_t][  apixels_a[ct_n] * A[inarea].d_b + apixels_b[ct_n]  ];

                 } else { /**average between both Methods**/
                     cmd->S_target[ct_t][ct_n] += cmd->S_from1[0][ct_t][  apixels_a[ct_n] * A[inarea].d_b + apixels_b[ct_n]  ];
                     cmd->S_target[ct_t][ct_n] *= 0.5;
                 }
             }
         }
     }

    return (DOUBLE)(0);
}



/****************************** total_population_motor_row *******************/
/* Get population response from each row and send to one target neuron each. */
/* Hence, #rows must be #target neurons.                                     */
/* q[0][i] - q[1][i] = motor- (output-) ranges for each row.                 */

DOUBLE total_population_motor_row (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {
    int X, Y;
    int inarea = cmd->n_from1[0];
    int d_n = A[cmd->area].d_n;
    if  ((cmd->anz_quantums != 2) && (cmd->anz_quantums != 3))
        fprintf (stderr, "\ntotal_population_motor_row: wrong no of quant ");
    if  (d_n != cmd->anz_quant[0])
        fprintf (stderr, "\ntotal_population_motor_row: wrong no of quant[0] ");
    if  (d_n != cmd->anz_quant[1])
        fprintf (stderr, "\ntotal_population_motor_row: wrong no of quant[1] ");
    if  (cmd->anz_quantums == 2)
    if  (d_n != A[inarea].d_a)
        fprintf (stderr, "\ntotal_population_motor_row: no of inarea's rows must be no target neurons ");

    /**standard version -- every target neuron and every row**/
    if  (cmd->anz_quantums == 2)
    for (X = 0; X < d_n; X++) {
        DOUBLE response = 0.0;
        DOUBLE tot_act  = 0.0;
        for (Y = 0; Y < A[inarea].d_b; Y++) {
            float value = cmd->quantum[0][X] + (cmd->quantum[1][X] - cmd->quantum[0][X]) * (float)Y / (float)(A[inarea].d_b - 1);
            response += cmd->S_from1[0][ct_t][X * A[inarea].d_b + Y] * value;
            tot_act += fabs (cmd->S_from1[0][ct_t][Y]);
        }
        cmd->S_target[ct_t][X] = response / tot_act;
    }

    /**do this only for one target neuron and one row**/
    if  (cmd->anz_quantums == 3)
    {   X = (int)cmd->quantum[2][0];
        DOUBLE response = 0.0;
        DOUBLE tot_act  = 0.0;
        for (Y = 0; Y < A[inarea].d_b; Y++) {
            float value = cmd->quantum[0][X] + (cmd->quantum[1][X] - cmd->quantum[0][X]) * (float)Y / (float)(A[inarea].d_b - 1);
            response += cmd->S_from1[0][ct_t][X * A[inarea].d_b + Y] * value;
            tot_act += fabs (cmd->S_from1[0][ct_t][Y]);
        }
        cmd->S_target[ct_t][X] = response / tot_act;
    }


    return (DOUBLE)(0);
}


/****************************** total_population_motor_2D ********************/
/* */

DOUBLE total_population_motor_2D (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int X, Y;
    int inarea = cmd->n_from1[0];

    if  (cmd->anz_quantums != 2)
        fprintf (stderr, "\ntotal_population_motor_2D: wrong no of quant ");
    if  (cmd->anz_quant[0] != 2)
        fprintf (stderr, "\ntotal_population_motor_2D: wrong no of quant[0] ");
    if  (cmd->anz_quant[1] != 2)
        fprintf (stderr, "\ntotal_population_motor_2D: wrong no of quant[1] ");
    if  (A[cmd->area].d_n != 2)
        fprintf (stderr, "\ntotal_population_motor_2D: must be 2 target units ");

    DOUBLE response_x = 0.0;
    DOUBLE response_y = 0.0;
    DOUBLE tot_act    = 0.0;

    for (X = 0; X < A[inarea].d_a; X++)
        for (Y = 0; Y < A[inarea].d_b; Y++) {

            float value_x = cmd->quantum[0][0] + (cmd->quantum[1][0] - cmd->quantum[0][0]) * (float)X / (float)(A[inarea].d_a - 1);
            float value_y = cmd->quantum[0][1] + (cmd->quantum[1][1] - cmd->quantum[0][1]) * (float)Y / (float)(A[inarea].d_b - 1);

            response_x += cmd->S_from1[0][ct_t][X * A[inarea].d_b + Y] * value_x;
            response_y += cmd->S_from1[0][ct_t][X * A[inarea].d_b + Y] * value_y;
            tot_act    += cmd->S_from1[0][ct_t][X * A[inarea].d_b + Y];
        }

    if  (tot_act != 0.0) {
        cmd->S_target[ct_t][0] = response_x / tot_act;
        cmd->S_target[ct_t][1] = response_y / tot_act;
    } else {
        cmd->S_target[ct_t][0] = 0.0;
        cmd->S_target[ct_t][1] = 0.0;
    }

    return (DOUBLE)(0);
}


/**************************** total_scalar_mult **************************/
/* Scalar multiply (S_from1,S_from2) and write to unit q00 of S_target.  */

DOUBLE total_scalar_mult (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int inarea1 = cmd->n_from1[0];
    int inarea2 = cmd->n_from2[0];
    int area = cmd->area;
    int target_neuron = (int)cmd->quantum[0][0];

    if  (A[inarea1].d_n != A[inarea2].d_n)
        fprintf (stderr, "\ntotal_scalar_mult wants two input areas of same size!\n");
    if  (target_neuron >= A[area].d_n)
        fprintf (stderr, "\ntotal_scalar_mult: parameter exceeds target area size!\n");

    DOUBLE sp = 0.0;
    for (int i = 0; i < A[inarea1].d_n; ++i)
        sp += cmd->S_from1[0][ct_t][i] * cmd->S_from2[0][ct_t][i];
 
    cmd->S_target[ct_t][target_neuron] = sp;

    return (DOUBLE)(0);
}


/************************** total_change_active_by ***********************/
/* At the positions where S_from1 is active (should be blobs on SC area) */
/* do modifications of the saccadic gain/offset represented in S_target. */
/* S_target is for vertical if q00=1 and S_from2 unit 0 or 1 is active.  */
/* S_target is for horiz'al if q00=2 and S_from2 unit 2 or 3 is active.  */

DOUBLE total_change_active_by (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int inarea1 = cmd->n_from1[0];
    int inarea2 = cmd->n_from2[0];
    int area = cmd->area;
    float scale_factor = cmd->quantum[1][0]; 

    if  (A[area].d_n != A[inarea1].d_n)
        fprintf (stderr, "\ntotal_change_active_by: target area should be same as the 1st input area!\n");
    if  (A[inarea2].d_n != 4)
        fprintf (stderr, "\ntotal_change_active_by: 2nd input area must have 4 units!\n");
    if  (cmd->S_from2[0][ct_t][0] + cmd->S_from2[0][ct_t][1] + cmd->S_from2[0][ct_t][2] + cmd->S_from2[0][ct_t][3] != 1.0)
        fprintf (stderr, "\ntotal_change_active_by: wrong activation %.2f %.2f %.2f %.2f on second input!\n",
                          cmd->S_from2[0][ct_t][0], cmd->S_from2[0][ct_t][1], cmd->S_from2[0][ct_t][2], cmd->S_from2[0][ct_t][3]);

    /**change vertical saccade offset**/
    if  (cmd->quantum[0][0] == 1.0) {

        if  (cmd->S_from2[0][ct_t][0] == 1.0)
            for (int i = 0; i < A[area].d_n; ++i)
                cmd->S_target[ct_t][i] += cmd->S_from1[0][ct_t][i] * scale_factor;
 
        if  (cmd->S_from2[0][ct_t][1] == 1.0)
            for (int i = 0; i < A[area].d_n; ++i)
                cmd->S_target[ct_t][i] -= cmd->S_from1[0][ct_t][i] * scale_factor;
    }

    /**change horizontal saccade offset**/
    if  (cmd->quantum[0][0] == 2.0) {

        if  (cmd->S_from2[0][ct_t][2] == 1.0)
            for (int i = 0; i < A[area].d_n; ++i)
                cmd->S_target[ct_t][i] += cmd->S_from1[0][ct_t][i] * scale_factor;
 
        if  (cmd->S_from2[0][ct_t][3] == 1.0)
            for (int i = 0; i < A[area].d_n; ++i)
                cmd->S_target[ct_t][i] -= cmd->S_from1[0][ct_t][i] * scale_factor;
    }

    return (DOUBLE)(0);
}



/**************************** total_converge_to **************************/
/*   */

DOUBLE total_converge_to (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    int inarea1 = cmd->n_from1[0];
    int inarea2 = cmd->n_from2[0];
    int area = cmd->area;

    if  (A[area].d_n != A[inarea1].d_n)
        fprintf (stderr, "\ntotal_converge_to wants all areas of same size!\n");
    if  (A[area].d_n != A[inarea2].d_n)
        fprintf (stderr, "\ntotal_converge_to wants all areas of same size!\n");
    if  ((cmd->quantum[0][0] < 0) || (cmd->quantum[0][0] > 3))
        fprintf (stderr, "\ntotal_converge_to allows occurrences 0-3 only!\n");

    static int firsttime[4] = {1,1,1,1};
    static DOUBLE *values[4];
    static int *refresh[4];
    static int counter[4] = {0,0,0,0};

    int occurrence = (int)(cmd->quantum[0][0]);

    if  (firsttime[occurrence]) {

        values[occurrence] = (DOUBLE *)calloc(A[area].d_n, sizeof (DOUBLE));
        refresh[occurrence] = (int *)calloc(A[area].d_n, sizeof (int));

        firsttime[occurrence] = 0;
    }

    for (int i = 0; i < A[area].d_n; ++i)
        if  (refresh[occurrence][i] == 1) {
            /**first time per unit always update with 0 (formerly with the new activations)**/
            values[occurrence][i] = 0.0; /* !! cmd->S_from1[0][ct_t][i]; !! */
            refresh[occurrence][i] = 0;
        } else {
            /**store values if decreasing and smallest so far**/
            if  (cmd->quantum[1][0] == -1) {
                if  ((cmd->S_from1[0][ct_t][i] < cmd->S_from2[0][ct_t][i]) && (cmd->S_from1[0][ct_t][i] < values[occurrence][i]))
                    values[occurrence][i] = cmd->S_from1[0][ct_t][i];
            }
            /**store values if increasing and largest so far**/
            if  (cmd->quantum[1][0] == 1) {
                if  ((cmd->S_from1[0][ct_t][i] > cmd->S_from2[0][ct_t][i]) && (cmd->S_from1[0][ct_t][i] > values[occurrence][i]))
                    values[occurrence][i] = cmd->S_from1[0][ct_t][i];
            }
        }

    if  (counter[occurrence] % (int)(cmd->quantum[2][0]) == 0) {
        for (int i = 0; i < A[area].d_n; ++i)
            refresh[occurrence][i] = 1;
        fprintf (stderr, "\ntotal_converge_to: resetting counter %d    ", occurrence);
    }

    for (int i = 0; i < A[area].d_n; ++i)
        cmd->S_target[ct_t][i] = values[occurrence][i];

    counter[occurrence] ++;

    return (DOUBLE)(0);
}



/************************** total_mean_left_min_right ************************/
/* S_target will all be the mean of all neurons activations (fixed time).    */

DOUBLE total_mean_left_min_right (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    double mean_left = 0.0, mean_right = 0.0;
    int ct_left = 0, ct_right = 0;
    int inarea = cmd->n_from1[0];
    int d_b = A[inarea].d_b;

    for (int X = 0; X < A[inarea].d_a; X++) {
        for (int Y = 0; Y < d_b; Y++) {
            if  (Y < d_b / 2) {
                mean_left += cmd->S_from1[0][ct_t][X * d_b + Y];
                ct_left ++;
            }
            if  (Y > d_b / 2) {
                mean_right += cmd->S_from1[0][ct_t][X * d_b + Y];
                ct_right ++;
            }
        }
    }

    mean_left  /= (double)(ct_left);
    mean_right /= (double)(ct_right);

    for (int i = 0; i < A[cmd->area].d_n; ++i)
        cmd->S_target[ct_t][i] = mean_left - mean_right;

    return (DOUBLE)(0);
}



/****************************** total_topo_gibbs_01 **************************/
/* inits act's with gaussian distributed random ON's                         */
/* over time, center remains constant, but different draws                   */
/* q[0][0] sparseness; if < 0 then exact negative number of ON units         */
/* q[0][1] saturation                                                        */
/* q[1][0] a-axis of gaussian (ori is random)                                */
/* q[1][1] b-axis of gaussian                                                */
/* q[2][0] mean number of gaussians according to poiss; -1: exactly 1        */
/* q[3][0] radius at corners to limit act-blobs (retry until near center)    */

/**aux function 1**/
DOUBLE gauss_distr () {

    DOUBLE dist, rand_1, rand_2;
    do  {
        rand_1 = -1.0 + 2.0 * drand48();
        rand_2 = -1.0 + 2.0 * drand48();
        dist   = rand_1*rand_1 + rand_2*rand_2;
    } while ((dist >= 1.0) || (dist == 0.0));

    dist = sqrt (-2.0 * log(dist) / dist);
    return (dist * rand_1);
}

/**aux function 2**/
int is_at_corner (int x_pos, int y_pos, int radius, int d_a, int d_b, int RGB) {
    int d_a_eff = d_a;

    if  (RGB == 3) {
        d_a_eff = d_a / 3;
        while (x_pos >= d_a_eff)
            x_pos -= d_a_eff;
    }

    /**upper left corner**/
    if  ((x_pos < radius) && (y_pos < radius))
        if  ( (x_pos - radius) * (x_pos - radius)
            + (y_pos - radius) * (y_pos - radius) > radius * radius)
            return 1;

    /**upper right corner**/
    if  ((x_pos < radius) && (d_b - 1 - y_pos < radius))
        if  ( (x_pos - radius) * (x_pos - radius)
            + (d_b - 1 - radius - y_pos) * (d_b - 1 - radius - y_pos) > radius * radius)
            return 1;

    /**lower left corner**/
    if  ((d_a_eff - 1 - x_pos < radius) && (y_pos < radius))
        if  ( (d_a_eff - 1 - radius - x_pos) * (d_a_eff - 1 - radius - x_pos)
            + (y_pos - radius) * (y_pos - radius) > radius * radius)
            return 1;

    /**lower right corner**/
    if  ((d_a_eff - 1 - x_pos < radius) && (d_b - 1 - y_pos < radius))
        if  ( (d_a_eff - 1 - radius - x_pos) * (d_a_eff - 1 - radius - x_pos)
            + (d_b - 1 - radius - y_pos) * (d_b - 1 - radius - y_pos) > radius * radius)
            return 1;

    return 0;
}


DOUBLE total_topo_gibbs_01 (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

  int anz_pix, i, ct_g, pos_x, pos_y;

  /**determine center and angle**/
  static double cm_x[32];
  static double cm_y[32];
  static double angle[32];
  static int    anz_gauss;
  const  int    area = cmd->area;

  if  (  (cmd->anz_quantums < 3)
      || (cmd->anz_quant[0] != 2) || (cmd->anz_quant[1] != 2))
    fprintf (stderr, "\nwrong use of total_topo_gibbs_01!\n");

  int radius = 0;   /**this is to disallow act-blob centers to be at corners**/
  if  (cmd->anz_quantums == 4)
      radius = (int)(cmd->quantum[3][0]);

  /* if  (ct_t == 0) */  /**originally the same CM was used in all relaxation steps**/
  {

      anz_gauss = (cmd->quantum[2][0] == -1.0) ? 1 : (int)(cmd->quantum[2][0]);
      if  (anz_gauss > 32) {
          anz_gauss = 32;
          fprintf (stderr, "\ntotal_topo_gibbs_01: gaussians limited to 32 ");
      }

      for (ct_g = 0; ct_g < anz_gauss; ++ct_g) {

          int isinside = 0;
          do {
              cm_x[ct_g]  = drand48() * (double)(A[area].d_a);
              cm_y[ct_g]  = drand48() * (double)(A[area].d_b);

              /**while act-blob center is at corner try again**/
              if  (! is_at_corner ((int)(cm_x[ct_g]), (int)(cm_y[ct_g]), radius, A[area].d_a, A[area].d_b, 1))
                  isinside = 1;

          } while (! isinside);

          angle[ct_g] = drand48() * 2.0 * M_PI;
      }
  }

  /**init act's zero**/
  for (i = 0; i < A[area].d_n; ++i)
      cmd->S_target[ct_t][i] = 0.0;

  for (ct_g = 0; ct_g < anz_gauss; ++ct_g) {

      /**determine number of ON's**/
      anz_pix = 0;
      for (i = 0; i < A[area].d_n; ++i)
          if  (drand48() < (1.0 / (1.0 + cmd->quantum[0][0])))
              anz_pix += 1;

      if  (cmd->quantum[0][0] < 0.0)
          anz_pix = (int)(-cmd->quantum[0][0]);


      /**set act's ON independently for all times**/
      /**kind of reflecting (no absorbing or periodic) boundary conditions**/
      for (i = 0; i < anz_pix; ++i) {

          /* int one_was_out = -1; */

          do  {
              double dist_long  = gauss_distr() * cmd->quantum[1][0];
              double dist_short = gauss_distr() * cmd->quantum[1][1];

              pos_x = (int)(cm_x[ct_g] + dist_long * cos (angle[ct_g]) + dist_short * sin (angle[ct_g]));
              pos_y = (int)(cm_y[ct_g] + dist_long * sin (angle[ct_g]) + dist_short * cos (angle[ct_g]));

              /* one_was_out += 1; */

          } while ((pos_x < 0) || (pos_x >= A[area].d_a) || (pos_y < 0) || (pos_y >= A[area].d_b));

          /*
          if  (one_was_out > 0)
              if  (drand48() > 0.5)           **the edges now get only a little more pixels ON !!! heuristic !!!** did'nt make much difference
                  i += 1;
          */

          cmd->S_target[ct_t][pos_x * A[area].d_b + pos_y] = cmd->quantum[0][1];     /**experience says: this should not be additive for topo-HM**/
      }
  }

    return (DOUBLE)(0);
}
