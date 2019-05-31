#include <stdlib.h>
#include <stdio.h>
#include <math.h>    /**only for sqrt in feed_l_covar**/
#include "iter.h"
#include "series.h"  /**only for COMMAND**/
#include "local.h"   /**for feed_mean_01, (etc.)**/


/****************************** feed_in_W ************************************/
/* Function returns inner activation of neuron (5) ct_n  in (3) cmd->area    */
/* at time (4) ct_t  dependent on A[].R/T at time ct_t + cmd->offset.        */
/* Initialization with -Theta if ch_Theta.                                   */

double feed_in_W (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int inarea, ct_l, ct_in;
    double *S, *W;

    const int area = cmd->area;
    double val = z->ch_Theta ?  -1.0 * A[area].Theta[ct_n] : 0.0;
    const int oldtime = ct_t + (int)(cmd->quantum[0][0]) < 0 ? 0
                      : ct_t + (int)(cmd->quantum[0][0]);

ERR("\nfeed_in_W");

    for (ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {         /**list of areas**/

        inarea = cmd->n_from1[ct_l];       /**inarea index directly from cmd**/

        S = cmd->S_from1[ct_l][oldtime];                       /**set source**/

        W = A[area].W[inarea][ct_n];                      /**maybe faster so**/

        for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)

            val += W[ct_in] * S[ct_in];
    }

ERR("feed_W_end");

    return (val);
}

/****************************** feed_in_V ************************************/
/* Like feed_in_W but input weights are V.                                   */

double feed_in_V (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int inarea, ct_l, ct_in;
    double *S, *W;

    const int area = cmd->area;
    double val = z->ch_Theta ?  -1.0 * A[area].Theta[ct_n] : 0.0;
    const int oldtime = ct_t + (int)(cmd->quantum[0][0]) < 0 ? 0
                      : ct_t + (int)(cmd->quantum[0][0]);

ERR("\nfeed_in_V");

    for (ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {         /**list of areas**/

        inarea = cmd->n_from1[ct_l];       /**inarea index directly from cmd**/

        S = cmd->S_from1[ct_l][oldtime];                       /**set source**/

        W = A[area].V[inarea][ct_n];                      /**maybe faster so**/

        for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)

            val += W[ct_in] * S[ct_in];
    }

ERR("feed_V_end");

    return (val);
}

/****************************** feed_in_X ************************************/
/* Like feed_in_W but input weights are X.                                   */

double feed_in_X (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int inarea, ct_l, ct_in;
    double *S, *W;

    const int area = cmd->area;
    double val = z->ch_Theta ?  -1.0 * A[area].Theta[ct_n] : 0.0;
    const int oldtime = ct_t + (int)(cmd->quantum[0][0]) < 0 ? 0
                      : ct_t + (int)(cmd->quantum[0][0]);

ERR("\nfeed_in_X");

    for (ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {         /**list of areas**/

        inarea = cmd->n_from1[ct_l];       /**inarea index directly from cmd**/

        S = cmd->S_from1[ct_l][oldtime];                       /**set source**/

        W = A[area].X[inarea][ct_n];                      /**maybe faster so**/

        for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)

            val += W[ct_in] * S[ct_in];
    }

ERR("feed_X_end");

    return (val);
}

/****************************** feed_in_Y ************************************/
/* Like feed_in_W but input weights are Y.                                   */

double feed_in_Y (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int inarea, ct_l, ct_in;
    double *S, *W;

    const int area = cmd->area;
    double val = z->ch_Theta ?  -1.0 * A[area].Theta[ct_n] : 0.0;
    const int oldtime = ct_t + (int)(cmd->quantum[0][0]) < 0 ? 0
                      : ct_t + (int)(cmd->quantum[0][0]);

ERR("\nfeed_in_Y");

    for (ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {         /**list of areas**/

        inarea = cmd->n_from1[ct_l];       /**inarea index directly from cmd**/

        S = cmd->S_from1[ct_l][oldtime];                       /**set source**/

        W = A[area].Y[inarea][ct_n];                      /**maybe faster so**/

        for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)

            val += W[ct_in] * S[ct_in];
    }

ERR("feed_Y_end");

    return (val);
}


/****************************** feed_in_V_pm *********************************/
/* From feed_in_V; no mixing combination of + and -. Only the strongest wins.*/
/* For positive&negative data in a Helmholtz machine.                        */

double feed_in_V_pm (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int inarea, ct_l, ct_in;
    double *S, *W;

    const int area = cmd->area;
    double val_pos = z->ch_Theta ?  -1.0 * A[area].Theta[ct_n] : 0.0;
    double val_neg = z->ch_Theta ?  -1.0 * A[area].Theta[ct_n] : 0.0;
    const int oldtime = ct_t + (int)(cmd->quantum[0][0]) < 0 ? 0
                      : ct_t + (int)(cmd->quantum[0][0]);

ERR("\nfeed_in_V_pm");

    for (ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {         /**list of areas**/

        inarea = cmd->n_from1[ct_l];       /**inarea index directly from cmd**/

        S = cmd->S_from1[ct_l][oldtime];                       /**set source**/

        W = A[area].V[inarea][ct_n];                      /**maybe faster so**/

        for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)

            if  (W[ct_in] >= 0)
                val_pos += W[ct_in] * S[ct_in];
            else
                val_neg += W[ct_in] * S[ct_in];
    }

ERR("feed_V_pm_end");

    return (fabs(val_pos) > fabs(val_neg) ? val_pos : val_neg);
}


/****************************** feed_cheap ***********************************/
/* Sum of input areas1 modulated by sum of input areas2. One value z->cheap. */

double feed_cheap (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    double val = 0.0, val2 = 0.0, *W, *V, *S;
    int ct_l, inarea, ct_in;

    double retval = z->ch_Theta ? -1.0 * A[cmd->area].Theta[ct_n] : 0.0;
    const int oldtime = ct_t + (int)(cmd->quantum[0][0]) < 0 ? 0
                      : ct_t + (int)(cmd->quantum[0][0]);

    for (ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        inarea = cmd->n_from1[ct_l];

        S = cmd->S_from1[ct_l][oldtime];

        W = A[cmd->area].W[inarea][ct_n];

        for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)

            val += W[ct_in] * S[ct_in];
    }

    for (ct_l = 0; ct_l < cmd->anz_from2; ++ct_l) {       /**list of areas 2**/

        inarea = cmd->n_from2[ct_l];                               /**area 2**/

        S = cmd->S_from2[ct_l][oldtime];                     /**set source 2**/

        V = A[cmd->area].V[inarea][ct_n];                       /**weights 2**/

        for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)

            val2 += V[ct_in] * S[ct_in];                            /**val 2**/
    }

    return (retval + val - (val - val * val2) * z->cheap);
}


/****************************** feed_moderate ********************************/
/* Every input area1 modulated by every input area2. Matrix of z->full.      */
/* Probably don't need this.                                                 */

double feed_moderate (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    static int firsttime = 1;
    static double *val, *val2;
    double *W, *V, *S;
    int inarea, ct_l, ct_l2, ct_in;

    double retval = z->ch_Theta ? -1.0 * A[cmd->area].Theta[ct_n] : 0.0;
    const int oldtime = ct_t + (int)(cmd->quantum[0][0]) < 0 ? 0
                      : ct_t + (int)(cmd->quantum[0][0]);

    if  (firsttime) {
        val  = (double *)malloc (cmd->anz_from1 * sizeof (double));
        val2 = (double *)malloc (cmd->anz_from2 * sizeof (double));
        firsttime = 0;
    }


    for (ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        inarea = cmd->n_from1[ct_l];

        S = cmd->S_from1[ct_l][oldtime];

        W = A[cmd->area].W[inarea][ct_n];

        val[ct_l] = 0.0;

        for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)

            val[ct_l] += W[ct_in] * S[ct_in];
    }


    for (ct_l = 0; ct_l < cmd->anz_from2; ++ct_l) {       /**list of areas 2**/

        inarea = cmd->n_from2[ct_l];                            /**area 2**/

        S = cmd->S_from2[ct_l][oldtime];                     /**set source 2**/

        V = A[cmd->area].V[inarea][ct_n];                       /**weights 2**/

        val2[ct_l] = 0.0;

        for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)

            val2[ct_l] += V[ct_in] * S[ct_in];                      /**val 2**/
    }


    for (ct_l = 0; ct_l < cmd->anz_from1; ++ct_l)

        for (ct_l2 = 0; ct_l2 < cmd->anz_from2; ++ct_l2)

            retval += val[ct_l] - (val[ct_l] - val[ct_l] * val2[ct_l2])
                                  * z->full[ct_l][ct_l2];

    return (retval);
}


/****************************** feed_full ************************************/
/* Every input neuron1 modulated by every input neuron2. Matrix of z->full.  */

double feed_full (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    double **U, *S, *T;
    int inarea, inarea2, ct_l, ct_l2, ct_in, ct_in2;

    double retval = z->ch_Theta ? -1.0 * A[cmd->area].Theta[ct_n] : 0.0;
    const int oldtime = ct_t + (int)(cmd->quantum[0][0]) < 0 ? 0
                      : ct_t + (int)(cmd->quantum[0][0]);

    for (ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        for (ct_l2 = 0; ct_l2 < cmd->anz_from2; ++ct_l2) {/**list of areas 2**/

            inarea = cmd->n_from1[ct_l];
            inarea2 = cmd->n_from2[ct_l2];                         /**area 2**/

            S = cmd->S_from1[ct_l][oldtime];

            T = cmd->S_from2[ct_l2][oldtime];                /**set source 2**/

            U = A[cmd->area].U[inarea2][inarea][ct_n];

            if  (z->full[ct_l][ct_l2] == 1.0) {

                for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)

                    for (ct_in2 = 0; ct_in2 < A[inarea2].d_r; ++ct_in2)

                        /**standard: symmetry between both inputs**/
                        retval += U[ct_in2][ct_in] * S[ct_in] * T[ct_in2];

            } else {

                if  (z->full[ct_l][ct_l2] == 0.0)

                    for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)

                        for (ct_in2 = 0; ct_in2 < A[inarea2].d_r; ++ct_in2)

                            /**not modulated but 2nd order i.e. more weights**/
                            retval += U[ct_in2][ct_in] * S[ct_in];

                else

                    for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)

                        for (ct_in2 = 0; ct_in2 < A[inarea2].d_r; ++ct_in2)

                            /**unsymmetry depends on area parameters**/
                            retval += U[ct_in2][ct_in]
                                    * ( S[ct_in]
                                      - ( S[ct_in] - S[ct_in] * T[ct_in2] )
                                        * z->full[ct_l][ct_l2] );
            }

        }
    }

    return (retval);
}



/****************************** feed_back_W **********************************/
/* Analogeous to feed_in_W but instead  sum_j w_ij r_j  does  sum_i w_ij u_i */
/* May make obsolete symmetric weight_copy.            No Theta's considered.*/

/* No feed_back_V/X/Y yet. First test this. Then "copy"!                     */

double feed_back_W (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int inarea, ct_l, ct_in;
    double *S, **W;                              /**change1 ** instead of * **/

    const int area = cmd->area;
    double val = 0.0;
    const int oldtime = ct_t + (int)(cmd->quantum[0][0]) < 0 ? 0
                      : ct_t + (int)(cmd->quantum[0][0]);

ERR("\nfeed_back");

    for (ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        inarea = cmd->n_from1[ct_l];

        S = cmd->S_from1[ct_l][oldtime];

        W = A[inarea].W[area];                                    /**change2**/
	   /**^^   !    ^**/

        for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)

            val += W[ct_in][ct_n] * S[ct_in];                     /**change3**/
    }

ERR("feed_back_end");

    return (val);
}



/****************************** feed_euclid_W ********************************/
/* Function returns Euclidian distance of neuron (5) ct_n's weight vector    */
/* in (3) cmd->area to its input (data) vector at time (4) ct_t              */
/* dependent on A[].R/T at time ct_t + cmd->offset.                          */
/* Initialization with -Theta if ch_Theta.                                   */

double feed_euclid_W (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int inarea, ct_l, ct_in;
    double *S, *W;

    const int area = cmd->area;
    double val = z->ch_Theta ?  -1.0 * A[area].Theta[ct_n] : 0.0;
    const int oldtime = ct_t + (int)(cmd->quantum[0][0]) < 0 ? 0
                      : ct_t + (int)(cmd->quantum[0][0]);

ERR("\nfeed_euclid_W");

    for (ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {         /**list of areas**/

        inarea = cmd->n_from1[ct_l];       /**inarea index directly from cmd**/

        S = cmd->S_from1[ct_l][oldtime];                       /**set source**/

        W = A[area].W[inarea][ct_n];                      /**maybe faster so**/

        for (ct_in = 0; ct_in < A[inarea].d_r; ++ct_in)

            val += (S[ct_in] - W[ct_in]) * (S[ct_in] - W[ct_in]);
    }

ERR("feed_euclid_W_end");

    return (sqrt (val));
}



/****************************** feed_l_add ***********************************/
/* Updates old source with eps_r-correction of new estimation.               */
/* One source allowed only ("here" or another area of the same size).        */
/* Could be local (_l_) if there was no offset.                              */

double feed_l_add (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int t_old = (cmd->anz_quantums == 2) ? 0      /**functions as dummy here**/
              : fprintf (stderr, "\nfeed_l_add wants 2 quantums!\n");
    const double eps_r = cmd->quantum[1][0];

    t_old    = ct_t + (int)(cmd->quantum[0][0]);

    if  (t_old < 0)
        return ( cmd->S_target[  0  ][ct_n]);      /*was from1[0]-->target*/
    else
        return ( cmd->S_target[t_old][ct_n]     /*first was from1[0]-->target*/
               + cmd->S_from1[0][ct_t ][ct_n] * eps_r);
}



/****************************** feed_l_replace *******************************/
/* Like feed_l_add, but NOT solely additive.                                 */

double feed_l_replace (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int t_old = (cmd->anz_quantums == 2) ? 0      /**functions as dummy here**/
              : fprintf (stderr, "\nfeed_l_replace wants 2 quantums!\n");
    const double eps_r = cmd->quantum[1][0];

    t_old = ct_t + (int)(cmd->quantum[0][0]);

    if  (t_old < 0)
        return ( cmd->S_target[  0  ][ct_n]);      /*was from1[0]-->target*/
    else
        return ( cmd->S_target[t_old][ct_n] * (1.0 - eps_r)      /**!**/
               + cmd->S_from1[0][ct_t][ct_n] * eps_r);
}



/****************************** feed_copy ************************************/
/* Copys old source to new time.                                             */
/* One source allowed only ("here" or another area of the same size).        */
/* Previously local (feed_l_copy).                                           */

double feed_copy (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    const int t_old = ct_t + (int)(cmd->quantum[0][0]);

    if  (t_old < 0)
        return ( cmd->S_from1[0][  0  ][ct_n]);
    else
        return ( cmd->S_from1[0][t_old][ct_n]);
}



/****************************** feed_copy_limited ****************************/
/* Copys, but only the first q[0][0] neurons.                                */
/* One source (prob. another area). q[0][0] must NOT be larger than this d_r.*/

double feed_copy_limited (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    if  (ct_n < cmd->quantum[0][0])
        return (cmd->S_from1[0][ct_t][ct_n]);
    else
        return (0.0);
}



/****************************** feed_l_theta *********************************/
/* Returns A[area].Theta[ct_n].                                              */
/* Could be local (_l_) if it would not ask for Theta.                       */

double feed_l_theta (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    return (A[cmd->area].Theta[ct_n]);
}



/****************************** feed_l_difu **********************************/
/* Activity transfer from 4 neighbors, scaled by 0 <= q[0][0] <= 0.8!        */
/* Toroidal boundary conditions.                                             */
/* Warning: cmd->b_target must NOT be cmd->b_from1 !!                        */
/* Warning: use extra stay to use current input !!                           */

double feed_l_difu (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    const double dif  = cmd->quantum[0][0];
    const double *mid = cmd->S_from1[0][ct_t];

    const int d_r   = z->d_a * z->d_b;
    const int right = ((ct_n + 1) % z->d_b == 0)
                                           ? ct_n - z->d_b + 1 : ct_n + 1;
    const int left  = ((ct_n + z->d_b - 1) % z->d_b == z->d_b - 1) 
                                           ? ct_n + z->d_b - 1 : ct_n - 1;

    return ( mid[ct_n] * (1.0 - dif)
           + ( mid[right]
             + mid[(ct_n + z->d_b) % d_r]         /**under**/
             + mid[left]
             + mid[(ct_n + d_r - z->d_b) % d_r]   /**upper**/
             ) * dif * 0.25);
}


/****************************** feed_l_rand_from *****************************/
/* Set to random values taken from a given set and with given probabilities. */
/* No input. Local. Parameters like:                                         */
/* q[0][.]: -1+-0.5+0.5+1         the states                                 */
/* q[1][.]: 0.25+0.25+0.25+0.25   their probabilities (need not sum to one). */

double feed_l_rand_from (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    double part_sum = 0.0;
    double choose;
    int    i;

    for (i = 0; i < cmd->anz_quant[1]; ++i)
        part_sum += cmd->quantum[1][i];

    choose = drand48 () * part_sum;

    i = -1;
    do {
        i += 1;
        choose -= cmd->quantum[1][i];

    } while (choose > 0.0);

    return (cmd->quantum[0][i]);
}





/****************************** get_l_covar **********************************/
/* Get "diagonalized correlation matrix" for act based on 1000 iterations(!!)*/
/* One source allowed only ("here").                                         */
/* Could be local (_l_) if there was no average over time.                   */

double feed_l_covar (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int ilen = 1000; /**!!! cannot get from x->ilen ; how about z->ilen? **/

    static int listlen = 0;
    int sel = -1;
    int i;
    struct liste { int area ; char letter ; int zaehler ; int maxzaehl ;
                   double *values ; double *oldvals ; };
    static struct liste *list;
    int area = cmd->area;

    /** use area and target value to distinguish different computations **/
    for (i = 0; i < listlen; ++i)
        if  ((list[i].area == area) && (list[i].letter == cmd->b_target))
            sel = i;

    /**allocate and initialize**/
    if  (sel == -1) {
        list = (struct liste *)realloc (list,
                                        (listlen + 1) * sizeof (struct liste));
        sel = listlen;
        list[sel].area    = area;
        list[sel].letter  = cmd->b_target;
        list[sel].zaehler = 0;
        list[sel].maxzaehl= ilen * A[area].d_r * (int)(1.0/cmd->moment + 0.01);
        list[sel].values  = (double *)calloc(A[area].d_r, sizeof(double));
        list[sel].oldvals = (double *)calloc(A[area].d_r, sizeof(double));
        for (i = 0; i < A[area].d_r; ++i) {
            list[sel].values[i]  = 0.0;                     /*init with 0*/
            list[sel].oldvals[i] = 1.0;                     /*init with 1*/
        }
        fprintf (stderr, "\n adding %d %c to the %d-iteration covar list  ",
                             list[sel].area, list[sel].letter, ilen);
        listlen += 1;
    }


    /**compute**/
    list[sel].values[ct_n] += cmd->S_from1[0][ct_t][ct_n]
                            * cmd->S_from1[0][ct_t][ct_n];


    /**avarage, also over time, and initialize**/
    if  ((list[sel].zaehler % list[sel].maxzaehl < A[area].d_r) &&
         (list[sel].zaehler >= A[area].d_r)) {

        /* fprintf (stderr, "%d", area); */

        list[sel].oldvals[ct_n] = list[sel].values[ct_n]
                       / (double)(list[sel].maxzaehl) * (double)(A[area].d_r);

        list[sel].values[ct_n] = 0.0;
    }

    if  (list[sel].zaehler % list[sel].maxzaehl == A[area].d_r)
        list[sel].zaehler = A[area].d_r;

    list[sel].zaehler += 1;

    return (list[sel].oldvals[ct_n]);            /** =1 at first time**/
}





/****************************** feed_l_hack **********************************/
/* q[0][.] sparseness parameters                                             */
/* q[1][.] spontaneous activity parameters                                   */
/* q[2][0] transf-fkt:=0 oldmodel; =1 mean; =2 gibbs; =3 rand; =4; add_rand; */
/* q[3][0] =0 neg&pos values; =else only positive output                     */
/* q[4][0] beta                                                              */
/* q[5][0] sat_val                                                           */

double feed_l_hack (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    const double val1 = cmd->S_from1[0][ct_t][ct_n];
    const double sparseness = (ct_n < z->d_a * z->d_b / 2)            /**!!!**/
                            ? cmd->quantum[0][0] : cmd->quantum[0][1];
    const double spont_act  = (ct_n < z->d_a * z->d_b / 2)            /**!!!**/
                            ? cmd->quantum[1][0] : cmd->quantum[1][1];
    const double beta    = cmd->quantum[4][0];
    const double sat_val = cmd->quantum[5][0];

    /**linear gradient, overwrites upper assignment
    sparseness = cmd->quantum[0][0] + (cmd->quantum[0][1] - cmd->quantum[0][0])
               * (double)(int)(ct_n / z->d_b) / (double)(z->d_a - 1);
    spont_act  = cmd->quantum[1][0] + (cmd->quantum[1][1] - cmd->quantum[1][0])
               * (double)(int)(ct_n / z->d_b) / (double)(z->d_a - 1);
    **/

    if  (cmd->anz_quantums != 6)
        fprintf (stderr, "\n6 parameters for feed_l_hack, please!\n");
    if  (spont_act != 0.0)
        if  (cmd->quantum[3][0] == 0)
            fprintf (stderr, "\nfeed_l_hack: no spont_act if pos&neg act!\n");
    if  (sat_val == 0.0)
        fprintf (stderr, "\nfeed_l_hack: sat_val should not be zero\n");

    /**oldmodel**/
    if  (cmd->quantum[2][0] == 0) {
        if  (cmd->quantum[3][0] == 0) {
            return (val1 - sparseness * 2.0 * val1 / (1.0 + val1 * val1));
        } else {
            if  (val1 > 0.0)
                return (val1 - sparseness * 2.0 * val1 / (1.0 + val1 * val1));
            else
                return (0.0);
        }
    }

    /**mean**/
    if  (cmd->quantum[2][0] == 1) {

        double exp_h = exp (beta * val1);

        if  (cmd->quantum[3][0] == 0) {                           /**neg&pos**/
            double exp_min_h = exp (-beta * val1);

            return (sat_val * (exp_h - exp_min_h)
                           / (exp_h + sparseness + exp_min_h));
        } else {                                                 /**pos only**/
/***old version with unclean spontaneous activity
            if  (spont_act == 0.0) {
                return (sat_val * exp_h / (exp_h + sparseness));
            } else {
                double f1 = exp_h / (exp_h + sparseness);
                return (sat_val * (f1 + spont_act * (1.0 - f1)));
            }
***/
            return ((exp_h + spont_act) / (exp_h + spont_act + sparseness));
        }
    }

    /**gibbs**/
    if  (cmd->quantum[2][0] == 2) {
        double exp_h     = exp (beta * val1);
        double rd = drand48 ();

        if  (cmd->quantum[3][0] == 0) {                           /**neg&pos**/
            double exp_min_h = exp (-beta * val1);
            double part_sum  = exp_h + sparseness + exp_min_h;

            if  (rd < (exp_min_h / part_sum)) {
                return (-sat_val);
            } else {
                if  (rd > (1.0 - exp_h / part_sum))
                    return (sat_val);
                else
                    return (0.0);
            }
        } else {                                                 /**pos only**/
/***
            double part_sum  = exp_h + sparseness;
            if  (rd < exp_h / part_sum) {
                return (sat_val);
            } else {                       **add spont.act. to 0-values only**
                if  (spont_act != 0.0) {
                    if  (drand48 () < spont_act)
                        return (sat_val);
                    else
                        return (0.0);
                } else {
                    return (0.0);
                }
            }
***/
            double part_sum = exp_h + spont_act + sparseness;
            if  (rd < (exp_h + spont_act) / part_sum)
                return (sat_val);
            else
                return (0.0);
        }
    }

    /**init rand_gibbs**/
    if  (cmd->quantum[2][0] == 3) {
        if  (cmd->quantum[3][0] == 0)
            return (drand48() < (1.0 / (1.0 + sparseness))
                    ? (drand48() < 0.5 ? sat_val : -sat_val)
                    : 0.0);
        else
            return (drand48() < (1.0 / (1.0 + sparseness))
                    ? sat_val : 0.0);
    }

    fprintf (stderr, "\nfeed_l_hack should not come to here\n");
    return (-9.9);
}


/****************************** feed_const_hack ******************************/
/* Like local_const but takes two (different) values in q[0][0] and q[1][0]. */

double feed_const_hack (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    if  (cmd->anz_quantums != 2)
        fprintf (stderr, "\n2 parameters for feed_const_hack, please!\n");

    if  (ct_n < z->d_a * z->d_b / 2)
        return (cmd->quantum[0][0]);
    else
        return (cmd->quantum[1][0]);
}





/******************************* feed_mean_01 ********************************/
/* Like corresponding local function, but uses activations as parameters.    */
/* par[0](kurt) for sparseness.                                              */
/* par[1](beta) is inverse Temperature.                                      */
/* par[2](aux2) resting act (degeneracy of spontaneous ON-state!)            */
/* par[3](bias) negative threshold, comparable to local_zhang (new!).        */
/* par[4](sat) for max amplitude.                                            */

double feed_mean_01 (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    double par[5];
    double val1 = cmd->S_from1[0][ct_t][ct_n];

    if  (cmd->anz_from2 != 5)  fprintf (stderr, "\nfeed_mean_01 error\n");

    par[0] = cmd->S_from2[0][ct_t][ct_n];
    par[1] = cmd->S_from2[1][ct_t][ct_n];
    par[2] = cmd->S_from2[2][ct_t][ct_n];
    par[3] = cmd->S_from2[3][ct_t][ct_n];
    par[4] = cmd->S_from2[4][ct_t][ct_n];

    return (local_mean_01 (par, val1, 0.0));
}

/******************************* feed_mean_01_inv ****************************/
/* Inverse.                                                                  */

double feed_mean_01_inv (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    double par[5];
    double val1 = cmd->S_from1[0][ct_t][ct_n];

    if  (cmd->anz_from2 != 5)  fprintf (stderr, "\nfeed_mean_01_inv error\n");

    par[0] = cmd->S_from2[0][ct_t][ct_n];
    par[1] = cmd->S_from2[1][ct_t][ct_n];
    par[2] = cmd->S_from2[2][ct_t][ct_n];
    par[3] = cmd->S_from2[3][ct_t][ct_n];
    par[4] = cmd->S_from2[4][ct_t][ct_n];

    return (local_mean_01_inv (par, val1, 0.0));
}

/******************************* feed_mean_01_diff ***************************/
/* Derivative w.r.t. input.                                                  */

double feed_mean_01_diff (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    double par[5];
    double val1 = cmd->S_from1[0][ct_t][ct_n];

    if  (cmd->anz_from2 != 5)  fprintf (stderr, "\nfeed_mean_01_diff error\n");

    par[0] = cmd->S_from2[0][ct_t][ct_n];
    par[1] = cmd->S_from2[1][ct_t][ct_n];
    par[2] = cmd->S_from2[2][ct_t][ct_n];
    par[3] = cmd->S_from2[3][ct_t][ct_n];
    par[4] = cmd->S_from2[4][ct_t][ct_n];

    return (local_mean_01_diff (par, val1, 0.0));
}

/******************************* feed_mean_01_d0 *****************************/
/* Derivative w.r.t. parameter 0, kurt.                                      */

double feed_mean_01_d0 (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    double par[5];
    double val1 = cmd->S_from1[0][ct_t][ct_n];

    if  (cmd->anz_from2 != 5)  fprintf (stderr, "\nfeed_mean_01_d0 error\n");

    par[0] = cmd->S_from2[0][ct_t][ct_n];
    par[1] = cmd->S_from2[1][ct_t][ct_n];
    par[2] = cmd->S_from2[2][ct_t][ct_n];
    par[3] = cmd->S_from2[3][ct_t][ct_n];
    par[4] = cmd->S_from2[4][ct_t][ct_n];

    return (local_mean_01_d0 (par, val1, 0.0));
}

/******************************* feed_mean_01_d1 *****************************/
/* Derivative w.r.t. parameter 1, beta.                                      */

double feed_mean_01_d1 (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    double par[5];
    double val1 = cmd->S_from1[0][ct_t][ct_n];

    if  (cmd->anz_from2 != 5)  fprintf (stderr, "\nfeed_mean_01_d1 error\n");

    par[0] = cmd->S_from2[0][ct_t][ct_n];
    par[1] = cmd->S_from2[1][ct_t][ct_n];
    par[2] = cmd->S_from2[2][ct_t][ct_n];
    par[3] = cmd->S_from2[3][ct_t][ct_n];
    par[4] = cmd->S_from2[4][ct_t][ct_n];

    return (local_mean_01_d1 (par, val1, 0.0));
}

/******************************* feed_mean_01_d2 *****************************/
/* Derivative w.r.t. parameter 2, rest_act.                                  */

double feed_mean_01_d2 (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    double par[5];
    double val1 = cmd->S_from1[0][ct_t][ct_n];

    if  (cmd->anz_from2 != 5)  fprintf (stderr, "\nfeed_mean_01_d2 error\n");

    par[0] = cmd->S_from2[0][ct_t][ct_n];
    par[1] = cmd->S_from2[1][ct_t][ct_n];
    par[2] = cmd->S_from2[2][ct_t][ct_n];
    par[3] = cmd->S_from2[3][ct_t][ct_n];
    par[4] = cmd->S_from2[4][ct_t][ct_n];

    return (local_mean_01_d2 (par, val1, 0.0));
}

/******************************* feed_mean_01_d3 *****************************/
/* Derivative w.r.t. parameter 3, bias.                                      */

double feed_mean_01_d3 (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    double par[5];
    double val1 = cmd->S_from1[0][ct_t][ct_n];

    if  (cmd->anz_from2 != 5)  fprintf (stderr, "\nfeed_mean_01_d3 error\n");

    par[0] = cmd->S_from2[0][ct_t][ct_n];
    par[1] = cmd->S_from2[1][ct_t][ct_n];
    par[2] = cmd->S_from2[2][ct_t][ct_n];
    par[3] = cmd->S_from2[3][ct_t][ct_n];
    par[4] = cmd->S_from2[4][ct_t][ct_n];

    return (local_mean_01_d3 (par, val1, 0.0));
}

/******************************* feed_mean_01_d4 *****************************/
/* Derivative w.r.t. parameter 4, sat.                                       */

double feed_mean_01_d4 (AGENT *z, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    double par[5];
    double val1 = cmd->S_from1[0][ct_t][ct_n];

    if  (cmd->anz_from2 != 5)  fprintf (stderr, "\nfeed_mean_01_d4 error\n");

    par[0] = cmd->S_from2[0][ct_t][ct_n];
    par[1] = cmd->S_from2[1][ct_t][ct_n];
    par[2] = cmd->S_from2[2][ct_t][ct_n];
    par[3] = cmd->S_from2[3][ct_t][ct_n];
    par[4] = cmd->S_from2[4][ct_t][ct_n];

    return (local_mean_01_d4 (par, val1, 0.0));
}
