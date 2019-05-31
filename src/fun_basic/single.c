#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../kernel/coco.h"
#include "../kernel/series.h"  /**only for COMMAND**/
#include "local.h"             /**only for single_circ_gauss**/


/****************************** single_ptr2act *******************************/
/* Target act's will be set to cmd-pointer's float_val.                      */

DOUBLE single_ptr2act (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    return (DOUBLE)(cmd->pointers[0]->float_val);
}

/****************************** set_ptr_int_val ******************************/
/* Can be used in all contexts ("nevermind" status).                         */
/* Later move all ptr-related functions to an extra file.                    */

DOUBLE set_ptr_int_val (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    cmd->pointers[0]->int_val = (int)(cmd->quantum[0][0]);

    return (DOUBLE)(0);
}

/******************************* set_ptr_val *********************************/
/* Set int_val and float_val of the pointer to S_from1.                      */ 

DOUBLE set_ptr_val (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    cmd->pointers[0]->int_val   = (int)(cmd->S_from1[0][ct_t][ct_n]);
    cmd->pointers[0]->float_val = cmd->S_from1[0][ct_t][ct_n];

    fprintf (stderr, "\nset_ptr_val to %f   ", cmd->S_from1[0][ct_t][ct_n]);

    return (cmd->S_from1[0][ct_t][ct_n]);
}


/****************************** single_copy **********************************/
/* Copys old source to new time.                                             */
/* One source allowed only ("here" or another area of the same size).        */
/* Previously feed_l_copy or feed_copy.                                      */

DOUBLE single_copy (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    const int t_old = ct_t + (int)(cmd->quantum[0][0]);

    if  (t_old < 0)
        return (cmd->S_from1[0][  0  ][ct_n]);
    else
        return (cmd->S_from1[0][t_old][ct_n]);
}


/****************************** single_copy_limited **************************/
/* Copys, but only the first q[0][0] neurons.                                */
/* One source (prob. another area). q[0][0] must NOT be larger than this d_r.*/
/* Has only special functionality in hacked version for saccade learning ... */

double single_copy_limited (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    if  (ct_n < cmd->quantum[0][0])
        return (cmd->S_from1[0][ct_t][ct_n]);
    else
        return (0.0);
}



/****************************** single_mean_back *****************************/
/* Averages acts between rlen time q[0][0] and q[1][0] and returns mean.     */

DOUBLE single_mean_back (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    if  (ct_t != 0)
        fprintf (stderr, "\nsingle_mean_back wants to operate only at t=0!\n");
    if  (cmd->anz_quantums != 2)
        fprintf (stderr, "\nsingle_mean_back wants begin and end\n");

    int begin = (int)(cmd->quantum[0][0]);
    int end   = (int)(cmd->quantum[1][0]);

    if  (end <= begin)
        fprintf (stderr, "\nsingle_mean_back wants end > begin\n");

    DOUBLE average = 0.0;
    for (int t = begin; t < end; ++t)
        average += cmd->S_from1[0][t][ct_n];

    return (average / (DOUBLE)(end - begin));
}



/******************************* single_circ_gauss ***************************/
/* Uses local_circ_gauss, but shift the mean by mu = 2PI / rlen * ct_t.      */

DOUBLE single_circ_gauss (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    DOUBLE params[2];

    params[0] = 2.0 * M_PI * (DOUBLE)ct_t / cmd->quantum[0][0];
    params[1] = cmd->quantum[0][1];

    return local_circ_gauss (params, cmd->S_from1[0][ct_t][ct_n], 0.0);
}



/****************************** single_l_add *********************************/
/* Updates old source with eps_r-correction of new estimation.               */
/* One source allowed only ("here" or another area of the same size).        */
/* Could be local (_l_) if there was no offset.                              */

DOUBLE single_l_add (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int t_old = (cmd->anz_quantums == 2) ? 0      /**functions as dummy here**/
              : fprintf (stderr, "\nsingle_l_add wants 2 quantums!\n");
    const DOUBLE eps_r = cmd->quantum[1][0];

    t_old    = ct_t + (int)(cmd->quantum[0][0]);

    if  (t_old < 0)
        return ( cmd->S_target[  0  ][ct_n]);      /*was from1[0]-->target*/
    else
        return ( cmd->S_target[t_old][ct_n]     /*first was from1[0]-->target*/
               + cmd->S_from1[0][ct_t][ct_n] * eps_r);
}



/****************************** feed_l_replace *******************************/
/* Like feed_l_add, but NOT solely additive.                                 */

DOUBLE feed_l_replace (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int t_old = (cmd->anz_quantums == 2) ? 0      /**functions as dummy here**/
              : fprintf (stderr, "\nfeed_l_replace wants 2 quantums!\n");
    const DOUBLE eps_r = cmd->quantum[1][0];

    t_old = ct_t + (int)(cmd->quantum[0][0]);

    if  (t_old < 0)
        return ( cmd->S_target[  0  ][ct_n]);      /*was from1[0]-->target*/
    else
        return ( cmd->S_target[t_old][ct_n] * (1.0 - eps_r)      /**!**/
               + cmd->S_from1[0][ct_t][ct_n] * eps_r);
}



/****************************** feed_l_rand_from *****************************/
/* Set to random values taken from a given set and with given probabilities. */
/* No input. Local. Parameters like:                                         */
/* q[0][.]: -1+-0.5+0.5+1         the states                                 */
/* q[1][.]: 0.25+0.25+0.25+0.25   their probabilities (need not sum to one). */

DOUBLE feed_l_rand_from (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    DOUBLE part_sum = 0.0;
    DOUBLE choose;
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

DOUBLE feed_l_covar (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int ilen = 1000; /**!!! cannot get from g->ilen ; how about z->ilen? **/

    static int listlen = 0;
    int sel = -1;
    int i;
    struct liste { int area ; char letter ; int zaehler ; int maxzaehl ;
                   DOUBLE *values ; DOUBLE *oldvals ; };
    static struct liste *list;
    int area = cmd->area;

    /** use area and target value to distinguish different computations **/
    for (i = 0; i < listlen; ++i)
        if  ((list[i].area == area) && (list[i].letter == cmd->ch_target))
            sel = i;

    /**allocate and initialize**/
    if  (sel == -1) {
        list = (struct liste *)realloc (list,
                                        (listlen + 1) * sizeof (struct liste));
        sel = listlen;
        list[sel].area    = area;
        list[sel].letter  = cmd->ch_target;
        list[sel].zaehler = 0;
        list[sel].maxzaehl= ilen * A[area].d_n * (int)(1.0/cmd->moment + 0.01);
        list[sel].values  = (DOUBLE *)calloc(A[area].d_n, sizeof(DOUBLE));
        list[sel].oldvals = (DOUBLE *)calloc(A[area].d_n, sizeof(DOUBLE));
        for (i = 0; i < A[area].d_n; ++i) {
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
    if  ((list[sel].zaehler % list[sel].maxzaehl < A[area].d_n) &&
         (list[sel].zaehler >= A[area].d_n)) {

        /* fprintf (stderr, "%d", area); */

        list[sel].oldvals[ct_n] = list[sel].values[ct_n]
                       / (DOUBLE)(list[sel].maxzaehl) * (DOUBLE)(A[area].d_n);

        list[sel].values[ct_n] = 0.0;
    }

    if  (list[sel].zaehler % list[sel].maxzaehl == A[area].d_n)
        list[sel].zaehler = A[area].d_n;

    list[sel].zaehler += 1;

    return (list[sel].oldvals[ct_n]);            /** =1 at first time**/
}
