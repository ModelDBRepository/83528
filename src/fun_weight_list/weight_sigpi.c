#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "../kernel/coco.h"
#include "../kernel/series.h"
#include "../kernel/utils.h"
#include "weight_sigpi.h"
#include "weight_list.h"




/**function not seen in other file, thus can be same name as in weight_list.c**/
int count_connections (CONN_SIGPI *conn) {
    int count = 0;

    if  (conn->next) {
        do  {
            conn = conn->next;
            count += 1;

        } while (conn->next != NULL);

        return count;

    } else {
        return 0;
    }
}


/**function not seen in other file, thus can be same name as in weight_list.c**/
void insert_conn (CONN_SIGPI *conn, DOUBLE newval, int a1, int b1, int n1, int a2, int b2, int n2) {

    CONN_SIGPI *newconn = (CONN_SIGPI *) malloc (sizeof (CONN_SIGPI));
    CONN_SIGPI *overnext = conn->next;
    conn->next = newconn;
    newconn->next = overnext;

    newconn->val  = newval;
    newconn->a1    = a1;
    newconn->b1    = b1;
    newconn->n1    = n1;
    newconn->a2    = a2;
    newconn->b2    = b2;
    newconn->n2    = n2;
}


/**used by: weight_sigpi_alloc_full, weight_sigpi_alloc_import.
   Does the necessary init only once automatically and points cmd->pointers[0]->data to these weights.
   Every neuron gets a NULL pointer to every inarea pair.
**/

DOUBLE weight_sigpi_alloc_init_once (PARAMS *g, AREA *A, COMMAND *cmd, int, int) {

    static int firsttime_weightsNotYetAllocated = 1; /**to allow several functions to use this even though it has to be done only once**/

    if  (firsttime_weightsNotYetAllocated) {

        /**allocate pointer space**/
        ALL_WEIGHTS_SIGPI *W = (ALL_WEIGHTS_SIGPI *) malloc (sizeof (ALL_WEIGHTS_SIGPI));

        /**architecture**/
        W->areas = g->areas;

        W->all_conn = (CONN_SIGPI *****) malloc (W->areas * sizeof (CONN_SIGPI ****));

        W->d_a = (int *) malloc (W->areas * sizeof (int));
        W->d_b = (int *) malloc (W->areas * sizeof (int));
        W->d_n = (int *) malloc (W->areas * sizeof (int));

        /**area parameters**/
        for (int ct_ar = 0; ct_ar < W->areas; ++ct_ar) {

            W->d_a[ct_ar] = A[ct_ar].d_a;
            W->d_b[ct_ar] = A[ct_ar].d_b;
            W->d_n[ct_ar] = A[ct_ar].d_n;

            W->all_conn[ct_ar] = (CONN_SIGPI ****) malloc (W->d_n[ct_ar] * sizeof (CONN_SIGPI ***));
        }

        /**neuron's pointers**/
        for (int ct_ar = 0; ct_ar < W->areas; ++ct_ar) {

            for (int ct_n = 0; ct_n < W->d_n[ct_ar]; ++ct_n) {

                W->all_conn[ct_ar][ct_n] = (CONN_SIGPI ***) malloc (W->areas * sizeof (CONN_SIGPI **));

                /**init each neuron's pointer to all (input) areas to NULL**/
                for (int ct_ar_ar1 = 0; ct_ar_ar1 < W->areas; ++ct_ar_ar1) {

                    W->all_conn[ct_ar][ct_n][ct_ar_ar1] = (CONN_SIGPI **) malloc (W->areas * sizeof (CONN_SIGPI *));

                    for (int ct_ar_ar2 = 0; ct_ar_ar2 < W->areas; ++ct_ar_ar2)

                        W->all_conn[ct_ar][ct_n][ct_ar_ar1][ct_ar_ar2] = NULL;
                }
            }
        }

        /**make the cmd pointer show to the right place (cmd->pointers[0] must NOT be lost!)**/
        cmd->pointers[0]->data = W;

        firsttime_weightsNotYetAllocated = 0;
    }

    return (DOUBLE)0;
}


/************************ weight_sigpi_alloc_full ****************************/
/* Allocates a full sigmapi connectivity between target- and input areas.    */
/* ATTENTION: DOUBLE ALLOCATION IF AREA MODULATED BY ITSELF! ADD CODE TO     */
/* PREVENT THIS! CONSIDER ALSO SEVERAL ARAES AND MODULATING AREAS!            */
/* q[0][0/1]=min/max of random initial connection values.                    */

DOUBLE weight_sigpi_alloc_full (PARAMS *g, AREA *A, COMMAND *cmd, int, int) {

    if  (cmd->anz_quant[0] != 2)
        fprintf (stderr, "\nweight_sigpi_alloc_full: insufficient parameters!");

    const DOUBLE upper = cmd->quantum[0][1];
    const DOUBLE lower = cmd->quantum[0][0];

    weight_sigpi_alloc_init_once (g, A, cmd, 0, 0);

    ALL_WEIGHTS_SIGPI *W = (ALL_WEIGHTS_SIGPI *) cmd->pointers[0]->data;

    /**all input areas in command list**/
    for (int ct_l1 = 0; ct_l1 < cmd->anz_from1; ++ct_l1)
    for (int ct_l2 = 0; ct_l2 < cmd->anz_from2; ++ct_l2) {

        const int inarea1 = cmd->n_from1[ct_l1];                       /**area number in global list**/
        const int inarea2 = cmd->n_from2[ct_l2];                       /**area number in global list**/

        for (int ct_n = 0; ct_n < W->d_n[cmd->area]; ++ct_n) {

            W->all_conn[cmd->area][ct_n][inarea1][inarea2] = (CONN_SIGPI *) malloc (sizeof (CONN_SIGPI));      /**first connection**/
            CONN_SIGPI *conn = W->all_conn[cmd->area][ct_n][inarea1][inarea2];

            /**create the admin element**/
            conn->val  = 0;
            conn->a1   = -1;
            conn->b1   = -1;
            conn->n1   = -1;
            conn->a2   = -1;
            conn->b2   = -1;
            conn->n2   = -1;
            conn->next = NULL;

            /**create the linked list**/
            for (int ct_a1 = 0; ct_a1 < W->d_a[inarea1]; ++ct_a1)
            for (int ct_b1 = 0; ct_b1 < W->d_b[inarea1]; ++ct_b1)
            for (int ct_a2 = 0; ct_a2 < W->d_a[inarea2]; ++ct_a2)
            for (int ct_b2 = 0; ct_b2 < W->d_b[inarea2]; ++ct_b2) {

                    if  (conn->next == NULL)
                        conn->next = (CONN_SIGPI *) malloc (sizeof (CONN_SIGPI));

                    conn       = conn->next;

                    conn->val  = lower + drand48() * (upper - lower);
                    conn->a1   = ct_a1;
                    conn->b1   = ct_b1;
                    conn->n1   = ct_a1 * W->d_b[inarea1] + ct_b1;
                    conn->a2   = ct_a2;
                    conn->b2   = ct_b2;
                    conn->n2   = ct_a2 * W->d_b[inarea2] + ct_b2;

                    if  ((conn->n1 == W->d_n[inarea1] - 1) && (conn->n2 == W->d_n[inarea2] - 1))
                        conn->next = NULL;
                    else
                        conn->next = (CONN_SIGPI *) malloc (sizeof (CONN_SIGPI));
            }
        }
    }

    return (DOUBLE)(0);
}



/**************************** weight_sigpi_feed ******************************/

DOUBLE weight_sigpi_feed (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    ALL_WEIGHTS_SIGPI *W = (ALL_WEIGHTS_SIGPI *)cmd->pointers[0]->data;

    DOUBLE act = (DOUBLE)(0);

    int intime = (ct_t + (int)(cmd->quantum[0][0]) < 0)
               ? 0
               : ct_t + (int)(cmd->quantum[0][0]);

    /**all input areas in command list**/
    for (int ct_l1 = 0; ct_l1 < cmd->anz_from1; ++ct_l1)
    for (int ct_l2 = 0; ct_l2 < cmd->anz_from2; ++ct_l2) {

        const int inarea1 = cmd->n_from1[ct_l1];                         /**area number in global list**/
        const int inarea2 = cmd->n_from2[ct_l2];                         /**area number in global list**/

        CONN_SIGPI *conn = W->all_conn[cmd->area][ct_n][inarea1][inarea2];     /**[area][ct_n][glob inarea][glob inarea2](first inarea1/inarea2 neuron)**/

        DOUBLE *input1 = cmd->S_from1[ct_l1][intime];                    /**[cmd inarea list][ct_t](inarea1 neurons)**/
        DOUBLE *input2 = cmd->S_from2[ct_l2][intime];                    /**[cmd inarea list][ct_t](inarea2 neurons)**/

        /**sum over weighted inputs**/
        if  (conn->next)
        do  {
            conn = conn->next;
            act += conn->val * input1[conn->n1] * input2[conn->n2];

        } while (conn->next != NULL);
    }

    return act;
}

/**************************** weight_sigpi_euclid ****************************/
/* Returns Euclidian distance of neuron ct_n's weight vector to data point.  */
/* Data and weights in from2 list; Depends on time ct_t.                     */

DOUBLE weight_sigpi_euclid (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    double dist = 0.0;

    ALL_WEIGHTS_SIGPI *W = (ALL_WEIGHTS_SIGPI *)cmd->pointers[0]->data;

    int intime = (ct_t + (int)(cmd->quantum[0][0]) < 0)
               ? 0
               : ct_t + (int)(cmd->quantum[0][0]);


    /**all input areas in command list**/
    for (int ct_l1 = 0; ct_l1 < cmd->anz_from1; ++ct_l1)
    for (int ct_l2 = 0; ct_l2 < cmd->anz_from2; ++ct_l2) {

        const int inarea1 = cmd->n_from1[ct_l1];                         /**area number in global list**/
        const int inarea2 = cmd->n_from2[ct_l2];                         /**area number in global list**/

        CONN_SIGPI *conn = W->all_conn[cmd->area][ct_n][inarea1][inarea2];     /**[area][ct_n][glob inarea][glob inarea2](first inarea1/inarea2 neuron)**/

        DOUBLE *input1 = cmd->S_from1[ct_l1][intime];                    /**[cmd inarea list][ct_t](inarea1 neurons)**/
        DOUBLE *input2 = cmd->S_from2[ct_l2][intime];                    /**[cmd inarea list][ct_t](inarea2 neurons)**/

        /**sum over weighted inputs**/
        if  (conn->next)
        do  {
            conn = conn->next;
            dist += (conn->val - input1[conn->n1] * input2[conn->n2]) * (conn->val - input1[conn->n1] * input2[conn->n2]);

        } while (conn->next != NULL);
    }

    return (sqrt (dist));
}




/**************************** weight_sigpi_kohonen ***************************/
/* Like weight_hebb, but pre[] replaced by (pre[conn->n] - conn->val).       */
/* Expects neighborhood function as post.                                    */

DOUBLE weight_sigpi_kohonen (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    ALL_WEIGHTS_SIGPI *W = (ALL_WEIGHTS_SIGPI *)cmd->pointers[0]->data;

    const DOUBLE post = cmd->S_from1[0][ct_t][ct_n];                     /**[ct_l="here"][ct_t][ct_n]**/

    if  (post != 0.0) {

        /**check only**/
        if  ((cmd->anz_from1 != 1) || (cmd->n_from1[0] != cmd->area))
            fprintf (stderr, "wrong use of weight_sigpi_kohonen");

        /**all input areas in command list**/
        for (int ct_l2 = 0; ct_l2 < cmd->anz_from2; ++ct_l2)
        for (int ct_l3 = 0; ct_l3 < cmd->anz_from3; ++ct_l3) {

            const int inarea2 = cmd->n_from2[ct_l2];                           /**area number in global list**/
            const int inarea3 = cmd->n_from3[ct_l3];                           /**area number in global list**/

            CONN_SIGPI *conn = W->all_conn[cmd->area][ct_n][inarea2][inarea3];               /**[area][ct_n][glob inarea](inarea admin neuron)**/

            DOUBLE *pre2 = cmd->S_from2[ct_l2][ct_t];                          /**[cmd inarea list][ct_t](inarea neurons)**/
            DOUBLE *pre3 = cmd->S_from3[ct_l3][ct_t];                          /**[cmd inarea list][ct_t](inarea neurons)**/

            /**input connections**/
            if  (conn->next)
            do  {
                conn = conn->next;
                conn->val += cmd->moment * cmd->quantum[0][0] * post * (pre2[conn->n1] * pre3[conn->n2] - conn->val);

            } while (conn->next != NULL);
        }
    }

    return (DOUBLE)(0);
}




/************************ weight_sigpi_histogram ******************************/
/* q[0][0] number of bins                                                     */

DOUBLE weight_sigpi_histogram (PARAMS *g, AREA *A, COMMAND *cmd, int, int) {

    fprintf (stderr, "\nweight_sigpi_histogram:");

    ALL_WEIGHTS_SIGPI *W = (ALL_WEIGHTS_SIGPI *)cmd->pointers[0]->data;

    const int area   = cmd->area;

    for (int ct_l1 = 0; ct_l1 < cmd->anz_from1; ++ct_l1)
    for (int ct_l2 = 0; ct_l2 < cmd->anz_from2; ++ct_l2) {

        DOUBLE max = -99;    /**for min and max weight**/
        DOUBLE min = 99;

        int count_all = 0;
        int number_without = 0;

        const int inarea1 = cmd->n_from1[ct_l1];                         /**area number in global list**/
        const int inarea2 = cmd->n_from2[ct_l2];                         /**area number in global list**/

        for (int ct_n = 0; ct_n < A[area].d_n; ++ct_n) {

            CONN_SIGPI *conn = W->all_conn[area][ct_n][inarea1][inarea2];
            int count_here = 0;

            if  (conn->next)
            do  {
                conn = conn->next;
                max = conn->val > max ? conn->val : max;
                min = conn->val < min ? conn->val : min;

                count_all += 1;
                count_here += 1;
 
            } while (conn->next != NULL);

            if  (count_here == 0)
                number_without += 1;
        }


        int bins = (int)cmd->quantum[0][0];
        int bincount[bins];
        DOUBLE interval = (max - min) / bins;

        /**init bins zero**/
        for (int i = 0; i < bins; ++i)
            bincount[i] = 0;

        /**for each neuron**/
        for (int ct_n = 0; ct_n < A[area].d_n; ++ct_n) {

            /**for all its existing connections**/
            CONN_SIGPI *conn = W->all_conn[area][ct_n][inarea1][inarea2];
            if  (conn->next)
            do  {
                conn = conn->next;

                DOUBLE lower = min;
                DOUBLE upper = min;

                /**assign the connection into its bin**/
                for (int i = 0; i < bins; ++i) {
                    upper += interval;
                    if  ((conn->val >= lower) && (conn->val < upper))
                        bincount[i] += 1;
                    lower += interval;
                }

            } while (conn->next != NULL);
        }

        /**print out**/
        fprintf (stderr, "\nW %d<-%d*%d strengths:", area, inarea1, inarea2);
        DOUBLE lower = min;
        DOUBLE upper = min;
        for (int i = 0; i < bins; ++i) {
            upper += interval;
            fprintf (stderr, "\n[% f - % f]: %d ", lower, upper, bincount[i]);
            lower += interval;
        }

        /**test weights for consistent coverage and a,b,n consistent?**/
        int ****connected = i_4tensor(W->d_a[inarea1], W->d_b[inarea1], W->d_a[inarea2], W->d_b[inarea2]);

        for (int ct_n = 0; ct_n < A[area].d_n; ++ct_n) {

            for (int i = 0; i < W->d_a[inarea1]; ++i)
            for (int j = 0; j < W->d_b[inarea1]; ++j)
            for (int k = 0; k < W->d_a[inarea2]; ++k)
            for (int h = 0; h < W->d_b[inarea2]; ++h)
                connected[i][j][k][h] = 0;

            CONN_SIGPI *conn = W->all_conn[area][ct_n][inarea1][inarea2];

            if  (conn->next)
            do  {
                conn = conn->next;

                connected[conn->a1][conn->b1][conn->a2][conn->b2] += 1;

                if  (conn->a1 * A[inarea1].d_b + conn->b1 != conn->n1)
                    fprintf (stderr, "\n\n\nInconsistent d_a1, d_b1 with d_n1!\n\n\n");
                if  (conn->a2 * A[inarea2].d_b + conn->b2 != conn->n2)
                    fprintf (stderr, "\n\n\nInconsistent d_a2, d_b2 with d_n2!\n\n\n");

            } while (conn->next != NULL);

            for (int i = 0; i < W->d_a[inarea1]; ++i)
            for (int j = 0; j < W->d_b[inarea1]; ++j)
            for (int k = 0; k < W->d_a[inarea2]; ++k)
            for (int h = 0; h < W->d_b[inarea2]; ++h)
                if  (connected[i][j][k][h] != 0)
                    if  (connected[i][j][k][h] != 1) {
                        fprintf (stderr, "\n\nInconsistent connected matrix!\n");
                        fprintf (stderr, "inarea1=%d inarea2=%d n=%d  connected[%d][%d][%d][%d]=%d\n\n", inarea1, inarea2, ct_n, i, j, k, h, connected[i][j][k][h]);
                    }
        }

        free_i_4tensor (connected, W->d_a[inarea1], W->d_b[inarea1], W->d_a[inarea2]);

        fprintf (stderr, "\n%s_%d_%d_%d has", cmd->ch_pointers[0], area, inarea1, inarea2);
        fprintf (stderr, " in average %.1f connections of %d possible ", (double)count_all/(double)A[area].d_n, A[inarea1].d_n * A[inarea2].d_n);
        fprintf (stderr, " and %d units have no connection  ", number_without);
    }


    return (DOUBLE)(0);
}




/**************************** weight_sigpi_hebb ******************************/
/* Works on every neuron but no return val. Updates only if act[ct_n] != 0.0.*/
/* cmd->pointers[0] is a ALL_WEIGHTS_SIGPI *. Arguments *g and *A are not used.    */

DOUBLE weight_sigpi_hebb (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    ALL_WEIGHTS_SIGPI *W = (ALL_WEIGHTS_SIGPI *)cmd->pointers[0]->data;

    const DOUBLE post = cmd->S_from1[0][ct_t][ct_n];                  /**[ct_l="here"][ct_t][ct_n]**/

    if  (post != 0.0) {

        /**check only**/
        if  ((cmd->anz_from1 != 1) || (cmd->n_from1[0] != cmd->area))
            fprintf (stderr, "wrong use of weight_sigpi_hebb");

        /**all input areas in command list**/
        for (int ct_l2 = 0; ct_l2 < cmd->anz_from2; ++ct_l2)
        for (int ct_l3 = 0; ct_l3 < cmd->anz_from3; ++ct_l3) {

            const int inarea2 = cmd->n_from2[ct_l2];                           /**area number in global list**/
            const int inarea3 = cmd->n_from3[ct_l3];                           /**area number in global list**/

            CONN_SIGPI *conn = W->all_conn[cmd->area][ct_n][inarea2][inarea3];               /**[area][ct_n][glob inarea](inarea admin neuron)**/

            DOUBLE *pre2 = cmd->S_from2[ct_l2][ct_t];                          /**[cmd inarea list][ct_t](inarea neurons)**/
            DOUBLE *pre3 = cmd->S_from3[ct_l3][ct_t];                          /**[cmd inarea list][ct_t](inarea neurons)**/

            /**input connections**/
            if  (conn->next)
            do  {
                conn = conn->next;
                conn->val += cmd->moment * cmd->quantum[0][0] * post * pre2[conn->n1] * pre3[conn->n2];

            } while (conn->next != NULL);

        }
    }

    return (DOUBLE)(0);
}





/**commenting out the following code**/
#if 0





/**************************** weight_sigpi_decay *****************************/
/* Decays weights by -m * conn->val * q[0][0]. Returns weight vector length. */

DOUBLE weight_sigpi_decay (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    ALL_WEIGHTS_SIGPI *W = (ALL_WEIGHTS_SIGPI *)cmd->pointers[0]->data;
    DOUBLE quad_length = (DOUBLE)0;

    /**all input areas in command list**/
    for (int ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        const int inarea = cmd->n_from1[ct_l];                           /**area number in global list**/

        CONN_SIGPI *conn = W->all_conn[cmd->area][ct_n][inarea];               /**[area][ct_n][glob inarea](inarea admin neuron)**/

        /**sum over weighted inputs**/
        if  (conn->next)
        do  {
            conn = conn->next;
            conn->val -= cmd->moment * conn->val * cmd->quantum[0][0];
            quad_length += conn->val * conn->val;

        } while (conn->next != NULL);
    }

    return quad_length;
}




/**************************** weight_sigpi_rectify ***************************/

DOUBLE weight_sigpi_rectify (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    ALL_WEIGHTS_SIGPI *W = (ALL_WEIGHTS_SIGPI *)cmd->pointers[0]->data;

    /**all input areas in command list**/
    for (int ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        const int inarea = cmd->n_from1[ct_l];                           /**area number in global list**/

        CONN_SIGPI *conn = W->all_conn[cmd->area][ct_n][inarea];               /**[area][ct_n][glob inarea](first inarea neuron)**/

        /**sum over weighted inputs**/
        if  (conn->next)
        do  {
            conn = conn->next;

            if  (cmd->quantum[1][0] == 1.0)
                conn->val = conn->val < cmd->quantum[0][0] ? cmd->quantum[0][0] : conn->val;
            else
                conn->val = conn->val > cmd->quantum[0][0] ? cmd->quantum[0][0] : conn->val;

        } while (conn->next != NULL);
    }

    return (DOUBLE)0;
}



/**************************** weight_sigpi_cutself ***************************/

DOUBLE weight_sigpi_cutself (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int num_cut = 0;

    /**test whether inner area connections exist at all**/
    int OK = 0;
    for (int ct_l = 0; ct_l < cmd->anz_from1; ++ct_l)
        if  (cmd->area == cmd->n_from1[ct_l])
            OK = 1;
    if  (!OK)
        fprintf (stderr, "\n\nweight_sigpi_noself applied to non-self-connections!\n");


    ALL_WEIGHTS_SIGPI *W = (ALL_WEIGHTS_SIGPI *)cmd->pointers[0]->data;

    CONN_SIGPI *conn = W->all_conn[cmd->area][ct_n][cmd->area];
    int more = 1;
    while (more) {

        CONN_SIGPI *prev = conn;
        conn = conn->next;
        CONN_SIGPI *tofree = NULL;

        if  (tofree) {
            free (tofree);
            tofree = NULL;
        }

        more = (conn->next == NULL) ? 0 : 1;    /**necessary, because conn can be removed**/

        if  (conn->n == ct_n) {

            if  (conn->next == NULL)    /**last element**/
                prev->next = NULL;
            else
                prev->next = conn->next;

            tofree = conn;

            num_cut += 1;
        }
    }

    if  (num_cut != 1)
        fprintf (stderr, "\n\nweight_sigpi_cutself has cut %d connection(s) at neuron %d instead of one!\n\n", num_cut, ct_n);

    return (DOUBLE)(0);
}



#endif



/**************************** weight_sigpi_cutsmall **************************/
/* q[0][0] = 1: regard cutting threshold as absolute; = 0 in % of min/max.   */
/* q[1][0] = negative cutting threshold in % of minimum value                */
/* q[1][1] = positive cutting threshold in % of maximum value                */
/* "single" function!                                                        */

DOUBLE weight_sigpi_cutsmall (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int num_cut = 0;
    DOUBLE min = 99;
    DOUBLE max = -99;

    if  ((cmd->anz_quantums != 2) || (cmd->anz_quant[1] != 2)) {
        fprintf (stderr, "\ninsufficient arguments in weight_sigpi_cutsmall\n");
        exit (1);
    }

    ALL_WEIGHTS_SIGPI *W = (ALL_WEIGHTS_SIGPI *)cmd->pointers[0]->data;

    /**all input areas in command list**/
    for (int ct_l1 = 0; ct_l1 < cmd->anz_from1; ++ct_l1)
    for (int ct_l2 = 0; ct_l2 < cmd->anz_from2; ++ct_l2) {

        const int inarea1 = cmd->n_from1[ct_l1];                         /**area number in global list**/
        const int inarea2 = cmd->n_from2[ct_l2];                         /**area number in global list**/

        CONN_SIGPI *conn = W->all_conn[cmd->area][ct_n][inarea1][inarea2];     /**[area][ct_n][glob inarea][glob inarea2](first inarea1/inarea2 neuron)**/

        static int average_cut  = 0;                    /**to print out for info**/
        static int total_before = 0;                    /**to print out for info**/
        static int total_after  = 0;
        if  (ct_n == 0) {
            fprintf (stderr, "\nweight_sigpi_cutsmall ");
            average_cut  = 0;
            total_before = 0;
            total_after  = 0;
        }

        total_before += count_connections (W->all_conn[cmd->area][ct_n][inarea1][inarea2]);

        if  (conn->next)
        do  {
            conn = conn->next;
            min = conn->val < min ? conn->val : min;
            max = conn->val > max ? conn->val : max;
        } while (conn->next != NULL);

        /**absolute cut mode**/
        if  (cmd->quantum[0][0] == 1) {
            max = 1;
            min = -1;
        }

        DOUBLE thres_neg = cmd->quantum[1][0] * min;
        DOUBLE thres_pos = cmd->quantum[1][1] * max;

        conn = W->all_conn[cmd->area][ct_n][inarea1][inarea2];

        int last;
        if  (conn->next)
        do  {
            CONN_SIGPI *prev = conn;
            conn = conn->next;

            last = (conn->next == NULL) ? 1 : 0;

            if  ((conn->val >= thres_neg) && (conn->val <= thres_pos)) {

                 if  (last)
                     prev->next = NULL;
                 else
                     prev->next = conn->next;  

                 num_cut += 1;
                 free (conn);
                 conn = prev;
            }
        } while (! last);

        average_cut += num_cut;
        total_after += count_connections (W->all_conn[cmd->area][ct_n][inarea1][inarea2]);

        if  (ct_n == A[cmd->area].d_n - 1) {
            fprintf (stderr, " cut %.1f connections in average, or %d total at obs_%s_%d_%d_%d  ",
                              (double)average_cut/(double)A[cmd->area].d_n, average_cut, cmd->ch_pointers[0], cmd->area, inarea1, inarea2);
            fprintf (stderr, "\ntotal connections before: %d;  cut: %d;  remain: %d  ", total_before, total_before - total_after, total_after);
        }
    }

    return (DOUBLE)(0);
}





/**commenting out the following code**/
#if 0



/**************************** weight_sigpi_sprout ****************************/
/* q[0][0] = connection value in % of minimum value around which to sprout   */
/* q[0][1] = connection value in % of maximum value around which to sprout   */
/* q[1][0] = value of new connections in % of minimum value                  */
/* q[1][1] = value of new connections in % of maximum value                  */
/* q[2][0] = 4/8 for sprouting at N/S/W/E or also at diagonal directions     */
/* q[3][0] = 3 for RGB-inarea, so sprouting does not cross the boundaries    */
/* q[3][1] = 2 for separate ON/OFF-inarea,    "                              */

DOUBLE weight_sigpi_sprout (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    if  (ct_n == 0)
        fprintf (stderr, "\nweight_sigpi_sprout ");

    DOUBLE max = -99;
    DOUBLE min = 99;

    ALL_WEIGHTS_SIGPI *W = (ALL_WEIGHTS_SIGPI *)cmd->pointers[0]->data;

    const int inarea = cmd->n_from1[0];       /**use only 1st inarea in list**/

    static int average = 0;                         /**to print out for info**/
    static int total_before = 0;                    /**to print out for info**/
    static int total_after  = 0;
    if  (ct_n == 0) {
        average = 0;
        total_before    = 0;
        total_after     = 0;
    }

    int **connected = i_matrix(W->d_a[inarea], W->d_b[inarea]);

    total_before += count_connections (W->all_conn[cmd->area][ct_n][inarea]);

    for (int i = 0; i < W->d_a[inarea]; ++i)
    for (int j = 0; j < W->d_b[inarea]; ++j)
        connected[i][j] = 0;

    CONN_SIGPI *conn = W->all_conn[cmd->area][ct_n][inarea];

    if  (conn->next)
    do  {
        conn = conn->next;

        connected[conn->a][conn->b] = 1;

        max = conn->val > max ? conn->val : max;
        min = conn->val < min ? conn->val : min;

    } while (conn->next != NULL);

    DOUBLE thres_neg = cmd->quantum[0][0] * min;
    DOUBLE thres_pos = cmd->quantum[0][1] * max;

    int count = 0;

    conn = W->all_conn[cmd->area][ct_n][inarea];               /**[area][ct_n][glob inarea](inarea admin neuron)**/

    /**loop over inputs**/
    if  (conn->next)
    do  {
        conn = conn->next;

        /**around large connections ...**/
        if  ((conn->val > thres_pos) || (conn->val < thres_neg)) {

            /**... put 4 connections with new value**/
            DOUBLE newval = (conn->val < thres_neg) ? min * cmd->quantum[1][0] : max * cmd->quantum[1][1];

            CONN_SIGPI *here = conn;
            int num_new = 0;

            int over_RGBborder    = 0;   /**test here->a+1     to prevent sprouting crossing the borders between the three RGB sub-areas**/
            int under_RGBborder   = 0;   /**test here->a-1**/
            int over_ONOFFborder  = 0;   /**test here->b+1     to prevent sprouting crossing the borders between the two ON/OFF sub-areas**/
            int under_ONOFFborder = 0;   /**test here->b-1**/

            if  (cmd->quantum[3][0] == 3)
                if  ((here->a == W->d_a[inarea] / 3) || (here->a == W->d_a[inarea] / 3 * 2))
                    under_RGBborder = 1;

            if  (cmd->quantum[3][0] == 3)
                if  ((here->a == W->d_a[inarea] / 3 - 1) || (here->a == W->d_a[inarea] / 3 * 2 - 1))
                    over_RGBborder = 1;

            if  (cmd->quantum[3][1] == 2)
                if  (here->b == W->d_b[inarea] / 2)
                    under_ONOFFborder = 1;

            if  (cmd->quantum[3][1] == 2)
                if  (here->b == W->d_b[inarea] / 2 - 1)
                    over_ONOFFborder = 1;


            /**below a*/
            if  (! under_RGBborder)
                if  (here->a - 1 >= 0)
                    if  (connected[here->a - 1][here->b] == 0) {
                        insert_conn (here, newval, here->a - 1, here->b, (here->a - 1) * W->d_b[inarea] + here->b);
                        connected[here->a - 1][here->b] = 1;
                        num_new += 1;
                    }

            /**below a and below b**/
            if  ((! under_RGBborder) && (! under_ONOFFborder))
                if  ((here->a - 1 >= 0) && (here->b - 1 >= 0) && (cmd->quantum[2][0] == 8))
                    if  (connected[here->a - 1][here->b - 1] == 0) {
                        insert_conn (here, newval, here->a - 1, here->b - 1, (here->a - 1) * W->d_b[inarea] + here->b - 1);
                        connected[here->a - 1][here->b - 1] = 1;
                        num_new += 1;
                    }

            /**below a and over b**/
            if  ((! under_RGBborder) && (! over_ONOFFborder))
                if  ((here->a - 1 >= 0) && (here->b + 1 < W->d_b[inarea]) && (cmd->quantum[2][0] == 8))
                    if  (connected[here->a - 1][here->b + 1] == 0) {
                        insert_conn (here, newval, here->a - 1, here->b + 1, (here->a - 1) * W->d_b[inarea] + here->b + 1);
                        connected[here->a - 1][here->b + 1] = 1;
                        num_new += 1;
                    }

            /**over a**/
            if  (! over_RGBborder)
                if  (here->a + 1 < W->d_a[inarea])
                    if  (connected[here->a + 1][here->b] == 0) {
                        insert_conn (here, newval, here->a + 1, here->b, (here->a + 1) * W->d_b[inarea] + here->b);
                        connected[here->a + 1][here->b] = 1;
                        num_new += 1;
                    }

            /**over a and below b**/
            if  ((! over_RGBborder) && (! under_ONOFFborder))
                if  ((here->a + 1 < W->d_a[inarea]) && (here->b - 1 >= 0) && (cmd->quantum[2][0] == 8))
                    if  (connected[here->a + 1][here->b - 1] == 0) {
                        insert_conn (here, newval, here->a + 1, here->b - 1, (here->a + 1) * W->d_b[inarea] + here->b - 1);
                        connected[here->a + 1][here->b - 1] = 1;
                        num_new += 1;
                    }

            /**over a and over b**/
            if  ((! over_RGBborder) && (! over_ONOFFborder))
                if  ((here->a + 1 < W->d_a[inarea]) && (here->b + 1 < W->d_b[inarea]) && (cmd->quantum[2][0] == 8))
                    if  (connected[here->a + 1][here->b + 1] == 0) {
                        insert_conn (here, newval, here->a + 1, here->b + 1, (here->a + 1) * W->d_b[inarea] + here->b + 1);
                        connected[here->a + 1][here->b + 1] = 1;
                        num_new += 1;
                    }

            /**below b**/
            if  (! under_ONOFFborder)
                if  (here->b - 1 >= 0)
                    if  (connected[here->a][here->b - 1] == 0) {
                        insert_conn (here, newval, here->a, here->b - 1, here->a * W->d_b[inarea] + here->b - 1);
                        connected[here->a][here->b - 1] = 1;
                        num_new += 1;
                    }

            /**over b**/
            if  (! over_ONOFFborder)
                if  (here->b + 1 < W->d_b[inarea])
                    if  (connected[here->a][here->b + 1] == 0) {
                        insert_conn (here, newval, here->a, here->b + 1, here->a * W->d_b[inarea] + here->b + 1);
                        connected[here->a][here->b + 1] = 1;
                        num_new += 1;
                    }


            for (int i = 0; i < num_new; ++i)
                conn = conn->next;

            count += num_new;
        }

    } while (conn->next != NULL);

    free_i_matrix (connected, W->d_a[inarea]);

    average += count;
    total_after += count_connections (W->all_conn[cmd->area][ct_n][inarea]);

    if  (ct_n == A[cmd->area].d_n - 1) {
        fprintf (stderr, " added %.1f connections in average or total %d at obs_%s_%d_%d ", (double)average/(double)(A[cmd->area].d_n), average, cmd->ch_pointers[0], cmd->area, inarea);
        fprintf (stderr, "\ntotal connections after - before sprout: %d - %d = %d  ", total_after, total_before, total_after - total_before);
    }

    return (DOUBLE)(0);
}

#endif



void weightssigpi2rgbmatrix (ALL_WEIGHTS_SIGPI *W, int area, int inarea1, int inarea2, DOUBLE ***rgbmatrix, int d_a, int d_b, int d_in1, int d_in2, DOUBLE *min, DOUBLE *max, int format) {

    if  ((format != 3) && (format != 9))
        format = 6;

    DOUBLE maxcol;

    if  (format == 3)
        maxcol = MAXCOL_INT;
    else /**format == 6**/
        maxcol = MAXCOL_CHAR;

    /**background init with a light green tone**/
    for (int i = 0; i < d_a * d_b; ++i)
        for (int j = 0; j < d_in1 * d_in2; ++j) {
            rgbmatrix[0][i][j] = 0.8 * maxcol;    /**red**/
            rgbmatrix[1][i][j] = maxcol;          /**green**/
            rgbmatrix[2][i][j] = 0.8 * maxcol;    /**blue**/

            if  (format == 9)
                rgbmatrix[0][i][j] = 0.0;         /**necessary for non-existing connections**/
        }

 
    /**find maximum and minimum weight strength**/
    *max = -99;
    *min = 99;
    for (int ct_n = 0; ct_n < d_a * d_b; ++ct_n) {

        CONN_SIGPI *conn = W->all_conn[area][ct_n][inarea1][inarea2];

        if  (conn->next)
        do  {
            conn = conn->next;
            *max = conn->val > *max ? conn->val : *max;
            *min = conn->val < *min ? conn->val : *min;

        } while (conn->next != NULL);
    }


    /**set rgbmatrix values where weights exist; positive&negative are normalised together w.r.t. to largest of max/min**/
    DOUBLE largest = *max > -*min ? *max : -*min;

    for (int ct_n = 0; ct_n < d_a * d_b; ++ct_n) {

        CONN_SIGPI *conn = W->all_conn[area][ct_n][inarea1][inarea2];

        if  (conn->next)
        do  {
            conn = conn->next;

            /**positive (blue remains bright)**/
            if  (conn->val >= 0) {
                rgbmatrix[0][ct_n][conn->n1 * d_in2 + conn->n2] = maxcol * (1.0 - conn->val / largest);
                rgbmatrix[1][ct_n][conn->n1 * d_in2 + conn->n2] = maxcol * (1.0 - conn->val / largest);
                rgbmatrix[2][ct_n][conn->n1 * d_in2 + conn->n2] = maxcol;
            }

            /**negative (red remains bright)**/
            if  (conn->val < 0) {
                rgbmatrix[0][ct_n][conn->n1 * d_in2 + conn->n2] = maxcol;
                rgbmatrix[1][ct_n][conn->n1 * d_in2 + conn->n2] = maxcol * (1.0 + conn->val / largest);
                rgbmatrix[2][ct_n][conn->n1 * d_in2 + conn->n2] = maxcol * (1.0 + conn->val / largest);
            }

            /**in format 9 only write val into first component only**/
            if  (format == 9)
                rgbmatrix[0][ct_n][conn->n1 * d_in2 + conn->n2] = conn->val;

        } while (conn->next != NULL);
    }
}






/************************ weight_sigpi_export ********************************/
/* Export weights to target area from all inarea pairs. File looks like:     */
/* sigmapi                                                                   */
/* format = 1                                                                */
/* area = %d                                                                 */
/* d_a = %d, d_b = %d                                                        */
/* inarea1[0] = %d, inarea1[1] = %d ...                                      */
/* d_a = %d, d_b = %d, d_a = %d, d_b = %d ...                                */
/* inarea2[0] = %d, inarea2[1] = %d ...                                      */
/* d_a = %d, d_b = %d, d_a = %d, d_b = %d ...                                */
/* DATA                                                                      */
/* inarea1 = %d, inarea2 = %d                                                */
/* n = %d                                                                    */
/* %f %d %d (val n1 n2)                                                      */
/* %f %d %d                                                                  */
/* %f %d %d                                                                  */
/* ...                                                                       */
/* NULL                                                                      */
/* ... for all neurons at this area until the last ...                       */
/* n = %d                                                                    */
/* ...                                                                       */
/* NULL                                                                      */
/* LAST ("LAST" is better, if not every unit has weights and would appear)   */
/* inarea1 = %d, inarea2 = %d                                                */
/* ... for all inarea1/2 combinations until last ...                         */
/* NULL                                                                      */
/* LAST                                                                      */
/* END                                                                       */

DOUBLE weight_sigpi_export (PARAMS *g, AREA *A, COMMAND *cmd, int, int) {

    fprintf (stderr, "  weight_sigpi_export  ");

    DOUBLE min, max;
    int format = (int)(cmd->quantum[0][0]);
    int export_pnm = 0;
    if  (cmd->anz_quantums == 2)
        export_pnm = 1;

    ALL_WEIGHTS_SIGPI *W = (ALL_WEIGHTS_SIGPI *)cmd->pointers[0]->data;

    const int area   = cmd->area;
    const int inarea1 = cmd->n_from1[0];       /**use only 1st inarea in list**/
    const int inarea2 = cmd->n_from2[0];       /**use only 1st inarea in list**/

    int d_a    = A[area].d_a;
    int d_b    = A[area].d_b;

    /**export pnm file if any of the input areas' dimensions is 1**/
    if  ((A[inarea1].d_a == 1) || (A[inarea1].d_b == 1) || (A[inarea2].d_a == 1) || (A[inarea2].d_b == 1) || (export_pnm == 1)) {

        int d_in1  = A[inarea1].d_a * A[inarea1].d_b;
        int d_in2  = A[inarea2].d_a * A[inarea2].d_b;

        DOUBLE ***rgbmatrix = d_tensor (3, d_a * d_b, d_in1 * d_in2);

        weightssigpi2rgbmatrix (W, area, inarea1, inarea2, rgbmatrix, d_a, d_b, d_in1, d_in2, &min, &max, format);

        rgbmatrix2file (cmd, rgbmatrix, d_a, d_b, d_in1, d_in2, min, max, format, 2);
                                                                                /*weights*/
        free_d_tensor (rgbmatrix, 3, d_a * d_b);
    }

    /**export sigma pi weights**/

    if  (cmd->anz_pointers != 2)
        fprintf (stderr, "\n\n\nweight_sigpi_export needs two pointers: weights and directory!\n\n\n");

    char fullname[512];
    sprintf (fullname, "%s/obs_%s_%d_%d_%d.spi", cmd->pointers[1]->words[0], cmd->ch_pointers[0],
                                  cmd->area, inarea1, inarea2);    /**pointers[1] gives the directory!**/

    FILE *fp = fopen (fullname, "w");

    /**write header**/
    fprintf (fp, "sigpi\n");
    fprintf (fp, "format = 1\n");
    fprintf (fp, "area = %d\n", area);
    fprintf (fp, "d_a = %d, d_b = %d\n", d_a, d_b);
    for (int i = 0; i < cmd->anz_from1; ++i) {
        fprintf (fp, "inarea1[%d] = %d", i, cmd->n_from1[i]);
        if  (i < cmd->anz_from1 - 1)
            fprintf (fp, ", ");
    }
    fprintf (fp, "\n");
    for (int i = 0; i < cmd->anz_from1; ++i) {
        fprintf (fp, "d_a = %d, d_b = %d", A[cmd->n_from1[i]].d_a, A[cmd->n_from1[i]].d_b);
        if  (i < cmd->anz_from1 - 1)
            fprintf (fp, ", ");
    }
    fprintf (fp, "\n");
    for (int i = 0; i < cmd->anz_from2; ++i) {
        fprintf (fp, "inarea2[%d] = %d", i, cmd->n_from2[i]);
        if  (i < cmd->anz_from2 - 1)
            fprintf (fp, ", ");
    }
    fprintf (fp, "\n");
    for (int i = 0; i < cmd->anz_from2; ++i) {
        fprintf (fp, "d_a = %d, d_b = %d", A[cmd->n_from2[i]].d_a, A[cmd->n_from2[i]].d_b);
        if  (i < cmd->anz_from2 - 1)
            fprintf (fp, ", ");
    }
    fprintf (fp, "\n");
    fprintf (fp, "DATA\n");

    /**write data**/
    for (int i = 0; i < cmd->anz_from1; ++i) {
        for (int j = 0; j < cmd->anz_from2; ++j) {
            fprintf (fp, "inarea1 = %d, inarea2 = %d\n", cmd->n_from1[i], cmd->n_from2[j]);
            for (int ct_n = 0; ct_n < A[area].d_n; ++ct_n) {
                fprintf (fp, "ct_n = %d\n", ct_n);

                CONN_SIGPI *conn = W->all_conn[area][ct_n][inarea1][inarea2];

                if  (conn->next)
                do  {
                    conn = conn->next;
                    fprintf (fp, "%f %d %d\n", conn->val, conn->n1, conn->n2);
                } while (conn->next != NULL);
                fprintf (fp, "NULL\n");
            }
            fprintf (fp, "LAST\n");
        }
    }
    fprintf (fp, "END\n");

    fclose (fp);

    return (DOUBLE)0;
}


/************************ weight_sigpi_alloc_import **************************/

DOUBLE weight_sigpi_alloc_import (PARAMS *g, AREA *A, COMMAND *cmd, int, int) {

    fprintf (stderr, "\nweight_sigpi_import ");

    int format = 1;
    int area;
    int d_a;
    int d_b;

    /**import sigma pi weights**/

    if  (cmd->anz_pointers != 2)
        fprintf (stderr, "\n\n\nweight_sigpi_import needs two pointers: weights and directory!\n\n\n");

    char fullname[512];
    sprintf (fullname, "%s/obs_%s_%d_%d_%d.spi", cmd->pointers[1]->words[0], cmd->ch_pointers[0],
                                  cmd->area, cmd->n_from1[0], cmd->n_from2[0]);    /**pointers[1] gives the directory!**/

    fprintf (stderr, " reading file %s\n", fullname);

    FILE *fp = fopen (fullname, "r");

    /**read header**/
    fscanf (fp, "sigpi\n");
    fscanf (fp, "format = %d\n", &format);
    if  (format != 1)
        fprintf (stderr, "\n\nweight_sigpi_import currently takes only format=1!!\n\n\n");
    fscanf (fp, "area = %d\n", &area);
    if  (area != cmd->area)
        fprintf (stderr, "\n\nweight_sigpi_import: area conflict!!\n\n\n");
    fscanf (fp, "d_a = %d, d_b = %d\n", &d_a, &d_b);
    if  ((d_a != A[area].d_a) || (d_b != A[area].d_b))
        fprintf (stderr, "\n\nweight_sigpi_import: requested d_a=%d d_b=%d but in file d_a=%d d_b=%d!!\n\n\n", A[area].d_a, A[area].d_b, d_a, d_b);

    /**require here that inarea order is same in command as well as file**/
    for (int i = 0; i < cmd->anz_from1; ++i) {
        int ii, curr_inarea;
        fscanf (fp, "inarea1[%d] = %d", &ii, &curr_inarea);
        if  ((ii != i) || (curr_inarea != cmd->n_from1[i]))
            fprintf (stderr, "\n\nweight_sigpi_import: cmd inarea list doesn't match file inarea list!!\n\n\n");
        if  (i < cmd->anz_from1 - 1)
            fscanf (fp, ", ");
    }
    fscanf (fp, "\n");
    for (int i = 0; i < cmd->anz_from1; ++i) {
        int curr_d_a, curr_d_b;
        fscanf (fp, "d_a = %d, d_b = %d", &curr_d_a, &curr_d_b);
        if  ((curr_d_a != A[cmd->n_from1[i]].d_a) || (curr_d_b != A[cmd->n_from1[i]].d_b))
            fprintf (stderr, "\n\nweight_sigpi_import: cmd inarea resolutions don't match file inarea resolutions!!\n\n\n");
        if  (i < cmd->anz_from1 - 1)
            fscanf (fp, ", ");
    }
    fscanf (fp, "\n");

    for (int i = 0; i < cmd->anz_from2; ++i) {
        int ii, curr_inarea;
        fscanf (fp, "inarea2[%d] = %d", &ii, &curr_inarea);
        if  ((ii != i) || (curr_inarea != cmd->n_from2[i]))
            fprintf (stderr, "\n\nweight_sigpi_import: cmd inarea list doesn't match file inarea list!!\n\n\n");
        if  (i < cmd->anz_from2 - 1)
            fscanf (fp, ", ");
    }
    fscanf (fp, "\n");
    for (int i = 0; i < cmd->anz_from2; ++i) {
        int curr_d_a, curr_d_b;
        fscanf (fp, "d_a = %d, d_b = %d", &curr_d_a, &curr_d_b);
        if  ((curr_d_a != A[cmd->n_from2[i]].d_a) || (curr_d_b != A[cmd->n_from2[i]].d_b))
            fprintf (stderr, "\n\nweight_sigpi_import: cmd inarea resolutions don't match file inarea resolutions!!\n\n\n");
        if  (i < cmd->anz_from2 - 1)
            fscanf (fp, ", ");
    }
    fscanf (fp, "\n");
    fscanf (fp, "DATA\n");

    /**allocate NULL pointer for every neuron**/
    weight_sigpi_alloc_init_once (g, A, cmd, 0, 0);

    ALL_WEIGHTS_SIGPI *W = (ALL_WEIGHTS_SIGPI *) cmd->pointers[0]->data;

    /**read data**/
    for (int i = 0; i < cmd->anz_from1; ++i) {
        for (int j = 0; j < cmd->anz_from2; ++j) {
            int curr_inarea1, curr_inarea2;
            fscanf (fp, "inarea1 = %d, inarea2 = %d\n", &curr_inarea1, &curr_inarea2);
            if  ((curr_inarea1 != cmd->n_from1[i]) || (curr_inarea2 != cmd->n_from2[j]))
                fprintf (stderr, "\n\nweight_sigpi_import: file DATA inareas don't match cmd inareas!!\n\n\n");

            /**require that all neurons appear in order -- drop this requirement if not all neurons get weights exported -- would be easy to drop requirement here**/
            for (int ct_n = 0; ct_n < A[area].d_n; ++ct_n) {
                int curr_n;
                fscanf (fp, "ct_n = %d\n", &curr_n);

                if  (curr_n != ct_n)
                    fprintf (stderr, "\n\nweight_sigpi_import: neurons scrambled -- file's curr_n=%d, index' ct_n=%d!!\n\n", curr_n, ct_n);

                W->all_conn[area][curr_n][curr_inarea1][curr_inarea2] = (CONN_SIGPI *) malloc (sizeof (CONN_SIGPI));      /**first connection**/
                CONN_SIGPI *conn = W->all_conn[area][curr_n][curr_inarea1][curr_inarea2];

                /**create the admin element**/
                conn->val  = 0;
                conn->a1   = -1;
                conn->b1   = -1;
                conn->n1   = -1;
                conn->a2   = -1;
                conn->b2   = -1;
                conn->n2   = -1;
                conn->next = NULL;

                char kommentarzeile[512];
                fgets (kommentarzeile, 512, fp);

                while (strncmp (kommentarzeile, "NULL", 4)) {

                    if  (conn->next == NULL)
                        conn->next = (CONN_SIGPI *) malloc (sizeof (CONN_SIGPI));

                    conn       = conn->next;

                    sscanf (kommentarzeile, "%f %d %d\n", &(conn->val), &(conn->n1), &(conn->n2));

                    conn->a1   = conn->n1 / A[curr_inarea1].d_b;
                    conn->b1   = conn->n1 % A[curr_inarea1].d_b;
                    conn->a2   = conn->n2 / A[curr_inarea2].d_b;
                    conn->b2   = conn->n2 % A[curr_inarea2].d_b;

                    conn->next = (CONN_SIGPI *) malloc (sizeof (CONN_SIGPI));
    
                    fgets (kommentarzeile, 512, fp);
                }
                conn->next = NULL;

            }
            fscanf (fp, "LAST\n");
        }
    }
    fscanf (fp, "END\n");

    fclose (fp);

    return (DOUBLE)0;
}



/**************************** weight_sigpi_gnuplot_cm ************************/
/* Plots a grid into a 2-D data space as normally done for Kohonen weights.  */
/* However, considers CM of RF instead direct weights, for population codes. */
/* q[0][0] = 0/1: keep j/k of w_ijk constant.                                */
/* q[1][0]/[1] = j/k (the values of the selected index along hor and vert).  */
/* q[2][0] threshold value; only larger weights for CM; not implemented.     */

DOUBLE weight_sigpi_gnuplot_cm (PARAMS *g, AREA *A, COMMAND *cmd, int, int) {

    int index_a = (int)(cmd->quantum[1][0]);
    int index_b = (int)(cmd->quantum[1][1]);

    const int area   = cmd->area;
    const int inarea1 = cmd->n_from1[0];       /**use only 1st inarea in list**/
    const int inarea2 = cmd->n_from2[0];       /**use only 1st inarea in list**/

    int d_a     = A[area].d_a;
    int d_b     = A[area].d_b;
    int d_a_in1 = A[inarea1].d_a;
    int d_b_in1 = A[inarea1].d_b;
    int d_a_in2 = A[inarea2].d_a;
    int d_b_in2 = A[inarea2].d_b;

    float *CMs_a = (float *) malloc (A[area].d_n * sizeof (float));
    float *CMs_b = (float *) malloc (A[area].d_n * sizeof (float));

    ALL_WEIGHTS_SIGPI *W = (ALL_WEIGHTS_SIGPI *)cmd->pointers[0]->data;

    for (int ct_n = 0; ct_n < A[area].d_n; ++ct_n) {

        CONN_SIGPI *conn = W->all_conn[area][ct_n][inarea1][inarea2];

        float CM_a = 0.0;
        float CM_b = 0.0;
        float norm = 0.0;
        int   num  = 0;

        if  (conn->next)
        do  {
            conn = conn->next;

            /**keep in1 const and integrate over in2**/
            if  (cmd->quantum[0][0] == 0) {
                if  ((conn->a1 == index_a) && (conn->b1 == index_b)) {

                    CM_a += conn->val * (float)(conn->a2);
                    CM_b += conn->val * (float)(conn->b2);
                    norm += conn->val;
                    num += 1;
                }

            } else { /**keep in2 const and integrate over in1**/

                if  ((conn->a2 == index_a) && (conn->b2 == index_b)) {

                    CM_a += conn->val * (float)(conn->a1);
                    CM_b += conn->val * (float)(conn->b1);
                    norm += conn->val;
                    num += 1;
                }
            }

        } while (conn->next != NULL);

        if  (norm != 0.0) {
            CMs_a[ct_n] = CM_a / norm;
            CMs_b[ct_n] = CM_b / norm;
        } else {
            CMs_a[ct_n] = -1.0;
            CMs_b[ct_n] = -1.0;
        }

        fprintf (stderr, "\nct_n=%d num=%d norm=%f CM_a=%f CM_b=%f ", ct_n, num, norm, CM_a, CM_b);
    }

    for (int ct_a = 0; ct_a < d_a; ++ct_a)
        for (int ct_b = 0; ct_b < d_b; ++ct_b)
            fprintf (stderr, "%f %f\n", CMs_a[ct_a * d_b + ct_b], CMs_b[ct_a * d_b + ct_b]);

    char fullname[512];
    sprintf (fullname, "%s/obs_%s_%d__%d_%d.gnu", cmd->pointers[1]->words[0], cmd->ch_pointers[0], area, index_a, index_b);    /**pointers[1] gives the directory!**/
    FILE *fp = fopen (fullname, "w");

    if  (cmd->quantum[0][0] == 0) {
        fprintf (fp, "set xrange [0:%d]\n", d_a_in2 - 1);
        fprintf (fp, "set yrange [0:%d]\n", d_b_in2 - 1);
    } else {
        fprintf (fp, "set xrange [0:%d]\n", d_a_in1 - 1);
        fprintf (fp, "set yrange [0:%d]\n", d_b_in1 - 1);
    }

    fprintf (fp, "set size 0.721,1.0\n");

    for (int ct_a = 0; ct_a < d_a; ++ct_a)
        for (int ct_b = 1; ct_b < d_b; ++ct_b)
            if  ((CMs_a[ct_a * d_b + ct_b - 1] != -1.0) && (CMs_b[ct_a * d_b + ct_b - 1] != -1.0) &&  (CMs_a[ct_a * d_b + ct_b] != -1.0) && (CMs_b[ct_a * d_b + ct_b] != -1.0))
                fprintf (fp, "set arrow from %f,%f to %f,%f nohead\n", CMs_a[ct_a * d_b + ct_b - 1], CMs_b[ct_a * d_b + ct_b - 1], CMs_a[ct_a * d_b + ct_b], CMs_b[ct_a * d_b + ct_b]);

    for (int ct_a = 1; ct_a < d_a; ++ct_a)
        for (int ct_b = 0; ct_b < d_b; ++ct_b)
            if  ((CMs_a[(ct_a - 1) * d_b + ct_b] != -1.0) && (CMs_b[(ct_a - 1) * d_b + ct_b] != -1.0) && (CMs_a[ct_a * d_b + ct_b] != -1.0) && (CMs_b[ct_a * d_b + ct_b] != -1.0))
                fprintf (fp, "set arrow from %f,%f to %f,%f nohead\n", CMs_a[(ct_a - 1) * d_b + ct_b], CMs_b[(ct_a - 1) * d_b + ct_b], CMs_a[ct_a * d_b + ct_b], CMs_b[ct_a * d_b + ct_b]);

    fprintf (fp, "set nokey; plot -0.1\n");

    fprintf (stderr, "\n\ngnuplot\nload \"%s\"\n\n", fullname);

    fclose (fp);
    return (DOUBLE)(0);
}

