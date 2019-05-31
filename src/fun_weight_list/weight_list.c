#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "../kernel/coco.h"
#include "../kernel/series.h"
#include "../kernel/utils.h"
#include "weight_list.h"




int count_connections (CONN *conn) {
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



void insert_conn (CONN *conn, DOUBLE newval, int a, int b, int n) {

    CONN *newconn = (CONN *) malloc (sizeof (CONN));
    CONN *overnext = conn->next;
    conn->next = newconn;
    newconn->next = overnext;

    newconn->val  = newval;
    newconn->a    = a;
    newconn->b    = b;
    newconn->n    = n;
}


/**used by:
   weight_list_alloc_full,
   weight_list_alloc_topo,
   weight_list_alloc_import;
   does the necessary init only once automatically**/

DOUBLE weight_list_alloc_init_once (PARAMS *g, AREA *A, COMMAND *cmd, int, int) {

    static int firsttime_weightsNotYetAllocated = 1; /**to allow several functions to use this even though it has to be done only once**/

    if  (cmd->pointers[0]->data == NULL) {

        if  (firsttime_weightsNotYetAllocated == 0)
	  fprintf (stderr, "\nweight_list_alloc_init_once warning: this should only occur if there are two different weight pointers (like w and v)\n");

        /**allocate pointer space**/
        ALL_WEIGHTS *W = (ALL_WEIGHTS *) malloc (sizeof (ALL_WEIGHTS));

        /**architecture**/
        W->areas = g->areas;

        W->all_conn = (CONN ****) malloc (W->areas * sizeof (CONN ***));

        W->d_a = (int *) malloc (W->areas * sizeof (int));
        W->d_b = (int *) malloc (W->areas * sizeof (int));
        W->d_n = (int *) malloc (W->areas * sizeof (int));

        /**area parameters**/
        for (int ct_ar = 0; ct_ar < W->areas; ++ct_ar) {

            W->d_a[ct_ar] = A[ct_ar].d_a;
            W->d_b[ct_ar] = A[ct_ar].d_b;
            W->d_n[ct_ar] = A[ct_ar].d_n;

            W->all_conn[ct_ar] = (CONN ***) malloc (W->d_n[ct_ar] * sizeof (CONN **));
        }

        /**neuron's pointers**/
        for (int ct_ar = 0; ct_ar < W->areas; ++ct_ar) {

            for (int ct_n = 0; ct_n < W->d_n[ct_ar]; ++ct_n) {

                W->all_conn[ct_ar][ct_n] = (CONN **) malloc (W->areas * sizeof (CONN *));

                /**init each neuron's pointer to all (input) areas to NULL**/
                for (int ct_ar_ar = 0; ct_ar_ar < W->areas; ++ct_ar_ar)

                    W->all_conn[ct_ar][ct_n][ct_ar_ar] = NULL;     /**for all areas, for all neurons, to every area  there exists a NULL pointer**/
            }
        }

        /**make the cmd pointer show to the right place (cmd->pointers[0] must NOT be lost!)**/
        cmd->pointers[0]->data = W;

        firsttime_weightsNotYetAllocated = 0;
    }

    return (DOUBLE)0;
}


/************************ weight_list_alloc_full *****************************/
/* Allocates a full connectivity between target- and all input areas.        */
/* q[0][0/1]=min/max of random initial connection values.                    */

DOUBLE weight_list_alloc_full (PARAMS *g, AREA *A, COMMAND *cmd, int, int) {

    if  (cmd->anz_quant[0] != 2)
        fprintf (stderr, "\nweight_list_alloc_full: insufficient parameters!");

    const DOUBLE upper = cmd->quantum[0][1];
    const DOUBLE lower = cmd->quantum[0][0];

    weight_list_alloc_init_once (g, A, cmd, 0, 0);

    ALL_WEIGHTS *W = (ALL_WEIGHTS *) cmd->pointers[0]->data;

    /**all input areas in command list**/
    for (int ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        const int inarea = cmd->n_from1[ct_l];                       /**area number in global list**/

        for (int ct_n = 0; ct_n < W->d_n[cmd->area]; ++ct_n) {

            W->all_conn[cmd->area][ct_n][inarea] = (CONN *) malloc (sizeof (CONN));           /**first connection to that inarea**/

            CONN *conn = W->all_conn[cmd->area][ct_n][inarea];

            /**create the admin element**/
            conn->val  = 0;
            conn->a    = -1;
            conn->b    = -1;
            conn->n    = -1;
            conn->next = NULL;

            /**create the linked list**/
            for (int ct_a = 0; ct_a < W->d_a[inarea]; ++ct_a)
            for (int ct_b = 0; ct_b < W->d_b[inarea]; ++ct_b) {

                if  (conn->next == NULL)
                    conn->next = (CONN *) malloc (sizeof (CONN));

                conn       = conn->next;

                conn->val  = lower + drand48() * (upper - lower);
                conn->a    = ct_a;
                conn->b    = ct_b;
                conn->n    = ct_a * W->d_b[inarea] + ct_b;
                conn->next = (conn->n < W->d_n[inarea] - 1) ? (CONN *) malloc (sizeof (CONN)) : NULL;
            }
        }
    }

    return (DOUBLE)(0);
}



/************************ weight_list_alloc_topo *****************************/
/* Allocates limited topographic connections between target/all input areas. */
/* q[0][0/1] = max/sigma as proportion of each dim of Gaussian init values.  */
/* q[1][0] =3: init on rgb input; =0/1: init normally.                       */
/* q[1][1] =2: init separate ON/OFF cells; then also only positive weights.  */
/* q[2][0] ="radius" as proportion of each dimension; useful is < 1.         */
/* Note that initialised fields will rather be ovals!                        */

DOUBLE weight_list_alloc_topo (PARAMS *g, AREA *A, COMMAND *cmd, int, int) {

    fprintf (stderr, "\nweight_list_alloc_topo ");

    if  (cmd->anz_quant[0] != 2)
        fprintf (stderr, "\nweight_list_alloc_topo: insufficient parameters!");

    if  (cmd->anz_quantums != 3)
        fprintf (stderr, "\nweight_list_alloc_topo: wrong number of parameters!");

    const DOUBLE upper = cmd->quantum[0][0];
    const DOUBLE sigma = cmd->quantum[0][1];

    weight_list_alloc_init_once (g, A, cmd, 0, 0);

    ALL_WEIGHTS *W = (ALL_WEIGHTS *) cmd->pointers[0]->data;

    /**all input areas in command list**/
    for (int ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        const int inarea = cmd->n_from1[ct_l];                       /**area number in global list**/

        const int d_a = W->d_a[cmd->area];
        const int d_b = W->d_b[cmd->area];

        for (int ct_a = 0; ct_a < d_a; ++ct_a)
        for (int ct_b = 0; ct_b < d_b; ++ct_b) {

            int ct_n = ct_a * d_b + ct_b;

            W->all_conn[cmd->area][ct_n][inarea] = (CONN *) malloc (sizeof (CONN));           /**first connection to that inarea**/

            CONN *conn = W->all_conn[cmd->area][ct_n][inarea];

            /**create the admin element**/
            conn->val  = 0;
            conn->a    = -1;
            conn->b    = -1;
            conn->n    = -1;
            conn->next = NULL;

            const int d_a_in = W->d_a[inarea];
            const int d_b_in = W->d_b[inarea];

            /**create the linked list**/
            for (int ct_a_in = 0; ct_a_in < d_a_in; ++ct_a_in)
            for (int ct_b_in = 0; ct_b_in < d_b_in; ++ct_b_in) {

                double diffA = (double)ct_a/(d_a-1.0) - (double)ct_a_in/(d_a_in-1.0);
                double diffB = (double)ct_b/(d_b-1.0) - (double)ct_b_in/(d_b_in-1.0);

int toroid = 1;
                if  (toroid) {
                    diffA = (diffA > 0.5) ? (1.0 - diffA) : diffA;
                    diffB = (diffB > 0.5) ? (1.0 - diffB) : diffB;
                    diffA = (diffA < -0.5) ? (-1.0 - diffA) : diffA;
                    diffB = (diffB < -0.5) ? (-1.0 - diffB) : diffB;
                    if  ((diffA < -1.0) || (diffB < -1.0) || (diffA > 1.0) || (diffB > 1.0))
                        fprintf (stderr, "\n\n\n weight_list_alloc_topo: re-consider the toroid part !! \n\n\n");
                }

                if  (cmd->quantum[1][0] == 3) {

                    const int d_a_in_eff = d_a_in / 3;

                    if  (ct_a_in < d_a_in_eff)                                         /**red**/
                        diffA = (double)ct_a/(d_a-1.0) - (double)ct_a_in/(d_a_in_eff-1.0);
                    if  ((ct_a_in >= d_a_in_eff) && (ct_a_in < 2 * d_a_in_eff))        /**green**/
                        diffA = (double)ct_a/(d_a-1.0) - (double)(ct_a_in - d_a_in_eff)/(d_a_in_eff-1.0);
                    if  (ct_a_in >= 2 * d_a_in_eff)                                    /**blue**/
                        diffA = (double)ct_a/(d_a-1.0) - (double)(ct_a_in - 2 * d_a_in_eff)/(d_a_in_eff-1.0);
                }

                if  (cmd->quantum[1][1] == 2) {

                    const int d_b_in_eff = d_b_in / 2;

                    if  (ct_b_in < d_b_in_eff)                                         /**ON center cells**/
                        diffB = (double)ct_b/(d_b-1.0) - (double)ct_b_in/(d_b_in_eff-1.0);
                    if  (ct_b_in >= d_b_in_eff)                                        /**OFF center cells**/
                        diffB = (double)ct_b/(d_b-1.0) - (double)(ct_b_in - d_b_in_eff)/(d_b_in_eff-1.0);
                }

                double dist_sq = diffA * diffA + diffB * diffB;

                if  (dist_sq < cmd->quantum[2][0] * cmd->quantum[2][0]) {

                    if  (conn->next == NULL)
                        conn->next = (CONN *) malloc (sizeof (CONN));

                    conn       = conn->next;

                    conn->val  = drand48 () * upper * exp (-0.5 * dist_sq / (sigma * sigma));
                    if  (cmd->quantum[1][1] != 2)
                        conn->val  = drand48() > 0.5 ? -conn->val : conn->val;
                    conn->a    = ct_a_in;
                    conn->b    = ct_b_in;
                    conn->n    = ct_a_in * W->d_b[inarea] + ct_b_in;
                    conn->next = NULL;
                }
            }
        }
    }

    return (DOUBLE)(0);
}



/************************ weight_list_alloc_invert ***************************/
/* Allocates connections from input area where connections to it must exist. */
/* For intra-area (lateral) weights, use different name & 2nd pointer!       */

DOUBLE weight_list_alloc_invert (PARAMS *g, AREA *A, COMMAND *cmd, int, int) {

fprintf (stderr, "\nweight_list_alloc_invert ... ");

    ALL_WEIGHTS *W_in;

    if  (cmd->anz_pointers == 2) {
                 W_in = (ALL_WEIGHTS *) cmd->pointers[1]->data;
                 weight_list_alloc_init_once (g, A, cmd, 0, 0);  /**this is to allocate the cmd->pointers[0]->data weights, if first time use of the letter (e.g. "v")**/
    } else {
                 W_in = (ALL_WEIGHTS *) cmd->pointers[0]->data;
    }

    ALL_WEIGHTS *W_out = (ALL_WEIGHTS *) cmd->pointers[0]->data;

    if  (cmd->anz_from1 != 1)
        fprintf (stderr, "\nweight_list_alloc_invert wants only one input area!\n");

    const int inarea = cmd->n_from1[0];               /**just one input area**/

    /**on this area, intialise each neuron's list**/
    for (int ct_n = 0; ct_n < W_out->d_n[cmd->area]; ++ct_n) {

            W_out->all_conn[cmd->area][ct_n][inarea] = (CONN *) malloc (sizeof (CONN));          /**new first connection to that inarea**/

            CONN *conn = W_out->all_conn[cmd->area][ct_n][inarea];

            /**create the admin element**/
            conn->val  = 0;
            conn->a    = -1;
            conn->b    = -1;
            conn->n    = -1;
            conn->next = NULL;
    }

    /**loop through inarea**/
    for (int ct_a_in = 0; ct_a_in < W_in->d_a[inarea]; ++ct_a_in)
    for (int ct_b_in = 0; ct_b_in < W_in->d_b[inarea]; ++ct_b_in) {

        int ct_n_in = ct_a_in * W_in->d_b[inarea] + ct_b_in;

        CONN *conn_in = W_in->all_conn[inarea][ct_n_in][cmd->area];            /**Note that inarea & area are reversed here!**/

        /**loop through connections of the inarea neuron!**/
        if  (conn_in->next)
        do  {
            conn_in = conn_in->next;

            /**find the target neuron "conn_in->n" at this area**/
            CONN *conn = W_out->all_conn[cmd->area][conn_in->n][inarea];

            /**go to list end at this area neuron**/
            if  (conn->next)
            while (conn->next != NULL)  {
                conn = conn->next;
            }

            /**add a connection to the end of the list of this area neuron**/
            conn->next = (CONN *) malloc (sizeof (CONN));

            conn->next->val  = conn_in->val;
            conn->next->a    = ct_a_in;
            conn->next->b    = ct_b_in;
            conn->next->n    = ct_n_in;
            conn->next->next = NULL;

        } while (conn_in->next != NULL);
    }

fprintf (stderr, "... end  ");

    return (DOUBLE)(0);
}


/************************ weight_list_invert_quick ***************************/
/* Assumes full bidirectional connections, and inverts via static array.     */
/* Use this function only for one size in program!                           */
/* Easy to work around full connectivity!                                    */

DOUBLE weight_list_invert_quick (PARAMS *g, AREA *A, COMMAND *cmd, int, int) {

    /*fprintf (stderr, "\nweight_list_invert_quick ... ");*/

    /**only one input area allowed**/
    const int inarea = cmd->n_from1[0];                           /**area number in global list**/

    ALL_WEIGHTS *W = (ALL_WEIGHTS *) cmd->pointers[0]->data;

    static DOUBLE **Wtemp;
    static int    **Wtemp_i;
    static int firsttime = 1;

    if  (firsttime) {
        Wtemp     = d_matrix (A[cmd->area].d_n, A[inarea].d_n);
        Wtemp_i   = i_matrix (A[cmd->area].d_n, A[inarea].d_n);
        firsttime = 0;
    }

    /**loop through inarea to fill matrix**/
    for (int ct_n_in = 0; ct_n_in < A[inarea].d_n; ++ct_n_in) {

        CONN *conn_in = W->all_conn[inarea][ct_n_in][cmd->area];            /**Note that inarea & area are reversed here!**/

        /**loop through connections of the inarea neuron!**/
        if  (conn_in != NULL)
        if  (conn_in->next)
        do  {
            conn_in = conn_in->next;
            Wtemp[conn_in->n][ct_n_in] = conn_in->val;
            Wtemp_i[conn_in->n][ct_n_in] = 1;

        } while (conn_in->next != NULL);
    }


    /**loop through target area to read out matrix**/
    for (int ct_n = 0; ct_n < A[cmd->area].d_n; ++ct_n) {

        CONN *conn = W->all_conn[cmd->area][ct_n][inarea];               /**[area][ct_n][glob inarea](first inarea neuron)**/

        /**loop over input weights**/
        if  (conn != NULL)
        if  (conn->next)
        do  {
            conn = conn->next;
            conn->val = Wtemp[ct_n][conn->n];

            if  (Wtemp_i[ct_n][conn->n] != 1)
                fprintf (stderr, "\nWarning: inconsistency in weight_list_invert_quick!\n");

        } while (conn->next != NULL);
    }


    /*fprintf (stderr, "... end  ");*/
    return (DOUBLE)(0);
}



/**************************** weight_list_mult *******************************/
/* Use like:  2(n) {N,w; weight_list_mult; 0, 1; , ;   }                     */

DOUBLE weight_list_mult (PARAMS *g, AREA *A, COMMAND *cmd, int, int) {

    fprintf (stderr, "\nweight_list_mult ... ");

    const int area_low = cmd->n_from1[0];     /**overall input**/
    const int area_mid = cmd->n_from2[0];     /**middle area**/
    const int area_hig = cmd->area;           /**overall output**/

    const int d_n_low = A[area_low].d_n;
    const int d_n_mid = A[area_mid].d_n;
    const int d_n_hig = A[area_hig].d_n;

    ALL_WEIGHTS *W = (ALL_WEIGHTS *) cmd->pointers[0]->data;

    static DOUBLE **W_right, **W_left, **W_tot;
    static int    firsttime = 1;

    if  (firsttime) {
        W_right   = d_matrix (d_n_mid, d_n_low);
        W_left    = d_matrix (d_n_hig, d_n_mid);
        W_tot     = d_matrix (d_n_hig, d_n_low);
        firsttime = 0;
    }

    int ct_left_zero = 0, ct_left_vals = 0;
    int ct_right_zero = 0, ct_right_vals = 0;
    int ct_tot_zero = 0, ct_tot_vals = 0;

    /**init matrices zero**/
    for (int j = 0; j < d_n_mid; ++j)
        for (int k = 0; k < d_n_low; ++k)
            W_right[j][k] = 0.0;
    for (int i = 0; i < d_n_hig; ++i)
        for (int j = 0; j < d_n_mid; ++j)
            W_left[i][j] = 0.0;
    for (int i = 0; i < d_n_hig; ++i)
        for (int k = 0; k < d_n_low; ++k)
            W_tot[i][k] = 0.0;

    /**fill right matrix**/
    for (int ct_n = 0; ct_n < d_n_mid; ++ct_n) {

        CONN *conn = W->all_conn[area_mid][ct_n][area_low];

        if  (conn != NULL)
        if  (conn->next)
        do  {
            conn = conn->next;
            W_right[ct_n][conn->n] = conn->val;

            if  (conn->val == 0.0)
                ct_right_zero += 1;
            else
                ct_right_vals += 1;

        } while (conn->next != NULL);
    }

    fprintf (stderr, "\nright: %d x zero, %d x nonzero, %d exist ",
             ct_right_zero, ct_right_vals, d_n_mid * d_n_low);

    /**fill left matrix**/
    for (int ct_n = 0; ct_n < d_n_hig; ++ct_n) {

        CONN *conn = W->all_conn[area_hig][ct_n][area_mid];

        if  (conn != NULL)
        if  (conn->next)
        do  {
            conn = conn->next;
            W_left[ct_n][conn->n] = conn->val;

            if  (conn->val == 0.0)
                ct_left_zero += 1;
            else
                ct_left_vals += 1;

        } while (conn->next != NULL);
    }

    fprintf (stderr, "\nleft: %d x zero, %d x nonzero, %d exist ",
             ct_left_zero, ct_left_vals, d_n_hig * d_n_mid);

    /**do the multiplication**/
    for (int i = 0; i < d_n_hig; ++i)
        for (int k = 0; k < d_n_low; ++k)
            for (int j = 0; j < d_n_mid; ++j)
                W_tot[i][k] += W_left[i][j] * W_right[j][k];

    for (int i = 0; i < d_n_hig; ++i)
        for (int k = 0; k < d_n_low; ++k)
            if  (W_tot[i][k] == 0.0)
                ct_tot_zero += 1;
            else
                ct_tot_vals += 1;

    fprintf (stderr, "\ntot: %d x zero, %d x nonzero, %d exist ",
             ct_tot_zero, ct_tot_vals, d_n_hig * d_n_low);


    /**write to list**/
    for (int ct_n = 0; ct_n < d_n_hig; ++ct_n) {

        CONN *conn = W->all_conn[area_hig][ct_n][area_low];

        if  (conn != NULL)
        if  (conn->next)
        do  {
            conn = conn->next;
            conn->val = W_tot[ct_n][conn->n];

        } while (conn->next != NULL);
    }

    fprintf (stderr, "... end\n");
    return (DOUBLE)(0);
}



/************************ weight_list_free ***********************************/
/* De-allocates all connectivity between target- and all input areas.        */
/* Useful e.g. so that weight_list_alloc_invert can be called repeatedly ... */

DOUBLE weight_list_free (PARAMS *g, AREA *A, COMMAND *cmd, int, int) {

fprintf (stderr, "\nweight_list_free ... ");

    ALL_WEIGHTS *W = (ALL_WEIGHTS *) cmd->pointers[0]->data;

    /**all input areas in command list**/
    for (int ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        const int inarea = cmd->n_from1[ct_l];                       /**area number in global list**/

        for (int ct_n = 0; ct_n < W->d_n[cmd->area]; ++ct_n) {

            CONN *conn = W->all_conn[cmd->area][ct_n][inarea];               /**[area][ct_n][glob inarea](first inarea neuron)**/
            CONN *conn_old;

            if  (conn != NULL) {

                if  (conn->next)
                do  {

                    conn_old = conn;
                    conn = conn->next;
                    free (conn_old);

                } while (conn->next != NULL);

                free (conn);                                                /**last element**/

                W->all_conn[cmd->area][ct_n][inarea] = NULL;                /**first connection to that inarea (no admin element!)**/
            }
        }
    }

fprintf (stderr, "... end  ");

    return (DOUBLE)(0);
}



/**************************** weight_list_feed *******************************/

DOUBLE weight_list_feed (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    ALL_WEIGHTS *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;

    DOUBLE act = (DOUBLE)(0);

    int intime = (ct_t + (int)(cmd->quantum[0][0]) < 0)
               ? 0
               : ct_t + (int)(cmd->quantum[0][0]);


    /**all input areas in command list**/
    for (int ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        const int inarea = cmd->n_from1[ct_l];                           /**area number in global list**/

        CONN *conn = W->all_conn[cmd->area][ct_n][inarea];               /**[area][ct_n][glob inarea](first inarea neuron)**/

        DOUBLE *input = cmd->S_from1[ct_l][intime];                        /**[cmd inarea list][ct_t](inarea neurons)**/

        /**sum over weighted inputs**/
        if  (conn != NULL)
        if  (conn->next)
        do  {
            conn = conn->next;
            act += conn->val * input[conn->n];

        } while (conn->next != NULL);
    }

    return act;
}



/**************************** WEIGHT_LIST_PUSH *******************************/
/* WRITES TO cmd->S_from1 !! USE CAREFULLY !! Prefer weight_list_AXON_total !*/
/* NOTE: output should have been initialised to zero before!                 */
/* q[1][0]=act threshold above which to push out.                            */

DOUBLE WEIGHT_LIST_PUSH (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    ALL_WEIGHTS *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;

    int intime = (ct_t + (int)(cmd->quantum[0][0]) < 0)
               ? 0
               : ct_t + (int)(cmd->quantum[0][0]);

    static int firsttime = 1;
    if  (firsttime) {
        fprintf (stderr, "\n\nWarning: you are using WEIGHT_LIST_PUSH which OVERWRITES S_from1 !\n\n");
        firsttime = 0;
    }

    /**do the feed_out only if the neuron's activation is significant**/
    if  (fabs (cmd->S_target[intime][ct_n]) >= cmd->quantum[1][0]) {

    /**all input areas in command list**/
    for (int ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        const int inarea = cmd->n_from1[ct_l];

        CONN *conn = W->all_conn[cmd->area][ct_n][inarea];               /**same as in weight_list_feed**/

        DOUBLE *_output_ = cmd->S_from1[ct_l][ct_t];                     /**will WRITE to S_from1 at current time!**/

        /**write weighted values to connected _output_s**/
        if  (conn->next)
        do  {
            conn = conn->next;
            _output_[conn->n] += conn->val * cmd->S_target[intime][ct_n];

        } while (conn->next != NULL);
    }

    }

    return cmd->S_target[ct_t][ct_n]; /**leave as is**/
}



/**************************** weight_list_AXON_total *************************/
/* Feeds using weights as outgoing weights! conn-pointer is switched!        */
/* Thus hebb must be inverted and decay is also the other way around!        */
/* q[0][0]=intime relative to ct_t.                                          */
/* q[1][0]=1: see following act threshold absolute; =2: as relative to max.  */
/* q[2][0]=act threshold above which to push out.                            */

DOUBLE weight_list_AXON_total (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    ALL_WEIGHTS *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;

    int intime = (ct_t + (int)(cmd->quantum[0][0]) < 0)
               ?  0
               :  ct_t + (int)(cmd->quantum[0][0]);

    /**init target to zero**/
    for (int ct_n = 0; ct_n < A[cmd->area].d_n; ++ct_n)
        cmd->S_target[ct_t][ct_n] = 0.0;

    /**all input areas in command list**/
    for (int ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        const int inarea = cmd->n_from1[ct_l];
        DOUBLE max = 0.0;

        /**get max act**/
        if  (cmd->quantum[1][0] == 2)
            for (int ct_in = 0; ct_in < A[inarea].d_n; ++ct_in)
                max = fabs (cmd->S_from1[ct_l][intime][ct_in]) > max ? fabs (cmd->S_from1[ct_l][intime][ct_in]) : max;


        for (int ct_in = 0; ct_in < A[inarea].d_n; ++ct_in) {

            int significant = 1;

            /**do the feed_out only if the neuron's activation is significant**/
            if  (cmd->quantum[1][0] == 1)
                if  (fabs (cmd->S_from1[ct_l][intime][ct_in]) < cmd->quantum[2][0])
                    significant = 0;

            if  (cmd->quantum[1][0] == 2)
                if  (fabs (cmd->S_from1[ct_l][intime][ct_in]) < cmd->quantum[2][0] * max)
                    significant = 0;

            if  (significant) {

                CONN *conn = W->all_conn[inarea][ct_in][cmd->area];               /**AXONS! This is switched! (because we use the list of the inarea neurons)**/

                /**write weighted values to connected _output_s**/
                if  (conn->next)
                do  {
                    conn = conn->next;
                    cmd->S_target[ct_t][conn->n] += conn->val * cmd->S_from1[ct_l][intime][ct_in];

                } while (conn->next != NULL);
            }
        }
    }
    return (DOUBLE)(0);
}



/**************************** weight_list_AXON_total_pos *********************/
/* Like weight_list_AXON_total but feeds only via pos or only neg weights.   */
/* q[3][0]=1: only via positive weights; =-1: only via negative weights.     */

DOUBLE weight_list_AXON_total_pos (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    ALL_WEIGHTS *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;

    int intime = (ct_t + (int)(cmd->quantum[0][0]) < 0)
               ?  0
               :  ct_t + (int)(cmd->quantum[0][0]);

    if  (cmd->anz_quantums != 4)
        fprintf (stderr, "\n\nwrong use 1 of weight_list_AXON_total_pos\n\n");
    if  ((cmd->quantum[3][0] != 1) && (cmd->quantum[3][0] != -1))
        fprintf (stderr, "\n\nwrong use 2 of weight_list_AXON_total_pos\n\n");

    /**init target to zero**/
    for (int ct_n = 0; ct_n < A[cmd->area].d_n; ++ct_n)
        cmd->S_target[ct_t][ct_n] = 0.0;

    /**all input areas in command list**/
    for (int ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        const int inarea = cmd->n_from1[ct_l];
        DOUBLE max = 0.0;

        /**get max act**/
        if  (cmd->quantum[1][0] == 2)
            for (int ct_in = 0; ct_in < A[inarea].d_n; ++ct_in)
                max = fabs (cmd->S_from1[ct_l][intime][ct_in]) > max ? fabs (cmd->S_from1[ct_l][intime][ct_in]) : max;


        for (int ct_in = 0; ct_in < A[inarea].d_n; ++ct_in) {

            int significant = 1;

            /**do the feed_out only if the neuron's activation is significant**/
            if  (cmd->quantum[1][0] == 1)
                if  (fabs (cmd->S_from1[ct_l][intime][ct_in]) < cmd->quantum[2][0])
                    significant = 0;

            if  (cmd->quantum[1][0] == 2)
                if  (fabs (cmd->S_from1[ct_l][intime][ct_in]) < cmd->quantum[2][0] * max)
                    significant = 0;

            if  (significant) {

                CONN *conn = W->all_conn[inarea][ct_in][cmd->area];               /**AXONS! This is switched! (because we use the list of the inarea neurons)**/

                /**write weighted values to connected _output_s**/
                if  (conn->next)
                do  {
                    conn = conn->next;

                    if  (cmd->quantum[3][0] == 1) {
                        if  (conn->val > 0.0)
                            cmd->S_target[ct_t][conn->n] += conn->val * cmd->S_from1[ct_l][intime][ct_in];
                    } else {
                        if  (conn->val < 0.0)
                            cmd->S_target[ct_t][conn->n] += conn->val * cmd->S_from1[ct_l][intime][ct_in];
                    }

                } while (conn->next != NULL);
            }
        }
    }
    return (DOUBLE)(0);
}



/**************************** weight_list_hebb *******************************/
/* Works on every neuron but no return val. Updates only if act[ct_n] != 0.0.*/
/* cmd->pointers[0] is a ALL_WEIGHTS *. Arguments *g and *A are not used.    */

DOUBLE weight_list_hebb (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    ALL_WEIGHTS *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;

    const DOUBLE post = cmd->S_from1[0][ct_t][ct_n];                  /**[ct_l="here"][ct_t][ct_n]**/

    //static DOUBLE average = 0.0;
    //if  (ct_n == 0)
    //    average = 0.0;

    if  (post != 0.0) {

        /**check only**/
        if  ((cmd->anz_from1 != 1) || (cmd->n_from1[0] != cmd->area))
            fprintf (stderr, "wrong use of weight_list_hebb");

        /**all input areas in command list**/
        for (int ct_l = 0; ct_l < cmd->anz_from2; ++ct_l) {

            const int inarea = cmd->n_from2[ct_l];                    /**area number in global list**/

            CONN *conn = W->all_conn[cmd->area][ct_n][inarea];        /**[area][ct_n][glob inarea](inarea admin neuron)**/

            DOUBLE *pre = cmd->S_from2[ct_l][ct_t];                   /**[cmd inarea list][ct_t](inarea neurons)**/

            /**input connections**/
            if  (conn->next)
            do  {
                conn = conn->next;
                conn->val += cmd->moment * cmd->quantum[0][0] * post * pre[conn->n];

                //average += cmd->moment * cmd->quantum[0][0] * post * pre[conn->n];

            } while (conn->next != NULL);
        }
    }

    //if  (ct_t == 0)
    //if  (ct_n == A[cmd->area].d_n - 1)
    //    fprintf (stdout, " %c:av=%f ", cmd->ch_target, average);

    return (DOUBLE)(0);
}



/**************************** weight_list_kohonen ****************************/
/* Like weight_hebb, but pre[] replaced by (pre[conn->n] - conn->val).       */
/* Expects neighborhood function as post.                                    */

DOUBLE weight_list_kohonen (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    ALL_WEIGHTS *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;

    const DOUBLE post = cmd->S_from1[0][ct_t][ct_n];                     /**[ct_l="here"][ct_t][ct_n]**/

    if  (post != 0.0) {

        /**check only**/
        if  ((cmd->anz_from1 != 1) || (cmd->n_from1[0] != cmd->area))
            fprintf (stderr, "wrong use of weight_list_kohonen");

        /**all input areas in command list**/
        for (int ct_l = 0; ct_l < cmd->anz_from2; ++ct_l) {

            const int inarea = cmd->n_from2[ct_l];                           /**area number in global list**/

            CONN *conn = W->all_conn[cmd->area][ct_n][inarea];               /**[area][ct_n][glob inarea](inarea admin neuron)**/

            DOUBLE *pre = cmd->S_from2[ct_l][ct_t];                          /**[cmd inarea list][ct_t](inarea neurons)**/

            /**input connections**/
            if  (conn->next)
            do  {
                conn = conn->next;
                conn->val += cmd->moment * cmd->quantum[0][0] * post * (pre[conn->n] - conn->val);

            } while (conn->next != NULL);
        }
    }

    return (DOUBLE)(0);
}



/**************************** weight_list_euclid *****************************/
/* Returns Euclidian distance of neuron ct_n's weight vector to data point.  */
/* Data and weights in from2 list; Depends on time ct_t.                     */

DOUBLE weight_list_euclid (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    double dist = 0.0;

    ALL_WEIGHTS *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;

    {
        /**check only**/
        if  ((cmd->anz_from1 != 1) || (cmd->n_from1[0] != cmd->area))
            fprintf (stderr, "wrong use of weight_list_euclid");

        /**all input areas in command list**/
        for (int ct_l = 0; ct_l < cmd->anz_from2; ++ct_l) {

            const int inarea = cmd->n_from2[ct_l];                           /**area number in global list**/

            CONN *conn = W->all_conn[cmd->area][ct_n][inarea];               /**[area][ct_n][glob inarea](inarea admin neuron)**/

            DOUBLE *pre = cmd->S_from2[ct_l][ct_t];                          /**[cmd inarea list][ct_t](inarea neurons)**/

            /**input connections**/
            if  (conn->next)
            do  {
                conn = conn->next;
                dist += (pre[conn->n] - conn->val) * (pre[conn->n] - conn->val);

            } while (conn->next != NULL);
        }
    }

    return (sqrt (dist));
}


/**************************** weight_list_euclid_ptr *************************/
/* Returns Euclidian dist of neuron ptr[1]->int_val's weight vector to data. */
/* Weights are between area S_from2 -> S_from1; Data point is in S_from2.    */
/* Note that output area should have only 1 neuron. Use as single-function.  */
/* q[0][0]=2 then take fabs(conn->val), e.g. to fit pos&neg RF with Gaussian */
/* Very unusual function .. don't copy from here!                            */

DOUBLE weight_list_euclid_ptr (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    double dist = 0.0;

    if  (ct_n != 0)
        fprintf (stderr, "\nweight_list_euclid_ptr: target area must have 1 unit!");

    ALL_WEIGHTS *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;

    {
        /**no command list; just one input area for simplicity**/

            const int inarea = cmd->n_from2[0];
            const int area   = cmd->n_from1[0]; /**! "target area of weights != target area of function!**/

            CONN *conn = W->all_conn[area][cmd->pointers[1]->int_val][inarea];

            DOUBLE *pre = cmd->S_from2[0][ct_t];

            /**input connections**/
            if  (conn->next)
            do  {
                conn = conn->next;

                if  (cmd->quantum[0][0] == 2)
                    dist += (pre[conn->n] - fabs (conn->val)) * (pre[conn->n] - fabs (conn->val));
                else
                    dist += (pre[conn->n] - conn->val) * (pre[conn->n] - conn->val);

            } while (conn->next != NULL);
    }

    return (sqrt (dist));
}


/**************************** weight_list_act2weight *************************/
/* Sets weights of neuron ct_n to activations in S_from1.                    */
/* For one selected neuron, use in a sweep like: sw (0; 1; $ct_n)            */

DOUBLE weight_list_act2weight (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    const int inarea = cmd->n_from1[0];

    ALL_WEIGHTS *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;

    /**check only**/
    if  (cmd->anz_from1 != 1)
        fprintf (stderr, "wrong use of weight_list_act2weight");

    CONN *conn = W->all_conn[cmd->area][ct_n][inarea];               /**[area][ct_n][glob inarea](inarea admin neuron)**/

    DOUBLE *pre = cmd->S_from1[0][ct_t];                          /**[cmd inarea list][ct_t](inarea neurons)**/

    /**input connections**/
    if  (conn->next)
        do  {
            conn = conn->next;
            conn->val = pre[conn->n];

        } while (conn->next != NULL);

    return (0.0);
}


/**************************** weight_list_act2weight_offset ******************/
/* Sets weights to S_from1: different neurons get acts dependent on location.*/
/* Assumes: inarea >= area. S_from1's coordinate (0,0) is placed onto center-*/
/* -field of RF (which is within S_from1) that has the size of area.         */
/* Toroidal treatment of S_from1 (assumes filter kernel, Gauss, etc.).       */
/* Tested so far ONLY for circularly symmetric S_from1, around (0,0)!        */
/* q00=0: offset act-vector faithfully by 1 unit in center area with frame   */
/* q00=N>0: inarea_size = N*area_size+1; centers spread N units; 1 pix frame */
/* q10 exists then non-periodic boundary .. cutoff of filter at 1/2 area dist*/

DOUBLE weight_list_act2weight_offset (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    fprintf (stderr, "\nweight_list_act2weight_offset ... ");

    const int area    = cmd->area;
    const int inarea  = cmd->n_from1[0];
    const int frame_a = (A[inarea].d_a - A[area].d_a) / 2;
    const int frame_b = (A[inarea].d_b - A[area].d_b) / 2;
    DOUBLE    *vec    = d_vector (A[inarea].d_n);

    ALL_WEIGHTS *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;

    /**check**/
    if  (cmd->anz_from1 != 1)
        fprintf (stderr, "wrong use of weight_list_act2weight");

    int N = (int)(cmd->quantum[0][0]);

    if  (cmd->quantum[0][0] != 0)
        if  ((A[inarea].d_a != A[area].d_a * N + 1) || (A[inarea].d_b != A[area].d_b * N + 1))
            fprintf (stderr, "\n\n\n\nweight_list_act2weight_offset: area size mismatch!\n\n\n\n");

    for (int ct_a = 0; ct_a < A[area].d_a; ++ct_a)
    for (int ct_b = 0; ct_b < A[area].d_b; ++ct_b) {

        CONN   *conn = W->all_conn[cmd->area][ct_a * A[area].d_b + ct_b][inarea];               /**[area][ct_n][glob inarea](inarea admin neuron)**/

        int at_a, at_b;

        if  (cmd->quantum[0][0] == 0) {
            at_a = ct_a + frame_a;
            at_b = ct_b + frame_b;
        } else {
            at_a = 1 + ct_a * N;
            at_b = 1 + ct_b * N;
        }

        /**check whether out of boundaries**/
        if  (at_a < 0) {
            fprintf (stderr, "\nweight_list_act2weight_offset: at_a=%d, set to zero  ", at_a);
            at_a = 0;
        }
        if  (at_b < 0) {
            fprintf (stderr, "\nweight_list_act2weight_offset: at_b=%d, set to zero  ", at_b);
            at_b = 0;
        }

        for (int X = 0; X < A[inarea].d_a; ++X) {
            for (int Y = 0; Y < A[inarea].d_b; ++Y) {

                int in_a = (X-at_a);
                if  (in_a < 0)
                    in_a += A[inarea].d_a;
                if  (in_a > A[inarea].d_a / 2)
                    in_a = A[inarea].d_a - in_a;

                int in_b = (Y-at_b);
                if  (in_b < 0)  in_b += A[inarea].d_b;
                if  (in_b > A[inarea].d_b / 2)
                    in_b = A[inarea].d_b - in_b;

                vec[X * A[inarea].d_b + Y] = cmd->S_from1[0][ct_t][in_a * A[inarea].d_b + in_b];

                if  (cmd->anz_quantums > 1) {
                    float dist = sqrt ((X-at_a)*(X-at_a)+(Y-at_b)*(Y-at_b));
                    if  (dist > (A[inarea].d_a + A[inarea].d_b) / 4.0)
                        vec[X * A[inarea].d_b + Y] = 0.0;
                }
            }
        }

        /**connections**/
        if  (conn->next)
        do  {
            conn = conn->next;
            conn->val = vec[conn->a * A[inarea].d_b + conn->b];

        } while (conn->next != NULL);
    }

    fprintf (stderr, " ... end   ");

    return (0.0);
}



/**************************** weight_list_hebb_turn **************************/
/* Like hebb, but leaves weight length as is. Returns squared weight length. */

DOUBLE weight_list_hebb_turn (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    ALL_WEIGHTS *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;

    const DOUBLE post = cmd->S_from1[0][ct_t][ct_n];                     /**[ct_l="here"][ct_t][ct_n]**/

    DOUBLE quad_length_before = (DOUBLE)0;
    DOUBLE quad_length_after  = (DOUBLE)0;

    //if  (post != 0.0) {

        /**check only**/
        if  ((cmd->anz_from1 != 1) || (cmd->n_from1[0] != cmd->area))
            fprintf (stderr, "wrong use of weight_hebb");

        /**all input areas in command list**/
        for (int ct_l = 0; ct_l < cmd->anz_from2; ++ct_l) {

            const int inarea = cmd->n_from2[ct_l];                           /**area number in global list**/

            CONN *conn = W->all_conn[cmd->area][ct_n][inarea];               /**[area][ct_n][glob inarea](inarea admin neuron)**/

            DOUBLE *pre = cmd->S_from2[ct_l][ct_t];                          /**[cmd inarea list][ct_t](inarea neurons)**/

            /**input connections**/
            if  (conn->next)
            do  {
                conn = conn->next;

                quad_length_before += conn->val * conn->val;
                conn->val += cmd->moment * cmd->quantum[0][0] * post * pre[conn->n];
                quad_length_after  += conn->val * conn->val;

            } while (conn->next != NULL);
        }

        DOUBLE factor = sqrt (quad_length_before) / sqrt (quad_length_after);

        /**all input areas in command list**/
        for (int ct_l = 0; ct_l < cmd->anz_from2; ++ct_l) {

            const int inarea = cmd->n_from2[ct_l];                           /**area number in global list**/

            CONN *conn = W->all_conn[cmd->area][ct_n][inarea];               /**[area][ct_n][glob inarea](inarea admin neuron)**/

            if  (conn->next)
            do  {
                conn = conn->next;

                conn->val *= factor;

            } while (conn->next != NULL);
        }
    //}

    return quad_length_before;
}



/**************************** weight_list_hebb_diff **************************/
/* Assumes Hebb-like use in Olshausen rule with pre=data-reconstruction.     */
/* Because only OTHER neuron's reconstruction should count,                  */
/* here this (ct_n's) neuron's reconstruction is added again.                */
/* cmd->pointers[1] used, because previous reconstruction divided by this val*/

DOUBLE weight_list_hebb_diff (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    ALL_WEIGHTS *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;

    const DOUBLE post = cmd->S_from1[0][ct_t][ct_n];                  /**[ct_l="here"][ct_t][ct_n]**/

    if  (post != 0.0) {

        /**check only**/
        if  ((cmd->anz_from1 != 1) || (cmd->n_from1[0] != cmd->area))
            fprintf (stderr, "wrong use of weight_list_hebb");

        /**all input areas in command list**/
        for (int ct_l = 0; ct_l < cmd->anz_from2; ++ct_l) {

            const int inarea = cmd->n_from2[ct_l];                    /**area number in global list**/

            CONN *conn = W->all_conn[cmd->area][ct_n][inarea];        /**[area][ct_n][glob inarea](inarea admin neuron)**/

            DOUBLE *pre = cmd->S_from2[ct_l][ct_t];                   /**[cmd inarea list][ct_t](inarea neurons)**/

            /**input connections**/
            if  (conn->next)
            do  {
                conn = conn->next;
                conn->val += cmd->moment * cmd->quantum[0][0] * post * (pre[conn->n]  +  post * conn->val / cmd->pointers[1]->float_val);

            } while (conn->next != NULL);
        }
    }

fprintf (stderr, "  backproj divided by %f     ", cmd->pointers[1]->float_val);

    return (DOUBLE)(0);
}


/****************************** weight_list_ica1 ****************************/
/* Formerly total_weight_ica1.                                               */
/* First part of the covariant ICA learning rule: Just add W ("anti-decay"). */

DOUBLE weight_list_ica1 (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    static int firsttime = 1;
    static int d_r_local;
    static int d_r_input;

//    double **W, **dW;
    ALL_WEIGHTS *W = (ALL_WEIGHTS *) cmd->pointers[0]->data;
    CONN *conn;

    if  (dummy != 0)
      fprintf (stderr, "\n\nweight_list_ica1 to be used in a total stay!\n");

    const int inarea   = cmd->n_from2[0];
    const int area     = cmd->area;

    DOUBLE eps  = cmd->quantum[0][0];

    if  (firsttime) {
        d_r_local = A[area].d_a * A[area].d_b;
        d_r_input = A[inarea].d_n;
        firsttime = 0;
    }

    if  (d_r_local != A[area].d_n)
        fprintf (stderr, "\ncan't use weight_ica1 twice\n");
    if  (d_r_local > A[inarea].d_n)
        fprintf (stderr, "\nno overcompleteness for weight_ica1, please!\n");
    if  (cmd->anz_from1 != 1)
        fprintf (stderr, "\nwrong use of weight_ica1\n");
    if  (cmd->anz_from2 != 1)
        fprintf (stderr, "\nwrong use of weight_ica1\n");
    if  (cmd->n_from1[0] != area)
        fprintf (stderr, "\nwrong areas in weight_ica1\n");

    // W  = cmd->W_target[inarea];
    // dW = cmd->dW_target[inarea];

    // for (i = 0; i < d_r_local; ++i)
    // for (j = 0; j < d_r_input; ++j)
    //     dW[i][j] += cmd->moment * W[i][j];

    for (int ct_n = 0; ct_n < A[area].d_n; ++ct_n) {
        conn = W->all_conn[cmd->area][ct_n][inarea];
        if  (conn->next)
        do  {
            conn = conn->next;
            conn->val += cmd->moment * conn->val * eps;
        } while (conn->next != NULL);
    }

    return (DOUBLE)(0);
}


/******************************* weight_list_ica2 ****************************/
/* Formerly total_weight_ica2.                                               */
/* Additional parts of the covariant ICA learning rule.                      */
/* Anti-Hebbian like terms which depend on the prior. Choose by quantum[0][0]*/

DOUBLE weight_list_ica2 (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    static int firsttime = 1;
    static int d_r_local;
    static int d_r_input;
    static DOUBLE **C;
    static DOUBLE **Wtemp, **Wsave;
    static DOUBLE *Y;
    static DOUBLE *FY;
    static DOUBLE *FYbeta, *FYzeta;

    // double **W, **dW;
    ALL_WEIGHTS *W = (ALL_WEIGHTS *) cmd->pointers[0]->data;
    CONN *conn;

    int function_mean_01       = 0;
    int function_tanh          = 0;
    int function_assymm_01     = 0;
    int function_gauss         = 0;
    int function_Lee_Sejnowski = 0;

         /****function_assymm_01:
         gnuplot
         f(x) = (exp (x))/(1 + exp(x))**2  
         g(x) = (exp (beta*x))/(1 + exp(beta*x))
         beta = 5.0
         plot f(x), g(x), f(x) * g(x)
         ****/

    const int inarea   = cmd->n_from2[0];
    const int area     = cmd->area;
    const DOUBLE *Pre  = cmd->S_from2[0][ct_t];

    if  (cmd->quantum[0][0] == 0)
        function_mean_01 = 1;
    if  (cmd->quantum[0][0] == 1)
        function_tanh = 1;
    if  (cmd->quantum[0][0] == 2)
        function_assymm_01 = 1;
    if  (cmd->quantum[0][0] == 3)
        function_gauss = 1;
    if  (cmd->quantum[0][0] == 4)
        function_Lee_Sejnowski = 1;

    if  (firsttime) {
        d_r_local = A[area].d_n;
        d_r_input = A[inarea].d_n;
        C         = d_matrix (d_r_local, d_r_local);
        Wtemp     = d_matrix (d_r_local, d_r_input);
        Wsave     = d_matrix (d_r_local, d_r_input);
        Y         = d_vector (d_r_local);
        FY        = d_vector (d_r_local);
        FYbeta    = d_vector (d_r_local);
        FYzeta    = d_vector (d_r_local);
        firsttime = 0;
    }

    if  (d_r_local != A[area].d_n)
        fprintf (stderr, "\ncan't use weight_ica2 twice\n");
    if  (d_r_local > A[inarea].d_n)
        fprintf (stderr, "\nno overcompleteness for weight_ica2, please!\n");
    if  (cmd->anz_from1 != 1)
        fprintf (stderr, "\nwrong use of weight_ica2\n");
    if  (cmd->anz_from2 != 1)
        fprintf (stderr, "\nwrong use of weight_ica2\n");
    if  (cmd->n_from1[0] != area)
        fprintf (stderr, "\nwrong areas in weight_ica2\n");
    if  (cmd->anz_quantums != 3)
        fprintf (stderr, "\nchoose function, beta and eps for weight_ica2\n");
    if  (cmd->moment != 1.0)
        fprintf (stderr, "\nweight_ica2: cmd->moment should be 1.0 !\n");

    DOUBLE beta = cmd->quantum[1][0];
    DOUBLE eps  = cmd->quantum[2][0];

//    W  = cmd->W_target[inarea];
//   dW = cmd->dW_target[inarea];
    // for (i = 0; i < d_r_local; ++i) {
    //     Y[i] = 0.0;
    //      for (j = 0; j < d_r_input; ++j)
    //          Y[i] += W[i][j] * Pre[j];
    //     /*
    //     if  (z->ch_Theta)
    //         Y[i] -= A[area].Theta[i];
    //     */
    //     if  (function_mean_01)
    //         FY[i] = 1.0 / (1.0 + exp (-Y[i]));
    //     if  (function_tanh)
    //         FY[i] = tanh (Y[i]);
    //     if  (function_assymm_01) {
    //         FY[i] = 1.0 / (1.0 + exp (-Y[i]));
    //         FYbeta[i] = 1.0 / (1.0 + exp (- beta * Y[i]));
    //     }
    // }

    /**get Y[i]**/
    for (int ct_n = 0; ct_n < A[area].d_n; ++ct_n) {
        Y[ct_n] = 0.0;
        conn = W->all_conn[cmd->area][ct_n][inarea];
        if  (conn->next)
        do  {
            conn = conn->next;
            Y[ct_n] += conn->val * Pre[conn->n];
        } while (conn->next != NULL);

        /*
        if  (z->ch_Theta)
            Y[i] -= A[area].Theta[i];
        */
        if  (function_mean_01)
            FY[ct_n] = 1.0 / (1.0 + exp (-Y[ct_n]));
        if  (function_tanh)
            FY[ct_n] = tanh (Y[ct_n]);
        if  (function_assymm_01) {
            FY[ct_n] = 1.0 / (1.0 + exp (-Y[ct_n]));
            FYbeta[ct_n] = 1.0 / (1.0 + exp (- beta * Y[ct_n]));
        }
        if  (function_Lee_Sejnowski) {
            FY[ct_n]     = tanh (Y[ct_n]);
            FYbeta[ct_n] = tanh (Y[ct_n] + beta);
            FYzeta[ct_n] = tanh (Y[ct_n] - beta);
        }
    }

    for (int i = 0; i < d_r_local; ++i)
    for (int k = 0; k < d_r_local; ++k) {
        if  (function_mean_01)
            C[i][k] = (1.0 - 2.0 * FY[i]) * Y[k];
        if  (function_tanh)
            C[i][k] = (    - 2.0 * FY[i]) * Y[k];
        if  (function_assymm_01)
            C[i][k] = (1.0 - 2.0 * FY[i] + beta - beta * FYbeta[i]) * Y[k];
        if  (function_gauss)
            C[i][k] = (    - Y[i] / (beta*beta)) * Y[k];
        if  (function_Lee_Sejnowski)
            C[i][k] = (      2.0*FY[i] - 2.0*FYbeta[i] - 2.0*FYzeta[i]) * Y[k];
        /**.../retina/SpatioChromaticRFs_Original-InfoMax_LeeSejnowski99.pdf**/
        /**page 10. beta=b=0 then Bell&Sejnowski; b=2 then subgaussian**/
    }

    /****this term is now NOT separate in total_weight_ica1****/
    for (int i = 0; i < d_r_local; ++i)
        C[i][i] += 1.0;
    

    /**fill weight matrix**/
    for (int ct_n = 0; ct_n < d_r_local; ++ct_n) {
        CONN *conn = W->all_conn[area][ct_n][inarea];
        if  (conn != NULL)
        if  (conn->next)
        do  {
            conn = conn->next;
            Wsave[ct_n][conn->n] = conn->val;
        } while (conn->next != NULL);
    }

    for (int i = 0; i < d_r_local; ++i)
    for (int j = 0; j < d_r_input; ++j) {
        Wtemp[i][j] = 0.0;
        for (int k = 0; k < d_r_local; ++k)
            Wtemp[i][j] += C[i][k] * Wsave[k][j];

            // Wtemp[i][j] += C[i][k] * W[k][j];
    }

    // for (i = 0; i < d_r_local; ++i)
    // for (j = 0; j < d_r_input; ++j)
    //     dW[i][j] += cmd->moment * Wtemp[i][j];

    for (int ct_n = 0; ct_n < A[area].d_n; ++ct_n) {
        conn = W->all_conn[cmd->area][ct_n][inarea];
        if  (conn->next)
        do  {
            conn = conn->next;
            conn->val += cmd->moment * Wtemp[ct_n][conn->n] * eps;
        } while (conn->next != NULL);
    }

   /*
   if  (z->ch_Theta)
       for (i = 0; i < d_r_local; ++i) {
           if  (function_mean_01)
               A[area].Theta_delta[i] -= cmd->moment * (1.0 - 2.0 * FY[i]);
                                    ** * A[area].Theta[i] * A[area].Theta[i] **
           if  (function_tanh)
               A[area].Theta_delta[i] -= cmd->moment * (    - 2.0 * FY[i]);
           if  (function_assymm_01)
               fprintf (stderr, "\nimplement Theta update for function_assymm!\n");
       }
   */

    for (int i = 0; i < d_r_local; ++i)
        cmd->S_target[ct_t][i] = FY[i];

    return (DOUBLE)(0);
}



/**************************** weight_list_decay ******************************/
/* Decays weights by -m * conn->val * q[0][0]. Returns weight vector length. */

DOUBLE weight_list_decay (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    ALL_WEIGHTS *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;
    DOUBLE quad_length = (DOUBLE)0;

    /**all input areas in command list**/
    for (int ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        const int inarea = cmd->n_from1[ct_l];                           /**area number in global list**/

        CONN *conn = W->all_conn[cmd->area][ct_n][inarea];               /**[area][ct_n][glob inarea](inarea admin neuron)**/

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


/**************************** weight_list_decay_const ************************/
/* Decays weights by -m * q[0][0].                                           */
/* (Used for Foldiak rule with q[0][0]<0 for constant growth.)               */

DOUBLE weight_list_decay_const (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    ALL_WEIGHTS *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;
    DOUBLE quad_length = (DOUBLE)0;

    /**all input areas in command list**/
    for (int ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        const int inarea = cmd->n_from1[ct_l];                           /**area number in global list**/

        CONN *conn = W->all_conn[cmd->area][ct_n][inarea];               /**[area][ct_n][glob inarea](inarea admin neuron)**/

        /**sum over weighted inputs**/
        if  (conn->next)
        do  {
            conn = conn->next;
            if  (conn->val > 0)
                conn->val -= cmd->moment * cmd->quantum[0][0];
            if  (conn->val < 0)
                conn->val += cmd->moment * cmd->quantum[0][0];
            quad_length += conn->val * conn->val;

        } while (conn->next != NULL);
    }

    return quad_length;
}



/**************************** weight_list_decay_quad *************************/
/* Decays weights by -m * conn->val * q[0][0] * quad_length.                 */

DOUBLE weight_list_decay_quad (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    ALL_WEIGHTS *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;
    DOUBLE quad_length = (DOUBLE)0;

    /**all input areas in command list**/
    for (int ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        const int inarea = cmd->n_from1[ct_l];                           /**area number in global list**/

        CONN *conn = W->all_conn[cmd->area][ct_n][inarea];               /**[area][ct_n][glob inarea](inarea admin neuron)**/

        /**sum over weighted inputs**/
        if  (conn->next)
        do  {
            conn = conn->next;
            quad_length += conn->val * conn->val;

        } while (conn->next != NULL);
    }

    /**all input areas in command list**/
    for (int ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        const int inarea = cmd->n_from1[ct_l];                           /**area number in global list**/

        CONN *conn = W->all_conn[cmd->area][ct_n][inarea];               /**[area][ct_n][glob inarea](inarea admin neuron)**/

        /**sum over weighted inputs**/
        if  (conn->next)
        do  {
            conn = conn->next;
            conn->val -= cmd->moment * conn->val * cmd->quantum[0][0] * quad_length;

        } while (conn->next != NULL);
    }

    return quad_length;
}



/**************************** weight_list_decay_post *************************/
/* Activity-dependent weight decay.                                          */
/* from1 must be here area to obtain post act,                               */
/* from2 are inareas used for selecting weights.                             */
/* q[1][0]=1: then funny decay (see below)                                   */

DOUBLE weight_list_decay_post (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    ALL_WEIGHTS *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;
/*
    if  (cmd->anz_quantums != 2)
        fprintf (stderr, "\nweight_list_decay_post wants 2 quantums!\n");
*/
    if  (cmd->anz_from2 != 1)
        fprintf (stderr, "\nweight_list_decay_post wants two from arguments!\n\n");
    if  (cmd->n_from1[0] != cmd->area)
        fprintf (stderr, "\nweight_list_decay_post: inconsistent target area!\n\n");

    const DOUBLE post = cmd->S_from1[0][ct_t][ct_n];                     /**[ct_l="here"][ct_t][ct_n]**/

    if  (post != 0.0)
    /**all input areas in command list**/
    for (int ct_l = 0; ct_l < cmd->anz_from2; ++ct_l) {

        const int inarea = cmd->n_from2[ct_l];                           /**area number in global list**/

        CONN *conn = W->all_conn[cmd->area][ct_n][inarea];               /**[area][ct_n][glob inarea](inarea admin neuron)**/

        /**sum over weighted inputs**/
        if  (conn->next)
        do  {
            conn = conn->next;

/*          if  (cmd->quantum[1][0] == 1.0) {
                if  (conn->val != 0.0)
                    conn->val -= cmd->moment * conn->val / fabs(conn->val) * cmd->quantum[0][0] * fabs (post);
            }
            else
*/
            {
                conn->val -= cmd->moment * conn->val * cmd->quantum[0][0] * post;    /** !!! fabs is now gone !!! **/
            }

        } while (conn->next != NULL);
    }

    return post;
}



/**************************** weight_list_decay_pre **************************/
/* Activity-dependent weight decay. from1 at inarea (as in normal decay)!    */
/* q[1][0]=1: then funny decay (see below)                                   */

DOUBLE weight_list_decay_pre (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    ALL_WEIGHTS *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;

    if  (cmd->anz_quantums != 2)
        fprintf (stderr, "\nweight_list_decay_pre wants 2 quantums!\n");

    /**all input areas in command list**/
    for (int ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        const int inarea = cmd->n_from1[ct_l];                           /**area number in global list**/

        CONN *conn = W->all_conn[cmd->area][ct_n][inarea];               /**[area][ct_n][glob inarea](inarea admin neuron)**/

        DOUBLE *pre = cmd->S_from1[ct_l][ct_t];                          /**[cmd inarea list][ct_t](inarea neurons)**/

        /**sum over weighted inputs**/
        if  (conn->next)
        do  {
            conn = conn->next;

            if  (cmd->quantum[1][0] == 1.0) {
                if  (conn->val != 0.0)
                    conn->val -= cmd->moment * conn->val / fabs(conn->val) * cmd->quantum[0][0] * fabs (pre[conn->n]);
            } else {
                conn->val -= cmd->moment * conn->val * cmd->quantum[0][0] * fabs (pre[conn->n]);    /** !!! fabs is new !!! **/
            }

        } while (conn->next != NULL);
    }

    return (DOUBLE)0;
}


/**************************** weight_list_decay_Baddeley *********************/
/* Like weight_list_decay_const, but acts only if weight length exceeds q10. */
/* Decays weights by -m * q[0][0].                                           */
/* q[1][0] threshold that has to be exceeded befor contraint is applied.     */

DOUBLE weight_list_decay_Baddeley (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    ALL_WEIGHTS *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;
    DOUBLE length = 0.0, length_pos = 0.0, length_neg = 0.0;

    static int monitor = 0;
    if  (ct_n == 0)
        monitor = 0;

    /**all input areas in command list**/
    for (int ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        const int inarea = cmd->n_from1[ct_l];                           /**area number in global list**/

        /**first, get weight vector length**/
        CONN *conn = W->all_conn[cmd->area][ct_n][inarea];               /**[area][ct_n][glob inarea](inarea admin neuron)**/

        if  (conn->next)
        do  {
            conn = conn->next;
            length += fabs(conn->val);
            if  (conn->val > 0.0)
                length_pos += conn->val;
            if  (conn->val < 0.0)
                length_neg += conn->val;

        } while (conn->next != NULL);

        /**then, apply constraint if vector too large**/
        conn = W->all_conn[cmd->area][ct_n][inarea];               /**[area][ct_n][glob inarea](inarea admin neuron)**/

        if  (cmd->anz_quant[0] == 1) {
            if  (length > cmd->quantum[1][0])
            if  (conn->next)
            do  {
                conn = conn->next;
                DOUBLE val_old = conn->val;
                if  (conn->val > 0)
                    conn->val -= cmd->moment * cmd->quantum[0][0];
                if  (conn->val < 0)
                    conn->val += cmd->moment * cmd->quantum[0][0];

                if  (conn->val * val_old < 0.0)  /**prevent sign-change!**/
                    conn->val = 0.0;

                monitor += 1;

            } while (conn->next != NULL);

        } else {

            if  ((cmd->anz_quant[0] != 2) || (cmd->anz_quant[1] != 2))
               fprintf (stderr, "\nweight_list_decay_Baddeley: check params!");

            if  (conn->next)
            do  {
                conn = conn->next;
                DOUBLE val_old = conn->val;
                if  (length_pos > cmd->quantum[1][0])
                    if  (conn->val > 0)
                        conn->val -= cmd->moment * cmd->quantum[0][0];
                if  (length_neg < -cmd->quantum[1][1])
                    if  (conn->val < 0)
                        conn->val += cmd->moment * cmd->quantum[0][1];

                if  (conn->val * val_old < 0.0)  /**prevent sign-change!**/
                    conn->val = 0.0;

            } while (conn->next != NULL);
        }
    }

    static int monitor_ct = 0;
    monitor_ct ++;

    if  (monitor_ct == 400) {
        if  (ct_n == A[cmd->area].d_n - 1)
            fprintf (stderr, " constraint applied to %d ", monitor);
        monitor_ct = 0;
    }

    return length;
}


/*********************** weight_list_decay_Baddeley_distance *****************/
/* Like weight_list_decay_const, but acts only if weight length exceeds q10. */
/* Decays weights by -m * q[0][0]                                            */
/* q[1][0] threshold that has to be exceeded before contraint is applied     */
/* q[2][0] param added to distance, so that at dist 0 there is punishment    */
/* q[3][0] exists then topographic                                           */

DOUBLE weight_list_decay_Baddeley_distance (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    ALL_WEIGHTS   *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;
    int           inarea = cmd->n_from1[0];               /**only one inarea**/
    static DOUBLE *vec;
    static int    firsttime = 1;

    if  (cmd->anz_quantums < 4)
        fprintf (stderr,"\nw_l_decay_Baddeley_distance: not enough params!\n");

    if  (firsttime) {
        vec   = d_vector (A[inarea].d_a * A[inarea].d_b);
        firsttime = 0;
    }

        static int monitor = 0;
        if  (ct_n == 0)
            monitor = 0;

        /**first, assign vec, and get weight vector length**/
        CONN *conn = W->all_conn[cmd->area][ct_n][inarea];

        if  (conn->next)
        do  {
            conn = conn->next;
            vec[conn->n] = conn->val;

        } while (conn->next != NULL);

        /**CM of RF**/
        float cm_a = 0.0;
        float cm_b = 0.0;
        if  (cmd->anz_quantums == 4) {             /**set CM topographically**/
            float area_a = (float)(ct_n / A[cmd->area].d_a);
            float area_b = (float)(ct_n % A[cmd->area].d_b);
            cm_a = area_a / (float)(A[cmd->area].d_a) * (float)(A[inarea].d_a);
            cm_b = area_b / (float)(A[cmd->area].d_b) * (float)(A[inarea].d_b);

        } else {                                /**take max val for CM of RF**/
            float max  = 0.0;
            for (int a = 0; a < A[inarea].d_a; ++a)
            for (int b = 0; b < A[inarea].d_b; ++b) {
                if  (fabs(vec[a * A[inarea].d_b + b]) > max) {
                    cm_a = (float)a;
                    cm_b = (float)b;
                    max  = fabs(vec[a * A[inarea].d_b + b]);
                }
            }
        }

        DOUBLE length = 0.0;
        for (int a = 0; a < A[inarea].d_a; ++a)
        for (int b = 0; b < A[inarea].d_b; ++b) {
            float dist = sqrt ( (float)(a - cm_a) * (a - cm_a)
                              + (float)(b - cm_b) * (b - cm_b));
            length += fabs (vec[a*A[inarea].d_b+b])
                      * (cmd->quantum[2][0] + dist);
        }

        /**then, apply constraint if vector too large**/
        conn = W->all_conn[cmd->area][ct_n][inarea];

        if  (length > cmd->quantum[1][0])
            if  (conn->next)
            do  {
                conn = conn->next;
                DOUBLE val_old = conn->val;

                float dist = sqrt ( (float)(conn->a - cm_a) * (conn->a - cm_a)
                                  + (float)(conn->b - cm_b) * (conn->b -cm_b));

                if  (conn->val > 0)
                    conn->val -= cmd->moment * cmd->quantum[0][0]
                               * (cmd->quantum[2][0] + dist);
                if  (conn->val < 0)
                    conn->val += cmd->moment * cmd->quantum[0][0]
                               * (cmd->quantum[2][0] + dist);

                if  (conn->val * val_old < 0.0)  /**prevent sign-change!**/
                    conn->val = 0.0;

                monitor += 1;

            } while (conn->next != NULL);

    static int monitor_ct = 0;
    monitor_ct ++;

    if  (monitor_ct == 800) {
        if  (ct_n == A[cmd->area].d_n - 1) {
            fprintf (stderr, " constraint applied to %d ", monitor);
            fprintf (stderr, " cm_a=%.2f cm_b=%.2f ", cm_a, cm_b);
        }
        monitor_ct = 0;
    }

    return 0.0;
}


/**************************** weight_list_normalize **************************/
/* Normalises weight vector. (Each input area separately!)                   */
/* q[0][0]: length after normalisation.                                      */
/* q[1][0]=1: use absolute values; =2: use root of sum of squares.           */

DOUBLE weight_list_normalize (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    ALL_WEIGHTS *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;

    /**all input areas in command list**/
    for (int ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        const int inarea = cmd->n_from1[ct_l];                           /**area number in global list**/
        DOUBLE length = 0.0;

        /**get length of weight vector**/
        CONN *conn = W->all_conn[cmd->area][ct_n][inarea];               /**[area][ct_n][glob inarea](inarea admin neuron)**/

        if  (conn->next)
        do  {
            conn = conn->next;

            if  (cmd->quantum[1][0] == 1)
                length += fabs (conn->val);

            if  (cmd->quantum[1][0] == 2)
                length += conn->val * conn->val;

        } while (conn->next != NULL);

        if  (cmd->quantum[1][0] == 2)
            length = sqrt (length);


        /**normalise**/
        conn = W->all_conn[cmd->area][ct_n][inarea];                     /**[area][ct_n][glob inarea](inarea admin neuron)**/

        if  (conn->next)
        do  {
            conn = conn->next;
            conn->val /= length;
            conn->val *= cmd->quantum[0][0];

        } while (conn->next != NULL);
    }

    return (DOUBLE)0;
}



/**************************** weight_list_rectify ****************************/

DOUBLE weight_list_rectify (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    ALL_WEIGHTS *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;

    /**all input areas in command list**/
    for (int ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

        const int inarea = cmd->n_from1[ct_l];                       /**area number in global list**/

        CONN *conn = W->all_conn[cmd->area][ct_n][inarea];           /**[area][ct_n][glob inarea](first inarea neuron)**/

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



/**************************** weight_list_cutself ****************************/

DOUBLE weight_list_cutself (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int num_cut = 0;

    /**test whether inner area connections exist at all**/
    int OK = 0;
    for (int ct_l = 0; ct_l < cmd->anz_from1; ++ct_l)
        if  (cmd->area == cmd->n_from1[ct_l])
            OK = 1;
    if  (!OK)
        fprintf (stderr, "\n\nweight_list_noself applied to non-self-connections!\n");


    ALL_WEIGHTS *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;

    CONN *conn = W->all_conn[cmd->area][ct_n][cmd->area];
    int more = 1;
    while (more) {

        CONN *prev = conn;
        conn = conn->next;
        CONN *tofree = NULL;

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
        fprintf (stderr, "\n\nweight_list_cutself has cut %d connection(s) at neuron %d instead of one!\n\n", num_cut, ct_n);

    return (DOUBLE)(0);
}



/**************************** weight_list_cutsmall ***************************/
/* q[0][0] = 1: regard cutting threshold as absolute; = 0 in % of min/max.   */
/* q[1][0] = negative cutting threshold in % of minimum value                */
/* q[1][1] = positive cutting threshold in % of maximum value                */

DOUBLE weight_list_cutsmall (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int num_cut = 0;
    DOUBLE min = 99;
    DOUBLE max = -99;

    if  ((cmd->anz_quantums != 2) || (cmd->anz_quant[1] != 2)) {
        fprintf (stderr, "\ninsufficient arguments in weight_list_cutsmall\n");
        exit (1);
    }

    static int average_cut  = 0;                    /**to print out for info**/
    static int total_before = 0;                    /**to print out for info**/
    static int total_after  = 0;
    if  (ct_n == 0) {
        fprintf (stderr, "\nweight_list_cutsmall ");
        average_cut  = 0;
        total_before = 0;
        total_after  = 0;
    }

    ALL_WEIGHTS *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;

    const int inarea = cmd->n_from1[0];       /**use only 1st inarea in list**/

    CONN *conn = W->all_conn[cmd->area][ct_n][inarea];

    total_before += count_connections (W->all_conn[cmd->area][ct_n][inarea]);

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

    conn = W->all_conn[cmd->area][ct_n][inarea];

    int last;
    if  (conn->next)
    do  {
        CONN *prev = conn;
        conn = conn->next;

        last = (conn->next == NULL) ? 1 : 0;

        if  ((conn->val > thres_neg) && (conn->val < thres_pos)) {

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
    total_after += count_connections (W->all_conn[cmd->area][ct_n][inarea]);

    if  (ct_n == A[cmd->area].d_n - 1) {
        fprintf (stderr, " cut %.1f connections in average, or %d total at obs_%s_%d_%d  ",
                          (double)average_cut/(double)A[cmd->area].d_n, average_cut, cmd->ch_pointers[0], cmd->area, inarea);
        fprintf (stderr, "\ntotal connections before - after = cut: %d - %d = %d  ", total_before, total_after, total_before - total_after);
    }

    return (DOUBLE)(0);
}



/**************************** weight_list_sprout *****************************/
/* q[0][0] = connection value in % of minimum value around which to sprout   */
/* q[0][1] = connection value in % of maximum value around which to sprout   */
/* q[1][0] = value of new connections in % of minimum value                  */
/* q[1][1] = value of new connections in % of maximum value                  */
/* q[2][0] = 4/8 for sprouting at N/S/W/E or also at diagonal directions     */
/* q[3][0] = 3 for RGB-inarea, so sprouting does not cross the boundaries    */
/* q[3][1] = 2 for separate ON/OFF-inarea,    "                              */

DOUBLE weight_list_sprout (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    if  (ct_n == 0)
        fprintf (stderr, "\nweight_list_sprout ");

    DOUBLE max = -99;
    DOUBLE min = 99;

    ALL_WEIGHTS *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;

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

    CONN *conn = W->all_conn[cmd->area][ct_n][inarea];

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

            CONN *here = conn;
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


void rgbmatrix2file (COMMAND *cmd, DOUBLE ***rgbmatrix, int d_a, int d_b, int d_a_in, int d_b_in, DOUBLE min, DOUBLE max, int format, int key) {

    char fullname[512];
    char ending[4];
    if  (format == 9)
        sprintf (ending, "ext");
    else
        sprintf (ending, "pnm");

    if  (key == 1) {
        /**for activations**/
        sprintf (fullname, "%s/obs_%c_%d.%s", cmd->pointers[0]->words[0], cmd->ch_target, cmd->area, ending);                             /**first pointer gives the directory!**/

    } else {
        /**for weights (key = 2)**/
        if  (cmd->anz_pointers != 2)
            fprintf (stderr, "\n\n\nrgbmatrix2file needs two pointers in cmd, because weights are exported!!\n\n\n");

        sprintf (fullname, "%s/obs_%s_%d_%d.%s", cmd->pointers[1]->words[0], cmd->ch_pointers[0], cmd->area, cmd->n_from1[0], ending);    /**second pointer gives the directory!**/
    }

    std::fstream fp(fullname, std::ios::out | std::ios::binary);  // Creates an ofstream object named fp

    if  (! fp) {            // Always test file open
        fprintf (stderr, "\nError opening output file %s\n", fullname);
        exit (1);
    }

    const int hcol_int = 666666;                  /**for format=3; int export: max values are MAXCOL_INT**/
    const /*unsigned*/ char hcol_char = 66;       /**for format=6; char export: values scaled between [0:255]**/

    /**write header 1st line**/
    if  (format == 3)
        fp << "P3" << std::endl;
    else
       if  (format == 6)
           fp << "P6" << std::endl;
       else
           if  (format == 9)
               fp << "P9" << std::endl;
           else
               fprintf (stderr, "\nexport format must be 3 or 6 or 9\n");

    /**write header comments: max, min**/
    fp << "# highS: " << max << "  lowS: " << min << std::endl;

    int h = 2;

    /**write header sizes**/
    fp << h + (d_b_in + h) * d_b << " " << h + (d_a_in + h) * d_a << std::endl;

    /**write grey values**/
    if  (format == 3)
        fp << MAXCOL_INT << std::endl;
    if  (format == 6)
        fp << MAXCOL_CHAR << std::endl;
    if  (format == 9)
        fp << 0 << std::endl; /**just because a line will be read in anyway**/

    for (int i = 0; i < h; ++i)                          /**1st white line**/
        for (int j = 0; j < h + (d_b_in + h) * d_b; ++j) {
            if  (format == 3)
                fp << hcol_int << " " << hcol_int << " " << hcol_int << " ";
            if  (format == 6) {
                // fp << hcol_char << hcol_char << hcol_char;
                fp.write ( &hcol_char, sizeof (/*unsigned*/ char));
                fp.write ( &hcol_char, sizeof (/*unsigned*/ char));
                fp.write ( &hcol_char, sizeof (/*unsigned*/ char));
            }
        }

    for (int i = 0; i < d_a; ++i) {
        for (int a = 0; a < d_a_in; ++a) {

            for (int k = 0; k < h; ++k) {
                if  (format == 3)
                    fp << hcol_int << " " << hcol_int << " " << hcol_int << " ";
                if  (format == 6) {
                    // fp << hcol_char << hcol_char << hcol_char;
                    fp.write ( &hcol_char, sizeof (/*unsigned*/ char));
                    fp.write ( &hcol_char, sizeof (/*unsigned*/ char));
                    fp.write ( &hcol_char, sizeof (/*unsigned*/ char));
                }
            }

            for (int j = 0; j < d_b; ++j) {
                for (int b = 0; b < d_b_in; ++b) {

                        if  (format == 3) {
                            int red   = (int)(rgbmatrix[0][i * d_b + j][a * d_b_in + b]);
                            int green = (int)(rgbmatrix[1][i * d_b + j][a * d_b_in + b]);
                            int blue  = (int)(rgbmatrix[2][i * d_b + j][a * d_b_in + b]);
                            fp << red << " " << green << " " << blue << " ";
                        }
                        if  (format == 6) {
                            const /*unsigned*/ char red   = (/*unsigned*/ char)(rgbmatrix[0][i * d_b + j][a * d_b_in + b]);
                            const /*unsigned*/ char green = (/*unsigned*/ char)(rgbmatrix[1][i * d_b + j][a * d_b_in + b]);
                            const /*unsigned*/ char blue  = (/*unsigned*/ char)(rgbmatrix[2][i * d_b + j][a * d_b_in + b]);
                            // fp << red << green << blue;
                            fp.write ( &red, sizeof (/*unsigned*/ char));
                            fp.write ( &green, sizeof (/*unsigned*/ char));
                            fp.write ( &blue, sizeof (/*unsigned*/ char));
                        }
                        if  (format == 9) {
                            DOUBLE value = rgbmatrix[0][i * d_b + j][a * d_b_in + b];
                            // fwrite (&value, sizeof (DOUBLE), 1 , fp);
                            fp.write ( (char *)&value, sizeof (DOUBLE));

                            if  (((value < min) || (value > max)) && (value != 0.0))
                                fprintf (stderr, "\nvalue out of bounds: %f while min=%f, max=%f\n", value, min, max);
                        }
                }

                for (int k = 0; k < h; ++k) {
                    if  (format == 3)
                        fp << hcol_int << " " << hcol_int << " " << hcol_int << " ";
                    if  (format == 6) {
                        // fp << hcol_char << hcol_char << hcol_char;
                        fp.write ( &hcol_char, sizeof (/*unsigned*/ char));
                        fp.write ( &hcol_char, sizeof (/*unsigned*/ char));
                        fp.write ( &hcol_char, sizeof (/*unsigned*/ char));
                    }
                }
	    }
        }

        for (int k = 0; k < h; ++k)         /**between + last white lines**/
            for (int j = 0; j < h + (d_b_in + h) * d_b; ++j) {
                if  (format == 3)
                    fp << hcol_int << " " << hcol_int << " " << hcol_int << " ";
                if  (format == 6) {
                    // fp << hcol_char << hcol_char << hcol_char;
                    fp.write ( &hcol_char, sizeof (/*unsigned*/ char));
                    fp.write ( &hcol_char, sizeof (/*unsigned*/ char));
                    fp.write ( &hcol_char, sizeof (/*unsigned*/ char));
                }
            }
    }

    fp.close();
}


void file2rgbmatrix (COMMAND *cmd, DOUBLE ***rgbmatrix, int d_a, int d_b, int d_a_in, int d_b_in, DOUBLE *min, DOUBLE *max, int format, int weight_or_act) {

    int width, height;
    int h = 2;  /**frame**/

    std::cerr << " file2rgbmatrix ";

    char fullname[512];
    char ending[4];
    if  (format == 9)
        sprintf (ending, "ext");
    else
        sprintf (ending, "pnm");

    if  (weight_or_act == 1)
        sprintf (fullname, "%s/obs_%s_%d_%d.%s", cmd->pointers[1]->words[0], cmd->ch_pointers[0], cmd->area, cmd->n_from1[0], ending);
    else
        sprintf (fullname, "%s/obs_%c_%d.%s", cmd->pointers[0]->words[0], cmd->ch_target, cmd->area, ending);

    std::cerr << " reads: " << fullname;

    std::fstream fp(fullname, std::ios::in | std::ios::binary);
    if  (! fp) {
        fprintf (stderr, "\nError opening input file %s\n", fullname);
        exit (1);
    }

    int hcol_int;                     /**for format=3; int export: max values are MAXCOL_INT**/
    /*unsigned*/ char hcol_char;      /**for format=6; char export: values scaled between [0:255]**/

    std::string kommentarzeile;
    const char *c_kommentarzeile;

    /**the P6**/
    getline (fp, kommentarzeile);
    std::cerr << std::endl << kommentarzeile.c_str() << std::endl;


    /**max/min**/
    getline (fp, kommentarzeile);
    c_kommentarzeile = kommentarzeile.c_str();

    double d_highS, d_lowS; /**must be double because use of %lf in next line**/
    sscanf  (c_kommentarzeile, "# highS: %lf lowS: %lf", &d_highS, &d_lowS);
    std::cerr << "highS=" << d_highS << ", lowS=" << d_lowS << std::endl;

    *min = (DOUBLE)d_lowS;
    *max = (DOUBLE)d_highS;


    /**dimensions**/
    getline (fp, kommentarzeile);
    c_kommentarzeile = kommentarzeile.c_str();

    sscanf (c_kommentarzeile, "%d %d\n", &width, &height);
    std::cerr << "width=" << width << " height=" << height << std::endl;

    /**check size consistency**/
    if  ((width != h + (d_b_in + h) * d_b) || (height != h + (d_a_in + h) * d_a)) {
        std::cerr << "inconsistent image size\n";
        std::cerr << "expected width " << h + (d_b_in + h) * d_b << ", and height " << h + (d_a_in + h) * d_a << std::endl;
        exit (0);
    }


    /**the 255**/
    getline (fp, kommentarzeile);
    std::cerr << kommentarzeile.c_str() << std::endl;


    for (int i = 0; i < h; ++i)                          /**1st white line**/
        for (int j = 0; j < h + (d_b_in + h) * d_b; ++j) {
            if  (format == 3)
                fp >> hcol_int >> hcol_int >> hcol_int;
            if  (format == 6) {
                // fp >> hcol_char >> hcol_char >> hcol_char;
                fp.read ( &hcol_char, sizeof (char));
                fp.read ( &hcol_char, sizeof (char));
                fp.read ( &hcol_char, sizeof (char));
            }
        }

    for (int i = 0; i < d_a; ++i) {
        for (int a = 0; a < d_a_in; ++a) {

            for (int k = 0; k < h; ++k) {
                if  (format == 3)
                    fp >> hcol_int >> hcol_int >> hcol_int;
                if  (format == 6) {
                    // fp >> hcol_char >> hcol_char >> hcol_char;
                    fp.read ( &hcol_char, sizeof (char));
                    fp.read ( &hcol_char, sizeof (char));
                    fp.read ( &hcol_char, sizeof (char));
                }
            }

            for (int j = 0; j < d_b; ++j) {
                for (int b = 0; b < d_b_in; ++b) {

                        if  (format == 3) {
                            int red, green, blue;
                            fp >> red >> green >> blue;
                            rgbmatrix[0][i * d_b + j][a * d_b_in + b] = (DOUBLE)red;
                            rgbmatrix[1][i * d_b + j][a * d_b_in + b] = (DOUBLE)green;
                            rgbmatrix[2][i * d_b + j][a * d_b_in + b] = (DOUBLE)blue;
                        }
                        if  (format == 6) {
                            /*unsigned*/ char red, green, blue;
                            // fp >> red >> green >> blue;
                            fp.read ( &red, sizeof (char));
                            fp.read ( &green, sizeof (char));
                            fp.read ( &blue, sizeof (char));
                            rgbmatrix[0][i * d_b + j][a * d_b_in + b] = (DOUBLE)red;
                            rgbmatrix[1][i * d_b + j][a * d_b_in + b] = (DOUBLE)green;
                            rgbmatrix[2][i * d_b + j][a * d_b_in + b] = (DOUBLE)blue;
                        }
                        if  (format == 9) {
                            DOUBLE value;
                            // fread (&value, sizeof (DOUBLE), 1, fp);
                            fp.read ( (char *)&value, sizeof (DOUBLE));
                            rgbmatrix[0][i * d_b + j][a * d_b_in + b] = value;
                        }
                }

                for (int k = 0; k < h; ++k) {
                    if  (format == 3)
                        fp >> hcol_int >> hcol_int >> hcol_int;
                    if  (format == 6) {
                        // fp >> hcol_char >> hcol_char >> hcol_char;
                        fp.read ( &hcol_char, sizeof (char));
                        fp.read ( &hcol_char, sizeof (char));
                        fp.read ( &hcol_char, sizeof (char));
                    }
                }
	    }
        }

        for (int k = 0; k < h; ++k)         /**between + last white lines**/
            for (int j = 0; j < h + (d_b_in + h) * d_b; ++j) {
                if  (format == 3)
                    fp >> hcol_int >> hcol_int >> hcol_int;
                if  (format == 6) {
                    // fp >> hcol_char >> hcol_char >> hcol_char;
                    fp.read ( &hcol_char, sizeof (char));
                    fp.read ( &hcol_char, sizeof (char));
                    fp.read ( &hcol_char, sizeof (char));
                }
            }
    }


    /*unsigned*/ char mycharacter;
    if (!fp.eof()) {
        std::cerr << "one dangling character in file: -->";
        fp >> mycharacter;
        std::cerr << mycharacter << "<--   ";
    }
    if (!fp.eof()) {
        std::cerr << "\n2nd dangling character in file: -->";
        fp >> mycharacter;
        std::cerr << mycharacter << "<--  this is bad!   " << std::endl;
    }
    if (!fp.eof()) {
        std::cerr << "\n3rd dangling character in file: -->";
        fp >> mycharacter;
        std::cerr << mycharacter << "<--  now this is really bad!  " << std::endl;
    }

    fp.close();

fprintf (stderr, " end of file2rgbmatrix ");
}


void weights2rgbmatrix (ALL_WEIGHTS *W, int area, int inarea, DOUBLE ***rgbmatrix, int d_a, int d_b, int d_a_in, int d_b_in, DOUBLE *min, DOUBLE *max, int format) {

    if  ((format != 3) && (format != 9))
        format = 6;

    DOUBLE maxcol;

    if  (format == 3)
        maxcol = MAXCOL_INT;
    else /**format == 6**/
        maxcol = MAXCOL_CHAR;

    /**background init with a light green tone**/
    for (int i = 0; i < d_a * d_b; ++i)
        for (int j = 0; j < d_a_in * d_b_in; ++j) {
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

        CONN *conn = W->all_conn[area][ct_n][inarea];

        if  (conn != NULL)
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

        CONN *conn = W->all_conn[area][ct_n][inarea];

        if  (conn != NULL)
        if  (conn->next)
        do  {
            conn = conn->next;

            /**positive (blue remains bright)**/
            if  (conn->val >= 0) {
                rgbmatrix[0][ct_n][conn->n] = maxcol * (1.0 - conn->val / largest);
                rgbmatrix[1][ct_n][conn->n] = maxcol * (1.0 - conn->val / largest);
                rgbmatrix[2][ct_n][conn->n] = maxcol;
            }

            /**negative (red remains bright)**/
            if  (conn->val < 0) {
                rgbmatrix[0][ct_n][conn->n] = maxcol;
                rgbmatrix[1][ct_n][conn->n] = maxcol * (1.0 + conn->val / largest);
                rgbmatrix[2][ct_n][conn->n] = maxcol * (1.0 + conn->val / largest);
            }

            /**in format 9 only write val into first component only**/
            if  (format == 9)
                rgbmatrix[0][ct_n][conn->n] = conn->val;

        } while (conn->next != NULL);
    }
}



/************************ weight_list_export *********************************/
/* Exports the weights to target area from only the first input area.        */
/* q[0][0]=3/9: export in integer / extended(.ext) format (else char).       */

DOUBLE weight_list_export (PARAMS *g, AREA *A, COMMAND *cmd, int, int) {

    fprintf (stderr, "  weight_list_export ... ");

    DOUBLE min, max;
    int format = ((int)(cmd->quantum[0][0]) == 0) ? 6 : (int)(cmd->quantum[0][0]);

    ALL_WEIGHTS *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;

    const int area   = cmd->area;
    const int inarea = cmd->n_from1[0];       /**use only 1st inarea in list**/

    int d_a    = A[area].d_a;
    int d_b    = A[area].d_b;
    int d_a_in = A[inarea].d_a;
    int d_b_in = A[inarea].d_b;

    DOUBLE ***rgbmatrix = d_tensor (3, d_a * d_b, d_a_in * d_b_in);

    weights2rgbmatrix (W, area, inarea, rgbmatrix, d_a, d_b, d_a_in, d_b_in, &min, &max, format);

    rgbmatrix2file (cmd, rgbmatrix, d_a, d_b, d_a_in, d_b_in, min, max, format, 2);
                                                                              /*weights*/
    free_d_tensor (rgbmatrix, 3, d_a * d_b);

    fprintf (stderr, " ... end ");

    return (DOUBLE)(0);
}


/************************ weight_list_export_col *****************************/
/* Exports weights assuming that the input area is RGB retina.               */
/* q[0][0] if > 0 scales the background greenish, but doesn't look so nice.  */

DOUBLE weight_list_export_col (PARAMS *g, AREA *A, COMMAND *cmd, int, int) {

    fprintf (stderr, "  weight_list_export_col  ");

    ALL_WEIGHTS *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;

    const int area   = cmd->area;
    const int inarea = cmd->n_from1[0];       /**use only 1st inarea in list**/

    const int d_a     = A[area].d_a;
    const int d_b     = A[area].d_b;
    const int d_a_eff = A[inarea].d_a / 3;
    const int d_b_in  = A[inarea].d_b;

    DOUBLE ***rgbmatrix = d_tensor (3, d_a * d_b, d_a_eff * d_b_in);
    int **occupymatrix  = i_matrix (d_a * d_b, d_a_eff * d_b_in);

    int format = 6;        /**put into cmd->quantum if adjustability desired!**/
    DOUBLE maxcol;

    if  (format == 3)
        maxcol = MAXCOL_INT;
    else /**format == 6**/
        maxcol = MAXCOL_CHAR;

    /**find maximum and minimum weight strength**/
    DOUBLE max = -99;
    DOUBLE min = 99;

    for (int ct_n = 0; ct_n < d_a * d_b; ++ct_n) {

        CONN *conn = W->all_conn[area][ct_n][inarea];

        if  (conn->next)
        do  {
            conn = conn->next;
            max = conn->val > max ? conn->val : max;
            min = conn->val < min ? conn->val : min;

        } while (conn->next != NULL);
    }


    /**set rgbmatrix values where weights exist; positive&negative are normalised together w.r.t. to largest of max/min**/
    DOUBLE largest = max > -min ? max : -min;

    /**only negative values disturb, but if all are positive, then subtraction distorts the true colors**/
    double shiftS = min < 0.0 ? min : 0.0;

    /**background init with zero**/
    for (int i = 0; i < d_a * d_b; ++i)
        for (int j = 0; j < d_a_eff * d_b_in; ++j) {
            rgbmatrix[0][i][j] = maxcol * (0.0 - shiftS) / (largest - shiftS);    /**red**/
            rgbmatrix[1][i][j] = maxcol * (0.0 - shiftS) / (largest - shiftS);    /**green**/
            rgbmatrix[2][i][j] = maxcol * (0.0 - shiftS) / (largest - shiftS);    /**blue**/

            occupymatrix[i][j] = 0;
        }

    for (int ct_n = 0; ct_n < d_a * d_b; ++ct_n) {

        CONN *conn = W->all_conn[area][ct_n][inarea];

        if  (conn->next)
        do  {
            conn = conn->next;

            /**red third of retina**/
            if  (conn->n < d_a_eff * d_b_in) {
                rgbmatrix[0][ct_n][conn->n]                        = maxcol * (conn->val - shiftS) / (largest - shiftS);
                occupymatrix[ct_n][conn->n]                        = 1;
            }
            /**green third of retina**/
            if  ((conn->n >= d_a_eff * d_b_in) && (conn->n < 2 * d_a_eff * d_b_in)) {
                rgbmatrix[1][ct_n][conn->n - d_a_eff * d_b_in]     = maxcol * (conn->val - shiftS) / (largest - shiftS);
                occupymatrix[ct_n][conn->n - d_a_eff * d_b_in]     = 1;
            }
            /**blue third of retina**/
            if  (conn->n >= 2 * d_a_eff * d_b_in) {
                rgbmatrix[2][ct_n][conn->n - 2 * d_a_eff * d_b_in] = maxcol * (conn->val - shiftS) / (largest - shiftS);
                occupymatrix[ct_n][conn->n - 2 * d_a_eff * d_b_in] = 1;
            }

        } while (conn->next != NULL);
    }

    /**background with no connections a bit greenish**/
    for (int i = 0; i < d_a * d_b; ++i)
        for (int j = 0; j < d_a_eff * d_b_in; ++j)
            if  (occupymatrix[i][j] == 0) {
                rgbmatrix[0][i][j] = maxcol * (0.0                          - shiftS) / (largest - shiftS);              /**red**/
                rgbmatrix[1][i][j] = maxcol * (cmd->quantum[0][0] * largest - shiftS) / (largest - shiftS);              /**green**/
                rgbmatrix[2][i][j] = maxcol * (0.0                          - shiftS) / (largest - shiftS);              /**blue**/
        }

    rgbmatrix2file (cmd, rgbmatrix, d_a, d_b, d_a_eff, d_b_in, min, max, format, 2);
                                                                               /*weights*/

    free_d_tensor (rgbmatrix, 3, d_a * d_b);
    free_i_matrix (occupymatrix, d_a * d_b);

    return (DOUBLE)(0);
}


/************************ weight_list_export_gnu *****************************/
/* Exports weights assuming that the input area has just 2 or 3 dims.        */

DOUBLE weight_list_export_gnu (PARAMS *g, AREA *A, COMMAND *cmd, int, int) {

    fprintf (stderr, "  weight_list_export_gnu  ");

    ALL_WEIGHTS *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;

    const int area   = cmd->area;
    const int inarea = cmd->n_from1[0];       /**use only 1st inarea in list**/
    const int d_a     = A[area].d_a;
    const int d_b     = A[area].d_b;
    const int d_in    = A[inarea].d_a * A[inarea].d_b;

    char fullname[512];
    sprintf (fullname, "%s/obs_%s_%d_%d.dat", cmd->pointers[1]->words[0], cmd->ch_pointers[0], cmd->area, cmd->n_from1[0]);

    FILE *fp = fopen (fullname, "w");

    for (int i = 0; i < d_a; ++i) {
        for (int j = 0; j < d_b; ++j) {

            CONN *conn = W->all_conn[area][i * d_b + j][inarea];

            if  (d_in == 2)
                fprintf (fp, "%f %f\n", conn->next->val, conn->next->next->val);
            if  (d_in == 3)
                fprintf (fp, "%f %f %f\n", conn->next->val, conn->next->next->val, conn->next->next->next->val);
        }
        fprintf (fp, "\n");
    }

    fclose (fp);
    return (DOUBLE)0;
}


/************************ weight_list_alloc_import ***************************/
/* Imports the weights to target area from only the first input area.        */
/* q[0][0]=3/9: import in integer / extended(.ext) format (else char).       */
/* q[1][0] existent, then assign connections also to zero-weights.           */

DOUBLE weight_list_alloc_import (PARAMS *g, AREA *A, COMMAND *cmd, int, int) {

    /**find maximum and minimum weight strength**/
    DOUBLE min, max;

    fprintf (stderr, "\nweight_list_alloc_import  ");

    weight_list_alloc_init_once (g, A, cmd, 0, 0);

    ALL_WEIGHTS *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;

    const int area   = cmd->area;
    const int inarea = cmd->n_from1[0];       /**use only 1st inarea in list**/

    int d_a    = A[area].d_a;
    int d_b    = A[area].d_b;
    int d_a_in = A[inarea].d_a;
    int d_b_in = A[inarea].d_b;

    DOUBLE ***rgbmatrix = d_tensor (3, d_a * d_b, d_a_in * d_b_in);

    int format = (int)(cmd->quantum[0][0]);

    if  ((format != 3) && (format != 9))
        format = 6;

    DOUBLE maxcol;

    if  (format == 3)
        maxcol = MAXCOL_INT;
    else /**format == 6**/
        maxcol = MAXCOL_CHAR;

    file2rgbmatrix (cmd, rgbmatrix, d_a, d_b, d_a_in, d_b_in, &min, &max, format, 1);

    DOUBLE largest = max > -min ? max : -min;

    for (int ct_n = 0; ct_n < d_a * d_b; ++ct_n) {

        W->all_conn[area][ct_n][inarea] = (CONN *) malloc (sizeof (CONN));

        CONN *conn = W->all_conn[area][ct_n][inarea];

        /**create the admin element**/
        conn->val  = 0;
        conn->a    = -1;
        conn->b    = -1;
        conn->n    = -1;
        conn->next = NULL;

        int j = -1;
        for (int ct_a_in = 0; ct_a_in < d_a_in; ++ct_a_in)
        for (int ct_b_in = 0; ct_b_in < d_b_in; ++ct_b_in) {
            j += 1;

            DOUBLE value = 0.0;
            int OK = 0;

            /**positive weight**/
            if  (rgbmatrix[2][ct_n][j] == maxcol) {
                value = largest * (1.0 - rgbmatrix[1][ct_n][j] / maxcol);
                OK = 1;
            }

            /**negative weight**/
            if  (rgbmatrix[0][ct_n][j] == maxcol) {
                value = largest * (rgbmatrix[1][ct_n][j] / maxcol - 1.0);
                OK = 1;
            }

            /**zero or no weight**/
            if  (rgbmatrix[1][ct_n][j] == maxcol) {
                value = 0.0;
                OK = 1;
            }

            /**ignore upper if this format is used**/
            if  (format == 9) {
                value = rgbmatrix[0][ct_n][j];
                OK = 1;
            }

/*            if  (((value < min - 0.000001) || (value > max + 0.000001)) && (value != 0.0))
                fprintf (stderr, "\nweight_list_alloc_import: value is out of bounds: %f while min=%f, max=%f  ", value, min, max);
*/
            if  (!OK) {
                fprintf (stderr, "\nNO VALUE ASSIGNED! ct_n=%d a_in=%d b_in=%d r,g,b = %.1f, %.1f, %.1f", ct_n, ct_a_in, ct_b_in, rgbmatrix[0][ct_n][j], rgbmatrix[1][ct_n][j], rgbmatrix[2][ct_n][j]);
            }

            int add_new_element = (value != 0.0) ? 1 : 0;
            if  (cmd->anz_quantums > 1)
                add_new_element = 1;

            if  (add_new_element) {

                conn->next = (CONN *) malloc (sizeof (CONN));

                conn       = conn->next;

                conn->val  = value;
                conn->a    = ct_a_in;
                conn->b    = ct_b_in;
                conn->n    = j;
                conn->next = NULL;
            }
        }
    }

fprintf (stderr, "\nweight_list_import  2 ");

    free_d_tensor (rgbmatrix, 3, d_a * d_b);

    fprintf (stderr, "weight_list_import end ");

    return (DOUBLE)(0);
}


/************************ weight_list_alloc_import_tiles *********************/
/* For huge map. Areas are 2x or 3x larger than imported (toroidal) weights. */
/* q[0][0]=3/9: import in integer / extended(.ext) format (else char).       */
/* q[1][0] threshold to assign weights ( <0 to assign also to zero-weights). */
/* q[2][0/1] number of tiles in x/y-direction (input as well as output layer)*/

DOUBLE weight_list_alloc_import_tiles (PARAMS *g, AREA *A, COMMAND *cmd, int, int) {

    /**find maximum and minimum weight strength**/
    DOUBLE min, max;

    fprintf (stderr, "\nweight_list_alloc_import_tiles  ");

    weight_list_alloc_init_once (g, A, cmd, 0, 0);

    ALL_WEIGHTS *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;

    const int area   = cmd->area;
    const int inarea = cmd->n_from1[0];       /**use only 1st inarea in list**/

    int format              = (int)(cmd->quantum[0][0]);
    double weight_threshold =       cmd->quantum[1][0];
    int num_tiles_a         = (int)(cmd->quantum[2][0]);
    int num_tiles_b         = (int)(cmd->quantum[2][1]);

    int d_a       = A[area].d_a;
    int d_b       = A[area].d_b;
    int d_a_in    = A[inarea].d_a;
    int d_b_in    = A[inarea].d_b;

    int small_d_a    = d_a / num_tiles_a;
    int small_d_b    = d_b / num_tiles_b;
    int small_d_n    = small_d_a * small_d_b;
    int small_d_a_in = d_a_in / num_tiles_a;
    int small_d_b_in = d_b_in / num_tiles_b;
    int small_d_n_in = small_d_a_in * small_d_b_in;

    DOUBLE ***rgbmatrix = d_tensor (3, small_d_n, small_d_n_in);

    if  ((format != 3) && (format != 9))
        format = 6;

    DOUBLE maxcol;

    if  (format == 3)
        maxcol = MAXCOL_INT;
    else /**format == 6**/
        maxcol = MAXCOL_CHAR;

    file2rgbmatrix (cmd, rgbmatrix, d_a / num_tiles_a, d_b / num_tiles_b, d_a_in / num_tiles_a, d_b_in / num_tiles_b, &min, &max, format, 1);

    DOUBLE largest = max > -min ? max : -min;


    for (int ct_a = 0; ct_a < d_a; ++ct_a)
    for (int ct_b = 0; ct_b < d_b; ++ct_b) {

        int ct_n = ct_a * d_b + ct_b;
        int small_ct_a = ct_a % small_d_a;
        int small_ct_b = ct_b % small_d_b;
        int small_ct_n = small_ct_a * small_d_b + small_ct_b;

        W->all_conn[area][ct_n][inarea] = (CONN *) malloc (sizeof (CONN));

        CONN *conn = W->all_conn[area][ct_n][inarea];

        /**create the admin element**/
        conn->val  = 0;
        conn->a    = -1;
        conn->b    = -1;
        conn->n    = -1;
        conn->next = NULL;

        int j = -1;
        for (int ct_a_in = 0; ct_a_in < d_a_in; ++ct_a_in)
        for (int ct_b_in = 0; ct_b_in < d_b_in; ++ct_b_in) {

          DOUBLE value = 0.0;
          int add_new_element = 0;

          if  ((ct_a_in < small_d_a_in) && (ct_b_in < small_d_b_in)) {

            j += 1;
            int OK = 0;

            /**positive weight**/
            if  (rgbmatrix[2][small_ct_n][j] == maxcol) {
                value = largest * (1.0 - rgbmatrix[1][small_ct_n][j] / maxcol);
                OK = 1;
            }

            /**negative weight**/
            if  (rgbmatrix[0][small_ct_n][j] == maxcol) {
                value = largest * (rgbmatrix[1][small_ct_n][j] / maxcol - 1.0);
                OK = 1;
            }

            /**zero or no weight**/
            if  (rgbmatrix[1][small_ct_n][j] == maxcol) {
                value = 0.0;
                OK = 1;
            }

            /**ignore upper if this format is used**/
            if  (format == 9) {
                value = rgbmatrix[0][small_ct_n][j];
                OK = 1;
            }

            if  (((value < min - 0.000001) || (value > max + 0.000001)) && (value != 0.0))
                fprintf (stderr, "\nweight_list_alloc_import_tiles: value is out of bounds: %f while min=%f, max=%f  ", value, min, max);

            if  (!OK) {
                fprintf (stderr, "\nNO VALUE ASSIGNED! ct_n=%d a_in=%d b_in=%d r,g,b = %.1f, %.1f, %.1f",
                                                       ct_n, ct_a_in, ct_b_in,
                                                       rgbmatrix[0][ct_n][j], rgbmatrix[1][ct_n][j], rgbmatrix[2][ct_n][j]);
            }

            if  (fabs (value) >= weight_threshold)
                add_new_element = 1;
          }

            if  (add_new_element) {

                conn->next = (CONN *) malloc (sizeof (CONN));

                conn       = conn->next;

                conn->val  = value;
                conn->a    = ct_a_in;
                conn->b    = ct_b_in;
                conn->n    = ct_a_in * d_b_in + ct_b_in  /*j*/;
                conn->next = NULL;
            }
        }
    }

    fprintf (stderr, "\nweight_list_alloc_import_tiles  imported ");

    free_d_tensor (rgbmatrix, 3, small_d_n);


    /**re-arrange the RF's**/
    for (int ct_a = 0; ct_a < d_a; ++ct_a)
    for (int ct_b = 0; ct_b < d_b; ++ct_b) {

        /**find maximum abs weight strength**/
        int ct_n = ct_a * d_b + ct_b;
        DOUBLE max = -99;
        int a_max = -1, b_max = -1;

        CONN *conn = W->all_conn[area][ct_n][inarea];

        if  (conn != NULL)
        if  (conn->next)
        do  {
            conn = conn->next;
            if  (fabs (conn->val) > max) {

                max   = fabs (conn->val);
                a_max = conn->a;
                b_max = conn->b;
            } 
        } while (conn->next != NULL);

        /**place the max into proper tile**/
        int tile_a = ct_a / small_d_a;
        int tile_b = ct_b / small_d_b;

        /**shift assumes that max is originally in tile (0,0)**/
        int a_max_new = a_max + tile_a * small_d_a_in;
        int b_max_new = b_max + tile_b * small_d_b_in;

        /**go through all connections and adjust their a and b positions**/
        conn = W->all_conn[area][ct_n][inarea];

        if  (conn != NULL)
        if  (conn->next)
        do  {
            CONN *prev = conn;
            conn = conn->next;
            int last = (conn->next == NULL) ? 1 : 0;

            if  (conn->a - a_max_new > small_d_a_in / 2) {
                conn->a -= small_d_a_in;
            } else {
                while (abs (conn->a - a_max_new) > small_d_a_in / 2) {
                    conn->a += small_d_a_in;
                }
            }
            if  (conn->b - b_max_new > small_d_b_in / 2) {
                conn->b -= small_d_b_in;
            } else {
                while (abs (conn->b - b_max_new) > small_d_b_in / 2) {
                    conn->b += small_d_b_in;
                }
            }
            conn->n = conn->a * d_b_in + conn->b;

            /**cut ectopic connections**/
            if  ((conn->a >= d_a_in) || (conn->a < 0) || (conn->b >= d_b_in) || (conn->b < 0)) {
                if  (last)
                    prev->next = NULL;
                else
                    prev->next = conn->next;  
                free (conn);
                conn = prev;
            }

        } while (conn->next != NULL);
    }


    fprintf (stderr, "weight_list_alloc_import_tiles end ");

    return (DOUBLE)(0);
}



/************************ weight_list_histogram ******************************/

DOUBLE weight_list_histogram (PARAMS *g, AREA *A, COMMAND *cmd, int, int) {

    fprintf (stderr, "\nweight_list_histogram:");

    ALL_WEIGHTS *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;

    const int area   = cmd->area;
    const int inarea = cmd->n_from1[0];       /**use only 1st inarea in list**/

    /**find maximum and minimum weight strength**/
    DOUBLE max = -99;
    DOUBLE min = 99;

    int count_all = 0;
    int number_without = 0;

    for (int ct_n = 0; ct_n < A[area].d_n; ++ct_n) {

        CONN *conn = W->all_conn[area][ct_n][inarea];
        int count_here = 0;

        if  (conn != NULL)
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
        CONN *conn = W->all_conn[area][ct_n][inarea];
        if  (conn != NULL)
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
    DOUBLE lower = min;
    DOUBLE upper = min;
    for (int i = 0; i < bins; ++i) {
        upper += interval;
        fprintf (stderr, "\n[%f - %f]: %d ", lower, upper, bincount[i]);
        lower += interval;
    }



    /**test weights for consistent coverage and a,b,n consistent?**/

    int **connected = i_matrix(W->d_a[inarea], W->d_b[inarea]);

    for (int ct_n = 0; ct_n < A[area].d_n; ++ct_n) {

        for (int i = 0; i < W->d_a[inarea]; ++i)
        for (int j = 0; j < W->d_b[inarea]; ++j)
            connected[i][j] = 0;

        CONN *conn = W->all_conn[area][ct_n][inarea];

        if  (conn != NULL)
        if  (conn->next)
        do  {
            conn = conn->next;

            connected[conn->a][conn->b] += 1;

            if  (conn->a * A[inarea].d_b + conn->b != conn->n)
                fprintf (stderr, "\n\n\nInconsistent d_a, d_b with d_n!\n\n\n");

        } while (conn->next != NULL);

        for (int i = 0; i < W->d_a[inarea]; ++i)
        for (int j = 0; j < W->d_b[inarea]; ++j)
            if  (connected[i][j] != 0)
                if  (connected[i][j] != 1)
                    fprintf (stderr, "\n\nInconsistent connected matrix with connected[%d][%d] = %d!\n\n", i, j, connected[i][j]);;
    }

    free_i_matrix (connected, W->d_a[inarea]);


    fprintf (stderr, " %s_%d_%d has", cmd->ch_pointers[0], area, inarea);
    fprintf (stderr, " in average %.1f connections of %d possible ", (double)count_all/(double)A[area].d_n, A[inarea].d_n);
    fprintf (stderr, " and %d units have no connection  ", number_without);

    return (DOUBLE)(0);
}




float pyramid (float x, float lower, float upper) {

   if  (lower < upper) {

       const float width = upper - lower;

       if  (x <= lower)
           return 0.0;
       if  (x >= upper)
           return 0.0;
       if  (x < lower + 0.25 * width)
           return (x-lower) / (0.25 * width);
       if  (x > lower + 0.75 * width)
           return ((upper-x) / (0.25 * width));

       // if  ((x >= lower + 0.25 * width) && (x <= lower + 0.75 * width))
       return 1.0;

   } else {

       const float width = 2.0 * (lower - upper);

       if  ((x >= upper) && (x <= lower))
           return 0.0;
       if  ((x > lower) && (x < lower + 0.25 * width))
           return (x-lower) / (0.25 * width);
       if  ((x < upper) && (x > upper - 0.25 * width))
           return (upper-x) / (0.25 * width);

       return 1.0;
   }
}

void give_color (DOUBLE *red_out, DOUBLE *green_out, DOUBLE *blue_out, DOUBLE color, DOUBLE color_range, DOUBLE bright, DOUBLE bright_range, int format) {

        float red   = pyramid ((float)color, 3.0/6.0*(float)color_range, 1.0/6.0*(float)color_range);
        float green = pyramid ((float)color, 5.0/6.0*(float)color_range, 3.0/6.0*(float)color_range);
        float blue  = pyramid ((float)color, 1.0/6.0*(float)color_range, 5.0/6.0*(float)color_range);

        red   = 1.0 - red   * bright / bright_range;
        green = 1.0 - green * bright / bright_range;
        blue  = 1.0 - blue  * bright / bright_range;

        DOUBLE maxcol;

        if  (format == 3)
            maxcol = MAXCOL_INT;
        else /**format == 6**/
            maxcol = MAXCOL_CHAR;

        *red_out   = (DOUBLE)(red * maxcol);
        *green_out = (DOUBLE)(green * maxcol);
        *blue_out  = (DOUBLE)(blue * maxcol);
}

/************************ observe_phase **************************************/
/* q[0][0]=3/9: export in integer/extended format (else char).               */

DOUBLE observe_phase (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {

    const int area   = cmd->area;

    int d_a    = A[area].d_a;
    int d_b    = A[area].d_b;

    DOUBLE ***rgbmatrix = d_tensor (3, (end - begin), d_a * d_b);

    int format = (int)(cmd->quantum[0][0]);

    if  ((format != 3) && (format != 9))
        format = 6;

    /**find maximum and minimum activation (and phase)**/
    DOUBLE max = -99;
    DOUBLE min = 99;
    DOUBLE phase_max = -99;
    DOUBLE phase_min = 99;
    DOUBLE used_phase_max = 1.0;   /**phases are between 0..1, i.e. NOT between 0..2PI**/
    for (int ct_t = begin; ct_t < end; ++ct_t)
        for (int ct_n = 0; ct_n < d_a * d_b; ++ct_n) {

            max = cmd->S_from1[0][ct_t][ct_n] > max ? cmd->S_from1[0][ct_t][ct_n] : max;
            min = cmd->S_from1[0][ct_t][ct_n] < min ? cmd->S_from1[0][ct_t][ct_n] : min;

            if  (cmd->anz_from2 == 1) {
                phase_max = cmd->S_from2[0][ct_t][ct_n] > phase_max ? cmd->S_from2[0][ct_t][ct_n] : phase_max;
                phase_min = cmd->S_from2[0][ct_t][ct_n] < phase_min ? cmd->S_from2[0][ct_t][ct_n] : phase_min;

                /* used_phase_max = cmd->quantum[0][0]; */ /* if used build in somehow again! ... */
            }
        }

    DOUBLE highest = max > -min ? max : -min;

    if  (phase_max > used_phase_max)
        fprintf (stderr, "\nobserve_phase warning: phase values larger than %f exist!\n", used_phase_max);
    if  (phase_min < 0.0)
        fprintf (stderr, "\nobserve_phase warning: phase values smaller than 0 exist!\n");

    for (int ct_t = begin; ct_t < end; ++ct_t)
        for (int ct_n = 0; ct_n < d_a * d_b; ++ct_n) {

            DOUBLE value = cmd->S_from1[0][ct_t][ct_n];
            DOUBLE phase;

            if  (cmd->anz_from2 == 1) {
                phase = cmd->S_from2[0][ct_t][ct_n];
            } else {
                if  (cmd->S_from1[0][ct_t][ct_n] >= 0.0) {   /** positive activations **/
                    phase = 0.0;                             /**will be displayed blue**/
                } else {                                     /** negative activations **/
                    phase = used_phase_max / 3.0;            /**will be displayed red **/
                    value = -cmd->S_from1[0][ct_t][ct_n];
                }
            }

            give_color (rgbmatrix[0][ct_t-begin]+ct_n, rgbmatrix[1][ct_t-begin]+ct_n, rgbmatrix[2][ct_t-begin]+ct_n,
                        phase, used_phase_max, value, highest, format);

            if  (format == 9)
                rgbmatrix[0][ct_t-begin][ct_n] = cmd->S_from1[0][ct_t][ct_n];   /**ignore rest and overwrite**/
        }

    rgbmatrix2file (cmd, rgbmatrix, 1, (end-begin), d_a, d_b, min, max, format, 1);
                                                                              /*act*/
    free_d_tensor (rgbmatrix, 3, (end-begin));

    return (DOUBLE)(0);
}


/*************************** import_phase ************************************/
/* To import activations from file. Make sure that rlen matches.             */
/* q[0][0]=3/9: import in integer/extended format (else char).               */

DOUBLE import_phase (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {

    fprintf (stderr, "\nimport_phase ");

    const int area   = cmd->area;
    DOUBLE min, max;

    int d_a    = A[area].d_a;
    int d_b    = A[area].d_b;

    DOUBLE ***rgbmatrix = d_tensor (3, (end - begin), d_a * d_b);

    int format = (int)(cmd->quantum[0][0]);

    if  ((format != 3) && (format != 9))
        format = 6;

    DOUBLE maxcol;

    if  (format == 3)
        maxcol = MAXCOL_INT;
    else /**format == 6**/
        maxcol = MAXCOL_CHAR;

    file2rgbmatrix (cmd, rgbmatrix, 1, (end-begin), d_a, d_b, &min, &max, format, 2);

    DOUBLE largest = max > -min ? max : -min;

    for (int ct_n = 0; ct_n < (end - begin); ++ct_n) {     /**was:  int ct_n = 1; (??)**/

        int j = -1;
        for (int ct_a_in = 0; ct_a_in < d_a; ++ct_a_in)
        for (int ct_b_in = 0; ct_b_in < d_b; ++ct_b_in) {
            j += 1;

            DOUBLE value = 0.0;
            int OK = 0;

            /**positive**/
            if  (rgbmatrix[2][ct_n][j] == maxcol) {
                value = largest * (1.0 - rgbmatrix[1][ct_n][j] / maxcol);
                OK = 1;
            }

            /**negative**/
            if  (rgbmatrix[0][ct_n][j] == maxcol) {
                value = largest * (rgbmatrix[1][ct_n][j] / maxcol - 1.0);
                OK = 1;
            }

            /**green pixels should only be if "no weight", but not here, except if everything is max (act=0)!**/
            if  ((rgbmatrix[1][ct_n][j] == maxcol) && (rgbmatrix[0][ct_n][j] != maxcol)) {
                value = 0.0;
                OK = 0;
            }

            /**ignore upper if this format is used**/
            if  (format == 9) {
                value = rgbmatrix[0][ct_n][j];
                OK = 1;
            }

            if  ((value < min) || (value > max))
                fprintf (stderr, "\nVALUE EXCEEDS BOUNDS: %f while min=%f, max=%f\n", value, min, max);

            if  (!OK) {
               fprintf(stderr,"\nNO VALUE ASSIGNED in import_phase! ct_n(=t)=%d a_in=%d b_in=%d r,g,b = %.1f, %.1f, %.1f",
                             ct_n, ct_a_in, ct_b_in, rgbmatrix[0][ct_n][j], rgbmatrix[1][ct_n][j], rgbmatrix[2][ct_n][j]);
            }

            cmd->S_target[ct_n][ct_a_in * d_b + ct_b_in] = value;
        }
    }

    free_d_tensor (rgbmatrix, 3, (end-begin));

    return (DOUBLE)(0);
}



int get_y_from_n (int n, int d_b) {

    return (n % d_b);
}


/**************************** weight_list_cuthalfinput ***********************/
/* q[0][0] = 0/1 for cut all input weights of left/right half of input area. */

DOUBLE weight_list_cuthalfinput (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    int last;

    ALL_WEIGHTS *W = (ALL_WEIGHTS *)cmd->pointers[0]->data;

    const int inarea = cmd->n_from1[0];       /**use only 1st inarea in list**/

    /**remove doomed connections**/
    CONN *conn = W->all_conn[cmd->area][ct_n][inarea];

    if  (conn->next)
    do  {
        CONN *prev = conn;
        conn = conn->next;

        last = (conn->next == NULL) ? 1 : 0;

      /*int x_pos = get_x_from_n (conn->n, A[inarea].d_a, A[inarea].d_b, 1);*/
        int y_pos = get_y_from_n (conn->n, A[inarea].d_b);            /**3 for RGB retina**/

        int cut_this_conn = 0;

        if  (cmd->quantum[0][0] == 0)
            if  (y_pos < A[inarea].d_b / 2)
                cut_this_conn = 1;

        if  (cmd->quantum[0][0] == 1)
            if  (y_pos >= A[inarea].d_b / 2)
                cut_this_conn = 1;

        if  (cut_this_conn == 1) {

            if  (last)
                prev->next = NULL;
            else
                prev->next = conn->next;  

            free (conn);
            conn = prev;
        }
    } while (! last);

    return (DOUBLE)(0);
}
