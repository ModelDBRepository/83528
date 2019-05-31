#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <strings.h>
#include <string.h>
#include <math.h>  /**for finite**/
/* #include <ieeefp.h> *for finite, on some computers*/

#include "coco.h"
#include "series.h"
#include "vehicle.h"      /**for print_command (debugging)**/
#include "relax.h" /*not needed here but so that compiler checks for consistency*/


/****************************** do_command ***********************************/
/* Calls appropriate cmd->function.                                          */

void do_command (PARAMS *x, AREA *A, COMMAND *cmd, int ct_t, int sel) {

#if  REPORT_ERR
fprintf (stderr, "%s  in area %d, neuron %d ", cmd->func_name, cmd->area, sel);
print_command (cmd);
#endif

    switch (cmd->func_key) {

        case 'l':
                  cmd->S_target[ct_t][sel] = (*cmd->localfunc) (cmd->quantum[0], cmd->S_from1[0][ct_t][sel], cmd->S_from2[0][ct_t][sel]);
                                                                                         /**  ^ the only "source" **/
                  if  (! finite (cmd->S_target[ct_t][sel])) {
                      std::cerr << std::endl << "local function " << cmd->func_name << " returns infinite number to " << cmd->ch_target;
                      std::cerr << " from " << cmd->ch_from1[0] << cmd->ch_from2[0] << " which now is " << cmd->S_from1[0][ct_t][sel] << std::endl;
                      exit (0);
                  }
                  break;

        case 's':
                  cmd->S_target[ct_t][sel] = (*cmd->func) (x, A, cmd, ct_t, sel);

                  if  (! finite (cmd->S_target[ct_t][sel])) {
                      fprintf (stderr, "\nfeed function %s returns infinite value.\n", cmd->func_name);
                      fprintf (stderr, "first input area=%d, in_activity=%c; target area=%d, target activity=%c; relaxation time=%d, neuron=%d\n",
                                            cmd->n_from1[0], cmd->ch_from1[0],     cmd->area,      cmd->ch_target,               ct_t,       sel);
                      exit (0);
                  }
                  break;

        case 'a': fprintf (stderr, "\nWARNING: attempt to invoke alltime_function: %s at neuron level!\n", cmd->func_name);
                  break;

        case 't': fprintf (stderr, "\nWARNING: attempt to invoke total_function: %s at neuron level (Argument \"(t)\" forgotten?)!\n", cmd->func_name);
                  break;

        default:  fprintf (stderr, "\ndo_command: no valid func_key \"%c\" at function: %s\n", cmd->func_key, cmd->func_name);
                  exit (0);
                  break;
    }

ERR("command_end");
}



/****************************** do_stay **************************************/
/* Counts neurons. Counts commands. Calls do_command.                        */
/* First, st->area is handed down to cmd->area.                              */

void do_stay (PARAMS *x, AREA *A, STAY *st, int ct_t) {

    int ct_n = 0, sel = 0, ct_c = 0;
    const char st_update = st->st_update;
    const int  area      = st->area;
    const int  d_n       = A[area].d_n;
    int        *shuffle  = A[area].shuffle;
    const int  anz_cmd   = st->anz_cmd;

ERR("\ndo_stay");

    /**hand down area to the commands (better do this in init!)**/
    for (ct_c = 0; ct_c < anz_cmd; ++ct_c)

        st->cmd[ct_c].area = area;


    /**only for total_functions**/
    if  (st_update == 't') {

        /**commands**/
        for (ct_c = 0; ct_c < anz_cmd; ++ct_c) {

            COMMAND *cmd = st->cmd+ct_c;

            (*cmd->func) (x, A, cmd, ct_t, 0);
        }

    } else {

        if  (st_update == 's') /*for shuffle (from http://remus.rutgers.edu/~rhoads/Code/random.shuffle)*/
            for (ct_n = 0; ct_n < d_n; ++ct_n)
                 shuffle[ct_n] = ct_n;

        /**neurons**/
        for (ct_n = 0; ct_n < d_n; ++ct_n) {

            if  (st_update == 'o')   /*order*/
                sel = ct_n;

            if  (st_update == 'r')   /*random*/
                sel = (int)(drand48 () * (double)(d_n));

            if  (st_update == 's') { /*shuffle*/
                int store     = shuffle[ct_n];
                int pick      = ct_n + (int)(drand48() * (d_n - ct_n));
                shuffle[ct_n] = shuffle[pick];
                shuffle[pick] = store;
                sel           = shuffle[ct_n];
            }

            /**commands**/
            for (ct_c = 0; ct_c < anz_cmd; ++ct_c)

                do_command (x, A, st->cmd+ct_c, ct_t, sel);
        }
    }

ERR("stay_end");
}


/****************************** do_sweep *************************************/
/* Counts relaxations. Counts stays.                                         */

void do_sweep (PARAMS *x, AREA *A, SWEEP *sw) {

    int ct_t, ct_st, sel_st, ct_cmd;
    COMMAND *cmd;
    STAY *st;

    int        area = 0;
    int        *shuffle  = A[area].shuffle;

ERR("\ndo_sweep");

    if  (!strcmp (sw->sw_update, "order"))

        /**relaxations**/
        for (ct_t = sw->begin; ct_t < sw->end; ++ct_t)

            /**stays**/
            for (ct_st = 0; ct_st < sw->anz_st; ++ct_st)

                /**condition**/
                if  ( (*sw->cond[ct_st]->test)(sw->cond[ct_st]->left, sw->cond[ct_st]->right) )

                    do_stay (x, A, sw->st+ct_st, ct_t);

ERR("sweep1");

    if  (!strcmp (sw->sw_update, "random"))

        /**relaxations**/
        for (ct_t = sw->begin; ct_t < sw->end; ++ct_t)

            /**stays**/
            for (ct_st = 0; ct_st < sw->anz_st; ++ct_st) {

                sel_st = (int)(drand48() * (double)(sw->anz_st));

                /**condition**/
                if  ( (*sw->cond[sel_st]->test)(sw->cond[sel_st]->left, sw->cond[sel_st]->right) )

                    do_stay (x, A, sw->st+sel_st, ct_t);
            }

ERR("sweep2");

    if  (!strcmp (sw->sw_update, "propto"))        /**invokes directly do_command without invoking do_stay -- uses 'shuffle'-update in order not to leave any neuron of all areas out**/

        /**relaxations**/
        for (ct_t = sw->begin; ct_t < sw->end; ++ct_t) {

            int sum_d_n, upper_sum_d_n, lower_sum_d_n,
                ct_st_n, ch_st = -1, ch_ar = -1, sel_n, ch_n = -1;

            /**copy activations from last time ** shouldn't be necessary, because of shuffled update **
            if  (ct_t)
                for (ct_st = 0; ct_st < sw->anz_st; ++ct_st) {
                    memcpy ((sw->st+ct_st)->cmd[0].S_target[ct_t],
                            (sw->st+ct_st)->cmd[0].S_target[ct_t - 1],
                            A[(sw->st+ct_st)->area].d_n * sizeof (double));
                }
            **/

            /**how many neurons at all?**/
            sum_d_n = 0;
            for (ct_st = 0; ct_st < sw->anz_st; ++ct_st)
                sum_d_n += A[(sw->st+ct_st)->area].d_n;

            /*init for shuffle (from http://remus.rutgers.edu/~rhoads/Code/random.shuffle)*/
            for (ct_st_n = 0; ct_st_n < sum_d_n; ++ct_st_n)
                shuffle[ct_st_n] = ct_st_n;

            /**stays and neurons**/
            for (ct_st_n = 0; ct_st_n < sum_d_n; ++ct_st_n) {

                /**select neuron from all of them
                sel_n = (int)(drand48 () * (double)(sum_d_n));
                **/

                { /**selection must be like "shuffle", because elsewise some neurons are not updated, and their activations cannot be retreived in the next iteration with time delay '-1'**/
                    int store        = shuffle[ct_st_n];
                    int pick         = ct_st_n + (int)(drand48() * (sum_d_n - ct_st_n));
                    shuffle[ct_st_n] = shuffle[pick];
                    shuffle[pick]    = store;
                    sel_n            = shuffle[ct_st_n];
                }

                /**get stay and area of this neuron, and its number**/
                lower_sum_d_n = 0;
                for (ct_st = 0; ct_st < sw->anz_st; ++ct_st) {

                    upper_sum_d_n = lower_sum_d_n
                                  + A[(sw->st+ct_st)->area].d_n;

                    /**is the neuron within the current stay/area?**/
                    if  ((sel_n >= lower_sum_d_n) && (sel_n < upper_sum_d_n)) {

                        /**actual stay**/
                        ch_st = ct_st;

                        /**actual area**/
                        ch_ar = (sw->st+ct_st)->area;

                        /**actual neuron within area**/
                        ch_n = sel_n - lower_sum_d_n;
                    }

                    lower_sum_d_n = upper_sum_d_n;
                }

                if  ((ch_st == -1) || (ch_ar == -1) || (ch_n == -1))
                    fprintf (stderr, "\n\nincorrect area or neuron\n\n");

                /**condition**/
                if  ( (*sw->cond[ch_st]->test)(sw->cond[ch_st]->left, sw->cond[ch_st]->right) )

                /**commands**/
                for (ct_cmd = 0; ct_cmd < (sw->st+ch_st)->anz_cmd; ++ct_cmd)

                    do_command (x, A, (sw->st+ch_st)->cmd+ct_cmd, ct_t, ch_n);
            }
        }

ERR("sweep3");

    if  (  (!strcmp (sw->sw_update, "alltime")   ) )

        for (ct_st = 0; ct_st < sw->anz_st; ++ct_st) {

            st = sw->st + ct_st;

            /**condition**/
            if  ( (*sw->cond[ct_st]->test)(sw->cond[ct_st]->left, sw->cond[ct_st]->right) )

            /**commands (outer loop; relaxations will be inside)**/
            for (ct_cmd = 0; ct_cmd < st->anz_cmd; ++ct_cmd) {

                cmd = st->cmd + ct_cmd;

                #if REPORT_ERR
                    fprintf (stderr, " do %s ", cmd->func_name);
                #endif

                (*cmd->func) (x, A, cmd, sw->begin, sw->end);
            }
        }

ERR("sweep4");

    int sel;
    if  (1 == sscanf (sw->sw_update, "%d", &sel))        /**if the sw_update word is a number, update that neuron (use local/single functions!)**/

        /**relaxations**/
        for (ct_t = sw->begin; ct_t < sw->end; ++ct_t)

            /**stays**/
            for (ct_st = 0; ct_st < sw->anz_st; ++ct_st) {

                st = sw->st + ct_st;

                /**commands**/
                for (int ct_c = 0; ct_c < st->anz_cmd; ++ct_c)

                    do_command (x, A, st->cmd+ct_c, ct_t, sel);
            }

ERR("sweep_end");
}




void do_series (PARAMS *x, AREA *A, SERIES *se) {

ERR("\ndo_series");

    x->iter = 1;        /** => header has no effect; necessary as a reset for following series**/
                        /**starts at 1, because "<=" in the "while" below, for convenient last iter at e.g. 1000**/

    fprintf (stderr, " at iter=%d for ilen=%d  ", x->iter, se->ilen);

    /**iterations**/
    while (x->iter <= se->ilen) {

        if  (x->iter % 50 == 0)
            fprintf (stderr, "\niter %d", x->iter);
        else
            fprintf (stderr, ".");

        /**sweeps**/
        for (int ct_sw = 0; ct_sw < se->anz_sw; ++ct_sw) {

            /**condition**/
            if  ( (*se->cond[ct_sw]->test)(se->cond[ct_sw]->left, se->cond[ct_sw]->right) )

                do_sweep (x, A, se->sw+ct_sw);
        }

        x->iter++;
    }

ERR("series_end");
}



void do_simulation (PARAMS *x, AREA *A, SIMULATION *si) {

    int ct_se;

ERR("\ndo_simulation");

    /**series**/
    for (ct_se = 0; ct_se < si->anz_se; ++ct_se) {

        fprintf (stderr, "\nconsidering series %d ", ct_se);

        /**condition**/
        if  ( (*si->cond[ct_se]->test)(si->cond[ct_se]->left, si->cond[ct_se]->right) )

             do_series (x, A, si->se+ct_se);
    }

ERR("simulation_end");
}
