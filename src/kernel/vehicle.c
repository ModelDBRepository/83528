#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <strings.h>
#include <string.h>

#include "coco.h"                   /**PARAMS, AGENT, AREA, STATE**/
#include "series.h"                 /**SERIES, SWEEP, STAY, COMMAND**/

#include "utils.h"                  /**for floatprint**/
#include "../fun_basic/local.h"     /**for cmd->local function pointers**/
#include "../fun_basic/single.h"    /**for cmd->feed function pointers**/
#include "../fun_basic/total.h"     /**for cmd->total function ptr's**/
#include "../fun_basic/observe.h"   /**for observe function pointers and importP36...**/
#include "../fun_special/reinforcement.h"
#include "../fun_weight_list/weight_list.h"
#include "../fun_weight_list/weight_sigpi.h"
#include "../fun_data/data.h"       /**for import function pointers**/

#include "vehicle.h"                /**for print_command at the bottom**/

#if  GAZEBO
#include "../fun_gazebo/gazebo.io.h"
#endif

#if  MIRO
#include "../fun_miro/miro.io.h"
#include "../fun_miro/mirror.io.h"
#endif

/**into file named one of:  agent.c aid.c assist.c transmit.c vehicle.c **/


/************************** assign_conditions ***************************/
/* While conditions are numbered according to occurrence in file in si, */
/* they need to be assigned to each series, sweep and stay.             */

void assign_conditions (SIMULATION *si) {

    int ct_all_cond = 0;

    fprintf (stderr, "\nassigning conditions:");

    /**simulation gets conditions for series**/
    si->cond = (CONDITION **) malloc (100 * sizeof (CONDITION *));

    for (int ct_se = 0; ct_se < si->anz_se; ++ct_se) {

        si->cond[ct_se] = si->all_cond[ct_all_cond ++ ];

        /*
        fprintf (stderr, "\nsi->cond[%d]             = si->all_cond[%d]   ", ct_se, ct_all_cond - 1);
        fprintf (stderr, "   op_tags = %d %d", si->cond[ct_se]->op_tag, si->all_cond[ct_all_cond-1]->op_tag);
        */

        /**series get conditions for sweeps**/
        si->se[ct_se].cond = (CONDITION **) malloc (10000 * sizeof (CONDITION *));

        for (int ct_sw = 0; ct_sw < si->se[ct_se].anz_sw; ++ct_sw) {

            si->se[ct_se].cond[ct_sw] = si->all_cond[ct_all_cond ++ ];

            /*
            fprintf (stderr, "\nsi->se[%d].cond[%d]       = si->all_cond[%d]   ", ct_se, ct_sw, ct_all_cond - 1);
            fprintf (stderr, "   op_tags = %d %d", si->se[ct_se].cond[ct_sw]->op_tag, si->all_cond[ct_all_cond-1]->op_tag);
            */

            /**sweeps get conditions for stays**/
            si->se[ct_se].sw[ct_sw].cond = (CONDITION **) malloc (10000 * sizeof (CONDITION *));

            for (int ct_st = 0; ct_st < si->se[ct_se].sw[ct_sw].anz_st; ++ct_st) {

                si->se[ct_se].sw[ct_sw].cond[ct_st] = si->all_cond[ct_all_cond ++ ];

                /*
                fprintf (stderr, "\nsi->se[%d].sw[%d].cond[%d] = si->all_cond[%d]   ", ct_se, ct_sw, ct_st, ct_all_cond - 1);
                fprintf (stderr, "   op_tags = %d %d", si->se[ct_se].sw[ct_sw].cond[ct_st]->op_tag, si->all_cond[ct_all_cond-1]->op_tag);
                */
            }
        }
    }
}

void printcondition (PARAMS *x, CONDITION *cond, FILE *fp, const int indent) {

            if  (cond->op_tag > 0)
                for (int i = 0; i < indent; ++i)
                    fprintf (fp, " ");

            if  (cond->op_tag == 0)
                ; /**nothing**/

            if  (cond->op_tag == 1) {
                if  (cond->left == &(x->iter))
                    fprintf (fp, "if (iter)\n");
                else
                    fprintf (fp, "if (%d)\n", cond->left[0]);
            }

            if  (cond->op_tag > 1) {
                fprintf (fp, "if (");
                if  (cond->left == &(x->iter))
                    fprintf (fp, "iter");
                else
                    fprintf (fp, "%d", cond->left[0]);
                if  (cond->op_tag == 2)
                    fprintf (fp, " = ");
                if  (cond->op_tag == 3)
                    fprintf (fp, " %% ");
                if  (cond->op_tag == 4)
                    fprintf (fp, " < ");
                if  (cond->op_tag == 5)
                    fprintf (fp, " > ");
                if  (cond->right == &(x->iter))
                    fprintf (fp, "iter");
                else
                    fprintf (fp, "%d", cond->right[0]);
                fprintf (fp, ")\n");
            }
}

void printparams (PARAMS *x, AREA *A, SIMULATION *si, FILE *fp) {
    int i, ii, ct_ar, ct_se, ct_sw, ct_st, ct_cmd, letters;
    SERIES  *se;
    COMMAND *cmd;

    #undef CONVERT_GLO
    #undef CONVERT_SCA
    #undef CONVERT_VEC
    #undef CONVERT_ARR
    #define CONVERT_GLO(a,b)  ;
    #define CONVERT_SCA(a,b)  ;
    #define CONVERT_VEC(a,b)  ;
    #define CONVERT_ARR(a,b)  ;

    /**global parameters**/
    #undef CONVERT_GLO
    #define CONVERT_GLO(a,b) \
        if  (!strcmp (#a, "int")) \
            fprintf (fp, " " #b " %d\n", (int)(x->b)); \
        if  (!strcmp (#a, "double")) \
            fprintf (fp, " " #b " %f\n", (double)(x->b));
    fprintf (fp, "\nglobal {\n");
    #include "parameters.h"
    #undef CONVERT_GLO
    #define CONVERT_GLO(a,b)  ;

    for (i = 0; i < x->num_pointers_at_global; ++i) {
        GLOBPTR *globptr = (GLOBPTR *)x->pointers[i];
        fprintf (fp, " %s \"", x->ch_pointers[i]);
        for (ii = 0; ii < globptr->num_words - 1; ++ii)
            fprintf (fp, "%s ", globptr->words[ii]);
        if  (globptr->num_words > 0)
            fprintf (fp, "%s", globptr->words[globptr->num_words - 1]);
        fprintf (fp, "\"\n");
    }

    fprintf (fp, "}\n");


    /**all area parameters**/
    #undef CONVERT_SCA
    #define CONVERT_SCA(a,b) \
        if  (!strcmp (#a, "int")) \
            fprintf (fp, " " #b " %d\n", (int)(A[x->areas].b)); \
        if  (!strcmp (#a, "double")) { \
            fprintf (fp, " " #b " "); \
            floatprint (fp, (double)(A[x->areas].b)); \
            fprintf (fp, "\n"); \
        }
    fprintf (fp, "\nall {\n");
    #include "parameters.h"
    fprintf (fp, "}\n");
    #undef CONVERT_SCA
    #define CONVERT_SCA(a,b)  ;


    /**single area parameters**/
    #undef CONVERT_SCA
    #define CONVERT_SCA(a,b) \
        if  (!strcmp (#a, "int")) \
            if  ((int)(A[x->areas].b) != (int)(A[ct_ar].b)) /*speudoarea*/ \
                fprintf (fp, " " #b " %d\n", (int)(A[ct_ar].b)); \
        if  (!strcmp (#a, "double")) \
            if  ((double)(A[x->areas].b) != (double)(A[ct_ar].b)) { \
                fprintf (fp, " " #b " "); \
                floatprint (fp, (double)(A[ct_ar].b)); \
                fprintf (fp, "\n"); \
            }

    for (ct_ar = 0; ct_ar < x->areas; ++ct_ar) {
        fprintf (fp, "\narea %d {\n", ct_ar);
        #include "parameters.h"
        fprintf (fp, "}\n");
    }
    #undef CONVERT_SCA
    #define CONVERT_SCA(a,b)  ;

    fprintf (fp, "\n");

    /**series**/
    for (ct_se = 0; ct_se < si->anz_se; ++ct_se) {

        /**condition at series**/
        printcondition (x, si->cond[ct_se], fp, 0);

        se = si->se+ct_se;
        fprintf (fp, "series (%d) {\n", se->ilen);

        /**sweeps**/
        for (ct_sw = 0; ct_sw < se->anz_sw; ++ct_sw) {

            /**condition at sweep**/
            printcondition (x, se->cond[ct_sw], fp, 1);

            /**control of each sweep**/
            fprintf (fp, " sw (%d; %d; %s) {\n", se->sw[ct_sw].begin, se->sw[ct_sw].end, se->sw[ct_sw].sw_update);

            /**stays**/
            for (ct_st = 0; ct_st < se->sw[ct_sw].anz_st; ++ct_st) {

                /**condition at stay**/
                printcondition (x, (se->sw+ct_sw)->cond[ct_st], fp, 2);

                fprintf (fp, "  %d(%c) ", se->sw[ct_sw].st[ct_st].area, se->sw[ct_sw].st[ct_st].st_update);

                /**commands**/
                for (ct_cmd=0; ct_cmd < se->sw[ct_sw].st[ct_st].anz_cmd; ++ct_cmd) {

                    cmd = se->sw[ct_sw].st[ct_st].cmd + ct_cmd;

                    /**print space if not 1st command**/
                    if  (ct_cmd)
                        fprintf (fp, "       ");

                    /**ch_target**/
                    fprintf(fp, "{%c", cmd->ch_target);

                    /**pointers**/
                    if (cmd->anz_pointers) {
                       fprintf (stderr, ", ");

                       for (i = 0; i < cmd->anz_pointers - 1; ++i) {
                           fprintf (stderr, "%s", cmd->ch_pointers[i]);
                           fprintf (stderr, "+");
                       }
                       fprintf (stderr, "%s", cmd->ch_pointers[cmd->anz_pointers - 1]);
                    }

                    /**func_name**/
                    fprintf(fp, "; %16s ; ", cmd->func_name);

                    /**n_from1**/
                    letters = 0;
                    for (i = 0; i < cmd->anz_from1; ++i) {
                        if  (i) {
                            fprintf (fp, "+");
                            letters += 1;
                        }
                        fprintf (fp, "%d", cmd->n_from1[i]);
                        letters += 1;
                    }

                    /**n_modu: consider later!**/

                    /**n_from2**/
                    if  (cmd->anz_arguments >= 2) {
                        fprintf (fp, ", ");
                        letters += 2;
                        for (i = 0; i < cmd->anz_from2; ++i) {
                            if  (i) {
                                fprintf (fp, "+");
                                letters += 1;
                            }
                            fprintf (fp, "%d", cmd->n_from2[i]);
                            letters += 1;
                        }
                    }

                    /**n_from3**/
                    if  (cmd->anz_arguments == 3) {
                        fprintf (fp, ", ");
                        letters += 2;
                        for (i = 0; i < cmd->anz_from3; ++i) {
                            if  (i) {
                                fprintf (fp, "+");
                                letters += 1;
                            }
                            fprintf (fp, "%d", cmd->n_from3[i]);
                            letters += 1;
                        }
                    }

                    if  (letters < 8)
                        for (i = 8; i > letters; --i)
                            fprintf (fp, " ");
                    fprintf (fp, " ; ");

                    /**ch_from1**/
                    letters = 0;
                    for (i = 0; i < cmd->anz_from1; ++i) {
                        if  (i) {
                            fprintf (fp, "+");
                            letters += 1;
                        }
                        fprintf (fp, "%c", cmd->ch_from1[i]);
                        letters += 1;
                    }

                    /**ch_from2**/
                    if  (cmd->anz_arguments >= 2) {
                        fprintf (fp, ", ");
                        letters += 2;
                        for (i = 0; i < cmd->anz_from2; ++i) {
                            if  (i) {
                                fprintf (fp, "+");
                                letters += 1;
                            }
                            fprintf (fp, "%c", cmd->ch_from2[i]);
                            letters += 1;
                        }
                    }

                    /**ch_from3**/
                    if  (cmd->anz_arguments == 3) {
                        fprintf (fp, ", ");
                        letters += 2;
                        for (i = 0; i < cmd->anz_from3; ++i) {
                            if  (i) {
                                fprintf (fp, "+");
                                letters += 1;
                            }
                            fprintf (fp, "%c", cmd->ch_from3[i]);
                            letters += 1;
                        }
                    }

                    if  (letters < 8)
                        for (i = 8; i > letters; --i)
                            fprintf (fp, " ");
                    fprintf (fp, " ; ");


                    for (i = 0; i < cmd->anz_quantums; ++i) {
                        if  (i) fprintf (fp, ", ");
                        for (ii = 0; ii < cmd->anz_quant[i]; ++ii) {
                            if  (ii)  fprintf (fp, "+");
                            if  (cmd->quantum[i][ii] != 0.0)
                                floatprint (fp, cmd->quantum[i][ii]);
                        }
                    }
                    fprintf (fp, "}\n");
                }
            }
            fprintf (fp, " }\n");
        }
        fprintf (fp, "}\n");
    }
}


void free_se (SIMULATION *si) {
    int i, ct_se, ct_sw, ct_st, ct_cmd;
    SERIES  *se;
    COMMAND *cmd;

    /**series**/
    for (ct_se = 0; ct_se < si->anz_se; ++ct_se) {

        se = si->se+ct_se;

        /**sweeps**/
        for (ct_sw = 0; ct_sw < se->anz_sw; ++ct_sw) {
            /**stays**/
            for (ct_st = 0; ct_st < se->sw[ct_sw].anz_st; ++ct_st) {
                /**commands**/
                for (ct_cmd=0; ct_cmd < se->sw[ct_sw].st[ct_st].anz_cmd; ++ct_cmd) {

                    cmd = se->sw[ct_sw].st[ct_st].cmd + ct_cmd;

                    free (cmd->func_name);

                    free (cmd->n_from1);
                    free (cmd->n_from2);
                    free (cmd->n_from3);

                    free (cmd->ch_from1);
                    free (cmd->ch_from2);
                    free (cmd->ch_from3);

                    for (i = 0; i < cmd->anz_quantums; ++i)
                        free (cmd->quantum[i]);
                    free (cmd->quantum);
                }
                free (se->sw[ct_sw].st[ct_st].cmd);
            }
            free (se->sw[ct_sw].st);
            free (se->sw[ct_sw].sw_update);
        }
        free (se->sw);
    }
    free (si->se);
}





/********************************* init_a ************************************/
/* Fills in contents of Areas A.                                             */

void init_a (AREA *A, PARAMS *x, SIMULATION *si) {
    AREA *a;
    int ct_ar;

    /*int sum_d_n = 0;*/
    /*int offset  = 0;*/

    /**essential within-one area values (must be ready for next loop)**/
    for (ct_ar = 0; ct_ar < x->areas; ++ct_ar) {

        a = A + ct_ar;

        a->d_n = a->d_a * a->d_b;

        a->shuffle = i_vector (a->d_n); /**the old version**reactivated, because new didn't work**/
        /*sum_d_n += a->d_n;*/
    }

    /*
    a = A + 0;
    a->shuffle = (int *)malloc (sum_d_n * sizeof (int));  ** the shuffle arrays  **
    for (ct_ar = 0; ct_ar < x->areas; ++ct_ar) {          **  of all areas are   **
        a = A + ct_ar;                                    **continuously arranged**
        a->shuffle = a->shuffle + offset;                 **  after each other   **
        offset += a->d_n;
    }
    */

    int max_rlen = 0;
    for (int ct_se = 0; ct_se < si->anz_se; ++ct_se) {
        SERIES *se = si->se + ct_se;
        for (int ct_sw = 0; ct_sw < se->anz_sw; ++ct_sw)
            max_rlen = se->sw[ct_sw].end > max_rlen
                     ? se->sw[ct_sw].end : max_rlen;
    }
    x->max_rlen = max_rlen;

    /**activations**/
    for (ct_ar = 0; ct_ar < x->areas; ++ct_ar) {

        a = A + ct_ar;

        a->A = d_matrix (x->mult * x->max_rlen, a->d_n);
        a->B = d_matrix (x->mult * x->max_rlen, a->d_n);
        a->C = d_matrix (x->mult * x->max_rlen, a->d_n);
        a->D = d_matrix (x->mult * x->max_rlen, a->d_n);
        a->E = d_matrix (x->mult * x->max_rlen, a->d_n);
        a->F = d_matrix (x->mult * x->max_rlen, a->d_n);
        a->G = d_matrix (x->mult * x->max_rlen, a->d_n);
        a->H = d_matrix (x->mult * x->max_rlen, a->d_n);
        a->I = d_matrix (x->mult * x->max_rlen, a->d_n);
        a->J = d_matrix (x->mult * x->max_rlen, a->d_n);
        a->K = d_matrix (x->mult * x->max_rlen, a->d_n);
        a->L = d_matrix (x->mult * x->max_rlen, a->d_n);
        a->M = d_matrix (x->mult * x->max_rlen, a->d_n);
        a->N = d_matrix (x->mult * x->max_rlen, a->d_n);
        a->O = d_matrix (x->mult * x->max_rlen, a->d_n);
        a->P = d_matrix (x->mult * x->max_rlen, a->d_n);
        a->Q = d_matrix (x->mult * x->max_rlen, a->d_n);
        a->R = d_matrix (x->mult * x->max_rlen, a->d_n);
        a->S = d_matrix (x->mult * x->max_rlen, a->d_n);
        a->T = d_matrix (x->mult * x->max_rlen, a->d_n);
        a->U = d_matrix (x->mult * x->max_rlen, a->d_n);
        a->V = d_matrix (x->mult * x->max_rlen, a->d_n);
        a->W = d_matrix (x->mult * x->max_rlen, a->d_n);
        a->X = d_matrix (x->mult * x->max_rlen, a->d_n);
        a->Y = d_matrix (x->mult * x->max_rlen, a->d_n);
        a->Z = d_matrix (x->mult * x->max_rlen, a->d_n);
    }
}


/******************************** free_a *************************************/
/* Frees memory for Areas A and contents.                                    */

void free_a (AREA *A, PARAMS *x, SIMULATION *si) {
    AREA  *a;
    int ct_ar;

    for (ct_ar = 0; ct_ar < x->areas; ++ct_ar) {

        a = A + ct_ar;

        /**free all activations here if cleanliness requested**/
    }

    free (A);
}



/***************************** init_weights_arbor ****************************/
void init_weights_arbor (double **S, int d_x, int d_y, int d_a, int d_b) {
    int a, b, x, y;

    double sigma = 0.2; /**relative to area dim's -> elliptic Gauss!                     changed from 0.25**/
    double diffA, diffB;
    double normfact = 1.0 / (2.0 * M_PI * sigma * sigma);

    for (x = 0; x < d_x; ++x)
    for (y = 0; y < d_y; ++y)
    for (a = 0; a < d_a; ++a)
    for (b = 0; b < d_b; ++b) {

        /**simple gradient (all below not needed!)**
          S[y + x * d_y][b + a * d_b]
        = (double)(a * x + b * y)
        / (double)((d_a-1) * (d_x-1) + (d_b-1) * (d_y-1)) * 0.1;
        */
        /**Gauss, periodic boundary, from p/c/difff7.c**/
        diffA = (a/(d_a-1.0) - x/(d_x-1.0)) < 0.5
             ? (a/(d_a-1.0) - x/(d_x-1.0)) : 1.0 - (a/(d_a-1.0) - x/(d_x-1.0));
        diffB = (b/(d_b-1.0) - y/(d_y-1.0)) < 0.5
             ? (b/(d_b-1.0) - y/(d_y-1.0)) : 1.0 - (b/(d_b-1.0) - y/(d_y-1.0));
        if  (diffA < -0.5)
              diffA = -1.0 - diffA;
        if  (diffB < -0.5)
              diffB = -1.0 - diffB;

        /**Gauss, clipped boundary***/
        diffA = (a/(d_a-1.0) - x/(d_x-1.0));
        diffB = (b/(d_b-1.0) - y/(d_y-1.0));
        

        if  ((d_a == 1) || (d_x == 1)) {
            diffA = 0;
            normfact = 1.0 / (sqrt (2.0 * M_PI) * sigma);
        }
        if  ((d_b == 1) || (d_y == 1)) {
            diffB = 0;
            normfact = 1.0 / (sqrt (2.0 * M_PI) * sigma);
        }

        /**needed only for Gauss!**/
        S[y + x * d_y][b + a * d_b] = diffA*diffA + diffB*diffB;
        S[y + x * d_y][b + a * d_b] = 0.02 * normfact                              /**changed from 0.01**/
                                    * exp (-0.5 * S[y + x * d_y][b + a * d_b]
                                                / (sigma * sigma))
                                  /** (1.0 - 2.0 * drand48())*/;                   /**commented out, to get smooth, positive-only weights**/
    }
     
}

/***************************** init_weights_arbor_color **********************/
void init_weights_arbor_color (double **S, int d_x, int d_y, int d_a, int d_b) {
    int a, b, x, y;

    double sigma = 0.2; /**relative to area dim's -> elliptic Gauss!                     changed from 0.25**/
    double diffA, diffB;
    double normfact = 1.0 / (2.0 * M_PI * sigma * sigma);

    d_a /= 3;

    for (x = 0; x < d_x; ++x)
    for (y = 0; y < d_y; ++y)
    for (a = 0; a < d_a; ++a)
    for (b = 0; b < d_b; ++b) {

        /**simple gradient (all below not needed!)**
          S[y + x * d_y][b + a * d_b]
        = (double)(a * x + b * y)
        / (double)((d_a-1) * (d_x-1) + (d_b-1) * (d_y-1)) * 0.1;
        */
        /**Gauss, periodic boundary, from p/c/difff7.c**/
        diffA = (a/(d_a-1.0) - x/(d_x-1.0)) < 0.5
             ? (a/(d_a-1.0) - x/(d_x-1.0)) : 1.0 - (a/(d_a-1.0) - x/(d_x-1.0));
        diffB = (b/(d_b-1.0) - y/(d_y-1.0)) < 0.5
             ? (b/(d_b-1.0) - y/(d_y-1.0)) : 1.0 - (b/(d_b-1.0) - y/(d_y-1.0));
        if  (diffA < -0.5)
              diffA = -1.0 - diffA;
        if  (diffB < -0.5)
              diffB = -1.0 - diffB;

        /**Gauss, clipped boundary***/
        diffA = (a/(d_a-1.0) - x/(d_x-1.0));
        diffB = (b/(d_b-1.0) - y/(d_y-1.0));
        

        if  ((d_a == 1) || (d_x == 1)) {
            diffA = 0;
            normfact = 1.0 / (sqrt (2.0 * M_PI) * sigma);
        }
        if  ((d_b == 1) || (d_y == 1)) {
            diffB = 0;
            normfact = 1.0 / (sqrt (2.0 * M_PI) * sigma);
        }

        /**needed only for Gauss!**/
        S[y + x * d_y][b + a * d_b] = diffA*diffA + diffB*diffB;
        S[y + x * d_y][b + a * d_b] = 0.02 * normfact                              /**changed from 0.01**/
                                    * exp (-0.5 * S[y + x * d_y][b + a * d_b]
                                                / (sigma * sigma))
                                    * (1.0 - 2.0 * drand48());
    }
     

    for (x = 0; x < d_x; ++x)
    for (y = 0; y < d_y; ++y)
    for (a = 0; a < d_a; ++a)
    for (b = 0; b < d_b; ++b) {
        S[y + x * d_y][b + (d_a + a) * d_b] = S[y + x * d_y][b + a * d_b];            /**green**/
        S[y + x * d_y][b + (2 * d_a + a) * d_b] = S[y + x * d_y][b + a * d_b];        /**blue**/
    }
}


/************************* assign_local_func *********************************/
/* Assigns specific function pointer; used by choose_functions.              */

void assign_local_func (COMMAND *cmd, DOUBLE (*localfunc)(DOUBLE *par, DOUBLE val1, DOUBLE val2), char *name) {

    if  (!strcmp (cmd->func_name, name)) {
        cmd->localfunc = localfunc;
        cmd->func_key = 'l';
    }
}

/************************* assign_func ***************************************/

void assign_func (COMMAND *cmd, DOUBLE (*func)(PARAMS *x, AREA *A, COMMAND *cmd, int ival1, int ival2), char *name, char key) {

    if  (!strcmp (cmd->func_name, name)) {
        cmd->func = func;
        cmd->func_key = key;

        /* fprintf (stderr, "\nassigning %s with key %c  ", name, key); */
    }
}


/*************************** choose_functions ********************************/
/* Assigns function pointers of x and s according to their names.            */

void choose_functions (PARAMS *x, SIMULATION *si) {

    int ct_se, ct_sw, ct_st, ct_cmd, dummy;
    SERIES  *se;
    COMMAND *cmd;

ERR("\nchoose_functions");

    /**data import functions**
    int OK = 0;
    if  (x->data_import_name == NULL) {
        x->data_import = NULL;
        OK = 1;
    } else {
        if  (!strcmp (x->data_import_name, "import_images")) {
            x->data_import = import_images;
            OK = 1;
        }
        if  (!strcmp (x->data_import_name, "import_points")) {
            x->data_import = import_points;
            OK = 1;
        }
        if  (!strcmp (x->data_import_name, "import_motor")) {
            x->data_import = import_motor;
            OK = 1;
        }
        if  (!strcmp (x->data_import_name, "import_lang_assoc")) {
            x->data_import = import_lang_assoc;
            OK = 1;
        }
    }
    if  (!OK)
        fprintf (stderr, "\n\ndata import function could not be assigned!\n");
    **/


ERR("functions_middle");

    /**series**/
    for (ct_se = 0; ct_se < si->anz_se; ++ct_se) {

      se = si->se + ct_se;

      /**sweeps**/
      for (ct_sw = 0; ct_sw < se->anz_sw; ++ct_sw) {

        /**stays**/
        for (ct_st = 0; ct_st < se->sw[ct_sw].anz_st; ++ct_st) {

            /**commands**/
            for (ct_cmd=0; ct_cmd < se->sw[ct_sw].st[ct_st].anz_cmd; ++ct_cmd) {

                cmd = se->sw[ct_sw].st[ct_st].cmd + ct_cmd;

                cmd->func_key = '0';

                /**assign local functions -- func_key is always 'l'**/

                assign_local_func (cmd, local_sum, "local_sum");
                assign_local_func (cmd, local_sum_const, "local_sum_const");
                assign_local_func (cmd, local_sub, "local_sub");
                assign_local_func (cmd, local_mult, "local_mult");
                assign_local_func (cmd, local_mult_const, "local_mult_const");
                assign_local_func (cmd, local_div, "local_div");
                assign_local_func (cmd, local_div_const, "local_div_const");
                assign_local_func (cmd, local_modulo_const, "local_modulo_const");
                assign_local_func (cmd, local_cyclic, "local_cyclic");
                assign_local_func (cmd, local_average, "local_average");
                assign_local_func (cmd, local_copy, "local_copy");
                assign_local_func (cmd, local_persist, "local_persist");
                assign_local_func (cmd, local_threshold, "local_threshold");
                assign_local_func (cmd, local_isbigger, "local_isbigger");
                assign_local_func (cmd, local_biggerof, "local_biggerof");
                assign_local_func (cmd, local_biggerfabsof, "local_biggerfabsof");
                assign_local_func (cmd, local_const, "local_const");
                assign_local_func (cmd, local_rectify, "local_rectify");
                assign_local_func (cmd, local_invert, "local_invert");
                assign_local_func (cmd, local_abs, "local_abs");
                assign_local_func (cmd, local_sign, "local_sign");
                assign_local_func (cmd, local_round, "local_round");
                assign_local_func (cmd, local_sqrt, "local_sqrt");
                assign_local_func (cmd, local_power, "local_power");
                assign_local_func (cmd, local_exp, "local_exp");
                assign_local_func (cmd, local_log_pos, "local_log_pos");
                assign_local_func (cmd, local_cos, "local_cos");
                assign_local_func (cmd, local_sin, "local_sin");
                assign_local_func (cmd, local_tan, "local_tan");
                assign_local_func (cmd, local_atan, "local_atan");
                assign_local_func (cmd, local_atan_to_2pi, "local_atan_to_2pi");
                assign_local_func (cmd, local_lin_01, "local_lin_01");
                assign_local_func (cmd, local_rand, "local_rand");
                assign_local_func (cmd, local_rand_pos_neg, "local_rand_pos_neg");
                assign_local_func (cmd, local_as_rand, "local_as_rand");
                assign_local_func (cmd, local_circ_gauss, "local_circ_gauss");
                assign_local_func (cmd, local_rand_gibbs, "local_rand_gibbs");
                assign_local_func (cmd, local_rand_gibbs_01, "local_rand_gibbs_01");
                assign_local_func (cmd, local_rand_exp, "local_rand_exp");
                assign_local_func (cmd, local_tanh, "local_tanh");
                assign_local_func (cmd, local_bp_tanh, "local_bp_tanh");
                assign_local_func (cmd, local_sparse, "local_sparse");
                assign_local_func (cmd, local_sparse_diff, "local_sparse_diff");
                assign_local_func (cmd, local_sparse_01, "local_sparse_01");
                assign_local_func (cmd, local_zhang, "local_zhang");
                assign_local_func (cmd, local_zhang_scale, "local_zhang_scale");
                assign_local_func (cmd, local_zhang2, "local_zhang2");
                assign_local_func (cmd, local_zhang_inv, "local_zhang_inv");
                assign_local_func (cmd, local_zhang2_inv, "local_zhang2_inv");
                assign_local_func (cmd, local_zhang_diff, "local_zhang_diff");
                assign_local_func (cmd, local_zhang2_diff, "local_zhang2_diff");
                assign_local_func (cmd, local_mean, "local_mean");
                assign_local_func (cmd, local_mean_diff, "local_mean_diff");
                assign_local_func (cmd, local_mean_01, "local_mean_01");
                assign_local_func (cmd, local_mean_scale, "local_mean_scale");
                assign_local_func (cmd, local_mean_scale_diff, "local_mean_scale_diff");
                assign_local_func (cmd, local_bp_mean_01, "local_bp_mean_01");
                assign_local_func (cmd, local_mean_01_inv, "local_mean_01_inv");
                assign_local_func (cmd, local_mean_01_diff, "local_mean_01_diff");
                assign_local_func (cmd, local_mean_01_d0, "local_mean_01_d0");
                assign_local_func (cmd, local_mean_01_d1, "local_mean_01_d1");
                assign_local_func (cmd, local_mean_01_d2, "local_mean_01_d2");
                assign_local_func (cmd, local_mean_01_d3, "local_mean_01_d3");
                assign_local_func (cmd, local_mean_01_d4, "local_mean_01_d4");
                assign_local_func (cmd, local_gibbs, "local_gibbs");
                assign_local_func (cmd, local_gibbs_01, "local_gibbs_01");
                assign_local_func (cmd, local_gibbs_01_2, "local_gibbs_01_2");
                assign_local_func (cmd, local_odds_01, "local_odds_01");

                /**assign other functions -- func_key: 's'=single, 't'=total, 'a'=alltime;  'v'=versatile (total or alltime), 'n'=nevermind (single or total or alltime)**/

                assign_func  (cmd, single_ptr2act, "single_ptr2act", 's');
                assign_func  (cmd, set_ptr_int_val, "set_ptr_int_val", 's');
                assign_func  (cmd, set_ptr_val, "set_ptr_val", 's');
                assign_func  (cmd, single_copy, "single_copy", 's');
                assign_func  (cmd, single_copy_limited, "single_copy_limited", 's');
                assign_func  (cmd, single_mean_back, "single_mean_back", 's');
                assign_func  (cmd, single_circ_gauss, "single_circ_gauss", 's');
                assign_func  (cmd, single_l_add, "single_l_add", 's');

                assign_func  (cmd, weight_list_alloc_full, "weight_list_alloc_full", 'a');
                assign_func  (cmd, weight_list_alloc_topo, "weight_list_alloc_topo", 'a');
                assign_func  (cmd, weight_list_alloc_invert, "weight_list_alloc_invert", 'a');
                assign_func  (cmd, weight_list_invert_quick, "weight_list_invert_quick", 'a');
                assign_func  (cmd, weight_list_mult, "weight_list_mult", 'a');
                assign_func  (cmd, weight_list_free, "weight_list_free", 'a');
                assign_func  (cmd, weight_list_feed, "weight_list_feed", 's');
                assign_func  (cmd, WEIGHT_LIST_PUSH, "WEIGHT_LIST_PUSH", 's');
                assign_func  (cmd, weight_list_AXON_total, "weight_list_AXON_total", 't');
                assign_func  (cmd, weight_list_AXON_total_pos, "weight_list_AXON_total_pos", 't');
                assign_func  (cmd, weight_list_hebb, "weight_list_hebb", 's');
                assign_func  (cmd, weight_list_hebb_diff, "weight_list_hebb_diff", 's');
                assign_func  (cmd, weight_list_ica1, "weight_list_ica1", 't');
                assign_func  (cmd, weight_list_ica2, "weight_list_ica2", 't');
                assign_func  (cmd, weight_list_kohonen, "weight_list_kohonen", 's');
                assign_func  (cmd, weight_list_euclid, "weight_list_euclid", 's');
                assign_func  (cmd, weight_list_euclid_ptr, "weight_list_euclid_ptr", 's');
                assign_func  (cmd, weight_list_act2weight, "weight_list_act2weight", 's');
                assign_func  (cmd, weight_list_act2weight_offset, "weight_list_act2weight_offset", 'v');
                assign_func  (cmd, weight_list_hebb_turn, "weight_list_hebb_turn", 's');
                assign_func  (cmd, weight_list_decay, "weight_list_decay", 's');
                assign_func  (cmd, weight_list_decay_const, "weight_list_decay_const", 's');
                assign_func  (cmd, weight_list_decay_quad, "weight_list_decay_quad", 's');
                assign_func  (cmd, weight_list_decay_post, "weight_list_decay_post", 's');
                assign_func  (cmd, weight_list_decay_pre, "weight_list_decay_pre", 's');
                assign_func  (cmd, weight_list_decay_Baddeley, "weight_list_decay_Baddeley", 's');
                assign_func  (cmd, weight_list_decay_Baddeley_distance, "weight_list_decay_Baddeley_distance", 's');
                assign_func  (cmd, weight_list_normalize, "weight_list_normalize", 's');
                assign_func  (cmd, weight_list_rectify, "weight_list_rectify", 's');
                assign_func  (cmd, weight_list_cutself, "weight_list_cutself", 's');
                assign_func  (cmd, weight_list_cutsmall, "weight_list_cutsmall", 's');
                assign_func  (cmd, weight_list_sprout, "weight_list_sprout", 's');
                assign_func  (cmd, weight_list_export, "weight_list_export", 'a');
                assign_func  (cmd, weight_list_export_col, "weight_list_export_col", 'a');
                assign_func  (cmd, weight_list_export_gnu, "weight_list_export_gnu", 'a');
                assign_func  (cmd, weight_list_alloc_import, "weight_list_alloc_import", 'a');
                assign_func  (cmd, weight_list_alloc_import_tiles, "weight_list_alloc_import_tiles", 'a');
                assign_func  (cmd, weight_list_histogram, "weight_list_histogram", 'a');
                assign_func  (cmd, weight_list_cuthalfinput, "weight_list_cuthalfinput", 's');

                assign_func  (cmd, observe_phase, "observe_phase", 'a');
                assign_func  (cmd, import_phase, "import_phase", 'a');

                assign_func  (cmd, observe_act, "observe_act", 'a');
                assign_func  (cmd, observe_col, "observe_col", 'a');
                assign_func  (cmd, observe_animgif, "observe_animgif", 'a');
                assign_func  (cmd, observe_animgif_iter, "observe_animgif_iter", 'a');
                assign_func  (cmd, observe_getc, "observe_getc", 'n');
                assign_func  (cmd, observe_gnu_act, "observe_gnu_act", 'n');
                assign_func  (cmd, observe_gnu_two_acts, "observe_gnu_two_acts", 'n');
                assign_func  (cmd, observe_act_hist, "observe_act_hist", 'v');
                assign_func  (cmd, observe_phase_hist, "observe_phase_hist", 'v');

                assign_func  (cmd, weight_sigpi_alloc_full, "weight_sigpi_alloc_full", 'a');
                assign_func  (cmd, weight_sigpi_feed, "weight_sigpi_feed", 's');
                assign_func  (cmd, weight_sigpi_euclid, "weight_sigpi_euclid", 's');
                assign_func  (cmd, weight_sigpi_kohonen, "weight_sigpi_kohonen", 's');
                assign_func  (cmd, weight_sigpi_histogram, "weight_sigpi_histogram", 'a');
                assign_func  (cmd, weight_sigpi_hebb, "weight_sigpi_hebb", 's');
                assign_func  (cmd, weight_sigpi_cutsmall, "weight_sigpi_cutsmall", 's');
                assign_func  (cmd, weight_sigpi_export, "weight_sigpi_export", 'a');
                assign_func  (cmd, weight_sigpi_alloc_import, "weight_sigpi_alloc_import", 'a');
                assign_func  (cmd, weight_sigpi_gnuplot_cm, "weight_sigpi_gnuplot_cm", 'a');

                assign_func  (cmd, data_gauss_move, "data_gauss_move", 'a');
                assign_func  (cmd, data_gauss_three, "data_gauss_three", 'a');
                assign_func  (cmd, data_gauss_3areas, "data_gauss_3areas", 'a');
                assign_func  (cmd, data_gauss_3areas_2D, "data_gauss_3areas_2D", 'a');
                assign_func  (cmd, data_gauss_fromfile, "data_gauss_fromfile", 'a');
                assign_func  (cmd, data_gauss_3areas_2D_anim, "data_gauss_3areas_2D_anim", 'a');
                assign_func  (cmd, data_gauss_3areas_2D_anim2, "data_gauss_3areas_2D_anim2", 'a');
                assign_func  (cmd, data_gauss_2D_anim3, "data_gauss_2D_anim3", 'a');
                assign_func  (cmd, data_half_circle, "data_half_circle", 'a');
                assign_func  (cmd, write_act_file, "write_act_file", 'v');
                assign_func  (cmd, read_act_file, "read_act_file", 'v');
                assign_func  (cmd, dense_act_file, "dense_act_file", 'v');
                assign_func  (cmd, data_act_file, "data_act_file", 'v');
                assign_func  (cmd, init_lines_hier, "init_lines_hier", 'a');
                assign_func  (cmd, file_flag, "file_flag", 'v');
                assign_func  (cmd, int_val_change, "int_val_change", 'v');
                assign_func  (cmd, data_gnu_connect_max, "data_gnu_connect_max", 'a');

                assign_func  (cmd, import_images, "import_images", 'a');
                assign_func  (cmd, init_image, "init_image", 'a');
                assign_func  (cmd, init_whole_image, "init_whole_image", 'a');
                assign_func  (cmd, init_onoff_image, "init_onoff_image", 'a');
                assign_func  (cmd, cut_image, "cut_image", 'a');
                assign_func  (cmd, cut_image_pantilt, "cut_image_pantilt", 'a');
                assign_func  (cmd, init_orange, "init_orange", 'a');
                assign_func  (cmd, image_color_blob, "image_color_blob", 'v');
                assign_func  (cmd, init_image_cosinus, "init_image_cosinus", 'v');


                assign_func  (cmd, total_cut_peaks, "total_cut_peaks", 't');
                assign_func  (cmd, total_rectify_save, "total_rectify_save", 't');
                assign_func  (cmd, total_sub_mean, "total_sub_mean", 't');
                assign_func  (cmd, total_sub_mean_col, "total_sub_mean_col", 't');
                assign_func  (cmd, total_average, "total_average", 't');
                assign_func  (cmd, total_sub_min, "total_sub_min", 't');
                assign_func  (cmd, total_set_mean, "total_set_mean", 't');
                assign_func  (cmd, total_time_mean, "total_time_mean", 't');
                assign_func  (cmd, all_time_mean, "all_time_mean", 'a');
                assign_func  (cmd, total_set_attime, "total_set_attime", 't');
                assign_func  (cmd, total_epoche_mean, "total_epoche_mean", 't');
                assign_func  (cmd, total_epoche_x_square, "total_epoche_x_square", 't');
                assign_func  (cmd, total_epoche_x_fourth, "total_epoche_x_fourth", 't');
                assign_func  (cmd, total_epoche_kurt, "total_epoche_kurt", 't');
              /*assign_func  (cmd, total_import, "total_import", 't');*/
              /*assign_func  (cmd, total_back_pipe, "total_back_pipe", 't');*/
                assign_func  (cmd, total_copy_back, "total_copy_back", 't');
                assign_func  (cmd, total_embed, "total_embed", 't');
                assign_func  (cmd, total_embed_tiles, "total_embed_tiles", 't');
                assign_func  (cmd, total_spread_time, "total_spread_time", 't');
                assign_func  (cmd, total_spread_xy, "total_spread_xy", 't');
                assign_func  (cmd, total_set_all, "total_set_all", 't');
                assign_func  (cmd, total_collapse_xy, "total_collapse_xy", 't');
                assign_func  (cmd, total_mean_to_one, "total_mean_to_one", 't');
                assign_func  (cmd, total_spherize, "total_spherize", 't');
                assign_func  (cmd, total_normalize, "total_normalize", 't');
                assign_func  (cmd, total_scale_to_max, "total_scale_to_max", 't');
                assign_func  (cmd, total_threshold, "total_threshold", 't');
                assign_func  (cmd, total_true_softmax, "total_true_softmax", 't');
                assign_func  (cmd, total_true_softmax_row, "total_true_softmax_row", 't');
                assign_func  (cmd, total_four, "total_four", 't');
                assign_func  (cmd, total_mult_compl, "total_mult_compl", 't');
                assign_func  (cmd, total_absq_compl, "total_absq_compl", 't');
                assign_func  (cmd, total_geom_vector_sum, "total_geom_vector_sum", 't');
              /*assign_func  (cmd, total_zhang_m_fix, "total_zhang_m_fix", 't');*/
                assign_func  (cmd, total_orange, "total_orange", 't');
                assign_func  (cmd, total_cosinus, "total_cosinus", 't');
                assign_func  (cmd, total_gabor, "total_gabor", 't');
                assign_func  (cmd, total_gauss_elliptic, "total_gauss_elliptic", 't');
                assign_func  (cmd, total_cos, "total_cos", 't');
                assign_func  (cmd, total_cut_at, "total_cut_at", 't');
                assign_func  (cmd, total_shift_torus, "total_shift_torus", 't');
                assign_func  (cmd, total_winner, "total_winner", 't');
                assign_func  (cmd, total_winner_per_row, "total_winner_per_row", 't');
                assign_func  (cmd, total_neigh_winner, "total_neigh_winner", 't');
                assign_func  (cmd, total_neigh_winner_3D, "total_neigh_winner_3D", 't');
                assign_func  (cmd, total_fit_gauss, "total_fit_gauss", 't');
                assign_func  (cmd, total_gauss_at, "total_gauss_at", 't');
                assign_func  (cmd, total_softmax, "total_softmax", 't');
                assign_func  (cmd, total_error_at, "total_error_at", 't');
                assign_func  (cmd, gnuplot_error_at, "gnuplot_error_at", 'v');
                assign_func  (cmd, total_dist_xy_middle, "total_dist_xy_middle", 't');
                assign_func  (cmd, total_rand, "total_rand", 't');
                assign_func  (cmd, total_as_rand, "total_as_rand", 't');
                assign_func  (cmd, total_any_nonzero, "total_any_nonzero", 't');
                assign_func  (cmd, total_pause, "total_pause", 't');

                assign_func  (cmd, total_init_behaviour, "total_init_behaviour", 't');
                assign_func  (cmd, total_high_vision, "total_high_vision", 't');
                assign_func  (cmd, total_mult_const_if, "total_mult_const_if", 't');
                assign_func  (cmd, total_set_nth, "total_set_nth", 't');
                assign_func  (cmd, total_swap_nth, "total_swap_nth", 't');
                assign_func  (cmd, total_copy_nth, "total_copy_nth", 't');
                assign_func  (cmd, total_rand_nth, "total_rand_nth", 't');
                assign_func  (cmd, total_compare_with, "total_compare_with", 't');
                assign_func  (cmd, total_compare_count, "total_compare_count", 't');
                assign_func  (cmd, total_retina_to_SC, "total_retina_to_SC", 't');
                assign_func  (cmd, total_population_motor_row, "total_population_motor_row", 't');
                assign_func  (cmd, total_population_motor_2D, "total_population_motor_2D", 't');
                assign_func  (cmd, total_scalar_mult, "total_scalar_mult", 't');
                assign_func  (cmd, total_change_active_by, "total_change_active_by", 't');
                assign_func  (cmd, total_converge_to, "total_converge_to", 't');
                assign_func  (cmd, total_mean_left_min_right, "total_mean_left_min_right", 't');

                assign_func  (cmd, total_topo_gibbs_01, "total_topo_gibbs_01", 't');

                /**special/reinforcement.c**/
                assign_func  (cmd, feed_reward, "feed_reward", 's');
                assign_func  (cmd, feed_angl_reward, "feed_angl_reward", 's');
                assign_func  (cmd, feed_dock_reward, "feed_dock_reward", 's');
                assign_func  (cmd, total_rand_winner, "total_rand_winner", 't');
                assign_func  (cmd, total_init_place, "total_init_place", 't');
                assign_func  (cmd, total_move_place, "total_move_place", 't');
                assign_func  (cmd, total_init_coord, "total_init_coord", 't');
                assign_func  (cmd, total_scan_coord, "total_scan_coord", 't');
                assign_func  (cmd, total_imag_coord, "total_imag_coord", 't');
                assign_func  (cmd, total_coord_imag, "total_coord_imag", 't');
                assign_func  (cmd, total_angl_coord, "total_angl_coord", 't');
                assign_func  (cmd, total_angl_where, "total_angl_where", 't');
                assign_func  (cmd, total_xyp_to_dtp, "total_xyp_to_dtp", 't');
                assign_func  (cmd, total_move_coord, "total_move_coord", 't');
                assign_func  (cmd, total_behaviour, "total_behaviour", 't');
                assign_func  (cmd, total_real_world, "total_real_world", 't');
#if GAZEBO
                /**fun_gazebo/gazebo.io.c**/
                assign_func  (cmd, gazebo_d_image, "gazebo_d_image", 'v');
                assign_func  (cmd, gazebo_t_odometry, "gazebo_t_odometry", 't');
                assign_func  (cmd, gazebo_t_move, "gazebo_t_move", 't');
                assign_func  (cmd, gazebo_t_pantilt, "gazebo_t_pantilt", 't');
                assign_func  (cmd, gazebo_t_gripper, "gazebo_t_gripper", 't');
                assign_func  (cmd, gazebo_t_reset, "gazebo_t_reset", 't');
                assign_func  (cmd, gazebo_t_reward, "gazebo_t_reward", 't');
                assign_func  (cmd, gazebo_t_table_touch, "gazebo_t_table_touch", 't');
#endif

#if MIRO
                /**fun_miro/miro.io.c**/
                assign_func  (cmd, miro_d_sonar, "miro_d_sonar", 'a');
                assign_func  (cmd, miro_d_compass, "miro_d_compass", 'a');
                assign_func  (cmd, miro_d_image, "miro_d_image", 'a');
                assign_func  (cmd, miro_t_pantilt, "miro_t_pantilt", 'v');
                assign_func  (cmd, miro_do_motor, "miro_do_motor", 'a');
                assign_func  (cmd, miro_t_move, "miro_t_move", 'v');
                assign_func  (cmd, miro_d_limp, "miro_d_limp", 'a');
                assign_func  (cmd, miro_t_limp_if, "miro_t_limp_if", 'v');
                assign_func  (cmd, miro_t_openGrip, "miro_t_openGrip", 'v');
                assign_func  (cmd, miro_t_odometry, "miro_t_odometry", 't');
                assign_func  (cmd, miro_t_close, "miro_t_close", 'v');
                /**fun_miro/mirror.io.c**/
                assign_func  (cmd, setObjectVector_xy, "setObjectVector_xy", 'v');
                assign_func  (cmd, getObjectVector_xy_angle, "getObjectVector_xy_angle", 'v');
#endif

                /**something not OK**/
                if  (cmd->func_key == '0') {
                    fprintf (stderr, "\nno cmd->func ");
                    fprintf (stderr, "\"%s\" assigned!\n", cmd->func_name);
                }

                /**test usage of local and single functions**/
                if  ((cmd->func_key == 'l') || (cmd->func_key == 's'))
                    if  ( !(  (!strcmp (se->sw[ct_sw].sw_update, "order"))
                           || (!strcmp (se->sw[ct_sw].sw_update, "random"))
                           || (!strcmp (se->sw[ct_sw].sw_update, "propto"))
                           || (1 == sscanf (se->sw[ct_sw].sw_update, "%d", &dummy))))
                        fprintf (stderr, "\nWARNING: local/single function %s not in any of order/random/propto sweep!\n", cmd->func_name);

                if  ((cmd->func_key == 'l') || (cmd->func_key == 's'))
                    if  (se->sw[ct_sw].st[ct_st].st_update == 't')
                        fprintf (stderr, "\nWARNING: local/single function %s is in total stay!\n", cmd->func_name);

                /**test usage of total functions**/
                if  (cmd->func_key == 't')
                    if  (se->sw[ct_sw].st[ct_st].st_update != 't')
                        fprintf (stderr, "\nWARNING: total function %s not in total stay!\n", cmd->func_name);

                /**test usage of alltime functions (they are also allowed in a total stay if only one relaxation is done)**/
                if  (cmd->func_key == 'a')
                    if  (! ((!strcmp (se->sw[ct_sw].sw_update, "alltime") || ((se->sw[ct_sw].st[ct_st].st_update == 't') && (se->sw[ct_sw].end - se->sw[ct_sw].begin == 1)))))
                        fprintf (stderr, "\nWARNING: alltime function %s not in alltime sweep and not exactly 1 relaxation done!\n", cmd->func_name);

            } /**commands**/
        } /**stays**/
      } /**sweeps**/
    } /**series**/

ERR("functions_end");
}


/*************************** choose_pointers *********************************/
/* Assign pointers to values within cmd.                                     */
/* Also: allocate mem for *** ptrs S_from(2).                                */

void choose_pointers (SIMULATION *si, AREA *A) {

    int ct_se, ct_sw, ct_st, ct_cmd, ct_l, inarea;
    SERIES  *se;
    COMMAND *cmd;
    int OK = 0;

    /**series**/
    for (ct_se = 0; ct_se < si->anz_se; ++ct_se) {

      se = si->se + ct_se;

      /**sweeps**/
      for (ct_sw = 0; ct_sw < se->anz_sw; ++ct_sw) {
        /**stays**/
        for (ct_st = 0; ct_st < se->sw[ct_sw].anz_st; ++ct_st) {
            /**commands**/
            for (ct_cmd=0; ct_cmd < se->sw[ct_sw].st[ct_st].anz_cmd; ++ct_cmd) {

                cmd = se->sw[ct_sw].st[ct_st].cmd + ct_cmd;

                cmd->area   = se->sw[ct_sw].st[ct_st].area;
                cmd->moment = 1.0 / (double)(se->sw[ct_sw].end - se->sw[ct_sw].begin);

                OK = 0;
                /**target activations**/
                if  (cmd->ch_target == 'A') {
                    cmd->S_target = A[cmd->area].A;
                    OK = 1;
                }
                if  (cmd->ch_target == 'B') {
                    cmd->S_target = A[cmd->area].B;
                    OK = 1;
                }
                if  (cmd->ch_target == 'C') {
                    cmd->S_target = A[cmd->area].C;
                    OK = 1;
                }
                if  (cmd->ch_target == 'D') {
                    cmd->S_target = A[cmd->area].D;
                    OK = 1;
                }
                if  (cmd->ch_target == 'E') {
                    cmd->S_target = A[cmd->area].E;
                    OK = 1;
                }
                if  (cmd->ch_target == 'F') {
                    cmd->S_target = A[cmd->area].F;
                    OK = 1;
                }
                if  (cmd->ch_target == 'G') {
                    cmd->S_target = A[cmd->area].G;
                    OK = 1;
                }
                if  (cmd->ch_target == 'H') {
                    cmd->S_target = A[cmd->area].H;
                    OK = 1;
                }
                if  (cmd->ch_target == 'I') {
                    cmd->S_target = A[cmd->area].I;
                    OK = 1;
                }
                if  (cmd->ch_target == 'J') {
                    cmd->S_target = A[cmd->area].J;
                    OK = 1;
                }
                if  (cmd->ch_target == 'K') {
                    cmd->S_target = A[cmd->area].K;
                    OK = 1;
                }
                if  (cmd->ch_target == 'L') {
                    cmd->S_target = A[cmd->area].L;
                    OK = 1;
                }
                if  (cmd->ch_target == 'M') {
                    cmd->S_target = A[cmd->area].M;
                    OK = 1;
                }
                if  (cmd->ch_target == 'N') {
                    cmd->S_target = A[cmd->area].N;
                    OK = 1;
                }
                if  (cmd->ch_target == 'O') {
                    cmd->S_target = A[cmd->area].O;
                    OK = 1;
                }
                if  (cmd->ch_target == 'P') {
                    cmd->S_target = A[cmd->area].P;
                    OK = 1;
                }
                if  (cmd->ch_target == 'Q') {
                    cmd->S_target = A[cmd->area].Q;
                    OK = 1;
                }
                if  (cmd->ch_target == 'R') {
                    cmd->S_target = A[cmd->area].R;
                    OK = 1;
                }
                if  (cmd->ch_target == 'S') {
                    cmd->S_target = A[cmd->area].S;
                    OK = 1;
                }
                if  (cmd->ch_target == 'T') {
                    cmd->S_target = A[cmd->area].T;
                    OK = 1;
                }
                if  (cmd->ch_target == 'U') {
                    cmd->S_target = A[cmd->area].U;
                    OK = 1;
                }
                if  (cmd->ch_target == 'V') {
                    cmd->S_target = A[cmd->area].V;
                    OK = 1;
                }
                if  (cmd->ch_target == 'W') {
                    cmd->S_target = A[cmd->area].W;
                    OK = 1;
                }
                if  (cmd->ch_target == 'X') {
                    cmd->S_target = A[cmd->area].X;
                    OK = 1;
                }
                if  (cmd->ch_target == 'Y') {
                    cmd->S_target = A[cmd->area].Y;
                    OK = 1;
                }
                if  (cmd->ch_target == 'Z') {
                    cmd->S_target = A[cmd->area].Z;
                    OK = 1;
                }
                if  (!OK)
                    fprintf (stderr, "\n\nch_target could not be assigned!\n");

                /**alloc mem for ptrs that are assigned dereferenced**/
                cmd->S_from1 = (DOUBLE ***) malloc (cmd->anz_from1 * sizeof (DOUBLE **));
                cmd->S_from2 = (DOUBLE ***) malloc (cmd->anz_from2 * sizeof (DOUBLE **));
                cmd->S_from3 = (DOUBLE ***) malloc (cmd->anz_from3 * sizeof (DOUBLE **));

                /**sources**/
                for (ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

                    inarea = cmd->n_from1[ct_l];

                    OK = 0;
                    /**source activations; each shows to inarea-activations**/
                    if  (cmd->ch_from1[ct_l] == 'A') {
                        cmd->S_from1[ct_l] = A[inarea].A;
                        OK = 1;
                    }
                    if  (cmd->ch_from1[ct_l] == 'B') {
                        cmd->S_from1[ct_l] = A[inarea].B;
                        OK = 1;
                    }
                    if  (cmd->ch_from1[ct_l] == 'C') {
                        cmd->S_from1[ct_l] = A[inarea].C;
                        OK = 1;
                    }
                    if  (cmd->ch_from1[ct_l] == 'D') {
                        cmd->S_from1[ct_l] = A[inarea].D;
                        OK = 1;
                    }
                    if  (cmd->ch_from1[ct_l] == 'E') {
                        cmd->S_from1[ct_l] = A[inarea].E;
                        OK = 1;
                    }
                    if  (cmd->ch_from1[ct_l] == 'F') {
                        cmd->S_from1[ct_l] = A[inarea].F;
                        OK = 1;
                    }
                    if  (cmd->ch_from1[ct_l] == 'G') {
                        cmd->S_from1[ct_l] = A[inarea].G;
                        OK = 1;
                    }
                    if  (cmd->ch_from1[ct_l] == 'H') {
                        cmd->S_from1[ct_l] = A[inarea].H;
                        OK = 1;
                    }
                    if  (cmd->ch_from1[ct_l] == 'I') {
                        cmd->S_from1[ct_l] = A[inarea].I;
                        OK = 1;
                    }
                    if  (cmd->ch_from1[ct_l] == 'J') {
                        cmd->S_from1[ct_l] = A[inarea].J;
                        OK = 1;
                    }
                    if  (cmd->ch_from1[ct_l] == 'K') {
                        cmd->S_from1[ct_l] = A[inarea].K;
                        OK = 1;
                    }
                    if  (cmd->ch_from1[ct_l] == 'L') {
                        cmd->S_from1[ct_l] = A[inarea].L;
                        OK = 1;
                    }
                    if  (cmd->ch_from1[ct_l] == 'M') {
                        cmd->S_from1[ct_l] = A[inarea].M;
                        OK = 1;
                    }
                    if  (cmd->ch_from1[ct_l] == 'N') {
                        cmd->S_from1[ct_l] = A[inarea].N;
                        OK = 1;
                    }
                    if  (cmd->ch_from1[ct_l] == 'O') {
                        cmd->S_from1[ct_l] = A[inarea].O;
                        OK = 1;
                    }
                    if  (cmd->ch_from1[ct_l] == 'P') {
                        cmd->S_from1[ct_l] = A[inarea].P;
                        OK = 1;
                    }
                    if  (cmd->ch_from1[ct_l] == 'Q') {
                        cmd->S_from1[ct_l] = A[inarea].Q;
                        OK = 1;
                    }
                    if  (cmd->ch_from1[ct_l] == 'R') {
                        cmd->S_from1[ct_l] = A[inarea].R;
                        OK = 1;
                    }
                    if  (cmd->ch_from1[ct_l] == 'S') {
                        cmd->S_from1[ct_l] = A[inarea].S;
                        OK = 1;
                    }
                    if  (cmd->ch_from1[ct_l] == 'T') {
                        cmd->S_from1[ct_l] = A[inarea].T;
                        OK = 1;
                    }
                    if  (cmd->ch_from1[ct_l] == 'U') {
                        cmd->S_from1[ct_l] = A[inarea].U;
                        OK = 1;
                    }
                    if  (cmd->ch_from1[ct_l] == 'V') {
                        cmd->S_from1[ct_l] = A[inarea].V;
                        OK = 1;
                    }
                    if  (cmd->ch_from1[ct_l] == 'W') {
                        cmd->S_from1[ct_l] = A[inarea].W;
                        OK = 1;
                    }
                    if  (cmd->ch_from1[ct_l] == 'X') {
                        cmd->S_from1[ct_l] = A[inarea].X;
                        OK = 1;
                    }
                    if  (cmd->ch_from1[ct_l] == 'Y') {
                        cmd->S_from1[ct_l] = A[inarea].Y;
                        OK = 1;
                    }
                    if  (cmd->ch_from1[ct_l] == 'Z') {
                        cmd->S_from1[ct_l] = A[inarea].Z;
                        OK = 1;
                    }

                    if  (!OK)
                        fprintf(stderr,"\ncouldnt assign ch_from[%d]!\n", ct_l);
                }

                /**sources 2  ...  as above**/
                for (ct_l = 0; ct_l < cmd->anz_from2; ++ct_l) {

                    inarea = cmd->n_from2[ct_l];

                    OK = 0;

                    if  (cmd->ch_from2[ct_l] == 'A') {
                        cmd->S_from2[ct_l] = A[inarea].A;
                        OK = 1;
                    }
                    if  (cmd->ch_from2[ct_l] == 'B') {
                        cmd->S_from2[ct_l] = A[inarea].B;
                        OK = 1;
                    }
                    if  (cmd->ch_from2[ct_l] == 'C') {
                        cmd->S_from2[ct_l] = A[inarea].C;
                        OK = 1;
                    }
                    if  (cmd->ch_from2[ct_l] == 'D') {
                        cmd->S_from2[ct_l] = A[inarea].D;
                        OK = 1;
                    }
                    if  (cmd->ch_from2[ct_l] == 'E') {
                        cmd->S_from2[ct_l] = A[inarea].E;
                        OK = 1;
                    }
                    if  (cmd->ch_from2[ct_l] == 'F') {
                        cmd->S_from2[ct_l] = A[inarea].F;
                        OK = 1;
                    }
                    if  (cmd->ch_from2[ct_l] == 'G') {
                        cmd->S_from2[ct_l] = A[inarea].G;
                        OK = 1;
                    }
                    if  (cmd->ch_from2[ct_l] == 'H') {
                        cmd->S_from2[ct_l] = A[inarea].H;
                        OK = 1;
                    }
                    if  (cmd->ch_from2[ct_l] == 'I') {
                        cmd->S_from2[ct_l] = A[inarea].I;
                        OK = 1;
                    }
                    if  (cmd->ch_from2[ct_l] == 'J') {
                        cmd->S_from2[ct_l] = A[inarea].J;
                        OK = 1;
                    }
                    if  (cmd->ch_from2[ct_l] == 'K') {
                        cmd->S_from2[ct_l] = A[inarea].K;
                        OK = 1;
                    }
                    if  (cmd->ch_from2[ct_l] == 'L') {
                        cmd->S_from2[ct_l] = A[inarea].L;
                        OK = 1;
                    }
                    if  (cmd->ch_from2[ct_l] == 'M') {
                        cmd->S_from2[ct_l] = A[inarea].M;
                        OK = 1;
                    }
                    if  (cmd->ch_from2[ct_l] == 'N') {
                        cmd->S_from2[ct_l] = A[inarea].N;
                        OK = 1;
                    }
                    if  (cmd->ch_from2[ct_l] == 'O') {
                        cmd->S_from2[ct_l] = A[inarea].O;
                        OK = 1;
                    }
                    if  (cmd->ch_from2[ct_l] == 'P') {
                        cmd->S_from2[ct_l] = A[inarea].P;
                        OK = 1;
                    }
                    if  (cmd->ch_from2[ct_l] == 'Q') {
                        cmd->S_from2[ct_l] = A[inarea].Q;
                        OK = 1;
                    }
                    if  (cmd->ch_from2[ct_l] == 'R') {
                        cmd->S_from2[ct_l] = A[inarea].R;
                        OK = 1;
                    }
                    if  (cmd->ch_from2[ct_l] == 'S') {
                        cmd->S_from2[ct_l] = A[inarea].S;
                        OK = 1;
                    }
                    if  (cmd->ch_from2[ct_l] == 'T') {
                        cmd->S_from2[ct_l] = A[inarea].T;
                        OK = 1;
                    }
                    if  (cmd->ch_from2[ct_l] == 'U') {
                        cmd->S_from2[ct_l] = A[inarea].U;
                        OK = 1;
                    }
                    if  (cmd->ch_from2[ct_l] == 'V') {
                        cmd->S_from2[ct_l] = A[inarea].V;
                        OK = 1;
                    }
                    if  (cmd->ch_from2[ct_l] == 'W') {
                        cmd->S_from2[ct_l] = A[inarea].W;
                        OK = 1;
                    }
                    if  (cmd->ch_from2[ct_l] == 'X') {
                        cmd->S_from2[ct_l] = A[inarea].X;
                        OK = 1;
                    }
                    if  (cmd->ch_from2[ct_l] == 'Y') {
                        cmd->S_from2[ct_l] = A[inarea].Y;
                        OK = 1;
                    }
                    if  (cmd->ch_from2[ct_l] == 'Z') {
                        cmd->S_from2[ct_l] = A[inarea].Z;
                        OK = 1;
                    }

                    if  (!OK)
                        fprintf(stderr,"\ncouldnt assign ch_from2[%d]!\n",ct_l);
                }

                /**sources 3  ...  as above**/
                for (ct_l = 0; ct_l < cmd->anz_from3; ++ct_l) {

                    inarea = cmd->n_from3[ct_l];

                    OK = 0;

                    if  (cmd->ch_from3[ct_l] == 'A') {
                        cmd->S_from3[ct_l] = A[inarea].A;
                        OK = 1;
                    }
                    if  (cmd->ch_from3[ct_l] == 'B') {
                        cmd->S_from3[ct_l] = A[inarea].B;
                        OK = 1;
                    }
                    if  (cmd->ch_from3[ct_l] == 'C') {
                        cmd->S_from3[ct_l] = A[inarea].C;
                        OK = 1;
                    }
                    if  (cmd->ch_from3[ct_l] == 'D') {
                        cmd->S_from3[ct_l] = A[inarea].D;
                        OK = 1;
                    }
                    if  (cmd->ch_from3[ct_l] == 'E') {
                        cmd->S_from3[ct_l] = A[inarea].E;
                        OK = 1;
                    }
                    if  (cmd->ch_from3[ct_l] == 'F') {
                        cmd->S_from3[ct_l] = A[inarea].F;
                        OK = 1;
                    }
                    if  (cmd->ch_from3[ct_l] == 'G') {
                        cmd->S_from3[ct_l] = A[inarea].G;
                        OK = 1;
                    }
                    if  (cmd->ch_from3[ct_l] == 'H') {
                        cmd->S_from3[ct_l] = A[inarea].H;
                        OK = 1;
                    }
                    if  (cmd->ch_from3[ct_l] == 'I') {
                        cmd->S_from3[ct_l] = A[inarea].I;
                        OK = 1;
                    }
                    if  (cmd->ch_from3[ct_l] == 'J') {
                        cmd->S_from3[ct_l] = A[inarea].J;
                        OK = 1;
                    }
                    if  (cmd->ch_from3[ct_l] == 'K') {
                        cmd->S_from3[ct_l] = A[inarea].K;
                        OK = 1;
                    }
                    if  (cmd->ch_from3[ct_l] == 'L') {
                        cmd->S_from3[ct_l] = A[inarea].L;
                        OK = 1;
                    }
                    if  (cmd->ch_from3[ct_l] == 'M') {
                        cmd->S_from3[ct_l] = A[inarea].M;
                        OK = 1;
                    }
                    if  (cmd->ch_from3[ct_l] == 'N') {
                        cmd->S_from3[ct_l] = A[inarea].N;
                        OK = 1;
                    }
                    if  (cmd->ch_from3[ct_l] == 'O') {
                        cmd->S_from3[ct_l] = A[inarea].O;
                        OK = 1;
                    }
                    if  (cmd->ch_from3[ct_l] == 'P') {
                        cmd->S_from3[ct_l] = A[inarea].P;
                        OK = 1;
                    }
                    if  (cmd->ch_from3[ct_l] == 'Q') {
                        cmd->S_from3[ct_l] = A[inarea].Q;
                        OK = 1;
                    }
                    if  (cmd->ch_from3[ct_l] == 'R') {
                        cmd->S_from3[ct_l] = A[inarea].R;
                        OK = 1;
                    }
                    if  (cmd->ch_from3[ct_l] == 'S') {
                        cmd->S_from3[ct_l] = A[inarea].S;
                        OK = 1;
                    }
                    if  (cmd->ch_from3[ct_l] == 'T') {
                        cmd->S_from3[ct_l] = A[inarea].T;
                        OK = 1;
                    }
                    if  (cmd->ch_from3[ct_l] == 'U') {
                        cmd->S_from3[ct_l] = A[inarea].U;
                        OK = 1;
                    }
                    if  (cmd->ch_from3[ct_l] == 'V') {
                        cmd->S_from3[ct_l] = A[inarea].V;
                        OK = 1;
                    }
                    if  (cmd->ch_from3[ct_l] == 'W') {
                        cmd->S_from3[ct_l] = A[inarea].W;
                        OK = 1;
                    }
                    if  (cmd->ch_from3[ct_l] == 'X') {
                        cmd->S_from3[ct_l] = A[inarea].X;
                        OK = 1;
                    }
                    if  (cmd->ch_from3[ct_l] == 'Y') {
                        cmd->S_from3[ct_l] = A[inarea].Y;
                        OK = 1;
                    }
                    if  (cmd->ch_from3[ct_l] == 'Z') {
                        cmd->S_from3[ct_l] = A[inarea].Z;
                        OK = 1;
                    }

                    if  (!OK)
                        fprintf(stderr,"\ncouldnt assign ch_from3[%d]!\n",ct_l);
                }

            } /**commands**/
        } /**stays**/
      } /**sweeps**/
    } /**series**/
}


/*************************** check_local *************************************/

void check_local (COMMAND *cmd, char *name, int anz_par) {

    if  (!strcmp (cmd->func_name, name)) {

        if  (cmd->anz_quantums != 1)
            fprintf (stderr, "\ncheck_local: %d instead 1 quantums in %s\n",
                             cmd->anz_quantums, cmd->func_name);

        if  (anz_par == 0) {
            if  ((cmd->anz_quant[0] != 1) || (cmd->quantum[0][0] != 0.0))
                fprintf (stderr, "\ncheck_local: err with quantums of %s\n",
                                 cmd->func_name);
        } else {
            if  (cmd->anz_quant[0] != anz_par)
                fprintf (stderr, "\ncheck_local prefers %d to %d quant's in %s\n",
                                 anz_par, cmd->anz_quant[0], cmd->func_name);
        }
    }
}


/*************************** check_local_quantums ****************************/

void check_local_quantums (COMMAND *cmd) {

    check_local (cmd, "local_sum"          , 0);
    check_local (cmd, "local_sub"          , 0);
    check_local (cmd, "local_mult"         , 0);
    check_local (cmd, "local_divide"       , 0);
    check_local (cmd, "local_average"      , 1);
    check_local (cmd, "local_copy"         , 0);
    check_local (cmd, "local_persist"      , 0);
    check_local (cmd, "local_const"        , 1);
    check_local (cmd, "local_rectify"      , 2);
    check_local (cmd, "local_invert "      , 0);
    check_local (cmd, "local_abs"          , 0);
    check_local (cmd, "local_sqrt"         , 0);
    check_local (cmd, "local_power"        , 1);
    check_local (cmd, "local_lin_01"       , 0);
    check_local (cmd, "local_rand"         , 2);
    check_local (cmd, "local_rand_pos_neg" , 0);
    check_local (cmd, "local_as_rand"      , 2);
    check_local (cmd, "local_gauss"        , 1);
    check_local (cmd, "local_gauss_pos"    , 1);
    check_local (cmd, "local_poiss"        , 1);
    check_local (cmd, "local_rand_gibbs"   , 2);
    check_local (cmd, "local_rand_gibbs_01", 2);
    check_local (cmd, "local_tanh"         , 0);
    check_local (cmd, "local_bp_tanh"      , 0);
    check_local (cmd, "local_sparse"       , 1);
    check_local (cmd, "local_sparse_diff"  , 1);
    check_local (cmd, "local_sparse_01"    , 1);
    check_local (cmd, "local_zhang"        , 4);
    check_local (cmd, "local_zhang2"       , 3);
    check_local (cmd, "local_zhang_inv"    , 4);
    check_local (cmd, "local_zhang2_inv"   , 3);
    check_local (cmd, "local_zhang_diff"   , 4);
    check_local (cmd, "local_zhang2_diff"  , 3);
    check_local (cmd, "local_mean"         , 3);
    check_local (cmd, "local_mean_diff"    , 2);
    check_local (cmd, "local_mean_01"      , 5);
    check_local (cmd, "local_bp_mean_01"   , 0);
    check_local (cmd, "local_mean_01_inv"  , 5);
    check_local (cmd, "local_mean_01_diff" , 5);
    check_local (cmd, "local_mean_01_d0"   , 5);
    check_local (cmd, "local_mean_01_d1"   , 5);
    check_local (cmd, "local_mean_01_d2"   , 5);
    check_local (cmd, "local_mean_01_d3"   , 5);
    check_local (cmd, "local_mean_01_d4"   , 5);
    check_local (cmd, "local_gibbs"        , 3);
    check_local (cmd, "local_gibbs_01"     , 5);
    check_local (cmd, "local_gibbs_01_2"   , 2);
    check_local (cmd, "local_odds_01"      , 1);
}


void check_vehicle (SIMULATION *si) {

    int ct_se, ct_sw, ct_st, ct_cmd;
    SERIES  *se;
    COMMAND *cmd;

    for (ct_se = 0; ct_se < si->anz_se; ++ct_se) {
        se = si->se + ct_se;
        /**sweeps**/
        for (ct_sw = 0; ct_sw < se->anz_sw; ++ct_sw) {
            /**stays**/
            for (ct_st = 0; ct_st < se->sw[ct_sw].anz_st; ++ct_st) {
                /**commands**/
                for (ct_cmd=0; ct_cmd < se->sw[ct_sw].st[ct_st].anz_cmd; ++ct_cmd) {

                    cmd = se->sw[ct_sw].st[ct_st].cmd + ct_cmd;

                    /**ch_target fits to func_key?**/
                    if  ((cmd->ch_target == 'W') || (cmd->ch_target == 'V') || (cmd->ch_target == 'X') || (cmd->ch_target == 'Y'))
                        if  ((cmd->func_key != 'w') && (cmd->func_key != 'o') && (cmd->func_key != 't'))
                            fprintf(stderr, "\nfunc_key does'nt fit to target!\n");

                    /**local function has 2 proper n_from, n_from_2?**/
                    if  (cmd->func_key == 'l') {
                        if  (cmd->anz_arguments != 2)
                            fprintf (stderr, "\nlocal func must have 2 args!\n");
                        if  ( (cmd->area != cmd->n_from1[0])
                           || (cmd->area != cmd->n_from2[0])
                           || (cmd->area != cmd->n_from3[0]))
                            fprintf (stderr, "\nlocal func must be local!\n");
                        check_local_quantums (cmd);
                    }
                }
            }
        }
    }
}



/******************************* decr_eps ************************************
 * Decrease all Z[i].eps_W,V[j] linearly to zero after x->elen / 2.          *
void decr_eps (PARAMS *x, AGENT *Z, int epoche, int elen) {
    int i, j;
    static int firsttime = 1;
    static double **save_eps_W;
    if  (firsttime) {
        save_eps_W = d_matrix (x->areas, x->areas);
        for (i = 0; i < x->areas; ++i)
            for (j = 0; j < x->areas; ++j) {
                save_eps_W[i][j] = Z[i].eps_W[j];
            }
        firsttime = 0;
    }
    if  (epoche > elen * 3 / 4) {

        for (i = 0; i < x->areas; ++i)
            for (j = 0; j < x->areas; ++j) {

                Z[i].eps_W[j] = save_eps_W[i][j] * (4.0 - 4.0*epoche/elen);

                if  (epoche % 25 == 0) {
                    if  (Z[i].eps_W[j] != 0.0)
                        fprintf (stderr, " eps_W[%d][%d]=%f ",
                                                  i, j, Z[i].eps_W[j]);
                }
            }
    }
    **reset values for a second, ..., series**
    if  (epoche == 0) {

        fprintf (stderr, "\neps are: ");

        for (i = 0; i < x->areas; ++i)
            for (j = 0; j < x->areas; ++j) {

                Z[i].eps_W[j] = save_eps_W[i][j];

                if  (Z[i].eps_W[j] != 0.0)
                    fprintf (stderr, "eps_W[%d][%d]=%f ", i, j, Z[i].eps_W[j]);
            }

        fprintf (stderr, "\n");
    }
}
**/


/******************************* print_command *******************************/
/* Print contents of a command for debugging only.                           */

void print_command (COMMAND *cmd) {

    int ct_l, ii;

    fprintf (stderr, "\ncmd->area=%d, ", cmd->area);
    fprintf (stderr, "moment=%f, ", cmd->moment);
    fprintf (stderr, "func_name=%s, ", cmd->func_name);
    fprintf (stderr, "ch_target=%c, ", cmd->ch_target);
    fprintf (stderr, "func_key=%c, ", cmd->func_key);
    fprintf (stderr, "anz_arguments=%d, ", cmd->anz_arguments);

    fprintf (stderr, "anz_from1=%d, ", cmd->anz_from1);
    fprintf (stderr, "anz_from2=%d, ", cmd->anz_from2);
    fprintf (stderr, "anz_from3=%d, ", cmd->anz_from3);

    fprintf (stderr, "n_from1= ");
    for (ct_l = 0; ct_l < cmd->anz_from1; ++ct_l)
        fprintf (stderr, "%d, ", cmd->n_from1[ct_l]);
    fprintf (stderr, "n_from2= ");
    for (ct_l = 0; ct_l < cmd->anz_from2; ++ct_l)
        fprintf (stderr, "%d, ", cmd->n_from2[ct_l]);
    fprintf (stderr, "n_from3= ");
    for (ct_l = 0; ct_l < cmd->anz_from3; ++ct_l)
        fprintf (stderr, "%d, ", cmd->n_from3[ct_l]);

    fprintf (stderr, "ch_from1= ");
    for (ct_l = 0; ct_l < cmd->anz_from1; ++ct_l)
        fprintf (stderr, "%c, ", cmd->ch_from1[ct_l]);
    fprintf (stderr, "ch_from2= ");
    for (ct_l = 0; ct_l < cmd->anz_from2; ++ct_l)
        fprintf (stderr, "%c, ", cmd->ch_from2[ct_l]);
    fprintf (stderr, "ch_from3= ");
    for (ct_l = 0; ct_l < cmd->anz_from3; ++ct_l)
        fprintf (stderr, "%c, ", cmd->ch_from3[ct_l]);

     for (ct_l = 0; ct_l < cmd->anz_quantums; ++ct_l) {
         if  (ct_l) fprintf (stderr, ", ");
         for (ii = 0; ii < cmd->anz_quant[ct_l]; ++ii) {
             if  (ii)  fprintf (stderr, "+");
             if  (cmd->quantum[ct_l][ii] != 0.0)
                 floatprint (stderr, cmd->quantum[ct_l][ii]);
         }
     }

    for (ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {
        if  ((cmd->ch_from1[ct_l] == 'K') || (cmd->ch_from1[ct_l] == 'L') ||
             (cmd->ch_from1[ct_l] == 'M') || (cmd->ch_from1[ct_l] == 'N') ||
             (cmd->ch_from1[ct_l] == 'O') || (cmd->ch_from1[ct_l] == 'P') ||
             (cmd->ch_from1[ct_l] == 'Q') || (cmd->ch_from1[ct_l] == 'R') ||
             (cmd->ch_from1[ct_l] == 'S') || (cmd->ch_from1[ct_l] == 'T')) {
            fprintf (stderr, " S_from1[%d][0][0]=", ct_l);
            fprintf (stderr, "%f, ", cmd->S_from1[ct_l][0][0]);
        }
    }

    for (ct_l = 0; ct_l < cmd->anz_from2; ++ct_l) {
        if  ((cmd->ch_from2[ct_l] == 'K') || (cmd->ch_from2[ct_l] == 'L') ||
             (cmd->ch_from2[ct_l] == 'M') || (cmd->ch_from2[ct_l] == 'N') ||
             (cmd->ch_from2[ct_l] == 'O') || (cmd->ch_from2[ct_l] == 'P') ||
             (cmd->ch_from2[ct_l] == 'Q') || (cmd->ch_from2[ct_l] == 'R') ||
             (cmd->ch_from2[ct_l] == 'S') || (cmd->ch_from2[ct_l] == 'T')) {
            fprintf (stderr, " S_from2[%d][0][0]=", ct_l);
            fprintf (stderr, "%f, ", cmd->S_from2[ct_l][0][0]);
        }
    }
	fprintf (stderr, "\n");

}
