#ifndef _weight_H
#define _weight_H

void weight_hebb        (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
void weight_antihebb    (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
void weight_kohonen     (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
void weight_heuristic   (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
void weight_log         (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
void weight_ica         (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
void weight_rec_mean_01 (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
void weight_copy_topo   (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
void weight_add_diag    (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
void weight_self_zero   (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
void dweight_self_set   (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
void weight_seung       (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
void weight_respcorr    (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
void weight_antirespcorr
                        (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
void weight_decay       (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
void weight_decay_topo  (AREA *A, COMMAND *cmd, int ct_t, int ct_n);

void weight_incr_glob   (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
void weight_incr_pos_glob
                        (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
void weight_act_decay   (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
void weight_norm_glob   (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
void weight_norm_each   (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
void weight_rectify     (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
void weight_rect_eps    (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
void weight_add_rand    (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
void weight_bound       (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
void weight_copy        (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
void weight_hack_init   (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
void weight_hack_slow   (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
void weight_hack_slow_inv
                        (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
void init_theta_const   (AREA *A, COMMAND *cmd, int ct_t, int ct_n);

/**used in sweep named weight**/
void weight_ica2       (AREA *A, COMMAND *cmd, int, int);
void weight_invert     (AREA *A, COMMAND *cmd, int, int);

/**used in iter.c**/
void weight_initdelta (PARAMS *x, AREA *A);
void weight_update    (PARAMS *x, AREA *A);

#endif
