#ifndef _feed_H
#define _feed_H

double feed_in_W         (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
double feed_in_V         (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
double feed_in_X         (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
double feed_in_Y         (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
double feed_in_V_pm      (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
double feed_cheap        (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
double feed_moderate     (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
double feed_full         (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
double feed_l_replace    (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
double feed_back_W       (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
double feed_euclid_W     (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
double feed_l_add        (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
double feed_copy         (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
double feed_copy_limited (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
double feed_l_theta      (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
double feed_l_difu       (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
double feed_l_rand_from  (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
double feed_l_covar      (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
double feed_l_hack       (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
double feed_const_hack   (AREA *A, COMMAND *cmd, int ct_t, int ct_n);

double feed_mean_01      (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
double feed_mean_01_inv  (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
double feed_mean_01_diff (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
double feed_mean_01_d0   (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
double feed_mean_01_d1   (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
double feed_mean_01_d2   (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
double feed_mean_01_d3   (AREA *A, COMMAND *cmd, int ct_t, int ct_n);
double feed_mean_01_d4   (AREA *A, COMMAND *cmd, int ct_t, int ct_n);

#endif
