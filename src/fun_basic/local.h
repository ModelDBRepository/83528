#ifndef _local_H
#define _local_H

DOUBLE local_sum             (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_sum_const       (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_sub             (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_mult            (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_mult_const      (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_div             (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_div_const       (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_modulo_const    (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_cyclic          (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_average         (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_copy            (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_persist         (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_threshold       (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_isbigger        (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_biggerof        (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_biggerfabsof    (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_const           (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_rectify         (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_invert          (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_abs             (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_sign            (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_round           (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_sqrt            (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_power           (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_exp             (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_log_pos         (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_cos             (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_sin             (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_tan             (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_atan            (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_atan_to_2pi     (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_lin_01          (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_rand            (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_rand_pos_neg    (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_as_rand         (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_gauss           (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_circ_gauss      (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_gauss_pos       (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_rand_gibbs      (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_rand_gibbs_01   (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_rand_exp        (DOUBLE *par, DOUBLE val1, DOUBLE val2);

DOUBLE local_tanh            (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_bp_tanh         (DOUBLE *par, DOUBLE val1, DOUBLE val2);

DOUBLE local_sparse          (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_sparse_diff     (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_sparse_01       (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_zhang           (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_zhang_scale     (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_zhang2          (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_zhang_inv       (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_zhang2_inv      (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_zhang_diff      (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_zhang2_diff     (DOUBLE *par, DOUBLE val1, DOUBLE val2);

DOUBLE local_mean            (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_mean_diff       (DOUBLE *par, DOUBLE val1, DOUBLE val2);

DOUBLE local_mean_01         (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_mean_scale      (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_mean_scale_diff (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_bp_mean_01      (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_mean_01_inv     (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_mean_01_diff    (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_mean_01_d0      (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_mean_01_d1      (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_mean_01_d2      (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_mean_01_d3      (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_mean_01_d4      (DOUBLE *par, DOUBLE val1, DOUBLE val2);

DOUBLE local_gibbs           (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_gibbs_01        (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_gibbs_01_2      (DOUBLE *par, DOUBLE val1, DOUBLE val2);
DOUBLE local_odds_01         (DOUBLE *par, DOUBLE val1, DOUBLE val2);

#endif
