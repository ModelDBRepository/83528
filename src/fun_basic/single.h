#ifndef _single_H
#define _single_H

DOUBLE single_ptr2act      (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE set_ptr_int_val     (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE set_ptr_val         (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE single_copy         (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE single_copy_limited (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE single_mean_back    (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE single_circ_gauss   (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE single_l_add        (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);

#endif
