#ifndef _reinforcement_H
#define _reinforcement_H

DOUBLE feed_reward          (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE feed_angl_reward     (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE feed_dock_reward     (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);

DOUBLE total_rand_winner    (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy);
DOUBLE total_init_place     (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy);
DOUBLE total_move_place     (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy);
DOUBLE total_init_coord     (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy);
DOUBLE total_scan_coord     (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy);
DOUBLE total_imag_coord     (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy);
DOUBLE total_coord_imag     (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy);
DOUBLE total_angl_coord     (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy);
DOUBLE total_angl_where     (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy);
DOUBLE total_xyp_to_dtp     (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy);
DOUBLE total_move_coord     (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy);

#endif
