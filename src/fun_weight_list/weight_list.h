#ifndef _weight_list_H
#define _weight_list_H


typedef struct CONN {

 DOUBLE val;            /**connection strength**/

 int    a, b;           /**index of prae-synaptic cell**/
 int    n;              /**n = a * d_b + b**/

 CONN   * next;         /**next connection**/

} CONN;


typedef struct ALL_WEIGHTS {

 int areas;             /**number of areas = x->areas**/
 int *d_a, *d_b, *d_n;  /**area sizes, d_n = d_a * d_b**/

 CONN ****all_conn;     /**all_conn[area][ct_n][inarea](is a conn * to the (first->admin) connection)**/

} ALL_WEIGHTS;


const int MAXCOL_INT  = 999999;
const int MAXCOL_CHAR = 127;


DOUBLE weight_list_alloc_full         (PARAMS *g, AREA *A, COMMAND *cmd, int, int);
DOUBLE weight_list_alloc_topo         (PARAMS *g, AREA *A, COMMAND *cmd, int, int);
DOUBLE weight_list_alloc_invert       (PARAMS *g, AREA *A, COMMAND *cmd, int, int);
DOUBLE weight_list_invert_quick       (PARAMS *g, AREA *A, COMMAND *cmd, int, int);
DOUBLE weight_list_mult               (PARAMS *g, AREA *A, COMMAND *cmd, int, int);
DOUBLE weight_list_free               (PARAMS *g, AREA *A, COMMAND *cmd, int, int);
DOUBLE weight_list_feed               (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE WEIGHT_LIST_PUSH               (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_list_AXON_total         (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy);
DOUBLE weight_list_AXON_total_pos     (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy);
DOUBLE weight_list_hebb               (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_list_hebb_diff          (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_list_ica1               (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy);
DOUBLE weight_list_ica2               (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy);
DOUBLE weight_list_kohonen            (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_list_euclid             (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_list_euclid_ptr         (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_list_act2weight         (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_list_act2weight_offset  (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_list_hebb_turn          (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_list_decay              (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_list_decay_const        (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_list_decay_quad         (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_list_decay_post         (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_list_decay_pre          (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_list_decay_Baddeley     (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_list_decay_Baddeley_distance  (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_list_normalize          (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_list_rectify            (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_list_cutself            (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_list_cutsmall           (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_list_sprout             (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_list_export             (PARAMS *g, AREA *A, COMMAND *cmd, int, int);
DOUBLE weight_list_export_col         (PARAMS *g, AREA *A, COMMAND *cmd, int, int);
DOUBLE weight_list_export_gnu         (PARAMS *g, AREA *A, COMMAND *cmd, int, int);
DOUBLE weight_list_alloc_import       (PARAMS *g, AREA *A, COMMAND *cmd, int, int);
DOUBLE weight_list_alloc_import_tiles (PARAMS *g, AREA *A, COMMAND *cmd, int, int);
DOUBLE weight_list_histogram          (PARAMS *g, AREA *A, COMMAND *cmd, int, int);
DOUBLE observe_phase                  (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE import_phase                   (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE weight_list_cuthalfinput       (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);

void rgbmatrix2file (COMMAND *cmd, DOUBLE ***rgbmatrix, int d_a, int d_b, int d_a_in, int d_b_in, DOUBLE min, DOUBLE max, int format, int key);

#endif
