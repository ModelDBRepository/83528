#ifndef _weight_sigpi_H
#define _weight_sigpi_H


typedef struct CONN_SIGPI {

 DOUBLE val;            /**connection strength**/

 int    a1, b1;         /**index of first prae-synaptic cell      **/
 int    n1;             /**n = a * d_b + b; for faster computation**/
 int    a2, b2;         /**    ~    second       ~                **/
 int    n2;             /**                                       **/

 CONN_SIGPI   * next;         /**next connection**/

} CONN_SIGPI;


typedef struct ALL_WEIGHTS_SIGPI {

 int areas;                    /**number of areas = x->areas**/
 int *d_a, *d_b, *d_n;         /**area sizes, d_n = d_a * d_b**/

 CONN_SIGPI *****all_conn;     /**all_conn[area][ct_n][inarea1][inarea2](is a conn * to the (first->admin) connection)**/

} ALL_WEIGHTS_SIGPI;



DOUBLE weight_sigpi_alloc_full     (PARAMS *g, AREA *A, COMMAND *cmd, int, int);
DOUBLE weight_sigpi_feed           (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_sigpi_euclid         (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_sigpi_kohonen        (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_sigpi_histogram      (PARAMS *g, AREA *A, COMMAND *cmd, int, int);

DOUBLE weight_sigpi_hebb           (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);

/*
DOUBLE weight_sigpi_alloc_topo     (PARAMS *g, AREA *A, COMMAND *cmd, int, int);
DOUBLE weight_sigpi_alloc_invert   (PARAMS *g, AREA *A, COMMAND *cmd, int, int);
DOUBLE WEIGHT_SIGPI_PUSH           (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_sigpi_AXON_total     (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy);
DOUBLE weight_sigpi_AXON_total_pos (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy);
DOUBLE weight_sigpi_hebb_turn      (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_sigpi_decay          (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_sigpi_decay_quad     (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_sigpi_decay_post     (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_sigpi_decay_pre      (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_sigpi_normalize      (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_sigpi_rectify        (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_sigpi_cutself        (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
*/
DOUBLE weight_sigpi_cutsmall       (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
/*
DOUBLE weight_sigpi_sprout         (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n);
DOUBLE weight_sigpi_export         (PARAMS *g, AREA *A, COMMAND *cmd, int, int);
DOUBLE weight_sigpi_export_col     (PARAMS *g, AREA *A, COMMAND *cmd, int, int);
DOUBLE weight_sigpi_export_gnu     (PARAMS *g, AREA *A, COMMAND *cmd, int, int);
*/

DOUBLE weight_sigpi_export         (PARAMS *g, AREA *A, COMMAND *cmd, int, int);
DOUBLE weight_sigpi_alloc_import   (PARAMS *g, AREA *A, COMMAND *cmd, int, int);
DOUBLE weight_sigpi_gnuplot_cm     (PARAMS *g, AREA *A, COMMAND *cmd, int, int);

#endif
