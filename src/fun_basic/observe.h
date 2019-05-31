#ifndef _observe_H
#define _observe_H

DOUBLE observe_act          (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE observe_col          (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE observe_animgif      (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE observe_animgif_iter (PARAMS *g, AREA *A, COMMAND *cmd, int, int);
DOUBLE observe_getc         (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE observe_gnu_act      (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE observe_gnu_two_acts (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE observe_act_hist     (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE observe_phase_hist   (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);

/*
void exportP36_matrix (double **S, int d_x, int d_y, int d_a, int d_b, int format, char datei[256]);
void importP36_matrix (double **S, int d_x, int d_y, int d_a, int d_b, char datei[256]);
void export_weights (PARAMS *g, AREA *A, char dir[256]);
void export_Theta   (PARAMS *g, AREA *A, char dir[256]);
void observe_onoff_weights (AREA *A, COMMAND *cmd, int begin, int end, int ilen);
void observe_col      (AREA *A, COMMAND *cmd, int begin, int end, int ilen);
void observe_act_gnu  (AREA *A, COMMAND *cmd, int begin, int end, int ilen);
void observe_quad     (AREA *A, COMMAND *cmd, int begin, int end, int ilen);
void observe_quad2    (AREA *A, COMMAND *cmd, int begin, int end, int ilen);
void observe_quad3    (AREA *A, COMMAND *cmd, int begin, int end, int ilen);
void observe_quad4    (AREA *A, COMMAND *cmd, int begin, int end, int ilen);
void observe_act_stat (AREA *A, COMMAND *cmd, int begin, int end, int ilen);
void observe_act_time (AREA *A, COMMAND *cmd, int begin, int end, int ilen);
void observe_sophie1d (AREA *A, COMMAND *cmd, int begin, int end, int ilen);
*/

#endif
