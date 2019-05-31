#ifndef _vehicle_H
#define _vehicle_H

void assign_conditions (SIMULATION *si);
void printparams (PARAMS *x, AREA *A, SIMULATION *si, FILE *fp);
void free_se (SIMULATION *si);
AREA *alloc_a (PARAMS *x, SIMULATION *si);
void free_a (AREA *A, PARAMS *x, SIMULATION *si);
void init_a (AREA *A, PARAMS *x, SIMULATION *si);
void choose_functions (PARAMS *x, SIMULATION *si);
void choose_pointers (SIMULATION *si, AREA *A);
void check_vehicle (SIMULATION *si);
void decr_eps (PARAMS *x, int epoche, int elen);
void print_command (COMMAND *cmd);                          /**for debugging**/

#endif
