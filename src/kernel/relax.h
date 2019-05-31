#ifndef _relax_H
#define _relax_H

void do_simulation (PARAMS *x, AREA *A, SIMULATION *si);
void do_series     (PARAMS *x, AREA *A, SERIES *se);              /** needed directly by NNsimObject.cpp, but NOT by iter.c**/
void do_epoche     (PARAMS *x, AREA *A, SERIES *se);              /**dosn't exist any more!!!**/

#endif
