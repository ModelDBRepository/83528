#ifndef _series_H
#define _series_H

typedef struct COMMAND {      /**instructions to update sth of an area   cmd**/

  int     area;           /**handed down from stay**/
  double  moment;         /**handed down from sweep**/

/*1*/
  char    ch_target;      /**char of what to compute**/
  DOUBLE  **S_target;    /**points to A[area].R/S/T**/

  int     anz_pointers;   /**number of pointers in the cmd, separated by "+"**/
  char    **ch_pointers;  /**char list of pointers**/
  GLOBPTR **pointers;


/*2*/
  char    *func_name;    /**name of function**/

  DOUBLE  (*localfunc)  (DOUBLE *par, DOUBLE val1, DOUBLE val2);                                /**simple**/
  DOUBLE  (*func)       (PARAMS *x, AREA *A, struct COMMAND *cmd, int ival1, int ival2);        /**ival1,2 = (feedfunc: ct_t,ct_n); (totalfunc: ct_t,dummy); (alltimefunc: begin,end)**/
  char    func_key;                                                                             /**see definitions in vehicle.c**/

/*3*/
  int     anz_arguments; /*check with function desire*/

  int     anz_from1;     /**input areas separated by "+"**/
  int     *n_from1;      /**list of input areas**/
  int     anz_from2;     /**for 2nd_cheap and local fkts (behind the ",")**/
  int     *n_from2;      /**list of modulating areas**/
  int     anz_from3;     /**for sigpi weight fkts (behind the second ",")**/
  int     *n_from3;      /**list of modulating areas**/

/*4*/
  char    *ch_from1;     /**char list of first input sources**/
  char    *ch_from2;     /**char list of second input sources**/
  char    *ch_from3;     /**char list of third input sources**/

  DOUBLE  ***S_from1;    /**pointers to A[inarea].R/T (first index like:ct_l) -> [ct_l][mult][ct_t][ct_n]**/
  DOUBLE  ***S_from2;    /**   "  ; inarea derived from first index**/
  DOUBLE  ***S_from3;    /**   "  ; inarea derived from first index**/

/*5*/
  int     anz_quantums;
  int     *anz_quant;
  DOUBLE  **quantum;

} COMMAND;


typedef struct CONDITION {

  int       op_tag; /**0=no test; 1=only left arg; 2="="; 3="%"; 4="<"; 5=">"**/
  int       *left;
  int       *right;
  int       (*test)(int *left, int *right);

} CONDITION;


typedef struct STAY {            /**stack of anz_cmd commands at a stay   st**/

  int       area;          /**where to stay**/

  char      st_update;     /**n(one) / o(rder) / r(andom) / t(otal)  no shuffle**/

  int       anz_cmd;
  COMMAND   *cmd;

} STAY;


typedef struct SWEEP {                  /**set of anz_st stays to relax   sw**/

  int       begin;         /**relaxation begin**/
  int       end;           /**relaxation end**/
  char      *sw_update;    /**string: order/random/propto / alltime**/

  int       anz_st;
  CONDITION **cond;
  STAY      *st;

} SWEEP;


typedef struct SERIES {                      /**pointer to anz_sw sweeps  se**/

  int       ilen;                                    /**number of iterations**/

  int       anz_sw;
  CONDITION **cond;
  SWEEP     *sw;

} SERIES;


typedef struct SIMULATION {                 /**pointer to anz_se series   si**/

  int       max_anz_se;

  int       anz_se;
  CONDITION **cond;
  SERIES    *se;

  int       num_all_cond;
  CONDITION **all_cond;

} SIMULATION;


#endif
