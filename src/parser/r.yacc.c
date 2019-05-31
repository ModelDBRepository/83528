%{
   #include <stdlib.h>
   #include <stdio.h>
   #include <string.h>

   extern int yylex();
   extern int yyparse();

   void yyerror(const char *msg);
   int test_none (int *left, int *right);
   int test_one (int *left, int *right);
   int test_equal (int *left, int *right);
   int test_modulo (int *left, int *right);
   int test_smaller (int *left, int *right);
   int test_larger (int *left, int *right);

   #include "../kernel/coco.h"
   #include "../kernel/series.h"
   #include "../kernel/relax.h"

   static PARAMS     *g;
   static AREA       *A;
   static SIMULATION *si;

   int     i, ii, ct_ar;
   int     ct_c = 0;
   int     ct_p = 0;
   GLOBPTR *globptr;
   int     op_tag;

   int     ch_from1_missing = 0, ch_from2_missing = 0, ch_from3_missing = 0;
   int     anz_ch_from1 = 0, anz_ch_from2 = 0, anz_ch_from3 = 0;
   SERIES  *se;
   SWEEP   *curr_sw;
   COMMAND *cmd;
   char    names[128][128];
   char    name[128];
   double  numbers[128];

   COMMAND oldcmd;
%}

%union {
  char *string; /* string buffer */
  int  cmd;     /* command value */
}

%token <string> ID NUMBER
%token <cmd> GLOBAL ALL AREAS SER SWE IFWORD ITERWORD EQ MOD SM LG
%token <cmd> GI_ _IG RI_ _IR
%token <cmd> COLON SEP KOMMA QUOTE PLUS MULT

%type <string> id num

%start blocks

%%

blocks  : block
        | blocks block
        ;

block   : global GI_  glines  _IG
        | all    GI_  alines  _IG
        | area   GI_  nlines  _IG
        | series RI_ ilen _IR GI_ sweeps _IG
        ;

        /**********global**********/

global  : GLOBAL           { fprintf (stderr, "\n0. global  ");

                             g = (PARAMS *) malloc (sizeof (PARAMS));
                             g->num_pointers = 0;
                             g->ch_pointers  = (char **) malloc (100 * sizeof (char *));                          /**max num pointers = 100!**/
                             g->pointers     = (GLOBPTR **) malloc (100 * sizeof (void *));

                             si = (SIMULATION *) malloc (sizeof (SIMULATION));
                             si->max_anz_se  = 30;                                                                /**max num series = 100 -- see also assign_conditions in vehicle.c!!**/
                             si->se          = (SERIES *) malloc (si->max_anz_se * sizeof (SERIES));
                             si->anz_se      = 0;

                             si->all_cond        = (CONDITION **) malloc (10000 * sizeof (CONDITION *));          /**max num conditions = 10000 -- see also assign_conditions in vehicle.c!**/
                           }

glines  : gline            { g->num_pointers_at_global = g->num_pointers;
                           }
        | gline glines
        ;

gline   : id num           { /**reading GLObal scalar parameters -- such as iter, areas, mult**/
                             #undef CONVERT_GLO
                             #define CONVERT_GLO(a,b) \
                                     if  (!strcmp ($1, #b)) \
                                         g->b = (a)(atof($2));
                             #include "../kernel/parameters.h"
                             #undef CONVERT_GLO
                             #define CONVERT_GLO(a,b) ;
                           }
        | ITERWORD num     { /**this hack is needed for "iter" which is defined in r.lex, so it isn't accepted any more as id**/
                             #undef CONVERT_GLO
                             #define CONVERT_GLO(a,b) \
                                     if  (!strcmp ("iter", #b)) \
                                         g->b = (a)(atof($2));
                             #include "../kernel/parameters.h"
                             #undef CONVERT_GLO
                             #define CONVERT_GLO(a,b) ;
                           }
        | id QUOTE glists QUOTE
                           { /**read & alloc pointers (which point at a GLOBPTR* if introduced here)**/
                             g->num_pointers += 1;

                             sprintf (name, "%s", $1);
                             g->ch_pointers[ct_p] = strdup (name);

                             for (i = 0; i < g->num_pointers - 1; ++i)
                                 if  (! strcmp (g->ch_pointers[i], g->ch_pointers[ct_p]))
                                     fprintf (stderr, "\n\nWARNING: twice the same global pointer name %s!\n", g->ch_pointers[ct_p]);

                             g->pointers[ct_p]  = (GLOBPTR *) malloc (sizeof (GLOBPTR));
                             globptr = (GLOBPTR *)g->pointers[ct_p];
                             globptr->num_words = ct_c;
                             globptr->words     = (char **) malloc (ct_c * sizeof (char *));
                             for (int in_c = 0; in_c < ct_c; ++in_c)
                                 globptr->words[in_c] = strdup (names[in_c]);
                             globptr->data      = NULL;
                             globptr->int_val   = 0;

                             fprintf (stderr, "\nassigning global pointer %d: %s   ", ct_p, g->ch_pointers[ct_p]);

                             ct_c = 0;
                             ++ct_p;
                           }
        | id QUOTE QUOTE   {
                             g->num_pointers += 1;

                             sprintf (name, "%s", $1);
                             g->ch_pointers[ct_p] = strdup (name);

                             for (i = 0; i < g->num_pointers - 1; ++i)
                                 if  (! strcmp (g->ch_pointers[i], g->ch_pointers[ct_p]))
                                     fprintf (stderr, "\n\nWARNING: twice the same global pointer name %s!\n", g->ch_pointers[ct_p]);

                             g->pointers[ct_p]  = (GLOBPTR *) malloc (sizeof (GLOBPTR));
                             globptr = (GLOBPTR *)g->pointers[ct_p];
                             globptr->num_words = 0;
                             globptr->words     = (char **) malloc (1 * sizeof (char *));
                             globptr->words[0]  = NULL;
                             globptr->data      = NULL;
                             globptr->int_val   = 0;

                             fprintf (stderr, "\nassigning global pointer %d: %s   ", ct_p, g->ch_pointers[ct_p]);

                             ct_c = 0;
                             ++ct_p;
                           }
        ;

glists  : ilists
     /* | nlists */
        ;

ilists  : ilist
        | ilist ilists
        ;

ilist   : id               { sprintf (names[ct_c], "%s", $1);
                             ++ct_c;
                           }
        ;

/** no numbers allowed together with id's because yacc would complain: reduce/reduce conflicts
nlists  : nlist
        | nlist nlists
        ;

nlist   : num              { sprintf (names[ct_c], "%s", $1);
                             ++ct_c;
                           }
        ;
**/


        /**********all areas default: in pseudo-area A[areas]**********/

all     : ALL              {
                             A = (AREA *) malloc ((g->areas + 1) * sizeof (AREA));
                           }

alines  : aline
        | aline alines
        ;

aline:   id num            { /**read all area SCAlar default parameters (inner loop for all ct_ar)**/
                             #undef CONVERT_SCA
                             #define CONVERT_SCA(a,b) \
                             if  (!strcmp ($1, #b)) \
                                 for (ct_ar = 0; ct_ar < g->areas + 1; ++ct_ar) \
                                     A[ct_ar].b = (a)(atof($2));
                             #include "../kernel/parameters.h"
                             #undef CONVERT_SCA
                             #define CONVERT_SCA(a,b) ;
                           }
        ;


        /**********numbered areas**********/

area    : AREAS num        { ct_ar = atoi($2);
                           }
nlines  : nline
        | nline nlines
        ;

nline   : id num           { /**read destinct area SCAlar parameters (no loop for ct_ar; specified above)**/
                             if  (ct_ar >= g->areas)
                                 yyerror ("ct_ar overrun");
                             #undef CONVERT_SCA
                             #define CONVERT_SCA(a,b) \
                             if  (!strcmp ($1, #b)) \
                                 A[ct_ar].b = (a)(atof($2));
                             #include "../kernel/parameters.h"
                             #undef CONVERT_SCA
                             #define CONVERT_SCA(a,b) ;
                           }
        ;


        /**********series**********/

series  : condition SER    { fprintf (stderr, "\n1. series ");

                             se = si->se + si->anz_se;

                             si->anz_se += 1;
                             if  (si->anz_se > si->max_anz_se)
                                 fprintf (stderr, "\n\ntoo many series!!\n\n");

                             /**init of sweep number and pointer for realloc**/
                             se->sw = NULL;
                             se->anz_sw = 0;

                             oldcmd.anz_pointers = 0; /*init for first time read*/
                             oldcmd.ch_pointers  = (char **) malloc (100 * sizeof (char *));
                             oldcmd.pointers     = (GLOBPTR **) malloc (100 * sizeof (void *));

                             oldcmd.anz_from1    = 0; /*init for first time read*/
                             anz_ch_from1        = 0;
                             oldcmd.anz_from2    = 0; /*init for first time read*/
                             anz_ch_from2        = 0;
                             oldcmd.anz_from3    = 0; /*init for first time read*/
                             anz_ch_from3        = 0;

                             oldcmd.n_from1      = (int *) malloc (128 * sizeof(int));
                             oldcmd.n_from2      = (int *) malloc (128 * sizeof(int));
                             oldcmd.n_from3      = (int *) malloc (128 * sizeof(int));
                             oldcmd.ch_from1     = (char *) malloc (128*sizeof(char));
                             oldcmd.ch_from2     = (char *) malloc (128*sizeof(char));
                             oldcmd.ch_from3     = (char *) malloc (128*sizeof(char));

                             oldcmd.anz_quantums = 0;
                             oldcmd.anz_quant=(int *) malloc (128 * sizeof(int));
                             oldcmd.quantum = (DOUBLE **) malloc (32*sizeof(DOUBLE *));
                             for (i = 0; i < 32; ++i) {
                                 oldcmd.anz_quant[i] = 0;
                                 oldcmd.quantum[i] = (DOUBLE *) malloc (32*sizeof(DOUBLE));
                             }
                           }
        ;


ilen    : num              { fprintf (stderr, " ilen ");
                             se->ilen = atoi ($1);
                           }
        ;


        /**********sweeps**********/

sweeps  : sweep
        | sweeps sweep
        ;

sweep   : condition sweepword RI_ control _IR GI_ stays _IG {
                           }
        ;

sweepword: SWE             { fprintf (stderr, "\n2. sweepword  ");
                             /**allocate the current sweep**/
                             se->anz_sw += 1;

                             /*if  (se->sw != NULL) free (se->sw);*/ /**program and computer will crash if free&malloc used!*/
                             /*se->sw = (SWEEP *)malloc (se->anz_sw * sizeof(SWEEP));*/
                             se->sw = (SWEEP *)realloc (se->sw, se->anz_sw * sizeof(SWEEP));

                             fprintf (stderr, "   anz_sw=%d  ", se->anz_sw);

                             /**init stays**/
                             curr_sw = se->sw + se->anz_sw - 1;
                             curr_sw->st = NULL;
                             curr_sw->anz_st = 0;
                           }
        ;

control : num SEP num SEP id {fprintf (stderr, "\n3. control  ");
                             curr_sw->begin = atoi ($1);
                             curr_sw->end   = atoi ($3);
                             curr_sw->sw_update = strdup ($5);
                           }
        | num SEP num SEP num {fprintf (stderr, "\n3. control (single unit)  ");
                             curr_sw->begin = atoi ($1);
                             curr_sw->end   = atoi ($3);
                             curr_sw->sw_update = strdup ($5);
                           }
        ;

stays   : stay
        | stays stay
        ;

stay    : condition stayword commands
        ;

stayword: num RI_ id _IR   { fprintf (stderr, "\n4. stayword  ");
                             /**allocate the current stay**/
                             curr_sw->anz_st += 1;
                             i = curr_sw->anz_st;

                             /*if  (curr_sw->st != NULL) free (curr_sw->st);*/ /**program and computer will crash if free&malloc used!*/
                             /*curr_sw->st = (STAY *)malloc (i * sizeof (STAY));*/
                             curr_sw->st = (STAY *)realloc (curr_sw->st, i * sizeof (STAY));

                             /**fill in area, st_update**/
                             curr_sw->st[i-1].area      = atoi ($1);
                             curr_sw->st[i-1].st_update = $3[0];

                             /**allocate first command**/
                             /**at the end stay will have one cmd too much**/
                             curr_sw->st[i-1].cmd = (COMMAND *) malloc (sizeof (COMMAND));
                             curr_sw->st[i-1].anz_cmd = 0;
                           }
        ;

condition :                { /**empty condition**/
                             si->all_cond[si->num_all_cond] = (CONDITION *) malloc (sizeof (CONDITION));
                             si->all_cond[si->num_all_cond]->op_tag   = 0;
                             si->all_cond[si->num_all_cond]->left     = NULL;
                             si->all_cond[si->num_all_cond]->right    = NULL;
                             si->all_cond[si->num_all_cond]->test     = test_none;
                             si->num_all_cond += 1;
                           }
        | IFWORD RI_ num _IR { /**condition with one variable -- here a number**/
                             si->all_cond[si->num_all_cond] = (CONDITION *) malloc (sizeof (CONDITION));
                             si->all_cond[si->num_all_cond]->op_tag   = 1;
                             si->all_cond[si->num_all_cond]->left     = (int *) malloc (sizeof (int));
                             si->all_cond[si->num_all_cond]->left[0]  = atoi ($3);
                             si->all_cond[si->num_all_cond]->right    = NULL;
                             si->all_cond[si->num_all_cond]->test     = test_one;
                             si->num_all_cond += 1;
                           }
        | IFWORD RI_ id _IR { /**condition with one variable -- here pointer's number**/
                             int sel = -1;
                             for (i = 0; i < g->num_pointers /*- 1*/; ++i)
                                 if  (! strcmp (g->ch_pointers[i], $3))
                                     sel = i;
                             if  (sel == -1) {
                                 fprintf (stderr, "\nid %s not recognised as pointer\n", $3);
                                 exit (1);
                             } else {
                                 fprintf (stderr, "\npointer %s used as number", $3);
                                 si->all_cond[si->num_all_cond] = (CONDITION *) malloc (sizeof (CONDITION));
                                 si->all_cond[si->num_all_cond]->op_tag   = 1;
                                 si->all_cond[si->num_all_cond]->left     = &(g->pointers[sel]->int_val);
                                 si->all_cond[si->num_all_cond]->right    = NULL;
                                 si->all_cond[si->num_all_cond]->test     = test_one;
                                 si->num_all_cond += 1;
                             }
                           }
        | IFWORD RI_ ITERWORD _IR { /**condition with one variable -- here special word iter**/
                             si->all_cond[si->num_all_cond] = (CONDITION *) malloc (sizeof (CONDITION));
                             si->all_cond[si->num_all_cond]->op_tag   = 1;
                             si->all_cond[si->num_all_cond]->left     = &(g->iter);
                             si->all_cond[si->num_all_cond]->right    = NULL;
                             si->all_cond[si->num_all_cond]->test     = test_one;
                             si->num_all_cond += 1;
                           }
        | IFWORD RI_ id OP num _IR { /**condition with two variables**/
                             int sel = -1;
                             for (i = 0; i < g->num_pointers - 1; ++i)
                                 if  (! strcmp (g->ch_pointers[i], $3))
                                     sel = i;
                             if  (sel == -1) {
                                 fprintf (stderr, "\nid %s not recognised as pointer\n", $3);
                                 exit (1);
                             } else {
                                 fprintf (stderr, "\npointer %s used as number", $3);

                                 si->all_cond[si->num_all_cond] = (CONDITION *) malloc (sizeof (CONDITION));
                                 si->all_cond[si->num_all_cond]->left     = &(g->pointers[sel]->int_val);
                                 si->all_cond[si->num_all_cond]->right    = (int *) malloc (sizeof (int));
                                 si->all_cond[si->num_all_cond]->right[0] = atoi ($5);

                                 if  (op_tag == 2) {
                                     si->all_cond[si->num_all_cond]->op_tag = 2;
                                     si->all_cond[si->num_all_cond]->test   = test_equal;
                                 }
                                 if  (op_tag == 3) {
                                     si->all_cond[si->num_all_cond]->op_tag = 3;
                                     si->all_cond[si->num_all_cond]->test   = test_modulo;
                                 }
                                 if  (op_tag == 4) {
                                     si->all_cond[si->num_all_cond]->op_tag = 4;
                                     si->all_cond[si->num_all_cond]->test   = test_smaller;
                                 }
                                 if  (op_tag == 5) {
                                     si->all_cond[si->num_all_cond]->op_tag = 5;
                                     si->all_cond[si->num_all_cond]->test   = test_larger;
                                 }
                                 si->num_all_cond += 1;
                             }
                           }
        | IFWORD RI_ num OP id _IR { /**condition with two variables**/
                             int sel = -1;
                             for (i = 0; i < g->num_pointers - 1; ++i)
                                 if  (! strcmp (g->ch_pointers[i], $5))
                                     sel = i;
                             if  (sel == -1) {
                                 fprintf (stderr, "\nid %s not recognised as pointer\n", $5);
                                 exit (1);
                             } else {
                                 fprintf (stderr, "\npointer %s used as number", $5);

                                 si->all_cond[si->num_all_cond] = (CONDITION *) malloc (sizeof (CONDITION));
                                 si->all_cond[si->num_all_cond]->left     = (int *) malloc (sizeof (int));
                                 si->all_cond[si->num_all_cond]->left[0]  = atoi ($3);
                                 si->all_cond[si->num_all_cond]->right    = &(g->pointers[sel]->int_val);

                                 if  (op_tag == 2) {
                                     si->all_cond[si->num_all_cond]->op_tag = 2;
                                     si->all_cond[si->num_all_cond]->test   = test_equal;
                                 }
                                 if  (op_tag == 3) {
                                     si->all_cond[si->num_all_cond]->op_tag = 3;
                                     si->all_cond[si->num_all_cond]->test   = test_modulo;
                                 }
                                 if  (op_tag == 4) {
                                     si->all_cond[si->num_all_cond]->op_tag = 4;
                                     si->all_cond[si->num_all_cond]->test   = test_smaller;
                                 }
                                 if  (op_tag == 5) {
                                     si->all_cond[si->num_all_cond]->op_tag = 5;
                                     si->all_cond[si->num_all_cond]->test   = test_larger;
                                 }
                                 si->num_all_cond += 1;
                             }
                           }
        | IFWORD RI_ id OP ITERWORD _IR { /**condition with two variables**/
                             int sel = -1;
                             for (i = 0; i < g->num_pointers - 1; ++i)
                                 if  (! strcmp (g->ch_pointers[i], $3))
                                     sel = i;
                             if  (sel == -1) {
                                 fprintf (stderr, "\nid %s not recognised as pointer\n", $3);
                                 exit (1);
                             } else {
                                 fprintf (stderr, "\npointer %s used as number", $3);

                                 si->all_cond[si->num_all_cond] = (CONDITION *) malloc (sizeof (CONDITION));
                                 si->all_cond[si->num_all_cond]->left     = &(g->pointers[sel]->int_val);
                                 si->all_cond[si->num_all_cond]->right    = &(g->iter);

                                 if  (op_tag == 2) {
                                     si->all_cond[si->num_all_cond]->op_tag = 2;
                                     si->all_cond[si->num_all_cond]->test   = test_equal;
                                 }
                                 if  (op_tag == 3) {
                                     si->all_cond[si->num_all_cond]->op_tag = 3;
                                     si->all_cond[si->num_all_cond]->test   = test_modulo;
                                 }
                                 if  (op_tag == 4) {
                                     si->all_cond[si->num_all_cond]->op_tag = 4;
                                     si->all_cond[si->num_all_cond]->test   = test_smaller;
                                 }
                                 if  (op_tag == 5) {
                                     si->all_cond[si->num_all_cond]->op_tag = 5;
                                     si->all_cond[si->num_all_cond]->test   = test_larger;
                                 }
                                 si->num_all_cond += 1;
                             }
                           }
        | IFWORD RI_ ITERWORD OP id _IR { /**condition with two variables**/
                             int sel = -1;
                             for (i = 0; i < g->num_pointers - 1; ++i)
                                 if  (! strcmp (g->ch_pointers[i], $5))
                                     sel = i;
                             if  (sel == -1) {
                                 fprintf (stderr, "\nid %s not recognised as pointer\n", $5);
                                 exit (1);
                             } else {
                                 fprintf (stderr, "\npointer %s used as number", $5);

                                 si->all_cond[si->num_all_cond] = (CONDITION *) malloc (sizeof (CONDITION));
                                 si->all_cond[si->num_all_cond]->left     = &(g->iter);
                                 si->all_cond[si->num_all_cond]->right    = &(g->pointers[sel]->int_val);

                                 if  (op_tag == 2) {
                                     si->all_cond[si->num_all_cond]->op_tag = 2;
                                     si->all_cond[si->num_all_cond]->test   = test_equal;
                                 }
                                 if  (op_tag == 3) {
                                     si->all_cond[si->num_all_cond]->op_tag = 3;
                                     si->all_cond[si->num_all_cond]->test   = test_modulo;
                                 }
                                 if  (op_tag == 4) {
                                     si->all_cond[si->num_all_cond]->op_tag = 4;
                                     si->all_cond[si->num_all_cond]->test   = test_smaller;
                                 }
                                 if  (op_tag == 5) {
                                     si->all_cond[si->num_all_cond]->op_tag = 5;
                                     si->all_cond[si->num_all_cond]->test   = test_larger;
                                 }
                                 si->num_all_cond += 1;
                             }
                           }
        | IFWORD RI_ ITERWORD OP num _IR { /**condition with two variables**/
                             si->all_cond[si->num_all_cond] = (CONDITION *) malloc (sizeof (CONDITION));
                             si->all_cond[si->num_all_cond]->left     = &(g->iter);
                             si->all_cond[si->num_all_cond]->right    = (int *) malloc (sizeof (int));
                             si->all_cond[si->num_all_cond]->right[0] = atoi ($5);

                             if  (op_tag == 2) {
                                 si->all_cond[si->num_all_cond]->op_tag = 2;
                                 si->all_cond[si->num_all_cond]->test   = test_equal;
                             }
                             if  (op_tag == 3) {
                                 si->all_cond[si->num_all_cond]->op_tag = 3;
                                 si->all_cond[si->num_all_cond]->test   = test_modulo;
                             }
                             if  (op_tag == 4) {
                                 si->all_cond[si->num_all_cond]->op_tag = 4;
                                 si->all_cond[si->num_all_cond]->test   = test_smaller;
                             }
                             if  (op_tag == 5) {
                                 si->all_cond[si->num_all_cond]->op_tag = 5;
                                 si->all_cond[si->num_all_cond]->test   = test_larger;
                             }
                             si->num_all_cond += 1;
                           }
        | IFWORD RI_ num OP ITERWORD _IR { /**condition with two variables**/
                             si->all_cond[si->num_all_cond] = (CONDITION *) malloc (sizeof (CONDITION));
                             si->all_cond[si->num_all_cond]->left     = (int *) malloc (sizeof (int));
                             si->all_cond[si->num_all_cond]->left[0]  = atoi ($3);
                             si->all_cond[si->num_all_cond]->right    = &(g->iter);

                             if  (op_tag == 2) {
                                 si->all_cond[si->num_all_cond]->op_tag = 2;
                                 si->all_cond[si->num_all_cond]->test   = test_equal;
                             }
                             if  (op_tag == 3) {
                                 si->all_cond[si->num_all_cond]->op_tag = 3;
                                 si->all_cond[si->num_all_cond]->test   = test_modulo;
                             }
                             if  (op_tag == 4) {
                                 si->all_cond[si->num_all_cond]->op_tag = 4;
                                 si->all_cond[si->num_all_cond]->test   = test_smaller;
                             }
                             if  (op_tag == 5) {
                                 si->all_cond[si->num_all_cond]->op_tag = 5;
                                 si->all_cond[si->num_all_cond]->test   = test_larger;
                             }
                             si->num_all_cond += 1;
                           }
        | IFWORD RI_ num OP num _IR { /**condition with two variables**/
                             si->all_cond[si->num_all_cond] = (CONDITION *) malloc (sizeof (CONDITION));
                             si->all_cond[si->num_all_cond]->left     = (int *) malloc (sizeof (int));
                             si->all_cond[si->num_all_cond]->left[0]  = atoi ($3);
                             si->all_cond[si->num_all_cond]->right    = (int *) malloc (sizeof (int));
                             si->all_cond[si->num_all_cond]->right[0] = atoi ($5);

                             if  (op_tag == 2) {
                                 si->all_cond[si->num_all_cond]->op_tag = 2;
                                 si->all_cond[si->num_all_cond]->test   = test_equal;
                             }
                             if  (op_tag == 3) {
                                 si->all_cond[si->num_all_cond]->op_tag = 3;
                                 si->all_cond[si->num_all_cond]->test   = test_modulo;
                             }
                             if  (op_tag == 4) {
                                 si->all_cond[si->num_all_cond]->op_tag = 4;
                                 si->all_cond[si->num_all_cond]->test   = test_smaller;
                             }
                             if  (op_tag == 5) {
                                 si->all_cond[si->num_all_cond]->op_tag = 5;
                                 si->all_cond[si->num_all_cond]->test   = test_larger;
                             }
                             si->num_all_cond += 1;
                           }
        ;

OP      : EQ               { op_tag = 2; }
        | MOD              { op_tag = 3; }
        | SM               { op_tag = 4; }
        | LG               { op_tag = 5; }
        ;

commands: command
        | commands command
        ;

command : GI_ target SEP id SEP n_source SEP b_source SEP n_quantums _IG {
                             fprintf (stderr, "\n5. command  ");

                             /**allocate a new next command**/
                             i = curr_sw->anz_st - 1;

                             /*if  (curr_sw->st[i].cmd != NULL) free (curr_sw->st[i].cmd);*/ /**program and computer will crash if free&malloc used!*/
                             /*curr_sw->st[i].cmd = (COMMAND *) malloc ((curr_sw->st[i].anz_cmd+2) * sizeof (COMMAND));*/

fprintf (stderr, "\n reallocating %d times %d   ", (curr_sw->st[i].anz_cmd+2), sizeof (COMMAND));

                             curr_sw->st[i].cmd = (COMMAND *) realloc (curr_sw->st[i].cmd, (curr_sw->st[i].anz_cmd+2) * sizeof (COMMAND));

                             curr_sw->st[i].anz_cmd += 1;

                             /**set pointer to cmd which is now filled**/
                             cmd = curr_sw->st[i].cmd + curr_sw->st[i].anz_cmd - 1;

                             /**fill in target**/
                             cmd->ch_target = oldcmd.ch_target;

                             cmd->anz_pointers = oldcmd.anz_pointers;
                             cmd->ch_pointers  = (char **) malloc (100 * sizeof (char *));
                             cmd->pointers     = (GLOBPTR **) malloc (100 * sizeof (void *));
                             for (i = 0; i < cmd->anz_pointers; ++i) {
                                 cmd->ch_pointers[i] = strdup (oldcmd.ch_pointers[i]);
                                 cmd->pointers[i] = oldcmd.pointers[i];
                             }

                             /**fill in func_name**/
                             cmd->func_name = strdup ($4);

                             /**source data from oldcmd obtained below**/
                             cmd->anz_arguments = oldcmd.anz_arguments;

                             cmd->anz_from1 = oldcmd.anz_from1;

                             /**check consistency**/
                             if  ((cmd->anz_from1 != anz_ch_from1) && ( ! ch_from1_missing))
                                 fprintf(stderr,"\n\nch_from1 inconsistent!\n");

                             /**alloc mem for n_from1 and fill-in**/
                             cmd->n_from1 = (int *) malloc (cmd->anz_from1 * sizeof (int));
                             for (i = 0; i < cmd->anz_from1; ++i)
                                 cmd->n_from1[i] = oldcmd.n_from1[i];

                             /**repeat for anz_from2**/
                             cmd->anz_from2 = oldcmd.anz_from2;

                             /**check consistency**/
                             if  ((cmd->anz_from2 != anz_ch_from2) && ( ! ch_from2_missing))
                                 fprintf(stderr,"\n\nch_from2 inconsistent!\n");

                             /**alloc mem for n_from2 and fill-in**/
                             cmd->n_from2 = (int *) malloc (cmd->anz_from2 * sizeof (int));
                             for (i = 0; i < cmd->anz_from2; ++i)
                                 cmd->n_from2[i] = oldcmd.n_from2[i];

                             /**repeat for anz_from3**/
                             cmd->anz_from3 = oldcmd.anz_from3;

                             /**check consistency**/
                             if  ((cmd->anz_from3 != anz_ch_from3) && ( ! ch_from3_missing))
                                 fprintf(stderr,"\n\nch_from3 inconsistent!\n");

                             /**alloc mem for n_from3 and fill-in**/
                             cmd->n_from3 = (int *) malloc (cmd->anz_from3 * sizeof (int));
                             for (i = 0; i < cmd->anz_from3; ++i)
                                 cmd->n_from3[i] = oldcmd.n_from3[i];

                             /**alloc mem for ch_from1 and fill-in**/
                             cmd->ch_from1 = (char *) malloc ((cmd->anz_from1) * sizeof (char));
                             for (i = 0; i < cmd->anz_from1; ++i) {
                                 cmd->ch_from1[i] = oldcmd.ch_from1[i];

                                 /**default: take target as b_source**/
                                 if  (ch_from1_missing)
                                     cmd->ch_from1[i] = cmd->ch_target;
                             }

                             /**alloc mem for ch_from2 and fill-in**/
                             cmd->ch_from2 = (char *) malloc (cmd->anz_from2 * sizeof (char));
                             for (i = 0; i < cmd->anz_from2; ++i) {
                                 cmd->ch_from2[i] = oldcmd.ch_from2[i];

                                 /**default: take target as b_source**/
                                 if  (ch_from2_missing)
                                     cmd->ch_from2[i] = cmd->ch_target;
                             }

                             /**alloc mem for ch_from3 and fill-in**/
                             cmd->ch_from3 = (char *) malloc (cmd->anz_from3 * sizeof (char));
                             for (i = 0; i < cmd->anz_from3; ++i) {
                                 cmd->ch_from3[i] = oldcmd.ch_from3[i];

                                 /**default: take target as b_source**/
                                 if  (ch_from3_missing)
                                     cmd->ch_from3[i] = cmd->ch_target;
                             }

                             /**no last Arg: anz_quantums=1, anz_quant[0]= 1**/
                             cmd->anz_quantums = oldcmd.anz_quantums;
                             cmd->anz_quant = (int *) malloc (cmd->anz_quantums * sizeof (int));
                             cmd->quantum = (DOUBLE **) malloc (cmd->anz_quantums * sizeof(DOUBLE *));
                             for (i = 0; i < cmd->anz_quantums; ++i) {
                                 /**fill in 0.0 if argument missing**/
                                 if  (oldcmd.anz_quant[i] == 0) {
                                     oldcmd.anz_quant[i] = 1;
                                     oldcmd.quantum[i][0] = 0.0;
                                 }
                                 cmd->anz_quant[i] = oldcmd.anz_quant[i];
                                 cmd->quantum[i] = (DOUBLE *) malloc (cmd->anz_quant[i] * sizeof(DOUBLE));
                                 for (ii = 0; ii < cmd->anz_quant[i]; ++ii)
                                     cmd->quantum[i][ii] = oldcmd.quantum[i][ii];
                             }

                             /*init for next command which is read*/
                             oldcmd.anz_pointers = 0;
                             oldcmd.anz_from1    = 0;
                             anz_ch_from1        = 0;
                             oldcmd.anz_from2    = 0;
                             anz_ch_from2        = 0;
                             oldcmd.anz_from3    = 0;
                             anz_ch_from3        = 0;
                             ch_from1_missing    = 0;
                             ch_from2_missing    = 0;
                             ch_from3_missing    = 0;

                             oldcmd.anz_quantums = 0;
                             for (i = 0; i < 32; ++i)
                                 oldcmd.anz_quant[i] = 0;
                           }
        ;

target  : id               {
                             oldcmd.ch_target = $1[0];
                           }
        | id KOMMA p_terms {
                             oldcmd.ch_target = $1[0];
                           }
        ;

p_terms : p_term
        | p_term PLUS p_terms              /*****is this the same as  p_terms PLUS p_term ??? Find out !!! **/
        ;

p_term  : id               {
                             int newpointer = 1;
                             int thepointer;

                             sprintf (name, "%s", $1);
                             for (int i = 0; i < g->num_pointers; ++i)
                                 if  (! strcmp (name, g->ch_pointers[i])) {
                                     newpointer = 0;
                                     thepointer = i;

                                     fprintf (stderr, "\nre-visiting pointer %d: %s   ", i, g->ch_pointers[i]);
                                 }

                             if  (newpointer) {
                                 g->ch_pointers[g->num_pointers] = strdup (name);

                                 /**
                                 g->pointers[g->num_pointers]    = (void *) malloc (sizeof (int));            ** <-- allocationg pointer; sizeof int == sizeof address (hopefully)**
                                 **/
                                 g->pointers[g->num_pointers]  = (GLOBPTR *) malloc (sizeof (GLOBPTR));
                                 globptr = (GLOBPTR *)g->pointers[g->num_pointers];
                                 globptr->num_words = 0;
                                 globptr->words     = (char **) malloc (1 * sizeof (char *));
                                 globptr->words[0]  = NULL;
                                 globptr->data      = NULL;
                                 globptr->int_val   = 0;

                                 thepointer = g->num_pointers;

                                 fprintf (stderr, "\nassigning command pointer %d: %s   ", g->num_pointers, g->ch_pointers[g->num_pointers]);

                                 g->num_pointers += 1;
                             }

                             oldcmd.ch_pointers[oldcmd.anz_pointers] = g->ch_pointers[thepointer];
                             oldcmd.pointers[oldcmd.anz_pointers]    = g->pointers[thepointer];                /** <-- assigning pointer to command **/
                             oldcmd.anz_pointers += 1;
                           }
        ;

n_source: n_terms KOMMA n_terms2 KOMMA n_terms3 { oldcmd.anz_arguments = 3; }
        | n_terms KOMMA n_terms2 { oldcmd.anz_arguments = 2; }
        | n_terms          { oldcmd.anz_arguments = 1; }
        ;

n_terms : n_term
        | n_terms PLUS n_term
        ;

n_term  : num MULT RI_ n_modu _IR  { fprintf (stderr, "\nconsider later!"); }
        | num              { oldcmd.n_from1[oldcmd.anz_from1] = atoi ($1);
                             oldcmd.anz_from1 += 1;
                           }
        |                  { /**get area from above stay if not given**/
                             if  (oldcmd.anz_from1 == 0) {
                                 oldcmd.anz_from1 = 1;
                                 oldcmd.n_from1[0] = curr_sw->st[curr_sw->anz_st - 1].area;
                             } else {
                                 fprintf(stderr,"\nparser: anz_from1 wrong\n\n");
                             }
                           }
        ;

n_terms2: n_term2
        | n_terms2 PLUS n_term2
        ;

n_term2 : num MULT RI_ n_modu _IR  { fprintf (stderr, "\ntoo large term2"); }
        | num              { oldcmd.n_from2[oldcmd.anz_from2] = atoi ($1);
                             oldcmd.anz_from2 += 1;
                           }
        |                  { /**get area from above stay (as above)**/
                             if  (oldcmd.anz_from2 == 0) {
                                 oldcmd.anz_from2 = 1;
                                 oldcmd.n_from2[0] = curr_sw->st[curr_sw->anz_st - 1].area;
                             } else {
                                 fprintf(stderr,"\nparser: anz_fr2 wrong\n\n");
                             }
                           }
        ;

n_terms3: n_term3
        | n_terms3 PLUS n_term3
        ;

n_term3 : num MULT RI_ n_modu _IR  { fprintf (stderr, "\ntoo large term3"); }
        | num              { oldcmd.n_from3[oldcmd.anz_from3] = atoi ($1);
                             oldcmd.anz_from3 += 1;
                           }
        |                  { /**get area from above stay (as above)**/
                             if  (oldcmd.anz_from3 == 0) {
                                 oldcmd.anz_from3 = 1;
                                 oldcmd.n_from3[0] = curr_sw->st[curr_sw->anz_st - 1].area;
                             } else {
                                 fprintf(stderr,"\nparser: anz_fr3 wrong\n\n");
                             }
                           }
        ;

n_modu  : num              { }
        | n_modu PLUS num  { }
        ;




b_source: b_terms KOMMA b_terms2 KOMMA b_terms3 { if  (oldcmd.anz_arguments != 3)
                                      fprintf (stderr, "\nb_term has komma in excess\n");
                           }
        | b_terms KOMMA b_terms2 { if  (oldcmd.anz_arguments != 2)
                                      fprintf (stderr, "\nb_term has komma\n");
                           }
        | b_terms          { if  (oldcmd.anz_arguments != 1)
                                 fprintf (stderr, "\nb_term has no komma\n");
                           }
        ;

b_terms: b_term
        | b_terms PLUS b_term
        ;

b_term  : id MULT RI_ b_modu _IR  { fprintf (stderr, "\nconsider later!"); }
        | id               { oldcmd.ch_from1[anz_ch_from1] = $1[0];
                             anz_ch_from1 += 1;
                           }
        |                  { /**copy from ch_target later, when cmd completed**/
                             ch_from1_missing = 1;
                           }
        ;

b_terms2: b_term2
        | b_terms2 PLUS b_term2
        ;

b_term2: id MULT RI_ b_modu _IR  { fprintf (stderr, "\nconsider later!"); }
        | id               { oldcmd.ch_from2[anz_ch_from2] = $1[0];
                             anz_ch_from2 += 1;
                           }
        |                  { /**copy from ch_target later, when cmd completed**/
                             ch_from2_missing = 1;
                           }
        ;

b_terms3: b_term3
        | b_terms3 PLUS b_term3
        ;

b_term3: id MULT RI_ b_modu _IR  { fprintf (stderr, "\nconsider later!"); }
        | id               { oldcmd.ch_from3[anz_ch_from3] = $1[0];
                             anz_ch_from3 += 1;
                           }
        |                  { /**copy from ch_target later, when cmd completed**/
                             ch_from3_missing = 1;
                           }
        ;

b_modu: id                 { }
        | b_modu PLUS id   { }
        ;


n_quantums: n_quantums KOMMA n_quants { oldcmd.anz_quantums += 1; }
        | n_quants         { oldcmd.anz_quantums += 1; }
        ;

n_quants: n_quants PLUS num {
                              oldcmd.quantum[oldcmd.anz_quantums][oldcmd.anz_quant[oldcmd.anz_quantums]] = atof ($3);
                              oldcmd.anz_quant[oldcmd.anz_quantums] += 1; 
                           }
        | n_quants PLUS    {
                              oldcmd.quantum[oldcmd.anz_quantums][oldcmd.anz_quant[oldcmd.anz_quantums]] = 0.0;
                              oldcmd.anz_quant[oldcmd.anz_quantums] += 1; 
                           }
        | num              {
                              oldcmd.quantum[oldcmd.anz_quantums][oldcmd.anz_quant[oldcmd.anz_quantums]] = atof ($1);
                              oldcmd.anz_quant[oldcmd.anz_quantums] += 1;
                           }
        |                  {
                              oldcmd.quantum[oldcmd.anz_quantums][oldcmd.anz_quant[oldcmd.anz_quantums]] = 0.0;
                              oldcmd.anz_quant[oldcmd.anz_quantums] += 1;
                           }
        ;


id      : ID               { $$ = $1;
                           }
        ;


num     : NUMBER           { $$ = $1;
                           }
/**     | id               { fprintf(stderr, "line %d: cannot use '%s' as number.\n", lineno, $1);
                             yyerror ("non valid number");
                             exit(-1);
                           }
** had to be out-commented, because "RI_ num _IR" and "RI_ id _IR" are both used for condition => shift reduce conflict **
**/
        ;


%%





int lineno = 1;




/*****************************************************************************/


int test_none (int *left, int *right) {

    return 1;
}

int test_one (int *left, int *right) {

    return (*left > 0 ? 1 : 0);
}

int test_equal (int *left, int *right) {

    return (*left == *right ? 1 : 0);
}

int test_modulo (int *left, int *right) {

    return (*left % *right == 0 ? 1 : 0);
}

int test_smaller (int *left, int *right) {

    return (*left < *right ? 1 : 0);
}

int test_larger (int *left, int *right) {

    return (*left > *right ? 1 : 0);
}



/********************************* yyerror ***********************************/
/* my version of yyerror; will be invoked automatically.                     */

void yyerror(const char *msg){
  extern char *yytext;
  
  fprintf(stderr, "\nparser error in line ");
  fprintf(stderr, "%d: %s at '%s'\n", lineno, msg, yytext);

  exit (1);
}



/*
int main (int argc, char **argv) {
  extern FILE *yyin;
  if((int)argc > 2){
    fprintf(stderr, "%s: usage [infile]\n", argv[0]);
    exit(-1);
  }
  if((int)argc > 1){
    ** set argument 1 as yyin (scanner should read from that) **
    yyin = fopen(argv[1], "r");
    if(yyin == NULL){
      fprintf(stderr, "%s: cannot open %s\n", argv[0], argv[1]);
      exit(-1);
    }
  }
  ** parse yyin **
  yyparse();
  return 0;
}
*/



/********************************** alloc_g **********************************/
/* Allocates and initializes g, Z and s dependent on parameters in 1) fp.    */

PARAMS *alloc_g (FILE *fp) {
  extern FILE *yyin;

  char kommentarzeile[256];

  yyin = fp;
  if  (yyin == NULL) {
      fprintf(stderr, "\nalloc_g: cannot open file\n");
      exit(-1);
  }

  #undef CONVERT_GLO
  #undef CONVERT_SCA
  #define CONVERT_GLO(a,b)  ;
  #define CONVERT_SCA(a,b)  ;

  fgets (kommentarzeile, 256, yyin);
  fprintf (stderr, "\nalloc_g reads first line from yyin: %s", kommentarzeile);
  fseek (yyin, 0, SEEK_SET);

/*
#define YYDEBUG 1
yydebug = 1;

if there happens an error in line 1 with yytext = '', then think about libraries:  parser/y.tab.o  and  parser/lex.yy.o  have to be included!
*/

  yyparse();

  fprintf (stderr, "\nyyparse finished\n");

  fflush (stdin);
  fflush (fp);
  fflush (yyin);

  free (oldcmd.n_from1);
  free (oldcmd.n_from2);
  free (oldcmd.n_from3);
  free (oldcmd.ch_from1);
  free (oldcmd.ch_from2);
  free (oldcmd.ch_from3);

  free (oldcmd.anz_quant);
  for (i = 0; i < 32; ++i)
      free (oldcmd.quantum[i]);
  free (oldcmd.quantum);

  return g;
}


AREA *alloc_a () {

  return A;
}



SIMULATION *alloc_si () {

  return si;
}
