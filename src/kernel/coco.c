#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "coco.h"
#include "series.h"
#include "../parser/r.yacc.h"
#include "relax.h"                       /**for do_series**/
#include "vehicle.h"                     /**for printparams**/
#include "utils.h"                       /**for get_tmp_uname_filename**/


int main (int argc, char **argv) {

    PARAMS     *g;
    SIMULATION *si;

    AREA   *A;
    //DATA   *d;

    int  analyze        = 0;
    int  seed           = 1;
    char file[512];
    char *directory;
    char sets[1024]     = "\0";
    char directory_filename[512];
    char wovon[512];

    FILE *fp;
    int  pos_xy[2]      = {0, 0};
    int  i;

    directory = get_pwd ();  /**used only in this file to access (i) vlink and (ii) prae.exe and (iii) prae.file**/


    if  (argc == 1) {
        time_t now = time (NULL);

        fprintf (stderr, "%s [-analyze <0|1|2|3|5|6|9|10>] [-seed <1|n>]" , argv[0]);
        fprintf (stderr, " [-file <paramfile>] [-set <search> <change>]\n");

        sprintf (file, "%s/vlink", directory);

        fprintf (stderr, "\n... continuing, using link %s to the parameter file    ", file);

        do  {
            /* nothing */
        } while (time (NULL) < now + 2);
    }

    fprintf (stderr, "argc=%d, ", argc);

    for (i = 1; i < argc - 1; ) {
        int OK = 0;
        if  (! strcmp (argv[i], "-analyze")) {
            sscanf (argv[i+1], "%d", &analyze);
            i += 2;
            OK = 1;
        } else {
            if  (! strcmp (argv[i], "-seed")) {
                sscanf (argv[i+1], "%d", &seed);
                i += 2;
                OK = 1;
            } else {
                if  (! strcmp (argv[i], "-file")) {
                    sscanf (argv[i+1], "%s", file);
                    i += 2;
                    OK = 1;
                } else {
                    if  (! strcmp (argv[i], "-pos_xy")) {
                        sscanf (argv[i+1], "%d", &pos_xy[0]);
                        sscanf (argv[i+2], "%d", &pos_xy[1]);
                        i += 3;
                        OK = 1;
                    } else {
                        if  (! strcmp (argv[i], "-set")) {
                            strcat (sets, "-set ");
                            strcat (sets, argv[i+1]);
                            strcat (sets, " ");
                            strcat (sets, argv[i+2]);
                            strcat (sets, " ");
                            i += 3;
                            OK = 1;
                        }
                    }
                }
            }
        }
        if  (!OK)
            fprintf (stderr, "%s: wrong arguments!\n", argv[0]);
    }

    fprintf (stderr, "analyze=%d, seed=%d, file=%s, sets=%s", analyze, seed, file, sets);

    srand48 (seed);

    fprintf (stderr, "\nsystem call (prae-processor): ");
    sprintf (wovon, "%s/prae.exe %s %s > %s/prae.file", directory, file, sets, directory);
    fprintf (stderr, "%s\n", wovon);

/*hope this is not necessary anymore and LINUX will never crash again with an overly long command line ...
    if  ((strlen (wovon) == 0) || (strlen (wovon) > 127)) {
        FILE *fp = fopen ("/dev/tty", "w");
        fprintf (fp, "\07");
        fclose (fp);
        fprintf (stderr, "I refuse to execute a command with %d characters\n", strlen (wovon));
        exit (0);
    }
*/
    system (wovon);

    /**import and allocate g, s using yacc**/
    sprintf (directory_filename, "%s/prae.file", directory);
    fp = fopen (directory_filename, "r");
    if  (fp == NULL)
        fprintf (stderr, "\ncould not open %s!\n", directory_filename);

    g = alloc_g (fp);
    fclose (fp);

    si = alloc_si ();

    A = alloc_a ();
    init_a (A, g, si);

    assign_conditions (si);

    printparams (g, A, si, stderr);

    choose_functions (g, si);

    choose_pointers (si, A);

/*
    check_vehicle (si);

    if  ((d = (DATA *)malloc (sizeof (DATA))) == NULL)
        fprintf (stderr, "\nallocation failure for d\n");

    if  (g->data_import != NULL) {
        fprintf (stderr, "\nimport data: ");
        (*g->data_import)(g, d);
    } else {
        fprintf (stderr, "\nno data imported; hope you init data within series\n");
    }
*/

fprintf(stderr,"\n----------------------------------------------------------");

    do_simulation (g, A, si);

/*
    if  ((analyze >= 2) && (analyze <10))
        contrast_curves (g, Z, A, si, d, analyze, tmp_uname, pos_xy);
    if  (analyze == 10)
        curves_cuecomb (g, Z, A, si, d, analyze, tmp_uname);
*/

    printparams (g, A, si, stderr);

    free_a (A, g, si);
    free_se (si);
    free (si);
    free (g);

    return (0) ;
}
