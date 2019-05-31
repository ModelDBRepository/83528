#include <stdio.h>
#include <stdlib.h>
#include <math.h>     /**fuer fabs**/
#include <string.h>   /**fuer strstr**/

#include "../kernel/coco.h"
#include "../kernel/series.h"

/************************* exportP36_matrix **********************************/
/* Export 1) a (weight) matrix of sizes d_x * d_y * d_a * d_b to 7) "datei". */
/* with 6) format: 3 (integer), 6 (char), 9 (double with fwrite)             */

void exportP36_matrix (DOUBLE **S, int d_x, int d_y, int d_a, int d_b,
                      int format, char datei[256]) {

//fprintf (stderr, "\nexportP36_matrix begin ...");

    FILE *fp;
    int h = 2, i, j, a, b, k, hcol = 0, red, green, blue;
    double value, maxS, highS = -1000000, lowS = 1000000;

    /**compute maxS**/
    for (i = 0; i < d_x * d_y; ++i)
        for (j = 0; j < d_a * d_b; ++j) {
            highS = highS > S[i][j] ? highS : S[i][j];
            lowS  = lowS  < S[i][j] ? lowS  : S[i][j];
        }
    maxS  = fabs (highS) > fabs (lowS) ? fabs (highS) : fabs (lowS);

    if  (format == 3)      /**integer export: all values multiplied by 10000**/
         hcol = (int)(4500.0 * maxS);
    if  (format == 6)     /**character export: values scaled between [0:256]**/
         hcol = 60;
    if  (format == 9)     /**double export: no frames**/
         h = 0;

    if  (  (datei[0] != '/') || (datei[1] != 't') || (datei[2] != 'm')
        || (datei[3] != 'p') || (datei[4] != '/'))
        fprintf(stderr, "\n--> %s", datei);

    if  ((fp = fopen (datei, "w")) == NULL) {
        fprintf(stderr, "\nno write file %s\n", datei);  exit(0);
    }

    /**write the Px**/
    if  (format == 3)  fprintf (fp, "P3\n");
    else
       if  (format == 6)  fprintf (fp, "P6\n");
       else
          if  (format == 9)  fprintf (fp, "P9\n");
          else
              fprintf (stderr, "\nexport format must be 3 or 6 or 9\n");

    /**write maxS**/
    fprintf (fp, "# highS: %f  lowS: %f\n", highS, lowS);

    /**no parameters written yet**/

    /**write sizes**/
    fprintf (fp, "%d %d\n", h + (d_b + h) * d_y, h + (d_a + h) * d_x);

    /**write grey values**/
    if  (format == 3)
        fprintf (fp, "%d\n", (int)(10000.0 * maxS));
    if  (format == 6)
        fprintf (fp, "127\n");
    if  (format == 9)
        fprintf (fp, "1\n"); /**must write an int because import reads int!!**/

    for (i = 0; i < h; ++i)                          /**1st white line**/
        for (j = 0; j < h + (d_b + h) * d_y; ++j) {
            if  (format == 3)
                fprintf (fp, "%d %d %d ", hcol, hcol, hcol);
            if  (format == 6)
                fprintf (fp, "%c%c%c", hcol, hcol, hcol);
        }

    for (i = 0; i < d_x; ++i) {
         for (a = 0; a < d_a; ++a) {

             for (k = 0; k < h; ++k) {
                 if  (format == 3)
                     fprintf (fp, "%d %d %d ", hcol, hcol, hcol);
                 if  (format == 6)
                     fprintf (fp, "%c%c%c", hcol, hcol, hcol);
             }

             for (j = 0; j < d_y; ++j) {
                 for (b = 0; b < d_b; ++b) {

                         value = S[j + i * d_y][b + a * d_b];

                         if  (value < 0.0) {
                             red = (format == 3)
                                 ? (int)(-10000.0 * value)
                                 : (int)(-127.0 * value / maxS);
                             green = 0;
                             blue = 0;
                         } else {
                             red = 0;
                             green = (format == 3)
                                   ? (int)(10000.0 * value)
                                   : (int)(127.0 * value / maxS);
                             blue = 0;
                         }

                         if  (format == 3)
                             fprintf (fp, "%d %d %d ", red, green, blue);
                         if  (format == 6)
                             fprintf (fp, "%c%c%c", red, green, blue);
                         if  (format == 9)
                             fwrite (&value, sizeof (double), 1 , fp);
                 }

                 for (k = 0; k < h; ++k) {
                     if  (format == 3)
                         fprintf (fp, "%d %d %d ", hcol, hcol, hcol);
                     if  (format == 6)
                         fprintf (fp, "%c%c%c", hcol, hcol, hcol);
                 }
	     }
         }

         for (k = 0; k < h; ++k)         /**between + last white lines**/
             for (j = 0; j < h + (d_b + h) * d_y; ++j) {
                 if  (format == 3)
                     fprintf (fp, "%d %d %d ", hcol, hcol, hcol);
                 if  (format == 6)
                     fprintf (fp, "%c%c%c", hcol, hcol, hcol);
             }
    }

    fclose (fp);

//fprintf (stderr, "...exportP36_matrix end\n");
}

/************************* exportP36_col *************************************/
/* Export weight/act matrix of sizes d_x * d_y * d_a * d_b to "datei".       */
/* Format: 3 (integer) or 6 (char).  d_a must be as RGB (3x grey).           */
/* Used in observe_col.                                                      */

void exportP36_col (DOUBLE **S, int d_x, int d_y, int d_a, int d_b,
                   int format, char datei[256]) {
    FILE *fp;
    int h = 2, i, j, a, b, k, hcol = 0, red, green, blue;
    double maxS, highS = -1000000, lowS = 1000000;
    int height;

    height = d_a / 3;
    if  (d_a % 3) {
        fprintf (stderr, "\nUse exportP36_col only if d_a represents RGB.\n");
        exit (0);
    }

    /**compute maxS**/
    for (i = 0; i < d_x * d_y; ++i)
        for (j = 0; j < d_a * d_b; ++j) {
            highS = highS > S[i][j] ? highS : S[i][j];
            lowS  = lowS  < S[i][j] ? lowS  : S[i][j];
        }
    maxS  = fabs (highS) > fabs (lowS) ? fabs (highS) : fabs (lowS);

    if  (format == 3)      /**integer export: all values multiplied by 10000**/
         hcol = (int)(4500.0 * maxS);
    if  (format == 6)     /**character export: values scaled between [0:256]**/
         hcol = 80;

    if  (  (datei[0] != '/') || (datei[1] != 't') || (datei[2] != 'm')
        || (datei[3] != 'p') || (datei[4] != '/'))
        fprintf(stderr, "\n--> %s", datei);

    if  ((fp = fopen (datei, "w")) == NULL) {
        fprintf(stderr, "\nno write file %s\n", datei);  exit(0);
    }

    /**write the Px**/
    if  (format == 3)  fprintf (fp, "P3\n");
    else
       if  (format == 6)  fprintf (fp, "P6\n");
       else
          fprintf (stderr, "\nexport format must be 3 or 6\n");

    /**write maxS**/
    fprintf (fp, "# highS: %f  lowS: %f\n", highS, lowS);

    /**no parameters written yet**/

    /**write sizes**/
    fprintf (fp, "%d %d\n", h + (d_b + h) * d_y, h + (height + h) * d_x);

    /**write grey values**/
    if  (format == 3)
        fprintf (fp, "%d\n", (int)(10000.0 * maxS));
    if  (format == 6)
        fprintf (fp, "255\n");

    for (i = 0; i < h; ++i)                          /**1st white line**/
        for (j = 0; j < h + (d_b + h) * d_y; ++j) {
            if  (format == 3)
                fprintf (fp, "%d %d %d ", hcol, hcol, hcol);
            if  (format == 6)
                fprintf (fp, "%c%c%c", hcol, hcol, hcol);
        }

    for (i = 0; i < d_x; ++i) {
         for (a = 0; a < height; ++a) {

             for (k = 0; k < h; ++k) {
                 if  (format == 3)
                     fprintf (fp, "%d %d %d ", hcol, hcol, hcol);
                 if  (format == 6)
                     fprintf (fp, "%c%c%c", hcol, hcol, hcol);
             }

             for (j = 0; j < d_y; ++j) {
                 for (b = 0; b < d_b; ++b) {

                         double valueR = S[j + i * d_y][b +  a                    * d_b];
                         double valueG = S[j + i * d_y][b + (a + height)          * d_b];
                         double valueB = S[j + i * d_y][b + (a + height + height) * d_b];

                         /**only negative values disturb, but if all are positive, then subtraction distorts the true colors**/
                         double shiftS = lowS < 0.0 ? lowS : 0.0;

                         if  (format == 3) {
                             red   = (int)((valueR - shiftS) * 10000.0);
                             green = (int)((valueG - shiftS) * 10000.0);
                             blue  = (int)((valueB - shiftS) * 10000.0);
                         } else {
                             red   = (int)((valueR - shiftS) / (highS - shiftS) * 255.0);
                             green = (int)((valueG - shiftS) / (highS - shiftS) * 255.0);
                             blue  = (int)((valueB - shiftS) / (highS - shiftS) * 255.0);
                         }

                         if  (format == 3)
                             fprintf (fp, "%d %d %d ", red, green, blue);
                         if  (format == 6)
                             fprintf (fp, "%c%c%c", red, green, blue);
                 }

                 for (k = 0; k < h; ++k) {
                     if  (format == 3)
                         fprintf (fp, "%d %d %d ", hcol, hcol, hcol);
                     if  (format == 6)
                         fprintf (fp, "%c%c%c", hcol, hcol, hcol);
                 }
	     }
         }

         for (k = 0; k < h; ++k)         /**between + last white lines**/
             for (j = 0; j < h + (d_b + h) * d_y; ++j) {
                 if  (format == 3)
                     fprintf (fp, "%d %d %d ", hcol, hcol, hcol);
                 if  (format == 6)
                     fprintf (fp, "%c%c%c", hcol, hcol, hcol);
             }
    }

    fclose (fp);
}


/************************* importP36_matrix **********************************/
/* Import 1) a (weight) matrix of sizes d_x * d_y * d_a * d_b                */
/* from filename 7) "datei". Automatic format detection (3 or 6 or 9).       */

void importP36_matrix (double **S, int d_x, int d_y, int d_a, int d_b,
                       char datei[256]) {
    FILE *fp;
    char kommentarzeile[256];

    int format = 0, hcol1, hcol2, hcol3;
    int h = 2, i, j, a, b, k, red = 0, green = 0, blue = 0, maxgrey;
    char chcol1, chcol2, chcol3, cred = 0, cgreen = 0, cblue = 0;
    double value = 0.0, maxS, highS, lowS;

    fprintf(stderr, "\nreading %s ", datei);

    if  ((fp = fopen (datei, "r")) == NULL) {
        fprintf(stderr, "\nno read file %s\n", datei);  exit(0);
    }

    /**read the Px**/
    fgets (kommentarzeile, 256, fp);
    if  (! strcmp (kommentarzeile, "P3\n"))
        format = 3;
    else {
        if  (! strcmp (kommentarzeile, "P6\n"))
            format = 6;
        else {
           if  (! strcmp (kommentarzeile, "P9\n")) {
               format = 9;
               h = 0;        /**no frames in format P9**/
           } else
               fprintf (stderr, " wrong header: %s", kommentarzeile);
        }
    }
    fprintf (stderr, " format=P%d ", format);

    /**read maxS**/
    fgets (kommentarzeile, 256, fp);
    sscanf  (kommentarzeile, "# highS: %lf lowS: %lf", &highS, &lowS);
    maxS  = fabs (highS) > fabs (lowS) ? fabs (highS) : fabs (lowS);
    fprintf(stderr, "got maxS: highS=%f, lowS=%f\n", highS, lowS);

    /**no other comment lines can be read**/

    /**read sizes**/
    fscanf (fp, "%d %d\n", &a, &b);

    if  ((a != h + (d_b + h) * d_y) || (b != h + (d_a + h) * d_x)) {
        fprintf (stderr, " sizes don't match in file %s\n", datei);
        exit (1);
    }

    /**read grey values ** this is an int, so write an int in export only !!**/
    fscanf (fp, "%d\n", &maxgrey);

    if  ((format == 3) && (fabs (10000.0 * maxS - maxgrey) > 1.0))
        fprintf (stderr, "\n\n\ninconsistent maximum grey value\n\n\n");
    if  ((format == 6) && (maxgrey != 127))
        fprintf (stderr, "\n\n\ninconsistent maximum grey value\n\n\n");

    fprintf (stderr, "pgm-size: %dx%d  grey: %d  ", a, b, maxgrey);
    fprintf (stderr, "weight-size: %dx%d x %dx%d ", d_x, d_y, d_a, d_b);

    for (i = 0; i < h; ++i)                          /**1st white line**/
        for (j = 0; j < h + (d_b + h) * d_y; ++j) {
            if  (format == 3)
                fscanf (fp, "%d %d %d ", &hcol1, &hcol2, &hcol3);
            if  (format == 6)
                fscanf (fp, "%c%c%c", &chcol1, &chcol2, &chcol3);
        }

    for (i = 0; i < d_x; ++i) {
         for (a = 0; a < d_a; ++a) {

             for (k = 0; k < h; ++k) {
                 if  (format == 3)
                     fscanf (fp, "%d %d %d ", &hcol1, &hcol2, &hcol3);
                 if  (format == 6)
                     fscanf (fp, "%c%c%c", &chcol1, &chcol2, &chcol3);
             }

             for (j = 0; j < d_y; ++j) {
                 for (b = 0; b < d_b; ++b) {

                     if  (feof (fp))
                         fprintf (stderr, "\n file end reached at x=%d a=%d y=%d b=%d\n", i, a, j, b);

                     if  (format == 3) {
                         fscanf (fp, "%d %d %d ", &red, &green, &blue);

                         if  ((red != 0) && (green == 0))
                             value = - (double)(red) / 10000.0;
                         if  ((red == 0) && (green != 0))
                             value = - (double)(red) / 10000.0;
                         if  ((red == 0) && (green == 0))
                             value = 0.0;
                         if  ((red != 0) && (green != 0))
                             fprintf (stderr, "warning: weights inconsistent");
                     }

                     if  (format == 6) {
                         fscanf (fp, "%c%c%c", &cred, &cgreen, &cblue);

                         if  ((cred != 0) && (cgreen == 0))
                             value = - (double)(cred) / 127.0 * maxS;
                         if  ((cred == 0) && (cgreen != 0))
                             value =   (double)(cgreen) / 127.0 * maxS;
                         if  ((cred == 0) && (cgreen == 0))
                             value = 0.0;
                         if  ((cred != 0) && (cgreen != 0))
                             fprintf (stderr, "warning: weights inconsistent");
                     }

                     if  (format == 9)
                         fread (&value, sizeof (double), 1, fp);

                     if  ((blue != 0) || (cblue != 0))
                         fprintf (stderr, " warning: hidden info in weights ");

                     S[j + i * d_y][b + a * d_b] = value;

                     if  (i == 0 && j == 0 && a == 0 && b == 0)
                         fprintf (stderr, "\nfirst value read is: %f", value);

                     if  (fabs (value) > 10000.0) {
                         fprintf (stderr, "\nWarning: large value at index [%d][%d][%d][%d]: %f\n", i, j, a, b, value);
                         exit (0);
                     }
                 }

                 for (k = 0; k < h; ++k) {
                     if  (format == 3)
                         fscanf (fp, "%d %d %d ", &hcol1, &hcol2, &hcol3);
                     if  (format == 6)
                         fscanf (fp, "%c%c%c", &chcol1, &chcol2, &chcol3);
                 }
	     }
         }

         for (k = 0; k < h; ++k)         /**between + last white lines**/
             for (j = 0; j < h + (d_b + h) * d_y; ++j) {
                 if  (format == 3)
                     fscanf (fp, "%d %d %d ", &hcol1, &hcol2, &hcol3);
                 if  (format == 6)
                     fscanf (fp, "%c%c%c", &chcol1, &chcol2, &chcol3);
             }
    }

    { /**just a test**/
      double max = -999.9;
      double min =  999.9;
      for (i = 0; i < d_x * d_y; ++i)
          for (j = 0; j < d_a * d_b; ++j) {
              max = S[i][j] > max ? S[i][j] : max;
              min = S[i][j] < min ? S[i][j] : min;
          }
      fprintf (stderr, " max=%f min=%f ", max, min);
      if  ((fabs (highS - max) > 0.01) || (fabs (lowS  - min) > 0.01)) {
          fprintf (stderr, "\n\nWarning: inconsistent max at weight import in %s:", datei);
          fprintf (stderr, "\n         Header: highS=%f lowS=%f  Values: max=%f min=%f\n\n", highS, lowS, max, min);
      }
    }

    fprintf (stderr, " file read\n");

    fclose (fp);
}




/**************************** observe_act ************************************/
/* Exports cmd->S_target as weight-matrix-imagefile.                         */

DOUBLE observe_act (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {

    int timespan;
    char fullname[256];
    int inarea = cmd->n_from1[0];

    sprintf (fullname, "%s/obs_%c_%d.pnm", cmd->pointers[0]->words[0], cmd->ch_target, inarea);

    if  ((timespan = end - begin) <= 0)
        fprintf (stderr, "\nso sense for observe_act without time\n");

    exportP36_matrix (cmd->S_target + begin, 1, timespan, A[inarea].d_a, A[inarea].d_b, 6, fullname);

    return (DOUBLE)(0);
}




/**************************** observe_col ************************************/
/* From but different from observe_act:                                      */
/* (i) takes S/W_from but writes file using B_target                         */
/* (ii) assumes RGB input (3x size) and calls exportP3_col                   */

void observe_col (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {

    int timespan = end - begin;
    char fullname[256];
    int inarea = cmd->n_from1[0];

/**is it a weight? (must be W or V)
        if  ((cmd->b_from1[0] == 'W') || (cmd->b_from1[0] == 'V')) {
            sprintf (fullname, "%s/obs_%c_%d_%d.pnm", tmp_uname, cmd->b_target, cmd->area, inarea);
            if  (cmd->b_from1[0] == 'W')
                exportP36_col (A[cmd->area].W[inarea], z->d_a, z->d_b, Z[inarea].d_a, Z[inarea].d_b, 6, fullname);
            else
                exportP36_col (A[cmd->area].V[inarea], z->d_a, z->d_b, Z[inarea].d_a, Z[inarea].d_b, 6, fullname);
        } else
**/

        {
            sprintf (fullname, "%s/obs_%c_%d.pnm", cmd->pointers[0]->words[0], cmd->ch_target, inarea);

            if  (timespan <= 0)
                fprintf (stderr, "\nso sense for observe_col without time\n");

            exportP36_col (cmd->S_from1[0] + begin, 1, timespan, A[inarea].d_a, A[inarea].d_b, 6, fullname);
        }
}



/**************************** observe_animgif ********************************/
/* Writes the activations from ct_t = begin to end as animated gif.          */
/* Creates file with ending .gif and auxiliary subdirectory called animgif.  */
/* q[0][0]=3 for RGB retina, else traditional act export like observe_act.   */
/* Calls gifsicle. View with /usr/local/bin/gifview -a animfile.gif          */

void observe_animgif (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {

    char fullname[256];
    char fullnamegif[256];
    int inarea = cmd->n_from1[0];
    char wovon[512];

    static int firsttime = 1;
    if  (firsttime) {
        sprintf (wovon, "mkdir %s/animgif", cmd->pointers[0]->words[0]);
        fprintf (stderr, "%s\n", wovon);
        system (wovon);
        sprintf (wovon, "rm %s/animgif/*", cmd->pointers[0]->words[0]);
        fprintf (stderr, "%s\n", wovon);
        system (wovon);
        firsttime = 0;
    }

    for (int ct_t = begin; ct_t < end; ++ct_t) {

        sprintf (fullname, "%s/animgif/obs_%c_%d__%3d.pnm", cmd->pointers[0]->words[0], cmd->ch_target, inarea, ct_t);
        sprintf (fullnamegif, "%s/animgif/obs_%c_%d__%3d.gif", cmd->pointers[0]->words[0], cmd->ch_target, inarea, ct_t);

        /**replace the spaces (occurred through %3d) by 0's**/
        for (unsigned int i = 0; i < strlen (fullname); ++i)
            if  (fullname[i] == ' ')
                fullname[i] = '0';
        for (unsigned int i = 0; i < strlen (fullnamegif); ++i)
            if  (fullnamegif[i] == ' ')
                fullnamegif[i] = '0';


        if  (cmd->quantum[0][0] == 3)
            exportP36_col    (cmd->S_from1[0] + ct_t, 1, 1, A[inarea].d_a, A[inarea].d_b, 6, fullname);
        else
            exportP36_matrix (cmd->S_from1[0] + ct_t, 1, 1, A[inarea].d_a, A[inarea].d_b, 6, fullname);


        sprintf (wovon, "convert %s %s", fullname, fullnamegif);
        //fprintf (stderr, "%s\n", wovon);
        system (wovon);
    }

    sprintf (wovon, "/usr/local/bin/gifsicle --delay 10 --loopcount=forever --colors 255 -o %s/obs_%c_%d.gif %s/animgif/obs_%c_%d__???.gif",
                                                                    cmd->pointers[0]->words[0], cmd->ch_target, inarea, cmd->pointers[0]->words[0], cmd->ch_target, inarea);
    //fprintf (stderr, "%s\n", wovon);
    fprintf (stderr, " gifsicle ");
    system (wovon);
}



/**************************** observe_animgif_iter ***************************/
/* Copies & "convert"s an exported pnm file to gif in animgif subdir.        */
/* Files get ending with ascending number for later creation of an animgif.  */
/* q[0][0] =0 for act; =1 for normal weights; =2 for topo weights            */

void observe_animgif_iter (PARAMS *g, AREA *A, COMMAND *cmd, int, int) {

    char fullname[256];
    char fullnamegif[256];
    int inarea = cmd->n_from1[0];
    char wovon[512];

    static int firsttime = 1;
    if  (firsttime) {
        sprintf (wovon, "mkdir %s/animgif", cmd->pointers[0]->words[0]);
        fprintf (stderr, "%s\n", wovon);
        system (wovon);
        sprintf (wovon, "rm %s/animgif/*", cmd->pointers[0]->words[0]);
        fprintf (stderr, "%s\n", wovon);
        system (wovon);
        firsttime = 0;
    }

    static int count = 0;
    count += 1;

    /**for activations**/
    if  (cmd->quantum[0][0] == 0)
        sprintf (fullname, "%s/obs_%c_%d.pnm", cmd->pointers[0]->words[0], cmd->ch_target, inarea);

    /**for normal weights**/
    if  (cmd->quantum[0][0] == 1)
        sprintf (fullname, "%s/obs_w_%d_%d.pnm", cmd->pointers[0]->words[0], cmd->area, inarea);

    /**for topo weights**/
    if  (cmd->quantum[0][0] == 2)
        sprintf (fullname, "%s/topo_w_%d_%d.pnm", cmd->pointers[0]->words[0], cmd->area, inarea);

    sprintf (fullnamegif, "%s/animgif/itergif_%c_%d_%d__%3d.gif", cmd->pointers[0]->words[0], cmd->ch_target, cmd->area, inarea, count);

    /**replace the spaces (occurred through %3d) by 0's**/
    for (unsigned int i = 0; i < strlen (fullnamegif); ++i)
        if  (fullnamegif[i] == ' ')
            fullnamegif[i] = '0';

    sprintf (wovon, "convert %s %s", fullname, fullnamegif);
    fprintf (stderr, "\n%s    ", wovon);
    system (wovon);

    fprintf (stderr, "\n\nYou need to execute:\n");
    fprintf (stderr, "/usr/local/bin/gifsicle --delay 10 --loopcount=forever --colors 255 -o %s/itergif_%c_%d_%d.gif %s/animgif/itergif_%c_%d_%d__???.gif",
                                                                                             cmd->pointers[0]->words[0], cmd->ch_target, cmd->area, inarea,
                                                                                             cmd->pointers[0]->words[0], cmd->ch_target, cmd->area, inarea);
    fprintf (stderr, "\nAnd then:  rm %s/animgif/itergif_%c_%d_%d__???.gif  ",            cmd->pointers[0]->words[0], cmd->ch_target, cmd->area, inarea);
}



/**************************** observe_getc ***********************************/

DOUBLE observe_getc (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {

    fprintf (stderr, "\npress RETURN");
    getc (stdin);
    return (DOUBLE)(0);
}



/**************************** observe_gnu_act ********************************/
/* Writes activation values into file. Modus depends on area geometry:       */
/* For plot with gnuplot if 1-dim area, or for splot if 2-dim; at one time.  */
/* If 1-neuron-area, then writes time course of acts till rlen.              */
/* q[0][0]=time at which to write activations (ignored if 1-neuron-area)     */
/* q[1][0]=1 (opt) writes also gnuplot commands to stdout for easy plotting  */
/* q[2][0] (optional) additional file name natural number (for batch script) */

DOUBLE observe_gnu_act (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {

    char fullname[256];
    FILE *fp;

    /**for activations only**/
    if  (cmd->anz_quantums > 2)
        sprintf (fullname, "%s/obs_%c_%d__%d.dat", cmd->pointers[0]->words[0], cmd->ch_target, cmd->area, (int)(cmd->quantum[2][0]));
    else
        sprintf (fullname, "%s/obs_%c_%d.dat", cmd->pointers[0]->words[0], cmd->ch_target, cmd->area);

    if  ((fp = fopen (fullname, "w")) == NULL) {
        fprintf(stderr, "\nobserve_gnu_act: no write file %s\n", fullname);
        exit(0);
    }

    int d_a = A[cmd->area].d_a;
    int d_b = A[cmd->area].d_b;

    int time = (int)cmd->quantum[0][0];

    if  ((d_a == 1) || (d_b == 1)) {
        if  ((d_a == 1) && (d_b == 1))
            for (int t = begin; t < end; ++t)
                fprintf (fp, "%f\n", cmd->S_from1[0][t][0]);                /**area has 1 neuron:   write time course over rlen**/
        else
            for (int i = 0; i < d_a * d_b; ++i)
                fprintf (fp, "%f\n", cmd->S_from1[0][time][i]);             /**area has 1 col/row:  writes activations in one column**/
    } else {
        for (int X = 0; X < d_a; X++) {
            for (int Y = 0; Y < d_b; Y++)
                fprintf (fp, "%f ", cmd->S_from1[0][time][X*d_b + Y]);      /**area has 2dim shape: writes activations in several columns, for splot**/
            fprintf (fp, "\n");
        }
    }

    fclose (fp);

    if  (cmd->anz_quantums > 1)
    if  (cmd->quantum[1][0] == 1) {
        if  ((d_a == 1) || (d_b == 1))
            fprintf (stdout, "\ngnuplot\nplot \"%s\"\n", fullname);
        else
            fprintf (stdout, "\ngnuplot\nsplot \"%s\"\n", fullname);
    }

    return (DOUBLE)0;
}

/**************************** observe_gnu_two_acts ***************************/
/* Writes acts of two areas (S_from1[0/1]) into file; for each neuron a row. */
/* For splot with gnuplot, can be polar if 1st value is an angle.            */
/* q[0][0]=1 (opt) writes also gnuplot commands to stdout for easy plotting  */
/* if q[1][0], then omit values if S_from2[0] area values are higher.        */

DOUBLE observe_gnu_two_acts (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {

    fprintf (stderr, "\nobserve_gnu_two_acts ");

    char fullname[256];
    FILE *fp;

    /**for activations only**/
    sprintf (fullname, "%s/obs_%c-%c_%d.dat", cmd->pointers[0]->words[0], cmd->ch_from1[0], cmd->ch_from1[1], cmd->area);

    if  ((fp = fopen (fullname, "w")) == NULL) {
        fprintf(stderr, "\nobserve_gnu_act: no write file %s\n", fullname);
        exit(0);
    }

    int d_n = A[cmd->area].d_a * A[cmd->area].d_b;

    /**exclude scaffold for discarding certain values**/
    int  *exclude = (int *)calloc (d_n, sizeof (int));
    if  (cmd->anz_quantums > 1) {
        if  (cmd->anz_from2 != 1)
            fprintf (stderr, "\nobserve_gnu_two_acts needs an S_from2 here!\n");
        fprintf (stderr, "\nDiscarding the following neurons: ");
        int ct = 0;
        for (int i = 0; i < d_n; ++i)
            if  (cmd->S_from2[0][begin][i] > cmd->quantum[1][0]) {
                exclude[i] = 1;
                ct ++;
                fprintf (stderr, "%d ", i);
            }
        fprintf (stderr, " hence, discarding %d of %d values ", ct, d_n);
    }

    for (int i = 0; i < d_n; ++i)
        if  (!exclude[i])
            fprintf (fp, "%f %f\n", cmd->S_from1[0][begin][i], cmd->S_from1[1][begin][i]);   /**writes both activations in one row**/
           /*fprintf (fp, "%f %f\n", (cmd->S_from1[0][begin][i]>cmd->S_from1[1][begin][i]?cmd->S_from1[0][begin][i]:cmd->S_from1[1][begin][i]),
                                     (cmd->S_from1[0][begin][i]>cmd->S_from1[1][begin][i]?cmd->S_from1[1][begin][i]:cmd->S_from1[0][begin][i]));*/

    fclose (fp);

    if  (cmd->quantum[0][0] == 1) {
        fprintf (stdout, "\ngnuplot\nplot \"%s\"\n", fullname);
    }

    return (DOUBLE)0;
}


/**************************** observe_act_hist *******************************/
/* Writes histogram of S_from1 act's into data file for gnuplot.             */
/* If S_from2 given, it can determine, which neurons' acts to discard.       */
/* q00: number of bins                                                       */
/* q10: threshold of S_from2 values over which S_from1 values discarded      */

DOUBLE observe_act_hist (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {

    char fullname[256];
    FILE *fp;

    fprintf (stderr, "\nobserve_act_hist ");

    /**for activations only**/
    sprintf (fullname, "%s/hist_%c_%d.dat", cmd->pointers[0]->words[0], cmd->ch_target, cmd->area);

    if  ((fp = fopen (fullname, "w")) == NULL) {
        fprintf(stderr, "\nobserve_gnu_hist: no write file %s\n", fullname);
        exit(0);
    }

    int bins      = (int)cmd->quantum[0][0];

    int      time = begin; /**re-program if other needed!**/
    int       d_n = A[cmd->area].d_a * A[cmd->area].d_b;
    int  *exclude = (int *)calloc (d_n, sizeof (int));
    int *bincount = (int *)calloc (bins, sizeof (int));

    /**exclude scaffold for discarding certain values**/
    if  (cmd->anz_quantums > 1) {
        if  (cmd->anz_from2 != 1)
            fprintf (stderr, "\nobserve_act_hist needs an S_from2 here!\n");
        fprintf (stderr, "\nDiscarding the following neurons: ");
        int ct = 0;
        for (int i = 0; i < d_n; ++i)
            if  (cmd->S_from2[0][time][i] > cmd->quantum[1][0]) {
                exclude[i] = 1;
                ct ++;
                fprintf (stderr, "%d ", i);
            }
        fprintf (stderr, " hence, discarding %d of %d values ", ct, d_n);
    }

    /**get Smax and Smin**/
    DOUBLE Smax = -99.99, Smin = 99.99;
    int index_Smin = -1, index_Smax = -1;
    for (int i = 0; i < d_n; ++i)
        if  (!exclude[i]) {
            if  (cmd->S_from1[0][time][i] > Smax) {
                Smax = cmd->S_from1[0][time][i];
                index_Smax = i;
            }
            if  (cmd->S_from1[0][time][i] < Smin) {
                Smin = cmd->S_from1[0][time][i];
                index_Smin = i;
            }
        }

/*Smin = 0.0;*/
    
    fprintf (stderr, "\nSmax=%f at neuron %d, Smin=%f at neuron %d ", Smax, index_Smax, Smin, index_Smin);

    DOUBLE  interval = (Smax - Smin) / bins;

    /**count values in bins**/
    for (int i = 0; i < d_n; ++i) {
        DOUBLE lower = Smin;
        DOUBLE upper = Smin;

        for (int b = 0; b < bins; ++b) {
            upper += interval;
            if  (b < bins - 1)
                if  ((cmd->S_from1[0][time][i] >= lower) && (cmd->S_from1[0][time][i] < upper) && (!exclude[i]))
                    bincount[b] += 1;
            if  (b == bins - 1)
                if  ((cmd->S_from1[0][time][i] >= lower) && (cmd->S_from1[0][time][i] <= Smax) && (!exclude[i]))
                    bincount[b] += 1;
            lower += interval;
        }
    }

    /**print out**/
    for (int i = 0; i < bins; ++i)
        fprintf (fp, "\n%d ", bincount[i]);

    fclose (fp);

    /**gnuplot commands for a help**/
    fprintf (stderr, "\ngnuplot");
    fprintf (stderr, "\nset style data histeps");
    fprintf (stderr, "\nset xrange [-0.5:%.1f]", bins - 0.5);
    fprintf (stderr, "\nset xtics (\"%.2f\" -0.5, \"%.2f\" %f, \"%.2f\" %f)", Smin, 0.5*(Smin+Smax), 0.5*bins-0.5, Smax, bins-0.5);
    fprintf (stderr, "\nplot \"%s\" title \"hist_%c_%d\"\n", fullname, cmd->ch_target, cmd->area);

    free (exclude);
    free (bincount);
    return (DOUBLE)0;
}



/**************************** observe_phase_hist *****************************/
/* Writes histogram of S_from1 act's overS_from2 into data file for gnuplot. */
/* q00/01: interval boundaries of phase                                      */
/* q10: phase which to center                                                */
/* q20: number of bins                                                       */
/* q30: optional int value used for file name                                */

DOUBLE observe_phase_hist (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {

    char fullname[256];
    FILE *fp;

    fprintf (stderr, "\nobserve_phase_hist ");

    /**for activations only**/
    if  (cmd->anz_quantums > 3)
        sprintf (fullname, "%s/phase_%c_%d__%d.dat", cmd->pointers[0]->words[0], cmd->ch_target, cmd->area, (int)(cmd->quantum[3][0]));
    else
        sprintf (fullname, "%s/phase_%c_%d.dat", cmd->pointers[0]->words[0], cmd->ch_target, cmd->area);

    if  ((fp = fopen (fullname, "w")) == NULL) {
        fprintf(stderr, "\nobserve_phase_hist: no write file %s\n", fullname);
        exit(0);
    }

    int bins         = (int)cmd->quantum[2][0];

    int         time = begin; /**re-program if other needed!**/
    int          d_n = A[cmd->area].d_a * A[cmd->area].d_b;
    DOUBLE *bincount = (DOUBLE *)calloc (bins, sizeof (DOUBLE));

    /**get PHmax and PHmin, only for warning if out of interval**/
    DOUBLE PHmax = -99.99, PHmin = 99.99;
    for (int i = 0; i < d_n; ++i) {
        if  (cmd->S_from2[0][time][i] > PHmax)
            PHmax = cmd->S_from1[0][time][i];
        if  (cmd->S_from2[0][time][i] < PHmin)
            PHmin = cmd->S_from1[0][time][i];
    }

    fprintf (stderr, "  PHmax=%f, PHmin=%f ", PHmax, PHmin);
    if  ((PHmin < cmd->quantum[0][0]) || (PHmax > cmd->quantum[0][1]))
        fprintf (stderr, "\n\n\nWARNING: phase out of interval boundaries\n\n");

    PHmin = cmd->quantum[0][0];
    PHmax = cmd->quantum[0][1];

    DOUBLE  interval = (PHmax - PHmin) / bins;

    /**count values in bins**/
    for (int i = 0; i < d_n; ++i) {
        DOUBLE lower = PHmin;
        DOUBLE upper = PHmin;

        for (int b = 0; b < bins; ++b) {
            upper += interval;
            if  (b < bins - 1)
                if  ((cmd->S_from2[0][time][i] >= lower) && (cmd->S_from2[0][time][i] < upper))
                    bincount[b] += cmd->S_from1[0][time][i];
            if  (b == bins - 1)
                if  ((cmd->S_from2[0][time][i] >= lower) && (cmd->S_from2[0][time][i] <= PHmax))
                    bincount[b] += cmd->S_from1[0][time][i];
            lower += interval;
        }
    }

    /**print out**/
    for (int i = 0; i < bins; ++i)
        fprintf (fp, "\n%f ", bincount[i]);

    fclose (fp);

    /**gnuplot commands for a help**/
    fprintf (stderr, "\ngnuplot");
    fprintf (stderr, "\nset style data histeps");
    fprintf (stderr, "\nset xrange [-0.5:%.1f]", bins - 0.5);
    fprintf (stderr, "\nset xtics (\"%.2f\" -0.5, \"%.2f\" %f, \"%.2f\" %f)", PHmin, 0.5*(PHmin+PHmax), 0.5*bins-0.5, PHmax, bins-0.5);
    fprintf (stderr, "\nplot \"%s\" title \"phase_%c_%d\"\n", fullname, cmd->ch_target, cmd->area);

    free (bincount);
    return (DOUBLE)0;
}
