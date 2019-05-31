#include <stdio.h>
#include <stdlib.h>    /**for total_four**/
#include <math.h>
#include <string.h>    /**for phoneme2array**/
#include <time.h>      /**for total_pause**/
#include <unistd.h>    /**for total_pause; usleep**/
#include "../kernel/coco.h"
#include "../kernel/series.h"
#include "../kernel/utils.h"

/******************************* feed_reward *********************************/
/* For reinforcement learning. One unit on target area!                      */
/* Output is set to reward=1 if max of input is at position x/y=q[0][0]/[1]. */
/* q[1][0]/[1] d_a/d_b of input area (not known in A or z)!                  */
/* q[2][0] = negative "reward" value at the borders.                         */

DOUBLE feed_reward (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    DOUBLE max = 0.0;
    int X_winner = 0;
    int Y_winner = 0;
    int d_a = (int)(cmd->quantum[1][0]);
    int d_b = (int)(cmd->quantum[1][1]);

    /**get the location on the input area (winner)**/
    for (int X = 0; X < d_a; X++)
        for (int Y = 0; Y < d_b; Y++)
            if  (cmd->S_from1[0][ct_t][X * d_b + Y] > max) {

                X_winner = X;
                Y_winner = Y;
                max = cmd->S_from1[0][ct_t][X * d_b + Y];
            }

    /**if target**/
    if  ((abs (X_winner - (int)(cmd->quantum[0][0])) <= 1) && (abs (Y_winner - (int)(cmd->quantum[0][1])) <= 1)) {
        fprintf (stderr, "+");
        return (1.0);
    }

    /**if border**/
    if  ((X_winner == 0) || (X_winner == d_a - 1) || (Y_winner == 0) || (Y_winner == d_b - 1)) {
        fprintf (stderr, "-");
        return (cmd->quantum[2][0]);
    }

    /**else**/
    return (0.0);
}



/******************************* feed_angl_reward ****************************/
/* For reinforcement learning. One unit on target area!                      */
/* q[0][0]/[1] x/y position at which output is set to reward=1 if max there. */
/* q[1][0]/[1]/[2] mult_factor of angle / in_d_a_orig / in_d_b of input area.*/
/* q[2][0]/[1] padding around reward position in x and y / in angle direction*/
/* q[3][0] = negative "reward" value at the borders.                         */

DOUBLE feed_angl_reward (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    DOUBLE max       = 0.0;
    int X_winner     = 0;
    int Y_winner     = 0;
    int Ang_winner   = 0;
    int angle_factor = (int)(cmd->quantum[1][0]);
    int in_d_a_orig     = (int)(cmd->quantum[1][1]);
    int in_d_b          = (int)(cmd->quantum[1][2]);


    /**get the location on the input area (winner)**/
    for (int angle_count = 0; angle_count < angle_factor; ++angle_count) {

        int offset = angle_count * in_d_a_orig * in_d_b;

        for (int X = 0; X < in_d_a_orig; X++)
            for (int Y = 0; Y < in_d_b; Y++)
                if  (cmd->S_from1[0][ct_t][X * in_d_b + Y + offset] > max) {

                    X_winner = X;
                    Y_winner = Y;
                    Ang_winner = angle_count;
                    max = cmd->S_from1[0][ct_t][X * in_d_b + Y + offset];
                }
    }


    /**if target**/
    if  (  (abs (X_winner - (int)(cmd->quantum[0][0])) <= (int)(cmd->quantum[2][0]))
        && (abs (Y_winner - (int)(cmd->quantum[0][1])) <= (int)(cmd->quantum[2][0]))
        && (abs (Ang_winner - angle_factor / 2)        <= (int)(cmd->quantum[2][1]))) {
        fprintf (stderr, " + X=%d Y=%d Ang=%d ", X_winner, Y_winner, Ang_winner);
        return (1.0);
    }

    /**if border**/
    if  ((X_winner == 0) || (X_winner == in_d_a_orig - 1) || (Y_winner == 0) || (Y_winner == in_d_b - 1)) {
        fprintf (stderr, " - ");
        return (cmd->quantum[3][0]);
    }

    /**else**/
    return (0.0);
}



/******************************* feed_dock_reward ****************************/
/* From feed_angl_reward. Needs real (x,y,phi) coordinates from other area!  */
/* Here negative "reward", if robot is near table AND has a large angle.     */
/* q[0][0]/[1] x/y position at which output is set to reward=1 if max there. */
/* q[1][0]/[1]/[2] mult_factor of angle / in_d_a_orig / in_d_b of input area.*/
/* q[2][0]/[1] padding around reward position in x and y / in angle direction*/
/* q[3][0] = negative "reward" value at the borders.                         */
/* q[4][0] = width of the robot (relevant for colliding with the table).     */

DOUBLE feed_dock_reward (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int ct_n) {

    DOUBLE max       = 0.0;
    int X_winner     = 0;
    int Y_winner     = 0;
    int Ang_winner   = 0;
    int angle_factor = (int)(cmd->quantum[1][0]);
    int in_d_a_orig     = (int)(cmd->quantum[1][1]);
    int in_d_b          = (int)(cmd->quantum[1][2]);


    /**get the location on the input area (winner)**/
    for (int angle_count = 0; angle_count < angle_factor; ++angle_count) {

        int offset = angle_count * in_d_a_orig * in_d_b;

        for (int X = 0; X < in_d_a_orig; X++)
            for (int Y = 0; Y < in_d_b; Y++)
                if  (cmd->S_from1[0][ct_t][X * in_d_b + Y + offset] > max) {

                    X_winner = X;
                    Y_winner = Y;
                    Ang_winner = angle_count;
                    max = cmd->S_from1[0][ct_t][X * in_d_b + Y + offset];
                }
    }


    /**if target**/
    if  (  (abs (X_winner - (int)(cmd->quantum[0][0])) <= (int)(cmd->quantum[2][0]))
        && (abs (Y_winner - (int)(cmd->quantum[0][1])) <= (int)(cmd->quantum[2][0]))
        && (abs (Ang_winner - angle_factor / 2)        <= (int)(cmd->quantum[2][1]))) {
        fprintf (stderr, " + X=%d Y=%d Ang=%d ", X_winner, Y_winner, Ang_winner);
        return (1.0);
    }

    /**if border**/
    if  ((X_winner == 0) || (X_winner == in_d_a_orig - 1) || (Y_winner == 0) || (Y_winner == in_d_b - 1)) {
        fprintf (stderr, " - ");
        return (cmd->quantum[3][0]);
    }

    /**check**/
    if  (cmd->anz_quantums != 5)
        fprintf (stderr, "\n\nfeed_dock_reward wants 5 quantums!\n\n");
    if  (cmd->anz_from2 != 1)
        fprintf (stderr, "\n\nfeed_dock_reward wants 2 input areas!\n\n");
    if  (A[cmd->n_from2[0]].d_n != 3)
        fprintf (stderr, "\n\nfeed_dock_reward assumes 3 neurons in n_from2[0]!\n\n");

    /**if x-distance from table smaller than side-intrusion of turned robot**/
    if  (cmd->S_from2[0][ct_t][0] < cmd->quantum[4][0] * fabs (sin (cmd->S_from2[0][ct_t][2]))) {
        fprintf (stderr, " ¬ ");
        return (cmd->quantum[3][0]);
    }

    /**else**/
    return (0.0);
}





/****************************** total_rand_winner ****************************/
/* Determine probabilistically the winner (which will then do an action).    */
/* Prob for activity of neuron i is:  exp(2*input_i)/(sum_j exp(2*input_j)   */

DOUBLE total_rand_winner (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

  double sum = 0.0;
  double p_i = 0.0;
  double rnd = drand48();
  int area = cmd->area;

  int d_r = A[area].d_n;

  for (int i = 0; i < d_r; ++i)
      sum += exp (2.0 * cmd->S_from1[0][ct_t][i]);

  for (int i = 0; i < d_r; ++i)
      cmd->S_target[ct_t][i] = 0.0;

  for (int i = 0; i < d_r; ++i) {
      p_i += exp (2.0 * cmd->S_from1[0][ct_t][i]) / sum;

      if  (p_i > rnd) {
          cmd->S_target[ct_t][i] = 1.0;
          rnd = 1.1; /**out of reach, so the next will not be turned ON**/
      }
  }

    return (DOUBLE)(0);
}



/****************************** total_init_place *****************************/
/* Set a Gaussian to a random location, if S_from1[0][ct_t][0] is large.     */
/* Else do nothing.                                                          */

DOUBLE total_init_place (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {
    int area = cmd->area;

  if  (cmd->S_from1[0][ct_t][0] > 0.9) {

    int X_winner = (int)(drand48() * A[area].d_a);
    int Y_winner = (int)(drand48() * A[area].d_b);

    double sigma = cmd->quantum[0][0] == 0.0 ? 1.0 : cmd->quantum[0][0];
    double normfactor = 1.0; /** / (sqrt (2.0 * M_PI) * sigma); **/

    for (int X = 0; X < A[area].d_a; X++)
        for (int Y = 0; Y < A[area].d_b; Y++) {

            /**non-periodic boundary**/
            double diffA = (double)(X - X_winner);
            double diffB = (double)(Y - Y_winner);

            cmd->S_target[ct_t][X * A[area].d_b + Y]
            = normfactor * exp (-0.5 * (diffA*diffA + diffB*diffB) / (sigma*sigma));
        }
  }

    return (DOUBLE)(0);
}



/****************************** total_move_place *****************************/
/* Move the present Gaussian one pixel, determined by S_from1[][][1] (motor).*/
/* This is the "forward model"!                                              */

DOUBLE total_move_place (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    double max = 0.0;
    int X_winner = 0;
    int Y_winner = 0;
    int area = cmd->area;

    /**get the present location of the Gaussian**/
    for (int X = 0; X < A[area].d_a; X++)
        for (int Y = 0; Y < A[area].d_b; Y++)
            if  (cmd->S_from1[1][ct_t][X * A[area].d_b + Y] > max) {

                X_winner = X;
                Y_winner = Y;
                max = cmd->S_from1[1][ct_t][X * A[area].d_b + Y];
            }


    /**move one pixel**/
    if  (cmd->S_from1[0][ct_t][0] == 1.0) {
        X_winner += 1;
        if  (X_winner == A[area].d_a)
            X_winner = A[area].d_a - 1;
    }
    if  (cmd->S_from1[0][ct_t][1] == 1.0) {
        X_winner -= 1;
        if  (X_winner == -1)
            X_winner = 0;
    }
    if  (cmd->S_from1[0][ct_t][2] == 1.0) {
        Y_winner += 1;
        if  (Y_winner == A[area].d_b)
            Y_winner = A[area].d_b - 1;
    }
    if  (cmd->S_from1[0][ct_t][3] == 1.0) {
        Y_winner -= 1;
        if  (Y_winner == -1)
            Y_winner = 0;
    }

    double sigma = cmd->quantum[0][0] == 0.0 ? 1.0 : cmd->quantum[0][0];
    double normfactor = 1.0; /** / (sqrt (2.0 * M_PI) * sigma); **/

    for (int X = 0; X < A[area].d_a; X++)
        for (int Y = 0; Y < A[area].d_b; Y++) {

            /**non-periodic boundary**/
            double diffA = (double)(X - X_winner);
            double diffB = (double)(Y - Y_winner);

            cmd->S_target[ct_t][X * A[area].d_b + Y]
            = normfactor * exp (-0.5 * (diffA*diffA + diffB*diffB) / (sigma*sigma));
        }

    return (DOUBLE)(0);
}



/****************************** total_init_coord *****************************/
/* Set to a random (x,y,phi=0) position, if S_from1[0][ct_t][0] != 0.        */
/* Else do nothing.                                                          */
/* Coordinates (x,y) scale like the pixels of the image.                     */
/* q[0][0] = max x distance from front wall for init                         */
/*     [1] =  " y     "                                                      */
/*     [2] = min x, if exists!                                               */
/* q[1][0] = maxcount; init after these many calls                           */
/* q[2][0] = init angle, but:                                                */
/*     [1], if exists, then q[2][0] .. q[2][1] is interval for random angle  */

DOUBLE total_init_coord (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

  static int counter = 0;

  if  (cmd->anz_quantums != 3)
      fprintf (stderr, "\n\ntotal_init_coord wants new init angle q[2][0]!\n\n");

  if  ((cmd->S_from1[0][ct_t][0] != 0.0) || (counter == (int)(cmd->quantum[1][0]))) {

      fprintf (stderr, " total_init_coord active .. ");

      if  (cmd->anz_quant[0] == 2) {
          cmd->S_target[ct_t][0] = drand48() * cmd->quantum[0][0]; /**distance x**/
          cmd->S_target[ct_t][1] = - cmd->quantum[0][1] / 2
                                 + drand48() * cmd->quantum[0][1]; /**lateral  y**/
      } else { /**cmd->anz_quant == 3 assumed**/
          cmd->S_target[ct_t][0] = cmd->quantum[0][2] + drand48() * (cmd->quantum[0][0] - cmd->quantum[0][2]); /**distance x**/
          cmd->S_target[ct_t][1] = - cmd->quantum[0][1] / 2
                                 + drand48() * cmd->quantum[0][1]; /**lateral  y**/
      }

      if  (cmd->anz_quant[2] == 1)
          cmd->S_target[ct_t][2] = cmd->quantum[2][0];         /**angle  phi**/
      else /**cmd->anz_quant == 2 assumed**/
          cmd->S_target[ct_t][2] = cmd->quantum[2][0] + drand48() * (cmd->quantum[2][1] - cmd->quantum[2][0]);

      if  (counter == (int)(cmd->quantum[1][0]))
          fprintf (stderr, " newpos ");

      counter = 0;

      fprintf (stderr, " .. done ");
  }

  counter += 1;

    return (DOUBLE)(0);
}



/****************************** total_scan_coord *****************************/
/* For flowfield.tcl display via file: /tmp/obs_flowfield.dat                */
/* q[0][0]==1: set x,y along grid and export x,y,phi to file                 */
/*             q[1][0/1]== area dimensions                                   */
/*             q[2][0/1]== grid x/y-dimensions ("resolutions" on the grid)   */
/*                         function exits the program after full scan        */
/*             q[3][0]  == phi                                               */
/* q[0][0]==2: export output values to (same) file                           */
/* Coordinates (x,y) scale like the pixels of the image.                     */
/* Flaw: cannot test whether data point is IN or OUTside of visual field!    */

DOUBLE total_scan_coord (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

  static int x_ct = 0;
  static int y_ct = 0;
  static double sumsq_0 = 0.0;
  static double sumsq_1 = 0.0;
  static int ctr_sum = 0;
    int area = cmd->area;

  /**first time write header**/
  if  ((x_ct == 0) && (y_ct == 0)) {

      FILE *fp = fopen ("/tmp/obs_flowfield.dat", "w");                           /**write (NOT append)**/
      fprintf (fp, "%d %d\n", (int)cmd->quantum[2][0], (int)cmd->quantum[2][1]);  /**slow x-dimension, fast y-dimension**/
      fprintf (fp, "%f %f\n", cmd->quantum[1][0], cmd->quantum[1][1]);            /**x-size, y-size (note that lateral y-distance is from: -ysize/2 to: +ysize/2 !)**/
      fprintf (fp, "%f\n", cmd->quantum[3][0]);                                   /**angle**/
      fclose (fp);
  }


  /**write x,y**/
  if  (cmd->quantum[0][0] == 1) {

      if  (cmd->anz_quantums != 4)
          fprintf (stderr, "\n\ntotal_scan_coord wants 4 quantums here!\n\n");

      double x_advance = (double)x_ct / (cmd->quantum[2][0] - 1.0);  /*- 1.0, because x_ct doesn't reach upper resolution value*/
      double y_advance = (double)y_ct / (cmd->quantum[2][1] - 1.0);

      cmd->S_target[ct_t][0] = x_advance * cmd->quantum[1][0]; /**distance x**/
      cmd->S_target[ct_t][1] = - cmd->quantum[1][1] / 2
                             + y_advance * cmd->quantum[1][1]; /**lateral  y**/
      cmd->S_target[ct_t][2] = cmd->quantum[3][0];             /**angle  phi**/

      FILE *fp = fopen ("/tmp/obs_flowfield.dat", "a");
      fprintf (fp, "%f %f", cmd->S_target[ct_t][0], cmd->S_target[ct_t][1]);
      fclose (fp);

      if  (x_ct == (int)cmd->quantum[2][0]) {
          fprintf (stderr, "\ntotal_scan_coord: end of loop; exiting.\n");
          exit (0);
      }

      y_ct += 1;
      if  (y_ct == (int)(cmd->quantum[2][1])) {
          x_ct += 1;
          y_ct = 0;
      }
  }


  /**write motor output (should be 4 values)**/
  if  (cmd->quantum[0][0] == 2) {

      FILE *fp = fopen ("/tmp/obs_flowfield.dat", "a");
      for (int i = 0; i < A[area].d_n; ++i)
          fprintf (fp, " %f", cmd->S_from1[0][ct_t][i]);
      fprintf (fp, "\n");
      fclose (fp);
  }


  /**sum up the errors and write to stderr**/
  if  (cmd->quantum[0][0] == 3) {
      int index = (int)(cmd->quantum[1][0]);

      if  (index == 0)
	for (int i = 0; i < A[area].d_n; ++i) {
              sumsq_0 += fabs (cmd->S_from1[0][ct_t][i]);
              ctr_sum += 1;
        }

      if  (index == 1)
          for (int i = 0; i < A[area].d_n; ++i)
              sumsq_1 += fabs (cmd->S_from1[0][ct_t][i]);

      fprintf (stderr, "\nsumsq_0 = %f, sumsq_1 = %f, ctr_sum = %d  ", sumsq_0, sumsq_1, ctr_sum);
  }

    return (DOUBLE)(0);
}


/****************************** total_imag_coord *****************************/
/* Compute perceived Gaussian or HD angle, given (x,y,phi)-position.         */
/* This is the "sensory model"!                                              */
/* q[0][0] == 1: perceived target location; == 2: HD cell activation profile */
/*         == 3: phi relative to nearest wall                                */
/* q[1][0] == sigma of Gaussian (for both cases)                             */
/* q[2][0] == multiplicity of angle (if HD cell act)                         */
/* q[2][1] == angle assigned to outermost neurons (if HD cell act)           */

DOUBLE total_imag_coord (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    double x   = cmd->S_from1[0][ct_t][0];
    double y   = cmd->S_from1[0][ct_t][1];
    double phi = cmd->S_from1[0][ct_t][2];
    int area = cmd->area;

    if  (cmd->anz_quantums < 2)
        fprintf (stderr, "\n\ntotal_imag_coord wants 2 or more quantums!!!\n");

    double sigma = cmd->quantum[1][0] == 0.0 ? 1.0 : cmd->quantum[1][0];

    /**perceived target location**/
    if  (cmd->quantum[0][0] == 1) {

        double d     = sqrt (x*x + y*y);
        double theta =   atan (y/x) - phi;

        double a = (double)(A[area].d_a) - d * cos (theta);
        double b = (double)(A[area].d_b) / 2.0 + d * sin (theta);

        /**correction (of perception only!) at edges**/
        /**rely on negative reinforcement here to initialize at new position**/
        if  (a < 0.0)
            a = 0.0;
        if  (a > (double)(A[area].d_a - 1))
            a = (double)(A[area].d_a - 1);
        if  (b < 0.0)
            b = 0.0;
        if  (b > (double)(A[area].d_b - 1))
            b = (double)(A[area].d_b - 1);


        double normfactor = 1.0; /** / (sqrt (2.0 * M_PI) * sigma); **/

        for (int X = 0; X < A[area].d_a; X++)
            for (int Y = 0; Y < A[area].d_b; Y++) {

                /**non-periodic boundary**/
                double diffA = (double)(X) - a;
                double diffB = (double)(Y) - b;

                cmd->S_target[ct_t][X * A[area].d_b + Y]
                = normfactor * exp (-0.5 * (diffA*diffA + diffB*diffB) / (sigma*sigma));
            }
    }

    /**make phi relative to the nearest wall**/
    if  (cmd->quantum[0][0] == 3) {

        const double ywall = 28 - 1;  /*16;*/  /**if values change then**/
        const double xwall = 56 + 10; /*48;*/  /**check also below!**/

        double relphi = 999;

        /**front wall is closest**/
        if  ((x < ywall - y) && (x < y + ywall) && (x < xwall / 2)) {
            relphi = phi;
        }

        /**end wall is closest**/
        if  ((xwall - x < ywall - y) && (xwall - x < y + ywall) && (x > xwall / 2)) {
            if  (phi > 0.0)
                relphi = phi - M_PI;
            else
                relphi = phi + M_PI;
        }

        /**right wall is closest**/
        if  ((x > y + ywall) && (xwall - x > y + ywall) && (y < 0.0)) {
            if  (phi > -M_PI / 2.0)
                relphi = phi - M_PI / 2.0;          /**subtract 90^o**/
            else
                relphi = phi + M_PI * 3.0 / 2.0;    /**add 270^o avoiding underflow**/
        }

        /**left wall is closest**/
        if ((x > ywall - y) && (xwall - x > ywall - y) && (y > 0.0)) {
            if  (phi <= M_PI / 2.0)
                relphi = phi + M_PI / 2.0;          /**add 90^o**/
            else
                relphi = phi - M_PI * 3.0 / 2.0;    /**subtract 270^o avoiding overflow**/
        }

        if  (relphi == 999)
            fprintf (stderr, "\ntotal_imag_coord warning: relphi not set!\n");

        phi = relphi;
    }


    /**rotation angle**/
    if  ((cmd->quantum[0][0] == 2) || (cmd->quantum[0][0] == 3)) {

        if  (cmd->anz_quantums != 3)
            fprintf (stderr, "\n\ntotal_imag_coord wants 3 quantums for angle!\n");

        int angle_factor = (int)(cmd->quantum[2][0]);     /**number of units along the array**/
        double max_phi   = cmd->quantum[2][1];            /**max angle -- angle will be between -max_phi and +max_phi**/
        double scaledAngle = phi / max_phi;

        if  (angle_factor != A[cmd->area].d_n)
            fprintf (stderr, "\ntotal_imag_coord: wrong dimensions for angle: angle_factor=%d, d_r=%d\n", angle_factor, A[cmd->area].d_n);

        if  (scaledAngle > 1.0)
            scaledAngle = 1.0;
        if  (scaledAngle < -1.0)
            scaledAngle = -1.0;
        scaledAngle *= (double)(angle_factor / 2);

        for (int angle_count = 0; angle_count < angle_factor; ++angle_count) {

            double diffAngle = scaledAngle - (double)(angle_count - angle_factor / 2);

            cmd->S_target[ct_t][angle_count] = exp (-0.5 * diffAngle*diffAngle / (sigma*sigma));
        }
    }

    return (DOUBLE)(0);
}



/****************************** total_coord_imag *****************************/
/* "Inverse" of total_image_coord, used to visualize the mental state.       */
/* Given perceived target-where and phi activations, compute robot (x,y,phi) */
/* q[0][0] == dim a of S_from1[0] input (where area)                         */
/* q[0][1] ==  "  b           "                                              */
/* q[1][0] == multiplicity of angle == dim a*b of S_from1[1] input           */
/* q[1][1] == angle assigned to outermost neurons                            */
/* S_from1[0] is perceived target pos.;  S_from1[1] is perceived robot angle */
/* Code partly from total_angl_where.                                        */

DOUBLE total_coord_imag (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    /***first, get the angle***/
    int angle_factor = (int)(cmd->quantum[1][0]);
    double max_phi   = cmd->quantum[1][1];

    double max = 0.0;
    double p_winner = 0;
    for (int p = 0; p < angle_factor; ++p)
        if  (cmd->S_from1[1][ct_t][p] > max) {
            p_winner = p;
            max = cmd->S_from1[1][ct_t][p];
        }

    p_winner -= (double)(angle_factor) / 2.0;
    double phi = p_winner * max_phi / (double)(angle_factor) * 2.0;
    cmd->S_target[ct_t][2] = phi;


fprintf (stderr, "\ntotal_coord_imag, ct_t=%d, phi=%f, ", ct_t, cmd->S_target[ct_t][2]);


    /***then, get the robot x,y position from perceived orange and angle***/
    int d_a = (int)(cmd->quantum[0][0]);
    int d_b = (int)(cmd->quantum[0][1]);

    max = 0.0;
    int X_winner = 0;
    int Y_winner = 0;

    /**get the present location of the ?Gaussian?**/
    for (int X = 0; X < d_a; X++)
        for (int Y = 0; Y < d_b; Y++)
            if  (cmd->S_from1[0][ct_t][X * d_b + Y] > max) {
                X_winner = X;
                Y_winner = Y;
                max = cmd->S_from1[0][ct_t][X * d_b + Y];
            }

    double a = (double)(d_a - X_winner);                 /**perceived x-distance (e.g. 0 .. 16-1)      **/
    double b = (double)Y_winner - (double)(d_b) / 2.0;   /**perceived y-distance (e.g. -12 .. 11 or so)**/

fprintf (stderr, "a=%.1f, b=%.1f, ", a, b);

    if  (a == 0.0)   /**just to avoid division by zero**/
        a += 0.0001;
    double theta = atan (b / a);  /**-?**/
    double T = tan (theta + phi); /**-?**/
    double d_sq = a*a + b*b;

    cmd->S_target[ct_t][0] = sqrt (d_sq / (1.0 + T*T));   /**x-position**/
    cmd->S_target[ct_t][1] = T * cmd->S_target[ct_t][0];  /**y-position**/

fprintf (stderr, "x=%.2f, y=%.2f     ", cmd->S_target[ct_t][0], cmd->S_target[ct_t][1]);

    return (DOUBLE)(0);
}



/****************************** total_angl_coord *****************************/
/* Compute perceived Gaussian, given (x,y,phi)  -- on "angle-blown-up input. */
/* Like total_imag_coord; output area increased (~ color images) for angles. */
/* q[0][0]=sigmas for x,y pos of orange                                      */
/* q[0][1]=sigma for angle-depth                                             */
/* q[1][0]=multiplicity of angle (how much larger is d_a than actual field)  */
/* q[1][1]=angle assigned to outermost neurons (as all larger angles)        */

DOUBLE total_angl_coord (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {
    int area = cmd->area;

    double x   = cmd->S_from1[0][ct_t][0];
    double y   = cmd->S_from1[0][ct_t][1];
    double phi = cmd->S_from1[0][ct_t][2];

    double d     = sqrt (x*x + y*y);  /**distance to object**/
    double theta = atan (y/x) - phi;  /**perceived angle to object**/

    int angle_factor = (int)(cmd->quantum[1][0]);
    int d_a_orig     = A[area].d_a / angle_factor;
    double max_phi   = cmd->quantum[1][1];

    double a = (double)(d_a_orig) - d * cos (theta);
    double b = (double)(A[area].d_b) / 2.0 + d * sin (theta);

    /**correction (of perception only!) at edges**/
    /**rely on negative reinforcement here to initialize at a new position**/
    if  (a < 0.0)
        a = 0.0;
    if  (a > (double)(d_a_orig - 1))
        a = (double)(d_a_orig - 1);
    if  (b < 0.0)
        b = 0.0;
    if  (b > (double)(A[area].d_b - 1))
        b = (double)(A[area].d_b - 1);

    double sigma       = cmd->quantum[0][0] == 0.0 ? 1.0 : cmd->quantum[0][0];
    double sigma_angle = cmd->quantum[0][1] == 0.0 ? 0.00001 : cmd->quantum[0][1];

    double normfactors[101];
    double normfactor = 0.0;
    for (int angle_count = 0; angle_count < angle_factor; ++angle_count) {

        double scaledAngle = phi / max_phi;
               if  (scaledAngle > 1.0)
                   scaledAngle = 1.0;
               if  (scaledAngle < -1.0)
                   scaledAngle = -1.0;
               scaledAngle *= (double)(angle_factor / 2);

        double diffAngle = scaledAngle - (double)(angle_count - angle_factor / 2);

        normfactors[angle_count] = exp (-0.5 * diffAngle*diffAngle / (sigma_angle*sigma_angle));
        normfactor += normfactors[angle_count];
    }

    for (int angle_count = 0; angle_count < angle_factor; ++angle_count) {

        int offset = angle_count * d_a_orig * A[area].d_b;
        normfactors[angle_count] /= normfactor;

        for (int X = 0; X < d_a_orig; X++)
            for (int Y = 0; Y < A[area].d_b; Y++) {

                /**non-periodic boundary**/
                double diffA = (double)(X) - a;
                double diffB = (double)(Y) - b;

                cmd->S_target[ct_t][X * A[area].d_b + Y + offset]
                = normfactors[angle_count] * exp (-0.5 * (diffA*diffA + diffB*diffB) / (sigma*sigma));
            }
    }

    return (DOUBLE)(0);
}



/****************************** total_angl_where *****************************/
/* Compute perceived Gaussian on angle-blown-up input, given "where" and phi.*/
/* S_from1 is the seen "where" area, S_from2[][][2] delivers the angle.      */
/* (Needed for real robot to get angle-blown-up where from "where".)         */
/* q[0][0]=sigmas for x,y pos of orange                                      */
/* q[0][1]=sigma for angle-depth                                             */
/* q[1][0]=multiplicity of angle (how much larger is d_a than actual field)  */
/* q[1][1]=angle assigned to outermost neurons (as all larger angles)        */

DOUBLE total_angl_where (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {
    int area = cmd->area;

    double phi = cmd->S_from2[0][ct_t][2];

    int angle_factor = (int)(cmd->quantum[1][0]);
    int d_a_orig     = A[area].d_a / angle_factor;
    double max_phi   = cmd->quantum[1][1];

    double max = 0.0;
    int X_winner = 0;
    int Y_winner = 0;

    /**get the present location of the Gaussian**/
    for (int X = 0; X < d_a_orig; X++)
        for (int Y = 0; Y < A[area].d_b; Y++)
            if  (cmd->S_from1[0][ct_t][X * A[area].d_b + Y] > max) {

                X_winner = X;
                Y_winner = Y;
                max = cmd->S_from1[0][ct_t][X * A[area].d_b + Y];
            }

    double a = (double)X_winner;
    double b = (double)Y_winner;

    /**below unchanged from total_angl_coord**/

    /**correction (of perception only!) at edges**/
    /**rely on negative reinforcement here to initialize at a new position**/
    if  (a < 0.0)
        a = 0.0;
    if  (a > (double)(d_a_orig - 1))
        a = (double)(d_a_orig - 1);
    if  (b < 0.0)
        b = 0.0;
    if  (b > (double)(A[area].d_b - 1))
        b = (double)(A[area].d_b - 1);

    double sigma       = cmd->quantum[0][0] == 0.0 ? 1.0 : cmd->quantum[0][0];
    double sigma_angle = cmd->quantum[0][1] == 0.0 ? 0.00001 : cmd->quantum[0][1];

    double normfactors[101];
    double normfactor = 0.0;
    for (int angle_count = 0; angle_count < angle_factor; ++angle_count) {

        double scaledAngle = phi / max_phi;
               if  (scaledAngle > 1.0)
                   scaledAngle = 1.0;
               if  (scaledAngle < -1.0)
                   scaledAngle = -1.0;
               scaledAngle *= (double)(angle_factor / 2);

        double diffAngle = scaledAngle - (double)(angle_count - angle_factor / 2);

        normfactors[angle_count] = exp (-0.5 * diffAngle*diffAngle / (sigma_angle*sigma_angle));
        normfactor += normfactors[angle_count];
    }

    for (int angle_count = 0; angle_count < angle_factor; ++angle_count) {

        int offset = angle_count * d_a_orig * A[area].d_b;
        normfactors[angle_count] /= normfactor;

        for (int X = 0; X < d_a_orig; X++)
            for (int Y = 0; Y < A[area].d_b; Y++) {

                /**non-periodic boundary**/
                double diffA = (double)(X) - a;
                double diffB = (double)(Y) - b;

                cmd->S_target[ct_t][X * A[area].d_b + Y + offset]
                = normfactors[angle_count] * exp (-0.5 * (diffA*diffA + diffB*diffB) / (sigma*sigma));
            }
    }

    return (DOUBLE)(0);
}


/****************************** total_xyp_to_dtp *****************************/
/* Transforms odometry coordinates (x,y,phi) to (sqrt(diff),theta,phi).      */
/* Former given by gazebo, latter also from neural coord. transformations.   */
/* Use in conjunction with total_gauss_at, replacing total_angl/coord_where! */

DOUBLE total_xyp_to_dtp (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    double x   = cmd->S_from1[0][ct_t][0];
    double y   = cmd->S_from1[0][ct_t][1];
    double phi = cmd->S_from1[0][ct_t][2];

    double dist    = sqrt (x*x + y*y);
    double sq_dist = sqrt (dist);
    double alpha   = - atan (y/x);
    double theta   = alpha + phi;            /**body-centred angle of orange**/

    cmd->S_target[ct_t][0] = sq_dist;
    cmd->S_target[ct_t][1] = theta;
    cmd->S_target[ct_t][2] = phi;

    fprintf (stderr, "\n x=%.2f y=%.2f phi=%.2f  d=%.2f t=%.2f a=%.2f ", x, y, phi, sq_dist, theta, alpha);

    return (DOUBLE)(0);
}



/****************************** total_move_coord *****************************/
/* Output (x,y,phi) new from old (S_from1[1][][]) and motor (S_from[0][][]). */
/* This is the "motor model"!                                                */

DOUBLE total_move_coord (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

    double x   = cmd->S_from1[1][ct_t][0];
    double y   = cmd->S_from1[1][ct_t][1];
    double phi = cmd->S_from1[1][ct_t][2];

    double advance = 0.0;
    double phi_advance = 0.0;

    if  (cmd->S_from1[0][ct_t][0] == 1.0)
        advance = cmd->quantum[0][0];
    if  (cmd->S_from1[0][ct_t][1] == 1.0)
        advance = -1.0 * cmd->quantum[0][0];
    if  (cmd->S_from1[0][ct_t][2] == 1.0)
        phi_advance = cmd->quantum[0][1];
    if  (cmd->S_from1[0][ct_t][3] == 1.0)
        phi_advance = -1.0 * cmd->quantum[0][1];

    cmd->S_target[ct_t][0] = x - cos (phi) * advance;   /**distance to table**/
    cmd->S_target[ct_t][1] = y - sin (phi) * advance;   /**lateral offset   **/
    cmd->S_target[ct_t][2] = phi + phi_advance;         /**angle            **/

    if  (cmd->S_target[ct_t][2] > M_PI)
        cmd->S_target[ct_t][2] -= 2.0 * M_PI;
    if  (cmd->S_target[ct_t][2] < -M_PI)
        cmd->S_target[ct_t][2] += 2.0 * M_PI;

if  (phi == 0.0)
if  (cmd->S_from1[1][ct_t][1] != cmd->S_target[ct_t][1])
    fprintf (stderr, "\n\nlook closer at total_move_coord!\n\n");

/*
fprintf (stderr, "\ntotal_move_coord: ct_t=%d, x=%f y=%f phi=%f  ", ct_t, cmd->S_target[ct_t][0], cmd->S_target[ct_t][1], cmd->S_target[ct_t][2]);
*/

    return (DOUBLE)(0);
}
