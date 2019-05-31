#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "../kernel/coco.h"
#include "../kernel/series.h"
#include "../kernel/utils.h"



/******************************** add_a_line *********************************/
/* Used for init_lines_hier/mixed. Puts one line into a square vector.       */

void add_a_line (DOUBLE *vec, int ho_a, int br_b, int ori, int offset,
                 int linestyle, double pixel_prob, double lum) {
    int j;

    /**draw line along fast counting index **/
    if  (ori == 0)
        for (j = 0; j < br_b; ++j)
        if  (drand48() < pixel_prob) {
            vec[offset * br_b + j] += lum;
            if  (linestyle >= 2)
             vec[((offset+1) * br_b) % (ho_a*br_b) + j] -= lum;
        }

    /**draw line along slow counting index **/
    if  (ori == 1)
        for (j = 0; j < ho_a; ++j)
        if  (drand48() < pixel_prob) {
            vec[j * br_b + offset] += lum;
            if  (linestyle >= 2)
                vec[j * br_b + (offset+1) % br_b] -= lum;
        }

    if  (ori == 2)
        for (j = 0; j < ho_a; ++j)
        if  (drand48() < pixel_prob) {
            vec[(offset * br_b + j * (1 + br_b))
                % (ho_a * br_b)] += lum;
            if  (linestyle >= 2)
                vec[((offset + 1) * br_b + j * (1 + br_b))
                    % (ho_a * br_b)] -= lum;
        }

    if  (ori == 3)
        for (j = 0; j < ho_a; ++j)
        if  (drand48() < pixel_prob) {
            double sign = linestyle == 3 ? -1.0 : 1.0;
            vec[(offset * br_b + (-j + br_b) % br_b
                 + j * br_b) % (ho_a * br_b)] += lum * sign;
            if  (linestyle >= 2)
                vec[((offset + 1) * br_b + (-j + br_b) % br_b
                    + j * br_b) % (ho_a * br_b)] -= lum * sign;
        }
}


/******************************** set_lines_scaffold *************************/
/* Used for init_lines_hier/mixed. Sets array 'line_is_set' to 0's and 1's.  */

void set_lines_scaffold (int *line_is_set, double mean_anz_lines, int ho_a) {

    int offset;

    /**set the scaffold vector 'line_is_set[]'**/
    if  (mean_anz_lines == -1.0) {         /**exactly one line is on**/
        for (offset = 0; offset < ho_a; ++offset)
            line_is_set[offset] = 0;
        line_is_set[(int)(drand48() * (double)(ho_a))] = 1;
    }
    if  (mean_anz_lines == -2.0) {           /**only if linestyle==3**/
        for (offset = 0; offset < ho_a; ++offset)
            line_is_set[offset] = 0;
        line_is_set[(int)(drand48() * (double)(ho_a / 2)) * 2] = 1;
    }
    if  (mean_anz_lines >= 0.0) {             /**usually**/
        for (offset = 0; offset < ho_a; ++offset)
            line_is_set[offset]
            = (drand48() < mean_anz_lines / (double)(ho_a)) ? 1 : 0;
    }
}


/******************************** init_lines_hier ****************************/
/* q[0][0] (cmd->offset) != 0.0 -> get new orientations                      */
/* q[1][0] (aux1) = no. of ori's                                             */
/* q[2][0] (aux2) = avrg no. of ori's; if -1 => exactly one ori is ON.       */
/* q[3][0] (aux3) = avrg no. of lines; if -1 => exactly one line is ON.      */
/* q[4][0] (aux4) =1: positive lines; =2: double +1/-1 line => mean zero     */
/*                =3: double +1/-1 line on every second pixel <- even pixels!*/
/*                !! double-lines are +1/-1 but not -1/+1 !!                 */
/* q[5][0] pixel probability                                                 */
/* q[6][0] fluctuation                                                       */
/* q[7][0] =0: quadratic input; =1: teacher values                           */

DOUBLE init_lines_hier (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {

    int ori, offset = 0, j, i, ct_t;
    DOUBLE *vec;
    double lum = 1.0;  /**brightness of the lines**/
    double rd;

    AREA *z = A + cmd->area;

    static int firsttime = 1;
    static int *ori_is_set;
    static int *line_is_set;
    static int *line_is_new;

    const int    dummy          = (cmd->anz_quantums == 8) ? 0 :
                 fprintf (stderr, "\n8 parameters for init_lines_hier, please! (second-last is new: fluctuation!) You have given %d params!\n", cmd->anz_quantums);

    const int    get_new_ori    =      cmd->quantum[0][0] == 0.0 ? 0 : 1;
    const int    max_anz_ori    = (int)cmd->quantum[1][0];
    const double mean_anz_ori   =      cmd->quantum[2][0];
          double mean_anz_lines =      cmd->quantum[3][0];
    const int    linestyle      = (int)cmd->quantum[4][0];
    const double pixel_prob     =      cmd->quantum[5][0];
    const double fluctuation    =      cmd->quantum[6][0];
    const int    teacher        = (int)cmd->quantum[7][0];

    int ho_a = z->d_a;
    int br_b = z->d_b;

    /**allocate mem for target values**/
    if  (firsttime) {
        ori_is_set  = i_vector (max_anz_ori);
        line_is_set = i_vector (ho_a >= br_b ? ho_a : br_b);
        line_is_new = i_vector (ho_a >= br_b ? ho_a : br_b);
        firsttime = 0;
    }

    /**get new orientation(s)**/
    if  (get_new_ori) {

        if  (mean_anz_ori == -1.0) {

            rd = drand48();

            /**one and only one orientation is on each time**/
            for (ori = 0; ori < max_anz_ori; ++ori) {
                ori_is_set[ori] = 0;
                if  (  (rd > (double)(ori)     / (double)max_anz_ori)
                    && (rd < (double)(ori + 1) / (double)max_anz_ori))
                    ori_is_set[ori] = 1;
            }
        } else {
            /**choose orientation (one of max_anz_ori = 2, 3 or 4)**/
            for (ori = 0; ori < max_anz_ori; ++ori)

                if  (drand48() < mean_anz_ori / (double)max_anz_ori)                /**each ori with prob mean_anz_ori/max_anz_ori (mean~1..<2)**/
                    ori_is_set[ori] = 1;
                else
                    ori_is_set[ori] = 0;
        }
    }



    /**paint the image**/
    if  (teacher == 0) {

        if  (ho_a != br_b)
            fprintf (stderr,"\nwrong input sizes 0 for init_lines!%d\n",dummy);

        if  (linestyle == 3)
            if  ((ho_a % 2 != 0) || (br_b % 2 != 0))
                fprintf (stderr, "\ninit_lines_hier wants even size here\n");

        if  (linestyle == 3)
            mean_anz_lines *= 2.0;     /**because lines only every 2nd pixel**/


        ct_t = begin;

            vec = cmd->S_target[ct_t];

            /**init with 0**/
            for (j = 0; j < ho_a * br_b; ++j)
                vec[j] = 0.0;

            for (ori = 0; ori < max_anz_ori; ++ori)
            if  (ori_is_set[ori]) {

        for (ct_t = begin; ct_t < end; ++ct_t) {

                if  (ct_t == begin)
                    set_lines_scaffold (line_is_set, mean_anz_lines, ho_a);

                /**slowly changing line scaffold**/
                if  (fluctuation != 0.0) {
                    set_lines_scaffold (line_is_new, mean_anz_lines, ho_a);

                    for (offset = 0; offset < ho_a; ++offset)
                        if  (drand48 () < fluctuation)
                            line_is_set[offset] = line_is_new[offset];
                }

                /**choose offset**/
                offset = 0;
                while (offset < ho_a) {

                    /**each offset with prob mean_anz_lines/Ho_a**/
                    if  (line_is_set[offset])

                        add_a_line (vec, ho_a, br_b, ori, offset, linestyle, pixel_prob, lum);

                    offset += 1;
                    if  (linestyle == 3)
                        offset += 1;
                }
            }
        }

	/**I pasted this in later, because image wasn't constant!!!**/
        /**copy to all times**/
        for (ct_t = begin; ct_t < end; ++ct_t)
            for (j = 0; j < ho_a * br_b; ++j)
                cmd->S_target[ct_t][j] = vec[j];
    }


    /**fill teacher values**/
    if  (teacher == 1) {

        vec = cmd->S_target[begin];

        if  (br_b != max_anz_ori)
            fprintf (stderr, "\nwrong input sizes 1 for init_lines!\n");

        for (j = 0; j < br_b; ++j)
            for (i = 0; i < ho_a; ++i) /**fill redundantly along ho_a**/
                if  (ori_is_set[j])
                    vec[br_b * i + j] = lum;
                else
                    vec[br_b * i + j] = 0.0;

        /**copy to all times**/
        for (ct_t = begin + 1; ct_t < end; ++ct_t)
            for (j = 0; j < ho_a * br_b; ++j)
                cmd->S_target[ct_t][j] = vec[j];
    }

    return (DOUBLE)0;
}






/******************************** data_gauss_move ****************************/
/* Places a Gaussian randomly on either upper or lower half and moves it.    */
/* Gaussian is Gaussian only along d_b, but horizontal along d_a.            */
/* q[0][0] = 0 .. 1 prob to choose a new position at the beginning           */
/* q[1][0] = sigma of Gaussian                                               */
/* q[2][0] = height of Gaussian                                              */
/* q[3][0]/[1] = velocity in upper/lower half (move-dist at each time step)  */

DOUBLE data_gauss_move (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {
     int X, Y, ct_t;
     static double cm_x = 0.0;
     static double cm_y = 0.0;
     double do_new, sigma, height, vel_upper, vel_lower;

     AREA *z = A + cmd->area;

     /**warnings of wrong usage**/
     if  (cmd->anz_quantums != 4)
         fprintf (stderr, "\nwrong number of quantums in data_gauss_move\n");
     if  (cmd->anz_quant[3] != 2)
         fprintf (stderr, "\nwrong number of quant's in data_gauss_move\n");

     do_new    = cmd->quantum[0][0];
     sigma     = cmd->quantum[1][0];
     height    = cmd->quantum[2][0];
     vel_upper = cmd->quantum[3][0];
     vel_lower = cmd->quantum[3][1];

     if  (do_new > drand48()) {
         cm_x = drand48() * (double)(z->d_a);  /**vertical starting position**/
         cm_y = drand48() * (double)(z->d_b); /*horizontal starting position**/
     }

     for (ct_t = begin; ct_t < end; ++ct_t) {

         if  (cm_x >= z->d_a / 2)                           /**in lower half**/
             cm_y += vel_lower;                                      /**move**/

         if  (cm_x < z->d_a / 2)                            /**in upper half**/
             cm_y += vel_upper;                                      /**move**/

         if  (cm_y >= z->d_b)                    /**if out of right boundary**/
             cm_y -= (double)z->d_b;                     /**start left again**/
         if  (cm_y < 0)                           /**if out of left boundary**/
             cm_y += (double)z->d_b;                    /**start right again**/

         for (X = 0; X < z->d_a; ++X)
             for (Y = 0; Y < z->d_b; ++Y) {
               /*double diffA = X - cm_x;*/
                 double diffB = Y - cm_y;

                 /**for periodic boundary**/
                 diffB = (fabs(diffB) <= fabs(z->d_b / 2.0)) ? diffB : - z->d_b + fabs(diffB);

                 if  (cm_x >= z->d_a / 2) {
                     if  (X >= z->d_a / 2)
                         cmd->S_target[ct_t][X * z->d_b + Y] = height * exp (-0.5 * (/*diffA*diffA+*/diffB*diffB)/(sigma*sigma));
                     else
                         cmd->S_target[ct_t][X * z->d_b + Y] = 0.0;
                 }

                 if  (cm_x < z->d_a / 2) {
                     if  (X < z->d_a / 2)
                         cmd->S_target[ct_t][X * z->d_b + Y] = height * exp (-0.5 * (/*diffA*diffA+*/diffB*diffB)/(sigma*sigma));
                     else
                         cmd->S_target[ct_t][X * z->d_b + Y] = 0.0;
                 }
             }
     }

     return (DOUBLE)(0);
}


/******************************** data_gauss_three ***************************/
/* Places three Gaussian, 1st & 2nd randomly and for 3rd: mu_3 = mu_1 + mu_2.*/
/* d_a must be divideable by three.                                          */
/* Gaussians are Gaussians only along d_b, but horizontal along d_a.         */
/* q[0][0] = sigma of Gaussians                                              */
/* q[1][0] = height of Gaussians                                             */
/* q[2][0] = row to set to 0; =0:none; =1 1st; =2 2nd; =3 3rd. (LATER!)      */

DOUBLE data_gauss_three (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {
     int X, B, ct_t;
     static double cm_x = 0.0; /**first Gaussian**/
     static double cm_y = 0.0; /**second Gaussian**/
     static double cm_z = 0.0; /**third Gaussian**/
     double sigma, height;

     AREA *z = A + cmd->area;

     /**warnings of wrong usage**/
     if  (cmd->anz_quantums != 3)
         fprintf (stderr, "\nwrong number of quantums in data_gauss_three\n");

     sigma     = cmd->quantum[0][0];
     height    = cmd->quantum[1][0];

     cm_x = drand48() * (double)(z->d_b / 2);
     cm_y = drand48() * (double)(z->d_b / 2);
     cm_z = cm_x + cm_y;

     for (ct_t = begin; ct_t < end; ++ct_t) {

         if  (cm_z >= z->d_b)                    /**if out of right boundary**/
             cm_z -= (double)z->d_b;                     /**start left again**/
         if  (cm_z < 0)                           /**if out of left boundary**/
             cm_z += (double)z->d_b;                    /**start right again**/

         for (X = 0; X < z->d_a; ++X)
             for (B = 0; B < z->d_b; ++B) {

                 double diffx = B - cm_x;
                 double diffy = B - cm_y;
                 double diffz = B - cm_z;

                 /**for periodic boundary**/
                 diffx = (fabs(diffx) <= fabs(z->d_b / 2.0)) ? diffx : - z->d_b + fabs(diffx);
                 diffy = (fabs(diffy) <= fabs(z->d_b / 2.0)) ? diffy : - z->d_b + fabs(diffy);
                 diffz = (fabs(diffz) <= fabs(z->d_b / 2.0)) ? diffz : - z->d_b + fabs(diffz);

                 if  (X < z->d_a / 3) {
                     if  (cmd->quantum[2][0] != 1)
                         cmd->S_target[ct_t][X * z->d_b + B] = height * exp (-0.5 * (diffx*diffx)/(sigma*sigma));
                     else
                         cmd->S_target[ct_t][X * z->d_b + B] = 0.0;
                 } else {
                     if  (X < z->d_a * 2 / 3) {
                         if  (cmd->quantum[2][0] != 2)
                             cmd->S_target[ct_t][X * z->d_b + B] = height * exp (-0.5 * (diffy*diffy)/(sigma*sigma));
                         else
                             cmd->S_target[ct_t][X * z->d_b + B] = 0.0;
                     } else {
                         if  (cmd->quantum[2][0] != 3)
                             cmd->S_target[ct_t][X * z->d_b + B] = height * exp (-0.5 * (diffz*diffz)/(sigma*sigma));
                         else
                             cmd->S_target[ct_t][X * z->d_b + B] = 0.0;
                     }
                 }
             }
     }

     return (DOUBLE)(0);
}


/******************************** data_gauss_3areas **************************/
/* Places 3 Gaussians. 0,1 randomly and for no.2: mu_2 = (mu_0+mu_1)*q[4][0] */
/* For 3 different areas.                                                    */
/* Gaussians are Gaussians only along d_b, but horizontal along d_a.         */
/* q[0][0] = initialise new if 1                                             */
/* q[1][0] = area 0, 1 or the "sum"-area 2                                   */
/* q[2][0] = sigma of Gaussian                                               */
/* q[3][0] = height of Gaussian                                              */
/* q[4][0] = should be 0.5 so that mu_2 doesn't "leave" the area boundary;   */
/*           if set to 1 then mu_2 folds, i.e. starts at left of area again. */

DOUBLE data_gauss_3areas (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {
     int X, B, ct_t;
     static double mu_0_x = 0.0; /**first Gaussian**/
     static double mu_1_x = 0.0; /**second Gaussian**/
     static double mu_2_x = 0.0; /**third Gaussian**/
     double cm_x = 0.0, sigma, height;

     AREA *z = A + cmd->area;

     /**warnings of wrong usage**/
     if  (cmd->anz_quantums < 4)
         fprintf (stderr, "\nwrong number of quantums in data_gauss_3areas\n");
     if  ((cmd->quantum[0][0] == 1) && (cmd->anz_quantums < 5))
         fprintf (stderr, "\nwrong number of quantums in data_gauss_3areas\n");

     sigma     = cmd->quantum[2][0];
     height    = cmd->quantum[3][0];

     if  (cmd->quantum[0][0] == 1) {
         mu_0_x = drand48();
         mu_1_x = drand48();
         mu_2_x = (mu_0_x + mu_1_x) * cmd->quantum[4][0];
         if  (mu_2_x > 1.0)                    /**if out of right boundary**/
             mu_2_x -= 1.0;                            /**start left again**/
     }

     if  (cmd->quantum[1][0] == 0)
         cm_x = mu_0_x * z->d_b; 
     if  (cmd->quantum[1][0] == 1)
         cm_x = mu_1_x * z->d_b; 
     if  (cmd->quantum[1][0] == 2)
         cm_x = mu_2_x * z->d_b; 

     for (ct_t = begin; ct_t < end; ++ct_t)
         for (X = 0; X < z->d_a; ++X)
             for (B = 0; B < z->d_b; ++B) {

                 double diffx = B - cm_x;

                 /**for periodic boundary**/
                 diffx = (fabs(diffx) <= fabs(z->d_b / 2.0)) ? diffx : - z->d_b + fabs(diffx);

                 cmd->S_target[ct_t][X * z->d_b + B] = height * exp (-0.5 * (diffx*diffx)/(sigma*sigma));
             }

     return (DOUBLE)(0);
}


/******************************** data_gauss_3areas_2D ***********************/
/* Places 3 Gaussians. 0,1 randomly and for no.2: mu_2 = (mu_0+mu_1)*q[4][0] */
/* For 3 different areas.                                                    */
/* Gaussians are Gaussians only along d_b, but horizontal along d_a.         */
/* NEW: avoids outer rim of units for placement of Gaussian centres !!!      */
/* q[0][0] =1: init all new; =2: init only areas 0 and 1 new (other "view"). */
/* q[1][0] = area 0, 1 or the "sum"-area 2                                   */
/* q[2][0] = sigma of Gauss; if q[2]1] exists then sigma=rand between these  */
/* q[3][0] = height of Gaussian                                              */
/* q[4][0] = mode: 1 or 2.                                                   */
/* q[5][0] = should be 0.5 so that mu_2 doesn't "leave" the area boundary;   */
/*           if set to 1 then mu_2 folds, i.e. starts at left of area again. */
/* q[6][0] = edge (distance to border) as variable away_from_boundary        */

DOUBLE data_gauss_3areas_2D (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {
     int X, B, ct_t;
     static double mu_0_x[100]; /**first Gaussian**/
     static double mu_0_y[100];
     static double mu_1_x[100]; /**second Gaussian**/
     static double mu_1_y[100];
     static double mu_2_x[100]; /**third Gaussian**/
     static double mu_2_y[100];
     double cm_x = 0.0, cm_y = 0.0;

     if  (end >= 100)
         fprintf (stderr, "\n\ndata_gauss_3areas_2D: please enlarge me for longer relax-times!\n\n");

     int periodic_boundary = 0;
     float away_from_boundary = 0.0;

     if  (cmd->anz_quantums == 7)
         away_from_boundary = cmd->quantum[6][0];   /**how far to place Gaussian centres away from outermost units**/

     AREA *z = A + cmd->area;

     /**warnings of wrong usage**/
     if  (cmd->anz_quantums < 4)
         fprintf (stderr, "\nwrong number of quantums in data_gauss_3areas; new is the distance to boundary\n");
     if  ((cmd->quantum[0][0] == 1) && (cmd->anz_quantums < 6))
         fprintf (stderr, "\nwrong number of quantums in data_gauss_3areas\n");

     double sigma     = cmd->quantum[2][0];

     if  (cmd->anz_quant[2] == 2)
         sigma = cmd->quantum[2][0] + drand48() * (cmd->quantum[2][1] - cmd->quantum[2][0]);

     double height    = cmd->quantum[3][0];

for (ct_t = begin; ct_t < end; ++ct_t) {

     if  (cmd->quantum[0][0] == 1) {

         int mode = (int)(cmd->quantum[4][0]);

         /**in this mode, determine first mu_0 and mu_1 and then mu_2; only this mode has a variable parameter q[5][0]**/
         if  (mode == 1) {
             mu_0_x[ct_t] = drand48();
             mu_0_y[ct_t] = drand48();
             mu_1_x[ct_t] = drand48();
             mu_1_y[ct_t] = drand48();

             mu_2_x[ct_t] = (mu_0_x[ct_t] + mu_1_x[ct_t]) * cmd->quantum[5][0];
             mu_2_y[ct_t] = (mu_0_y[ct_t] + mu_1_y[ct_t]) * cmd->quantum[5][0];

             if  (mu_2_x[ct_t] > 1.0)                      /**if out of right boundary**/
                 mu_2_x[ct_t] -= 1.0;                              /**start left again**/
             if  (mu_2_y[ct_t] > 1.0)                      /**if out of lower boundary**/
                 mu_2_y[ct_t] -= 1.0;                            /**start on top again**/

             /**produce positions shifted systematically, for an "anim", but within one rlen, i.e. displayable on one pnm**
             int num_blocks = 3;
             mu_0_y[ct_t] = (ct_t * num_blocks / (end - begin)) / (double)(num_blocks - 1);
             mu_1_y[ct_t] = (ct_t * num_blocks) % (end - begin) / (double)(end - begin);
             **/
         }

         /**in this mode, determine first mu_2 and then mu_1 and mu_0 in a way which is possible**/
         if  (mode == 2) {

#define DATA_GNU_CONNECT_MAX_POINTS_PER_NODE 5
#define DATA_GNU_CONNECT_MAX_COLS_AND_LINES 10

/** for test of precision / variation **
static int count_for_data_gnu_connect_max = -1;
count_for_data_gnu_connect_max += 1;
if  (count_for_data_gnu_connect_max % DATA_GNU_CONNECT_MAX_POINTS_PER_NODE == 0)
**/

             if  (ct_t == begin) {

                 mu_2_x[ct_t] = drand48();
                 mu_2_y[ct_t] = drand48();

/** for test of precision / variation **
static int count_shift = -1;
count_shift += 1;
float border = 0.05;
mu_2_x[ct_t] = border + (count_shift % DATA_GNU_CONNECT_MAX_COLS_AND_LINES) * (1.0 - 2.0 * border) / (DATA_GNU_CONNECT_MAX_COLS_AND_LINES - 1.0);
mu_2_y[ct_t] = border + (count_shift / DATA_GNU_CONNECT_MAX_COLS_AND_LINES) * (1.0 - 2.0 * border) / (DATA_GNU_CONNECT_MAX_COLS_AND_LINES - 1.0);
**/

                 /**produce sloped distribution (until "other" < "half" means that "other" will be preferentially small)**
                 double half = drand48();
                 double other;
                 do {
                     other = drand48();
                     mu_2_x[ct_t] = other;
                     mu_2_y[ct_t] = other;
                 } while (other > half);
                 **/

             } else {
                 mu_2_x[ct_t] = mu_2_x[ct_t - 1];
                 mu_2_y[ct_t] = mu_2_y[ct_t - 1];
             }

             /**produce systematically shifted positions in z (= mu_2) but random x and y (mu_0 and mu_1)**
             int tmpz = ct_t / 4;
             mu_2_x[ct_t] = drand48();
             mu_2_y[ct_t] = (double)tmpz / (double)(A[cmd->area].d_a * A[cmd->area].d_b);
             if  (mu_2_y[ct_t] <= 0.0)
                 mu_2_y[ct_t] = 0.000001;
             if  (mu_2_y[ct_t] >= 1.0)
                 mu_2_y[ct_t] = 0.999999;
             **/

             if  (drand48() < 0.5) {
                 do {
                    mu_0_x[ct_t] = drand48();
                    mu_1_x[ct_t] = 2.0 * mu_2_x[ct_t] - mu_0_x[ct_t];
                    /**produce: x + y^2 = z -- see also below again**
                    if  (mu_1_x[ct_t] >= 0)
                        mu_1_x[ct_t] = sqrt (mu_1_x[ct_t]);
                    **/
                 } while ((mu_1_x[ct_t] < 0.0) || (mu_1_x[ct_t] > 1.0));
             } else {
                 do {
                    mu_1_x[ct_t] = drand48();
                    mu_0_x[ct_t] = 2.0 * mu_2_x[ct_t] - mu_1_x[ct_t];
                    /**produce: x + y^2 = z -- see also below again**
                    if  (mu_1_x[ct_t] >= 0)
                        mu_1_x[ct_t] = sqrt (mu_1_x[ct_t]);
                    **/
                 } while ((mu_0_x[ct_t] < 0.0) || (mu_0_x[ct_t] > 1.0));
             }

             if  (drand48() < 0.5) {
                 do {
                    mu_0_y[ct_t] = drand48();
                    mu_1_y[ct_t] = 2.0 * mu_2_y[ct_t] - mu_0_y[ct_t];
                    /**produce: x + y^2 = z -- see also above again**
                    if  (mu_1_y[ct_t] >= 0)
                        mu_1_y[ct_t] = sqrt (mu_1_y[ct_t]);
                    **/
                 } while ((mu_1_y[ct_t] < 0.0) || (mu_1_y[ct_t] > 1.0));
             } else {
                 do {
                    mu_1_y[ct_t] = drand48();
                    mu_0_y[ct_t] = 2.0 * mu_2_y[ct_t] - mu_1_y[ct_t];
                    /**produce: x + y^2 = z -- see also above again**
                    if  (mu_1_y[ct_t] >= 0)
                        mu_1_y[ct_t] = sqrt (mu_1_y[ct_t]);
                    **/
                 } while ((mu_0_y[ct_t] < 0.0) || (mu_0_y[ct_t] > 1.0));
             }
         }
     }

     /**in this case, keep mu_2, and determine another mu_1 and mu_0 in a way which is possible (for repeated collection of the same data from different "viewing angle")**/
     if  (cmd->quantum[0][0] == 2) {

             if  (drand48() < 0.5) {
                 do {
                    mu_0_x[ct_t] = drand48();
                    mu_1_x[ct_t] = 2.0 * mu_2_x[ct_t] - mu_0_x[ct_t];
                 } while ((mu_1_x[ct_t] < 0.0) || (mu_1_x[ct_t] > 1.0));
             } else {
                 do {
                    mu_1_x[ct_t] = drand48();
                    mu_0_x[ct_t] = 2.0 * mu_2_x[ct_t] - mu_1_x[ct_t];
                 } while ((mu_0_x[ct_t] < 0.0) || (mu_0_x[ct_t] > 1.0));
             }

             if  (drand48() < 0.5) {
                 do {
                    mu_0_y[ct_t] = drand48();
                    mu_1_y[ct_t] = 2.0 * mu_2_y[ct_t] - mu_0_y[ct_t];
                 } while ((mu_1_y[ct_t] < 0.0) || (mu_1_y[ct_t] > 1.0));
             } else {
                 do {
                    mu_1_y[ct_t] = drand48();
                    mu_0_y[ct_t] = 2.0 * mu_2_y[ct_t] - mu_1_y[ct_t];
                 } while ((mu_0_y[ct_t] < 0.0) || (mu_0_y[ct_t] > 1.0));
             }
     }


     if  (cmd->quantum[1][0] == 0) {
         cm_x = away_from_boundary + mu_0_x[ct_t] * (z->d_b - 1 - 2*away_from_boundary);
         cm_y = away_from_boundary + mu_0_y[ct_t] * (z->d_a - 1 - 2*away_from_boundary);
     }
     if  (cmd->quantum[1][0] == 1) {
         cm_x = away_from_boundary + mu_1_x[ct_t] * (z->d_b - 1 - 2*away_from_boundary); 
         cm_y = away_from_boundary + mu_1_y[ct_t] * (z->d_a - 1 - 2*away_from_boundary);
     }
     if  (cmd->quantum[1][0] == 2) {
         cm_x = away_from_boundary + mu_2_x[ct_t] * (z->d_b - 1 - 2*away_from_boundary); 
         cm_y = away_from_boundary + mu_2_y[ct_t] * (z->d_a - 1 - 2*away_from_boundary);
     }

         for (X = 0; X < z->d_a; ++X)
             for (B = 0; B < z->d_b; ++B) {

                 double diffx = B - cm_x;
                 double diffy = X - cm_y;

                 /**for periodic boundary**/
                 if  (periodic_boundary) {
                     diffx = (fabs(diffx) <= fabs(z->d_b / 2.0)) ? diffx : - z->d_b + fabs(diffx);
                     diffy = (fabs(diffy) <= fabs(z->d_b / 2.0)) ? diffx : - z->d_b + fabs(diffx);
                 }

                 cmd->S_target[ct_t][X * z->d_b + B] = height * exp (-0.5 * (diffx*diffx+diffy*diffy)/(sigma*sigma));
             }
}
     return (DOUBLE)(0);
}


/******************************** data_gauss_fromfile ************************/
/* Places 2 Gaussians, centres from data read from file by read_act_file.    */
/* For 1D or 2D-areas, dependent anz_quant[1].                               */
/* For 2 different time steps ct_t=0 and ct_t=1.                             */
/* q[0][0] = 1: select new data point pair.                                  */
/* q[1][0(/1)] = column(s) from data file; if 2D area, 2 columns to be given */
/* q[2][0(/1)] = lower bound(s) of data value(s)                             */
/* q[3][0(/1)] = upper    "            "                                     */
/* q[4][0] = sigma of Gaussian                                               */

DOUBLE data_gauss_fromfile (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {
     int X, B, ct_t;
     int away_from_boundary = 0;   /**how far to place Gaussian centres away from outermost units**/
     static int sel[2];

     /**warnings of wrong usage**/
     if  (cmd->anz_quantums != 5)
         fprintf (stderr, "\nwrong number of quantums in data_gauss_fromfile\n");

     if  ((begin > 0) || (end > 2))
         fprintf (stderr, "\n\ndata_gauss_fromfile: relax-time should be [begin=0 .. end=2[   or just one step!\n\n");

     float ** data = (float **)(cmd->pointers[0]->data);
     int anz = (int)data[0][0];

     /*
     static int firsttime = 1;
     if  (firsttime) {
         fprintf (stderr, "\ndata_gauss_fromfile has %d data points\n", anz);
         for (int i = 1; i  <= anz; ++i)
             fprintf (stderr, "\n%d %f %f %f %f ", (int)data[i][0], data[i][1], data[i][2], data[i][3], data[i][4]);
         fprintf (stderr, "\ndata_gauss_fromfile has %d data points\n", anz);
         firsttime = 0;
     }
     */

     AREA *z = A + cmd->area;

     double sigma  = cmd->quantum[4][0];
     double height = 1.0;

     if  ((int)(cmd->quantum[0][0]) == 1) {

         /**select first data point sel1 for ct_t=0**/
         int num[2];
         ct_t = 0;
         sel[ct_t] = (int)(drand48() * (double)(anz)) + 1;  /**the sel'th data point**/
         num[ct_t] = (int)(data[sel[ct_t]][0]);                  /**diff robot poses (not diff views) have diff num's**/

         /**select second data point sel2 for ct_t=1**/
         ct_t = 1;
         do {
             sel[ct_t] = (int)(drand48() * (double)(anz)) + 1;
             num[ct_t] = (int)(data[sel[ct_t]][0]);
         } while (num[1] != num[0]); /**stop if num1==num2 (demanding sel1!=sel2 omitted as maybe just one data)**/
     }


     for (ct_t = begin; ct_t < end; ++ct_t) {

         int   column = (int)(cmd->quantum[1][0]);
         float lower  = cmd->quantum[2][0];
         float upper  = cmd->quantum[3][0];
         float mu     = data[sel[ct_t]][column];
         float mu_aux = (mu - lower) / (upper - lower);

         if  ((mu_aux < 0.0) || (mu_aux > 1.0)) {
             fprintf (stderr, "\n\nmu_aux=%f exceeds boundaries:\n", mu_aux);
             fprintf (stderr, "sel[%d]=%d, column=%d, data=%f, lower=%f, upper=%f\n",
                               ct_t, sel[ct_t], column, mu, lower, upper);
         }

         /**1D case**/
         if  (cmd->anz_quant[1] == 1) {

             if  ((A[cmd->area].d_a != 1) && (A[cmd->area].d_b != 1))
                 fprintf (stderr, "\ndata_gauss_fromfile: 1D case inconsistent!\n");

             float cm = away_from_boundary + mu_aux * (z->d_a*z->d_b - 1 - 2*away_from_boundary); 

             for (X = 0; X < z->d_a * z->d_b; ++X) {
                 float diff = X - cm;
                 cmd->S_target[ct_t][X] = height * exp (-0.5 * (diff*diff)/(sigma*sigma));
             }

         } else { /**2D case**/

             if  ((A[cmd->area].d_a == 1) || (A[cmd->area].d_b == 1))
                 fprintf (stderr, "\ndata_gauss_fromfile: 2D case inconsistent!\n");

             int   column2 = (int)(cmd->quantum[1][1]);
             float lower2  = cmd->quantum[2][1];
             float upper2  = cmd->quantum[3][1];
             float mu2     = data[sel[ct_t]][column2];
             float mu_aux2 = (mu2 - lower2) / (upper2 - lower2);

             if  ((mu_aux2 < 0.0) || (mu_aux2 > 1.0))
                 fprintf (stderr, "\n\nmu_aux2=%f exceeds boundaries!\n\n", mu_aux2);

             float cm_x = away_from_boundary + mu_aux  * (z->d_b - 1 - 2*away_from_boundary); 
             float cm_y = away_from_boundary + mu_aux2 * (z->d_a - 1 - 2*away_from_boundary);

             for (X = 0; X < z->d_a; ++X)
                 for (B = 0; B < z->d_b; ++B) {

                     float diffx = B - cm_x;
                     float diffy = X - cm_y;

                     cmd->S_target[ct_t][X * z->d_b + B] = height * exp (-0.5 * (diffx*diffx+diffy*diffy)/(sigma*sigma));
                 }
         }
     }

     return (DOUBLE)(0);
}



/******************************** data_gauss_3areas_2D_anim ******************/
/* Places 3 Gaussians. 0,1 randomly and for no.2: mu_2 = (mu_0+mu_1)*q[4][0] */
/* For 3 different areas.                                                    */
/* Gaussians are Gaussians only along d_b, but horizontal along d_a.         */
/* NEW: avoids outer rim of units for placement of Gaussian centres !!!      */
/* q[0][0] = initialise new if 1                                             */
/* q[1][0] = area 0, 1 or the "sum"-area 2                                   */
/* q[2][0] = sigma of Gaussian                                               */
/* q[3][0] = height of Gaussian                                              */
/* q[4][0] = mode: 1 or 2.                                                   */
/* q[5][0] = should be 0.5 so that mu_2 doesn't "leave" the area boundary;   */
/*           if set to 1 then mu_2 folds, i.e. starts at left of area again. */

DOUBLE data_gauss_3areas_2D_anim (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {
     int X, B, ct_t;
     static double mu_0_x = 0.0; /**first Gaussian**/
     static double mu_0_y = 0.0;
     static double mu_1_x = 0.0; /**second Gaussian**/
     static double mu_1_y = 0.0;
     static double mu_2_x = 0.0; /**third Gaussian**/
     static double mu_2_y = 0.0;
     double cm_x = 0.0, cm_y = 0.0;

     int periodic_boundary = 0;
     int away_from_boundary = 1;   /**how far to place Gaussian centres away from outermost units**/

     AREA *z = A + cmd->area;

     /**warnings of wrong usage**/
     if  (cmd->anz_quantums < 4)
         fprintf (stderr, "\nwrong number of quantums in data_gauss_3areas\n");
     if  ((cmd->quantum[0][0] == 1) && (cmd->anz_quantums < 6))
         fprintf (stderr, "\nwrong number of quantums in data_gauss_3areas\n");

     double sigma     = cmd->quantum[2][0];
     double height    = cmd->quantum[3][0];

     int numPoints = 14;
     static int time = 0;
     int smalltime = time % numPoints;
     static int smallcounter = 0;
     int bigtime = (time / 4) % numPoints;
     static int bigcounter = 0;

     if  (cmd->quantum[0][0] == 1) {

         int mode = (int)(cmd->quantum[4][0]);

         /**in this mode, determine first mu_0 and mu_1 and then mu_2; only this mode has a variable parameter q[5][0]**/
         if  (mode == 1) {

             /**draw diagonal 1**/
             if  (smallcounter == 0) {
                 mu_1_x = (double)smalltime / numPoints;
                 mu_1_y = (double)smalltime / numPoints;
             }

             /**draw opposite 1**/
             if  (smallcounter == 1) {
                 mu_1_x = 1.0 - (double)smalltime / numPoints;
                 mu_1_y = 1.0;
             }

             /**draw bottom 1**/
             if  (smallcounter == 2) {
                 mu_1_x = 0.0;
                 mu_1_y = 1.0 - (double)smalltime / numPoints;
             }

             /**draw diagonal 0**/
             if  (bigcounter == 0) {
                 mu_0_x = (double)bigtime / numPoints;
                 mu_0_y = (double)bigtime / numPoints;
             }

             /**draw opposite 0**/
             if  (bigcounter == 1) {
                 mu_0_x = 1.0 - (double)bigtime / numPoints;
                 mu_0_y = 1.0;
             }

             /**draw bottom 0**/
             if  (bigcounter == 2) {
                 mu_0_x = 0.0;
                 mu_0_y = 1.0 - (double)bigtime / numPoints;
             }

             mu_2_x = (mu_0_x + mu_1_x) * cmd->quantum[5][0]; /**this parameter must be 0.5**/
             mu_2_y = (mu_0_y + mu_1_y) * cmd->quantum[5][0];

             time += 1;

             if  (time % numPoints == 0)
                 smallcounter += 1;
             if  (smallcounter == 3)
                 smallcounter = 0;

             if  (time % (4 * numPoints) == 0)
                 bigcounter += 1;
             if  (bigcounter == 3) {
                 bigcounter = 0;
                 fprintf (stderr, "\n\nrepeating ... \n\n");
             }
         }
     }

     if  (cmd->quantum[1][0] == 0) {
         cm_x = away_from_boundary + mu_0_x * (z->d_b - 1 - 2*away_from_boundary);
         cm_y = away_from_boundary + mu_0_y * (z->d_a - 1 - 2*away_from_boundary);
     }
     if  (cmd->quantum[1][0] == 1) {
         cm_x = away_from_boundary + mu_1_x * (z->d_b - 1 - 2*away_from_boundary); 
         cm_y = away_from_boundary + mu_1_y * (z->d_a - 1 - 2*away_from_boundary);
     }
     if  (cmd->quantum[1][0] == 2) {
         cm_x = away_from_boundary + mu_2_x * (z->d_b - 1 - 2*away_from_boundary); 
         cm_y = away_from_boundary + mu_2_y * (z->d_a - 1 - 2*away_from_boundary);
     }

     for (ct_t = begin; ct_t < end; ++ct_t)
         for (X = 0; X < z->d_a; ++X)
             for (B = 0; B < z->d_b; ++B) {

                 double diffx = B - cm_x;
                 double diffy = X - cm_y;

                 /**for periodic boundary**/
                 if  (periodic_boundary) {
                     diffx = (fabs(diffx) <= fabs(z->d_b / 2.0)) ? diffx : - z->d_b + fabs(diffx);
                     diffy = (fabs(diffy) <= fabs(z->d_b / 2.0)) ? diffx : - z->d_b + fabs(diffx);
                 }

                 cmd->S_target[ct_t][X * z->d_b + B] = height * exp (-0.5 * (diffx*diffx+diffy*diffy)/(sigma*sigma));
             }

     return (DOUBLE)(0);
}

/******************************** data_gauss_3areas_2D_anim2 *****************/
/* Places 3 Gaussians. 0,1 randomly and for no.2: mu_2 = (mu_0+mu_1)*q[4][0] */
/* For 3 different areas.                                                    */
/* Gaussians are Gaussians only along d_b, but horizontal along d_a.         */
/* NEW: avoids outer rim of units for placement of Gaussian centres !!!      */
/* q[0][0] = initialise new if 1                                             */
/* q[1][0] = area 0, 1 or the "sum"-area 2                                   */
/* q[...] others, possibly from data_gauss_3areas_2D_anim are ignored here   */
/* Use in a skript (export with rlen=1):                                     */
/* convert +append file_0.pnm file_1.pnm file_2.pnm file_app.pnm             */
/* convert -shave 1x1 obs_app.pnm obs_app_shave.pnm                          */
/* if also  convert -border 10x0 -bordercolor #858585 in.pnm out.pnm  then   */
/*    also -crop would have to be used, so leave this.                       */

DOUBLE data_gauss_3areas_2D_anim2 (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {
     int X, B;
     static double mu_0_x = 0.0; /**first Gaussian**/
     static double mu_0_y = 0.0;
     static double mu_1_x = 0.0; /**second Gaussian**/
     static double mu_1_y = 0.0;
     static double mu_2_x = 0.0; /**third Gaussian**/
     static double mu_2_y = 0.0;
     double cm_x = 0.0, cm_y = 0.0;

     int periodic_boundary = 0;
     int away_from_boundary = 1;   /**how far to place Gaussian centres away from outermost units**/

     AREA *z = A + cmd->area;

     /**warnings of wrong usage**/
     if  (cmd->anz_quantums < 2)
         fprintf (stderr, "\nwrong number of quantums in data_gauss_3areas_2D_anim2\n");

     double sigma     = 1.0;
     double height    = 1.0;

     int numPoints = 14;
     static int time = 0;
     int smalltime = time % numPoints;
     static int smallcounter = 0;
     int bigtime = (time / 4) % numPoints;
     static int bigcounter = 0;

     if  (cmd->quantum[0][0] == 1) {

             /**draw diagonal 1**/
             if  (smallcounter == 0) {
                 mu_1_x = (double)smalltime / numPoints;
                 mu_1_y = (double)smalltime / numPoints;
             }

             /**draw opposite 1**/
             if  (smallcounter == 1) {
                 mu_1_x = 1.0 - (double)smalltime / numPoints;
                 mu_1_y = 1.0;
             }

             /**draw bottom 1**/
             if  (smallcounter == 2) {
                 mu_1_x = 0.0;
                 mu_1_y = 1.0 - (double)smalltime / numPoints;
             }

             /**draw diagonal 0**/
             if  (bigcounter == 0) {
                 mu_0_x = (double)bigtime / numPoints;
                 mu_0_y = (double)bigtime / numPoints;
             }

             /**draw opposite 0**/
             if  (bigcounter == 1) {
                 mu_0_x = 1.0 - (double)bigtime / numPoints;
                 mu_0_y = 1.0;
             }

             /**draw bottom 0**/
             if  (bigcounter == 2) {
                 mu_0_x = 0.0;
                 mu_0_y = 1.0 - (double)bigtime / numPoints;
             }

             mu_2_x = (mu_0_x + mu_1_x) * 0.5;
             mu_2_y = (mu_0_y + mu_1_y) * 0.5;

             time += 1;

             if  (time % numPoints == 0)
                 smallcounter += 1;
             if  (smallcounter == 3)
                 smallcounter = 0;

             if  (time % (4 * numPoints) == 0)
                 bigcounter += 1;
             if  (bigcounter == 3) {
                 bigcounter = 0;
                 fprintf (stderr, "\n\nrepeating ... \n\n");
             }
     }

     if  (cmd->quantum[1][0] == 0) {
         cm_x = away_from_boundary + mu_0_x * (z->d_b - 1 - 2*away_from_boundary);
         cm_y = away_from_boundary + mu_0_y * (z->d_a - 1 - 2*away_from_boundary);
     }
     if  (cmd->quantum[1][0] == 1) {
         cm_x = away_from_boundary + mu_1_x * (z->d_b - 1 - 2*away_from_boundary); 
         cm_y = away_from_boundary + mu_1_y * (z->d_a - 1 - 2*away_from_boundary);
     }
     if  (cmd->quantum[1][0] == 2) {
         cm_x = away_from_boundary + mu_2_x * (z->d_b - 1 - 2*away_from_boundary); 
         cm_y = away_from_boundary + mu_2_y * (z->d_a - 1 - 2*away_from_boundary);
     }

         for (X = 0; X < z->d_a; ++X)
             for (B = 0; B < z->d_b; ++B) {

                 double diffx = B - cm_x;
                 double diffy = X - cm_y;

                 /**for periodic boundary**/
                 if  (periodic_boundary) {
                     diffx = (fabs(diffx) <= fabs(z->d_b / 2.0)) ? diffx : - z->d_b + fabs(diffx);
                     diffy = (fabs(diffy) <= fabs(z->d_b / 2.0)) ? diffx : - z->d_b + fabs(diffx);
                 }

                 cmd->S_target[0][X * z->d_b + B] = height * exp (-0.5 * (diffx*diffx+diffy*diffy)/(sigma*sigma));
             }

     if  (cmd->quantum[1][0] == 0) {

        for (X = 0; X < z->d_a; ++X)
             for (B = 0; B < z->d_b; ++B)
                 cmd->S_target[0][X * z->d_b + B] = 0.0;

         double radius = z->d_b * 0.5 * (-1.0 + time / 200.0);
         /* if  (fabs(radius) > z->d_b * 0.5 - 1) */
             radius = z->d_b * 0.5 - 1;
	 double inc = (double)(time / 100.0 * 6.283);
         B = (int)(z->d_b / 2.0 - 0.0 + radius * cos(inc));
         X = (int)(z->d_a / 2.0 - 0.0 + radius * sin(inc));
         cmd->S_target[0][X * z->d_b + B] = 1.0;
     }

     if  (cmd->quantum[1][0] == 1) {

         for (X = 0; X < z->d_a; ++X)
             for (B = 0; B < z->d_b; ++B)
                 cmd->S_target[0][X * z->d_b + B] = 0.0;

         /**draw an "F"**
         for (X = 1; X < z->d_a - 1; ++X) {
             B = z->d_b / 4;
             cmd->S_target[0][X * z->d_b + B] = 1.0;
         }
         X = 1;
         for (B = z->d_b / 4; B < 3*z->d_b/4+1; ++B)
             cmd->S_target[0][X * z->d_b + B] = 1.0;
         X = z->d_a / 2 - 1;
         for (B = z->d_b / 4; B < 3*z->d_b/4-1; ++B)
             cmd->S_target[0][X * z->d_b + B] = 1.0;
         **/

         /**draw a "C"**
         double half_axis_hor = z->d_b * 0.3;
         double half_axis_ver = z->d_a * 0.4;
	 for (double inc = 0.7; inc < 5.5; inc += 0.01) {
             B = (int)(z->d_b / 2.0 - 0.0 + half_axis_hor * cos(inc));
             X = (int)(z->d_a / 2.0 - 0.0 + half_axis_ver * sin(inc));
             cmd->S_target[0][X * z->d_b + B] = 1.0;
         }
         **/

         /**draw an "L"**
         for (X = 1; X < z->d_a - 1; ++X) {
             B = z->d_b / 3;
             cmd->S_target[0][X * z->d_b + B] = 1.0;
         }
         X = z->d_a - 2;
         for (B = z->d_b / 3; B < 3*z->d_b/4+1; ++B)
             cmd->S_target[0][X * z->d_b + B] = 1.0;
         **/

         /**draw an "X" (assuming square area)**
         for (X = 1; X < z->d_a - 1; ++X) {
             B = X;
             cmd->S_target[0][X * z->d_b + B] = 1.0;
         }
         for (X = 1; X < z->d_a - 1; ++X) {
             B = z->d_b - 1 - X;
             cmd->S_target[0][X * z->d_b + B] = 1.0;
         }
         **/

         /**draw a "+" (assuming square area)**/
         for (X = 1; X < z->d_a - 1; ++X) {
             B = z->d_b / 2;
             cmd->S_target[0][X * z->d_b + B] = 1.0;
         }
         for (B = 1; B < z->d_b - 1; ++B) {
             X = z->d_a / 2;
             cmd->S_target[0][X * z->d_b + B] = 1.0;
         }
         /**/
     }

     return (DOUBLE)(0);
}

/******************************** data_gauss_2D_anim3 ************************/
/* Anim leads a Gaussian around a circle and along the mid vert&horiz.       */
/* For one 2D area -- used for visualization of log-polar retina-SC mapping. */
/* From data_gauss_3areas_2D_anim2. Old comments:                            */
/* NEW: avoids outer rim of units for placement of Gaussian centres !!!      */
/* q[0][0] = initialise new if 1                                             */
/* q[1][0] = area 0, 1 or the "sum"-area 2                                   */
/* q[...] others, possibly from data_gauss_3areas_2D_anim are ignored here   */
/* Use in a skript (export with rlen=1):                                     */
/* convert +append file_0.pnm file_1.pnm file_2.pnm file_app.pnm             */
/* convert -shave 1x1 obs_app.pnm obs_app_shave.pnm                          */
/* if also  convert -border 10x0 -bordercolor #858585 in.pnm out.pnm  then   */
/*    also -crop would have to be used, so leave this.                       */

DOUBLE data_gauss_2D_anim3 (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {
     int X, B;
     static double mu_0_x = 0.0; /**first Gaussian**/
     static double mu_0_y = 0.0;
     double cm_x = 0.0, cm_y = 0.0;

     int periodic_boundary = 0;
     int away_from_boundary = 1;   /**how far to place Gaussian centres away from outermost units**/

     AREA *z = A + cmd->area;

     /**warnings of wrong usage**/
     if  (cmd->anz_quantums < 2)
         fprintf (stderr, "\nwrong number of quantums in data_gauss_2D_anim3\n");

     double sigma     = 1.0;
     double height    = 1.0;

     int numPoints = 50;
     static int time = 0;
     int bigtime = (time) % numPoints;
     static int bigcounter = 0;

     if  (cmd->quantum[0][0] == 1) {

             /**draw horizontal**/
             if  (bigcounter == 0) {
                 mu_0_x = (double)bigtime / numPoints;
                 mu_0_y = 0.5;
             }
             /**draw circle part**/
             if  (bigcounter == 1) {
                 double radius = 0.5;
                 double inc = (double)bigtime / numPoints * 3.1415927 / 2.0 + 3.1415927 / 2.0;
                 mu_0_y = 0.5 + radius * cos(inc);
                 mu_0_x = 0.5 + radius * sin(inc);
	     }
             /**draw vertical**/
             if  (bigcounter == 2) {
                 mu_0_x = 0.5;
                 mu_0_y = (double)bigtime / numPoints;
             }
             /**draw circle part**/
             if  (bigcounter == 3) {
                 double radius = 0.5;
                 double inc = -(double)bigtime / numPoints * 3.1415927 / 2.0;
                 mu_0_y = 0.5 + radius * cos(inc);
                 mu_0_x = 0.5 + radius * sin(inc);
	     }

             time += 1;

             if  (time % (numPoints) == 0)
                 bigcounter += 1;
             if  (bigcounter == 4) {
                 bigcounter = 0;
                 fprintf (stderr, "\n\nrepeating ... \n\n");
             }
     }

     if  (cmd->quantum[1][0] == 0) {
         cm_x = away_from_boundary + mu_0_x * (z->d_b - 1 - 2*away_from_boundary);
         cm_y = away_from_boundary + mu_0_y * (z->d_a - 1 - 2*away_from_boundary);
     }

         for (X = 0; X < z->d_a; ++X)
             for (B = 0; B < z->d_b; ++B) {

                 double diffx = B - cm_x;
                 double diffy = X - cm_y;

                 /**for periodic boundary**/
                 if  (periodic_boundary) {
                     diffx = (fabs(diffx) <= fabs(z->d_b / 2.0)) ? diffx : - z->d_b + fabs(diffx);
                     diffy = (fabs(diffy) <= fabs(z->d_b / 2.0)) ? diffx : - z->d_b + fabs(diffx);
                 }

                 cmd->S_target[0][X * z->d_b + B] = height * exp (-0.5 * (diffx*diffx+diffy*diffy)/(sigma*sigma));
             }

     return (DOUBLE)(0);
}



/******************************** data_half_circle ***************************/
/* Creates 3D data points: 2D for half-circle and 1D for rotation angle.     */
/* Used for self-organisation of long-range docking state space on a SOM.    */
/* q[0][0] = max distance from center (fruit location).                      */
/* Note: important for SOM are relative scales of the 3 dimensions.          */

DOUBLE data_half_circle (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {

     double max_dist = cmd->quantum[0][0];  /**circle radius**/
     double max_ang  = M_PI / 2.0;          /**cake-slice width**/
     double max_rot  = M_PI / 2.0;          /**robot facing away from center -- this is not the real rotation angle of the robot**/

     double dist = drand48() * max_dist;
     double ang  = -max_ang + drand48() * max_ang * 2.0;
     double rot  = -max_rot + drand48() * max_rot * 2.0;

     for (int ct_t = begin; ct_t < end; ++ct_t) {
         cmd->S_target[ct_t][0] = dist;
         cmd->S_target[ct_t][1] = ang;
         cmd->S_target[ct_t][2] = rot;
     }

     return (DOUBLE)(0);
}


/****************************** write_act_file *******************************/
/* Writes act at time step ct_t. Versatile function.                         */
/* q[0][0] = 1: write newline AFTER writing the data point.                  */
/* Previosly, data were written as INT;  Modify function as needed !!        */

DOUBLE write_act_file (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

  char fullname[512];
  FILE *fp;

  sprintf (fullname, "%s", cmd->pointers[0]->words[0]);  /**first pointer gives path/file-name**/
if  (cmd->anz_from1 > 1) {

  if  ((fp = fopen (fullname, "a")) == 0)
      fprintf (stderr, "\n\n\nError opening file %s for append\n", cmd->pointers[0]->words[0]);

  for (int ct_l = 0; ct_l < cmd->anz_from1; ++ct_l) {

      if  (ct_l == 0)
          fprintf (fp, "%d ", (int)(cmd->S_from1[ct_l][ct_t][0]));
      if  ((ct_l > 0) && (ct_l < 3))
          for (int i = 0; i < 2 /*A[inarea].d_n*/ ; ++i)
              fprintf (fp, "%f ", cmd->S_from1[ct_l][ct_t][i]);
      if  (ct_l == 3)
          for (int i = 0; i < 3 /*A[inarea].d_n*/ ; ++i)
              fprintf (fp, "%f ", cmd->S_from1[ct_l][ct_t][i]);
  }

} else {

  if  ((fp = fopen (fullname, "w")) == 0)
      fprintf (stderr, "\n\n\nError opening file %s for append\n", cmd->pointers[0]->words[0]);

      int inarea = cmd->n_from1[0];
      for (int i = 0; i < A[inarea].d_n; ++i)
          fprintf (fp, "%f ", cmd->S_from1[0][ct_t][i]);
}


  if  (cmd->quantum[0][0] == 1.0)
      fprintf (fp, "\n");
  else
      fprintf (fp, " ");

  fclose (fp);

  return (DOUBLE)0;
}


/****************************** read_act_file ********************************/
/* Reads act written by write_act_file. Versatile function.                  */
/* q[0][0] = data dim (length of a line in file). NEW!                       */
/* Previously, data were int;  Modify function as needed !!                  */

DOUBLE read_act_file (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

  char fullname[512];
  char kommentarzeile[512];
  FILE *fp;
  int anz;
  int data_dim = (int)(cmd->quantum[0][0]);

  if  (data_dim == 0) {
      fprintf (stderr, "\n\nread_act_file: data_dim is zero !!!\n");
      exit (1);
  }

  sprintf (fullname, "%s", cmd->pointers[0]->words[0]);  /**first pointer gives path/file-name**/
  if  ((fp = fopen (fullname, "r")) == 0)
      fprintf (stderr, "\n\n\nError opening file %s for read\n", cmd->pointers[0]->words[0]);

  for (anz = 0; ! feof (fp); ++anz)
      fgets (kommentarzeile, 512, fp);
  fprintf (stderr, "\nLine length of file (=num of data): %d\n", anz);

  anz--;  /**because the loop above has shot one over the goal**/

  fseek (fp, 0, SEEK_SET);

  float ** data = f_matrix (anz + 1, data_dim);

  data[0][0] = (float)anz;

  for (int i = 1; i <= anz; ++i)
      for (int j = 0; j < data_dim; ++j) {
          float val;
          fscanf (fp, "%f", &val);
          data[i][j] = val;
      }

  fclose (fp);

  cmd->pointers[0]->data = (void *)data;

  return (DOUBLE)0;
}



/****************************** dense_act_file *******************************/
/* Writes density of data onto each area. Versatile function.                */
/* q[0][0]=0/1/2 for 1st, 2nd and 3rd area (each 2dim; 2 int's denote loc).  */

DOUBLE dense_act_file (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

  int **data = (int **)cmd->pointers[0]->data;
  int anz = data[0][0];
  int thisarea = (int)cmd->quantum[0][0];    /*must be 0, 1 or 2*/

  for (int n = 0; n < A[cmd->area].d_n; ++n)
      cmd->S_target[ct_t][n] = 0.0;

  for (int i = 1; i <= anz; ++i) {
      int neuron_x = data[i][thisarea*2];
      int neuron_y = data[i][thisarea*2 + 1];

      if  (neuron_x * A[cmd->area].d_a + neuron_y >= A[cmd->area].d_n)
          fprintf (stderr, "\n\nWARNING: area burst in dense_act_file!\n\n");

      cmd->S_target[ct_t][neuron_x * A[cmd->area].d_b + neuron_y] += 1.0;
  }

  return (DOUBLE)0;
}


/****************************** data_act_file ********************************/
/* Writes x,y-component of data onto each area (2 units). Versatile function.*/
/* q[0][0] triggers new data point if != 0. =1: random data; =2: even prob.  */
/* q[1][0]=0/1/2 for 1st, 2nd and 3rd area (each 2dim; 2 int's denote loc).  */
/* Currently not yet consider data density ...                               */

DOUBLE data_act_file (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

  int **data = (int **)cmd->pointers[0]->data;
  int anz = data[0][0];
  int thisarea = (int)cmd->quantum[1][0];    /*must be 0, 1 or 2*/

  static int sel = 1;
  static DOUBLE * data_invfreq;
  static DOUBLE * data_upperinterval;
  static DOUBLE scale = 0.0;
  static int firsttime = 1;

  if  (firsttime) {

      data_invfreq = d_vector (anz + 1);
      data_upperinterval = d_vector (anz + 1);

      int inarea = cmd->n_from1[0];

      /**assign each data point its frequency of occurrence according to S_from**/
      for (int i = 1; i <= anz; ++i) {
          int neuron_x = data[i][4];
          int neuron_y = data[i][5];
          data_invfreq[i] = 1.0 / cmd->S_from1[0][ct_t][neuron_x * A[inarea].d_b + neuron_y];
      }

      /**add up inverse frequencies on a continuous scale and assign each data point an upper interval on this scale**/
      data_upperinterval[0] = 0.0;
      for (int i = 1; i <= anz; ++i) {
          scale += data_invfreq[i];
          data_upperinterval[i] = scale;
      }

      firsttime = 0;
  }

  /**even-probability mode**/
  if  (cmd->quantum[0][0] == 2.0) {
      DOUBLE scale_sel = scale * drand48();
      for (int i = 1; i <= anz; ++i)
          if  ((data_upperinterval[i-1] < scale_sel) && (data_upperinterval[i] >= scale_sel))
              sel = i;
  }


  /**normal random data mode**/
  if  (cmd->quantum[0][0] == 1.0)
      sel = (int)(drand48() * (anz-1.0) + 1.0);

  int neuron_x = data[sel][thisarea*2];
  int neuron_y = data[sel][thisarea*2 + 1];

  cmd->S_target[ct_t][0] = (DOUBLE)neuron_x;
  cmd->S_target[ct_t][1] = (DOUBLE)neuron_y;

  return (DOUBLE)0;
}



/****************************** file_flag ************************************/
/* Write/read an int to/from a file for communication between programs.      */
/* Pointer has double use for path/filename and int_val which to be returned.*/
/* q[0][0] = 0: write 0.                                                     */
/* q[0][0] = 1: write 1.                                                     */
/* q[0][0] = 2: return 0/1 if read 0/1.                                      */
/* q[0][0] = 3: return 1/0 if read 0/1.                                      */

DOUBLE file_flag (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

  char fullname[512];
  FILE *fp;

  sprintf (fullname, "%s", cmd->pointers[0]->words[0]);  /**first pointer gives path/file-name**/

  if  ((cmd->quantum[0][0] == 0) || (cmd->quantum[0][0] == 1)) {

      int val = (int)(cmd->quantum[0][0]);

      if  ((fp = fopen (fullname, "w")) == 0)
          fprintf (stderr, "\n\n\nError opening file %s for write\n", cmd->pointers[0]->words[0]);

      fprintf (fp, "%d", val);

      fclose (fp);

      fprintf (stderr, "\nwrote %d to file   ", val);
  }


  if  ((cmd->quantum[0][0] == 2) || (cmd->quantum[0][0] == 3)) {
      int val;

      if  ((fp = fopen (fullname, "r")) == 0)
          fprintf (stderr, "\n\n\nError opening file %s for read\n", cmd->pointers[0]->words[0]);

      fscanf (fp, "%d", &val);

      fclose (fp);

      if  (cmd->quantum[0][0] == 2)
          cmd->pointers[0]->int_val = val;

      if  (cmd->quantum[0][0] == 3)
          cmd->pointers[0]->int_val = !val;

      fprintf (stderr, "\nread %d from file ", val);
      if  (val != cmd->pointers[0]->int_val)
          fprintf (stderr, "but returning %d  !!  ", cmd->pointers[0]->int_val);
  }

  return (DOUBLE)0;
}



/****************************** int_val_change *******************************/
/* Returns 1 to pointers[0]->int_val                                         */
/* if pointers[1]->int_val changes from q[0][0] to q[1][0],                  */
/* else returns 0 to pointers[0]->int_val.                                   */
/* q[2][0] = initialisation value.                                           */

DOUBLE int_val_change (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy) {

  static int old_val = (int)(cmd->quantum[2][0]);  /**init**/

  if  (  (old_val                   == (int)(cmd->quantum[0][0]))
      && (cmd->pointers[1]->int_val == (int)(cmd->quantum[1][0]))) {

          cmd->pointers[0]->int_val = 1;
          fprintf (stderr, "\n\n int_val_change from %d to %d !  \n\n   ", old_val, cmd->pointers[1]->int_val);
  } else {
          cmd->pointers[0]->int_val = 0;
          fprintf (stderr, "\n int_val_change: old=%d new=%d     ", old_val, cmd->pointers[1]->int_val);
  }

  old_val = cmd->pointers[1]->int_val;

  return (DOUBLE)0;
}


/**************************** data_gnu_connect_max ***************************/
/* On a 2-Dim layer, connects the maxima of consecutive act patterns.        */
/* Used the variations in transformations with different inputs.             */
/* Works together with inset in  data_gauss_3areas_2D.                       */

DOUBLE data_gnu_connect_max (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {

 float max_act = 0.0;
 float X_winner = -1;
 float Y_winner = -1;
 static float X_old = -1;
 static float Y_old = -1;

        /**find winner: highest activation**
        for (int X = 0; X < A[cmd->area].d_a; X++)
            for (int Y = 0; Y < A[cmd->area].d_b; Y++)
                if  (cmd->S_from1[0][begin][Y + A[cmd->area].d_b * X] > max_act) {
                    max_act = cmd->S_from1[0][begin][Y + A[cmd->area].d_b * X];
                    X_winner = (float)X;
                    Y_winner = (float)Y;
                }
        **/

        /**find winner on "continuous" plane**/
        float inc = 0.5;
        for (float X = 0.0; X < A[cmd->area].d_a; X += inc)
        for (float Y = 0.0; Y < A[cmd->area].d_b; Y += inc) {

            float curr_act = 0.0;
            for (int XX = 0; XX < A[cmd->area].d_a; XX++)
            for (int YY = 0; YY < A[cmd->area].d_b; YY++) {
                float diffsq = (X-(float)XX)*(X-(float)XX) + (Y-(float)YY)*(Y-(float)YY);
                curr_act += cmd->S_from1[0][begin][YY + A[cmd->area].d_b * XX] * exp (-0.5 * diffsq);
            }

            if  (curr_act > max_act) {
                max_act = curr_act;
                X_winner = X;
                Y_winner = Y;
            }
        }
        /**/


 static int count_for_data_gnu_connect_max = -1;
 count_for_data_gnu_connect_max += 1;
 char fullname[512];
 sprintf (fullname, "%s/connect_max.gnu", cmd->pointers[0]->words[0]);    /**pointers[1] gives the directory!**/

 static int gridcounter = -1;
 static float CMs_a[15*15]; /*even though only 8*8 used in display*/
 static float CMs_b[15*15];

 if  (count_for_data_gnu_connect_max == 0)
     for (int i = 0; i < 15*15; i++) {
         CMs_a[i] = 0.0;
         CMs_b[i] = 0.0;
     }

 static int linestyle = 0;
 if  (count_for_data_gnu_connect_max % DATA_GNU_CONNECT_MAX_POINTS_PER_NODE == 0) {

         gridcounter += 1;
         linestyle = (int)(100.0 * drand48()) + 1;

         FILE *fp = fopen (fullname, "a");
         fprintf (fp, "set linestyle %d\n", linestyle);
         fclose (fp);

 } else {

         FILE *fp = fopen (fullname, "a");
         /* fprintf (fp, "set label \"x\" at %f,%f center\n", X_winner, Y_winner); */
         fprintf (fp, "set arrow from %f,%f to %f,%f nohead ls %d\n", X_old, Y_old, X_winner, Y_winner, linestyle);
         fclose (fp);
 }

 X_old = X_winner;
 Y_old = Y_winner;

 CMs_a[gridcounter] += X_winner / (float)(DATA_GNU_CONNECT_MAX_POINTS_PER_NODE);
 CMs_b[gridcounter] += Y_winner / (float)(DATA_GNU_CONNECT_MAX_POINTS_PER_NODE);

 fprintf (stderr, " %d ", gridcounter);
 if  (  (gridcounter == DATA_GNU_CONNECT_MAX_COLS_AND_LINES * DATA_GNU_CONNECT_MAX_COLS_AND_LINES - 1)
     && (count_for_data_gnu_connect_max % DATA_GNU_CONNECT_MAX_POINTS_PER_NODE == DATA_GNU_CONNECT_MAX_POINTS_PER_NODE - 1)) {

         fprintf (stderr, "\nplotting final grid\n");

         FILE *fp = fopen (fullname, "a");

         fprintf (fp, "set linestyle 1\n");

         for (int ct_a = 0; ct_a < DATA_GNU_CONNECT_MAX_COLS_AND_LINES; ++ct_a)
             for (int ct_b = 1; ct_b < DATA_GNU_CONNECT_MAX_COLS_AND_LINES; ++ct_b)
                 fprintf (fp, "set arrow from %f,%f to %f,%f nohead ls 1\n", CMs_a[ct_a * DATA_GNU_CONNECT_MAX_COLS_AND_LINES + ct_b - 1], CMs_b[ct_a * DATA_GNU_CONNECT_MAX_COLS_AND_LINES + ct_b - 1], CMs_a[ct_a * DATA_GNU_CONNECT_MAX_COLS_AND_LINES + ct_b], CMs_b[ct_a * DATA_GNU_CONNECT_MAX_COLS_AND_LINES + ct_b]);

         for (int ct_a = 1; ct_a < DATA_GNU_CONNECT_MAX_COLS_AND_LINES; ++ct_a)
             for (int ct_b = 0; ct_b < DATA_GNU_CONNECT_MAX_COLS_AND_LINES; ++ct_b)
                 fprintf (fp, "set arrow from %f,%f to %f,%f nohead ls 1\n", CMs_a[(ct_a - 1) * DATA_GNU_CONNECT_MAX_COLS_AND_LINES + ct_b], CMs_b[(ct_a - 1) * DATA_GNU_CONNECT_MAX_COLS_AND_LINES + ct_b], CMs_a[ct_a * DATA_GNU_CONNECT_MAX_COLS_AND_LINES + ct_b], CMs_b[ct_a * DATA_GNU_CONNECT_MAX_COLS_AND_LINES + ct_b]);

         fclose (fp);
 }

 return (DOUBLE)0;
}
