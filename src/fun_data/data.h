#ifndef _data_H
#define _data_H

/**
  was originally in PARAMS:
  int (*data_import)(struct PARAMS *x);
  char *data_import_name;
  char **data_files;
  int  anz_data_files;
**/

typedef struct DATA {
  int    anzahl;
  int    *Ho_a;                      /**langsam zaehlender index**/
  int    *Br_b;                      /**schnell zaehlender index**/
  /**for image files**/
  DOUBLE **Bilder;
  int    *Bild_grau;
  double *Bild_max;
  double *Bild_min;
  double *Bild_mindiff;
} DATA;


/**from data.images.c**/
DOUBLE import_images        (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
int import_points           (PARAMS *g, DATA *d);
DOUBLE init_image           (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE init_whole_image     (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE init_onoff_image     (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE cut_image            (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE cut_image_pantilt    (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE init_orange          (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE image_color_blob     (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE init_image_cosinus   (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);

/**from data.digits.c**/
DOUBLE init_digit           (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE init_mnist           (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);


/**from data.made.c**/
DOUBLE init_lines_hier      (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE init_lines_mixed     (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE init_lines_sparse    (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE init_points          (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE init_burglar         (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE init_cuecomb         (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE data_sub_mean        (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);

DOUBLE data_rand_gibbs_01   (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE data_hack_init_lines (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE data_zeppelin        (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE data_gauss           (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE data_gauss_nontorus  (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE data_gauss_move      (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE data_gauss_three     (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE data_gauss_3areas    (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE data_gauss_3areas_2D (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE data_gauss_fromfile  (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE data_gauss_3areas_2D_anim (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE data_gauss_3areas_2D_anim2 (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE data_gauss_2D_anim3  (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE data_half_circle     (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE write_act_file       (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy);
DOUBLE read_act_file        (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy);
DOUBLE dense_act_file       (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy);
DOUBLE data_act_file        (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy);
DOUBLE file_flag            (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy);
DOUBLE int_val_change       (PARAMS *g, AREA *A, COMMAND *cmd, int ct_t, int dummy);
DOUBLE data_gnu_connect_max (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);

DOUBLE data_gauss_move3     (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE data_gauss_move4     (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE data_gauss_motor     (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE data_zhang           (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE data_block           (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE data_trichter        (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE data_gauss123        (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);


/**from data.webots.c**/
#if WEBOTS
DOUBLE webot_supervisor   (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE webot_supervisor_2 (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE webot_rotate       (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE webot_pioneer      (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
#endif

/**from data.bttv.c**/
DOUBLE bttv_image         (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);


/**from data.fixed.c**/
/** DATA *import_snns_data (PARAMS *g);   **later maybe**/


/**from data.motor.c**/
int import_motor          (PARAMS *g);
int import_lang_assoc     (PARAMS *g);
DOUBLE data_motor         (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE data_lang_assoc    (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
#if  USE_CORBA
DOUBLE data_read_understood(PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
DOUBLE data_write_string  (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end);
#endif

#endif
