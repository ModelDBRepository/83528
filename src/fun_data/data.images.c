/**taken from /home/ni/cweber/p/co/bm/src/images.c**/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>  /**for M_PI in ...gabor**/
#include "../kernel/coco.h"
#include "../kernel/series.h"
#include "../kernel/utils.h"
#include "data.h"
#include "../fun_basic/total.h"     /**for init_image_cosinus**/
#include "../fun_basic/local.h"

/************************** importP5_img *************************************/
/* Allocates double vector d->Bilder[.] and reads images from 3) "datei".    */
/* New: gets the format (P5 for raw bw or P6 for raw color).                 */

void importP5_img (DATA *d, int bild_nr, char *datei, int *format) {
    int i, num_colors = 1;
    FILE *fp;
    char kommentarzeile[256];

    fprintf(stderr, "\n%d<--%18s", bild_nr, datei);
    if  ((fp = fopen (datei, "r")) == NULL) {
        fprintf(stderr, "\nno file %s in\n", datei);  exit(0);
    }

    /**get the P5 or the P6**/
    fscanf (fp, "%s\n", kommentarzeile);
    fprintf (stderr, " %s ", kommentarzeile);

    /**get the format**/
    if  (kommentarzeile[1] == '5') {
        *format = 5;
        num_colors = 1;
        fprintf (stderr, " format=5 ");
    } else {
        if  (kommentarzeile[1] == '6') {
            *format = 6;
            num_colors = 3;
            fprintf (stderr, " format=6 ");
	} else
	  fprintf (stderr, "\n\n%s: image format not known!\n\n", datei);
    }


    /**forget comments; read the 2 lines with Br_b, Ho_a, Bild_grau**/
    for (i = 0; i < 2; ) {
        fgets (kommentarzeile, 256, fp);
        if  ((kommentarzeile[0] != '#') && (strlen(kommentarzeile) > 1)) {
            if  (i == 0)
                sscanf (kommentarzeile, "%d %d",
                        &d->Br_b[bild_nr], &d->Ho_a[bild_nr]);
            if  (i == 1)
                sscanf (kommentarzeile, "%d", &d->Bild_grau[bild_nr]); 
            ++i;
        }
    }

    fprintf (stderr, " br:%d ho:%d grau:%d",
             d->Br_b[bild_nr], d->Ho_a[bild_nr], d->Bild_grau[bild_nr]);

    if  (d->Bild_grau[bild_nr] != 255)
        fprintf (stderr, "\n\nWarning: not 255 grey levels!\n\n");

    d->Bilder[bild_nr] = d_vector (d->Br_b[bild_nr] * d->Ho_a[bild_nr] * num_colors);

    /**read the data**/
       for (i = 0; i < d->Br_b[bild_nr] * d->Ho_a[bild_nr] * num_colors; ++i)
            d->Bilder[bild_nr][i] = (double)(fgetc (fp));

    fclose (fp);

    fprintf (stderr, " ... read");
}



/**************************** import_images **********************************/
/* Allocates memory for struct d and points cmd->pointers[0]->data to it.    */
/* Calls importP5_img for images as contents.                                */
/* Use_d_ sub_mean_vector and spherize_vector for each image file.           */
/* Now: NO preprocessing. Only d->Bild_max,min,mindiff[bild_nr] retreived.   */
/* If color image, then d->Bilder[bild_nr] is 3=num_colors times as long.    */

DOUBLE import_images (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {

    int bild_nr, i;
    int format, num_colors = 1;

    DATA *d = (DATA *)malloc (sizeof (DATA));
    cmd->pointers[0]->data = d;

    d->anzahl = cmd->pointers[0]->num_words;

    /** get memory for image information **/
    if  ((d->Br_b = (int *)malloc (d->anzahl * sizeof (int))) == NULL)
        fprintf (stderr, "\nallocation failure for d->Br_b\n");
    if  ((d->Ho_a = (int *)malloc (d->anzahl * sizeof (int))) == NULL)
        fprintf (stderr, "\nallocation failure for d->Ho_a\n");
    if  ((d->Bild_grau = (int *)malloc (d->anzahl * sizeof (int))) == NULL)
        fprintf (stderr, "\nallocation failure for d->Bild_grau\n");
    if  ((d->Bild_max = (double *)malloc(d->anzahl * sizeof (double))) == NULL)
        fprintf (stderr, "\nallocation failure for d->Bild_\n");
    if  ((d->Bild_min = (double *)malloc(d->anzahl * sizeof (double))) == NULL)
        fprintf (stderr, "\nallocation failure for d->Bild_\n");
    if  ((d->Bild_mindiff = (double *)malloc (d->anzahl * sizeof (double)))
         == NULL)  fprintf (stderr, "\nallocation failure for d->Bild_\n");

    /** get memory for our pointers to pictures **/
    if  ((d->Bilder = (DOUBLE **)malloc(d->anzahl*sizeof(DOUBLE *))) == NULL)
        fprintf (stderr, "\nout of memory\n");

    for (bild_nr = 0; bild_nr < d->anzahl; ++bild_nr) {

        importP5_img (d, bild_nr, cmd->pointers[0]->words[bild_nr], &format);

        if  (format == 5)
            num_colors = 1;
        else {
            if  (format == 6)
                num_colors = 3;
            else
                fprintf (stderr, "\n\nimport_images: format not known!\n\n");
        }

        /** !!! 
        sub_mean_vector (d->Bilder[bild_nr], d->Br_b[bild_nr] * d->Ho_a[bild_nr] * num_colors);
        spherize_vector (d->Bilder[bild_nr], d->Br_b[bild_nr] * d->Ho_a[bild_nr] * num_colors);
        **/

        /**get maximum and minimum of each image ** from end of importP5_img**/
        d->Bild_max[bild_nr] = 0.0;
        d->Bild_min[bild_nr] = 255.0;
        for (i = 0; i < d->Br_b[bild_nr] * d->Ho_a[bild_nr] * num_colors; ++i) {
            d->Bild_max[bild_nr] = d->Bilder[bild_nr][i] > d->Bild_max[bild_nr]
                                 ? d->Bilder[bild_nr][i]
                                 : d->Bild_max[bild_nr];
            d->Bild_min[bild_nr] = d->Bilder[bild_nr][i] < d->Bild_min[bild_nr]
                                 ? d->Bilder[bild_nr][i]
                                 : d->Bild_min[bild_nr];
        }
       d->Bild_mindiff[bild_nr] = (d->Bild_max[bild_nr] - d->Bild_min[bild_nr])
                                 / 10.0;

        fprintf (stderr, " max:%.2f min:%.2f",
                 d->Bild_max[bild_nr], d->Bild_min[bild_nr]);

    }

    return (DOUBLE)(0);
}


/**************************** import_points **********************************/
/* From import_images. Writes just one image using init_points.              */

int import_points (PARAMS *g, DATA *d) {
    int bild_nr = 0, i;

    d->anzahl = 1;

    /** get memory for image information **/
    if  ((d->Br_b = (int *)malloc (d->anzahl * sizeof (int))) == NULL)
        fprintf (stderr, "\nallocation failure for d->Br_b\n");
    if  ((d->Ho_a = (int *)malloc (d->anzahl * sizeof (int))) == NULL)
        fprintf (stderr, "\nallocation failure for d->Ho_a\n");
    if  ((d->Bild_grau = (int *)malloc (d->anzahl * sizeof (int))) == NULL)
        fprintf (stderr, "\nallocation failure for d->Bild_grau\n");
    if  ((d->Bild_max = (double *)malloc(d->anzahl * sizeof (double))) == NULL)
        fprintf (stderr, "\nallocation failure for d->Bild_\n");
    if  ((d->Bild_min = (double *)malloc(d->anzahl * sizeof (double))) == NULL)
        fprintf (stderr, "\nallocation failure for d->Bild_\n");
    if  ((d->Bild_mindiff = (double *)malloc (d->anzahl * sizeof (double)))
         == NULL)  fprintf (stderr, "\nallocation failure for d->Bild_\n");

    /** get memory for our pointers to pictures **/
    if  ((d->Bilder = (DOUBLE **)malloc(d->anzahl*sizeof(DOUBLE *))) == NULL)
        fprintf (stderr, "\nout of memory\n");

    d->Ho_a[bild_nr] = 500;
    d->Br_b[bild_nr] = 500;

    d->Bilder[bild_nr] = d_vector (d->Br_b[bild_nr] * d->Ho_a[bild_nr]);

    {
        COMMAND aux_cmd;
        AREA aux_z;

        aux_z.d_a = d->Ho_a[bild_nr];
        aux_z.d_b = d->Br_b[bild_nr];
        aux_cmd.S_target      = (DOUBLE **)malloc (sizeof (DOUBLE *));
        aux_cmd.S_target[0]   = (DOUBLE *) malloc (d->Ho_a[bild_nr] * d->Br_b[bild_nr] * sizeof (DOUBLE));
        aux_cmd.quantum       = (DOUBLE **)malloc (2 * sizeof (DOUBLE *));
        aux_cmd.quantum[0]    = (DOUBLE *) malloc (sizeof (DOUBLE));
        aux_cmd.quantum[1]    = (DOUBLE *) malloc (sizeof (DOUBLE));
        aux_cmd.anz_quantums  = 2;
        aux_cmd.quantum[0][0] = 1;
        aux_cmd.quantum[1][0] = 1;

        fprintf (stderr, "\ncalling init_points once for the big %dx%d image", d->Ho_a[bild_nr], d->Br_b[bild_nr]);

// in data.made.c        init_points (g, &aux_z, &aux_cmd, d, 0, 1);

        fprintf (stderr, "...done ");

        d->Bilder[bild_nr] = aux_cmd.S_target[0];
    }

    /**get maximum and minimum of each image ** from end of importP5_img**/
    d->Bild_max[bild_nr] = 0.0;
    d->Bild_min[bild_nr] = 255.0;
    for (i = 0; i < d->Br_b[bild_nr] * d->Ho_a[bild_nr]; ++i) {
        d->Bild_max[bild_nr] = d->Bilder[bild_nr][i] > d->Bild_max[bild_nr]
                             ? d->Bilder[bild_nr][i]
                             : d->Bild_max[bild_nr];
        d->Bild_min[bild_nr] = d->Bilder[bild_nr][i] < d->Bild_min[bild_nr]
                             ? d->Bilder[bild_nr][i]
                             : d->Bild_min[bild_nr];
    }
    d->Bild_mindiff[bild_nr] = (d->Bild_max[bild_nr] - d->Bild_min[bild_nr])
                             / 10.0;

    fprintf (stderr, " max:%.2f min:%.2f",
             d->Bild_max[bild_nr], d->Bild_min[bild_nr]);

    return 0;
}


/******************************** get_image_color ****************************/
/* Cut out an image part at ho_pos, br_pos in size of d_a / 3, d_b.          */
/* Returns 1 if part is over the border of the whole image.                  */

int get_image_color (DATA *d, int bild_nr, DOUBLE *Image, int d_a, int d_b,
                                                    int ho_pos, int br_pos) {
    int k, i, j;
    int border = 0;

    /**test if input fits for color**/
    if  (d_a % 3 != 0)
        fprintf (stderr, "\n\nd_a must be a multiple of 3 for color images\n");

    /**move ectopic locations onto the edge of the image**/
    if  (ho_pos < 0) {
        ho_pos = 0;
        border = 1;
    }
    if  (br_pos < 0) {
        br_pos = 0;
        border = 1;
    }
    if  (ho_pos > d->Ho_a[bild_nr] - d_a / 3) {
        ho_pos = d->Ho_a[bild_nr] - d_a / 3;
        border = 2;
    }
    if  (br_pos > d->Br_b[bild_nr] - d_b) {
        br_pos = d->Br_b[bild_nr] - d_b;
        border = 2;
    }

    /**red**/
    for (i = 0; i < d_a / 3; ++i)
        for (j = 0; j < d_b; ++j)
            Image[j + i * d_b]
            = d->Bilder[bild_nr][(br_pos + ho_pos * d->Br_b[bild_nr]
                                  + j    +      i * d->Br_b[bild_nr]) * 3];
    /**green**/
    for (k = 0; k < d_a / 3; ++k) {
        i = k + d_a / 3;
        for (j = 0; j < d_b; ++j)
            Image[j + i * d_b]
            = d->Bilder[bild_nr][(br_pos + ho_pos * d->Br_b[bild_nr]
                                  + j    +      k * d->Br_b[bild_nr]) * 3 + 1];
    }
    /**blue**/
    for (k = 0; k < d_a / 3; ++k) {
        i = k + 2 * d_a / 3;
        for (j = 0; j < d_b; ++j)
            Image[j + i * d_b]
            = d->Bilder[bild_nr][(br_pos + ho_pos * d->Br_b[bild_nr]
                                  + j    +      k * d->Br_b[bild_nr]) * 3 + 2];
    }

    return (border);
}

/******************************** get_image **********************************/
/* Cut out an image part at ho_pos, br_pos in size of d_a, d_b.              */
/* Returns 1 if part is over the border of the whole image.                  */

int get_image (DATA *d, int bild_nr, DOUBLE *Image, int d_a, int d_b,
                                                    int ho_pos, int br_pos) {
    int i, j;
    int border = 0;

    /**move ectopic locations onto the edge of the image**/
    if  (ho_pos > d->Ho_a[bild_nr] - d_a) {
        ho_pos = d->Ho_a[bild_nr] - d_a;       /**<-- this could be too far, with the "-1", but (see below) won't harm**/
        border = 2;
    }
    if  (br_pos > d->Br_b[bild_nr] - d_b) {
        br_pos = d->Br_b[bild_nr] - d_b;       /**<-- this could be too far, with the "-1", but (see below) won't harm**/
        border = 2;
    }
    if  (ho_pos < 0) {
        ho_pos = 0;
        border = 1;
    }
    if  (br_pos < 0) {
        br_pos = 0;
        border = 1;
    }

    for (i = 0; i < d_a; ++i)
        for (j = 0; j < d_b; ++j)
            Image[j + i * d_b]
            = d->Bilder[bild_nr][br_pos + ho_pos * d->Br_b[bild_nr]
                                 + j    +      i * d->Br_b[bild_nr]];

    /**just a debugging rest ...**/
    for (j = 0; j < d_a * d_b; ++j)
        if  (Image[j] > 255) {
            fprintf (stderr, "\nget_image: bildnr=%d d_a=%d d_b=%d ho_pos=%d br_pos=%d", bild_nr, d_a, d_b, ho_pos, br_pos);
            fprintf (stderr, "\nd->Br_b[bild_nr]= %d  ", d->Br_b[bild_nr]);
            fprintf (stderr, "\ntoo large Image[%d]: %.2f ", j, Image[j]);
            fprintf (stderr, "\nwhile d->Bilder[%d][%d] = %f ", bild_nr, j, d->Bilder[bild_nr][j]);
            exit (1);
        }


/**  if  (drand48 () <= 0.5) { ** images not rotated anymore **
	... **upper part**
     } else { **rotated by 90'**
        **hier: d_b || hoehe**
        ho_pos = (int)(drand48 () * (double)(d->Ho_a[bild_nr] - d_b));
        **hier: d_a || breite**
        br_pos = (int)(drand48 () * (double)(d->Br_b[bild_nr] - d_a));

        for (i = 0; i < d_a; ++i)
            for (j = 0; j < d_b; ++j)
                Image[j + d_b * i]
                = d->Bilder[bild_nr][br_pos + ho_pos * d->Br_b[bild_nr]
                                     + i    +      j * d->Br_b[bild_nr]];
    }
**/

    return (border);
}

/******************************** saccade ************************************/
/* Go to a new image and a new starting point there.                         */

void saccade (AREA *z, DATA *d, int edge, int *bild_nr,
              int *ho_pos, int *br_pos, double *dist_ho, double *dist_br, int color) {

        int height = z->d_a;

        if  (color == 1)
            height = z->d_a / 3;

        *bild_nr = (int)(drand48 () * (double)(d->anzahl));
        *ho_pos = edge + (int)(drand48 () * (double)(d->Ho_a[*bild_nr] - height - edge - edge));         /**d_a || hoehe:  langsam zaehlender index**/
        *br_pos = edge + (int)(drand48 () * (double)(d->Br_b[*bild_nr] - z->d_b - edge - edge));         /**d_b || breite: schnell zaehlender index**/
        *dist_ho = 0.0;
        *dist_br = 0.0;
}


/******************************** init_image *********************************/
/* cmd->S_target is initialized with an image part of size z->d_a * z->d_b.  */
/* Calls get_image.                                                          */
/* quantum[0][0] each q[0][0]'th time get new image (also at border!).       */
/* quantum[1][0] scales the image pixels (for bounded transfkt of backprop). */
/* quantum[2][0] == 1: sub_mean_vector of each datapoint.                    */
/* quantum[3][0] == 1: spherize; == 2: normalize.                            */
/* quantum[4][0]/[4][1]: sigma of image shift vert/horiz at each time step.  */
/* quantum[5][0] == 1: color, else black&white.    new!                      */

DOUBLE init_image (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {

    int j, ct_t;
    static int firsttime = 1;
    static DOUBLE *image;
    int d_r = A[cmd->area].d_a * A[cmd->area].d_b;;
    // static long int idum[1] = {-1};
    int edge = (cmd->quantum[4][0] == 0 && cmd->quantum[4][1] == 0) ? 0 : 5;   /**changed edge from 2 to 0 !!!**/

    static int count_me = 0;
    static int border = 1;
    static int bild_nr;
    static int ho_pos;
    static int br_pos;
    static double dist_ho;
    static double dist_br;

    DATA *d = (DATA *)cmd->pointers[0]->data;

    if  (d->anzahl == 0)
        fprintf (stderr, "\nthere are no images in sight!\n");
    if  (cmd->anz_quantums != 6)
        fprintf (stderr, "\n6 parameters for init_image, please (color new)!\n");
    if  (cmd->quantum[1][0] == 0.0)
        fprintf (stderr, "\n\ninit_image wants quantum[1][0] set!\n\n");
    if  (cmd->anz_quant[4] != 2)
        fprintf (stderr, "\n\ninit_image wants 2 quant[4] (new!)!\n\n");

    /**allocate mem to store image**/
    if  (firsttime) {
        image = d_vector (d_r);
        firsttime = 0;
    }

        /****take only images with some intensity differences**
             double bild_min, bild_max; int badimage = 1; do  { ...
             bild_min = 255.0 ; bild_max = 0.0; for (j = 0; j < d_r; ++j) {
                 bild_min = image[j] < bild_min ? image[j] : bild_min;
                 bild_max = image[j] > bild_max ? image[j] : bild_max;
             } if  ((bild_max - bild_min) > d->Bild_mindiff[bild_nr])
             badimage = 0; } while (badimage);
        ****/

    if  (count_me % (int)(cmd->quantum[0][0]) == 0)
        saccade (A + cmd->area, d, edge, &bild_nr, &ho_pos, &br_pos, &dist_ho, &dist_br, (int)(cmd->quantum[5][0]));
                /*risk here!!*/                                                          /*color*/

    ct_t = begin;
    do {
        if  (cmd->quantum[5][0] == 0)
            border = get_image (d, bild_nr, image, A[cmd->area].d_a, A[cmd->area].d_b, ho_pos + (int)dist_ho, br_pos + (int)dist_br);
        else
            border = get_image_color (d, bild_nr, image, A[cmd->area].d_a, A[cmd->area].d_b, ho_pos + (int)dist_ho, br_pos + (int)dist_br);

        if  (cmd->quantum[2][0])                  /**must NOT be set for ICA**/
            sub_mean_vector (image, d_r);         /**but should for BM / old**/

        if  (cmd->quantum[3][0] == 1)
            spherize_vector (image, d_r);

        if  (cmd->quantum[3][0] == 2)
            normalize (image, d_r, 1.0);


        if  (cmd->quantum[4][0] > 0.0)    dist_ho += cmd->quantum[4][0] * (-0.5 + drand48()); /**not tested!!**/

        if  (cmd->quantum[4][0] < 0.0)    dist_ho -= cmd->quantum[4][0];

        if  (cmd->quantum[4][1] > 0.0)    dist_br += cmd->quantum[4][1] * (-0.5 + drand48()); /**not tested!!**/

        if  (cmd->quantum[4][1] < 0.0)    dist_br -= cmd->quantum[4][1];
 
        for (j = 0; j < d_r; ++j)
            cmd->S_target[ct_t][j] = image[j] * cmd->quantum[1][0]; /**scale**/

        ct_t ++;

        if  (border)
            saccade (A + cmd->area, d, edge, &bild_nr, &ho_pos, &br_pos, &dist_ho, &dist_br, (int)(cmd->quantum[5][0]));
                    /*risk here!!*/                                                          /*color*/

    } while (ct_t < end);

    count_me += 1;

    return (DOUBLE)(0);
}



/******************************** init_whole_image ***************************/
/* cmd->S_target is initialized with an image part of size z->d_a * z->d_b.  */
/* Calls get_image.                                                          */
/* quantum[0][0] each q[0][0]'th time get new image (also at border!).       */
/* quantum[1][0] scales the image pixels (for bounded transfkt of backprop). */
/* quantum[2][0] == 1: sub_mean_vector of each datapoint.                    */
/* quantum[3][0] == 1: spherize; == 2: normalize.                            */
/* quantum[4][0]/[4][1]: dummy for compatibility with init_image.         .  */

DOUBLE init_whole_image (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {

    int j, ct_t, bild_nr;
    static int firsttime = 1;
    static DOUBLE *image;
    int d_r = A[cmd->area].d_a * A[cmd->area].d_b;;

    DATA *d = (DATA *)cmd->pointers[0]->data;

    if  (d->anzahl == 0)
        fprintf (stderr, "\nthere are no images in sight!\n");
    if  (cmd->anz_quantums != 5)
        fprintf (stderr, "\n5 parameters for init_whole_image, please!\n");
    if  (cmd->quantum[1][0] == 0.0)
        fprintf (stderr, "\n\ninit_whole_image wants quantum[1][0] set!\n\n");
    if  (cmd->anz_quant[4] != 2)
        fprintf (stderr, "\n\ninit_whole_image wants 2 quant[4] (dummy!)!\n\n");


    /**allocate mem to store image**/
    if  (firsttime) {
        image = d_vector (d_r);
        firsttime = 0;
    }

    bild_nr = (int)(drand48 () * (double)(d->anzahl));

    ct_t = begin;
    do {
        get_image (d, bild_nr, image, A[cmd->area].d_a, A[cmd->area].d_b, 0, 0);

        if  (cmd->quantum[2][0])                  /**must NOT be set for ICA**/
            sub_mean_vector (image, d_r);         /**but should for BM / old**/

        if  (cmd->quantum[3][0] == 1)
            spherize_vector (image, d_r);

        if  (cmd->quantum[3][0] == 2)
            normalize (image, d_r, 1.0);

        for (j = 0; j < d_r; ++j)
            cmd->S_target[ct_t][j] = image[j] * cmd->quantum[1][0]; /**scale**/

        ct_t ++;

    } while (ct_t < end);

    return (DOUBLE)(0);
}





/******************************** init_onoff_image ***************************/
/* Calls init_image. Needs a double-zise area (d_a=2*d_b) for ON and OFF.    */

DOUBLE init_onoff_image (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {

    int d_r, d_r2, ct_t, j, max_shift = 0;

    if  (A[cmd->area].d_a != 2 * A[cmd->area].d_b)
        fprintf (stderr, "\n\ninit_onoff_image wants d_a = 2 * d_b!\n\n");

    /**directly from init_image**/
    if  (cmd->quantum[4][0] != 0.0)
        max_shift = 5;
    d_r = (A[cmd->area].d_a + (2*max_shift)) * (A[cmd->area].d_b + (2*max_shift));
    d_r2 = d_r / 2;



//    zz.d_a = A[cmd->area].d_a / 2;

// The area size passed to init_image must be halved !!!
// Do this by passing another "A"-argument!
    init_image (g, A, cmd, begin, end);

    for (ct_t = begin; ct_t < end; ++ct_t) {

        /**first(!) set 2nd-half, formerly negative, pixels**/
        for (j = 0; j < d_r2; ++j)
            if  (cmd->S_target[ct_t][j] > 0.0)
                cmd->S_target[ct_t][j + d_r2] = 0.0;
            else
                cmd->S_target[ct_t][j + d_r2] = -cmd->S_target[ct_t][j];

        /**then(!) cut first-half pixels to be positive**/
        for (j = 0; j < d_r2; ++j)
            if  (cmd->S_target[ct_t][j] < 0.0)
                cmd->S_target[ct_t][j] = 0.0;
    }

    return (DOUBLE)(0);
}


/******************************** cut_image **********************************/
/* Replaces init_image. No preprocessing. No edges. No shift. No memory.     */
/* Later: do filtering on a larger area and also use both for shift.         */
/* cmd->S_target is set to an image part of size A[area].d_a * A[area].d_b.  */
/* Calls get_image(_color).                  Value range is still 0 .. 255.  */
/* q[0][0] >=0 <1: minimal pixel range as fraction of range of whole image.  */
/* q[1][0] == 3: color, else black&white.  (NEW!!!)                          */
/* q[1][1] == 2: ON/OFF separate.          (NEW!!!)                          */
/* q[2][0], q[2][1] == br_pos, ho_pos; random if == -1; middle if -2.        */

/* Use  import_images  first!                                                */

DOUBLE cut_image (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {

    int j, ct_t;
    static int firsttime = 1;
    static DOUBLE *image;

    int bild_nr, ho_pos = 0, br_pos = 0;

    int badimage, border;

    DATA *d = (DATA *)cmd->pointers[0]->data;

    if  (d->anzahl == 0)
        fprintf (stderr, "\nthere are no images in sight!\n");

    if  (cmd->anz_quantums != 3)
        fprintf (stderr, "\n3 params for cut_image, please! New: pos_y+x\n");
    if  (cmd->anz_quant[1] != 2)
        fprintf (stderr, "\n2 parameters for cut_image, q[1], please!\n");
    if  (cmd->anz_quant[2] != 2)
        fprintf (stderr, "\n2 parameters for cut_image, q[2], please!\n");

    const int d_a = A[cmd->area].d_a;
    const int d_b = A[cmd->area].d_b;

    int height = d_a;
    int width  = d_b;
    if  (cmd->quantum[1][0] == 3)   /**color**/
        height = d_a / 3;
    if  (cmd->quantum[1][1] == 2)   /**ON/OFF**/
        width = d_b / 2;

    /**allocate mem to store image**/
    if  (firsttime) {
        image = d_vector (d_a * width);  /**d_a contains 3 colors**/
        firsttime = 0;
    }

    /**take only images with some intensity differences**/
    badimage = 1;
    do  {

        /**choose random pic and patch; should also work for patch_size == image_size**/
        bild_nr = (int)(drand48 () * (double)(d->anzahl));

        /**random position**/
        if  (cmd->quantum[2][0] == -1)
            br_pos  = (int)(floor)(drand48 () * (double)(d->Br_b[bild_nr] - width));     /**d_b || breite: schneller index**/
        if  (cmd->quantum[2][1] == -1)
            ho_pos  = (int)(floor)(drand48 () * (double)(d->Ho_a[bild_nr] - height));    /**d_a || hoehe:  langsamer index**/

        /**defined position**/
        if  (cmd->quantum[2][0] >= 0)
            br_pos  = (int)(cmd->quantum[2][0]);
        if  (cmd->quantum[2][1] >= 0)
            ho_pos  = (int)(cmd->quantum[2][1]);

        /**middle position**/
        if  (cmd->quantum[2][0] == -2)
            br_pos  = (d->Br_b[bild_nr] - width) / 2;
        if  (cmd->quantum[2][1] == -2)
            ho_pos  = (d->Ho_a[bild_nr] - height) / 2;


        /**cut chosen patch out**/
        if  (cmd->quantum[1][0] == 3)
            border = get_image_color (d, bild_nr, image, d_a, width, ho_pos, br_pos);
        else
            border = get_image (d, bild_nr, image, d_a, width, ho_pos, br_pos);

        if  (border)
            fprintf (stderr, "\nWarning: cut image patch is outside of original image (border=%d)!\n", border);

        /**b/w images NOT YET TESTED !!!**/
        if  (cmd->quantum[1][0] == 1) {

            /**OLD! accept only if minimal intensity differences compared to minimal pixel range ... OLD!
            double bild_min = 255.0, bild_max = 0.0;
            for (j = 0; j < d_a * width; ++j) {
                bild_min = image[j] < bild_min ? image[j] : bild_min;
                bild_max = image[j] > bild_max ? image[j] : bild_max;
            }
            if  ((bild_max - bild_min) >= cmd->quantum[0][0] * (d->Bild_max[bild_nr] - d->Bild_min[bild_nr]))
                badimage = 0;
            **/

            /**like below for color but here using "Red" as the only, b/w channel**/
            double meanR = 0.0;
            for (j = 0; j < height * width; ++j)
                meanR += image[j];
            meanR /= (double)(height * width);

            double sig_R = 0.0;
            for (j = 0; j < height * width; ++j)
                sig_R += (image[j]-meanR) * (image[j]-meanR);
            sig_R /= (double)(height * width);
            sig_R = sqrt(sig_R);

            if  (sig_R > cmd->quantum[0][0])
                badimage = 0;
        }

        /**color: accept only patches where the variance in each color channel is larger than q[0][0]**/
        if  (cmd->quantum[1][0] == 3) {
            double meanR = 0.0, meanG = 0.0, meanB = 0.0;
            for (j = 0; j < height * width; ++j)
                meanR += image[j];
            for (j = height * width; j < 2 * height * width; ++j)
                meanG += image[j];
            for (j = 2 * height * width; j < 3 * height * width; ++j)
                meanB += image[j];
            meanR /= (double)(height * width);
            meanG /= (double)(height * width);
            meanB /= (double)(height * width);

            double sig_R = 0.0, sig_G = 0.0, sig_B = 0.0;
            for (j = 0; j < height * width; ++j)
                sig_R += (image[j]-meanR) * (image[j]-meanR);
            for (j = height * width; j < 2 * height * width; ++j)
                sig_G += (image[j]-meanG) * (image[j]-meanG);
            for (j = 2 * height * width; j < 3 * height * width; ++j)
                sig_B += (image[j]-meanB) * (image[j]-meanB);
            sig_R /= (double)(height * width);
            sig_G /= (double)(height * width);
            sig_B /= (double)(height * width);
            sig_R = sqrt(sig_R);
            sig_G = sqrt(sig_G);
            sig_B = sqrt(sig_B);

            if  ((sig_R > cmd->quantum[0][0]) && (sig_G > cmd->quantum[0][0]) && (sig_B > cmd->quantum[0][0]))
                badimage = 0;
	}

        if  (badimage)
            fprintf (stderr, "r");

    } while (badimage);


    /**copy to all times**/
    if  (cmd->quantum[1][1] == 2) {

        sub_mean_vector (image, d_a * width);

        for (ct_t = begin; ct_t < end; ++ct_t)
            for (int a = 0; a < d_a; ++a)
                for (int b = 0; b < d_b; ++b)
                    cmd->S_target[ct_t][a * d_b + b] = 0.0;

        for (ct_t = begin; ct_t < end; ++ct_t)
            for (int a = 0; a < d_a; ++a)
                for (int b = 0; b < width; ++b)
                    if  (image[a * width + b] >= 0)
                        cmd->S_target[ct_t][a * d_b + b] = image[a * width + b];
                    else
                        cmd->S_target[ct_t][a * d_b + b + width] = -image[a * width + b];
    } else {
        for (ct_t = begin; ct_t < end; ++ct_t)
            for (j = 0; j < d_a * width; ++j)
                cmd->S_target[ct_t][j] = image[j];
    }

    return (DOUBLE)(0);
}


/******************************** cut_patch_at *******************************/
/* From cut_image. br_pos and ho_pos are read from from S_from1              */
/* cmd->S_target is set to an image part of size A[area].d_a * A[area].d_b.  */
/* Calls get_image(_color).                  Value range is still 0 .. 255.  */
/* q[0][0] == 3: color, else black&white.                                    */

/* Use  import_images  first!                                                */

DOUBLE cut_image_at (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {

    int j, ct_t;
    static int firsttime = 1;
    static DOUBLE *image;
    int bild_nr, ho_pos = 0, br_pos = 0;
    int border;
    DATA *d = (DATA *)cmd->pointers[0]->data;
    int inarea    = cmd->n_from1[0];

    if  (d->anzahl == 0)
        fprintf (stderr, "\nthere are no images in sight!\n");

    if  (cmd->anz_quantums != 1)
        fprintf (stderr, "\n1 parameter for cut_image_at, please!\n");
    if  (cmd->anz_quant[0] != 1)
        fprintf (stderr, "\n1 parameter for cut_image_at, q[0], please!\n");
    if  (A[inarea].d_n != 2)
        fprintf (stderr, "\ncut_image_at wants inarea with 2 units!\n");

    const int d_a = A[cmd->area].d_a;
    const int d_b = A[cmd->area].d_b;

    int height = d_a;
    int width  = d_b;
    if  (cmd->quantum[0][0] == 3)   /**color**/
        height = d_a / 3;

    /**allocate mem to store cut-out image**/
    if  (firsttime) {
        image = d_vector (d_a * width);  /**d_a contains 3 colors**/
        firsttime = 0;
    }

    /**copy to all times**/
    for (ct_t = begin; ct_t < end; ++ct_t) {

        /**choose random pic and patch; should also work for patch_size == image_size**/
        bild_nr = (int)(drand48 () * (double)(d->anzahl));

        /**position**/
        br_pos  = (d->Br_b[bild_nr] - width)  / 2 + (int)cmd->S_from1[0][ct_t][1];
        ho_pos  = (d->Ho_a[bild_nr] - height) / 2 + (int)cmd->S_from1[0][ct_t][0];

        /**cut chosen patch out**/
        if  (cmd->quantum[0][0] == 3)
            border = get_image_color (d, bild_nr, image, d_a, width, ho_pos, br_pos);
        else
            border = get_image (d, bild_nr, image, d_a, width, ho_pos, br_pos);

        if  (border)
            fprintf (stderr, "\nWarning: cut image patch is outside of original image (border=%d)!\n", border);

        for (j = 0; j < d_a * width; ++j)
            cmd->S_target[ct_t][j] = image[j];
    }

    return (DOUBLE)(0);
}


/******************************** cut_image_pantilt **************************/
/* As cut_image: No preprocess/edges/shift/memory. Here also: No bad_image.  */
/* cmd->S_target is initialized with an image part of size z->d_a * z->d_b.  */
/* Calls get_image(_color).                  Value range is still 0 .. 255.  */
/* q[0][0] == -1 new image, else same as before.                             */
/* q[1][0] == 1: color, else black&white.                                    */
/* q[2][0],[1] == ho_pos, br_pos; == -1 random; == -2 middle; == -3 "move".  */
/* "move" defined by first 2 pixels of S_target! Uniquely dirty mechanism!!  */
/* Time t="end" of relax counts!                                             */

DOUBLE cut_image_pantilt (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {

    int j, ct_t, d_r;
    static int firsttime = 1;
    static DOUBLE *image;
    int border;
    static int bild_nr = 0;
    static int ho_pos  = 0;
    static int br_pos  = 0;
    int height = A[cmd->area].d_a;

    DATA *d = (DATA *)cmd->pointers[0]->data;

    if  (d->anzahl == 0)
        fprintf (stderr, "\nthere are no images in sight!\n");
    if  (cmd->anz_quantums != 3)
        fprintf (stderr, "\n3 parameters for cut_image_pantilt, please!\n");
    if  (cmd->anz_quant[2] != 2)
        fprintf (stderr, "\n2 parameters for cut_image_pantilt, q[2], please!\n");

    d_r = A[cmd->area].d_a * A[cmd->area].d_b; /**true also for color**/

    /**allocate mem to store image**/
    if  (firsttime) {
        image = d_vector (d_r);
        firsttime = 0;
    }

    if  (cmd->quantum[1][0] == 1)   /**color**/
        height = A[cmd->area].d_a / 3;

    /**new random pic**/
    if  (cmd->quantum[0][0] == -1)
        bild_nr = (int)(drand48 () * (double)(d->anzahl));

    /**new random position**/
    if  (cmd->quantum[2][0] == -1)
        ho_pos  = (int)(drand48 () * (double)(d->Ho_a[bild_nr] - height));    /**d_a || hoehe:  langsamer index**/
    if  (cmd->quantum[2][1] == -1)
        br_pos  = (int)(drand48 () * (double)(d->Br_b[bild_nr] - A[cmd->area].d_b));    /**d_b || breite: schneller index**/

    /**new defined position**/
    if  (cmd->quantum[2][0] > 0)
        ho_pos  = (int)(cmd->quantum[2][1]);
    if  (cmd->quantum[2][1] > 0)
        br_pos  = (int)(cmd->quantum[2][0]);

    /**new middle position**/
    if  (cmd->quantum[2][0] == -2)
        ho_pos  = (d->Ho_a[bild_nr] - height) / 2;
    if  (cmd->quantum[2][1] == -2)
        br_pos  = (d->Br_b[bild_nr] - A[cmd->area].d_b) / 2;

    /**Uniquely dirty mechanism(!!): new pos according to S_target(!!) at time "end" of relaxation.**/
    /**Position may be out of data image! (border warning below, then)**/
    if  (cmd->quantum[2][0] == -3) {
        ho_pos  += (int)(cmd->S_target[end - 1][0]);
fprintf (stderr, "\ncut_image_pantilt: moving h+=%.2f", cmd->S_target[end - 1][0]);
    }
    if  (cmd->quantum[2][1] == -3) {
        br_pos  += (int)(cmd->S_target[end - 1][1]);
fprintf (stderr, " b+=+=%.2f\n", cmd->S_target[end - 1][1]);
    }

    /**old position stays if cmd->quantum[2][0]/[1] == 0**/

    /**cut chosen patch out**/
    if  (cmd->quantum[1][0] == 0)
        border = get_image (d, bild_nr, image, A[cmd->area].d_a, A[cmd->area].d_b, ho_pos, br_pos);
    else
        border = get_image_color (d, bild_nr, image, A[cmd->area].d_a, A[cmd->area].d_b, ho_pos, br_pos);
    if  (border)
        fprintf (stderr, "\nWarning: cut image patch left the original image (border=%d)!\n", border);

    /**copy to all times**/
    for (ct_t = begin; ct_t < end; ++ct_t)
        for (j = 0; j < d_r; ++j)
            cmd->S_target[ct_t][j] = image[j];

    return (DOUBLE)(0);
}


/******************************** init_orange ********************************/
/* cmd->S_target initialized with orange, zero surround. Range: 0 .. 255.    */
/* q[0][0] == 1: new rand location and -- if q[2[0]=3 -- new fruit           */
/* q[1][0] == diameter of orange / sigma of Gaussian.                        */
/* q[2][0] == 0: Gaussian (greyscale);                                       */
/*            1: orange (RGB format); 2:apple; 3: orange or apple (50% prob).*/

/* ATTENTION:  q[2][0] changed (was before "2" for Gaussian, now "0")        */
/* ATTENTION2: apple color manipulated                                       */
/* ATTENTION3: a fruit must be selected/placed first BEFORE language area set*/

DOUBLE init_orange (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {

    int a, b, ct_t, d_r, height;
    static double ho_center;
    static double br_center;
    static int fruit_select;
    double radius = (double)cmd->quantum[1][0];
    int proto1, proto2;
    double colmix, red = 0.0, green = 0.0, blue = 0.0;

    double O_col[39+5][3] = {/*R G    B        orange-likeliness              from file   **/

                       {241, 109,  24},   /**prototype**/          /**only the best from below**/
                       {230,  88,  16},   /**nice prototype**/
                       {227,  86,  14},   /**nice medium**/
                       {221,  92,  26},   /**nice prototype**/
                       {234, 146,  36},   /**prototype**/

                       {253, 133,  80},   /**brighter                       fruit_02.jpg**/
                       {247, 122,  30},   /**brighter**/
                       {244, 126,  39},   /**brighter**/
                       {241, 109,  24},   /**prototype**/
                       {238, 144, 106},   /**bright shiny spot**/
                       {230,  88,  16},   /**nice prototype**/
                       {227,  86,  14},   /**nice medium**/
                       {221,  92,  26},   /**nice prototype**/
                       {206,  69,  14},   /**darker**/
                       {203,  59,  25},   /**darker, reddish**/
                       {175,  66,  33},   /**darker**/
                       {161,  39,  34},   /**dark red-brown**/
                       {254, 182,  72},   /**pale                           orangegrapefruitgiftbox.jpg**/
                       {242, 147,  55},   /**brighter, prototype region**/
                       {240, 175,   0},   /**grapefruit?**/
                       {234, 146,  36},   /**prototype**/
                       {199,  86,   8},   /**darker region**/
                       {195, 124,  20},   /**darker pale**/

                       { 89,  48,  37},   /**shadow      snap1.ppm**/
                       {136,  57,  27},   /**  to  **/
                       {203, 102,  36},   /**  ..  **/
                       {238, 130,  41},   /**bright**/
                       {115,  54,  34},   /**shadow      snap2.ppm**/
                       {129,  57,  12},   /**  to  **/
                       {180,  85,  35},   /**  ..  **/
                       {231, 124,  53},   /**bright**/
                       {110,  47,  30},   /**shadow      snap3.ppm**/
                       {144,  66,  28},   /**  to  **/
                       {158,  77,  40},   /**  ..  **/
                       {202,  91,  46},   /**bright (not brightest)**/
                       { 85,  40,  19},   /**shadow      snap6.ppm**/
                       {153, 105,  39},   /**  to  **/
                       {198, 151,  60},   /**bright**/
                       { 77,  35,  15},   /**shadow      snap7.ppm**/
                       {125,  64,  35},   /**  to  **/
                       {194, 126,  69},   /**bright**/
                       {106,  56,  28},   /**shadow, other orange**/
                       {170, 119,  70},   /**  to  **/
                       {222, 177, 128}    /**bright**/
                      };

    double A_col[19][3] = {/*R G    B        apple-likeliness             from file   **/
                       {102, 108,  34},   /**medium**/
                       {112, 131,  76},   /**medium**/
                       {111, 120,  41},   /****/
                       {110, 118,  43},   /****/
                       {111, 135,  57},   /****/
                       {107, 117,  28},   /****/

                       {220, 224, 104},   /****/
                       {215, 217,  94},   /****/
                       {213, 217,  94},   /****/

                       {168, 170,  63},   /**shade**/
                       {173, 175,  68},   /****/


                       {182, 173,  52},   /****/
                       {179, 178,  86},    /****/

                       {234, 245, 215},   /**shiny spot**/
                       {189, 204, 175},   /**shiny spot**/
                       {234, 241, 171},   /**brightest spot      apple_green_large.jpg**/
                       {139, 137,  28},   /**yellowish spot**/
                       { 97,  84,  39},   /**dark                         fruit_02.jpg**/
                       { 90,  83,  37}    /**dark**/
                      };

    if  (cmd->anz_quantums != 3)
        fprintf (stderr, "\n3 parameters for init_orange, please!\n");

    d_r = A[cmd->area].d_a * A[cmd->area].d_b; /**true also for color**/

    if  (cmd->quantum[2][0] == 0)
        height = A[cmd->area].d_a;       /**in this case only used for choose location, here below**/
    else
        height = A[cmd->area].d_a / 3;   /**because color**/

    /**choose new random location and -- if random -- fruit**/
    if  (cmd->quantum[0][0] == 1.0) {
        ho_center = radius + drand48 () * ((double)(height) - radius - radius);    /**d_a || hoehe:  langsamer index**/
        br_center = radius + drand48 () * ((double)(A[cmd->area].d_b) - radius - radius);    /**d_b || breite: schneller index**/

        /**randomly orange or apple**/
        if  (cmd->quantum[2][0] == 3.0) {
            if  (drand48() < 0.5)
                fruit_select = 1;
            else
                fruit_select = 2;
        }
    }

    /**choose location according to input area <-- CAN'T BE IMPLEMANTED BECAUSE INAREA NOT IN ARGUMENTS => MOVE THIS TO A TOTAL_FUNCTION
    if  (cmd->quantum[0][0] == 2.0) {
        int inarea = cmd->n_from1[0];
        ho_center = radius + drand48 () * ((double)(height) - radius - radius);
        br_center = radius + drand48 () * ((double)(A[cmd->area].d_b) - radius - radius);
    }
    **/

    /**draw orange or apple**/
    if  ((cmd->quantum[2][0] == 1.0) || (cmd->quantum[2][0] == 2.0) || (cmd->quantum[2][0] == 3.0)) {

        /**choose fruit**/

        if  (cmd->quantum[2][0] == 1.0)
            fruit_select = 1;

        if  (cmd->quantum[2][0] == 2.0)
            fruit_select = 2;

        /**choose color**/

        colmix = drand48();

        /**orange**/
        if  (fruit_select == 1) {
            proto1 = (int)(drand48() * 5); /**only the 5 nicest of 39 + 5**/
            proto2 = (int)(drand48() * 5);
            red    = colmix * O_col[proto1][0] + (1.0 - colmix) * O_col[proto2][0];
            green  = colmix * O_col[proto1][1] + (1.0 - colmix) * O_col[proto2][1];
            blue   = colmix * O_col[proto1][2] + (1.0 - colmix) * O_col[proto2][2];
        }
        /**apple**/
        if  (fruit_select == 2) {
            proto1 = (int)(drand48() * 11); /**only the 11 nicest of 19**/
            proto2 = (int)(drand48() * 11);
            red    = colmix * A_col[proto1][0] + (1.0 - colmix) * A_col[proto2][0];
            green  = colmix * A_col[proto1][1] + (1.0 - colmix) * A_col[proto2][1];
            blue   = colmix * A_col[proto1][2] + (1.0 - colmix) * A_col[proto2][2];

            red   -= 10.0;    /**!!!**/
            green += 10.0;    /**!!!**/
        }

        /**draw orange/apple disc (for all times) and maintain other background**/
        for (ct_t = begin; ct_t < end; ++ct_t) {
            for (a = 0; a < height; ++a)
            for (b = 0; b < A[cmd->area].d_b; ++b)
            if  (((double)(a)-ho_center)*((double)(a)-ho_center) + ((double)(b)-br_center)*((double)(b)-br_center) < (radius*radius))
                cmd->S_target[ct_t][a*A[cmd->area].d_b + b] = red;
            else
                cmd->S_target[ct_t][a*A[cmd->area].d_b + b] = cmd->S_from1[0][ct_t][a*A[cmd->area].d_b + b];

            for (a = height; a < 2 * height; ++a)
            for (b = 0; b < A[cmd->area].d_b; ++b)
            if  (((double)(a-height)-ho_center)*((double)(a-height)-ho_center) + ((double)(b)-br_center)*((double)(b)-br_center)
                < (radius*radius))
                cmd->S_target[ct_t][a*A[cmd->area].d_b + b] = green;
            else
                cmd->S_target[ct_t][a*A[cmd->area].d_b + b] = cmd->S_from1[0][ct_t][a*A[cmd->area].d_b + b];

            for (a = 2 * height; a < 3 * height; ++a)
            for (b = 0; b < A[cmd->area].d_b; ++b)
            if  (((double)(a-2*height)-ho_center)*((double)(a-2*height)-ho_center) + ((double)(b)-br_center)*((double)(b)-br_center)
                < (radius*radius))
                cmd->S_target[ct_t][a*A[cmd->area].d_b + b] = blue;
            else
                cmd->S_target[ct_t][a*A[cmd->area].d_b + b] = cmd->S_from1[0][ct_t][a*A[cmd->area].d_b + b];
        }

    }


    if  (cmd->quantum[2][0] == 0.0) {  /**draw greyscale(!) Gaussian**/

        for (ct_t = begin; ct_t < end; ++ct_t) {
            for (a = 0; a < A[cmd->area].d_a; ++a)
            for (b = 0; b < A[cmd->area].d_b; ++b)
                    cmd->S_target[ct_t][a*A[cmd->area].d_b + b]
                = 1.0 * exp (- (((double)a-ho_center)*((double)a-ho_center)+((double)b-br_center)*((double)b-br_center))
                            /  (2.0 * radius * radius));
        }
    }


    if  (cmd->quantum[2][0] == 4.0) { /**language area: one row or column ON**/

        for (ct_t = begin; ct_t < end; ++ct_t) {
            for (a = 0; a < A[cmd->area].d_a; ++a)
                for (b = 0; b < A[cmd->area].d_b; ++b)
                    cmd->S_target[ct_t][a*A[cmd->area].d_b + b] = 0.0;
            if  (A[cmd->area].d_a == 2)
                for (b = 0; b < A[cmd->area].d_b; ++b)
                    cmd->S_target[ct_t][(fruit_select-1)*A[cmd->area].d_b + b] = 1.0;
            if  (A[cmd->area].d_b == 2)
                for (a = 0; a < A[cmd->area].d_a; ++a)
                    cmd->S_target[ct_t][a*A[cmd->area].d_b + (fruit_select-1)] = 1.0;
        }
    }

    return (DOUBLE)(0);
}


/******************************** image_color_blob ***************************/
/* Select region from RGB S_from1 image based on color. Range: 0 .. 255.     */
/* Gaussian placed on CM on ANOTHER (target) area (without RGB layers).      */
/* q[0][0]>0: sigma for Gaussian; <=0: set blob region pixels to 1 else 0.   */
/* q[1][0/1/2]: lower R/G/B limits for color interval                        */
/* q[2][0/1/2]: upper R/G/B limits for color interval                        */
/* orange: 245+182+61, 255+212+87 divided by 2 yields 122+91+30, 127+107+44   */

DOUBLE image_color_blob (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {

    int a, b, ct_t;
    double CM_x = 0.0;
    double CM_y = 0.0;
    int count = 0;
    static int failures = 0;

    double sigma  = cmd->quantum[0][0];
    int inarea    = cmd->n_from1[0];
    int height_in = A[inarea].d_a / 3;   /**because color**/
    double R_low  = cmd->quantum[1][0];
    double G_low  = cmd->quantum[1][1];
    double B_low  = cmd->quantum[1][2];
    double R_high = cmd->quantum[2][0];
    double G_high = cmd->quantum[2][1];
    double B_high = cmd->quantum[2][2];

    if  (end == 0)
        end = begin + 1;    /**if used as total function, because end would then be 0**/

    /**init**/
    for (ct_t = begin; ct_t < end; ++ct_t)
    for (a = 0; a < A[cmd->area].d_a; ++a)
    for (b = 0; b < A[cmd->area].d_b; ++b)
        cmd->S_target[ct_t][a*A[cmd->area].d_b + b] = 0.0;

    for (ct_t = begin; ct_t < end; ++ct_t) {
        for (a = 0; a < height_in; ++a)
        for (b = 0; b < A[inarea].d_b; ++b) {

            int match = 0;

            double red = cmd->S_from1[0][ct_t][a * A[inarea].d_b + b];
            if  ((red >= R_low) && (red <= R_high))
                match += 1;

            double green = cmd->S_from1[0][ct_t][(height_in + a) * A[inarea].d_b + b];
            if  ((green >= G_low) && (green <= G_high))
                match += 1;

            double blue = cmd->S_from1[0][ct_t][(2 * height_in + a) * A[inarea].d_b + b];
            if  ((blue >= B_low) && (blue <= B_high))
                match += 1;

            if  (match == 3) {
                CM_x += a;
                CM_y += b;
                count += 1;

                if  (cmd->quantum[0][0] <= 0.0) {
                    int x_pos = (int)((float)a / (float)height_in * (float)A[cmd->area].d_a);
                    int y_pos = (int)((float)b / (float)(A[inarea].d_b) * (float)A[cmd->area].d_b);
                    cmd->S_target[ct_t][x_pos*A[cmd->area].d_b + y_pos] = 1.0;
                }
            }
        }
    }

    if  (count == 0) {
        failures += 1;
        fprintf (stderr, "\nimage_color_blob: %d-th time no blob  ", failures);
    }

    if  (cmd->quantum[0][0] > 0.0) {

        if  (count > 0) {

            CM_x /= (double)count;
            CM_y /= (double)count;

            double ho_center = CM_x * A[cmd->area].d_a / (double)height_in;
            double br_center = CM_y * A[cmd->area].d_b / (double)A[inarea].d_b;

            for (ct_t = begin; ct_t < end; ++ct_t)
            for (a = 0; a < A[cmd->area].d_a; ++a)
            for (b = 0; b < A[cmd->area].d_b; ++b) {
                double diffA = (double)a-ho_center;
                double diffB = (double)b-br_center;

                cmd->S_target[ct_t][a*A[cmd->area].d_b + b]
                = 1.0 * exp (- 0.5 * (diffA*diffA + diffB*diffB) / (sigma * sigma));
            }
        }
    }

    return (DOUBLE)(0);
}




/*****************************************************************************/
/*****************************************************************************/
/*****************************************************************************/




/************************** choose_active_lines ******************************/
/**if number_active == -1: stochasticalle some of vec_active are switched on**/
/**                        # active cells is returned                       **/
/**if number_active > -1:              as many of vec_active are switched on**/
/**                        0 is returned                                    **/
int choose_active_lines (int number_active, int *vec_active, int dim) {
    int i_n;

    for (i_n = 0; i_n < dim; ++i_n)
        vec_active[i_n] = 0;

    if  (number_active == -1) {
        number_active = 0;
        for (i_n = 0; i_n < dim; ++i_n)
            if  (drand48() < 1.0 / (double)(dim)) {
                                            /**average: one stimulus per eye**/
                                                   /**2 out of 30 neurons on**/
                vec_active[i_n] = 1;
                number_active += 1;
            }
    } else {
        while (number_active > 0) {
            i_n = (int)(drand48() * (double)(dim));
            if  (vec_active[i_n] == 0) {
                vec_active[i_n] = 1;
                number_active -= 1;
            }
        }
    }

    return (number_active);
}


/*
int R_init_offset (double *vec, int d_a, int d_b, DATA *d, PARAMS *g) {

    int i_r, i_n, i_ns, sel_m, number_active = 0;

    int n  = g->d_b[0];       ** instances (upper pixel) **
    int m  = g->d_y[0];       ** offsets (classes) **
                              ** n * m = 15 possible events for each eye **
    int ns = g->d_a[0] / 2;   ** for training, harder if ns > 1 **

    static int firsttime = 1;
    static int *vec_active;

    if  (firsttime) {
        vec_active = i_vector (n);
        firsttime = 0;
    }

    if  (g->anz_r > 2)             fprintf (stderr, "\nanz_r?\n");
    if  (g->anz_t > g->anz_r + 1)  fprintf (stderr, "\nanz_t?\n");

    if  (g->d_x[0] != n)           fprintf (stderr, "\nd_b .. d_x?\n");
    if  (g->d_a[0] % 2 != 0)       fprintf (stderr, "\nd_a?\n");


    **choose exactly one class**
    sel_m = (int)(drand48() * (double)(m));

    **middle layer: two neurons of 2 * n*m are on**
    **bottom layer: 2 * aux3 neurons of ns * 2 * n neurons are on**
    **look at update dynamics -- stabilize somehow (?)**


    for (i_r = 0; i_r < g->anz_r; ++i_r) {                           **eyes**
        if  (i_r == 0)
            number_active = choose_active_lines (-1, vec_active, n);
        else
            number_active = choose_active_lines (number_active, vec_active, n);
        for (i_n = 0; i_n < n; ++i_n) {                         **instances**
            if  (vec_active[i_n] == 1) {
                for (i_ns = 0; i_ns < ns; ++i_ns) {            **along line**
                    if  (drand48() < g->aux3[0] / (double)(ns)) **old staff**
                        vec[i_r * (2 * ns * n) + i_ns * n + i_n] = 1.0;
                     **    count eye(s)  **  lines  ** instance**
                    if  (drand48() < g->aux3[0] / (double)(ns)) **old staff**
                        vec[i_r * (2 * ns * n) + i_ns * n
                                   + ns * n + (i_n + sel_m) % n] = 1.0;
                                **half-eye ** class offset  **
                }
            }
        }
    }

    **"noise" according to zero-probabilities kurt_t**
    for (i_r = 0; i_r < g->anz_r; ++i_r)                             **eyes**
        for (i_n = 0; i_n < n; ++i_n)                           **instances**
            for (sel_m = 0; sel_m < m; ++sel_m)                   **classes**
                if  (drand48() < 1.0 / g->kurt_t[0])
                    for (i_ns = 0; i_ns < ns; ++i_ns) {        **along line**
                        vec[i_r * (2 * ns * n) + i_ns * n + i_n] = 1.0;
                        vec[i_r * (2 * ns * n) + i_ns * n
                                   + ns * n + (i_n + sel_m) % n] = 1.0;
                    }
    does not make a difference ...
    **

    return (0);
}
*/

/*
int R_init (double *inputvec, PARAMS *g, DATA *d) {
    int anfj, ari;

    anfj = 0;

    if  (g->ch_look == 1) { **natural images (on 1st retina)**
        R_init_data (inputvec, g->d_a[0], g->d_b[0], d);
        sub_mean_vector (inputvec, g->d_a[0] * g->d_b[0]);
        ** do not spherize here, because korners become large **
        ** spherize_vector (inputvec, g->d_a[0] * g->d_b[0]); **

        if  (g->anz_r == 2) {
            anfj += g->d_a[0] * g->d_b[0];
            R_init_gabor3 (inputvec + anfj, g->d_a[1], g->d_b[1], d);
            sub_mean_vector (inputvec + anfj, g->d_a[1] * g->d_b[1]);
            spherize_vector (inputvec + anfj, g->d_a[1] * g->d_b[1]);
        }

        if  (g->anz_r > 2)
            fprintf (stderr, "\ntoo many retinae\n");
    }

    if  (g->ch_look == 2) { **different statistics on 1st / 2nd retina**
        R_init_gabor3 (inputvec, g->d_a[0], g->d_b[0], d);
        sub_mean_vector (inputvec, g->d_a[0] * g->d_b[0]);
        spherize_vector (inputvec, g->d_a[0] * g->d_b[0]);

        anfj += g->d_a[0] * g->d_b[0];
        for (ari = 1; ari < g->anz_r; ++ari) {
            R_init_gabor4 (inputvec + anfj, g->d_a[ari], g->d_b[ari], d);
            sub_mean_vector (inputvec + anfj, g->d_a[ari] * g->d_b[ari]);
            spherize_vector (inputvec + anfj, g->d_a[ari] * g->d_b[ari]);
            anfj += g->d_a[ari] * g->d_b[ari];
        }
    }

    if  (g->ch_look == 3) { **2nd image somehow dependent on 1st image**
        if  ((g->anz_r != 2) && (g->anz_r != 4)) {
            fprintf (stderr, "\nstereo vision with %d eye(s)?\n", g->anz_r);
            exit (13);
        }

        anfj = 0;
        for (ari = 0; ari < g->anz_r; ++ari) {
            if  ((g->d_a[ari] != g->d_a[0]) || (g->d_b[ari] != g->d_b[0]))
                fprintf (stderr, "\nstereo vision with unequal eyes?\n");

            R_init_stereo (inputvec + anfj, g->d_a[ari], g->d_b[ari], d, g);
            sub_mean_vector (inputvec + anfj, g->d_a[ari] * g->d_b[ari]);
            spherize_vector (inputvec + anfj, g->d_a[ari] * g->d_b[ari]);
            anfj += g->d_a[ari] * g->d_b[ari];
        }
    }

    if  (g->ch_look == 4) { **one pixel on**
        int i, ch;

        if  (g->anz_r > 1) {
            fprintf (stderr, "\n%d are too many eyes\n", g->anz_r);
            exit (13);
        }

        for (i = 0; i < g->d_a[0] * g->d_b[0]; ++i)
            inputvec[i] = 0.0;
        ch = (int)(drand48() * g->d_a[0] * g->d_b[0]);
        inputvec[ch] = 1.0;
    }

    if  (g->ch_look == 5) { **3 pixels in a row**
        if  (g->anz_r > 1) {
            fprintf (stderr, "\n%d are too many eyes\n", g->anz_r);
            exit (13);
        }

        R_init_3row (inputvec, g->d_a[0], g->d_b[0], d, g);
    }

    if  (g->ch_look == 6) **instances of a chosen class**
        R_init_offset (inputvec, g->d_a[0], g->d_b[0], d, g);


    if  ((g->ch_look < 1) || (g->ch_look > 6))
        fprintf (stderr, "\nch_look out of limits\n");

    return (0);
}
*/



/******************************** init_image_cosinus *************************/
/* Wrapper function for total_cosinus. Modify as convenient.                 */
/* Versatile function, hence "int end" may not be used or zero or something. */

DOUBLE init_image_cosinus (PARAMS *g, AREA *A, COMMAND *cmd, int begin, int end) {

    static COMMAND *give_cmd;
    static int firsttime = 1;

    if  (firsttime) {
        fprintf (stderr, "\nfirst time init_image_cosinus   ");

        give_cmd = (COMMAND *)malloc (sizeof (COMMAND));
        give_cmd->quantum = (DOUBLE **)malloc (10 * sizeof (DOUBLE *));
        for (int i = 0; i < 10; ++i)
            give_cmd->quantum[i] = (DOUBLE *)malloc (2 * sizeof (DOUBLE));
        firsttime = 0;
    }

    give_cmd->S_target = cmd->S_target;
    give_cmd->area     = cmd->area;

/*  -- this is from total_cosinus --
    cm_x          = cmd->quantum[0][0];
    cm_y          = cmd->quantum[0][1];
    sq_diameter   = cmd->quantum[1][0] * cmd->quantum[1][0];
    amplitude_in  = cmd->quantum[2][0];
    amplitude_out = cmd->quantum[2][1];
    frequency_in  = cmd->quantum[3][0];
    frequency_out = cmd->quantum[3][1];
    phase_in      = cmd->quantum[4][0];
    phase_out     = cmd->quantum[4][1];
    angle_in      = cmd->quantum[5][0];
    angle_out     = cmd->quantum[5][1];
    velocity_in   = cmd->quantum[6][0];
    velocity_out  = cmd->quantum[6][1];
    sphere        = (int)(cmd->quantum[7][0]);  --not implemented--
*/

  /*cm_x         */   give_cmd->quantum[0][0] = 0.0;
  /*cm_y         */   give_cmd->quantum[0][1] = 0.0;
  /*sq_diameter  */   give_cmd->quantum[1][0] = 100;
  /*amplitude_in */   give_cmd->quantum[2][0] = cmd->quantum[0][0] + (cmd->quantum[0][1] - cmd->quantum[0][0]) * drand48 ();
  /*amplitude_out*/   give_cmd->quantum[2][1] = 0.0;
  /*frequency_in */   give_cmd->quantum[3][0] = cmd->quantum[1][0] + (cmd->quantum[1][1] - cmd->quantum[1][0]) * drand48 ();
  /*frequency_out*/   give_cmd->quantum[3][1] = 0.0;
  /*phase_in     */   give_cmd->quantum[4][0] = cmd->quantum[2][0] + (cmd->quantum[2][1] - cmd->quantum[2][0]) * drand48 () * 2.0 * M_PI;
  /*phase_out    */   give_cmd->quantum[4][1] = 0.0;
  /*angle_in     */   give_cmd->quantum[5][0] = 0.5 * M_PI; /* 2.0 * M_PI * drand48 (); */
  /*angle_out    */   give_cmd->quantum[5][1] = 0.0;
  /*velocity_in  */   give_cmd->quantum[6][0] = 0.0;
  /*velocity_out */   give_cmd->quantum[6][1] = 0.0;
  /*sphere       */   give_cmd->quantum[7][0] = 0.0;
  /*admin inform */   give_cmd->anz_quantums  = 8;

    int ct_t = begin;
    do  {
        total_cosinus (g, A, give_cmd, ct_t, 0);
    } while (ct_t < end);

    return (DOUBLE)(0);
}
