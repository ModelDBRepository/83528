#if BTTV


// Source and Documentation:
// http://www.exploits.org/v4l/    -> The v4l mailing list archives
// http://www.linuxhq.com/kernel/v2.2/doc/video4linux/API.html.html

// TCap v0.1
// (C) Thomas Hargrove
// http://toonarchive.com/tcap

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include <stdlib.h>
//#include <vga.h>
//#include <vgagl.h>
#include <stdio.h>
#include <fcntl.h>
#include <sys/ioctl.h>
#include <sys/mman.h>
#include <errno.h>
#include <time.h>

#include "iter.h"
#include "series.h"
#include "utils.h"
#include "data.h"
#include "data.videodev.h"   /**original file: videodev.h**/


// Only 640x480 and 320x240 work here
//const int 	mywidth = 64;
//const int 	myheight = 48;

static struct	video_capability  capability;
static int 	fd = -1;
static struct 	video_mbuf gb_buffers = { 2*0x151000, 0, {0,0x151000 }};
static char 	*map = NULL;
static struct 	video_mmap my_buf;

//GraphicsContext *physicalscreen;
//GraphicsContext *virtualscreen;
clock_t		cstart=0;
clock_t		clast=0;
float		fps;
int		frame_count=0;

// ------------------------------------------------------------
// This function copies the image from the capture card (tmap)
// and draws each pixel on the off screen buffer.  Then it
// copies the off screen buffer to the actual screen.
// ------------------------------------------------------------
// This function is quite slow, so if you need a higher fps
// just optimize this routine
// ------------------------------------------------------------
/*
**************************************************************

Modified the above comment to cop it to a file

*************************************************************
*/


void copytofile(char* tmap) {

/*   int x,y;
   unsigned char r,g,b;
   for(y=0; y<myheight; y++) {
	for(x=0; x<mywidth; x++) {
		r = *(tmap+2);
		g = *(tmap+1);
		b = *tmap;
		gl_setpixelrgb(x, y, r, g, b); 
		tmap+=3;
	}
   }
   gl_printf(10,5,"TCap v0.1");
   gl_printf(10,15,"Tom Hargrove");

   gl_printf(10,25,"Fps %0.1f", fps);
   if ( ((cstart=clock())-clast) / CLOCKS_PER_SEC > 1) {
	fps = (float) frame_count / ((float) (cstart-clast) / CLOCKS_PER_SEC);
	frame_count=0;
	clast=cstart;
   }
   frame_count++;

   gl_copyscreen(physicalscreen);	// This copies the image to the screen!
*/
/*
FILE *fp;
int x,y;
unsigned char r,g,b;

int maxr = 0, maxg = 0, maxb = 0;

  fp = fopen ("image2.pnm", "w");
  fprintf (fp, "P6\n");
  fprintf (fp, "%d %d\n", mywidth, myheight);
  fprintf (fp, "255\n");

  for(y=0; y<myheight; y++) {
	for(x=0; x<mywidth; x++) {
		r = *(tmap+2);
		g = *(tmap+1);
		b = *tmap;
                maxr = r > maxr ? r : maxr;
                maxg = g > maxg ? g : maxg;
                maxb = b > maxb ? b : maxb;
		fprintf (fp, "%c%c%c", r, g, b);
		tmap+=3;
	}
   }

  fclose (fp);

   frame_count++;

  fprintf (stderr, "maxr=%d  maxg=%d  maxb=%d  frame=%d\n",
                    maxr, maxg, maxb, frame_count);
*/
}


/******************************** bttv_image *********************************/
/* Fille d->Bilder, similar to import_images, but as functon pointer.        */
/* q[0][0] = -1 init, else refresh; reads 2 frames but discards 1st => fresh */
/* q[1][0] = d->Br_b[bild_nr]                                                */
/* q[1][1] = d->Ho_a[bild_nr]                                                */

void bttv_image (AGENT *z, COMMAND *cmd, DATA *d, int begin, int end) {

   char *my_video_dev = "/dev/video0";
   const int bild_nr = 0;
   unsigned char r,g,b;
   int i, refresh;

   /**initialize**/
   if  (cmd->quantum[0][0] == -1) {

       // ----------------------------------------------------------------------
       // Get the v4l capture set up
       // ----------------------------------------------------------------------
       if  (-1 == (fd = open(my_video_dev, O_RDWR))) {
           printf("Error opening device: %s\n", my_video_dev);
           exit (1);
       }

       if  (-1 == ioctl(fd,VIDIOCGCAP,&capability)) {
           printf("Error: ioctl(fd,VIDIOCGCAP,&capability)\n");
           exit (1);
       }

       /*Set the close-on-exec flag as specified by FD_CLOEXEC*/
       fcntl(fd,F_SETFD,FD_CLOEXEC);

       if  (-1 == ioctl(fd,VIDIOCGMBUF,&gb_buffers)) {
           printf("Error: Error getting buffers, I think\n");
           exit (1);
       }

       map = mmap(0,gb_buffers.size,PROT_READ|PROT_WRITE,MAP_SHARED,fd,0); 
       if  (map == NULL) {
           printf("Error: Mmap returned NULL\n");
           exit (1);
       }

       // Set up out capture to use the correct resolution
       my_buf.width = cmd->quantum[1][0];
       my_buf.height = cmd->quantum[1][1];
       my_buf.format = VIDEO_PALETTE_RGB24;


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
      if  ((d->Bilder = (double **)malloc(d->anzahl*sizeof(double *))) == NULL)
          fprintf (stderr, "\nout of memory\n");

      /**(from importP5_img)**/
      d->Br_b[bild_nr]      = cmd->quantum[1][0];
      d->Ho_a[bild_nr]      = cmd->quantum[1][1];
      d->Bild_grau[bild_nr] = 255;
      d->Bilder[bild_nr]    = d_vector (d->Br_b[bild_nr] * d->Ho_a[bild_nr] * 3);
      d->Bild_min[bild_nr]  = 0.0;
      d->Bild_max[bild_nr]  = 255.0;
      d->Bild_mindiff[bild_nr] = 0.0;
   }

   // Tell the capture card to fill frame 0
/*
   my_buf.frame = 0;
   if (-1 == ioctl(fd, VIDIOCMCAPTURE, &my_buf)) { 
	printf("Error: Grabber chip can't sync (no station tuned in?)\n"); 
	goto err;
   }
*/   

   // --------------------------------------------------------------------------
   // This is the infinate loop
   // We basically:
   //	capture frame 1
   //   sync frame 0
   //   process frame 0
   //	capture frame 0 
   //   sync frame 1
   //   process frame 1
   // For more information, read the programming how-to that came with xawtv
   // --------------------------------------------------------------------------

   /**do the whole thing 2ce so that the newest frame will be taken**/
   for (refresh = 0; refresh < 2; ++refresh) {

       my_buf.frame = 0; /**changed this from 1 !!**/
       if  (-1 == ioctl(fd, VIDIOCMCAPTURE, &my_buf)) {
           printf("Error: Grabber chip can't sync (no station tuned in?)\n"); 
           printf("my_buf.width=%d my_buf.height=%d ... too small?\n", my_buf.width, my_buf.height); 
           exit (1);
       }
 
       my_buf.frame = 0;
       if  (-1 == ioctl(fd, VIDIOCSYNC, &my_buf.frame)) {
            printf("Error on sync!\n"); 
            exit (1);
       }

       for (i = 0; i < d->Br_b[bild_nr] * d->Ho_a[bild_nr] * 3; i += 3) {
           r = *(map + i + 2);
           g = *(map + i + 1);
           b = *(map + i);
           d->Bilder[bild_nr][i]   = (double)(r);
           d->Bilder[bild_nr][i+1] = (double)(g);
           d->Bilder[bild_nr][i+2] = (double)(b);
       }
   }

/*
	my_buf.frame = 0;
	if (-1 == ioctl(fd, VIDIOCMCAPTURE, &my_buf)) {
		printf("Error: Grabber chip can't sync (no station tuned in?)\n"); 
		goto err;
	}

	my_buf.frame = 1;
	if (-1 == ioctl(fd, VIDIOCSYNC, &my_buf.frame)) {
		printf("Error on sync!\n"); 
		goto err;
	}

	copytofile(map + gb_buffers.offsets[1]);
*/

   // Return screen to text mode.
   //vga_setmode(TEXT);
}


#endif

