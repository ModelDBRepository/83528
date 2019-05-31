#include "NNsimObject.h"
#include <stdio.h> 
#include <stdlib.h> 
#include <string> 
#include <iostream>


NNsimObject::NNsimObject(const int _analyze,
                         const std::string & _tmp_uname,
			 const std::string & _language_file,
			 const std::string & _name) :
  analyze(_analyze),
  tmp_uname((char *)_tmp_uname.c_str()),
  language_file((char *)_language_file.c_str()),
  name((char *)_name.c_str())
{
  char wovon[512];
  char tmp_uname_filename[512];
  char *NNSIM_ROOT = getenv ("NNSIM_ROOT");
  FILE *fp;

  std::cout << "Constructing NNsimObject." << std::endl;

  fprintf (stderr, "\nParameters: tmp_uname = %s", tmp_uname);
  fprintf (stderr, "\n        language_file = %s\n", language_file);
  fprintf (stderr, "\n                 name = %s\n", name);

  if  (NNSIM_ROOT == NULL)
      fprintf (stderr, "\n\nyou have to set the NNSIM_ROOT environment variable!\n\n");

  // pre-process the language file
  sprintf (wovon, "%s/prae.exe %s > %s/vonprae", NNSIM_ROOT, language_file, tmp_uname);
  fprintf (stderr, "\nExecuting system call (prae-processor): %s\n", wovon);
  system (wovon);

  // open the pre-processed file
  sprintf (tmp_uname_filename, "%s/vonprae", tmp_uname);
  fprintf (stderr, "opening to read: %s\n", tmp_uname_filename);
  fp = fopen (tmp_uname_filename, "r");
  if  (fp == NULL)
      fprintf (stderr, "\ncould not open %s!\n", tmp_uname_filename);

  {/**for debugging**/
   char kommentarzeile[256];
   fgets (kommentarzeile, 256, fp);
   fprintf (stderr, "NNsimObject.cpp reads first line: %s\n", kommentarzeile);
   fseek (fp, 0, SEEK_SET);
  }

    // import and allocate x, Z, s using yacc
    x = alloc_x (fp);
    fclose (fp);
    Z = alloc_z ();
    si = alloc_si ();

    printparams (x, Z, si, stderr);

    A = alloc_a (x, Z, si);
    init_a (x, Z, A, tmp_uname);

    choose_functions (x, si);

    choose_pointers (si, A);

    check_vehicle (si);

    if  ((d = (DATA *)malloc (sizeof (DATA))) == NULL)
        fprintf (stderr, "\nallocation failure for d\n");

    if  (x->data_import != NULL) {
        fprintf (stderr, "\nimport data: ");
        (*x->data_import)(x, d);
    } else {
        fprintf (stderr, "\nno data imported; hope you init data within series\n");
    }

fprintf(stderr,"\n----------------------------------------------------------");
}

void
NNsimObject::FirstSeries()
{
    if  (si->anz_se > 1) {
        if  ((si->se->elen == 1) && (si->se->ilen == 1))
            do_series (x, Z, A, si->se+0, d, analyze, tmp_uname);
        else
            fprintf (stderr, "\n\nFirstSeries Warning: elen or ilen aren't 1 |||\n\n");
    }
}

int
NNsimObject::MainSeries()
{
    /*
    if  (si->anz_se == 1)  do_series (x, Z, A, si->se+0, d, analyze, tmp_uname);
    if  (si->anz_se > 1)   do_series (x, Z, A, si->se+1, d, analyze, tmp_uname);
    */

    static int epoche = 0;
    SERIES *series = (si->anz_se == 1) ? si->se+0 : si->se+1;

    fprintf (stderr, "\nMainSeries: epoche=%d  ", epoche);

    if  (epoche < series->elen)
        do_epoche (x, Z, A, series, d, analyze, tmp_uname, epoche);

    epoche += 1;

    if  (epoche % series->elen == 0) {
        ThirdSeries ();
        epoche = 0;
        return 1;
    }

    return 0;
}

void
NNsimObject::ThirdSeries()
{
    if  (si->anz_se > 2) {
        do_series (x, Z, A, si->se+2, d, analyze, tmp_uname);
    }
}
