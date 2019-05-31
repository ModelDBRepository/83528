#include <stdio.h>
#include <stdlib.h>
#include <string.h>

void sort_longest (char **ch_main, char **ch_with, int anz) {

     char temp_main[512], temp_with[512];
     int i, j, longest;

     for (i = 0; i < anz; ++i) {

         /**find the index of the longest word**/
         longest = i;
         for (j = i; j < anz; ++j)
             longest = (strlen (ch_main[j]) > strlen (ch_main[longest])) ? j : longest;

         /**this is because of a former crash if argument line > 128 bytes**/
         if  (strlen (ch_with[longest]) == 0) {
             fprintf (stderr, "\nprae: strlen of one argument is zero!");
             fprintf (stderr, "\nargument is %s", ch_with[longest]);
	     if  (longest-1 >= 0)
                 fprintf (stderr, "\nprevious argument is %s\n", ch_with[longest-1]);
             exit (0);
	 }

         /**put this longest word to the front position i**/
         strcpy (temp_main, ch_main[longest]);
         strcpy (ch_main[longest], ch_main[i]);
         strcpy (ch_main[i], temp_main);

         strcpy (temp_with, ch_with[longest]);
         strcpy (ch_with[longest], ch_with[i]);
         strcpy (ch_with[i], temp_with);
     }
}


int main (int argc, char *argv[]) {

    FILE *fp;
    char zeile[512], *searchword[512], *changeword[512];
    int i, j, ct_ch, sel, set_line;
    int ct_replace = 0;
    int block_this = 0;

    fprintf (stderr, "\nprogram prae starting");

    if  (argc % 3 != 2) {
        fprintf (stderr, "\nWrong use of %s\n\n", argv[0]);
        exit (0);
    }
    fprintf (stderr, "\n%s filename [-set searchword changeword] [...]     used as:\n", argv[0]);
    for (i = 0; i < argc; ++i)
        fprintf (stderr, "%s ", argv[i]);
    fprintf (stderr, "\n");

    /**for command line options to become search- changeword's**/
    for (i = 2; i < argc - 2; i+= 3)
        if  (! strcmp (argv[i], "-set")) {

                    /**memory allocate only for searchword**/
                    searchword[ct_replace] = (char *)malloc (512 * sizeof (char));

                    /**fill-in search- and changeword**/
                    searchword[ct_replace][0] = '$';
                    strcat (searchword[ct_replace], argv[i + 1]);
                    changeword[ct_replace] = argv[i + 2];

                    /**count searchwords**/
                    ct_replace += 1;
        }


    if  ((fp = fopen (argv[1], "r")) == NULL) {
        fprintf (stderr, "\n%s could not open %s!\n", argv[0], argv[1]);
        exit (0);
    }

    fgets (zeile, 512, fp);

    while (! feof (fp)) {
        set_line = 0;

        /**get "set searchword changeword"; "set " has to be at the beginning of a line**/
        if  (strlen (zeile) > 1)
                if  (! strncmp ("set ", zeile, 4)) {

                    /**allocate memory**/
                    searchword[ct_replace] = (char *)malloc (512 * sizeof (char));
                    changeword[ct_replace] = (char *)malloc (512 * sizeof (char));

                    /**fill-in search- and changeword**/
                    searchword[ct_replace][0] = '$';
                    sscanf  (zeile, "set %s %s", searchword[ct_replace] + 1, changeword[ct_replace]);

                    /**test if searchword already exists (block the new then, keep the old, e.g. from the command line)**/
                    for (j = 0; j < ct_replace; ++j)
                        if  (! strcmp (searchword[j], searchword[ct_replace])) {
                            block_this = 1;
                            fprintf (stderr, "%s is set to %s intead of %s by file %s.\n",
                                             searchword[j], changeword[j], changeword[ct_replace], argv[1]);
                        }

                    /**count searchwords**/
                    if  (block_this == 0)
                        ct_replace += 1;
                    block_this = 0;

                    /**mark line with a "set"**/
                    set_line = 1;
                }


        sort_longest (searchword, changeword, ct_replace);

        if  ((strlen (zeile) > 1) && (set_line == 0))

            /**for each letter in the text**/
            for (i = 0; i < (signed int)(strlen (zeile)); ++i) {

                sel = -1;
                for (ct_ch = 0; ((ct_ch < ct_replace) && (sel == -1)); ++ct_ch)
                    if  (! strncmp (searchword[ct_ch], zeile + i, strlen (searchword[ct_ch])))
                        sel = ct_ch;                              /**state that and which word to replace**/

                if  (sel != -1) {
                    i += strlen (searchword[sel]) - 1;            /**skip the searchword (i counted anyway)**/
                    fprintf (stdout, "%s", changeword[sel]);      /**print the changeword**/
                } else {
                    fprintf (stdout, "%c", zeile[i]);             /**print all regular characters**/
                }
        }

        fgets (zeile, 512, fp);
    }

    return (0);
}
