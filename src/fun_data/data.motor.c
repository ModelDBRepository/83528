#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <math.h> /**for data_lang_assoc**/

#include "iter.h"
#include "series.h"
#include "data.h"
#include "utils.h"


/**************************** import_motor ***********************************/
/* Assigned under global: data_import and reads files data_files.            */
/* Marks file format, as given by PDP. E.g. norm2_train_12.pat.              */

int import_motor (PARAMS *x, DATA *d) {

    const int num_sens = 12;
    const int num_time = 10;
    const int num_item = 14;
    const int max_data = 10000;
    const int num_zeros = 144;

    if  ((d->Bilder = (double **)malloc(max_data*sizeof(double *))) == NULL)
        fprintf (stderr, "\nout of memory\n");

    d->anzahl = 0;

    for (int file_nr = 0; file_nr < x->anz_data_files; ++file_nr) {

        fprintf (stderr, "\nimport_motor: reading %s", x->data_files[file_nr]);

        FILE *fp = fopen (x->data_files[file_nr], "r");

        while (!feof(fp)) {

            d->Bilder[d->anzahl] = d_vector (num_sens * num_time + num_item);

            char tag[1000];
            fscanf (fp, "%s", tag);
            fprintf (stderr, "\nfile: %s", x->data_files[file_nr]);
            fprintf (stderr, "  tag: %s", tag);
            fprintf (stderr, "  num: %d\n", d->anzahl);
            if  (strncmp(tag, "event", 5) != 0)
                fprintf (stderr,"\nerror in line %d", d->anzahl);

            {
                int OK = 0;

                char *item_word = strstr (tag + 6, "_") + 1;
                if  ((item_word > tag + 9) || (item_word <= tag + 7))
                    item_word = tag + 8;
                fprintf (stderr, "item_word: %s\n", item_word);

                for (int i = 0; i < num_item; ++i)
                    d->Bilder[d->anzahl][num_sens * num_time + i] = 0.0;

                if  (! strcmp (item_word, "headdown")) {
                    d->Bilder[d->anzahl][num_sens * num_time + 0] = 1.0;
                    OK = 1;
                }
                if  (! strcmp (item_word, "headleft") ||
                     ! strcmp (item_word, "eventleft") ||
                     ! strcmp (item_word, "headdleft")) {
                    d->Bilder[d->anzahl][num_sens * num_time + 1] = 1.0;
                    OK = 1;
                }
                if  (! strcmp (item_word, "headright") ||
                     ! strcmp (item_word, "heeadright")) {
                    d->Bilder[d->anzahl][num_sens * num_time + 2] = 1.0;
                    OK = 1;
                }
                if  (! strcmp (item_word, "headup") ||
                     ! strcmp (item_word, "headu")) {
                    d->Bilder[d->anzahl][num_sens * num_time + 3] = 1.0;
                    OK = 1;
                }
                if  (! strcmp (item_word, "put")) {
                    d->Bilder[d->anzahl][num_sens * num_time + 4] = 1.0;
                    OK = 1;
                }
                if  (! strcmp (item_word, "drop")) {
                    d->Bilder[d->anzahl][num_sens * num_time + 5] = 1.0;
                    OK = 1;
                }
                if  (! strcmp (item_word, "pick")) {
                    d->Bilder[d->anzahl][num_sens * num_time + 6] = 1.0;
                    OK = 1;
                }
                if  (! strcmp (item_word, "lift")) {
                    d->Bilder[d->anzahl][num_sens * num_time + 7] = 1.0;
                    OK = 1;
                }
                if  (! strcmp (item_word, "touch")) {
                    d->Bilder[d->anzahl][num_sens * num_time + 8] = 1.0;
                    OK = 1;
                }
                if  (! strcmp (item_word, "forward")) {
                    d->Bilder[d->anzahl][num_sens * num_time + 9] = 1.0;
                    OK = 1;
                }
                if  (! strcmp (item_word, "goto")) {
                    d->Bilder[d->anzahl][num_sens * num_time + 10] = 1.0;
                    OK = 1;
                }
                if  (! strcmp (item_word, "moveback")) {
                    d->Bilder[d->anzahl][num_sens * num_time + 11] = 1.0;
                    OK = 1;
                }
                if  (! strcmp (item_word, "turnleft") ||   /**because:**/
                     ! strcmp (item_word, "_turnleft")) {  /**"evewnt"**/
                    d->Bilder[d->anzahl][num_sens * num_time + 12] = 1.0;
                    OK = 1;
                }
                if  (! strcmp (item_word, "turnright")) {
                    d->Bilder[d->anzahl][num_sens * num_time + 13] = 1.0;
                    OK = 1;
                }
                if  (!OK) {
                    fprintf (stderr, "\nitem_word=%s\n", item_word);
                    exit (1);
                }
           }

            for (int i = 0; i < num_sens * num_time; ++i)
                fscanf (fp, "%lf", d->Bilder[d->anzahl] + i);

            for (int i = 0; i < num_sens * num_time; ++i)
                fprintf (stderr, "%f", d->Bilder[d->anzahl][i]);
            fprintf (stderr, "\n");

            for (int i = 0; i < num_zeros; ++i) {
                double zero = 1.0;
                fscanf (fp, "%lf", &zero);
		if  (zero != 0.0)
                    fprintf (stderr, "\nnon-zero in line %d", d->anzahl);
            }

            d->anzahl += 1;
        }

	fclose (fp);
        d->anzahl -= 1;
    }

    fprintf (stderr, "  import_motor end with %d data read\n", d->anzahl);

    return (0);
}


/**************************** import_lang_assoc ******************************/
/* Assigned under global: data_import and reads files data_files.            */
/* Like import_motor, but for Mark's data file: lang_assoc.txt.              */

int import_lang_assoc (PARAMS *x, DATA *d) {

    const int num_sens = 13;
    const int num_time = 10;
    const int num_item = 14;
    const int max_data = 10000;

    if  ((d->Bilder = (double **)malloc(max_data*sizeof(double *))) == NULL)
        fprintf (stderr, "\nout of memory\n");

    d->anzahl = 0;

    for (int file_nr = 0; file_nr < x->anz_data_files; ++file_nr) {

        fprintf (stderr, "\nimport_motor: reading %s", x->data_files[file_nr]);

        FILE *fp = fopen (x->data_files[file_nr], "r");

	/**read the first line to skip for data**/
        fprintf (stderr, "\nimport_motor: first line:\n");
        char ch = 'x';
        while (ch != '\n') {
              fscanf (fp, "%c", &ch);
              fprintf (stderr, "%c", ch);
        }


        while (!feof(fp)) {

            d->Bilder[d->anzahl] = d_vector (num_sens * num_time + num_item);

            char tag[1000];
            fscanf (fp, "%s", tag);
            fprintf (stderr, "\nfile: %s", x->data_files[file_nr]);
            fprintf (stderr, "  tag: %s", tag);
            fprintf (stderr, "  num: %d\n", d->anzahl);
            if  (strncmp(tag, "event", 5) != 0)
                fprintf (stderr,"\nerror in line %d", d->anzahl);

            {
                int OK = 0;

                char *item_word = strstr (tag + 6, "_") + 1;
                if  ((item_word > tag + 9) || (item_word <= tag + 7))
                    item_word = tag + 8;
                fprintf (stderr, "item_word: %s\n", item_word);

                for (int i = 0; i < num_item; ++i)
                    d->Bilder[d->anzahl][num_sens * num_time + i] = 0.0;

                if  (! strcmp (item_word, "headdown")) {
                    d->Bilder[d->anzahl][num_sens * num_time + 0] = 1.0;
                    OK = 1;
                }
                if  (! strcmp (item_word, "headleft") ||
                     ! strcmp (item_word, "eventleft") ||
                     ! strcmp (item_word, "headdleft")) {
                    d->Bilder[d->anzahl][num_sens * num_time + 1] = 1.0;
                    OK = 1;
                }
                if  (! strcmp (item_word, "headright") ||
                     ! strcmp (item_word, "heeadright")) {
                    d->Bilder[d->anzahl][num_sens * num_time + 2] = 1.0;
                    OK = 1;
                }
                if  (! strcmp (item_word, "headup") ||
                     ! strcmp (item_word, "headu")) {
                    d->Bilder[d->anzahl][num_sens * num_time + 3] = 1.0;
                    OK = 1;
                }
                if  (! strcmp (item_word, "put")) {
                    d->Bilder[d->anzahl][num_sens * num_time + 4] = 1.0;
                    OK = 1;
                }
                if  (! strcmp (item_word, "drop")) {
                    d->Bilder[d->anzahl][num_sens * num_time + 5] = 1.0;
                    OK = 1;
                }
                if  (! strcmp (item_word, "pick")) {
                    d->Bilder[d->anzahl][num_sens * num_time + 6] = 1.0;
                    OK = 1;
                }
                if  (! strcmp (item_word, "lift")) {
                    d->Bilder[d->anzahl][num_sens * num_time + 7] = 1.0;
                    OK = 1;
                }
                if  (! strcmp (item_word, "touch")) {
                    d->Bilder[d->anzahl][num_sens * num_time + 8] = 1.0;
                    OK = 1;
                }
                if  (! strcmp (item_word, "forward")) {
                    d->Bilder[d->anzahl][num_sens * num_time + 9] = 1.0;
                    OK = 1;
                }
                if  (! strcmp (item_word, "goto")) {
                    d->Bilder[d->anzahl][num_sens * num_time + 10] = 1.0;
                    OK = 1;
                }
                if  (! strcmp (item_word, "moveback")) {
                    d->Bilder[d->anzahl][num_sens * num_time + 11] = 1.0;
                    OK = 1;
                }
                if  (! strcmp (item_word, "turnleft") ||   /**because:**/
                     ! strcmp (item_word, "_turnleft")) {  /**"evewnt"**/
                    d->Bilder[d->anzahl][num_sens * num_time + 12] = 1.0;
                    OK = 1;
                }
                if  (! strcmp (item_word, "turnright")) {
                    d->Bilder[d->anzahl][num_sens * num_time + 13] = 1.0;
                    OK = 1;
                }
                if  (!OK) {
                    fprintf (stderr, "\nitem_word=%s\n", item_word);
                    exit (1);
                }
           }

            for (int i = 0; i < num_sens * num_time; ++i)
                fscanf (fp, "%lf", d->Bilder[d->anzahl] + i);

            for (int i = 0; i < num_sens * num_time; ++i)
                fprintf (stderr, "%f", d->Bilder[d->anzahl][i]);
            fprintf (stderr, "\n");

	    /**no discarding of zero's any more**/

            d->anzahl += 1;
        }

	fclose (fp);
        d->anzahl -= 1;
    }

    fprintf (stderr, "  import_lang_assoc end with %d data read\n", d->anzahl);

    return (0);
}


/******************************** data_motor *********************************/
/* Selects a random data item stored in d->Bilder.                           */

void data_motor (AGENT *z, COMMAND *cmd, DATA *d, int begin, int end) {

    const int num_sens = 12;
    const int num_time = 10;
    const int num_item = 14;

    static int select = 0;

    if  (cmd->quantum[0][0] == 1)
        select = (int)(drand48() * d->anzahl);

    if  ((z->d_b == num_sens) || (z->d_a == num_time)) {
        for (int ct_t = begin; ct_t < end; ++ct_t)
            for (int j = 0; j < z->d_a * z->d_b; ++j)
                cmd->S_target[ct_t][j] = d->Bilder[select][j];
        return;
    }

    if  (z->d_a * z->d_b == num_item) {
        for (int ct_t = begin; ct_t < end; ++ct_t)
            for (int j = 0; j < z->d_a * z->d_b; ++j)
                cmd->S_target[ct_t][j]
                                  = d->Bilder[select][num_sens * num_time + j];
        return;
    }


    fprintf (stderr, "\ndata_motor: dimensions must be ");
    fprintf (stderr, " d_a=%d, d_b=%d ", num_time, num_sens);
    fprintf (stderr, " or  d_a*d_b=%d!\n", num_item);
    exit(0);
}



/******************************** data_lang_assoc ****************************/
/* Selects a random data item stored in d->Bilder.                           */
/* q[0][0] = 1: select a new data point                                      */
/*         < 0: select (positive) chosen example (only useful for Gaussian)  */
/* q[1][0] = 1: draw a Gaussian which moves along d_b in time                */
/*         = 2: set target to data values which also change in time          */
/* q[2][0] = sigma of Gaussian                                               */
/* q[3][0] = height    "                                                     */
/* q[4][0] = velocity  "                                                     */

void data_lang_assoc (AGENT *z, COMMAND *cmd, DATA *d, int begin, int end) {

    const int num_sens = 13;
    const int num_time = 10;
    // const int num_item = 14;

    static int select = 0;

    if  (cmd->anz_quantums != 5)
      fprintf (stderr, "\ndata_lang_assoc wants 5 parameters!\n\n");

    if  (cmd->quantum[0][0] == 1)
        select = (int)(drand48() * d->anzahl);


    /** !!! !!! !!! !!! !!! !!! used for 04.motor.simplified.enr!!! !!! !!! !!! !!! !!! !!! !!! !!! !!! *

    if  (cmd->anz_quant[2] == 2)
        if  (cmd->quantum[2][1] == 4)
            if  (cmd->quantum[0][0] == 1)
                select = (int)(drand48() * 4 * 20);

    * !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! **/



    if  (cmd->quantum[1][0] == 1) {

        double sigma  = cmd->quantum[2][0];
        double height = cmd->quantum[3][0];
        double vel    = cmd->quantum[4][0];
        double cm     = 0.0;

        for (int ct_t = begin; ct_t < end; ++ct_t) {

            cm += vel;                                      /**move**/

            for (int i = 0; i < z->d_a; ++i)
            for (int j = 0; j < z->d_b; ++j) {

                double diffB = j - cm;

                if  (d->Bilder[select][num_sens * num_time + i] == 1.0)
                    cmd->S_target[ct_t][i * z->d_b + j] = height * exp (-0.5 * (diffB*diffB)/(sigma*sigma));
                else
                    cmd->S_target[ct_t][i * z->d_b + j] = 0.0;


                /** !!! !!! !!! !!! !!! !!! used for 04.motor.simplified.enr!!! !!! !!! !!! !!! !!! !!! !!! !!! !!! *

                if  (cmd->anz_quant[2] == 2)
                    if  (cmd->quantum[2][1] == 9)
                        cmd->S_target[ct_t][i * z->d_b + j] = height * exp (-0.5 * (diffB*diffB)/(sigma*sigma));

                  * !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! !!! **/


                /**in this case, determine the selected Gaussian**/
                if  (cmd->quantum[0][0] < 0) {
                    if  (i == (int)(-cmd->quantum[0][0]))
                        cmd->S_target[ct_t][i * z->d_b + j] = height * exp (-0.5 * (diffB*diffB)/(sigma*sigma));
                    else
                        cmd->S_target[ct_t][i * z->d_b + j] = 0.0;
                }
            }
        }
        return;
    }

    if  (cmd->quantum[1][0] == 2) {

        if  (z->d_a * z->d_b > num_sens)
            fprintf (stderr, "\ndata_lang_assoc wants small (1dim) area for max %d sensors", num_sens);

        for (int ct_t = 0; ct_t < end - begin; ++ct_t) /**during interval (always start with data from time 0)**/
            for (int j = 0; j < z->d_a * z->d_b; ++j)
                cmd->S_target[ct_t][j] = d->Bilder[select][j + num_sens * ct_t];

        return;
    }

    fprintf (stderr, "\ndata_lang_assoc: should not reach here ");
    exit(0);
}



/******************************** data_read_understood ***********************/
/* Reads an int via CORBA read_understood; writes it to binary target vector.*/

#if USE_CORBA

#include "HandleitC.h"
#include "HandleitI.h"
#include <CosNamingC.h>

void data_read_understood (AGENT *z, COMMAND *cmd, DATA *d, int begin, int end) {
    const int num_item = 14;
    CORBA::Long understood_value = -1;

   int argc = 1;
   char *dummy = "dummy";
   char **argv = &dummy;
   static Handleit_var f;
   static int firsttime = 1;

   fprintf (stderr, "\ndata_read_understood: ");

   if  (firsttime) {
       // initialize the ORB
       CORBA::ORB_var orb = CORBA::ORB_init (argc, argv);

       try {
           // resolve the naming service
            CORBA::Object_var nsobj =
                 orb->resolve_initial_references ("NameService");

           if (CORBA::is_nil(nsobj.in())) {
               cerr << "can't resolve NameService\n";
               exit(1);
           }

           // narrow the root naming context
           CosNaming::NamingContext_var nc
                       = CosNaming::NamingContext::_narrow(nsobj.in());

           // create a name component
           CosNaming::Name name;
           name.length (1);
           name[0].id = CORBA::string_dup ("Handleit");
           name[0].kind = CORBA::string_dup ("");

           // resolve the name component with the naming service
           CORBA::Object_var obj = nc->resolve(name);

           // downcast this object to Squareit
           f = Handleit::_narrow (obj.in ());

       }
       catch(...)
       {
         cerr << "\ndata_read_understood: CORBA exception 1 raised!" << endl;
       }

       firsttime = 0;
   }

   fprintf (stderr, " read ");

   try {
       // the call to the handler
       f->read_understood (understood_value);
       // orb->destroy(); // (not done)
   }
   catch(...)
   {
     cerr << "\ndata_read_understood: CORBA exception 2 raised!" << endl;
   }


   fprintf (stderr, "value is %d ", understood_value);

    if  (z->d_a /**z->d_b**/ != num_item)                                                    /**must match: first column of the area**/
        fprintf (stderr, "\ndata_read_understood: dimensions don't fit!\n");

    if  ((understood_value < 0) || (understood_value >= num_item))
        fprintf (stderr, "\ndata_read_understood: value %d out of range!\n",
                          understood_value);

    /*write to target*/
    for (int ct_t = begin; ct_t < end; ++ct_t)
        for (int i = 0; i < z->d_a; ++i) {                                                    /**counts only along columns**/
            for (int j = 0; j < z->d_b; ++j) {                                                    /**counts only along columns**/
                double diffB = j - 1;
                double height = 1.2;
                double sigma = 1.0;
                if  (i == understood_value)
                    cmd->S_target[ct_t][i * z->d_b + j] = height * exp (-0.5 * (diffB*diffB)/(sigma*sigma));
                else
                    cmd->S_target[ct_t][i * z->d_b + j] = 0.0;
            }
        }
}



void data_write_string (AGENT *z, COMMAND *cmd, DATA *d, int begin, int end) {

   int argc = 1;
   char *dummy = "dummy";
   char **argv = &dummy;
   static Handleit_var f;
   static int firsttime = 1;

   if  (firsttime) {
       // initialize the ORB
       CORBA::ORB_var orb = CORBA::ORB_init (argc, argv);

       try {
           // resolve the naming service
            CORBA::Object_var nsobj =
                 orb->resolve_initial_references ("NameService");

           if (CORBA::is_nil(nsobj.in())) {
               cerr << "can't resolve NameService\n";
               exit(1);
           }

           // narrow the root naming context
           CosNaming::NamingContext_var nc
                       = CosNaming::NamingContext::_narrow(nsobj.in());

           // create a name component
           CosNaming::Name name;
           name.length (1);
           name[0].id = CORBA::string_dup ("Handleit");
           name[0].kind = CORBA::string_dup ("");

           // resolve the name component with the naming service
           CORBA::Object_var obj = nc->resolve(name);

           // downcast this object to Squareit
           f = Handleit::_narrow (obj.in ());

       }
       catch(...)
       {
         cerr << "\ndata_write_string: CORBA exception 1 raised!" << endl;
       }

       firsttime = 0;
   }

   try {

       double max_act = 0.0;
       int X_winner = -1;
       int Y_winner = -1;
       for (int X = 0; X < z->d_a; X++)
           for (int Y = 0; Y < z->d_b; Y++)
               if  (cmd->S_from1[0][begin][Y + z->d_b * X] > max_act) {
                   max_act = cmd->S_from1[0][begin][Y + z->d_b * X];
                   X_winner = X;
                   Y_winner = Y;
               }

       char sentence[1024];
       sprintf (sentence, "Winning unit has coordinates (a=%d b=%d)", X_winner, Y_winner);

       CORBA::Long sentence_length = strlen (sentence);
       unsigned char * buffer = new unsigned char[sentence_length];

       for (int i = 0; i < sentence_length; ++i)
           buffer[i] = sentence[i];

       OctetSequenceIDL * togive = new OctetSequenceIDL(sentence_length, sentence_length, buffer, 1);

       // the call to the handler
       f->write_sonar_location (sentence_length, * togive);
       // orb->destroy(); // (not done)
   }
   catch(...)
   {
     cerr << "\ndata_write_string: CORBA exception 2 raised!" << endl;
   }
}

#endif /**#if USE_CORBA**/
