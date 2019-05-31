#ifndef NNsimObject_h
#define NNsimObject_h

#include <stdlib.h> 
#include <iostream>
#include <string>

#include "iter.h"
#include "series.h"
#include "data.h"
#include "parser/r.yacc.h"
#include "relax.h"               /**for do_series**/
#include "vehicle.h"             /**for printparams**/
#include "observe.h"             /**for export_weights**/
#include "utils.h"               /**for get_tmp_uname_filename**/

class NNsimObject
{
public:

  // constructor reads the language file and allocates memory  (used firsttime in Miro's open())
  NNsimObject (const int _analyze,
               const std::string & _tmp_uname,
               const std::string & _language_file,
               const std::string & _name);

  // if there are 2 series and the first has ilen=elen=1, then do it  (used firsttime in Miro's open())
  void FirstSeries ();

  // ilen=repeat within one action(); elen=periodicity of weight export; never ends  (used in Miro's action())
  int MainSeries ();

  // if there are more than 2 series, then do third (used in ...)
  void ThirdSeries ();

private:

  PARAMS     *x;
  AGENT      *Z;
  SIMULATION *si;
  AREA       *A;
  DATA       *d;

  int analyze;         // a behaviourParameter
  char *tmp_uname;     // a behaviourParameter
  char *language_file; // a behaviourParameter
  char *name;          // the behaviour name
};

#endif
