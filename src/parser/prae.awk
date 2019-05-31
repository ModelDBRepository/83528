#!/bin/awk -f

BEGIN { i=0; }

function alength(a){
  i = 0; for(v in a) i++; return i;
}

# insertion sort found in traceroute
function sortthisarrays(ar1,ar2){
   n = alength(ar1);
   for (i = 2; i <= n; ++i) {
     v = ar1[i]; v2 = ar2[i]; j = i - 1;
     while (ar1[j] > v) {
             ar1[j+1] = ar1[j];
             ar2[j+1] = ar2[j];
             j = j - 1;
             if (j < 0)
                break;
     }
     ar1[j+1] = v;
     ar2[j+1] = v2
   }
}


/^set/ { n = split($0,erg);
         searchword[i] = erg[2];
         changeword[i] = erg[3];
         i++;
         next;         # Stop processing the current input record (line).
       }

  {
    sortthisarrays(searchword,changeword);      # sort such that longest words are replaced first (else, beginning parts of longer words could be replaced)
    for (u in searchword)
      if (length(searchword[u])>0)
        gsub(sprintf("\\\$%s", searchword[u]),
             sprintf("%s", changeword[u]),$0);
    print $0
  }

END { }
