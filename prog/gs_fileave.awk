#!/usr/bin/gawk -f
#
# gs_fileave.awk --- average over file with same structure
#
# This script takes several files, groups all rows having an identical
# first column and averages the rest of the columns across files.
#

BEGIN {ncols=0;}

/^#/ {print $0; next;}

{
  if (NF>ncols) ncols=NF;  
  for (f=2; f<=NF; f++) {
      key[$1]=1;
      count[$1,f]++;
      Q=$f-ave[$1,f];
      ave[$1,f]+=Q/count[$1,f];
  }
}

END {
  for (k in key) {
    printf "%s",k;
    for (f=2; f<=ncols; f++)
      if ((k,f) in count) printf " %g",ave[k,f];
      else printf " --";
    printf "\n";
  }
}
