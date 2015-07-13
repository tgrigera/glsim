#!/usr/bin/gawk -f
#
# gs_columnave.awk: simple averaging of column-formatted files. Lines starting
#                   with # are ignored.
#

BEGIN {
  if (ARGV[ARGC-1]=="-sd") {
    sd=1;
    ARGC--;
  }
  N=0;
  print "# gs_columnave.awk: averages using the following files:";
  for (i=1; i<ARGC; i++)
    printf "# %s\n",ARGV[i];
}

/^#/ {print $0; next;}

{
  for (i=1; i<=NF; i++) {
    sum[i]+=$i;
    sumsq[i]+=$i*$i;
  }
  N++;
}

END {
  if (sd) printf "ave: ";
  for (i=1; i<=NF; i++) printf "%g ",sum[i]/N;
  print "";
  if (sd) {
    printf "sd: "
      for (i=1; i<=NF; i++) printf "%g ", sqrt(sumsq[i]/N-sum[i]*sum[i]/(N*N));
    print "";
  }
    
}
