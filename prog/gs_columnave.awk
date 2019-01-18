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
	Q=$i-ave[i];
	R=Q/(N+1);
	ave[i]+=R;
	var[i]+=Q*R*N;
    }
    if (NF>0) N++;
}


END {
  print "#\n# Nsamples =",N;
  if (sd) printf "ave: ";
  for (i=1; i<=NF; i++) printf "%.10g ",ave[i];
  print "";
  if (sd) {
      printf "sd: "
      for (i=1; i<=NF; i++) printf "%.10g ", sqrt(var[i]/(N-1));
      print "";
  }
}
