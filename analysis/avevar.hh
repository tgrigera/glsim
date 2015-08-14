/*
 * avevar.hh -- Class for average and variance with West recurrence
 *
 * This file is part of glsim, a numerical simulation class library and
 * helper programs.
 *
 * glsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015
 * by Tomas S. Grigera.
 * 
 * glsim is free software; you can redistribute it and/or modify it under
 * the terms of the GNU General Public License (GPL) as published by the
 * Free Software Foundation, with the additional requirements of
 * attribution and nonmisrepresentation. You can use either version 3, or
 * (at your option) any later version.
 * 
 * Additional terms under GNU GPL version 3 section 7:
 * 
 * When you redistribute this software, you are required to preserve its
 * author attributions. If you distribute verbatim copies, you must not
 * alter the AUTHORS file or attributions inserted in the source files,
 * and you must not change the software's name. If you distribute a
 * modified copy, then you must give clear notice that your work is
 * different from but based on glsim. You must distribute it under a
 * different name, but include a prominent notice specifying that "(your
 * package) is based on glsim version x.x", and provide a pointer to the
 * glsim distribution.
 *
 * glsim is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.
 * 
 * For details see the file LICENSE in the home directory. If the file is
 * missing, contact the maintainers.
 *
 */

#ifndef AVEVAR_HH
#define AVEVAR_HH

#include <limits>

namespace glsim {

/** \class AveVar
    \ingroup Analysis
    \brief Compute average and variance using West recurrence
*/
template <bool maxmin=true>
class AveVar {
public:
  AveVar() {clear();};
  AveVar& push(double);
  AveVar& clear();
  double  ave() const {return ave_;}
  double  var() const {return var_/(N_-1);}
  long    N() const {return N_;}
  double  max() const {return max_;}
  double  min() const {return min_;}
  
private:
  double ave_,var_;
  long   N_;
  double max_,min_;

} ;

template <bool maxmin>
inline AveVar<maxmin>& AveVar<maxmin>::clear()
{
  ave_=0;
  var_=0;
  N_=0;
  max_=std::numeric_limits<double>::min();
  min_=std::numeric_limits<double>::max();
}
  

template <bool maxmin>
inline AveVar<maxmin>& AveVar<maxmin>::push(double x)
{
  double Q,R;

  N_++;
  Q=x-ave_;
  R=Q/N_;
  ave_+=R;
  var_+=Q*R*(N_-1);

  if (maxmin) {
    if (x<min_) min_=x;
    if (x>max_) max_=x;
  }

  return *this;
}


} /* namespace */

#endif /* AVEVAR_HH */
