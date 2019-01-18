/*
 * configuration.hh -- declaration of class Configuration
 * This file is part of glsim, a numerical simulation class library and
 * helper programs.
 *
 * glsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016,
 *                        2017, 2018, 2019
 * by Tomas S. Grigera.
 * 
 * glsim is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License (GPL) as
 * published by the Free Software Foundation. You can use either
 * version 3, or (at your option) any later version.
 * 
 * If you use glsim to produced published work, or if you redistribute a
 * modified version of glsim, or code based on glsim, please cite us
 *
 *  - T. S. Grigera, /glsim: A general library for numerical
 *    simulation/.  Computer Physics Communications *182*, 2122-2131
 *    (2011).
 *
 * glsim is developed as part of academic work.  As such author
 * attributions are not only a way to give recognition to the original
 * authors, but also help to ensure continued development, since
 * citations show that the work has been useful to others and thus
 * help convince funding agencies that it is worth funding.  * glsim
 * distribution.
 *
 * glsim is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 * or FITNESS FOR A PARTICULAR PURPOSE.
 * 
 * For details see the file LICENSE in the home directory. If the file is
 * missing, contact the maintainers.
 *
 */

#ifndef CONFIGURATION_HH
#define CONFIGURATION_HH

#include <string>

namespace glsim {

/** \class Configuration
    \ingroup Simulation
*/
class Configuration {
public:
  std::string name;

  Configuration(const std::string& name_="[no name]") : name(name_) {}
  virtual ~Configuration() {}

  void init(const std::string& fname) 
            {fname.empty() ? this->init() : this->init(fname.c_str());}
  void load(const std::string& fname) {this->load(fname.c_str());}
  void save(const std::string& fname) {this->save(fname.c_str());}

  virtual void init() {init(0);}
  virtual void init(const char* fname)=0;
  virtual void load(const char* fname)=0;
  virtual void save(const char* fname)=0;

} ;


} /* namespace */

#endif /* CONFIGURATION_HH */
