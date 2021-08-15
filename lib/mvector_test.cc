/*
 * mvector_test.cc -- Test for math vectors
 *
 * This file is part of glsim, a numerical simulation class library and
 * helper programs.
 *
 * glsim is copyright (C) 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016,
 *                        2017, 2018, 2019, 2020, 2021
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

#include "blib/mvector.hh"

int main(int argc, char *argv[])
{
  std::cout << "***** Testing mvector_2D" << "\n";

  glsim::mvector_2D a2,b2;
  a2={2,2};
  b2={3,3};

  std::cout << "Loaded vectors, dimension " << a2.dim() << " " << a2 << " " << b2 << "\n";
  b2=a2;
  std::cout << "Assigned b2 " << b2 << "\n";

  glsim::mvector_2D c2;
  c2=a2+b2;

  std::cout << "a2+b2" << c2 << "\n";
  std::cout << "a2-b2" << a2-b2 << "\n";
  a2+={1,1};
  std::cout << "a2+=(1,1):  " << a2 << "\n";


  std::cout << "\n\n***** Testing mvector_3D" << "\n";

  glsim::mvector_3D a3,b3;
  a3={1,1,1};

  std::cout << "Loaded vectors, dimension " << a3.dim() << " " << a3 << " " << b3 << "\n";
  b3=a3;
  std::cout << "Assigned b3 " << b3 << "\n";
  b3={10,20,5};
  std::cout << "a3 + " << b3 << " =  " << a3+b3 << "\n";
  a3+=b3;
  std::cout << "a3 + b3 - b3 = " << a3-b3 << "\n";

  std::cout << "Normalize, modsq " << a3.normalize().modsq() << "\n";
  std::cout << "Scalar mult b3*(-1) = " << -1*b3 << "\n";
  std::cout << "a3 cross b3 " << cross(a3,b3) << "\n";
  std::cout << "b3 cross a3 " << cross(b3,a3) << "\n";
  std::cout << "crossp dotp " << a3.dot(cross(a3,b3)) << "\n";




  std::cout << "\n\n***** Testing mvector<10>" << "\n";

  glsim::mvector<10> a10,b10;
  a10={2,5,80};
  b10=a10;
  std::cout << "Loaded vectors, dimension " << a10.dim() << " " << a10 << " " << b10 << "\n";
  std::cout << "b10 + a10 = " << b10+a10 << "\n";


  return 0;
}
