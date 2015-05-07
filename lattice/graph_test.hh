/*
 * graph_test.hh -- testing for classes defined in graph.hh
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

#ifndef GRAPH_TEST_HH
#define GRAPH_TEST_HH

template <typename GraphT>
class Graph_test {
public:
  virtual void test()=0;

  virtual void test_info(GraphT& graph,int size,int connectivity);
  virtual void test_plain_access(GraphT& graph);
  virtual void test_iterator(GraphT& graph);
} ;

template <typename GraphT>
void Graph_test<GraphT>::test_info(GraphT &graph,int size,int connectivity)
{
  check_result("size",graph.size(),size);
  check_result("coordination_number",graph.coordination_number(),connectivity);
}

template <typename GraphT>
void Graph_test<GraphT>::test_plain_access(GraphT &graph)
{
  double *d=graph.data();
  int N=graph.size();
  for (int i=0; i<N; ++i)
    *(d+i)=22;

  check_result("plain access",graph[N-3],22.0);
  check_result("plain access w/bound check",graph.at(N-3),22.0);

  std::cout << "Testing plain access w/boudn check 2: ";
  bool ok=false;
  try {
    check_result("plain access w/bound check",graph.at(N+3),22.0);
  } catch (glsim::Out_of_range& e) {
    ok=true;
  }
  if (ok)
    std::cout << "OK\n";
  else
    throw Test_failure("out of range exception","nothing");

  check_result("plain access: node id from pointer",graph.id(d+2),(ptrdiff_t) 2);
}

//
// Test bidirectional iterator operations
//
template <typename GraphT>
void Graph_test<GraphT>::test_iterator(GraphT &graph)
{
  typename GraphT::node_iterator ni1(graph);
  typename GraphT::node_iterator ni2(graph);

  ni1=graph.begin();
  ni2=typename GraphT::node_iterator(graph,graph.data()+1);

  check_result("node_iterator::operator==",ni1==ni2,false);
  ni2=ni1;
  check_result("node_iterator::operator==",ni1==ni2,true);
  check_result("node_iterator::operator!=",ni1!=ni2,false);
  ++ni1;
  check_result("node_iterator::operator!=",ni1!=ni2,true);
  check_result("node_iterator::operator!=",ni1!=graph.data(),true);
  check_result("node_iterator::operator==",ni2==graph.data(),true);

  for (int i=0; i<graph.size(); ++i)
    *(graph.data()+i)=i;

  ni1.to(graph.data());
  ni2.to(graph.data()+2);
  check_result("node iterator move and access",*ni1,0.);
  check_result("node iterator move and access",*ni2,2.);
  --ni2;
  check_result("node iterator move and access",*ni2,1.);
  double *d=ni1;
  check_result("node iterator cast to pointer",*d,*ni1);
}

#endif /* GRAPH_TEST_HH */
