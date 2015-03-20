/*
 * graph.hh -- classes for simulations on graphs/lattices
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

#ifndef _GRAPH_HH_
#define _GRAPH_HH_

#include <fstream>

#include <glsim/exception.hh>
#include <glsim/random.hh>

/*
 @T
 \section{Graph Base}

 [[GraphBase]] is not intended to be instantiated.  It defines the
 interface common to all graphs, and some internal data structures.

 we could add a consistency check function here, for the benefit of
 Graph's children, such as check all nodes have the right number of
 neighbours, etc

 These classes conceptually abstratct.  GraphBase is the base for all
 graph or lattice classes, stores nodes and optionally a neighbour for
 each node.  GraphBondBase add functionality to list and walk through
 bonds.

 Types: 
 \begin{itemize}
 \item nodes will have a unique id, of type [[ptrdiff_t]]
 \item nodeT (anything "simple")
 \item [[bondw_t]] (currently double)
 \end{itemize}

 GraphBase guarantees access through methods like vector.

 The topology is hidden here; neighbours must be made available by children.

 The implementation currently is simply a C-style array.  Conversion
 between pointers and ids is just a matter of substracting; but this
 is not something we commit to.  Children should convert through
 methods: pointer-->id with  [[id(nodeT*)]], id-->pointer and begin()+n

 For operations involving bulk operations which do not depend on
 topology, the fastest way is to acquire a pointer to the contiguous
 storage of nodes through the data() method (like C++11 vectors).

 The iterators provided are potentially expensive
 
 \begin{itemize}

 \item [[node_iterator]] iterates through nodes, but with knowledge of
 topology.  Capable of locating neighbours (w/specializations for
 fastest approach)

 \item [[neighbour_iterator]] iterates through all of a given node's
 neighbours.  Might be substantially slower than using the provided
 [[for_each...]] methods

 \end{itemize}

 @c */

typedef double bondw_t;

template <typename nodeT>
class bond {
public:
  nodeT   *n[2];
  bondw_t weight;

  bond(nodeT *node1,nodeT *node2,bondw_t w) :
    weight(w)
  {n[0]=node1; n[1]=node2;}
} ;

template <typename nodeT>
class GraphBase {
public:

  // Types and constants

  typedef      nodeT      node_t;
  static const ptrdiff_t  nilnode=-1;
  
  // Graph info

  int  size() const;
  int  coordination_number() const;
  int  neighbour_size(int i) const {return connectivity>0 ? connectivity : Nneighbours[i];}

  // Node access

  nodeT&       operator[](ptrdiff_t i);
  const nodeT& operator[](ptrdiff_t i) const;
  nodeT&       at(ptrdiff_t i);              // (w/bound check)
  const nodeT& at(ptrdiff_t i) const;
  ptrdiff_t    id(node_t *n) const;
  nodeT*       data() noexcept {return nodes;}
  const nodeT* data() const noexcept {return nodes;}

  // Input/output

  void read(std::istream&);
  void write(std::ostream&);

  // Neighbour access : only through iterators to avoid exposing the
  // possibly void neighbour list

  // Iterators

  class   node_iterator;
  class   neighbour_iterator;
  // class bond_iterator;

  node_iterator      begin();
  nodeT*             end();
  neighbour_iterator neighbour_begin(ptrdiff_t);
  neighbour_iterator neighbour_begin(nodeT* n) {return neigbour_begin(n-nodes);}
  nodeT**            neighbour_end(ptrdiff_t);
  nodeT**            neighbour_end(nodeT* n) {return neighbour_end(n-nodes);}

protected:
  GraphBase();
  GraphBase(int size, int connectivity, bool use_neighbour_list);
  ~GraphBase() {cleanup();}
  GraphBase(const GraphBase&);
  GraphBase& operator=(const GraphBase&);

  void add_bond(int n1,int n2);
  nodeT **node_neighbours(nodeT *n) {return neighbours[n-nodes];}
  nodeT **node_neighbours(ptrdiff_t i) {return neighbours[i];}

  int    Nnodes;
  nodeT* nodes;
  int    connectivity,*Nneighbours; // connectivity<0 means fluctuating
  nodeT  ***neighbours;

private:
  void cleanup();
  void copy_into_empty(const GraphBase&);
} ;

// @q

// @T
// \subsection{Construction and destruction}
// @c 

template <typename nodeT> inline
GraphBase<nodeT>::GraphBase() :
  nodes(0), Nnodes(0),
  neighbours(0), Nneighbours(0)
{}

template <typename nodeT>
GraphBase<nodeT>::GraphBase(int sizen,int connect,bool use_neighbour_list) : 
  Nnodes(sizen),
  connectivity(connect),
  neighbours(0), 
  Nneighbours(0)
{
  nodes=new nodeT[Nnodes];
  if (use_neighbour_list) {
    Nneighbours=new int[Nnodes+1];  // to allow a pointer 1-past-end
    for (int i=0; i<Nnodes; i++) Nneighbours[i]=0;
    neighbours=new nodeT**[Nnodes];
    for (int i=0; i<Nnodes; i++) neighbours[i]=0;
  }
}

// Copy construction and assignment

template <typename nodeT> inline
GraphBase<nodeT>::GraphBase(const GraphBase& g) :
  nodes(0), neighbours(0), Nneighbours(0)
{
  copy_into_empty(g);
}

template <typename nodeT> inline
GraphBase<nodeT>& GraphBase<nodeT>::operator=(const GraphBase& g) 
{
  if (this==&g) return *this;
  cleanup();
  copy_into_empty(g);
  return *this;
}

/* @t
 Private methods for copying and cleanup.
 @c */

template <typename nodeT>
void GraphBase<nodeT>::copy_into_empty(const GraphBase& g)
{
  Nnodes=g.Nnodes;
  nodes=new nodeT[Nnodes];
  memcpy(nodes,g.nodes,Nnodes*sizeof(nodeT));
  connectivity=g.connectivity;
  if (g.Nneighbours!=0) {
    Nneighbours=new int[Nnodes];
    memcpy(Nneighbours,g.Nneighbours,Nnodes*sizeof(int));
    neighbours=new nodeT**[Nnodes];
    for (int n=0; n<Nnodes; n++) {
      neighbours[n]=(nodeT **) malloc(sizeof(nodeT*)*Nneighbours[n]);
      for (int j=0; j<Nneighbours[n]; j++) {
	int nn=g.neighbours[n][j]-g.nodes;
	neighbours[n][j]=nodes+nn;
      }
    }
  } else Nneighbours=0;
}

template <typename node>
void GraphBase<node>::cleanup()
{
  if (Nneighbours) delete[] Nneighbours;
  if (neighbours) {
    for (int i=0; i<Nnodes; i++) 
      if (neighbours[i]) free(neighbours[i]);
    delete[] neighbours;
  }
  if (nodes) delete[] nodes;
}

/* @t
  This is never called from the class, it is to be used by children to add links.  GraphBondBase must override this to add the bond to its list of bonds.
  @c */

template <typename nodeT>
void GraphBase<nodeT>::add_bond(int n1,int n2)
{
  Nneighbours[n1]++;
  neighbours[n1]=(nodeT **) realloc(neighbours[n1],
				   sizeof(nodeT*)*Nneighbours[n1]);
  neighbours[n1][Nneighbours[n1]-1]=nodes+n2;

  Nneighbours[n2]++;
  neighbours[n2]=(nodeT **) realloc(neighbours[n2],
				   sizeof(nodeT*)*Nneighbours[n2]);
  neighbours[n2][Nneighbours[n2]-1]=nodes+n1;
}  

/* @t
   \subsection{Element access}
   @c */

template <typename nodeT>
inline int GraphBase<nodeT>::size() const
{
  return Nnodes;
}

template <typename nodeT> inline
int GraphBase<nodeT>::coordination_number() const
{
  return connectivity;
}

template <typename nodeT> inline
ptrdiff_t GraphBase<nodeT>::id(nodeT *n) const
{
  return n-nodes;
}

template <typename nodeT> inline
const nodeT& GraphBase<nodeT>::operator[](ptrdiff_t i) const
{
  return nodes[i];
}

template <typename nodeT> inline
nodeT& GraphBase<nodeT>::operator[](ptrdiff_t i)
{
  return nodes[i];
}

template <typename nodeT> inline
const nodeT& GraphBase<nodeT>::at(ptrdiff_t i) const
{
  if (i<0 || i>=Nnodes)
    throw glsim::Out_of_range();
  return nodes[i];
}

template <typename nodeT> inline
nodeT& GraphBase<nodeT>::at(ptrdiff_t i)
{
  if (i<0 || i>=Nnodes)
    throw glsim::Out_of_range();
  return nodes[i];
}

/* @t
   \subsection{In/out}
   @c */

template <typename nodeT>
void GraphBase<nodeT>::write(std::ostream &of)
{
  of.write((char*) &Nnodes,sizeof(Nnodes));
  of.write((char*) nodes,Nnodes*sizeof(nodeT));
  of.write((char*) &connectivity,sizeof(int));
  bool use_neighbour_list=Nneighbours!=0;
  of.write((char*) &use_neighbour_list,sizeof(bool));
  if (use_neighbour_list) {
    of.write((char*) Nneighbours,Nnodes*sizeof(int));
    for (int n=0; n<Nnodes; n++)
      for (int j=0; j<Nneighbours[n]; j++) {
	int nn=neighbours[n][j]-nodes;
	of.write((char*) &nn,sizeof(nn));
      }
  }
}

template <typename nodeT>
void GraphBase<nodeT>::read(std::istream &ifs)
{
  cleanup();
  ifs.read((char*) &Nnodes,sizeof(Nnodes));
  nodes=new nodeT[Nnodes];
  ifs.read((char*) nodes,Nnodes*sizeof(nodeT));
  ifs.read((char*) &connectivity,sizeof(int));
  bool use_neighbour_list;
  ifs.read((char*) &use_neighbour_list,sizeof(bool));
  if (use_neighbour_list) {
    Nneighbours=new int[Nnodes];
    ifs.read((char*) Nneighbours,Nnodes*sizeof(int));
    neighbours=new nodeT**[Nnodes];
    for (int n=0; n<Nnodes; n++) {
      neighbours[n]=(nodeT **) malloc(sizeof(nodeT*)*Nneighbours[n]);
      for (int j=0; j<Nneighbours[n]; j++) {
 	int nn;
 	ifs.read((char*) &nn,sizeof(nn));
 	neighbours[n][j]=nodes+nn;
      }
    }
  } else Nneighbours=0;
}

/* @t
   \subsection{Iterators}
   @c */

template <typename nodeT> inline
typename GraphBase<nodeT>::neighbour_iterator
GraphBase<nodeT>::neighbour_begin(ptrdiff_t i)
{
  return neighbour_iterator(neighbours[i]);
}

template <typename nodeT> inline
nodeT** GraphBase<nodeT>::neighbour_end(ptrdiff_t i)
{
  return neighbours[i]+Nneighbours[i];
}

template <typename nodeT>
class GraphBase<nodeT>::neighbour_iterator :
    public std::iterator<std::bidirectional_iterator_tag,nodeT>
{
public:
  neighbour_iterator(nodeT** first_n) :
    nn(first_n) {}

  bool operator==(const GraphBase<nodeT>::neighbour_iterator &i)
  {return nn==i.nn;}

  bool operator!=(const GraphBase<nodeT>::neighbour_iterator &i)
  {return nn!=i.nn;}

  bool operator==(nodeT** p)
  {return nn==p;}

  nodeT& operator*() const {return **nn;}
  nodeT* operator->() const {return *nn;}

  neighbour_iterator& operator++() {++nn; return *this;}
  neighbour_iterator& operator++(int) {neighbour_iterator c=*this; ++nn; return *this;}
  neighbour_iterator& operator--() {nn--; return *this;}
  neighbour_iterator& operator--(int) {neighbour_iterator c=*this; --nn; return *this;}

protected:
  nodeT **nn;
} ;

/* @t
   \subsection{Complex iterators}
   @c */

template <typename nodeT> inline
typename GraphBase<nodeT>::node_iterator
GraphBase<nodeT>::begin()
{
  return node_iterator(*this);
}

template <typename nodeT> inline
nodeT* GraphBase<nodeT>::end()
{
  return nodes+Nnodes;
}

template <typename nodeT>
class GraphBase<nodeT>::node_iterator :
  public std::iterator<std::bidirectional_iterator_tag,nodeT>
{
public:
  node_iterator(GraphBase<nodeT> &g,int ini_offset=0) :
    graph(g), nn(g.nodes+ini_offset) {}

  node_iterator(const GraphBase<nodeT>::node_iterator &i) :
    graph(i.graph), nn(i.nn) {}

  node_iterator& operator=(const GraphBase<nodeT>::node_iterator &i)
  {graph=i.graph; nn=i.nn; return *this;}

  bool operator==(const node_iterator &i)
  {return nn==i.nn;}

  bool operator!=(const node_iterator &i)
  {return nn!=i.nn;}

  bool operator==(nodeT* p)
  {return nn==p;}

  bool operator!=(nodeT* p)
  {return nn!=p;}

  nodeT& operator*() const {return *nn;}
  nodeT* operator->() const {return nn;}

  node_iterator& operator++() {++nn; return *this;}
  node_iterator& operator++(int) {node_iterator c=*this; ++nn; return c;}

  node_iterator& operator--() {nn--; return *this;}
  node_iterator& operator--(int) {node_iterator c=*this; --nn; return c;}

  // Above methods are required by a standard iterator; we add the following

  node_iterator& to(int n) {nn=graph.nodes+nn; return *this;}
  node_iterator& to_neighbour(int n) {nn=graph.node_neighbours(nn)+n; return *this;}
  nodeT* neighbour(int i) const {return graph.node_neighbours(nn)+i;}

protected:
  GraphBase<nodeT>& graph;
  nodeT *nn;
};

#endif /* _GRAPH_HH_ */
