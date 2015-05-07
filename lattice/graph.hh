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

namespace glsim {

/**\defgroup lattice Graph and lattice
 *
 * \brief Classes for simulations on graphs and lattices
 * 
 * This is a collection of classes for lattice or graph based
 * simulations.  They aim to provide a common, STL-container-like
 * interface, that allows to formulate the algorithms in a generic
 * way, independent of the graph topology.
 *
 * Graphs are made of course of _nodes_ and _bonds_.  The data type of
 * the nodes is a template parameter.  Any type with a default
 * constructor is in principle possilbe, but the classes have been
 * designed under the assumption that nodes will be relatively simple
 * objects.  Bonds are not necessarily be represented explicitly,
 * since in many regular graphs it is more efficient to enumerate or
 * walk through them by following some simple rule.
 *
 * The base class for all graphs is GraphBase. This provides the
 * common interface, but it should not be insantiated.  The classes
 * implement specific kinds of graph (like :BetheLattice) add
 * functionality as necessary.  Acces to graph elements is best done
 * through GraphBase for bulk operations on nodes (GraphBase::data()
 * gives a pointer to the nodes, or you can use subscripting with
 * `GraphBase::operator[]()` or GraphBase::at()).  For operations that
 * need knowledge of the topology (neighbour structure), each graph
 * class provides an appropriate iterator (a bidirectional iterator).
 * Iterators have acommon interface (those operations required by a
 * bidirectional iterator) plus graph-dependent methdos to access
 * neighbours.  However in generic programming these are best avoided
 * in favor of algorithm functions such as for_each_neighbour().  Such
 * algorithm functions are defined generically for all GraphBase
 * descendants but optimized specializations are provided many
 * particular graphs, e.g. to allow loop unrolling.
 *
 * Finally, note that at present [April 2015] the classes are designed
 * to be used as templates, and avoid `virtual` functions: thus to be
 * used generically, use the template instantiation mechanism rather
 * than base class pointers.
 *
 */

/// \ingroup lattice
typedef ptrdiff_t id_t;  ///< \brief Type of node ids


/** \class GraphBase 
 * \ingroup lattice
 * \brief Base with common interface for all graph topologies
 *
 * This class defines the interface common to all graphs, and some
 * internal data structures, but it is not intended to be instantiated
 * (its constructor is protected).
 *
 * The template parameter is the data type of the nodes.  Any type
 * with a default constructor is in principle possilbe, but the
 * classes have been designed under the assumption that nodes will be
 * relatively simple objects.  Bonds are not necessarily represented
 * explicitly, as in many regular graphs it is more efficient to
 * enumerate them by following some rule rather than to store a list
 * of them.  Accordingly, GraphBase will hold the nodes (in a simple
 * array), but not a list of bonds (which we leave to `GraphBondBase`.
 * GraphBase can store (at the request of the descendant class), for
 * each node, a list of its neighbours (nodes joined to it by a bond).
 * A list of this kind is probably inefficient for many regular graphs
 * (such as Bravais lattices), and for this reason its construction is
 * left as an option.
 *    
 * GraphBase provides an interface to access nodes, acting as far as
 * possible as a standard container.  Nodes will be assigned a unique
 * id, of type `ptrdiff_t` (an integer type, `typedef`'d as `id_t`),
 * and GraphBase provides access to nodes by id through
 * `operator[](id_t)` or the method `at(id_t)` (which does
 * bound checking).  Of course the whole point of a graph class is to
 * give easy access to the topology of the graph: its neighbour
 * structure.  This is done through the graph's iterator (the most
 * general version is described below).  The graph iterator is
 * necessarily more complex than a pointer, so for bulk operations on
 * the nodes it is better to use `operator[](id_t id)` or to acquire a
 * pointer to the node storage (guaranteed contiguous) through the
 * `data()` method.
 *
 * The method `id(nodeT*)` returns the id of a node identified by a
 * pointer.  This is useful, for example, to obtain an iterator from a
 * simple pointer.
 *
 * The constructors are protected, so they can only be called from
 * derived classes.  The default constructor creates an empty graph;
 * this is only useful for later reading the graph from file,
 * otherwise the `init_storage()` method must be called explicitly at
 * some point, before the graph is used.
 *
 * # Specialization
 *
 * This class is not intended to be instantiated by the user.  It is
 * to be used to inherit from to define specific graphs (where a list
 * of bonds is not needed).  To inherit one must provide
 *
 * - New constructors with appropriate parameters, which should at
 *   some point call GraphBase(N,z,use_neighbour_list) or
 *   init_storage().
 *
 * - If the neighbour list is to be used, it must be constructed by
 *   calling add_bond().
 *
 * - If new data members are added, override read() and write() (but
 *   call the GraphBase versions to ensure the base data members are
 *   read or written).
 *
 * - If you don't use the neighbour list (or if you can provide a more
 *   efficient version), overrite begin() and end() to return a more
 *   specialized iterator.
 *
 * - If you can, write a specialized for_each_neighbour() function
 *   that relies on the particular topology of your graph to provide a
 *   more efficient implementation.
 *
 */
template <typename nodeT> class GraphBase {
public:

  /// \name Types and constants
  /// @{

  typedef      nodeT  node_t;       ///< Give acces to the type of nodes
  static const id_t   nilnode=-1;   ///< An invalid node id (to be used as null pointer)

  ///@}

  /// \name Graph information
  ///@{
  int  size() const;                 ///< Number of nodes
  int  coordination_number() const;  ///< Connectivity (neighbours per node); negative means fluctuating connectivity
  int  neighbour_size(id_t id) const ///Connectivity of node id
  {return connectivity>0 ? connectivity : Nneighbours[id];}

  ///@}
  /// \name Node access
  ///@{
  nodeT&       operator[](id_t i);            ///< Get reference to node `i`
  const nodeT& operator[](id_t i) const; ///< Get constant reference to node `i`
  nodeT&       at(id_t i);               ///< Reference to node, with bound check
  const nodeT& at(id_t i) const;         ///< Constant eference to node, with bound check
  id_t         id(node_t *n) const;           ///< Get node id from pointer to node
  nodeT*       data() noexcept {return nodes;}///< Access raw data through pointer
  const nodeT* data() const noexcept          /// Access raw data through `const` pointer
  {return nodes;}

  ///@}
  
  ///\name Disk i/o
  ///@{

  void read(std::istream&);  ///< Read graph from stream
  void write(std::ostream&); ///< Write graph to stream

  ///@}
  
  /// \name Iterators
  ///@{
  class   node_iterator;

  node_iterator      begin();
  nodeT*             end();   ///< node_iterator objects can be compared to pointers

  ///@}
  
protected:
  /// \name Construction, copy and destruction
  ///@{
  GraphBase();
  GraphBase(int N, int connectivity, bool use_neighbour_list);
  ~GraphBase() {cleanup();}
  GraphBase(const GraphBase&);
  GraphBase& operator=(const GraphBase&);

  ///@}

  /// \name For children
  ///@{
  void init_storage(int N,int connect,bool use_neighbour_list);
  void add_bond(id_t n1,id_t n2);
  /// Returns neighbour list for node (identified by pointer)
  nodeT **node_neighbours(nodeT *n) {return neighbours[n-nodes];}
  /// Returns neighbour list for node (identified by node id)
  nodeT **node_neighbours(id_t i) {return neighbours[i];}

  ///@}

private:
  int    Nnodes;
  nodeT* nodes;
  int    connectivity,*Nneighbours; // connectivity<0 means fluctuating
  nodeT  ***neighbours;

  void cleanup();
  void copy_into_empty(const GraphBase&);
} ;


/********************************************************************************
 *
 * Construction, copy and destruction
 * 
 */

/** \brief Constructs empty graph
 */
template <typename nodeT> inline
GraphBase<nodeT>::GraphBase() :
  nodes(0), Nnodes(0),
  neighbours(0), Nneighbours(0)
{}

/** \brief Constructs and reserves space for nodes
 *
 * \param N        Number of nodes
 *
 * \param connect Connectivity (number of neighbours per node).
 * Give -1 to indicate that connectivity fluctuates from node to node
 * 
 * \param use_neighbour_list If true, the neighbour list will be
 * initialized.  Use this if you plan to call add_bond(id_t,id_t).  If
 * false, the list is not initialized and add_bond() cannot be
 * called.  This is for regular graphs where neighbours are most
 * efficiently found through some mathematical rule, or for more
 * specialized data structures (as in periodic lattices).
 */
template <typename nodeT> inline
GraphBase<nodeT>::GraphBase(int N,int connect,bool use_neighbour_list) : 
  Nnodes(N),
  connectivity(connect),
  neighbours(0), 
  Nneighbours(0)
{
  init_storage(N,connect,use_neighbour_list);
}

/** \brief Initialize graph storage; arguments as in GraphBase::GraphBase(int,int,bool)
 */
template <typename nodeT>
void GraphBase<nodeT>::init_storage(int N,int connect,bool use_neighbour_list)
{
  Nnodes=N;
  connectivity=connect;
  nodes=new nodeT[Nnodes];
  if (use_neighbour_list) {
    Nneighbours=new int[Nnodes+1];  // to allow a pointer 1-past-end
    for (int i=0; i<Nnodes; i++) Nneighbours[i]=0;
    neighbours=new nodeT**[Nnodes];
    for (int i=0; i<Nnodes; i++) neighbours[i]=0;
  }
}
  
/** \brief Copy construction
*/
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

/* This implements the actual copy construction */
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

/* Implementation of destruction */
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

/******************************************************************************
 *
 * Methods for children
 *
 */
  
/** \brief Add a link between the named nodes
 *
 * This method is never called from GraphBond, it is to be used by
 * children to add links.  Note however that it updates the respective
 * nodes' neighbour lists, but no bond list is maintained.
 * GraphBondBase must override this to update its list of bonds.
 *
 */
template <typename nodeT>
void GraphBase<nodeT>::add_bond(id_t n1,id_t n2)
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

/*******************************************************************************
 *
 * Element access
 *
 */

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

/*******************************************************************************
 *
 * Disk i/o
 *
 */

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

/******************************************************************************
 *
 * Iterator
 *
 */

/** \brief Return a node_iterator (a topology-aware iterator).

This method returns an iterator that knows about neighbours, which as
a result may perform slightly worse than a pointer.  Use only if you
need to use the graph's topology; for bulk operations on the nodes
without regard to structure, it's better to acquire a pointer to the
raw storage through the data() method.

 */
template <typename nodeT> inline
typename GraphBase<nodeT>::node_iterator
GraphBase<nodeT>::begin()
{
  return node_iterator(*this);
}

/** \brief A one-past-end pointer, suitable to compare against node_iterator.
 */
template <typename nodeT> inline
nodeT* GraphBase<nodeT>::end()
{
  return nodes+Nnodes;
}

/**
 \brief A bidirectional iterator that knows about a site's neighbours.

 This is an iterator that knows about neighbours (see methods
 neighbour(), to_neighbour(), neighbour_size()).

 To define a specialized iterator it may be convenient to inherit from
 here, so we have kept the graph reference and node pointer protected.

 */
template <typename nodeT>
class GraphBase<nodeT>::node_iterator :
  public std::iterator<std::bidirectional_iterator_tag,nodeT>
{
public:
  ///\name Construction, copy and comparison
  //@{
  node_iterator(GraphBase<nodeT> &g,id_t ini_site=0) :
    graph(g), nn(g.nodes+ini_site) {}

  node_iterator(GraphBase<nodeT> &g,nodeT *n_) :
    graph(g), nn(n_) {}

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

  //@}
  
  ///\name Operators required by standard bidirectional iterators
  ///@{
  nodeT& operator*() const {return *nn;}
  nodeT* operator->() const {return nn;}

  node_iterator& operator++() {++nn; return *this;}
  node_iterator& operator++(int) {node_iterator c=*this; ++nn; return c;}

  node_iterator& operator--() {nn--; return *this;}
  node_iterator& operator--(int) {node_iterator c=*this; --nn; return c;}

  ///@}

  ///\name Extra methods (jumps and neighbour access)
  ///@{

  node_iterator& to(int n) {nn=graph.nodes+nn; return *this;}
  node_iterator& to(nodeT* n_) {nn=n_; return *this;}

  int            neighbour_size() const
  {return graph.neighbour_size(graph.id(nn));}
  node_iterator& to_neighbour(int n)
  {nn=graph.node_neighbours(nn)[n]; return *this;}

  nodeT& neighbour(int i) const {return *(graph.node_neighbours(nn)[i]);}

  operator nodeT*() const {return nn;}

  ///@}

protected:
  GraphBase<nodeT>& graph;
  nodeT *nn;
};

/******************************************************************************
 *
 * Algorithm functions
 *
 */


/* \internal Implementation of for_each_neighbour
 *
 * This is the actual code for_each_neighbour will run in the general
 * case.  It is placed somewhat awkwardly in a struct because we will
 * use partial template specialization to write code for certain types
 * of graphs, and partial specialization is only allowed for classes,
 * not function templates.
 *
 */
template <typename nodeT,typename Function,typename Iterator>
struct implement_for_each_neighbour {
  inline static Function fen(Iterator n,Function f)
  {
    for (int i=0; i<n.neighbour_size(); ++i)
      f(n.neighbour(i));
    return f;
  }
} ;

/** \ingroup lattice

 This is a STL-style `for_each` function that will visit all
 neighbours of a given node (identified through an iterator).  Use of
 this function allows optimizations (such as unrolling the loop over
 neighbours) for certain regular classes (see e.g. the periodic square
 lattice implementation).

 */
template <typename Iterator,typename Function>
inline Function for_each_neighbour(Iterator n,Function f)
{
  return implement_for_each_neighbour<
    typename std::iterator_traits<Iterator>::value_type,
    Function,Iterator>::fen(n,f);
}

/******************************************************************************
 *
 * Bonds and weights (not implemented)
 *
 */

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

} /* namespace glsim */

#endif /* _GRAPH_HH_ */
