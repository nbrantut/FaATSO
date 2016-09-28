// safeguard library
#ifndef FMMLIB_INCLUDED
#define FMMLIB_INCLUDED
//

# include <array>
# include <vector>
# include <boost/heap/fibonacci_heap.hpp>

//forward declare class Data:
class Data;

/// The class Node corresponds to a node in a three dimensional orthonormal grid.
class Node
{
    
public:
  int index; ///< the index of the node in the grid, a unique identifier. It corresponds to the linear index in a 3D matrix, similar to Matlab linear indices.
  int i; ///< the subscript index along the x direction.
  int j; ///< the subscript index along the y direction.
  int k; ///< the subscript index along the z direction.

  /// an array of the 6 indices corresponding to the 6 neighbours of the node in the grid.
  /*! If the index is -1, it means that the neighbour does not exist (point outside the grid).
   * nhb[0] corresponds to the neighbour in +x
   * nhb[1] corresponds to the neighbour in -x
   * nhb[2] corresponds to the neighbour in +y
   * nhb[3] corresponds to the neighbour in -y
   * nhb[4] corresponds to the neighbour in +z
   * nhb[5] corresponds to the neighbour in -z
   */
  int nhb[6]; 
  double T; ///< the arrival time of the wavefront at the node.
  double V; ///< the horizontal wavespeed at this node.
  double E; ///< the anisotropy parameter (ratio of vertical to horizontal wavespeed).
  bool trial; ///< boolean tag for "trial" nodes.
  bool known; ///< boolean tag for "known" nodes.
  bool unknown; ///< boolean tag for "unknown nodes.

  Node();

  /// define comparison operator > for class "node" (compares the arrival time T)
  bool operator > (const Node& pt) const
  {
    return ( T > pt.T);
  }

  /// define comparison operator < for class "node" (compares the arrival time T)
  bool operator < (const Node& pt) const
  {
    return ( T < pt.T);
  }

};


/// The class Tetra embeds the information (velocity, anisotropy, arrival time of neighbours, etc) and methods required to compute the arrival time at a given node. It is named "Tetra" in reference to the tetrahedron that is formed by a given node of interest and three of its neighbours (each in x, y and z direction, resp.).
class Tetra
{
    
public:
  double vel; ///< the horizontal wavespeed at the node of interest (i,j,k).
  double eps; ///< the aniotropy parameter at the node of interest (i,j,k).
  double h; ///< grid spacing.
  bool test_a; ///< true is the x-neighbour (i+/-1,j,k) exists and is known.
  bool test_b; ///< true if the y-neighbour (i,j+/-1,k) exists and is known.
  bool test_c; ///< true if the z-neighbour (i,j,k+/-1) exists and is known. 
  double a; ///< arrival time at the x-neighbour.
  double ap; ///< arrival time at the x-neighbour of the x-neighbour (use for second order FD).
  double sw_a; ///< switch: equal to 1 if the x-neighbour of the x-neighbour exists and is knwon, 0 other wise. 
  double b; ///< same as for a, corresponding to the y-neighbour.
  double bp;  ///< same as for ap, corresponding to the y-neighbour.
  double sw_b; ///< same as for sw_a, corresponding to the y-neighbour.
  double c;  ///< same as for a, corresponding to the z-neighbour.
  double cp;  ///< same as for ap, corresponding to the z-neighbour.
  double sw_c;  ///< same as for sw_a, corresponding to the z-neighbour.
    
  Tetra(const double&, const double&, const double&);
  Tetra();
  double compute_T();
    
private:

  double threeD_aniso();
  double twoD_aniso_a();
  double twoD_aniso_b();

};


//type declaration:
typedef boost::heap::fibonacci_heap<Node, boost::heap::compare<std::greater<Node> > >::handle_type handle_t;


/// The class Grid corresponds to the set of nodes at which arrival time is to be computed, and contains the methods to do that.
class Grid
{

public:
  int Nx; ///< the number of nodes in the x direction
  int Ny; ///< the number of nodes in the y direction
  int Nz; ///< the number of nodes in the z direction
  int Ntot; ///< the total number of nodes
  double h; ///< the node spacing
  std::vector<Node> node; ///< the vector of nodes
  std::vector<handle_t> tab_handle_trial; ///< vector of handles in the trial heap
  int ind_source; ///< linear index of source node
  

  Grid(const int&, const int&, const int&, const double&);
  ~Grid();

  void initialise(const std::vector<double>& , const std::vector<double>& );
  void initialise_log(const std::vector<double>& , const std::vector<double>& );
  void reset_T();
  void reset(const std::vector<double>& , const std::vector<double>&);
  int march(const int& , const bool& );
  std::vector<double> export_Tvector();
  void import_T(const Data& );
  int sub2ind(const int& ,const int& , const int& );
  void ind2sub(const int& , int& , int& , int& );    
  
private:

  boost::heap::fibonacci_heap<Node, boost::heap::compare<std::greater<Node> > > trial_heap; ///< the heap structure in which "trial" nodes are pushed.

  Tetra buildtetra(const Node& , const int& ,const int&, const int&);
  std::array<Tetra,4> findalltetra(const int& ,const int& );
  double update(const int& , const int& );
  void init_box();
    
};

#endif
