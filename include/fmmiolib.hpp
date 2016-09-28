// safeguard library
#ifndef FMMIOLIB_INCLUDED
#define FMMIOLIB_INCLUDED
//

#include <vector>

//forward declare class grid
class Grid;

/// The class Data is used to load and save the binaries used for input/output in the fmm main code. It stores <double> values at every grid point, as well as the grid dimensions and spacing.

class Data
{
public:
  std::vector<double> dat; ///< a vector of data (will be of size Nx*Ny*Nz)
  int Nx; ///< number of nodes in x direction
  int Ny; ///< number of nodes in y direction
  int Nz; ///< number of nodes in z direction
  double h; ///< grid spacing

  Data();
  Data(const int&, const int&, const int&, const double&);
  void load(const char *) ;
  void save(const char *) ;
  Data refine(int);
  void get_interp_ind_coef(int*, double*, const int&, const int&, const int&, const int&, const int&);

private:
  void ind2sub(const int&, int&, int&, int&);
  void ind2sub(const int&, int&, int&, int&, const int&, const int&);
  int sub2ind(const int&, const int&, const int&);
};

namespace locations
{
  /// The class Source is used to store source points, with their coordinates and their corresponding linear index in a Grid.
  
  class Location
  {
  public:
    double x; ///< the x coordinate of the source
    double y; ///< the y coordinate of the source
    double z; ///< the z coordinate of the source
    double t0; ///< the origin time of the source
    int ind;  ///< the linear index of the source (i.e., the index of the closest node) in a Grid

    Location();

    Location(const double&, const double&, const double&, const double&);

    void findind(Grid&);
    void findind(const int&, const int&, const int&, const double&);
    void findxyz(Grid&);
  };

  std::vector<Location> readlocations(const char *);
  void writelocations(const char *, const std::vector<Location>);

}

namespace times
{
  std::vector<double> readtimes(const char*);
  void writetimes(const char*, const std::vector<double>);
}

#endif
