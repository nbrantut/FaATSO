//safeguard library
#ifndef RAYLIB_INCLUDED
#define RAYLIB_INCLUDED
//

#include <vector>

//forward declare classes Grid and Data:
class Grid;
class Data;

///The class Ray corresponds to a seismic ray.
/*! A ray is traced a posteriori based on a regular Grid with computed arrival times. The tracing starts from the receiver (given as argument to the method Ray::trace), and progresses with steps equal to the grid spacing along the ray direction. 
 *In isotropic media, the ray direction at each piont is normal to the wavefront, i.e., along the gradient of the arrival time.
 *In anisotropic media, the ray direction is NOT normal to the wavefront, and at each position the ray angle is found based on a conversion between the phase angle and the group angle.
 */
class Ray
{

public:

  std::vector<double> x; ///< array of \f$x\f$ coordinates of ray path
  std::vector<double> y; ///< array of \f$y\f$ coordinates of ray path
  std::vector<double> z; ///< array of \f$z\f$ coordinates of ray path
  std::vector<int> ind; ///< array of indices of nearest neighbour nodes along the ray path
  std::vector<double> tan_theta; ///< array storing the tangent of the phase angle along the ray path
  int Nx; ///< Nx dimension of the grid in which the ray is traced.
  int Ny; ///< Ny dimension of the grid in which the ray is traced.
  int Nz; ///< Nz dimension of the grid in which the ray is traced.
  double h; ///< grid spacing.
  double der_x0; ///< derivative of arrival time with respect to x position of first ray point (i.e. target point)
  double der_y0; ///< derivative of arrival time with respect to y position of first ray point
  double der_z0; ///< derivative of arrival time with respect to z position of first ray point
  std::vector<double> der_v; ///< array of derivatives with respect to v=ln(V) (V being horizontal velocity)
  std::vector<double> der_E; ///< array of derivatives with respect to E

  Ray();
  void trace(Grid&, const double&, const double&, const double&);
  void export_ascii(const char *);
  void compute_interp_der(int, Data&, Data&);
  double compute_t(const Data&, const Data&);
  double compute_t_and_interp_der(const Data&, const Data&, int, Data&, Data&, double&, double&, double&);
 

private:

  double tan_group_angle(const double&, const double&);
  double F(const double&, const double&);
  double G(const double&, const double&);

};

#endif
