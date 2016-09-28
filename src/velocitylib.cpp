#include <iostream>
#include <cmath>
#include <cstdlib>
#include <functional>

#include "numtools.hpp"
#include "velocitylib.hpp"
///class constructor
/*!\param v the horizontal wavespeed
 *\param e the anisotropy parameter
 */
Velocity::Velocity(const double& v, const double& e)
{
  vh = v;
  E = e;
}

///default class constructor
Velocity::Velocity()
{
  vh = 1.0;
  E = 0.0;
}

/// method to compute the phase velocity
/*!\param theta the phase angle
 *\return the phase velocity
 */
double Velocity::Vphase(const double& theta)
{
  double v;
  v = vh*(1 + E*std::pow(std::cos(theta), 2) );
  return v;
}

/// method to compute the tangent of the group angle from the tangent of the phase angle
/*!\param tan_theta the tangent of the phase angle
 *\return the tangent of the group angle
 */
double Velocity::tan_group_angle(const double& tan_theta)
{
  return tan_theta*(1 - E + tan_theta*tan_theta)/(1 + E + (1+2*E)*tan_theta*tan_theta);
}

/// method to compute the group angle from the phase angle
/*!\param theta the phase angle
 *\return the group angle
 */
double Velocity::group_angle(const double& theta)
{
  return std::atan(
		   tan_group_angle(std::tan(theta))
		   );
}

///method to compute the phase angle from the group angle
/*!There is no analytical formula to do the conversion, so a numerical approach is used. For simplicity and robustness I use Ridders' method. Important to note is that the input group angle \f$\phi\f$ is larger than \f$\pi/2\f$, a situation unphysical in principle but often met in real life because of the floating-point nature of the computation, this function will assume that \f$\phi=\pi/2\f$. Same thing for \f$\phi<-\pi/2\f$.
 *\param phi the group angle
 *\return the phase angle
 */
double Velocity::phase_angle(const double& phi)
{
  const double PI=3.141592653589793;
  double theta;
  double x1, x2;
  std::function<double(double)> f_PHI = std::bind( &Velocity::group_angle, this, std::placeholders::_1 );
    
  // initialise bounds
  if (phi>0)
    {
      x1 = 0.0;
      x2 = PI/2.0;
      // if the angle is larger than pi/2, assume it is equal to it
      // this is useful to interfacing with things like matlab or anything, because of the floating point nature of PI used here (nothing will ever be exactly equal to PI/2).
      if (phi>x2)
	return x2;
    }
  else if (phi<0)
    {
      x1 = -PI/2.0;
      x2 = 0.0;
      // same as above
      if (phi<x1)
	return x1;
    } 

  theta = numtools::fx_eq_x( f_PHI, phi, x1, x2 );

  return theta;

}

///method to compute the group velocity
/*!\param phi the group angle
 */
double Velocity::Vgroup(const double& phi)
{
  double V;
  double theta;

  // first compute theta from phi:
  theta = phase_angle(phi);

  // and then compute the group velocity
  V = Vphase(theta)*std::sqrt( 1+ E*E*std::sin(2*theta)*std::sin(2*theta)
			       /((1+E*std::cos(theta)*std::cos(theta))*(1+E*std::cos(theta)*std::cos(theta))) );
    
  return V;
}

