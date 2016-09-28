//safeguard library
#ifndef VELOCITYLIB_INCLUDED
#define VELOCITYLIB_INCLUDED
//

/// the class Velocity corresponds to an anisotropic seismic velocity model, characterised by two independent parameters.
/*! the model is very simple, and is more phenomenological than purely physical.
 * the model is for transversely isotropic medium, and the phase velocity is elliptical:
 * \f$v(\theta) = v_\mathrm{h}(1 + E \cos^2(\theta))\f$
 * where \f$\theta\f$ is the phase angle, \f$v_\mathrm{h}\f$ is the horizontal wavespeed, and \f$E\f$ is the anisotropy parameter.
 */
class Velocity
{

public:

  double vh; ///< the horizontal wavespeed
  double E;  ///< the aniostropy parameter, such that the vertical wavespeed is equal to vh*(1+E)

  Velocity(const double& , const double& );
  Velocity();

  double Vphase(const double&);
  double tan_group_angle(const double& );
  double group_angle(const double& );
  double phase_angle(const double& );
  double Vgroup(const double& );

};

#endif
