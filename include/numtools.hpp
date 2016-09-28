#ifndef NUMTOOLS_INCLUDED
#define NUMTOOLS_INCLUDED

#include <functional>

namespace numtools
{
  int sign(const double& );

  double fx_eq_x(const std::function<double(double)>& , const double&, const double& , const double&);

}

#endif
