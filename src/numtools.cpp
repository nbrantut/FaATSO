#include <functional>
#include <cmath>
#include <iostream>
#include <cstdlib>

#include "numtools.hpp"

///the sign function
/*!\param x a real number
 *\return the sign of x
 */
int numtools::sign(const double& x)
{
  if (x>0.0)
    return 1;
  if (x<0.0)
    return -1;
  return 0;
}

///function to find the x such that f(x)=x, employing Ridder's method.
/*!\param f the function \f$f(x)\f$
 *\param x the variable to look for in \f$f(x)=x\f$
 *\param xl value of \f$x\f$ at which \f$f(x)<x\f$
 *\param xh value of \f$x\f$ at which \f$f(x)>x\f$
 *\return the zero of \f$f(x)-x\f$
 */
double numtools::fx_eq_x(const std::function<double(double)>& f, const double& x, const double& xl, const double& xh)
{
  const double tol=1e-9;
  double x1, x2, xm;
  double fx1, fx2, fxm;
  double xnew;
  double fxnew;
  double s;
  int it=0;
  int maxit=25;

  // assing values of end points
  x1 = xl;
  x2 = xh;

  // evaluate at end points
  fx1 = f(x1)-x;
  fx2 = f(x2)-x;
    
  // initialise solution
  xnew = -9.99e99;
  
  if ((fx1<0.0 && fx2>0.0) || (fx1>0.0 && fx2<0.0)) 
    {
      do
	{
	  // iterate
	  it++;
	  
	  // take midpoint
	  xm = 0.5*(x1+x2);
	  
	  // evaluate function
	  fxm = f(xm) - x;
	  
	  //update forumla
	  s = std::sqrt( fxm*fxm - fx1*fx2 );
	  if (s==0)
	    return xnew;
	  
	  xnew = xm + (xm-x1)*( sign(fx1-fx2)*fxm )/s;
	  
	  //evaluate at updated point
	  fxnew = f(xnew) - x;
	  
	  if (sign(fxnew) != sign(fxm))
	    {
	      x1 = xm;
	      fx1 = fxm;
	      x2 = xnew;
	      fx2 = fxnew;
	    }
	  else if (sign(fxnew) != sign(fx1))
	    {
	      x2 = xnew;
	      fx2 = fxnew;
	    }
	  else if (sign(fxnew) != sign(fx2) )
	    {
	      x1 = xnew;
	      fx1 = fxnew;
	    }
	} while ( (std::abs(fxnew)>tol) && (it<maxit));
      return xnew;
    }
  else
    {
      if (fx1==0.0) return x1;
      if (fx2==0.0) return x2;
      std::cout << "Error in fzero: zero not bracketed\n";
      return xnew;
    }
}

