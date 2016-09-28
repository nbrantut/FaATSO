#include <cmath>
#include <cstdlib>
#include <iostream>
#include <array>
#include <vector>
#include <fstream>
#include <limits>



#include "fmmlib.hpp"
#include "raylib.hpp"
#include "fmmiolib.hpp"

/// method to compute the ray path based on a Grid.
/*!\param grd the Grid in which arrival times T have been computed
 *\param xr \f$x\f$ coordinate of receiver point
 *\param yr \f$y\f$ coordinate of receiver point
 *\param zr \f$z\f$ coordinate of receiver point
 */
void Ray::trace(Grid& grd, const double& xr, const double& yr, const double& zr)
{
  int i_near, j_near, k_near;
  int ind_near;
  int i_fwd, i_bck, j_fwd, j_bck, k_fwd, k_bck;
  int ind_fwd, ind_bck;
  double x0, y0, z0;
  double xr_l, yr_l, zr_l;
  double xmax, ymax, zmax;
  unsigned n, maxn;
  double dist;
  double grad_x, grad_y, grad_z;
  double g_xy, norm_g;
  double f, dfde;
  const double small=1e-23;
  double coef=0.0;

  // (re)initialise ray attributes:
  Nx = grd.Nx;
  Ny = grd.Ny;
  Nz = grd.Nz;
  h = grd.h;
  x.resize(0);
  y.resize(0);
  z.resize(0);
  ind.resize(0);
  tan_theta.resize(0);
  der_v.resize(0);
  der_E.resize(0);

  
  
  //std::cout << "calling trace...\n";
  // set maximum number of elements along ray:
  maxn = 3*std::max(Nx, std::max(Ny, Nz) );

  // set max coordinates
  xmax = (Nx-1)*h;
  ymax = (Ny-1)*h;
  zmax = (Nz-1)*h;

  // assign coordinates of source point
  x0 = grd.node[grd.ind_source].i *h;
  y0 = grd.node[grd.ind_source].j *h;
  z0 = grd.node[grd.ind_source].k *h;

  // assign first coordinates
  // first check that receiver is not outside the grid:
  if ( (xr<0) || (xr>xmax) || (yr<0) || (yr>ymax) || (zr<0) || (zr>zmax) )
    {
      std::cout << "Receiver coordinates out of grid. Choosing nearest grid point.\n";
      xr_l = std::max( std::min(xr, xmax) , 0.0);
      yr_l = std::max( std::min(yr, ymax) , 0.0);
      zr_l = std::max( std::min(zr, zmax) , 0.0);
    }
  else
    {
      xr_l = xr;
      yr_l = yr;
      zr_l = zr;
    }

  // assign the first element of each coordinate vector
  x.push_back(xr_l);
  y.push_back(yr_l);
  z.push_back(zr_l);

  // set counter:
  n = x.size()-1;

  // compute distance of first ray point to source point
  dist = std::sqrt( (x[n]-x0)*(x[n]-x0) 
		    + (y[n]-y0)*(y[n]-y0) 
		    + (z[n]-z0)*(z[n]-z0) 
		    );

  //std::cout << "dist is " << dist << "\n";

  // find nearest neighbour in the grid:
  i_near = (int) round(x[n]/h);
  j_near = (int) round(y[n]/h);
  k_near = (int) round(z[n]/h);
    
  //std::cout << "(i,j,k)_near is (" << i_near << "," << j_near << "," << k_near << ")\n";

  ind_near = grd.sub2ind(i_near, j_near, k_near);
  // assign first element of nearest neighbour indices
  ind.push_back(ind_near);

  //std::cout << "ind_near is " << ind_near << "\n";


  //std::cout << "start propagating ray...\n";
  while (dist>h && n<maxn)
    {
      // find the subscript indices of neighbouring points
      i_fwd = std::min(i_near+1, Nx-1);
      i_bck = std::max(i_near-1, 0);
	    
      j_fwd = std::min(j_near+1, Ny-1);
      j_bck = std::max(j_near-1, 0);

      k_fwd = std::min(k_near+1, Nz-1);
      k_bck = std::max(k_near-1, 0);

      //      std::cout << "i: (" <<i_fwd << "," << i_bck << ")\n";	    
      //      std::cout << "j: (" <<j_fwd << "," << j_bck << ")\n";
      //      std::cout << "k: (" <<k_fwd << "," << k_bck << ")\n";

      //std::cout << "compute gradient...\n";
      // compute gradient of arrival time in x direction
      ind_fwd = grd.sub2ind(i_fwd, j_near, k_near);
      ind_bck = grd.sub2ind(i_bck, j_near, k_near);
      grad_x = (grd.node[ind_fwd].T - grd.node[ind_bck].T)
	/( (i_fwd - i_bck)*h + small);

      // compute gradient of arrival time in y direction
      ind_fwd = grd.sub2ind(i_near, j_fwd, k_near);
      ind_bck = grd.sub2ind(i_near, j_bck, k_near);
      grad_y = (grd.node[ind_fwd].T - grd.node[ind_bck].T)
	/( (j_fwd - j_bck)*h + small);
	    
      // compute gradient of arrival time in z direction
      ind_fwd = grd.sub2ind(i_near, j_near, k_fwd);
      ind_bck = grd.sub2ind(i_near, j_near, k_bck);
      grad_z = (grd.node[ind_fwd].T - grd.node[ind_bck].T)
	/( (k_fwd - k_bck)*h + small);

      //std::cout << "compute ray angle...\n";
      // modify the z component of gradient to obtain the direction
      // of the ray propagation according to the group angle	
      if (grad_z!=0)
	{
	  // assign dimension of projected gradient on xy plane for easier use:
	  g_xy = std::sqrt(grad_x*grad_x + grad_y*grad_y);
	  
	  // tangent of phase angle:
	  tan_theta.push_back(g_xy/grad_z);

	  // compute the correct gradient using group angle:
	  if (g_xy>0.0)
	    {
	      grad_z = g_xy/tan_group_angle(tan_theta.back(), grd.node[ind_near].E);
	    }
	  

	  // compute useful values to get the derivatives of arrival time with respect to parameters:
	  f = F(tan_theta.back(),grd.node[ind_near].E);
	  dfde = G(tan_theta.back(),grd.node[ind_near].E)/(2.0*f);
	}
      else
	{
	  f = 1.0;
	  dfde = 0.0;
	  tan_theta.push_back(std::numeric_limits<double>::infinity());
	}
      
      // compute normalisation factor for the gradient
      norm_g = std::sqrt(grad_x*grad_x + grad_y*grad_y + grad_z*grad_z);
	    
      //      std::cout << "grad is (" << grad_x << "," << grad_y << "," << grad_z << ")\n";

      // update ray
      x.push_back( x[n] - h*grad_x/norm_g );
      y.push_back( y[n] - h*grad_y/norm_g );
      z.push_back( z[n] - h*grad_z/norm_g );

      //compute the derivatives with respect to the local v=ln(V) and E:
      //std::cout << "pushback der_v\n";
      der_v.push_back( -h/(grd.node[ind_near].V * f) );
      //std::cout << "der_v ok\n";
      der_E.push_back( der_v.back() *dfde/f);
      //std::cout << "der_E ok\n";
	    
      // update counter
      n = x.size()-1;

      // make sure the current ray point remains in the grid!
      if ( (x[n]<0) || (x[n]>xmax) || (y[n]<0) || (y[n]>ymax) || (z[n]<0) || (z[n]>zmax) )
	{
	  std::cout << "Current ray point outside grid: forcing it back in.\n";
	  x[n] = std::max( std::min( x[n], xmax), 0.0);
	  y[n] = std::max( std::min( y[n], ymax), 0.0);
	  z[n] = std::max( std::min( z[n], zmax), 0.0);
	}

      // compute distance of new ray point to the source point
      dist = std::sqrt( (x[n]-x0)*(x[n]-x0) 
			+ (y[n]-y0)*(y[n]-y0) 
			+ (z[n]-z0)*(z[n]-z0) 
			);

      // find nearest neighbour in the grid:
      i_near = (int) round(x[n]/h);
      j_near = (int) round(y[n]/h);
      k_near = (int) round(z[n]/h);

      // compute linear index
      ind_near = grd.sub2ind(i_near, j_near, k_near);

      // add this index to the ray.ind property
      ind.push_back(ind_near);    
    }

  // we have everything now, but the last pont on the ray is not associated with any derivative, so push back zeros to make der_v and der_E the same size as x,y,z etc:
  der_v.push_back(0.0);
  der_E.push_back(0.0);
  // same goes for the tangent of the phase angle:
  tan_theta.push_back(0.0);

  // now we have the full ray computed. last thing is to compute the derivatives of the arrival time with respect to teh x,y,z positions of the first ray point (the first ray point here is considered as the seismic source of the event: the arrival times are computed using the sensor as a virtual source, taking advantage of the reversibility of the wave equation).

  // first check that there is at least one point in the der_v vector. If yes, use that for the value of the coefficient -h/V. If no, there is no ray (empty) so the coefficient is anyways zero.
  if (der_v.size()>0)
    coef = der_v[0];
  
  der_x0 = coef*(x[1] - x[0])/(h*h);
  der_y0 = coef*(y[1] - y[0])/(h*h);
  der_z0 = coef*(z[1] - z[0])/(h*h);
  
}

void Ray::export_ascii(const char *filename)
{
  std::ofstream file;
  unsigned n=0;

  file.open(filename);

  file << "index,x,y,z,der_v,der_E\n";

  while ( n<x.size() )
    {
      file << ind[n]<< "," << x[n] << "," << y[n] << "," << z[n] << "," << der_v[n] << "," << der_E[n] << "\n";
      n++;
    }

  file.close();

}

/// method to extract the derivatives of the arrival time with respect to the (log(V)) or E at each node, on a coarser grid (by a factor n in all dimensions) than the one used to compute the ray initially. The values of log(V) and E on a fine grid are supposed to correspond to a trilinear interpolation of the values on the coarse grid.
/*!\param n the refinement factor.
 *\param ddv a Data structure which will contain the interpolated derivatives with respect to v=lnV.
 *\param ddE a Data structure which will contain the interpolated derivatives with respect to E.
 */
void Ray::compute_interp_der(int n, Data& ddv, Data& ddE)
{
  int I[8], i_000, i_100, i_010, i_110, i_001, i_101, i_011, i_111;
  double xyz[3], xx, omx, yy, omy, zz, omz;
  double a_000, a_100, a_010, a_110, a_001, a_101, a_011, a_111;

  // initialise derivatives data structures:
  ddv.Nx = (Nx-1)/n + 1;
  ddv.Ny = (Ny-1)/n + 1;
  ddv.Nz = (Nz-1)/n + 1;
  ddv.h = h*n;
  ddv.dat.assign(ddv.Nx*ddv.Ny*ddv.Nz, 0.0);

  ddE = ddv;
  
  for (int p=0; p<ind.size(); p++)
    {
      ddv.get_interp_ind_coef(I,xyz,ind[p],Nx,Ny,Nz,n);

      xx = xyz[0];
      omx = 1.0 - xx;
      yy = xyz[1];
      omy = 1.0 - yy;
      zz = xyz[2];
      omz = 1.0 - zz;

      i_000 = I[0];
      i_100 = I[1];
      i_010 = I[2];
      i_110 = I[3];
      i_001 = I[4];
      i_101 = I[5];
      i_011 = I[6];
      i_111 = I[7];

      a_000 = omx*omy*omz;
      a_100 = xx*omx*omy;
      a_010 = omx*yy*omz;
      a_110 = xx*yy*omz;
      a_001 = omx*omy*zz;
      a_101 = xx*omy*zz;
      a_011 = omx*yy*zz;
      a_111 = xx*yy*zz;

      ddv.dat[i_000] += a_000*der_v[p];
      ddv.dat[i_100] += a_100*der_v[p];
      ddv.dat[i_010] += a_010*der_v[p];
      ddv.dat[i_110] += a_110*der_v[p];
      ddv.dat[i_001] += a_001*der_v[p];
      ddv.dat[i_101] += a_101*der_v[p];
      ddv.dat[i_011] += a_011*der_v[p];
      ddv.dat[i_111] += a_111*der_v[p];

      ddE.dat[i_000] += a_000*der_E[p];
      ddE.dat[i_100] += a_100*der_E[p];
      ddE.dat[i_010] += a_010*der_E[p];
      ddE.dat[i_110] += a_110*der_E[p];
      ddE.dat[i_001] += a_001*der_E[p];
      ddE.dat[i_101] += a_101*der_E[p];
      ddE.dat[i_011] += a_011*der_E[p];
      ddE.dat[i_111] += a_111*der_E[p];
    }

  return;

}

/// method to compute a travel time along the ray, given a velocity and anisotropy models (given on the same grid as the one used to initialise the ray).
/*! This essentially computes \f$t = \int_{Ray} (1/V_{group})ds.\f$ This is useful to recompute arrival times in a new velocity model assuming the ray path is fixed (e.g., avoiding to recompute the whole times in a new grid).
 *\param tab_v Data structure containing the ln(V) data,
 *\param tab_E Data structure containing the E data.
 *\return the travel time.
 */
double Ray::compute_t(const Data& tab_v, const Data& tab_E)
{
  double tt=0.0;
  int p, k;
  double Vg;

  if (tab_v.h != h || tab_v.Nx!= Nx || tab_v.Ny != Ny || tab_v.Nz != Nz ||
      tab_E.h != h || tab_E.Nx!= Nx || tab_E.Ny != Ny || tab_E.Nz != Nz )
    {
      std::cout << "Trying to recompute arrival time based on incompatible grids... ABORT.\n";
      exit(1);
    }

  //integrate using trapezoidal rule. Take care of end points first:
  k = ind[0];
  Vg = std::exp(tab_v.dat[k]) * F(tan_theta[0], tab_E.dat[k]);
  tt += 0.5*h/Vg;

  k = ind[ind.size()-1];
  Vg = std::exp(tab_v.dat[k]) * F(tan_theta[ind.size()-1], tab_E.dat[k]);
  tt += 0.5*h/Vg;
  
  for (p=1; p<ind.size()-1; p++)
    {
      k = ind[p];
      Vg = std::exp(tab_v.dat[k]) * F(tan_theta[p], tab_E.dat[k]);
      tt += h/Vg;
    }
  return tt;
}

/// method to compute a travel time along the ray, given a velocity and anisotropy models (given on the same grid as the one used to initialise the ray), as well as the derivatives of the arrival time with respect to lnV and E, computed on a coarser grid (n times less nodes in each directions).
/*! This essentially computes \f$t = \int_{Ray} (1/V_{group})ds.\f$ and \f$dt/dv, dt/dE\f$. This is useful to recompute arrival times and frechet kernels in a new velocity model assuming the ray path is fixed (e.g., avoiding to recompute the whole times in a new grid).
 *\param tab_v Data structure containing the ln(V) data,
 *\param tab_E Data structure containing the E data.
 *\param n the grid refinement factor.
 *\param ddv Data structure to store the values of dt/dv (v on coraser grid).
 *\param ddE Data structure to store the values of dt/dE (E on coarser grid).
 *\param ddx the value (double) of dt/dx0.
 *\param ddy the value (double) of dt/dy0.
 *\param ddz the value (double) of dt/dz0.
 *\return the travel time.
 */
double Ray::compute_t_and_interp_der(const Data& tab_v, const Data& tab_E, int n, Data& ddv, Data& ddE, double& ddx, double& ddy, double& ddz)
{
  double tt=0.0, w;
  int p, k;
  double Vg,V,f,dfde,tmp_der_v,tmp_der_E;
  int I[8], i_000, i_100, i_010, i_110, i_001, i_101, i_011, i_111;
  double xyz[3], xx, omx, yy, omy, zz, omz;
  double a_000, a_100, a_010, a_110, a_001, a_101, a_011, a_111;

  // first check if grid sizes are OK:
  if (tab_v.h != h || tab_v.Nx!= Nx || tab_v.Ny != Ny || tab_v.Nz != Nz ||
      tab_E.h != h || tab_E.Nx!= Nx || tab_E.Ny != Ny || tab_E.Nz != Nz )
    {
      std::cout << "Trying to recompute arrival time based on incompatible grids... ABORT.\n";
      exit(1);
    }

  // now perform the computation of travel time and derivatives along the ray.
  // initialise derivatives data structures:
  ddv.Nx = (Nx-1)/n + 1;
  ddv.Ny = (Ny-1)/n + 1;
  ddv.Nz = (Nz-1)/n + 1;
  ddv.h = h*n;
  ddv.dat.assign(ddv.Nx*ddv.Ny*ddv.Nz, 0.0);

  ddE = ddv;
  
  for (p=0; p<ind.size(); p++)
    {
      // short hand notation for the index:
      k = ind[p];
      // and for horizontal velocity:
      V = std::exp(tab_v.dat[k]);
      // and for angular dependencies:
      f = F(tan_theta[p], tab_E.dat[k]);
      dfde = G(tan_theta[p], tab_E.dat[k])/(2.0*f);
      
      // group velocity:
      Vg = V * f;

      // take care of first and last points along the ray:
      if (p==0)
	{
	  // assign weight w to compute integral over raypath
	  w = 0.5;
	  // assign derivatives with respect to source position
	  ddx = -(x[1]-x[0])/(h*Vg);
	  ddy = -(y[1]-y[0])/(h*Vg);
	  ddz = -(z[1]-z[0])/(h*Vg);
	}
      else if (p==ind.size()-1)
	w= 0.5;
      else
	w = 1.0;
      
      // compute the travel time increment:
      tt += w*h/Vg;

      if (p<ind.size()-1)
	{
	  // compute derivative with respect to v:
	  tmp_der_v = -h/(V * f);
	  // and with respect to E
	  tmp_der_E = tmp_der_v *dfde/f;
	  
	  // now interpolate the derivatives on the coarse grid (refinement factor n):
	  ddv.get_interp_ind_coef(I,xyz,k,Nx,Ny,Nz,n);
	  
	  xx = xyz[0];
	  omx = 1.0 - xx;
	  yy = xyz[1];
	  omy = 1.0 - yy;
	  zz = xyz[2];
	  omz = 1.0 - zz;
	  
	  i_000 = I[0];
	  i_100 = I[1];
	  i_010 = I[2];
	  i_110 = I[3];
	  i_001 = I[4];
	  i_101 = I[5];
	  i_011 = I[6];
	  i_111 = I[7];
	  
	  a_000 = omx*omy*omz;
	  a_100 = xx*omx*omy;
	  a_010 = omx*yy*omz;
	  a_110 = xx*yy*omz;
	  a_001 = omx*omy*zz;
	  a_101 = xx*omy*zz;
	  a_011 = omx*yy*zz;
	  a_111 = xx*yy*zz;
	  
	  ddv.dat[i_000] += a_000*tmp_der_v;
	  ddv.dat[i_100] += a_100*tmp_der_v;
	  ddv.dat[i_010] += a_010*tmp_der_v;
	  ddv.dat[i_110] += a_110*tmp_der_v;
	  ddv.dat[i_001] += a_001*tmp_der_v;
	  ddv.dat[i_101] += a_101*tmp_der_v;
	  ddv.dat[i_011] += a_011*tmp_der_v;
	  ddv.dat[i_111] += a_111*tmp_der_v;
	  
	  ddE.dat[i_000] += a_000*tmp_der_E;
	  ddE.dat[i_100] += a_100*tmp_der_E;
	  ddE.dat[i_010] += a_010*tmp_der_E;
	  ddE.dat[i_110] += a_110*tmp_der_E;
	  ddE.dat[i_001] += a_001*tmp_der_E;
	  ddE.dat[i_101] += a_101*tmp_der_E;
	  ddE.dat[i_011] += a_011*tmp_der_E;
	  ddE.dat[i_111] += a_111*tmp_der_E;
	}
    }
  return tt;

}


/// default class constructor
Ray::Ray()
{
  Nx = Ny = Nz = 0;
  h = 0.0;
}

/// method to compute the tangent of the group angle from the tangent of the phase angle.
/*! this function is specific to the type of anisotropy dealt with here (alliptical).
 *\param tan_phase_angle the tangent of the phase angle
 *\param E the anisotropy parameter
 *\return the tangent of the group angle
 */
double Ray::tan_group_angle(const double& tan_phase_angle, const double& E)
{
  return tan_phase_angle*(1 - E + tan_phase_angle*tan_phase_angle)/(1 + E + (1+2*E)*tan_phase_angle*tan_phase_angle);
}

/// method to compute the angular dependency of the group velocity, based on the phase angle and the anisotropy parameter.
/*! The function F is given by: \f$V_g = V_h\times F(\theta,E)\f$, where \f$V_h\f$ is the horizontal velocity.
 *\param tan_phase_angle the tangent of the phase angle
 *\param E the anisotropy parameter
 *\return the value of \f$F(\theta,E)\f$
 */
double Ray::F(const double& tan_phase_angle, const double& E)
{
  double cossq, sinsq;
  //  std::cout << "compute F...\n";
  cossq = 1.0/(1.0 + tan_phase_angle*tan_phase_angle);
  sinsq = 1 - cossq;
  return std::sqrt(1 + 2.0*E*cossq + E*E*cossq*(1.0+3.0*sinsq));
}

/// method to compute an intermediate results required for the derivative of the arrival time with respect to the anisotropy parameter, based on the phase angle and the anisotropy parameter.
/*! The function G is given by: \f$dF/dE = G(\theta,E)/(2F)\f$.
 *\param tan_phase_angle the tangent of the phase angle
 *\param E the anisotropy parameter
 *\return the value of \f$G(\theta,E)\f$
 */
double Ray::G(const double& tan_phase_angle, const double& E)
{
  double cossq, sinsq;
  //  std::cout << "compute G...\n";
  cossq = 1.0/(1.0 + tan_phase_angle*tan_phase_angle);
  sinsq = 1 - cossq;
  return 2.0*cossq + 2.0*E*cossq*(1.0+3.0*sinsq); 
}
