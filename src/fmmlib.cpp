# include <iostream>
# include <cstdlib>
# include <cmath>
# include <limits>
# include <array>
# include <vector>
# include <algorithm>
# include <boost/heap/fibonacci_heap.hpp>

# include "fmmlib.hpp"
# include "fmmiolib.hpp"
# include "velocitylib.hpp"

/// constructor for class Node. Sets index and coordinates to 0, T to infinity, all neighbour indices to -1, unknown to true, and trial and known to false.
Node::Node()
{
  int p;
  index = 0;
  i = j = k = 0;
  for (p=0; p<6; p++)
    nhb[p]=-1;
  T = std::numeric_limits<double>::infinity();
  V = E = 1;
  trial = known = false;
  unknown = true;
}

    
/// Class constructor
/*!
 * \param hh the grid spacing (h)
 * \param vv the horizontal velocity (vel)
 * \param ee the anisotropy parameter (eps)
 */
Tetra::Tetra(const double& hh, const double& vv, const double& ee)
{
  h = hh;
  vel = vv;
  eps = ee;
  a = -1;
  ap = -1;
  sw_a = 0;
  b = -1;
  bp = -1;
  sw_b = 0;
  c = -1;
  cp = -1;
  sw_c = 0;
  test_a = false;
  test_b = false;
  test_c = false;
}
    
/// Class default constructor
Tetra::Tetra()
{
  h = 1;
  vel = 1;
  eps = 0;
  a = -1;
  ap = -1;
  sw_a = 0;
  b = -1;
  bp = -1;
  sw_b = 0;
  c = -1;
  cp = -1;
  sw_c = 0;
  test_a = false;
  test_b = false;
  test_c = false;
}

/// method to compute the arrival time T at the node of interest.
/*! It checks which neighbour is known and uses the appropriate method to compute the time based on that.
 */
double Tetra::compute_T()
{
  double Tbar;
  
  // initialise Tbar
  Tbar = std::numeric_limits<double>::infinity();
        
        
  // compute the time based on the known values
  if (test_a && test_b && test_c)
    {
      // compute T from the full 3d fd
      Tbar = threeD_aniso();
      return Tbar;
    }
        
  if (test_a && test_b)
    {
      // there is the analytical solution here since the medium is isotropic in the xy plane:
      double delta;
      delta = ( (1+sw_a*0.5)*(1+sw_a*0.5) + (1+sw_b*0.5)*(1+sw_b*0.5) )*h*h/(vel*vel)
	- std::pow(
		   (1+sw_b*0.5)*(a + sw_a*(a-ap*0.5)) - (1+sw_a*0.5)*(b + sw_b*(b-bp*0.5))
		   ,2);

      if (delta<0)
	{
	  std::cout << "Warning: Second order method failed. Revert ot first order.\n";
	  sw_a = 0;
	  sw_b = 0;
	  delta = 2.0*h*h/(vel*vel) - (a - b)*(a - b);
	}

      Tbar = ( (1+sw_a*0.5)*(a + sw_a*(a-ap*0.5)) + (1+sw_b*0.5)*(b + sw_b*(b-bp*0.5))
	       + std::sqrt(delta) )
	/ ( (1+sw_a*0.5)*(1+sw_a*0.5) + (1+sw_b*0.5)*(1+sw_b*0.5) );
	  
	
      //std::cout << "vel   = " << vel << "\n";
	

      if (delta<0)
	{
	  //std::cout << "delta = " << delta << "\n";
	  //std::cout << "Tbar = " << Tbar << "\n";
	  //exit(1);
	  std::cout << "merde\n";
	}
            
      return Tbar;
    }
        
  if (test_a && test_c)
    {
      // compute T from the 2d fd
      Tbar = twoD_aniso_a();
      return Tbar;
    }
        
  if (test_b && test_c)
    {
      // compute T from the 2d fd
      Tbar = twoD_aniso_b();
      return Tbar;
    }
        
  if (test_a)
    {
      // compute T from 1d fd
      Tbar = (1./(1.+sw_a*0.5))*(h/vel + a*(1.+sw_a) - ap*sw_a*0.5);
      return Tbar;
    }
        
  if (test_b)
    {
      // compute T from 1d fd
      Tbar = (1./(1.+sw_b*0.5))*(h/vel + b*(1.+sw_b) - bp*sw_b*0.5);
      return Tbar;
    }
        
  if (test_c)
    {
      // compute T from 1d fd
      Tbar = (1./(1.+sw_c*0.5))*(h/(vel*(1.+eps)) + c*(1.+sw_c) - cp*sw_c*0.5);
      return Tbar;
    }
  exit(1);
}
    

/// Method to compute the arrival time at the node of interest based on three known neighbours.
/*! No check of the tests (a, b or c) is made so this function is private.
 * Because of the anisotropy, the computation uses a Newton-Raphson method to solve for the arrival time. 
 * If the computation fails (too many iterations) then the returned value if infinite (so that it won't be picked as something OK by the fmm).
 */
double Tetra::threeD_aniso()
{
  int it=1;
  const int maxit=15;
  const double tol = 1e-9;
  double x, f, fprime, dx;
  double V2h2, xa, xb, xc, xa2, xb2, xc2;

  // compute useful shortcut values
  V2h2 = (vel*vel/(h*h));
        
  // initialise first guess
  x = std::max( a, std::max(b,c) ) + h/vel;
        
  // compute shortcuts
  xa = x-a + 0.5*sw_a*(x - 2*a + ap);
  xb = x-b + 0.5*sw_b*(x - 2*b + bp);
  xc = x-c + 0.5*sw_c*(x - 2*c + cp);
  xa2 = std::pow(xa,2);
  xb2 = std::pow(xb,2);
  xc2 = std::pow(xc,2);
        
  // residual (function to find the zero of) :
  f = V2h2*( xa2 +xb2 + xc2 ) * std::pow(1 + eps* xc2/(xa2 + xb2 + xc2 )  , 2)  - 1;
        

  while ( (std::abs(f)>tol) && it<maxit )
    {
      // iterate
      it = it+1;


      // if we reach max. number of iterations, it means that the second order method does not work, so switch back to first order and start again:
      if (it==maxit)
	{
	  std::cout << "Warning: 2nd order scheme failed. Reverting to first order.\n";
	  it = 1;
	  x = std::max( a, std::max(b,c) ) + h/vel;
	  sw_a = 0;
	  sw_b = 0;
	  sw_c = 0;
	}
      else
	{            
	  // derivative with respect to x:
	  fprime = V2h2*(
			 ( 2*xa*(1+sw_a*0.5) + 2*xb*(1+sw_b*0.5) + 2*xc*(1+sw_c*0.5) )*std::pow( 1 + eps* xc2 /( xa2 + xb2 + xc2 ) ,2 )
			 + (xa2 + xb2 + xc2 )* 2*(1 + eps* xc2 /( xa2 + xb2 + xc2 )  )
			 * (eps) *( 2*xc*(1+sw_c*0.5) / ( xa2 + xb2 + xc2 )
				    -  ( 2*xa*(1+sw_a*0.5) + 2*xb*(1+sw_b*0.5) + 2*xc*(1+sw_c*0.5) ) * xc2 / std::pow( xa2 + xb2 + xc2 , 2) )
			 );
            
	  // compute change in x
	  dx = -f/fprime;
            
	  // update x
	  x = x+dx;
	}

      // compute updated shortcuts
      xa = x-a + 0.5*sw_a*(x - 2*a + ap);
      xb = x-b + 0.5*sw_b*(x - 2*b + bp);
      xc = x-c + 0.5*sw_c*(x - 2*c + cp);
      xa2 = std::pow(xa,2);
      xb2 = std::pow(xb,2);
      xc2 = std::pow(xc,2);
            
      // compute updated residual:
      f = V2h2*( xa2 + xb2 + xc2 ) * std::pow(1 + eps* xc2/(xa2 + xb2 + xc2 )  , 2)  - 1;

    }
        
  if (it==maxit)
    {
      // it means the algorithm did not converge, so there is no zero of f
      std::cout << "merde\n";	
      return std::numeric_limits<double>::infinity();
      //exit(1);
    }
  else
    {
      return x;
    }
}
    
/// Method to compute the arrival time at the node of interest based on thw known neighbours (a and c).
/*! No check of the tests (a, b or c) is made so this function is private.
 * Because of the anisotropy, the computation uses a Newton-Raphson method to solve for the arrival time.
 * If the computation fails (too many iterations) then the returned value if infinite (so that it won't be picked as something OK by the fmm).
 */
double Tetra::twoD_aniso_a()
{
  int it=1;
  const int maxit=15;
  const double tol = 1e-9;
  double x, f, fprime, dx;
  double V2h2, xa, xc, xa2, xc2;
        
  // compute useful shortcut values
  V2h2 = (vel*vel/(h*h));
        
  // initialise first guess
  x = std::max( a, c ) + h/vel;
        
  // compute shortcuts
  xa = x-a + 0.5*sw_a*(x - 2*a + ap);
  xc = x-c + 0.5*sw_c*(x - 2*c + cp);
  xa2 = std::pow(xa,2);
  xc2 = std::pow(xc,2);
        
  // residual (function to find the zero of) :
  f = V2h2*( xa2 + xc2 ) * std::pow(1 + eps* xc2/(xa2 + xc2 )  , 2)  - 1;
    
  while ( (std::abs(f)>tol) && (it<maxit) )
    {
      // iterate
      it++;

      // if we reach max. number of iterations, it means that the second order method does not work, so switch back to first order and start again:
      if (it==maxit)
	{
	  std::cout << "Warning: 2nd order scheme failed. Reverting to first order.\n";
	  it = 1;
	  x = std::max( a,c ) + h/vel;
	  sw_a = 0;
	  sw_c = 0;
	}
      else
	{
            
	  // derivative with respect to x:
	  fprime = V2h2*(
			 ( 2*xa*(1+sw_a*0.5) + 2*xc*(1+sw_c*0.5) )*std::pow( 1 + eps* xc2 /( xa2 + xc2 ) ,2 )
			 + (xa2 + xc2 )* 2*(1 + eps* xc2 /( xa2 + xc2 )  )
			 * (eps) *( 2*xc*(1+sw_c*0.5) / ( xa2 + xc2 )
				    -  ( 2*xa*(1+sw_a*0.5) + 2*xc*(1+sw_c*0.5) ) * xc2 / std::pow( xa2 + xc2 , 2) )
			 );
            
	  // compute change in x
	  dx = -f/fprime;
            
	  // update x
	  x += dx;
	}

      // compute updated shortcuts
      xa = x-a + 0.5*sw_a*(x - 2*a + ap);
      xc = x-c + 0.5*sw_c*(x - 2*c + cp);
      xa2 = std::pow(xa,2);
      xc2 = std::pow(xc,2);
            
      // compute updated residual:
      f = V2h2*( xa2 + xc2 ) * std::pow(1 + eps* xc2/(xa2 + xc2 )  , 2)  - 1;
	
    }
        
  if (it==maxit)
    {
      // it means the algorithm did not converge, so there is no zero of f
      std::cout << "merde\n";	
      return std::numeric_limits<double>::infinity();
      //exit(1);
    }
  else
    return x;
}
    
/// Method to compute the arrival time at the node of interest based on two known neighbours (b and c).
/*! No check of the tests (a, b or c) is made so this function is private.
 * Because of the anisotropy, the computation uses a Newton-Raphson method to solve for the arrival time. 
 * If the computation fails (too many iterations) then the returned value if infinite (so that it won't be picked as something OK by the fmm).
 */
double Tetra::twoD_aniso_b()
{
  int it=1;
  const int maxit=15;
  const double tol = 1e-9;
  double x, f, fprime, dx;
  double V2h2, xb, xc, xb2, xc2;
        
  // compute useful shortcut values
  V2h2 = (vel*vel/(h*h));
        
  // initialise first guess
  x = std::max(b,c) + h/vel;
        
  // compute shortcuts
  xb = x-b + 0.5*sw_b*(x - 2*b + bp);
  xc = x-c + 0.5*sw_c*(x - 2*c + cp);
  xb2 = std::pow(xb,2);
  xc2 = std::pow(xc,2);
        
  // residual (function to find the zero of) :
  f = V2h2*( xb2 + xc2 ) * std::pow(1 + eps* xc2/(xb2 + xc2 )  , 2)  - 1;

  //std::cout << b << " " << sw_b  << " " << bp  << " " << c  << " " << sw_c  << " " << cp << "\n";
  //std::cout << f << "\n";
        
  while ( (std::abs(f)>tol) && (it<maxit) )
    {
      // iterate
      it++;

      // if we reach max. number of iterations, it means that the second order method does not work, so switch back to first order and start again:
      if (it==maxit)
	{
	  std::cout << "Warning: 2nd order scheme failed. Reverting to first order.\n";
	  it = 1;
	  x = std::max(b,c) + h/vel;
	  sw_b = 0;
	  sw_c = 0;
	}
      else
	{
            
	  // derivative with respect to x:
	  fprime = V2h2*(
			 ( 2*xb*(1+sw_b*0.5) + 2*xc*(1+sw_c*0.5) )*std::pow( 1 + eps* xc2 /( xb2 + xc2 ) ,2 )
			 + (xb2 + xc2 )* 2*(1 + eps* xc2 /( xb2 + xc2 )  )
			 * (eps) *( 2*xc*(1+sw_c*0.5) / (xb2 + xc2 )
				    -  ( 2*xb*(1+sw_b*0.5) + 2*xc*(1+sw_c*0.5) ) * xc2 / std::pow( xb2 + xc2 , 2) )
			 );
            
	  // compute change in x
	  dx = -f/fprime;
            
	  // update x
	  x += dx;
	}
	
      // compute updated shortcuts
      xb = x-b + 0.5*sw_b*(x - 2*b + bp);
      xc = x-c + 0.5*sw_c*(x - 2*c + cp);
      xb2 = std::pow(xb,2);
      xc2 = std::pow(xc,2);
            
      // compute updated residual:
      f = V2h2*( xb2 + xc2 ) * std::pow(1 + eps* xc2/( xb2 + xc2 )  , 2)  - 1;
      //std::cout << f << "\n";
    }
        
  if (it==maxit)
    // it means the algorithm did not converge, so there is no zero of f
    {
      std::cout << "merde\n";
      return std::numeric_limits<double>::infinity();
    }
  else
    return x;
}



/// Class constructor
/*! \param nx the number of nodes in x direction (Nx)
 * \param ny the number of nodes in y direction (Ny)
 * \param nz the number of nodes in z direction (Nz)
 * \param hh the grid spacing (h)
 */
Grid::Grid(const int& nx, const int& ny, const int& nz, const double& hh)
{
  Nx = nx;
  Ny = ny;
  Nz = nz;
  Ntot = nx*ny*nz;
  h = hh;
  ind_source=-1;
  tab_handle_trial.resize(Ntot);
}

/// Class destructor
Grid::~Grid(){} 

/// method to initialise grid.
/*! \param vel a vector of size Ntot which contains the horizontal wavespeed at each node.
 * \param eps a vector of size Ntot which contains the anisotropy parameter at each node.
 * Note that no check is made that vel and eps are indeed of size Ntot. The code will crash is they aren't.
 */
void Grid::initialise(const std::vector<double>& vel, const std::vector<double>& eps)
{
  int p, si, sj, sk;
        
  for (p=0; p<Ntot; p++)
    {
      node.push_back(Node());
      node[p].index = p;
      ind2sub(p, si, sj, sk);
      node[p].i = si;
      node[p].j = sj;
      node[p].k = sk;
      node[p].nhb[0] = sub2ind(node[p].i-1,node[p].j,node[p].k);
      node[p].nhb[1] = sub2ind(node[p].i+1,node[p].j,node[p].k);
      node[p].nhb[2] = sub2ind(node[p].i,node[p].j-1,node[p].k);
      node[p].nhb[3] = sub2ind(node[p].i,node[p].j+1,node[p].k);
      node[p].nhb[4] = sub2ind(node[p].i,node[p].j,node[p].k-1);
      node[p].nhb[5] = sub2ind(node[p].i,node[p].j,node[p].k+1);
      node[p].V = vel[p];
      node[p].E = eps[p];
      if (node[p].V<=0)
	{
	  node[p].unknown = false;
	  node[p].known = true;
	}
    }
}

/// method to initialise grid using a vecotr of logV instead of V.
/*! \param lvel a vector of size Ntot which contains the log horizontal wavespeed at each node.
 * \param eps a vector of size Ntot which contains the anisotropy parameter at each node.
 * Note that no check is made that vel and eps are indeed of size Ntot. The code will crash is they aren't.
 */
void Grid::initialise_log(const std::vector<double>& lvel, const std::vector<double>& eps)
{
  int p, si, sj, sk;
        
  for (p=0; p<Ntot; p++)
    {
      node.push_back(Node());
      node[p].index = p;
      ind2sub(p, si, sj, sk);
      node[p].i = si;
      node[p].j = sj;
      node[p].k = sk;
      node[p].nhb[0] = sub2ind(node[p].i-1,node[p].j,node[p].k);
      node[p].nhb[1] = sub2ind(node[p].i+1,node[p].j,node[p].k);
      node[p].nhb[2] = sub2ind(node[p].i,node[p].j-1,node[p].k);
      node[p].nhb[3] = sub2ind(node[p].i,node[p].j+1,node[p].k);
      node[p].nhb[4] = sub2ind(node[p].i,node[p].j,node[p].k-1);
      node[p].nhb[5] = sub2ind(node[p].i,node[p].j,node[p].k+1);
      node[p].V = std::exp(lvel[p]);
      node[p].E = eps[p];
      if (node[p].V<=0)
	{
	  node[p].unknown = false;
	  node[p].known = true;
	}
    }
}

/// method to reset all nodes of the grid to +infinity, so that the grid is ready for another computation using the same V and E etc.
void Grid::reset_T()
{
  int p;

  ind_source = -1;

  for (p=0; p<Ntot; p++)
    {
      node[p].T = std::numeric_limits<double>::infinity();
      if (node[p].V>0)
	{
	  node[p].trial = node[p].known = false;
	  node[p].unknown = true;
	}
    }
}


/// method to reset all nodes of the grid to +infinity, and changes the values of V and E at all nodes, so that the grid is ready for another computation using the new V and E.
/*\param vel a vector of size Ntot containing  horizontal velocities.
 *\param eps a vector of size Ntot containing anisotropy parameters.
 */
void Grid::reset(const std::vector<double>& vel, const std::vector<double>& eps)
{
  int p;

  if (vel.size()!=node.size() || eps.size()!=node.size())
    {
      std::cout << "Cannot reset grid: wrong size of input vectors.";
      exit(1);
    }
      
  ind_source = -1;

  for (p=0; p<Ntot; p++)
    {
      node[p].T = std::numeric_limits<double>::infinity();
      node[p].V = vel[p];
      node[p].E = eps[p];
      if (node[p].V>0)
	{
	  node[p].trial = node[p].known = false;
	  node[p].unknown = true;
	}
      else
	{
	  node[p].trial = node[p].unknown = false;
	  node[p].known = true;
	}
    }
}

/// method to perform the fast marching method and compute arrival times at each node.
/*! \param ind0 index of the source node, which has an arrival time of 0.0
 * \param box true if a sourcebox with analytical solutions is to be used, false otherwise.
 * \return the number of iterations.
 */
int Grid::march(const int& ind0, const bool& box)
{
  int n, ind_current, indn, iter;
  double Tbar;

  // set source point
  ind_source = ind0;

  node[ind0].T = 0.0;
  node[ind0].known = true;
  node[ind0].unknown = false;

  // update the neighbours
  if (box)
    {
      // initialise first trial arrivals by using exact solutions for all the direct neighbour and diagonal points around the source, using the source's velocity and anisotropy parameter:
      init_box();
    }
  else
    { 
      for (n=0; n<6; n++)
	{
	  indn = node[ind0].nhb[n]; // store index of current neighbour into indn for convenience
	  if (indn>-1)
	    {
	      Tbar = update( indn , n );
		
	      node[indn].T = Tbar;
	      node[indn].unknown = false;
	      node[indn].trial = true;
	      tab_handle_trial[indn] = trial_heap.push( node[indn] );
	    }
	}
    }

  iter = 1;

  // main loop !!
  while ( trial_heap.size()>0 )
    {
      //new iteration
      iter++;

      // find min in trial _heap
      ind_current = trial_heap.top().index;
      // remove point from heap top
      trial_heap.pop();
      // now it is known:
      node[ind_current].trial = false;
      node[ind_current].known = true;

      //std::cout << "current node: (" << node[ind_current].i << "," << node[ind_current].j << "," << node[ind_current].k << ")" << "\n";
      //std::cout << "time of current node is " <<  node[ind_current].T << "\n";

      // update the neighbours 
      for (n=0; n<6; n++)
	{
	  indn = node[ind_current].nhb[n]; // store index of current neighbour into indn for convenience
	  if (indn>-1 && ( node[indn].unknown || node[indn].trial) )
	    {
		
	      Tbar = update( indn , n );
	      if (node[indn].unknown)
		{
		  node[indn].T = Tbar;
		  node[indn].unknown = false;
		  node[indn].trial = true;
		  tab_handle_trial[indn] = trial_heap.push( node[indn] );
		}

	      if ( (node[indn].trial) && (Tbar<node[indn].T) )
		{
		  node[indn].T = Tbar;
		  trial_heap.increase( tab_handle_trial[indn], node[indn] );
		}
	    }
	}
    }
  return iter;
		
}

/// method to export the arrival times into a Data structure, to save it later.
/*! \return a vector of <double> containing T of all nodes.
 */
std::vector<double> Grid::export_Tvector()
{
  std::vector<double> expt;
  int k;
    
  for (k=0; k<Ntot; k++)
    expt.push_back(node[k].T);

  return expt;
}

///method to import arrival times from a Data structure.
/*!\param d a Data object (loaded with arrival times).
 */
void Grid::import_T(const Data& d)
{
  int p;

  if ( (h!=d.h) || (Nx!=d.Nx) || (Ny!=d.Ny) || (Nz!=d.Nz) )
    {
      std::cout << "Error in Grid.import_T: Data does not match the grid size or step size\n" ;
    }
  else
    {
      std::cout << "Importing arrival time data...\n";
      for (p=0; p<d.dat.size(); p++)
	{
	  node[p].T = d.dat[p];
	  if (node[p].T==0.0)
	    {
	      ind_source = p;
	      std::cout << "Found the source node at index " << p << "\n";
	    }
	}
      std::cout << "Done !\n";
    }
}

/// method to convert subscript to linear indices.
/*! \param i subcript index (i-th node in x direction in the grid)
 * \param j subcript index (j-th node in y direction in the grid)
 * \param k subcript index (k-th node in z direction in the grid)
 * \return the linear index of the node in the grid.
 */
int Grid::sub2ind(const int& i, const int& j, const int& k)
{
  int ind;
  if ( (i>=0) && (i<Nx) && (j>=0) && (j<Ny) && (k>=0) && (k<Nz) )
    ind = i + j*Nx + k*Nx*Ny;
  else
    ind = -1;
  return ind;
}
    
/// method to convert linear to subscript indices:
/*!\param ind the linear index to be computed 
 * \param i subcript index (i-th node in x direction in the grid)
 * \param j subcript index (j-th node in y direction in the grid)
 * \param k subcript index (k-th node in z direction in the grid)
 */
void Grid::ind2sub(const int& ind, int& i, int& j, int& k)
{
  int vk;
        
  vk = int ( std::fmod(ind, Nx*Ny) ) ;
        
  k = (ind - vk)/(Nx*Ny) ;
  i = int ( std::fmod(vk, Nx) ) ;
  j =  (vk - i)/Nx ;
        
  return;
}

    


/// method to build a Tetra (tetrahedron) around a given node.
/*!\param nd the node at which T is to be computed
 * \param nhbx an integer (-1, 0 or 1) corresponding to the index in the nhb array of the x-neighbour of node nd which should be considered for the x-vertex of the tetrahedron.
 * \param nhby an integer (-1, 2 or 3) corresponding to the index in the nhb array of the y-neighbour of node nd which should be considered for the y-vertex of the tetrahedron.
 * \param nhbz an integer (-1, 4 or 5) corresponding to the index in the nhb array of the z-neighbour of node nd which should be considered for the z-vertex of the tetrahedron.
 * \return the Tetra structure with all the attributes correctly set.
 */
Tetra Grid::buildtetra(const Node& nd, const int& nhbx, const int& nhby, const int& nhbz)
{
  Tetra tt(h, nd.V, nd.E);
  int ind_1, ind_2;

  // assign "a" values, corresponding to X (corner of tetrahedron in x direction)
  ind_1 = nd.nhb[nhbx];
  if ( (ind_1>-1) && (node[ind_1].known) )
    {
      tt.test_a = true;
      tt.a = node[ind_1].T;
      ind_2 = node[ind_1].nhb[nhbx];
      if ( (ind_2>-1) && (node[ind_2].known) && (node[ind_2].T< node[ind_1].T ) )
	{
	  tt.ap = node[ind_2].T;
	  tt.sw_a = 1;
	}
    }
        
  // assign "b" values, corresponding to Y (corner of tetrahedron in y direction)
  ind_1 = nd.nhb[nhby];
  if ( (ind_1>-1) && node[ind_1].known )
    {
      tt.test_b = true;
      tt.b = node[ind_1].T;
      ind_2 = node[ind_1].nhb[nhby];
      if  ( (ind_2>-1) && (node[ind_2].known) && (node[ind_2].T< node[ind_1].T ) )
	{
	  tt.bp = node[ind_2].T;
	  tt.sw_b = 1;
	}
    }
        
  // assign "c" values, corresponding to Z (corner of tetrahedron in z direction)
  ind_1 = nd.nhb[nhbz];
  if ( (ind_1>-1) && node[ind_1].known )
    {
      tt.test_c = true;
      tt.c = node[ind_1].T;
      ind_2 = node[ind_1].nhb[nhbz];
      if  ( (ind_2>-1) && (node[ind_2].known) && (node[ind_2].T< node[ind_1].T ) )
	{
	  tt.cp = node[ind_2].T;
	  tt.sw_c = 1;
	}
    }
  return tt;
}
    
/// method to find vertices of tetrahedron given one segment of it and its direction.
/*! \param indn the index of the "central" node (node of interest) which needs to be updated.
 * \param dir an integer between 0 and 5, corresponding to an index in the nhb array of node indn, pointing to the index of the node from which the neighbour indn is updated.
 * \return an array of 4 Tetra corresponding to all the possible tetrahedra sharing the node indn and its neighbour nhb[dir].
 */
std::array<Tetra,4> Grid::findalltetra(const int& indn, const int& dir)
{
  std::array<Tetra,4> vertices;

  switch (dir)
    {
    case 0:
      {
	vertices[0] = buildtetra(node[indn], 1, 2, 4);
	vertices[1] = buildtetra(node[indn], 1, 2, 5);
	vertices[2] = buildtetra(node[indn], 1, 3, 4);
	vertices[3] = buildtetra(node[indn], 1, 3, 5);
	break;
      }
    case 1:
      {
	vertices[0] = buildtetra(node[indn], 0, 2, 4);
	vertices[1] = buildtetra(node[indn], 0, 2, 5);
	vertices[2] = buildtetra(node[indn], 0, 3, 4);
	vertices[3] = buildtetra(node[indn], 0, 3, 5);
	break;
      }
    case 2:
      {
	vertices[0] = buildtetra(node[indn], 0, 3, 4);
	vertices[1] = buildtetra(node[indn], 0, 3, 5);
	vertices[2] = buildtetra(node[indn], 1, 3, 4);
	vertices[3] = buildtetra(node[indn], 1, 3, 5);
	break;
      }
    case 3:
      {
	vertices[0] = buildtetra(node[indn], 0, 2, 4);
	vertices[1] = buildtetra(node[indn], 0, 2, 5);
	vertices[2] = buildtetra(node[indn], 1, 2, 4);
	vertices[3] = buildtetra(node[indn], 1, 2, 5);
	break;
      }
    case 4:
      {
	vertices[0] = buildtetra(node[indn], 0, 2, 5);
	vertices[1] = buildtetra(node[indn], 0, 3, 5);
	vertices[2] = buildtetra(node[indn], 1, 2, 5);
	vertices[3] = buildtetra(node[indn], 1, 3, 5);
	break;
      }
    case 5:
      {
	vertices[0] = buildtetra(node[indn], 0, 2, 4);
	vertices[1] = buildtetra(node[indn], 0, 3, 4);
	vertices[2] = buildtetra(node[indn], 1, 2, 4);
	vertices[3] = buildtetra(node[indn], 1, 3, 4);
	break;
      }
    }
  return vertices;
}
    
/// method to update a node indn from a neighbouring node just accepted.
/*! \param indn the linear index of the node to be updated.
 * \param dir an integer between 0 and 5 corresponding to an index in the array nhb (of node indn), pointing to the linear index of the node which has just been accepted and from which the node indn should be updated.
 * \return the arrival time T at node indn
 */
double Grid::update(const int& indn, const int& dir)
{
  std::array<Tetra,4> vert;
  int k;
  double Tbar = std::numeric_limits<double>::infinity();

  if (node[indn].unknown || node[indn].trial)
    {
      vert = findalltetra(indn, dir);

      //std::cout << "E = " << vert[0].eps << "\n";
      //std::cout << "V = " << vert[0].vel << "\n";
	
      //std::cout << "(" << node[indn].i << "," << node[indn].j << "," << node[indn].k << ")" << "\n";
	
      for(k=0; k<4; k++)
	{
	  Tbar = std::min(Tbar, vert[k].compute_T());

	}
      //std::cout << "T = " << Tbar << "\n";
      //if ( (node[indn].i==30) && (node[indn].j==31) )
      //  exit(1);

    }
    
  return Tbar;
}

/// method to initialise the arrival times around the source node.
/*! \param ind0 the linear index of the source node.
 * This method initialises the arrival times around the source point dy using analytical solutions (to avoid large 1st order errors propagating everywhere), assuming a constant velocity and anisotropy around the source point. It basically updates all the points in a box of size 2h*2h*2h around the source. This trick is due to Rickett and Fomel, SEP Report 2000.
 */
void Grid::init_box()
{
  int i0, j0, k0, ind0, indn_box;
  int ai, aj, ak;
  double v, eps, v45;
  const double PI=3.14159265;


  ind0 = ind_source;

  // short hand notations:
  v = node[ind0].V;
  eps = node[ind0].E;
  
  // compute group velocity at 45degrees, using Velocity class:
  Velocity vel(v,eps);
  v45 = vel.Vgroup(PI/4.0);

  // get the subscript indices of the source point
  ind2sub(ind0, i0, j0, k0);

  // browse through all the points around the source, including diagonals:
  for (ai=-1; ai<2; ai++)
    {
      for (aj=-1; aj<2; aj++)
	{
	  for (ak=-1; ak<2; ak++)
	    {
	      // shorthand notation for the point index:
	      indn_box = sub2ind(i0+ai, j0+aj, k0+ak);
	      if (indn_box>-1) // if the point exist in the grid,
		{
		  if (ak!=0) // if there is a vertical component
		    {
		      if ( (ai==0) && (aj==0) )
			node[indn_box].T = h/(v*(1+eps)); // this is the direct point
		      else if ( ( (ai!=0) && (aj==0) ) || ( (ai==0) && (aj!=0) ) )
			node[indn_box].T = (h/v45)*std::sqrt(2); // this is a 2d diagonal point
		      else 
			node[indn_box].T = (h/v45)*std::sqrt(3); // this is a 3d diagonal point

		      // now update the status of the point and puch it the heap
		      node[indn_box].unknown = false;
		      node[indn_box].trial = true;
		      tab_handle_trial[indn_box] = trial_heap.push( node[indn_box] );
		    }
		  else if ( (ai!=0) || (aj!=0) ) // if there is no vertical component
		    {
		      if ( (ai==0) || (aj==0) )
			node[indn_box].T = h/v; // this is a direct neighbour
		      else
			node[indn_box].T = (h/v)*std::sqrt(2); // this is a 2d diagonal point in the plane of isotropy

		      // and update the status and push the point
		      node[indn_box].unknown = false;
		      node[indn_box].trial = true;
		      tab_handle_trial[indn_box] = trial_heap.push( node[indn_box] );
		    }
		}
	    }
	}
    }
}
