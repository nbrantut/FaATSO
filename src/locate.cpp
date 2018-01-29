#include<iostream>
#include<cstdlib>
#include<cmath>
#include<cfloat>
#include<armadillo>

# include "fmmlib.hpp"
# include "fmmiolib.hpp"

int main(int nargin, char *varargin[] )
{
  // read input file
  std::ifstream fis;
  std::string filename, line;
  std::vector<std::string> list;
  
  filename.assign(varargin[1]);

  fis.open(filename.c_str());
  while(!fis.eof())
    {
      std::getline(fis, line);
      std::cout << line << "\n";
          
      list.push_back(line);
     
    }  
  
  // parse input file
  std::string sensorfile, tobsfile, vfile, Efile, islog, locfile;
  int n;
  double sD;
  
  sensorfile = list[0];
  tobsfile = list[1];
  vfile = list[2];
  Efile = list[3];
  islog = list[4];
  n = std::stoi(list[5]);
  sD = std::stod(list[6]);
  locfile = list[7];

  
  // read velocity and anisotropy
  Data v0, E0;

  v0.load(vfile.c_str());
  E0.load(Efile.c_str());
  // check grid size consistency between v and E files:
  if ( (v0.Nx!=E0.Nx) || (v0.Ny!=E0.Ny) || (v0.Nz!=E0.Nz) || (v0.h!=E0.h) )
    {
      std::cout << "The two input files seem to correspond to different grids !\n";
      std::cout << "Grid size in " << vfile << " is (" << v0.Nx << "," << v0.Ny << "," << v0.Nz << "), spacing is " << v0.h << "\n";
      std::cout << "Grid size in " << Efile << " is (" << E0.Nx << "," << E0.Ny << "," << E0.Nz << "), spacing is " << E0.h << "\n";
      
      exit(1);
    }
  else
    {
      std::cout << "Grid size is (" << v0.Nx << "," << v0.Ny << "," << v0.Nz << "), spacing is " << v0.h << "\n";
    }
  
  // refine grid
  Data vf, Ef;
  
  std::cout << "Refine grids...\n";
  vf = v0.refine(n);
  Ef = E0.refine(n);

  // initialise grid
  Grid grid(vf.Nx,
	    vf.Ny,
	    vf.Nz,
	    vf.h);

  if (islog.compare("log") == 0)
    grid.initialise_log(vf.dat, Ef.dat);
  else
    grid.initialise(vf.dat, Ef.dat);

  // read sensor positions
  std::vector<locations::Location> sensors;
  
  std::cout << "Read locations of sensors from file " << sensorfile << " ...\n";
  sensors = locations::readlocations(sensorfile.c_str());

  std::cout << "Found " << sensors.size() << " sensors.\n";

  // find sensor indices within the grid:
  for (int k=0; k<sensors.size(); k++)
    {
      sensors[k].findind(vf.Nx, vf.Ny, vf.Nz, vf.h);
      std::cout << "Velocity at sensor " << k << " is: " << vf.dat[sensors[k].ind] << "\n";
    }


  // read arrival times
  arma::mat tobs;
   
  std::cout << "Read arrival times...\n";
  tobs.load(tobsfile, arma::raw_ascii);


  
  // compute arrival times from each sensor:
  arma::mat Tcalc(grid.Ntot, sensors.size());
  
  for (int k=0; k<sensors.size(); k++)
    {
      // reset times and march grid from the current virtual source
      grid.reset_T();
      std::cout << "March grid from virtual source " << k << " ...\n";
      grid.march(sensors[k].ind, true);

      // save arrival times into the big Tcalc matrix
      std::cout << "Store computed arrivals ...\n";
      Tcalc.col(k) = arma::colvec( grid.export_Tvector() );
    }
  
  // determine the locations using least L1 norm
  arma::mat loc(tobs.n_rows, 4);
  std::cout << "\nLocating events...\n";
  
  for (int i=0; i<tobs.n_rows; i++)
    {

      std::cout << "Event " << i << " : ";
      
      // indices of nonzero picks
      arma::rowvec tobsi = tobs.row(i);
      arma::uvec ind = arma::find(tobsi>0);
      int jmax = 0;
      double L;
      double Lmax = 0.0;
      double torigin;
      
      for (int j=0; j<grid.Ntot; j++)
	{
	  arma::rowvec Tcalcn = Tcalc.row(j);
	  // origin time of event
	  double to_trial = arma::median( tobsi(ind) - Tcalcn(ind));

	  // likelihood
	  if (std::abs(to_trial)>1/DBL_EPSILON)
	    L = 0.0;
	  else
	    L = std::exp( -arma::sum( arma::abs( tobsi(ind) - to_trial - Tcalcn(ind) )/sD ) );

	  // check if minimum
	  if (L > Lmax)
	    {
	      jmax = j;
	      Lmax = L;
	      torigin = to_trial;
	    }
	}
      // now I have the index of the best location in the grid
      loc(i, 0) = grid.node[jmax].i*grid.h;
      loc(i, 1) = grid.node[jmax].j*grid.h;
      loc(i, 2) = grid.node[jmax].k*grid.h;
      loc(i, 3) = torigin;

      std::cout << "(" << loc(i,0) << "," << loc(i,1) << "," << loc(i,2) << ";" << loc(i,3) << ")\n";      
    }

  std::cout <<"Done!\n";
  
  // export that stuff to the appropriate file
  std::cout <<"Saving file...\n";
  loc.save(locfile, arma::raw_ascii);

  return 0;
  
}
