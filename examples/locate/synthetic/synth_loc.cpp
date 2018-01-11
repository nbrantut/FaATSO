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

  std::cout << "Starting " << varargin[0] << "...\n\n";
  
  filename.assign(varargin[1]);

  std::cout << "Reading input file "  <<  varargin[1] << " ...\n";
  fis.open(filename.c_str());
  while(!fis.eof())
    {
      std::getline(fis, line);
      std::cout << line << "\n";
      
      list.push_back(line);
      
    }


  std::cout << "\nParsing the input file...\n";
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
  std::cout << "Reading velocity and aniostropy...\n";

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

  // read events positions
  std::vector<locations::Location> events;

  std::cout << "Read locations of events from file " << locfile << " ...\n";
  events = locations::readlocations(locfile.c_str());

  std::cout << "Found " << events.size() << " events.\n";

  // find event indices within the grid:
  for (int k=0; k<events.size(); k++)
    {
      events[k].findind(vf.Nx, vf.Ny, vf.Nz, vf.h);
    }
  
  // compute arrival times from each sensor:
  arma::mat tobs(events.size(), sensors.size());
  
  for (int k=0; k<sensors.size(); k++)
    {
      // reset times and march grid from the current virtual source
      grid.reset_T();
      std::cout << "March grid from virtual source " << k << " ...\n";
      grid.march(sensors[k].ind, true);

      // save arrival times
      std::cout << "Store arrivals ...\n";
      for (int i=0; i<events.size(); i++)
	{
	  tobs(i,k) = grid.node[events[i].ind].T + events[i].t0;
	}
    }

  // add noise
  tobs = tobs + sD*(2.0*arma::randu(arma::size(tobs)) - 1.0 );

  std::cout <<"Done!\n";
  
  // export that stuff to the appropriate file
  std::cout <<"Saving file...\n";
  tobs.save(tobsfile, arma::raw_ascii);

  return 0;
  
}
