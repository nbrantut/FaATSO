#include<cstdio>
#include<cstdlib>
#include<iostream>
#include<vector>
#include<armadillo>

#include "fmmiolib.hpp"
#include "fmmlib.hpp"

using namespace std;

int main(int nargin, char **varargin)
{
  vector<locations::Location> sensors;
  vector<locations::Location> events;
  arma::mat t_shots, t_events;
  arma::uvec I_shots, I_events;
  Data tab_V, tab_E, Vf, Ef;
  char *Vfile, *Efile, *sensfile, *evtfile, *tshotfile, *tevtfile;
  double noise;

  cout << "Running program " << varargin[0] << "\n";
  // assign command line parameters (file names):
  noise = atof(varargin[1]);
  Vfile = varargin[2];
  Efile = varargin[3];
  sensfile = varargin[4];
  evtfile = varargin[5];
  tshotfile = varargin[6];
  tevtfile = varargin[7];

  // say a few things about input files, to check:
  cout << "File with velocities is......... " << Vfile << "\n";
  cout << "File with anisotropy is......... " << Efile << "\n";
  cout << "File with sensor locations is... " << sensfile << "\n";
  cout << "File with event locations is.... " << evtfile << "\n";
  cout << "File with t_shots is............ " << tshotfile << "\n";
  cout << "File with t_events is........... " << tevtfile << "\n";

  // load data from files:
  tab_V.load(Vfile);
  tab_E.load(Efile);
  Vf = tab_V.refine(10);
  Ef = tab_E.refine(10);
  Grid grid(Vf.Nx,Vf.Ny,Vf.Nz,Vf.h);
  grid.initialise_log(Vf.dat,Ef.dat);
  
  sensors = locations::readlocations(sensfile);
  events = locations::readlocations(evtfile);

  cout << "Read " << sensors.size() << " sensors.\n";
  cout << "Read " << events.size() << " events.\n";
 
  // find indices of sensors and events in the grid
  I_shots.set_size(sensors.size());
  for (int k=0; k<sensors.size(); k++)
    {
      sensors[k].findind(grid);
      I_shots(k) = sensors[k].ind;
      //cout << grid.node[I_shots(k)].V << "\n";
    }

  //I_shots.print();

  I_events.set_size(events.size());
  for (int k=0; k<events.size(); k++)
    {
      events[k].findind(grid);
      I_events(k) = events[k].ind;
    }

  t_shots.set_size(sensors.size(), sensors.size());
  t_events.set_size(events.size(), sensors.size());
  
  // compute times
  cout << "Compute forward model...\n";
  for (int k=0; k<sensors.size(); k++)
    {
      cout << "March grid using sensor " << k << " as source...\n";
      grid.march(sensors[k].ind, true);

      // compute all the interesting times
      cout << "Extract times at locations of interest...\n";

      arma::vec tmp(grid.export_Tvector());
      
      t_shots.col(k) = tmp(I_shots);
      t_events.col(k) = tmp(I_events);

      grid.reset_T();
    }

  //correct for t0 for each event:
  for (int k=0; k<events.size(); k++)
    t_events.row(k) += events[k].t0;

    // add noise
  arma::arma_rng::set_seed_random();
  t_events += noise
    *(arma::randu<arma::mat>(t_events.n_rows,t_events.n_cols) - 0.5);

  cout << "Write times in text files...\n";
  t_shots.save(tshotfile, arma::raw_ascii);
  t_events.save(tevtfile, arma::raw_ascii);

  cout << "End !\n";
  return 0;
}
