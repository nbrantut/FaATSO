# include <iostream>
# include <cstdlib>

# include "fmmlib.hpp"
# include "fmmiolib.hpp"
# include "raylib.hpp"

using namespace std;


int main(int nargin, char *varargin[] )
{
  int ind0, n;
  clock_t tim = clock();
  Data tab_V, tab_E, tab_T;
  bool box=false;


  // use nargin and varargin to assign parameters from read files
  if (nargin<6)
    {
      cout << "Wrong number of input parameters: should be [(char) Vfilename] [(char) Efilename] [(int) Sourceindex] [(1/0) sourcebox] [(char) Tfilename]\n";
      exit(1);
    }
  else
    {
      cout << "Calling " << varargin[0] << " ...\n";
      cout << "Reading file " << varargin[1] << " for V...\n";

      tab_V.load(varargin[1]);

      cout << "Reading file " << varargin[2] << " for E...\n";
      
      tab_E.load(varargin[2]);

      ind0 = atoi( varargin[3] );

      if ( atoi( varargin[4] )!=0)
	box = true;

      if ( (tab_V.Nx!=tab_E.Nx) || (tab_V.Ny!=tab_E.Ny) || (tab_V.Nz!=tab_E.Nz) || (tab_V.h!=tab_E.h) )
	{
	  cout << "The two input files seem to correspond to different grids !\n";
	  cout << "Grid size for Vfilename is (" << tab_V.Nx << "," << tab_V.Ny << "," << tab_V.Nz << "), spacing is " << tab_V.h << "\n";
	  cout << "Grid size for Vfilename is (" << tab_E.Nx << "," << tab_E.Ny << "," << tab_E.Nz << "), spacing is " << tab_E.h << "\n";

	  exit(1);
	}
      else
	{
	  cout << "Grid size is (" << tab_V.Nx << "," << tab_V.Ny << "," << tab_V.Nz << "), spacing is " << tab_V.h << "\n";
	  cout << "Source is at index : "<< ind0 << "\n";
	}
    }


  // make the grid
  cout << "Make the grid...\n";
  Grid grid(tab_V.Nx,tab_V.Ny,tab_V.Nz,tab_V.h);

  // initialise the nodes of the grid using the arrays V and E
  grid.initialise(tab_V.dat,tab_E.dat);
 
  // set source point
  if ( (ind0<0) || (ind0>=grid.Ntot) ) // check that the source is inside the grid!
    {
      cout << "Source point is not in the grid ! \n";
      cout << "Number of grid points is " << grid.Ntot << "\n";
      exit(1);
    }

  cout << "Compute times...\n";

  n = grid.march(ind0,box);

  cout << "Completed the march in " << n << " iterations.\n";
  cout << "Export file with arrival times...\n" ;
  
  tab_T.dat = grid.export_Tvector();
  tab_T.Nx = grid.Nx;
  tab_T.Ny = grid.Ny;
  tab_T.Nz = grid.Nz;
  tab_T.h = grid.h;

  tab_T.save(varargin[5]);

  tim = clock() - tim;

  cout << "time elapsed: " << double(tim)/CLOCKS_PER_SEC << " sec\n";

  return 0;
}
