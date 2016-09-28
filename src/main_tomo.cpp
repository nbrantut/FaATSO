# include "tomolib.hpp"

using namespace std;

int main(int nargin, char *varargin[])
{
  Parameters param;
  string paramfilename;
  M model;
  D data;
  QuasiNewton qn;

  if (nargin!=2)
    {
      cout << "Wrong number of input arguments ! Needs a parameter file as input.\n";
	exit(1);
    }
  
  cout << "Running " << varargin[0] << "...\n";
  cout << "Using " << varargin[1] << " as input parameter file.\n";

  paramfilename.assign(varargin[1]);
  
  param.initialise(paramfilename);

  model.initialise(param);
  data.initialise(param, model);
  qn.initialise(param);

  qn.run(model, data);


  return 0;
}
