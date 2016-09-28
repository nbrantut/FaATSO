#include<iostream>
#include<cstdlib>
#include<cmath>

# include "tomolib.hpp"
# include "fmmlib.hpp"
# include "raylib.hpp"

/// basic class constructor. Sets default values. These are quite arbitrary so it is not very useful at this stage (was written as a test).
Parameters::Parameters()
{
  folder = "";
  paramfile = "none";
  sensorfile = "sensors.txt";
  eventfile = "events.txt";
  tsurveyfile = "tsurveys.txt";
  teventfile = "tevents.txt";
  vpriorfile = "v0.dat";
  Epriorfile = "E0.dat";
  vpostfile = "vpost.dat";
  Epostfile = "Epost.dat";
  eventpostfile = "evtpost.txt";
  dcalcfile = "dcalc.txt";
  Cmpostfile = "CMpost.txt";
  mpostfile = "mpost.txt";
  residualfile = "residuals.txt";
  n = 10;
  s_surveys = 0.1;
  s_events = 0.5;
  maxit = 12;
  tol = 1e-6;
  mu = 1.0;
  s_v = 2.0;
  s_E = 0.1;
  lambda = 100.0;
  s_x = 10.0;
  s_y = 10.0;
  s_z = 10.0;
  s_t0 = 10.0;
}

/// class constructor with call to a parameter file.
Parameters::Parameters(const std::string& filename)
{
  initialise(filename);
}

/// initialisation of parameters from an ascii parameter file.
void Parameters::initialise(const std::string& filename)
{
  std::ifstream fis;
  std::string line;
  std::vector<std::string> list;
  std::size_t found;
  
  paramfile = filename;

  fis.open(filename.c_str());
  while(!fis.eof())
    {
      std::getline(fis, line);
      found = line.find_first_of("\t\n");
      if (found<std::string::npos)
	{
	  list.push_back(line.substr(0,found));
	}
    }
  
  // at this point we have a vector of strings containing all parameter values
  // now extract and assign to members of PArameters class:
  if (list.size()!=27)
    {
      std::cout << "Invalid parameter file ! Needs to be 27 lines." << "\n";
      std::cout << "The file " << filename << " contains " <<list.size() << " lines.\n";
      std::cout << "Last line contains: "<< list.back() << " \n";
      exit(EXIT_FAILURE);
    }
  else
    {
      folder = list[0];
      sensorfile = folder + list[1];
      eventfile = folder + list[2];
      tsurveyfile = folder + list[3];
      teventfile = folder + list[4];
      vpriorfile = folder + list[5];
      Epriorfile = folder + list[6];
      vpostfile = folder + list[7];
      Epostfile = folder + list[8];
      eventpostfile = folder + list[9];
      dcalcfile = folder + list[10];
      Cmpostfile = folder + list[11];
      mpostfile = folder + list[12];
      residualfile = folder + list[13];
      n = std::stoi(list[14]);
      s_surveys = std::stod(list[15]);
      s_events = std::stod(list[16]);
      maxit = std::stoi(list[17]);
      tol = std::stod(list[18]);
      mu = std::stod(list[19]);
      s_v = std::stod(list[20]);
      s_E = std::stod(list[21]);
      lambda = std::stod(list[22]);
      s_x = std::stod(list[23]);
      s_y = std::stod(list[24]);
      s_z = std::stod(list[25]);
      s_t0 = std::stod(list[26]);
    }
  
}

/// Default constructor
M::M(){}

/// constructor with call to a parameters set
M::M(const Parameters& param)
{
  initialise(param);
}

/// Initialisation of M from parameters.
/*!\param param the structure containing input parameters.
 */
void M::initialise(const Parameters& param)
{
  arma::vec mv, mE;
  arma::vec tmp(4);

  // read filenames
  vpostfile = param.vpostfile;
  Epostfile = param.Epostfile;
  eventpostfile = param.eventpostfile;
  mpostfile = param.mpostfile;
  
  // read locations of events
  std::cout << "Read locations of events from file " << param.eventfile << " ...\n";
  events = locations::readlocations(param.eventfile.c_str());

  std::cout << "Found " << events.size() << " events.\n";

  // read velocity and anisotropy
  std::cout << "Read prior velocity and anisotropy from files " << param.vpriorfile << " and " << param.Epriorfile << " ...\n";
  v0.load(param.vpriorfile.c_str());
  E0.load(param.Epriorfile.c_str());

  // check grid size consistency between v and E files:
  if ( (v0.Nx!=E0.Nx) || (v0.Ny!=E0.Ny) || (v0.Nz!=E0.Nz) || (v0.h!=E0.h) )
    {
      std::cout << "The two input files seem to correspond to different grids !\n";
      std::cout << "Grid size in " << param.vpriorfile << " is (" << v0.Nx << "," << v0.Ny << "," << v0.Nz << "), spacing is " << v0.h << "\n";
      std::cout << "Grid size in " << param.Epriorfile << " is (" << E0.Nx << "," << E0.Ny << "," << E0.Nz << "), spacing is " << E0.h << "\n";
      
      exit(1);
    }
  else
    {
      std::cout << "Grid size is (" << v0.Nx << "," << v0.Ny << "," << v0.Nz << "), spacing is " << v0.h << "\n";
    }

  // covariance parameters
  std::cout << "Read covariance parameters:\n";
  
  lambda = param.lambda;
  s_v = param.s_v;
  s_E = param.s_E;
  s_x = param.s_x;
  s_y = param.s_y;
  s_z = param.s_z;
  s_t0 = param.s_t0;
  
  std::cout << "lambda = " << lambda << "\n";
  std::cout << "sigma_v = " << s_v << "\n";
  std::cout << "sigma_E = " << s_E << "\n";
  std::cout << "sigma_x = " << s_x << "\n";
  std::cout << "sigma_y = " << s_y << "\n";
  std::cout << "sigma_z = " << s_z << "\n";
  std::cout << "sigma_t0 = " << s_t0 << "\n";
  
  // now convert to vector m: start with velocity and anisotropy
  std::cout << "Start forming vector of model parameters...\n";
  mv = arma::vec(v0.dat);
  mE = arma::vec(E0.dat);
  m_prior = arma::join_cols(mv,mE); // make preliminary m_prior vector

  // find infinite values:
  ind_m_null = arma::find_nonfinite(m_prior);
  // store inifinity variable
  if (!ind_m_null.is_empty())
    {
      std::cout << "There are " << ind_m_null.size() << " infinte values.\n";
      
      // assign zeros to avoid having non finite values when doing linear algebra later on.
      m_prior.elem(ind_m_null).fill(V_SHADOW);
    }

  // replace elements of v0 and E0 with non infinite values:
  // extract good old pointer to vector m:
  int Ntot = v0.Nx*v0.Ny*v0.Nz;
  const double *pr = m_prior.memptr();

  v0.dat.assign(pr,pr+Ntot);
  E0.dat.assign(pr+Ntot,pr+2*Ntot);

  std::cout << "Refine grids...\n";
  n = param.n;
  vf = v0.refine(n);
  Ef = E0.refine(n);

  // then add event positions
  std::cout << "Determine events positions in grid, and append to model vector...\n";
  for (int k=0; k<events.size(); k++)
    {
      // determine event positions in the refined grid:
      events[k].findind(vf.Nx, vf.Ny, vf.Nz, vf.h);

      // extract (linear) vector of positions
      tmp(0) =  events[k].x;
      tmp(1) =  events[k].y;
      tmp(2) =  events[k].z;
      tmp(3) =  events[k].t0;

      // append vector to m_prior
      m_prior.insert_rows(m_prior.size(), tmp );
    }
  std::cout << "m_prior is formed !\n";
  
  // initialise m equal to m_prior
  m = m_prior;

  // build covariance matrix:
  std::cout << "Build (inverse of) covariance matrix...\n";
  build_inv_covariance();
  
}

/// method to compute covariance matrix.
void M::build_inv_covariance()
{
  int Nx = v0.Nx;
  int Ny = v0.Ny;
  int Nz = v0.Nz;
  double h = v0.h;
  int N = Nx*Ny*Nz;
  int i1, j1, k1, i2, j2, k2, I, J, tmp;
  arma::mat CvE(N,N,arma::fill::zeros);
  arma::mat CvEi, Csourcei;
  arma::vec tmp_vec(4);

  // initialise full CM^-1
  Cmi.zeros(m_prior.size(), m_prior.size());

  // compute part of CM^-1 corresponding to ln(v) and E:
  double hl = h/lambda;
  
  for (I=0; I<N; I++)
    {
      for (J=I; J<N; J++)
	{
	  tmp = int ( fmod(I,Nx*Ny));
	  k1 = (I-tmp)/(Nx*Ny);
	  i1 = int ( fmod(tmp, Nx));
	  j1 = (tmp -i1)/Nx;

	  tmp = int ( fmod(J,Nx*Ny));
	  k2 = (J-tmp)/(Nx*Ny);
	  i2 = int ( fmod(tmp, Nx));
	  j2 = (tmp -i2)/Nx;
	  
	  CvE(I,J) = exp(-hl*sqrt((i1-i2)*(i1-i2)
				 + (j1-j2)*(j1-j2)
				 + (k1-k2)*(k1-k2))
			);
	}
    }
  // assign temporary matrix:
  CvEi = arma::inv( arma::symmatu(CvE) );
  
  // assign the inverse covariance on v and E:
  Cmi.submat(0, 0, N-1, N-1) = (1/(s_v*s_v))*CvEi;
  Cmi.submat(N, N, N+N-1, N+N-1) = (1/(s_E*s_E))*CvEi;

  //compute part of CM^-1 corresponding to source locations:
  tmp_vec(0) = 1.0/(s_x*s_x);
  tmp_vec(1) = 1.0/(s_y*s_y);
  tmp_vec(2) = 1.0/(s_z*s_z);
  tmp_vec(3) = 1.0/(s_t0*s_t0);
  Csourcei = arma::diagmat( arma::repmat(tmp_vec,events.size(),1 ));
  // assign the inverse covariance on source parameters:
  Cmi.submat(N+N, N+N, N+N+4*events.size()-1, N+N+4*events.size()-1) = Csourcei;
}

/// method to convert m values (from the vector m) into data structures v0, E0, events. Important note: the values of m where m_prior should be -inf (i.e., zero speed, not computed) are imposed equal to a very negative number during export to v0 (and hence also for the refined grid vf). So when computing arrivla times in the next iterations, the rays will be prevented from going through V=0 regions.
void M::update_internals()
{
  int Ntot = v0.Nx*v0.Ny*v0.Nz;

  // define temporary m vector which will be masked:
  arma::vec m_tmp = m;

  // perform masking: enforce near-zero speed where needed
  m_tmp.elem(ind_m_null).fill(V_SHADOW);

  // extract good old pointer to vector m:
  const double *pr = m_tmp.memptr();

  v0.dat.assign(pr,pr+Ntot);
  E0.dat.assign(pr+Ntot,pr+2*Ntot);

  vf = v0.refine(n);
  Ef = E0.refine(n);

  for (int k=0; k<events.size(); k++)
    {
      events[k].x = m_tmp(2*Ntot+4*k);
      events[k].y = m_tmp(2*Ntot+4*k+1);
      events[k].z = m_tmp(2*Ntot+4*k+2);
      events[k].t0 = m_tmp(2*Ntot+4*k+3);
      events[k].findind(vf.Nx, vf.Ny, vf.Nz, vf.h);
    }
}

/// method to export velocity model and event locations into readable files.
/*!\param param the structre containing input parameters and file names.
 */
void M::export_data()
{
  v0.save(vpostfile.c_str());
  E0.save(Epostfile.c_str());
  writelocations(eventpostfile.c_str(), events);
}

/// method to export directly the vector m into a dedicated text file.
void M::export_m()
{
  m.save(mpostfile, arma::raw_ascii);
}

/// Default constructor
D::D(){}

/// constructor with call to a parameters set
D::D(const Parameters& param, const M& model)
{
  initialise(param, model);
}

/// Initialisation of D from parameters.
/*!\param param the object containing input parameters and file names.
 */
void D::initialise(const Parameters& param, const M& model)
{
   
    
  // read filenames
  dcalcfile = param.dcalcfile;

  
  
  std::cout << "Read locations of sensors from file " << param.sensorfile << " ...\n";
  sensors = locations::readlocations(param.sensorfile.c_str());

  std::cout << "Found " << sensors.size() << " sensors.\n";

  // find sensor indices within the grid:
  for (int k=0; k<sensors.size(); k++)
    {
      sensors[k].findind(model.vf.Nx, model.vf.Ny, model.vf.Nz, model.vf.h);
      std::cout << "ln(Velocity) at sensor " << k << " is: " << model.vf.dat[sensors[k].ind] << "\n";
    }

   // covariance parameters
  std::cout << "Read covariance parameters:\n";
  
  s_surveys = param.s_surveys;
  s_events = param.s_events;
  
  std::cout << "sigma_surveys = " << s_surveys << "\n";
  std::cout << "sigma_events = " << s_events << "\n";


  
  std::cout << "Make data vector...\n";
  build_data_vector(param.tsurveyfile,
		    param.teventfile,
		    sensors.size(),
		    model.events.size()
		    );
  

  std::cout << "Build data covariance matrix...\n";
  build_inv_covariance();
  
  std::cout << "CD^-1 done !\n";
}

/// This method builds the data vector based on input text files. It excludes arrival times taht are equal or less than zero. It also builds the field "path" containing, for each sensor, the list of indices of valid (i.e., tobs>0) shots and events.
/*!\param tsurveyfile the name of the file containing survey arrival times.
 *\param teventfile the name of the file containing events arrival times.
 *\param Nsensors the number of sensors.
 *\param Nevents the number of events.
 */
void D::build_data_vector(std::string tsurveyfile, std::string teventfile, int Nsensors, int Nevents)
{

   arma::mat t_shots, t_events;
   
    // arrival times
  std::cout << "Read arrival times...\n";
  t_shots.load(tsurveyfile, arma::raw_ascii);
  t_events.load(teventfile, arma::raw_ascii);
  
  // check that there is the correct number of arrival times (one for each source and sensor: Nsrc x Nrcv):
  if (t_shots.size() != Nsensors*Nsensors)
    {
      std::cout << "Wrong number of arrival times in file " << tsurveyfile << "\n";
      std::cout << "Found " << t_shots.size() << " whereas there should be " << Nsensors*Nsensors << "\n";
      exit(1);
    }
  if (t_events.size() != Nsensors*Nevents )
    {
      std::cout << "Wrong number of arrival times in file " << teventfile << "\n";
      std::cout << "Found " << t_events.size() << " whereas there should be " << Nsensors*Nevents << "\n";
      exit(1);
    }
  
  // only take the strictly positive arrival times: t<=0 means that the arrival time is not picked or not good. (That is my convention).

  // intiialise path vector of Sources:
  path.resize(Nsensors);
  // initialise d_obs with mx. number of elements:
  d_obs.set_size(Nsensors*(Nsensors + Nevents));
  // intialise counters:
  int cd=0; // that is for the d_obs vector
  int cp=0; // that is for the path vectors

  Ndshots = 0;
  Ndevents = 0;
  
  for (int k=0; k<Nsensors; k++)
    {
      path[k].active.set_size(Nsensors);
      cp=0;
      for (int p=0; p<Nsensors; p++)
	{
	  if (t_shots(p,k)>0.0)
	    {
	      d_obs(cd) = t_shots(p,k);
	      path[k].active(cp) = p;
	      cd++;
	      cp++;
	    }
	}
      // now resize path[k].active to eliminate zeros:
      path[k].active.set_size(cp);
      // update the number of good shots
      Ndshots += cp;
      
      path[k].passive.set_size(Nevents);
      cp=0;
      for (int p=0; p<Nevents; p++)
	{
	  if (t_events(p,k)>0.0)
	    {
	      d_obs(cd) = t_events(p,k);
	      path[k].passive(cp) = p;
	      cd++;
	      cp++;
	    }
	}
      // now resize path[k].passive to eliminate zeros:
      path[k].passive.set_size(cp);
      // update the number of good events
      Ndevents += cp;
    }
  // finally resize d_obs to elminate pending zeros:
  d_obs.set_size(cd);

  std::cout << "d_obs is of size " << arma::size(d_obs) << "\n";
    //  std::cout << "d_obs done !\n";

  // initialise d with correct size:
  d.set_size(d_obs.size());
}

void D::build_inv_covariance()
{
  // initialise with correct dimensions
  Cdi.zeros(d_obs.size(),d_obs.size());

  int c=0;

  for (int k=0; k<path.size(); k++)
    {
      for (int p=0; p<path[k].active.size(); p++)
	{
	  Cdi(c,c) = 1.0/(s_surveys*s_surveys);
	  c++;
	}
      for (int p=0; p<path[k].passive.size(); p++)
	{
	  Cdi(c,c) = 1.0/(s_events*s_events);
	  c++;
	}
    }
}
			     

/// method to export data (observed and computed) into files.
/*!\param param the Parameters structures which includes files names etc..
 */
void D::export_all()
{
  arma::mat d_th(d_obs.n_rows,4);

  int c=0;

  for (int k=0; k<path.size(); k++)
    {
      for (int p=0; p<path[k].active.size(); p++)
	{
	  d_th(c,0) = k;
	  d_th(c,1) = path[k].active(p);
	  d_th(c,2) = d_obs(c);
	  d_th(c,3) = d(c);
	  c++;
	}
      for (int p=0; p<path[k].passive.size(); p++)
	{
	  d_th(c,0) = k;
	  d_th(c,1) = path[k].passive(p);
	  d_th(c,2) = d_obs(c);
	  d_th(c,3) = d(c);
	  c++;
	}
    }

  d_th.save(dcalcfile, arma::raw_ascii);
}

/// method to compute theoretical arrival times from a model.
/*\param model an object M containing velocity, anisotropy fields and event locations.
 *\param G a sparse matrix of size (nb of data x nb of model parameters) where the jacobian of the forward model will be stored.
 */
void D::forward_model(M& model, arma::mat & G)
{
  // declarations
  Ray ray;
  Data der_v, der_E;
  //  clock_t tim;
  int Nm = model.m.size();
  int Nd = d.size();
  int Nv = model.v0.dat.size();
  int NE = model.E0.dat.size();
  int NvE = Nv+NE;

  // // // prepare matrix containing locations of (potentially) non zero elements in sparse matrix G:
  // // int Nnz = Ndshots*NvE + Ndevents*(NvE + 4);
  // arma::umat Gloc(2, Nnz);
  // arma::vec Gval(Nnz);

  // set zeros everywhere in G to start:
  G.zeros(Nd, Nm);

  // initialise grid
  Grid grid(model.vf.Nx,
	    model.vf.Ny,
	    model.vf.Nz,
	    model.vf.h);
  
  grid.initialise_log(model.vf.dat, model.Ef.dat);

  // start a big loop, for each sensor (i.e., each virtual source):
  //    - allocate grid
  //    - march grid
  //    - trace rays and compute along-ray derivatives
  //    - fill d and G line by line

  int i_data = 0;
  int i_s;
  
  for (int k=0; k<path.size(); k++)
    {
      // reset times and march grid from the current virtual source
      grid.reset_T();
      std::cout << "March grid from virtual source " << k << " ...\n";
      grid.march(sensors[k].ind, true);

      std::cout << "Trace rays to sensor " << k << " ...\n";

      for (int p=0; p<path[k].active.size(); p++)
	{
	  // find index of real source (vritual receiver)
	  i_s = path[k].active(p);
	  std::cout << "   ... from active source " << i_s <<" | ";
	  
	  // extract arrival time
	  d(i_data) = grid.node[sensors[i_s].ind].T + sensors[i_s].t0;
	  std::cout << "Dt(obs-calc) = " << d_obs(i_data) - d(i_data) << "\n";
	  // trace ray
	  ray.trace(grid, sensors[i_s].x, sensors[i_s].y, sensors[i_s].z);

	  // interpolate derivative on the coarse grid
	  ray.compute_interp_der(model.n, der_v, der_E);

	  // assign derivatives to elements of G:
	  //tim = clock();
	  
	  G(i_data, arma::span(0, der_v.dat.size()-1)) = arma::rowvec(der_v.dat);
	  G(i_data, arma::span(der_v.dat.size(), der_v.dat.size()+der_E.dat.size()-1)) = arma::rowvec(der_E.dat);

	  // for (int m=0; m<Nv; m++)
	  //   {
	  //     if (der_v.dat[m]!=0.0)
	  // 	{
	  // 	  Gloc(0,c) = i_data;
	  // 	  Gloc(1,c) = m;
	  // 	  Gval(c) = der_v.dat[m];
	  // 	  c++;
	  // 	}
	  //   }
	  // for (int m=Nv; m<NvE; m++)
	  //   {
	  //     if (der_E.dat[m]!=0.0)
	  // 	{
	  // 	  Gloc(0,c) = i_data;
	  // 	  Gloc(1,c) = m;
	  // 	  Gval(c) = der_E.dat[m];
	  // 	  c++;
	  // 	}
	  //   }
	  
	  // tim = clock() - tim;
	  // std::cout << "time elapsed: " << double(tim)/CLOCKS_PER_SEC << " sec\n";
	  // increment index
	  i_data++;
	}
      
      for (int p=0; p<path[k].passive.size(); p++)
	{
	  // find index of real source (vritual receiver)
	  i_s = path[k].passive(p);
	  std::cout << "   ... from passive source " << i_s <<" | ";
	  // extract arrival time
	  d(i_data) = grid.node[model.events[i_s].ind].T + model.events[i_s].t0;b
	  std::cout << "Dt(obs-calc) = " << d_obs(i_data) - d(i_data) << "\n";
	  // trace ray
	  ray.trace(grid, model.events[i_s].x, model.events[i_s].y, model.events[i_s].z);
	  // interpolate derivative on the coarse grid
	  ray.compute_interp_der(model.n, der_v, der_E);
	  // assign derivatives to elements of G:
	  
	  //tim = clock();
	  
	  G(i_data, arma::span(0, der_v.dat.size()-1)) = arma::rowvec(der_v.dat);
	  G(i_data, arma::span(der_v.dat.size(), der_v.dat.size()+der_E.dat.size()-1)) = arma::rowvec(der_E.dat);


	  // derivatives with respect to source positions	  
	  G(i_data, NvE+4*p)   = ray.der_x0;
	  G(i_data, NvE+4*p+1) = ray.der_y0;
	  G(i_data, NvE+4*p+2) = ray.der_z0;
	  G(i_data, NvE+4*p+3) = 1.0;

	  // for (int m=0; m<Nv; m++)
	  //   {
	  //     if (der_v.dat[m]!=0.0)
	  // 	{
	  // 	  Gloc(0,c) = i_data;
	  // 	  Gloc(1,c) = m;
	  // 	  Gval(c) = der_v.dat[m];
	  // 	  c++;
	  // 	}	      
	  //   }	  
	  // for (int m=Nv; m<NvE; m++)
	  //   {
	  //     if (der_E.dat[m]!=0.0)
	  // 	{
	  // 	  Gloc(0,c) = i_data;
	  // 	  Gloc(1,c) = m;
	  // 	  Gval(c) = der_E.dat[m];
	  // 	  c++;
	  // 	}
	  //   }

	  // Gloc(0,c) = i_data;
	  // Gloc(1,c) = NvE+4*p;
	  // Gval(c) = ray.der_x0;
	  // c++;

	  // Gloc(0,c) = i_data;
	  // Gloc(1,c) = NvE+4*p + 1;
	  // Gval(c) = ray.der_y0;
	  // c++;

	  // Gloc(0,c) = i_data;
	  // Gloc(1,c) = NvE+4*p + 2;
	  // Gval(c) = ray.der_z0;
	  // c++;

	  // Gloc(0,c) = i_data;
	  // Gloc(1,c) = NvE+4*p + 3;
	  // Gval(c) = 1.0;
	  // c++;

	  // tim = clock() - tim;
	  // std::cout << "time elapsed: " << double(tim)/CLOCKS_PER_SEC << " sec\n";
	  
	  // increment index
	  i_data++;
	  
	}
    }
  
  // //remove what is not needed:
  // Gloc.shed_cols(c, Gloc.n_cols-1);
  // Gval.shed_rows(c, Gval.n_rows-1);
  
  // assign G !!
  // std::cout << "Build G...\n";
  // G = arma::sp_mat(Gtmp);
  // G = arma::sp_mat(Gloc, Gval, Nd, Nm, true, true);

  // Gloc.save("Gloc.txt", arma::raw_ascii);
  // Gval.save("Gval.txt", arma::raw_ascii);
}

/// Method to compute least-square residual.
/*!\param model an object of class M containing prior and current model vecotr and inverse of coviarance matrix.
 *\return the least square residual.
 */
double D::get_residual(M& model)
{
  double r;
  r = arma::as_scalar(
		      (d.t()-d_obs.t())*Cdi*(d-d_obs)
		      + (model.m.t()-model.m_prior.t())
		      *model.Cmi
		      *(model.m-model.m_prior)
		      );
  return r;
}

/// Default constructor
QuasiNewton::QuasiNewton(){}

/// Constructor with assignment of parameter values.
/*!\param param an object of class Parameters containing input parameters and filenames.
 */
QuasiNewton::QuasiNewton(const Parameters& param)
{
  initialise(param);
}

/// Method to initialise parameter values.
/*!\param param an object of class Parameters containing input values.
 */
void QuasiNewton::initialise(const Parameters& param)
{
  residualfile = param.residualfile;
  Cmpostfile = param.Cmpostfile;
  maxit = param.maxit;
  tol = param.tol;
  stepsize = param.mu;
}

/// Method to perform one iteration of the QuasiNewton inversion procedure. During this process, the jacobian matrix G is computed and stored in the corresponding field of the class QuasiNewton, and the least-sqaure residual is also computed and stored in a vector. At the end of the iteration, the model and data are exported to external files, so that one can monitor the progress after each iteration.
/*!\param model an object of class M containing the model vector, covariances etc.
 *\param data an object of class D containing the data vector and method to compute forward model and jacobian.
 */
void QuasiNewton::iterate(M& model, D& data)
{
  
  
  // proceed with inversion step:
  std::cout << "Proceed with inversion step:\n";
  std::cout << "     (1) Build matrix...\n";
    
  arma::mat A =  G.t()*data.Cdi*G + model.Cmi;

  std::cout << "     (2) Build vector...\n";
    
  arma::vec y = G.t()*data.Cdi*(data.d-data.d_obs)
    + model.Cmi*(model.m-model.m_prior);
  
  std::cout << "     (3) Solve " << A.n_rows << " by " << A.n_cols << " linear system...\n";
    
  dm = arma::solve(A,y);
  
  std::cout << "Done !\n";
      
  model.m = model.m - stepsize*dm;
      
  // look at model update:
  std::cout << "mean absolute change in model= " << arma::mean(arma::abs(stepsize*dm)) << "\n";

  // now update model internal parameters: (i.e., convert vector m into structures usable to compute forward model again)
  std::cout << "Update model parameters from model vector...\n";
  model.update_internals();


}

/// Method to run a set number of QuasiNewton iterations. this method ends when either the max number of iterations has been reached, or when the residual becomes smaller than the set tolerance.
/*!\param model an object of class M corresponding to the model space and properties.
 *\param data an object of class D corresponding to the data space and properties.
 */
void QuasiNewton::run(M& model, D& data)
{

  // initialise residual
  double res;
  double Dres = -1.0;
  
  
  for (int it=0; it<=maxit; it++)
    {

      
      //bool success;
      std::cout << "Compute forward model...\n";
      data.forward_model(model, G);

      // compute residual
      res = data.get_residual(model);
      std::cout << "Residual is " << res << "\n";

      // append the new residual:
      residual.resize(residual.n_elem+1);
      residual(residual.n_elem-1) = res;

      //save residual into external file
      residual.save(residualfile, arma::raw_ascii);

      // compute change in residual if possible
      if (it>0)
	Dres = res - residual(it-1);

      // test if we continue:
      if (res<=tol || Dres>0.0)
	{
	  std::cout << "Stop now...\n";
	  break;
	}
      else
	{
	  // iterate QN step
	  iterate(model, data);
	  
	  // export posterior matrices and vectors
	  std::cout << "Export to files...\n";
	  export_posterior(model, data);
	}
    }

  // when finished, compute posterior covariance matrix (this takes some time because we are looking at a full matrix inversion):
  std::cout << "Computing posterior covariance matrix...\n";
  Cmpost = arma::inv_sympd( G.t()*data.Cdi*G + model.Cmi );
  Cmpost.save(Cmpostfile, arma::raw_binary);
  
}

/// Method to export all relevant vectors (m, d) and posterior covariance matrix into text files, for use elsewhere (notably, to compute the posterior movie).
/*!\param model an object of class M corresponding to the model space and properties.
 *\param data an object of class D corresponding to the data space and properties.
 */
void QuasiNewton::export_posterior(M& model, D& data)
{
  std::cout << "   ... m vector...\n";
  model.export_m();

  std::cout << "   ... (v,E) and events...\n";
  model.export_data();

  std::cout << "   ... data...\n";
  data.export_all();
}
