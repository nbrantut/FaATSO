//safeguard library
#ifndef TOMOLIB_INCLUDED
#define TOMOLIB_INCLUDED

# include<armadillo>
# include<vector>
# include<string>

# include "fmmiolib.hpp"

/// a global constant to use in place of -inf values in ln(v) (i.e., a very large negative number but not infinite...)
const double V_SHADOW = -100.0; // remember that the wavespeed will then be exp(-100) = 3.72e-44 ... close enough from zero !

/// the class Parameters is a structure containing inversion parameters, file names, etc.
class Parameters
{
public:
  std::string folder; ///< name fo folder to store all input/saved files.
  std::string paramfile; ///< name of parameter file (ascii)
  std::string sensorfile; ///< name of sensor file (ascii)
  std::string eventfile; ///< name of event position file (ascii)
  std::string tsurveyfile; ///< name of survey arrival times file (ascii)
  std::string teventfile; ///< name of event arrival times file (ascii)
  std::string vpriorfile; ///< name of prior ln(v) structure file (.dat)
  std::string Epriorfile; ///< name of prior E structure file (.dat)
  std::string vpostfile; ///< name of posterior ln(v) structure file (.dat)
  std::string Epostfile; ///< name of posterior E structure file (.dat)
  std::string eventpostfile; ///< name of posterior event positions file (ascii)
  std::string dcalcfile; ///< name of file for theoretical arrival times (ascii)
  std::string Cmpostfile; ///< name of file for posterior covariance matrix (ascii)
  std::string mpostfile; ///< name of file for posterior model values (ascii)
  std::string residualfile; ///< name of file with residuals
  int n; ///< grid refinement factor
  double s_surveys; ///< variance on arrival times from surveys
  double s_events; ///< variance on arrival times from events
  int maxit; ///< max number of QN iterations
  double tol; ///< tolerance to exit QN loop
  double mu; ///< quasi-Newton step size
  double s_v; ///< std dev on ln(v)
  double s_E; ///< std dev on E
  double lambda; ///< correlation length
  double s_x; ///< std dev on x position of events
  double s_y; ///< std dev on x position of events
  double s_z; ///< std dev on x position of events
  double s_t0; ///< std dev on x position of events

  Parameters();
  Parameters(const std::string& );
  void initialise(const std::string& );
};

/// The class M corresponds to the model parameters in the inverse problem. It contains not only the vector of (cartesian) model parameters, but all other informations (file names, prior values, etc) and methods to build covariance matrices etc.
class M
{
public:
  arma::vec  m; ///< vector of model parameters
  arma::vec  m_prior; ///< vector of prior values of m
  arma::mat  Cmi; ///< inverse of correlation matrix
  std::vector<locations::Location> events; ///< vector storing event positions
  Data v0; ///< data structure containing ln(v) on a grid
  Data E0; ///< data structure containing E on a grid
  Data vf; ///< data structure containing ln(v) on a refined grid
  Data Ef; ///< data structure containing E on a refined grid
  int n; ///< grid refinement factor
  
  M();
  M(const Parameters&);

  void initialise(const Parameters&);
  void update_internals();
  void export_data();
  void export_m();

private:
  double lambda; ///< correlation length
  double s_v; ///< std dev on ln(v)
  double s_E; ///< std dev on E
  double s_x; ///< std dev on x position of events
  double s_y; ///< std dev on x position of events
  double s_z; ///< std dev on x position of events
  double s_t0; ///< std dev on x position of events
  arma::uvec ind_m_null; ///< vector of indices where v=0 (m=inf)

  std::string vpostfile; ///< name of posterior ln(v) structure file (.dat)
  std::string Epostfile; ///< name of posterior E structure file (.dat)
  std::string eventpostfile; ///< name of posterior event positions file (ascii)
  std::string mpostfile; ///< name of file to export vector m_post

  void build_inv_covariance();
};

/// The class Sources is essentially a map containing the indices of the active and passive sources for which data exist.
class Sources
{
public:
  arma::uvec active; ///< an array of the sensors indices for which valid data exist.
  arma::uvec passive; ///< an array of the events indices for which valid data exist.
};


/// The class D corresponds to the data/observables and all the associated informations in the inverse problem.
class D
{
public:
  arma::vec d;
  arma::vec d_obs;
  arma::mat Cdi;
  std::vector<Sources> path; ///<this is a vector of Sources elements, of size equal to the number of sensors (virtual sources).
  std::vector<locations::Location> sensors; ///< vector storing sensor positions

  D();
  D(const Parameters&, const M&);

  void initialise(const Parameters&, const M&);
  void export_all();
  void forward_model(M &, arma::mat &);
  double get_residual(M&);
  

private:
  int Ndshots; ///< number of data corresponding to shots
  int Ndevents; ///< number of data corresponding to events
  std::string dcalcfile; ///< filename whare to save theoretical data
  double s_surveys; ///< std dev on arrival times from surveys
  double s_events; ///< std dev on arrival times from events

  void build_inv_covariance();
  void build_data_vector(std::string, std::string, int, int);
};

/// The class QuasiNewton corresponds to the method and associated tools to solve the inverse problem. It is intended as a high-level wrapper to run iterations and exports to files, using objetcs of class M and D as inputs.
class QuasiNewton
{
public:
  arma::mat G; ///< contains the current Jacobian
  arma::vec dm; ///< contains the current step in m
  arma::vec residual; ///< contains the sequence of residuals
  arma::mat Cmpost; ///< posterior covariance matrix on m

  QuasiNewton();
  QuasiNewton(const Parameters&);
 
  void initialise(const Parameters&);

  void run(M&, D&);
  void iterate(M&, D&);
  void export_posterior(M&, D&);
  
private:
  std::string residualfile; ///< name of file which stores residuals after each iteration
  std::string Cmpostfile; ///< name of file which stores posterior covariance matrix on model parameters.
  int maxit; ///< max number if iterations
  double tol; ///< tolerance of the QN algorithm
  double stepsize; ///< QN stepsize
    
};

#endif
