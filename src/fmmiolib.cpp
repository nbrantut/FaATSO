#include <cstdlib>
#include <cmath> 
#include <vector>
#include <iostream>
#include <cstdio>

#include "fmmiolib.hpp"
#include "fmmlib.hpp"


/// basic class constructor
Data::Data()
{
  Nx=Ny=Nz=0;
  h=0.0;
}

/// class constructor
/*!\param n1 number of grid points in the x direction.
 *\param n2 number of grid points in the y direction.
 *\param n3 number of grid points in the z direction.
 *\param s grid spacing.
 */
Data::Data(const int& n1, const int& n2, const int&n3, const double& s)
{
  Nx=n1;
  Ny=n2;
  Nz=n3;
  h=s;
  dat.resize(Nx*Ny*Nz);
}

/// member function that loads data from a formatted binary source.
/*! The file should be a formatted binary which contains the following sequence:
 * 3 integer types (number of grid points Nx, Ny, Nz),
 * 1 double type (grid spacing h),
 * Ntot double types (Ntot = Nx*Ny*Nz) which are the data attached to every grid point.
 *
 * \param filename a character string containing the name of file to read.
 */
void Data::load(const char *filename) 
{
  FILE * file;
  long lsize;
  int hsize;
  int Ntot;
  int k; 
  double tmp; 

  file = fopen(filename,"rb");

  if (file==NULL) 
    {
      std::cout << "File error!\n";
      exit (1);
    }

  hsize =  sizeof(int)*3 + sizeof(double);

  fread (&Nx, sizeof(int), 1, file);
  fread (&Ny, sizeof(int), 1, file);
  fread (&Nz, sizeof(int), 1, file);
  fread (&h, sizeof(double), 1, file);

  Ntot = Nx*Ny*Nz;

  // obtain total file size:
  fseek (file , 0 , SEEK_END);
  lsize = ftell (file) - hsize;
  fseek (file, hsize, SEEK_SET);

  if ( lsize/sizeof(double) == Ntot )
    {
      for (k=0; k<Ntot; k++)
	{
	  fread(&tmp, sizeof(double), 1, file);
	  dat.push_back(tmp);
	}
      fclose(file);
    }
  else
    {
      std::cout << "File error: wrong number of entries\n";
      std::cout << "File size after header is " << lsize/sizeof(double) << "\n";
      std::cout << "It should be " << Ntot << "\n";

      fclose(file);

      exit (1);
    }
}

/// member function that saves data into formatted binary source.
/*! The file generated is a formatted binary which contains the following sequence:
 * 3 integer types (number of grid points Nx, Ny, Nz),
 * 1 double type (grid spacing h),
 * Ntot double types (Ntot = Nx*Ny*Nz) which are the data attached to every grid point.
 *
 * \param filename a character string containing the name of file to write.
 */
void Data::save(const char *filename)
{
  FILE * file;
  int k;
  double tmp;

  file = fopen(filename,"wb");
  
  fwrite (&Nx, sizeof(int), 1, file);
  fwrite (&Ny, sizeof(int), 1, file);
  fwrite (&Nz, sizeof(int), 1, file);
  fwrite (&h, sizeof(double), 1, file);

  for (k=0; k<dat.size(); k++)
    {
      tmp = dat[k];
      fwrite (&tmp, sizeof(double), 1, file);
    }
  fclose(file);
}

/// member function that uses the data stored in dat, based on a given (coarse) grid specified by Nx, Ny, Nz and h, to generate a new Data object containing interpolated values of dat on a refined grid.
/*! The refinement factor is given by an integer, to be chosen so that the subdivision produces an integer number of new nodes, i.e., n should divide (Nx-1), (Ny-1) and (Nz-1).
 *\param n the refinement factor.
 *\return a Data object with the refined (linearly interpolated) values of dat.
 */ 
Data Data::refine(int n)
{
  Data fine;
  int Nxf, Nyf, Nzf;
  fine.Nx = Nxf = (Nx-1)*n + 1;
  fine.Ny = Nyf = (Ny-1)*n + 1;
  fine.Nz = Nzf = (Nz-1)*n + 1;
  fine.h = h/n;

  // this is the linear interpolation, looping through all indices on the new refined grid:
  for (int ind=0; ind<Nxf*Nyf*Nzf; ind++)
    {
      // Notations here follow those of the wikipedia page on trilinear interpolation.
      int I[8], i_000, i_001, i_010, i_011, i_100, i_101, i_110, i_111;
      double xyz[3], x, omx, y, omy, z, omz;
      double c00=0.0, c01=0.0, c10=0.0, c11=0.0;
      double c0, c1, c;
      
      get_interp_ind_coef(I,xyz,ind,Nxf,Nyf,Nzf,n);

      x = xyz[0];
      omx = 1.0 - x;
      y = xyz[1];
      omy = 1.0 - y;
      z = xyz[2];
      omz = 1.0 - z;

      i_000 = I[0];
      i_100 = I[1];
      i_010 = I[2];
      i_110 = I[3];
      i_001 = I[4];
      i_101 = I[5];
      i_011 = I[6];
      i_111 = I[7];

      // I have to separate the cases x=0 or (1-x)=0 here; otherwise, if I have a node where dat[i_xxx]=+/-Inf, the result is a nan (I get to compute things like 0*Inf). So the cases x=0 or 1 are best tackles in the simplest way by direct assignations to values at nodes.
      if (x==0)
	{
	  c00 = dat[i_000];
	  c10 = dat[i_010];
	  c01 = dat[i_001];
	  c11 = dat[i_011];
	}
      else if (omx==0)
	{
	  c00 = dat[i_100];
	  c10 = dat[i_110];
	  c01 = dat[i_101];
	  c11 = dat[i_111];
	}
      else
	{
	  c00 = omx*dat[i_000] + x*dat[i_100];
	  c10 = omx*dat[i_010] + x*dat[i_110];
	  c01 = omx*dat[i_001] + x*dat[i_101];
	  c11 = omx*dat[i_011] + x*dat[i_111];
	}
      
      if (y==0)
	{
	  c0 = c00;
	  c1 = c01;
	}
      else if (omy==0)
	{
	  c0 = c10;
	  c1 = c11;
	}
      else
	{
	  c0 = omy*c00 + y*c10;
	  c1 = omy*c01 + y*c11;
	}
      
      if (z==0)
	c = c0;
      else if (omz==0)
	c = c1;
      else
	c = omz*c0 + z*c1;
      
      fine.dat.push_back(c);
    }

  return fine;
}

/// method to obtain surrounding node indices and interpolation coefficients based on a linear index in a refined grid of size Nxf, Nyf, Nzf (refinement factor n).
/*!\param II an array of 8 integers to store the surrounding nodes indices. this is the output value (C style).
 *\param xyz an array of 3 doubles to store the relative (x,y,z) position of the node (in the fine grid) of interest within the cube formed by the 8 neighbouring nodes in the coarse grid. this is the output value (C style).
 *\param ind the linear index of node of interest.
 *\param Nxf the number of nodes along the x direction.
 *\param Nyf the number of nodes along the y direction.
 *\param Nzf the number of nodes along the z direction.
 *\param n the refinement factor.
 */
void Data::get_interp_ind_coef(int* II, double* xyz, const int& ind, const int& Nxf, const int& Nyf, const int& Nzf, const int& n)
{
  // Note to self: I is an array of size 8, xyz is an array of size 3.
  int i, j, k, I, J, K;
  int i_000, i_001, i_010, i_011, i_100, i_101, i_110, i_111;
  double x, y, z;

  // first find the subscript indices corresponding to the linear index ind IN THE FINE GRID.
  ind2sub(ind, i, j, k, Nxf, Nyf);

  //  std::cout << ind << "=(" << i << "," << j << "," << k << ")\n";

  // then by construction:
  I = i/n;
  J = j/n;
  K = k/n;

  x = ((double) (i%n))/n;
  y = ((double) (j%n))/n;
  z = ((double) (k%n))/n;

  xyz[0] = x;
  xyz[1] = y;
  xyz[2] = z;

  if ((I<Nx) && (J<Ny) && (K<Nz))
    {
      i_000 = sub2ind(I,J,K);
      i_100 = sub2ind(I+1,J,K);
      i_010 = sub2ind(I,J+1,K);
      i_001 = sub2ind(I,J,K+1);
      i_110 = sub2ind(I+1,J+1,K);
      i_101 = sub2ind(I+1,J,K+1);
      i_011 = sub2ind(I,J+1,K+1);
      i_111 = sub2ind(I+1,J+1,K+1);
    }
  else
    {
      int safeIplus = std::min(I+1,Nx);
      int safeJplus = std::min(J+1,Ny);
      int safeKplus = std::min(K+1,Nz);
      
      i_000 = sub2ind(I,J,K);
      i_100 = sub2ind(safeIplus,J,K);
      i_010 = sub2ind(I,safeJplus,K);
      i_001 = sub2ind(I,J,safeKplus);
      i_110 = sub2ind(safeIplus,safeJplus,K);
      i_101 = sub2ind(safeIplus,J,safeKplus);
      i_011 = sub2ind(I,safeJplus,safeKplus);
      i_111 = sub2ind(safeIplus,safeJplus,safeKplus);      
    }

  II[0] = i_000;
  II[1] = i_100;
  II[2] = i_010;
  II[3] = i_110;
  II[4] = i_001;
  II[5] = i_101;
  II[6] = i_011;
  II[7] = i_111;

  return;
  
}

/// method to convert linear to subscript indices, based on the grid parameters specified in the Data structure. Same as the method of the same name in the Grid object.
/*!\param ind the linear index to be converted.
 *\param i the resulting subscript index along the x direction.
 *\param j the resulting subscript index along the y direction.
 *\param k the resulting subscript index along the z direction.
 */
void Data::ind2sub(const int& ind, int& i, int& j, int& k)
{
  int vk;
        
  vk = int ( std::fmod(ind, Nx*Ny) ) ;
        
  k = (ind - vk)/(Nx*Ny) ;
  i = int ( std::fmod(vk, Nx) ) ;
  j =  (vk - i)/Nx ;
        
  return;
}

/// method to convert linear to subscript indices, based on the grid parameters specified in the arguments passed to the method. Otherwise same as the method of the same name in the Grid object.
/*!\param ind the linear index to be converted.
 *\param i the resulting subscript index along the x direction.
 *\param j the resulting subscript index along the y direction.
 *\param k the resulting subscript index along the z direction.
 *\param Nxf the number of nodes along the x direction.
 *\param Nyf the number of nodes along the y direction.
 */
void Data::ind2sub(const int& ind, int& i, int& j, int& k, const int& Nxf, const int& Nyf)
{
  int vk;
        
  vk = int ( std::fmod(ind, Nxf*Nyf) ) ;
        
  k = (ind - vk)/(Nxf*Nyf) ;
  i = int ( std::fmod(vk, Nxf) ) ;
  j =  (vk - i)/Nxf ;
        
  return;
}

/// method to convert subscript indices to linear indices. Same as the method of the same name in the Grid object.
/*!\param i subscript index along the x direction.
 *\param j subscript index along the y direction.
 *\param k subscript index along the z direction.
 *\return ind linear index.
 */
int Data::sub2ind(const int& i, const int& j, const int& k)
{
  int ind;
  if ( (i>=0) && (i<Nx) && (j>=0) && (j<Ny) && (k>=0) && (k<Nz) )
    ind = i + j*Nx + k*Nx*Ny;
  else
    ind = -1;
  return ind;
}

/// default constructor for class Source. Sets the (x,y,z) position to (0,0,0) and the linear index to 0.
locations::Location::Location()
{
  x = 0.0;
  y = 0.0;
  z = 0.0;
  t0 = 0.0;
  ind = 0;
}

/// constructor for class Location.
/*!\param xi x coordinate
 *\param yi y coordinate
 *\param zi z coordinate
 */
locations::Location::Location(const double& xi, const double& yi, const double& zi, const double& t0i)
{
  x=xi;
  y=yi;
  z=zi;
  t0=t0i;
}

/// method to determine the source index based on a given Grid.
/*!
 * \param grd the grid in which the source is.
 */
void locations::Location::findind(Grid& grd)
{
  double xmax = (grd.Nx-1)*grd.h;
  double ymax = (grd.Ny-1)*grd.h;
  double zmax = (grd.Nz-1)*grd.h;
  int i, j, k;


  // check that source is inside the computational grid
  if ( x>xmax || x<0.0 || y>ymax || y<0.0 || z>zmax || z<0.0)
    {
      std::cout << "Location is out of the grid; using closest grid point instead\n";
      x = std::max( std::min(x, xmax) , 0.0);
      y = std::max( std::min(y, ymax) , 0.0);
      z = std::max( std::min(z, zmax) , 0.0);
    }

  // find nearest neighbour in the grid:
  i = (int) std::round(x/grd.h);
  j = (int) std::round(y/grd.h);
  k = (int) std::round(z/grd.h);

  // assign the linear index of the source
  ind = grd.sub2ind(i, j, k);
}

/// method to determine the source index based on a given Grid dimensions.
/*!
 * \param Nx number of grid points in the x direction.
 * \param Ny number of grid points in the y direction.
 * \param Nz number of grid points in the z direction.
 * \param h grid spacing.
 */
void locations::Location::findind(const int &Nx, const int& Ny, const int& Nz, const double &h)
{
  double xmax = (Nx-1)*h;
  double ymax = (Ny-1)*h;
  double zmax = (Nz-1)*h;
  int i, j, k;


  // check that source is inside the computational grid
  if ( x>xmax || x<0.0 || y>ymax || y<0.0 || z>zmax || z<0.0)
    {
      std::cout << "Location is out of the grid; using closest grid point instead\n";
      x = std::max( std::min(x, xmax) , 0.0);
      y = std::max( std::min(y, ymax) , 0.0);
      z = std::max( std::min(z, zmax) , 0.0);
    }

  // find nearest neighbour in the grid:
  i = (int) std::round(x/h);
  j = (int) std::round(y/h);
  k = (int) std::round(z/h);

  // assign the linear index of the source
  ind = i + j*Nx + k*Nx*Ny;
}


/// method to determine the x,y,z position of an event based on a Grid.
/*!
 * \param grd the grid in which the source is.
 */
void locations::Location::findxyz(Grid& grd)
{
  x = grd.node[ind].i*grd.h;
  y = grd.node[ind].j*grd.h;
  z = grd.node[ind].k*grd.h;
}

/// method to read source positions in an ascii file (space-separated x y z coordinates, one line per source).
/*!
 * \param filename name of the file to read
 * \return a vector of Location
 */
std::vector<locations::Location> locations::readlocations(const char *filename)
{
  std::vector<locations::Location> loc_vec;
  double x, y, z, t0;
  FILE * file;

  file = fopen(filename,"r");

  if (file==NULL) 
    {
      std::cout << "File error!\n";
      exit (1);
    }

  while (fscanf(file,"%lf %lf %lf %lf", &x, &y, &z, &t0)==4)
    {
      //std::cout << x << "," << y << "," << z << "\n";
      loc_vec.push_back(Location(x,y,z,t0));
    }
  
  fclose(file);
  
  return loc_vec;

}

void locations::writelocations(const char *filename, const std::vector<Location> loc_vec)
{
  FILE * file;
  size_t k=0;
  
  file = fopen(filename, "w");

  if (file==NULL)
    {
      std::cout << "File error!\n";
      exit(1);
    }

  while (k<loc_vec.size())
    {
      fprintf(file, "%lf %lf %lf %lf\n", loc_vec[k].x, loc_vec[k].y, loc_vec[k].z, loc_vec[k].t0);
      k++;
    }

  fclose(file);
}

std::vector<double> times::readtimes(const char* filename)
{
  FILE * file;
  std::vector<double> t_vec;
  double t;
  
  file = fopen(filename,"r");

  if (file==NULL) 
    {
      std::cout << "File error!\n";
      exit (1);
    }

  while (fscanf(file,"%lf", &t)==1)
    {
      t_vec.push_back(t);
    }
  
  fclose(file);
  
  return t_vec;
  
}

void times::writetimes(const char* filename, std::vector<double> t_vec)
{
  FILE *file;
  size_t k=0;
  
  file = fopen(filename,"w");

  if (file==NULL)
    {
      std::cout << "file error!\n";
      exit(1);
    }

  while (k<t_vec.size())
    {
      fprintf(file, "%lf\n", t_vec[k]);
      k++;
    }

  fclose(file);
}
