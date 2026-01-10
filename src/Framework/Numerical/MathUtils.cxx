//____________________________________________________________________________
/*
 Copyright (c) 2003-2025, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org

 Costas Andreopoulos <c.andreopoulos \at cern.ch>
 University of Liverpool
*/
//____________________________________________________________________________

#include <float.h>

#include <TMath.h>

#include "Framework/Messenger/Messenger.h"
#include "Framework/Numerical/RandomGen.h"
#include "Framework/Numerical/MathUtils.h"

//using namespace genie::utils::math;

//____________________________________________________________________________
TMatrixD genie::utils::math::CholeskyDecomposition(const TMatrixD& cov_matrix)
{
// Perform a Cholesky decomposition of the input covariance matrix and
// return the lower triangular matrix
//
   const double epsilon = 1E-12;

   int ncols = cov_matrix.GetNcols();
   int nrows = cov_matrix.GetNrows();

   assert(ncols==nrows);

   int n = nrows;

   TMatrixD L(n, n);

   for (int i = 0; i < n; ++i) {

     // calculate the diagonal term first
     L(i,i) = cov_matrix(i,i);
     for (int k = 0; k < i; ++k) {
        double tmp = L(k,i);
        L(i,i) -= tmp*tmp;
     }//k

     if(L(i,i) <= 0) {
       if(fabs(L(i,i)) < epsilon){
         L(i,i)=epsilon;
         LOG("Cholesky", pINFO)
           << "Changed element (" << i << ", " << i << ") to " << L(i,i);
       }
       else{
         LOG("Cholesky", pERROR)
            << "Decomposed covariance matrix not positive-definite";
         LOG("Cholesky", pERROR)
            << "L(" << i << "," << i << ") = " << L(i,i);
         exit(1);
       }
     }
     L(i,i) = TMath::Sqrt(L(i,i));
     // then the off-diagonal terms
     for (int j = i+1; j < n; ++j) {
        L(i,j) = cov_matrix(i,j);
        for (int k = 0; k < i; ++k) {
           L(i,j) -= L(k,i)*L(k,j);
        }
        L(i,j) /= L(i,i);
     }//j
  }//i

  // create the transpose of L
  TMatrixD LT(TMatrixD::kTransposed,L);

  return LT;
}
//____________________________________________________________________________
TVectorD genie::utils::math::CholeskyGenerateCorrelatedParams (
    const TMatrixD& cholesky_triangular, TVectorD& mean_params)
{
// Generate a vector of correlated params

  int ncols = cholesky_triangular.GetNcols();
  int nrows = cholesky_triangular.GetNrows();
  int npars = mean_params.GetNrows();

  if(ncols != nrows) {
    LOG("Cholesky", pERROR)
        << "Mismatch between number of columns (" << ncols
        << ") & rows (" << nrows << ")";
    exit(1);
  }
  if(npars != nrows) {
    LOG("Cholesky", pERROR)
        << "Mismatch between number of parameters (" << npars
        << ") & array size (" << nrows << ")";
    exit(1);
  }

  int n = nrows;

  // create a vector of unit Gaussian variables
  // and multiply by Lt to introduce the appropriate correlations
  TVectorD g(n);
  for (int k = 0; k < n; ++k) {
    g(k) = RandomGen::Instance()->RndNum().Gaus();
  }
  g *= cholesky_triangular;

  // add the mean value offsets and store the results
  TVectorD correlated_params(n);
  for (int i = 0; i < n; ++i) {
     double v = mean_params[i];
     v += g(i);
     correlated_params[i] = v;
  }

  return correlated_params;
}
//____________________________________________________________________________
TVectorD genie::utils::math::CholeskyGenerateCorrelatedParams (
 const TMatrixD& cholesky_triangular, TVectorD& mean_params, TVectorD& g_uncorrelated)
{
// Generate a vector of correlated params

  int ncols = cholesky_triangular.GetNcols();
  int nrows = cholesky_triangular.GetNrows();
  int npars = mean_params.GetNrows();
  int nunco = g_uncorrelated.GetNrows();

  if(ncols != nrows) {
    LOG("Cholesky", pERROR)
        << "Mismatch between number of columns (" << ncols
        << ") & rows (" << nrows << ")";
    exit(1);
  }
  if(npars != nrows) {
    LOG("Cholesky", pERROR)
        << "Mismatch between number of parameters (" << npars
        << ") & array size (" << nrows << ")";
    exit(1);
  }
  if(nunco != nrows) {
    LOG("Cholesky", pERROR)
        << "Mismatch between size of uncorrelated parameter vector (" << nunco
        << ") & array size (" << nrows << ")";
    exit(1);
  }

  int n = nrows;

  // create a vector of unit Gaussian variables
  // and multiply by Lt to introduce the appropriate correlations
  g_uncorrelated *= cholesky_triangular;

  // add the mean value offsets and store the results
  TVectorD correlated_params(n);
  for (int i = 0; i < n; ++i) {
     double v = mean_params[i];
     v += g_uncorrelated(i);
     correlated_params[i] = v;
  }

  return correlated_params;
}
//____________________________________________________________________________
TVectorD genie::utils::math::CholeskyGenerateCorrelatedParamVariations (
    const TMatrixD& cholesky_triangular)
{
  int ncols = cholesky_triangular.GetNcols();
  int nrows = cholesky_triangular.GetNrows();

  assert(ncols==nrows);

  int n = nrows;

  // create a vector of unit Gaussian variables
  // and multiply by Lt to introduce the appropriate correlations
  TVectorD g(n);
  for (int k = 0; k < n; ++k) {
    g(k) = RandomGen::Instance()->RndNum().Gaus();
  }
  g *= cholesky_triangular;

  return g;
}
//____________________________________________________________________________
TVectorD genie::utils::math::CholeskyCalculateCorrelatedParamVariations (
    const TMatrixD& cholesky_triangular, TVectorD & g_uncorrelated)
{
  int ncols = cholesky_triangular.GetNcols();
  int nrows = cholesky_triangular.GetNrows();
  int npars = g_uncorrelated.GetNrows();

  assert(ncols==nrows);
  assert(npars==nrows);

  // create a vector of unit Gaussian variables
  // and multiply by Lt to introduce the appropriate correlations
  TVectorD g(g_uncorrelated);
  g *= cholesky_triangular;

  return g;
}
//____________________________________________________________________________
double genie::utils::math::KahanSummation(double x[], unsigned int n)
{
// the Kahan summation algorithm - minimizes the error when adding a sequence
// of finite precision floating point numbers (compensated summation)

  double sum = x[0];
  double c   = 0.0;
  for(unsigned int i=1; i<n; i++) {
    double y = x[i]-c;
    double t = sum+y;
    c   = (t-sum) - y;
    sum = t;
  }
  return sum;
}
//____________________________________________________________________________
double genie::utils::math::KahanSummation(const vector<double> & x)
{
// the Kahan summation algorithm - minimizes the error when adding a sequence
// of finite precision floating point numbers (compensated summation)

  double sum = x[0];
  double c   = 0.0;
  for(unsigned int i=1; i<x.size(); i++) {
    double y = x[i]-c;
    double t = sum+y;
    c   = (t-sum) - y;
    sum = t;
  }
  return sum;
}
//____________________________________________________________________________
bool genie::utils::math::AreEqual(double x1, double x2)
{
  double err = 0.001*DBL_EPSILON;
  double dx  = TMath::Abs(x1-x2);
  if(dx<err) {
    LOG("Math", pINFO) << x1 << " := " << x2;
    return true;
  }
  return false;;
}
//____________________________________________________________________________
bool genie::utils::math::AreEqual(float x1, float x2)
{
  float err = FLT_EPSILON;
  float dx  = TMath::Abs(x1-x2);
  if(dx<err) {
    LOG("Math", pINFO) << x1 << " := " << x2;
    return true;
  }
  return false;;
}
//____________________________________________________________________________
bool genie::utils::math::IsWithinLimits(double x, Range1D_t range)
{
  return ( x >= range.min && x <= range.max );
}
//____________________________________________________________________________
bool genie::utils::math::IsWithinLimits(float x, Range1F_t range)
{
  return ( x >= range.min && x <= range.max );
}
//____________________________________________________________________________
bool genie::utils::math::IsWithinLimits(int i, Range1I_t range)
{
  return ( i >= range.min && i <= range.max );
}
//____________________________________________________________________________
double genie::utils::math::NonNegative(double x)
{
// this is used to handle very small numbers in sqrts

  return TMath::Max(0., x);
}
//____________________________________________________________________________
double genie::utils::math::NonNegative(float x)
{
// this is used to handle very small numbers in sqrts

  return TMath::Max( (float)0., x);
}
//____________________________________________________________________________
genie::utils::math::GaussLegendreQuadrature * 
genie::utils::math::GaussLegendreQuadrature::fInstance = 0;
//____________________________________________________________________________
genie::utils::math::GaussLegendreQuadrature * 
genie::utils::math::GaussLegendreQuadrature::Instance() {
  if( fInstance == 0 ) fInstance = new genie::utils::math::GaussLegendreQuadrature;
  return fInstance;
}
//____________________________________________________________________________
genie::utils::math::GaussLegendreQuadrature::GaussLegendreQuadrature() : fDataPath("") {
  this->Init();
}
//____________________________________________________________________________
genie::utils::math::GaussLegendreQuadrature::~GaussLegendreQuadrature() {
  this->Clear();
}
//____________________________________________________________________________
// Change the data path of the GL coefficients. 
// This should *always* be accompanied by GaussLegendreQuadrature::ReadGLFile() !!
void genie::utils::math::GaussLegendreQuadrature::SetDataPath(std::string newPath) {
  LOG("Math", pWARN) << "Setting new path for Gauss-Legendre quadrature coefficients."
		     << " Please ensure GaussLegendreQuadrature::ReadGLFile() is called!";
  fDataPath = newPath;
}
//____________________________________________________________________________
void genie::utils::math::GaussLegendreQuadrature::Init() {
  // Hard coding to find this in 
  // $GENIE/data/numerical/Gauss-Legendre_quadrature_coefficients.csv
  std::string genie_path = gSystem->Getenv("GENIE");
  genie_path += "/";
  std::string relative_data_path = "data/numerical/Gauss-Legendre_quadrature_coefficients.csv";
  fDataPath = genie_path + relative_data_path;

  // Read the coefficients into `GaussLegQuad` objects and add these to the table
  this->ReadGLFile();
}
//____________________________________________________________________________
// Returns a row from fGLTable. This is dummy (n=0) if n was not found
const genie::utils::math::GaussLegQuad 
genie::utils::math::GaussLegendreQuadrature::GetGLQuad(int n) {
  if( fGLTable.find(n) == fGLTable.end() ) {
    LOG("Math", pERROR) << "Gauss-Legendre quadrature coefficients at order = " << n
			<< " not found in data path " << fDataPath
			<< " . Returning a dummy row. Any integrals will be wrong / break.";
    return fGLTable[0];
  }
  return fGLTable[n];
}
//____________________________________________________________________________
// Simple parser to handle streaming of the CSV-ified Gauss-Legendre coefficients
/* CSV structure:
 * n, N, x_1, ..., x_N, w_1, ..., w_N, K_n
 * n = order of approximation
 * N = floor(M/2) with M = number of control points
 * x_i = control points (if X is a ctrl-pt, so is -X)
 * w_i = weights (same for +X and -X)
 * K_n = coefficient of the error term, see CRC 33rd ed 8.3.1.7
 */
void genie::utils::math::GaussLegendreQuadrature::ReadGLFile() {
  // Add a dummy row to the table for error handling
  if( fGLTable.find(0) == fGLTable.end() ) {
    genie::utils::math::GaussLegQuad dummy = 
      { 0, std::vector<double>{}, std::vector<double>{}, -1.0 };
    this->AddGL(dummy);
  } // insert dummy if needed

  LOG("MATH", pNOTICE) << "Loading Gauss-Legendre quadrature coefficients from " 
		       << fDataPath;

  std::string line; std::ifstream fin(fDataPath.c_str());
  while( std::getline( fin, line ) ) {
    std::istringstream iss(line); char comma;
    if (line[0] == '#') continue;
    double nd; if( ! (iss >> nd >> comma) ) continue;
    double Nd; if( ! (iss >> Nd >> comma) ) continue;
    if( std::floor(nd) != nd || std::floor(Nd) != Nd ) {
      LOG("Math", pWARN) << "Skipping malformed line: \n\t" << line
			 << "\nbecause the first two fields should be integers, and they're not.";
      continue;
    }
    int n = static_cast<int>(nd); int N = static_cast<int>(Nd);

    // store the control points and weights
    // Assuming the control points are ordered, the vectors will be ordered too
    std::vector<double> nodes; std::vector<double> weights;
    double first_node = -999.9;
    for( int i = 0; i < N; i++ ) {
      double x; if( ! (iss >> x >> comma) ) break;
      first_node = (first_node < -1) ? x : first_node;
      nodes.emplace_back(x); if(x!=0.0) { nodes.insert(nodes.begin(), -x); }
    } // get the nodes
    for( int i = 0; i < N; i++ ) {
      double w; if( ! (iss >> w >> comma) ) break;
      weights.emplace_back(w); if( ! (first_node==0.0 && i==0) )
				 { weights.insert(weights.begin(), w); }
    } // get the weights

    double err; if( ! (iss >> err) ) continue;
    genie::utils::math::GaussLegQuad glq = { n, nodes, weights, err };
    this->AddGL(glq);
    
  } // read lines from fDataPath
}
