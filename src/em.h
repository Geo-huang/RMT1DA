/*****************************************************************************/
/*                                                                           */
/*  Copyright 2023                                                           */
/*  Huang Chen                                                               */
/*  chenhuang@cqu.edu.cn                                                     */
/*                                                                           */
/*****************************************************************************/



#ifndef _NAMESPACEEM_H
#define _NAMESPACEEM_H

#define EIGEN_DONT_VECTORIZE // remove issues with STL containers
#define EIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT
//#define EIGEN_USE_MKL_ALL

#include <complex>
#include <Eigen/Dense>  // Linear Algebra Lib.
#include <Eigen/Sparse> // Sparse Lib
#include <Eigen/StdVector>

class ComplexPoint;

//----------------------------------------------------------------------------
namespace EM
{  
  // Data-type 
  typedef double Real; 
  typedef std::complex<Real> Dcomplex;  
  typedef Eigen::SparseMatrix<Dcomplex, 
          Eigen::RowMajor>                      EigenSparseMatrix; 
  typedef Eigen::VectorXcd                      EigenDenseVector;
  typedef Eigen::Matrix<Dcomplex, 3, 1>         Grad;
  typedef Eigen::Matrix<Real,     6, 6>         Matrix6D;
  typedef Eigen::Matrix<Dcomplex, 6, 6>         Matrix6C;
  typedef Eigen::Matrix<Real, 2, 2>             Matrix2D;
  typedef Eigen::Matrix<Dcomplex, 2, 2>         Matrix2C;
  typedef Eigen::Matrix<Dcomplex, 2, 1>         Vector2C;
  typedef Eigen::Matrix<Real,   3,   3>         Matrix3D; 
  typedef Eigen::Matrix<Dcomplex,   3,   3>     Matrix3C; 
  typedef Eigen::Matrix<Dcomplex,   4,   4>     Matrix4C;
  typedef Eigen::Matrix<Dcomplex, 4, 1>         Vector4C;
  typedef Eigen::Matrix<Dcomplex, 6, 1>         Vector6C;


  // Pi=3.14159... 
  static const Real           pi   = 3.1415926535897932384626433832795;
  static const Dcomplex       zero = Dcomplex(0.,0.);
  static const Dcomplex       II   = Dcomplex(0.,1.);
  static const unsigned int   INVALID_UNIT=static_cast<unsigned int>(-1);
  // changed from 1e-14 to 1e-20 by Huang Chen on 2020.12.25
  static const Real           TOL = 1e-20;  
  static const Real           mu0       =4.0*pi*1e-7; 
  static const Real           epsilon0  =8.854187817*1e-12;
  // + or - signs 
  double sgn(const double m);
}

inline 
double EM::sgn(const double m) {
  if (std::abs(m)<EM::TOL) return 0.;
  else return m/std::abs(m);
}

//---------------------------------------------
using namespace EM;

#endif // #define _NAMESPACEEM_H
