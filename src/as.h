/***********************************************************************************
 //File Name	  : as.h
 //Author         : Huang Chen
 //Version        : 1.1
 //Copyright      : (C) 2023 Huang Chen, CQU Chongqing
 //Mail           : chenhuang@cqu.edu.cn
 //Created Time   : 22/03/18 21:36:02
 //Description    : Analytical solutions for 1D anisotropic MT/RMT problems
 **********************************************************************************/

// This class offers electricmagnetic responses (eg. rho, phase) at a given depth d 
// for a n-layers anisotropic Earth model
// with an vertical incident plance wave at a given frequency (f). 
// N is the number of layers below the air space (N>=1). 
// N=1, the half-space model.

// usage:  
// exp(-i*w*t)
// +z downward, z=0 is the flat air-earth-interface
// For isotropic layer:
// TE-mode,  E=(Ex,0,0),  H=(0,Hy,Hz);
// TM-model, E=(0,Ey,Ez), H=(Hx,0,0);
// For anistropic layer: the so-called TE-TM model are compled,
// and E=(Ex,Ey,Ez_-,Ez_+); H=(Hx, Hy, Hz_-, Hz_+)
// written by Huang Chen, 01/19/2018

#ifndef _AS_H
#define _AS_H

#include <vector>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include "em.h"

//----------------------------------------------------------------------
class AS
{
public:
  // Constructor 
  AS(const std::vector<std::vector<Real> >& electric_parameters, const Real frequency);
  // Destructor
  ~AS();

  // init
  void init();
  // transform electric_parameters to wavenumber, admittivity, impedivity etc.
  void inner_parameters(const std::vector<std::vector<Real> >& ep, const Real frequency);
 
  // Compute impedance at the the layered interface by impedance propagation.
  void compute_impedance();
  // Compute responses inclding apparent resistivity, phase, etc. at the ground.
  void post_process(std::ofstream &out);

 public:
  const Real                   f;        // observation frequencies
  unsigned int                 N;        // N-layer earth model
  std::vector<Matrix3D>        sigma;    // conductivity tensor of each layer
  std::vector<Matrix3D>        epsilon;  // epsilon tensor
  std::vector<Matrix3C>        y_hat;    // admittivity
  std::vector<Dcomplex>        z_hat;    // z_hat = i*omega*mu
  std::vector<Matrix2C>        A;        // auxiliary matrix A for sigma tensor 
  std::vector<double>          depth;    // depth of the bottom of every layer, depth = zi
                                         // to the air-Earth interface(z=0).
  std::vector<double>          h;        // thickness of each layer, h>0      
  std::vector<Dcomplex>        K1;       // K1 = sqrt(-i*omega*u*A1)
  std::vector<Dcomplex>        K2;       // K2 = sqrt(-i*omega*u*A2) 
  std::vector<Dcomplex>        Q1;       // Q1 = (K1^2+i*omega*u*Axx)/i*omega*u*Axy
  std::vector<Dcomplex>        Q2;       // Q2 = (K2^2+i*omega*u*Axx)/i*omega*u*Axy
  std::vector<Dcomplex>        gama1;    // gama1 = -K1/(i*omega*u)
  std::vector<Dcomplex>        gama2;    // gama2 = -K2/(i*omega*u)
  std::vector<Dcomplex>        D1;       // D1 = 1+e^-2(K1+K2)h-e^-2K1*h-e^-2K2*h
  std::vector<Dcomplex>        D2;       // D2 = 1+e^-2(K1+K2)h-e^+2K1*h-e^+2K2*h
  std::vector<Dcomplex>        D3;       // D3 = 1-e^-2(K1+K2)h-e^+2K1*h-e^-2K2*h
  std::vector<Dcomplex>        D4;       // D4 = 1-e^-2(K1+K2)h-e^-2K1*h-e^+2K2*h
  std::vector<Dcomplex>        D5;       // D4 = 4e^-(K1+K2)h  

  std::vector<Vector6C>        P;        // The following vectors are used for calculating of det and M   
  std::vector<Vector6C>        U1;     
  std::vector<Vector6C>        U2; 
  std::vector<Vector6C>        U3;     
  std::vector<Vector6C>        U4; 

  std::vector<Matrix2C>        Z;        // the surface impedance of l-th layer, Z[l-1], l = 1,...,N
                                         // l = 1, means observation point on the ground.
  Matrix2C                     ZN;       // the ground impedance of N-th layer.

};


//------------------------Implementation-------------------------------
inline
AS::AS(const std::vector<std::vector<Real> >& electric_parameters,
                 const Real frequency): f(frequency)
{
	
  assert(std::abs(frequency)>0);
  this->N = electric_parameters.size();
  assert(N>=1);
  this->sigma.resize(N+1);
  this->epsilon.resize(N+1);
  this->y_hat.resize(N+1);
  this->z_hat.resize(N+1);
  this->A.resize(N+1);
  this->depth.resize(N+1);
  this->h.resize(N+1);
  this->K1.resize(N+1);
  this->K2.resize(N+1);
  this->Q1.resize(N+1);
  this->Q2.resize(N+1);
  this->gama1.resize(N+1);
  this->gama2.resize(N+1);
  this->D1.resize(N+1);
  this->D2.resize(N+1);
  this->D3.resize(N+1);
  this->D4.resize(N+1);
  this->D5.resize(N+1);

  this->P.resize(N+1);
  this->U1.resize(N+1);
  this->U2.resize(N+1);
  this->U3.resize(N+1);
  this->U4.resize(N+1);

  // All the first element, e.g. U4[0], of the above vectors are 
  // not used in this program because we don't take the air layer into account.
  // The following two elements of h vector are not used as well.
  this->h[0] = 1E50;
  this->h[N] = 1E50;

  this->Z.resize(N);

  // Prepare k, z_hat, y_hat.
  this->inner_parameters(electric_parameters,frequency);


}

inline
AS::~AS()
{

}

#endif //_AS_H

