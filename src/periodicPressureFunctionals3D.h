/*
 * This file is part of the LBDEMcoupling software.
 *
 * LBDEMcoupling is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Copyright 2014 Johannes Kepler University Linz
 *
 * Author: Philippe Seil (philippe.seil@jku.at)
 */

#ifndef PERIODIC_PRESSURE_FUNCTIONALS_3D
#define PERIODIC_PRESSURE_FUNCTIONALS_3D

namespace plb {
  /*
   * implements a periodic pressure gradient across a periodic boundary
   * following 
   * Zhang and Kwok, Phys.Rev.E 73, 047702 (2006)
   */

  template<typename T, template<typename U> class Descriptor>
  class ZhangPeriodicPressureFunctional3D : public BoxProcessingFunctional3D_L<T,Descriptor> {
  public:
    ZhangPeriodicPressureFunctional3D(T rhoTarget_, T rhoAvg_, plint dimension_, plint direction_);
    void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice);
    virtual ZhangPeriodicPressureFunctional3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const;

  private:
    typedef std::vector<plint> IndexVec;
    IndexVec rescalePop;
    T rhoTarget,rhoAvg;
    plint dimension, direction;

  };

  /**
   * this derives directly from BoxProcessingFunctional
   * it creates a coupling between three objects: a lattice,
   * a scalar field and a tensor field.
   */

  template<typename T1, template<typename U> class Descriptor, typename T2, typename T3, int nDim>
  class BoxProcessingFunctional3D_LST : public BoxProcessingFunctional3D {
  public:
    virtual void process(Box3D domain, BlockLattice3D<T1,Descriptor> &lattice,
                         ScalarField3D<T2> &scalar, TensorField3D<T3,nDim> &tensor) =0;
    void processGenericBlocks (Box3D domain, std::vector<AtomicBlock3D*> atomicBlocks);

  };

  /*
   * implements a periodic pressure gradient across a periodic boundary
   * following 
   * Kim and Pitsch, Phys. Fluids 19, 108101(2007)
   */

  template<typename T1, template<typename U> class Descriptor, typename T2, typename T3, int nDim>
  class KimPeriodicPressureFunctional3D : public BoxProcessingFunctional3D_LST<T1,Descriptor,T2,T3,nDim> {
  public:
    KimPeriodicPressureFunctional3D(T1 rhoTarget_, T1 rhoAvg_, plint dimension_, plint direction_);
    void process(Box3D domain, BlockLattice3D<T1,Descriptor>& lattice,
                 ScalarField3D<T2> &rho, TensorField3D<T3,nDim> &u);
    virtual KimPeriodicPressureFunctional3D<T1,Descriptor,T2,T3,nDim>* clone() const;
    virtual void getTypeOfModification(std::vector<plb::modif::ModifT>& modified) const;

  private:
    typedef std::vector<plint> IndexVec;
    IndexVec rescalePop;
    T1 rhoTarget,rhoAvg;
    plint dimension, direction;
    T1 computeEquilibrium(Cell<T1,Descriptor>& cell,
                          plint iEq,
                          T2 rho,
                          Array<T3,Descriptor<T1>::d> u);  

  };

  /*
   * convenience wrapper class that contains instantiation and application 
   * of KimPeriodicPressureFunctional3D
   * usage:
   * ppm.preColl(lattice);
   * lattice.collideAndStream();
   * ppm.postColl(lattice);
   * dimension: 0 = x, 1 = y, 2 = z
   * in/out direction: +1 if boundary
   */

  template<typename T, template<typename U> class Descriptor>
  class PeriodicPressureManager {
  public:
    PeriodicPressureManager(MultiBlockLattice3D<T,Descriptor> &lattice, 
                            T rhoIn_, T rhoOut_, Box3D inlet_, Box3D outlet_, 
                            plint dimension_, plint inDirection_, plint outDirection_);
    
    void preColl(MultiBlockLattice3D<T,Descriptor> &lattice);
    void postColl(MultiBlockLattice3D<T,Descriptor> &lattice);
  private:
    T rhoIn, rhoOut, rhoAvgIn, rhoAvgOut;
    Box3D inlet,outlet,tmp;
    plint dimension, inDirection, outDirection;
    MultiScalarField3D<T> rho;
    MultiTensorField3D<T,Descriptor<T>::d> u;
  };
  
  /*
   * class PressureGradient
   * helper class to initialize a domain with a constant pressure gradient
   * usage:
   * initializeAtEquilibrium( lattice, lattice.getBoundingBox(), 
   *                          PressureGradient<T>(rhoHi,rhoLo, nGrid, dimension) );
   * nGrid: grid points in direction of gradient
   * dimension: 0 = x, 1 = y, 2 = z
   */

  template<typename T>
  class PressureGradient {
  public:
    PressureGradient(T pHi_, T pLo_, plint n_, plint dimension_);
    void operator() (plint iX, plint iY, plint iZ, T& density, Array<T,3>& velocity) const;
  private:
    T pHi, pLo;
    plint n,dimension;
  };

  template<typename T>
  class PoiseuilleProfileAndPressureGradient {
  public:
    PoiseuilleProfileAndPressureGradient(T pHi_, T pLo_,  T uMax, 
                                         plint nx_, plint ny_, plint nz_, plint dimension_);
    void operator() (plint iX, plint iY, plint iZ, T& density, Array<T,3>& velocity) const;
  private:
    T pHi,pLo,uMax;
    plint nx,ny,nz,dimension;
    T poiseuilleVelocity(T i, T n) const;
  };

}; /* namespace plb */

#include "periodicPressureFunctionals3D.hh"

#endif /* PERIODIC_PRESSURE_FUNCTIONALS_3D */
