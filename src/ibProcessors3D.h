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
#ifndef IB_PROCESSORS_3D_H
#define IB_PROCESSORS_3D_H

#include "liggghtsCouplingWrapper.h"
#include "physunits.h"

namespace plb{
  template<typename T>
  struct ParticleData {
    typedef typename std::vector<Array<T,3> > ParticleDataArrayVector;
    typedef typename std::vector<T> ParticleDataScalarVector;
  };

  template<typename T, template<typename U> class Descriptor>
  struct SetSingleSphere3D : public BoxProcessingFunctional3D_L<T,Descriptor> {
  public:
    SetSingleSphere3D(T *x_, T *v_, T *omega_, T r_, int id_, bool initVelFlag_)
      : x(x_),v(v_),omega(omega_),r(r_),id(id_), initVelFlag(initVelFlag_) {}
    SetSingleSphere3D(T *x_, T *v_, T *omega_, T r_, int id_)
      : x(x_),v(v_),omega(omega_),r(r_),id(id_), initVelFlag(false) {}
    SetSingleSphere3D(T *x_, T *v_, T r_, int id_)
      : x(x_),v(v_),omega(0),r(r_),id(id_), initVelFlag(false) {}
    SetSingleSphere3D(const SetSingleSphere3D &orig)
      : x(orig.x),v(orig.v),omega(orig.omega),r(orig.r),
      id(orig.id),initVelFlag(orig.initVelFlag) {}

    Box3D getBoundingBox() const { return Box3D(x[0]-r-2,x[0]+r+2,
                                                x[1]-r-2,x[1]+r+2,
                                                x[2]-r-2,x[2]+r+2);}

    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor> &lattice);
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    SetSingleSphere3D<T,Descriptor>* clone() const;
  private:
    T *x,*v,*omega,r;
    int id;
    bool initVelFlag;
    T calcSolidFraction(T const dx_, T const dy_, T const dz_, T const r_);
    T calcSolidFractionRec(T const dx_ ,T const dy_, T const dz_, T const r_,
                           plint const recursionLevel);
    void setValues(IBdynamicsParticleData<T,Descriptor> &p, 
                   T const sf, T const dx, T const dy, T const dz);
    void setToZero(IBdynamicsParticleData<T,Descriptor> &p);
  };
  
  /*
   * this functional does nothing, it is just used to get 
   * all the static variables right after looping over 
   * the spheres
   */
  template<typename T, template<typename U> class Descriptor>
  struct AttributeFunctional : public BoxProcessingFunctional3D_L<T,Descriptor> {
  public:
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor> &lattice) {}
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const
    {modified[0] = modif::staticVariables;}

    AttributeFunctional<T,Descriptor>* clone() const
    { return new AttributeFunctional<T,Descriptor>(*this);}
  };
  
  // template<typename T, template<typename U> class Descriptor>
  // void setSpheresOnLattice(MultiBlockLattice3D<T,Descriptor> &lattice,
  //                          T **x, T **v, T **omega, T *r, int **id, plint nSpheres, bool initVelFlag);
  template<typename T, template<typename U> class Descriptor>
  void setSpheresOnLattice(MultiBlockLattice3D<T,Descriptor> &lattice,
                           LiggghtsCouplingWrapper &wrapper,
                           PhysUnits3D<T> const &units,
                           bool initVelFlag);


  template<typename T, template<typename U> class Descriptor>
  void getForcesFromLattice(MultiBlockLattice3D<T,Descriptor> &lattice,
                            LiggghtsCouplingWrapper &wrapper,
                            PhysUnits3D<T> const &units);

  template<typename T, template<typename U> class Descriptor>
  void setSpheresOnLatticeNew(MultiBlockLattice3D<T,Descriptor> &lattice,
                              LiggghtsCouplingWrapper &wrapper,
                              PhysUnits3D<T> const &units,
                              std::vector<plint> &excludeType,
                              bool initVelFlag);

  template<typename T, template<typename U> class Descriptor>
  void setSpheresOnLatticeNew(MultiBlockLattice3D<T,Descriptor> &lattice,
                              LiggghtsCouplingWrapper &wrapper,
                              PhysUnits3D<T> const &units,
                              bool initVelFlag);

  template<typename T, template<typename U> class Descriptor>
  void getForcesFromLatticeNew(MultiBlockLattice3D<T,Descriptor> &lattice,
                               LiggghtsCouplingWrapper &wrapper,
                               PhysUnits3D<T> const &units);

  template<typename T, template<typename U> class Descriptor>
  struct SumForceTorque3D : public BoxProcessingFunctional3D_L<T,Descriptor> {
  public:
    SumForceTorque3D(typename ParticleData<T>::ParticleDataArrayVector &x,
                     T *force_, T *torque_, LiggghtsCouplingWrapper &wrapper_);
    SumForceTorque3D(SumForceTorque3D<T,Descriptor> const &orig);

    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice);

    SumForceTorque3D<T,Descriptor>* clone() const;
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const;

  private:
    typename ParticleData<T>::ParticleDataArrayVector &x;
    T *force,*torque;
    LiggghtsCouplingWrapper &wrapper;
    void addForce(plint const partId, plint const coord, T const value);
    void addTorque(plint const partId, plint const coord, T const value);
  };

  template<typename T1, template<typename U> class Descriptor, typename T2>
  struct GetScalarQuantityFromDynamicsFunctional : public BoxProcessingFunctional3D_LS<T1,Descriptor,T2> {
  public:

    enum Quantity { SolidFraction, ParticleId };
    Quantity which;
    GetScalarQuantityFromDynamicsFunctional(Quantity const which_)
      : which(which_) {}
    
    virtual void process(Box3D domain, BlockLattice3D<T1,Descriptor>& lattice,
                         ScalarField3D<T2>& solfrac) 
    {
      for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
          for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
            
            Cell<T1,Descriptor>& cell = lattice.get(iX,iY,iZ);
            Dynamics<T1,Descriptor> *dyn = &(cell.getDynamics());
            T1 val = 0;
            if(dyn->isComposite()){
              IBcompositeDynamics<T1,Descriptor> *cDyn = 
                static_cast< IBcompositeDynamics<T1,Descriptor>* >( dyn );
              switch(which){
              case SolidFraction:
                val = (T2) cDyn->particleData.solidFraction;
                break;
              case ParticleId:
                val = (T2) cDyn->particleData.partId;
                break;

              }

            }
            solfrac.get(iX,iY,iZ) = val;
          }
        }
      }
    }

    GetScalarQuantityFromDynamicsFunctional<T1,Descriptor,T2>* clone() const
    {
      return new GetScalarQuantityFromDynamicsFunctional<T1,Descriptor,T2>(*this);
    }
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const
    {
      modified[0] = modif::nothing;
      modified[1] = modif::staticVariables;
    }

    
  };

}; // plb

#include "ibProcessors3D.hh"

#endif // IB_PROCESSORS_3D_H
