/*
 * ibProcessors3D.hh
 */
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

#include "ibDef.h"
#include "ibDynamics3D.h"

#include "lammps.h"
#include "atom.h"
#include "modify.h"
#include "fix_lb_coupling_onetoone.h"

namespace plb{

  /*
   * implementation of SetSingleSphere3D
   */

  template<typename T, template<typename U> class Descriptor>
  void SetSingleSphere3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
  {
    //modified[0] = modif::staticVariables;
    modified[0] = modif::nothing; 
    // modifies stuff, but we don't want updates running
    // after each particle set
  }

  template<typename T, template<typename U> class Descriptor>
  SetSingleSphere3D<T,Descriptor>* SetSingleSphere3D<T,Descriptor>::clone() const 
  {
    return new SetSingleSphere3D<T,Descriptor>(*this);
  }

  template<typename T, template<typename U> class Descriptor>
  void SetSingleSphere3D<T,Descriptor>::process(Box3D domain, BlockLattice3D<T,Descriptor> &lattice)
  {
    Dot3D const relativePosition = lattice.getLocation();

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
      for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
        for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
          Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
          Dynamics<T,Descriptor> *dyn = &(cell.getDynamics());

          // no composite --> no IB
          if(!dyn->isComposite()) continue;
          if(dyn->isBoundary())
            dyn = &(static_cast<CompositeDynamics<T,Descriptor>* >(dyn))->getBaseDynamics();

          IBcompositeDynamics<T,Descriptor> *cDyn = 
            static_cast< IBcompositeDynamics<T,Descriptor>* >( dyn );
          
          IBdynamicsParticleData<T,Descriptor> &particleData = cDyn->particleData;
          
          T const xGlobal = (T) (relativePosition.x + iX);
          T const yGlobal = (T) (relativePosition.y + iY);
          T const zGlobal = (T) (relativePosition.z + iZ);
          
          T const sf_old = particleData.solidFraction;
          int const id_old = (int) particleData.partId;
          
          
          T const dx = xGlobal - x[0];
          T const dy = yGlobal - x[1];
          T const dz = zGlobal - x[2];
          T const sf = calcSolidFraction(dx,dy,dz,r);
          
          // pcout << iX << " " << iY << " " << iZ << " | "
          //       << id_old << " " << id << " | " << sf_old << " " << sf << std::endl;
          
          plint const decFlag = (sf > SOLFRAC_MIN) + 2*(sf_old > SOLFRAC_MIN);
          
          // pcout << "dec " << decFlag << std::endl;
          
          switch(decFlag){
          case 0: // sf == 0 && *sfPtr == 0
            setToZero(particleData);
            break; // do nothing
          case 1: // sf > 0 && *sfPtr == 0
            setValues(particleData,sf,dx,dy,dz);
            break;
          case 2: // sf == 0 && *sfPtr > 0
            if( id_old == id )
              setToZero(particleData);
            break; // else do nothing
          case 3: // sf > 0 && *sfPtr > 0
            if( sf > sf_old || id_old == id )
              setValues(particleData,sf,dx,dy,dz);
            break; // else do nothing
          }
          // if desired, initialize interior of sphere with sphere velocity
          if(initVelFlag && sf > SOLFRAC_MAX)
            cell.defineVelocity(particleData.uPart);
        }
      }
    }
    
  }

  template<typename T, template<typename U> class Descriptor>
  T SetSingleSphere3D<T,Descriptor>::calcSolidFraction(T const dx_, T const dy_, T const dz_, T const r_)
  {
    // return calcSolidFractionRec(dx_,dy_,dz_,r_);

    static plint const slicesPerDim = 5;
    static T const sliceWidth = 1./((T)slicesPerDim-1);
    static T const fraction = 1./((T)slicesPerDim*slicesPerDim*slicesPerDim);
    
    static const T sqrt3half = (T) sqrt(3)/2.;

    T const dist = dx_*dx_ + dy_*dy_ + dz_*dz_;

    if (dist > (r_+sqrt3half)*(r_+sqrt3half)) return 0;
    if (dist < (r_-sqrt3half)*(r_-sqrt3half)) return 1;

    T const r_sq = r_*r_;
    T dx_sq[slicesPerDim],dy_sq[slicesPerDim],dz_sq[slicesPerDim];

    for(plint i=0;i<slicesPerDim;i++){
      T delta = -0.5 + ((T)i)*sliceWidth;
      T dx = dx_+delta; dx_sq[i] = dx*dx;
      T dy = dy_+delta; dy_sq[i] = dy*dy;
      T dz = dz_+delta; dz_sq[i] = dz*dz;
    }

    plint n(0);
    for(plint i=0;i<slicesPerDim;i++){
      for(plint j=0;j<slicesPerDim;j++){
        for(plint k=0;k<slicesPerDim;k++){
          if(dx_sq[i] + dy_sq[j] + dz_sq[k] < r_sq) n++;
        }
      }
    }

    return fraction*((T)n);
  }

  template<typename T, template<typename U> class Descriptor>
  T SetSingleSphere3D<T,Descriptor>::calcSolidFractionRec(T const dx_ ,T const dy_, T const dz_,
                                                        T const r_, plint const recursionLevel = 1)
  {
    static const plint maxRecDepth = 3;
    static const T sqrt3half = (T) sqrt(3)/2.;
    
    const T len = 1./((T)recursionLevel);
    const T vol = len*len*len;
    const T cutoffLen = sqrt3half*len;

    T distSq = dx_*dx_ + dy_*dy_ + dz_*dz_;

    if (recursionLevel == maxRecDepth){
      T r_sq = r_*r_;
      if(distSq < r_sq) return vol;
      else return 0;
    }

    if (distSq > (r_+cutoffLen)*(r_+cutoffLen)) return 0;
    if (distSq < (r_-cutoffLen)*(r_-cutoffLen)) return vol;

    const T shift = len/4.;
    plint recursionLevelNew = recursionLevel + 1;

    T dx_plus = dx_+shift, dx_min = dx_+shift;
    T dy_plus = dy_+shift, dy_min = dy_+shift;
    T dz_plus = dz_+shift, dz_min = dz_+shift;

    return calcSolidFractionRec(dx_plus,dy_plus,dz_plus,r_,recursionLevelNew)
      + calcSolidFractionRec(dx_min,dy_plus,dz_plus,r_,recursionLevelNew)
      + calcSolidFractionRec(dx_plus,dy_min,dz_plus,r_,recursionLevelNew)
      + calcSolidFractionRec(dx_min,dy_min,dz_plus,r_,recursionLevelNew)
      + calcSolidFractionRec(dx_plus,dy_plus,dz_min,r_,recursionLevelNew)
      + calcSolidFractionRec(dx_min,dy_plus,dz_min,r_,recursionLevelNew)
      + calcSolidFractionRec(dx_plus,dy_min,dz_min,r_,recursionLevelNew)
      + calcSolidFractionRec(dx_min,dy_min,dz_min,r_,recursionLevelNew);

  }


  template<typename T, template<typename U> class Descriptor>
  void SetSingleSphere3D<T,Descriptor>::setValues(IBdynamicsParticleData<T,Descriptor> &p,
                                                  T const sf, T const dx, T const dy, T const dz)
  {    
    p.uPart.from_cArray(v);
    if(omega != 0){
      p.uPart[0] += omega[1]*dz - omega[2]*dy;
      p.uPart[1] += -omega[0]*dz + omega[2]*dx; 
      p.uPart[2] += omega[0]*dy - omega[1]*dx;
    }
    p.solidFraction = sf;
    p.partId = id;
  }
  
  template<typename T, template<typename U> class Descriptor>
  void SetSingleSphere3D<T,Descriptor>::setToZero(IBdynamicsParticleData<T,Descriptor> &p)
  {
    p.uPart[0] = 0;
    p.uPart[1] = 0;
    p.uPart[2] = 0;
    p.solidFraction = 0;
    p.partId = 0;
  }

  /*
   * implementation of SumForceTorque3D
   */

  /* --------------------------------------------- */
  template<typename T, template<typename U> class Descriptor>
  SumForceTorque3D<T,Descriptor>::SumForceTorque3D(typename ParticleData<T>::ParticleDataArrayVector &x_,
                                                   T *force_, T *torque_, LiggghtsCouplingWrapper &wrapper_)
    : x(x_), force(force_), torque(torque_), wrapper(wrapper_)
  {}

  template<typename T, template<typename U> class Descriptor>
  SumForceTorque3D<T,Descriptor>::SumForceTorque3D(SumForceTorque3D<T,Descriptor> const &orig)
    : BoxProcessingFunctional3D_L<T,Descriptor>(orig), x(orig.x),
      force(orig.force), torque(orig.torque), wrapper(orig.wrapper) {}
  
  template<typename T, template<typename U> class Descriptor>
  void SumForceTorque3D<T,Descriptor>::process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice)
  {
    Dot3D const relativePosition = lattice.getLocation();
    
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
      for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
        for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
          
          Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
          Dynamics<T,Descriptor> *dyn = &(cell.getDynamics());


          // no composite --> no IB
          if(!dyn->isComposite()) continue;
          if(dyn->isBoundary())
            dyn = &(static_cast<CompositeDynamics<T,Descriptor>* >(dyn))->getBaseDynamics();

          IBcompositeDynamics<T,Descriptor> *cDyn = 
            static_cast< IBcompositeDynamics<T,Descriptor>* >( dyn );

          IBdynamicsParticleData<T,Descriptor> particleData = cDyn->particleData;
          // LIGGGHTS indices start at 1
          plint const id = particleData.partId;

          if(id < 1) continue; // no particle here

          plint const ind = wrapper.lmp->atom->map(id);
          
          T const xGlobal = (T) (relativePosition.x + iX);
          T const yGlobal = (T) (relativePosition.y + iY);
          T const zGlobal = (T) (relativePosition.z + iZ);
          
          T const forceX = particleData.hydrodynamicForce[0];
          T const forceY = particleData.hydrodynamicForce[1];
          T const forceZ = particleData.hydrodynamicForce[2];
          
          addForce(ind,0,forceX);
          addForce(ind,1,forceY);
          addForce(ind,2,forceZ);
          
          //TODO: get torque evaluation right for periodic boundary conditions
          T const dx = xGlobal - x[ind][0];
          T const dy = yGlobal - x[ind][1];
          T const dz = zGlobal - x[ind][2];

          addTorque(ind,0,dy*forceZ - dz*forceY);
          addTorque(ind,1,-dx*forceZ + dz*forceX);
          addTorque(ind,2,dx*forceY - dy*forceX);
          
          cDyn->particleData = particleData;
        }
      }
    }
  }
        

  template<typename T, template<typename U> class Descriptor>
  void SumForceTorque3D<T,Descriptor>::addForce(plint const partId, plint const coord, T const value)
  {
    force[3*partId+coord] += value;
  }
  template<typename T, template<typename U> class Descriptor>
  void SumForceTorque3D<T,Descriptor>::addTorque(plint const partId, plint const coord, T const value)
  {
    torque[3*partId+coord] += value;
  }

  template<typename T, template<typename U> class Descriptor>
  SumForceTorque3D<T,Descriptor>* SumForceTorque3D<T,Descriptor>::clone() const
  { 
    return new SumForceTorque3D<T,Descriptor>(*this);
  }

  template<typename T, template<typename U> class Descriptor>
  void SumForceTorque3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
  {
    modified[0] = modif::nothing;
  }

  /* ------------------------------------------- */
  /*
   * implementation of setSpheresOnLattice
   */

  template<typename T, template<typename U> class Descriptor>
  void setSpheresOnLattice(MultiBlockLattice3D<T,Descriptor> &lattice,
                           LiggghtsCouplingWrapper &wrapper,
                           PhysUnits3D<T> const &units,
                           bool initVelFlag)
  {
    plint nx=lattice.getNx(), ny=lattice.getNy(), nz=lattice.getNz();
    for(plint iS=0;iS<wrapper.getNumParticles();iS++){
      T x[3],v[3],omega[3];
      T r; 
      plint id = (plint) round( (T)wrapper.id[iS][0] + 0.1 );

      for(plint i=0;i<3;i++){
        x[i] = units.getLbLength(wrapper.x[iS][i]);
        v[i] = units.getLbVel(wrapper.v[iS][i]);
        omega[i] = units.getLbFreq(wrapper.omega[iS][i]);
      }
      r = units.getLbLength(wrapper.r[iS][0]);

      SetSingleSphere3D<T,Descriptor> *sss 
        = new SetSingleSphere3D<T,Descriptor>(x,v,omega,r,id,initVelFlag);
      Box3D sss_box = sss->getBoundingBox();
      applyProcessingFunctional(sss,sss_box,lattice);

      
      if(lattice.periodicity().get(0)){
        T x_per[3] = {x[0],x[1],x[2]};
        bool addPerImg(true);
        if(sss_box.x0 <= 0)
          x_per[0] += (T)nx;
        else if(sss_box.x1 >= nx-1)
          x_per[0] -= (T)nx;
        else 
          addPerImg = false;

        if(addPerImg){
          SetSingleSphere3D<T,Descriptor> *sss_per 
            = new SetSingleSphere3D<T,Descriptor>(x_per,v,omega,r,id,initVelFlag);
          applyProcessingFunctional(sss_per,sss_per->getBoundingBox(),lattice);
        }        
      }
      if(lattice.periodicity().get(1)){
        T x_per[3] = {x[0],x[1],x[2]};
        bool addPerImg(true);
        if(sss_box.y0 <= 0)
          x_per[1] += (T)ny;
        else if(sss_box.y1 >= ny-1)
          x_per[1] -= (T)ny;
        else 
          addPerImg = false;

        if(addPerImg){
          SetSingleSphere3D<T,Descriptor> *sss_per 
            = new SetSingleSphere3D<T,Descriptor>(x_per,v,omega,r,id,initVelFlag);
          applyProcessingFunctional(sss_per,sss_per->getBoundingBox(),lattice);
        }
      }
      if(lattice.periodicity().get(2)){
        plint nz = lattice.getNz();
        T x_per[3] = {x[0],x[1],x[2]};
        bool addPerImg(true);
        if(sss_box.z0 <= 0)
          x_per[2] += (T)nz;
        else if(sss_box.z1 >= nz-1)
          x_per[2] -= (T)nz;
        else 
          addPerImg = false;

        if(addPerImg){
          SetSingleSphere3D<T,Descriptor> *sss_per 
            = new SetSingleSphere3D<T,Descriptor>(x_per,v,omega,r,id,initVelFlag);
          applyProcessingFunctional(sss_per,sss_per->getBoundingBox(),lattice);
        }
      }
      
      if(lattice.periodicity().get(0) && lattice.periodicity().get(1)){
        T x_per[3] = {x[0],x[1],x[2]};
        bool addPerImg(true);
        if(sss_box.x0 <= 0)
          x_per[0] += (T)nx;
        else if(sss_box.x1 >= nx-1)
          x_per[0] -= (T)nx;
        else 
          addPerImg = false;
        if(addPerImg){
          if(sss_box.y0 <= 0)
            x_per[1] += (T)ny;
          else if(sss_box.y1 >= ny-1)
            x_per[1] -= (T)ny;
          else 
            addPerImg = false;
        }

        if(addPerImg){
          SetSingleSphere3D<T,Descriptor> *sss_per 
            = new SetSingleSphere3D<T,Descriptor>(x_per,v,omega,r,id,initVelFlag);
          applyProcessingFunctional(sss_per,sss_per->getBoundingBox(),lattice);
        }
        
      }
      if(lattice.periodicity().get(1) && lattice.periodicity().get(2)){
        T x_per[3] = {x[0],x[1],x[2]};
        bool addPerImg(true);
        if(sss_box.y0 <= 0)
          x_per[1] += (T)ny;
        else if(sss_box.y1 >= ny-1)
          x_per[1] -= (T)ny;
        else 
          addPerImg = false;
        if(addPerImg){
          if(sss_box.z0 <= 0)
            x_per[2] += (T)nz;
          else if(sss_box.z1 >= nz-1)
            x_per[2] -= (T)nz;
          else 
            addPerImg = false;
        }
        if(addPerImg){
          SetSingleSphere3D<T,Descriptor> *sss_per 
            = new SetSingleSphere3D<T,Descriptor>(x_per,v,omega,r,id,initVelFlag);
          applyProcessingFunctional(sss_per,sss_per->getBoundingBox(),lattice);
        }
      }
      if(lattice.periodicity().get(2) && lattice.periodicity().get(0)){
        T x_per[3] = {x[0],x[1],x[2]};
        bool addPerImg(true);
        if(sss_box.z0 <= 0)
          x_per[2] += (T)nz;
        else if(sss_box.z1 >= nz-1)
          x_per[2] -= (T)nz;
        else 
          addPerImg = false;

        if(addPerImg){
          if(sss_box.x0 <= 0)
            x_per[0] += (T)nx;
          else if(sss_box.x1 >= nx-1)
            x_per[0] -= (T)nx;
          else 
            addPerImg = false;
        }
        if(addPerImg){
          SetSingleSphere3D<T,Descriptor> *sss_per 
            = new SetSingleSphere3D<T,Descriptor>(x_per,v,omega,r,id,initVelFlag);
          applyProcessingFunctional(sss_per,sss_per->getBoundingBox(),lattice);
        }
        
      }
      if(lattice.periodicity().get(0) && lattice.periodicity().get(1) && lattice.periodicity().get(2)){
        T x_per[3] = {x[0],x[1],x[2]};
        bool addPerImg(true);

        if(sss_box.x0 <= 0)
          x_per[0] += (T)nx;
        else if(sss_box.x1 >= nx-1)
          x_per[0] -= (T)nx;
        else 
          addPerImg = false;

        if(addPerImg){
          if(sss_box.y0 <= 0)
            x_per[1] += (T)ny;
          else if(sss_box.y1 >= ny-1)
            x_per[1] -= (T)ny;
          else 
            addPerImg = false;
        }
        if(addPerImg){
          if(sss_box.z0 <= 0)
            x_per[2] += (T)nz;
          else if(sss_box.z1 >= nz-1)
            x_per[2] -= (T)nz;
          else 
            addPerImg = false;
        }
        if(addPerImg){
          SetSingleSphere3D<T,Descriptor> *sss_per 
            = new SetSingleSphere3D<T,Descriptor>(x_per,v,omega,r,id,initVelFlag);
          applyProcessingFunctional(sss_per,sss_per->getBoundingBox(),lattice);
        }
      }
    }
    // this one returns modif::staticVariables and forces an update of those along processor
    // boundaries
    applyProcessingFunctional(new AttributeFunctional<T,Descriptor>(),lattice.getBoundingBox(),lattice);
  }

  template<typename T, template<typename U> class Descriptor>
  void getForcesFromLattice(MultiBlockLattice3D<T,Descriptor> &lattice,
                            LiggghtsCouplingWrapper &wrapper,
                            PhysUnits3D<T> const &units)
  {
    static typename ParticleData<T>::ParticleDataArrayVector x_lb;
    for(plint iPart=0;iPart<wrapper.getNumParticles();iPart++)
      x_lb.push_back( Array<T,3>( units.getLbLength(wrapper.x[iPart][0]),
                                  units.getLbLength(wrapper.x[iPart][1]),
                                  units.getLbLength(wrapper.x[iPart][2]) ) );
    
    plint const n_force = wrapper.getNumParticles()*3;
    
    double *force = new T[n_force];
    double *torque = new T[n_force];
    double *force_tmp = new T[n_force];
    double *torque_tmp = new T[n_force];
    
    for(plint i=0;i<n_force;i++){
      force_tmp[i] = 0;
      torque_tmp[i] = 0;
    }
 
    SumForceTorque3D<T,Descriptor> *sft = new SumForceTorque3D<T,Descriptor>(x_lb,force_tmp,torque_tmp);
    
    applyProcessingFunctional(sft,lattice.getBoundingBox(), lattice);
    
    MPI_Allreduce(force_tmp,force,n_force,MPI_DOUBLE,MPI_SUM,
                  global::mpi().getGlobalCommunicator());
    MPI_Allreduce(torque_tmp,torque,n_force,MPI_DOUBLE,MPI_SUM,
                  global::mpi().getGlobalCommunicator());

    T const nProc = (T) global::mpi().getSize();
    for(plint iPart=0;iPart<wrapper.getNumParticles();iPart++){
      for(int i=0;i<3;i++){
        wrapper.f[iPart][i] = units.getPhysForce(force[3*iPart+i])/nProc;
        wrapper.t[iPart][i] = units.getPhysTorque(torque[3*iPart+i])/nProc;
        // division grudge because of double MPI_allreduce
      }
    }
    
    
    delete[] force;
    delete[] torque;
  }



  template<typename T, template<typename U> class Descriptor>
  void setSpheresOnLatticeNew(MultiBlockLattice3D<T,Descriptor> &lattice,
                              LiggghtsCouplingWrapper &wrapper,
                              PhysUnits3D<T> const &units,
                              std::vector<plint> &excludeType,
                              bool initVelFlag)
  {
    plint r = global::mpi().getRank();

    // this relies on the fact that there is exactly one block on each lattice
    plint iBlock = lattice.getLocalInfo().getBlocks()[0];
    std::map<plint,Box3D> blockmap = lattice.getSparseBlockStructure().getBulks();
    Box3D localBB = blockmap[iBlock];

    plint nx=lattice.getNx(), ny=lattice.getNy(), nz=lattice.getNz();
    plint nPart = wrapper.lmp->atom->nlocal + wrapper.lmp->atom->nghost;

    for(plint iS=0;iS<nPart;iS++){
      plint type = (plint)wrapper.lmp->atom->type[iS];
      bool excludeFlag(false);
      for(plint iT=0;iT<excludeType.size();iT++){
        //        pcout << iS << " " << type << " " << excludeType[iT] << std::endl;
        if(type == excludeType[iT]){
          excludeFlag = true;
          break;
        }
      }

      if(excludeFlag) continue;

      T x[3],v[3],omega[3];
      T r;
      plint id = (plint) round( (T)wrapper.lmp->atom->tag[iS] + 0.1 );

      for(plint i=0;i<3;i++){
        x[i] = units.getLbLength(wrapper.lmp->atom->x[iS][i]);
        v[i] = units.getLbVel(wrapper.lmp->atom->v[iS][i]);
        omega[i] = units.getLbFreq(wrapper.lmp->atom->omega[iS][i]);
      }
      r = units.getLbLength(wrapper.lmp->atom->radius[iS]);
      SetSingleSphere3D<T,Descriptor> *sss 
        = new SetSingleSphere3D<T,Descriptor>(x,v,omega,r,id,initVelFlag);
      Box3D sss_box = sss->getBoundingBox();
      
      // only go over part that lies on local processor
      // to avoid unnecessary communication overhead
      Box3D sss_box_intersect(0,0,0,0,0,0);
      bool boxes_intersect = intersect(sss_box,localBB,sss_box_intersect);
      // std::cout << r << " intersect " << boxes_intersect << " "
      //           << sss_box_intersect.x0 << " " << sss_box_intersect.x1 << " "
      //           << sss_box_intersect.y0 << " " << sss_box_intersect.y1 << " "
      //           << sss_box_intersect.z0 << " " << sss_box_intersect.z1 << " "
      //           << std::endl;
      if(boxes_intersect)
        applyProcessingFunctional(sss,sss_box_intersect,lattice);
      else
        delete sss;
      // applyProcessingFunctional(sss,sss_box,lattice);
    }

    // this one returns modif::staticVariables and forces an update of those along processor
    // boundaries
    applyProcessingFunctional(new AttributeFunctional<T,Descriptor>(),lattice.getBoundingBox(),lattice);

  }

  template<typename T, template<typename U> class Descriptor>
  void setSpheresOnLatticeNew(MultiBlockLattice3D<T,Descriptor> &lattice,
                              LiggghtsCouplingWrapper &wrapper,
                              PhysUnits3D<T> const &units,
                              bool initVelFlag)
  {
    std::vector<plint> dummyExcludeType;
    setSpheresOnLatticeNew(lattice,wrapper,units,dummyExcludeType,initVelFlag);
  }



  template<typename T, template<typename U> class Descriptor>
  void getForcesFromLatticeNew(MultiBlockLattice3D<T,Descriptor> &lattice,
                               LiggghtsCouplingWrapper &wrapper,
                               PhysUnits3D<T> const &units)
  {
    // debug stuff
    plint r = global::mpi().getRank();

    static std::vector<T> force,torque;
    static typename ParticleData<T>::ParticleDataArrayVector x_lb;

    plint const nPart = wrapper.lmp->atom->nlocal + wrapper.lmp->atom->nghost;
    plint const n_force = nPart*3;

    if(nPart == 0) return; // no particles - no work

    if(nPart > x_lb.size()){
      for(plint iPart=0;iPart<x_lb.size();iPart++){
        x_lb[iPart][0] = units.getLbLength(wrapper.lmp->atom->x[iPart][0]);
        x_lb[iPart][1] = units.getLbLength(wrapper.lmp->atom->x[iPart][1]);
        x_lb[iPart][2] = units.getLbLength(wrapper.lmp->atom->x[iPart][2]);
      }
      for(plint iPart = x_lb.size();iPart < nPart; iPart++)
        x_lb.push_back( Array<T,3>( units.getLbLength(wrapper.lmp->atom->x[iPart][0]),
                                    units.getLbLength(wrapper.lmp->atom->x[iPart][1]),
                                    units.getLbLength(wrapper.lmp->atom->x[iPart][2]) ) );
    } else{
      for(plint iPart=0;iPart<nPart;iPart++){
        x_lb[iPart][0] = units.getLbLength(wrapper.lmp->atom->x[iPart][0]);
        x_lb[iPart][1] = units.getLbLength(wrapper.lmp->atom->x[iPart][1]);
        x_lb[iPart][2] = units.getLbLength(wrapper.lmp->atom->x[iPart][2]);
      }
    }

    if(n_force > force.size()){
      for(plint i=0;i<force.size();i++){
        force[i] = 0;
        torque[i] = 0;
      }
      for(plint i=force.size();i<n_force;i++){
        force.push_back(0.);
        torque.push_back(0.);
      }
    } else {
      for(plint i=0;i<n_force;i++){
        force[i] = 0;
        torque[i] = 0;
      }
    }

    SumForceTorque3D<T,Descriptor> *sft = new SumForceTorque3D<T,Descriptor>(x_lb,
                                                                             &force.front(),&torque.front(),
                                                                             wrapper);

    // this relies on the fact that there is exactly one block on each processor
    plint iBlock = lattice.getLocalInfo().getBlocks()[0];
    std::map<plint,Box3D> blockmap = lattice.getSparseBlockStructure().getBulks();
    Box3D localBB = blockmap[iBlock];
    applyProcessingFunctional(sft,localBB, lattice);

    // // experimental....
    // SparseBlockStructure3D sparseBlock = lattice.getSparseBlockStructure();
    // std::vector<Box3D> boxVec;
    // std::vector<plint> localBlocks = lattice.getLocalInfo().getBlocks();
    // sparseBlock.intersect(lattice.getBoundingBox(),localBlocks,boxVec);
    // std::cout << r << " boxvec " << boxVec.size() << std::endl;
    // for(plint i=0;i<boxVec.size();i++)
    //   applyProcessingFunctional(sft,boxVec[i],lattice);
    
    // // basic version with complete lattice
    // applyProcessingFunctional(sft,lattice.getBoundingBox(), lattice);

    LAMMPS_NS::FixLbCouplingOnetoone 
      *couplingFix 
      = dynamic_cast<LAMMPS_NS::FixLbCouplingOnetoone*>
      (wrapper.lmp->modify->find_fix_style("couple/lb/onetoone",0));

    double **f_liggghts = couplingFix->get_force_ptr();
    double **t_liggghts = couplingFix->get_torque_ptr();

    for(plint iPart=0;iPart<nPart;iPart++){
      int tag = wrapper.lmp->atom->tag[iPart];
      int liggghts_ind = wrapper.lmp->atom->map(tag);
      for(plint i=0;i<3;i++){
        f_liggghts[liggghts_ind][i] = units.getPhysForce(force[3*iPart+i]);
        t_liggghts[liggghts_ind][i] = units.getPhysTorque(torque[3*iPart+i]);
      }
    }
    couplingFix->comm_force_torque();

  }

  
};
