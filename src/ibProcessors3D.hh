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
      Dot3D relativePosition = lattice.getLocation();
     
      for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
        for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
          for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
            Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
	
            T xGlobal = (T) (relativePosition.x + iX);
            T yGlobal = (T) (relativePosition.y + iY);
            T zGlobal = (T) (relativePosition.z + iZ);
            
            T *sfPtr = cell.getExternal(Descriptor<T>::ExternalField::volumeFractionBeginsAt);
            T *idPtr = cell.getExternal(Descriptor<T>::ExternalField::particleIdBeginsAt);
          
            T dx = xGlobal - x[0];
            T dy = yGlobal - x[1];
            T dz = zGlobal - x[2];
            T sf = calcSolidFraction(dx,dy,dz,r);

            plint decFlag = (sf > SOLFRAC_MIN) + 2*(*sfPtr > SOLFRAC_MIN);

            switch(decFlag){
            case 0: // sf == 0 && *sfPtr == 0
              setToZero(cell);
              break; // do nothing
            case 1: // sf > 0 && *sfPtr == 0
              setValues(cell,sf,dx,dy,dz);
              break;
            case 2: // sf == 0 && *sfPtr > 0
              if( (plint)(*idPtr) == id )
                setToZero(cell);
              break; // else do nothing
            case 3: // sf > 0 && *sfPtr > 0
              if( sf > *sfPtr || (plint)(*idPtr) == id )
                setValues(cell,sf,dx,dy,dz);
              break; // else do nothing
            }

          }
        }
      }
      
  }

  template<typename T, template<typename U> class Descriptor>
  T SetSingleSphere3D<T,Descriptor>::calcSolidFraction(T dx_, T dy_, T dz_, T r_)
  {
    plint const slicesPerDim = 5;
    T const sliceWidth = 1./((T)slicesPerDim-1);
    T const fraction = 1./((T)slicesPerDim*slicesPerDim*slicesPerDim);
    
    if (dx_*dx_ + dy_*dy_ + dz_*dz_ > (r_+2)*(r_+2))
      return 0;

    if (dx_*dx_ + dy_*dy_ + dz_*dz_ < (r_-2)*(r_-2))
      return 1;

    T r_sq = r_*r_;

    dx_ = dx_ - 0.5;
    dy_ = dy_ - 0.5;
    dz_ = dz_ - 0.5;

    plint n(0);
    for(plint i=0;i<slicesPerDim;i++){
      T dx = dx_+i*sliceWidth;
      T dx_sq(dx*dx);
      for(plint j=0;j<slicesPerDim;j++){
        T dy = dy_+j*sliceWidth;
        T dy_sq(dy*dy);
        for(plint k=0;k<slicesPerDim;k++){
          T dz = dz_+k*sliceWidth;
          T dz_sq(dz*dz);
          n += (plint)(dx_sq + dy_sq + dz_sq < r_sq);
        }
      }
    }
    return fraction*((T)n);
  }

  template<typename T, template<typename U> class Descriptor>
  void SetSingleSphere3D<T,Descriptor>::setValues(Cell<T,Descriptor>& c, T sf, T dx, T dy, T dz)
  {
    T *sfPtr = c.getExternal(Descriptor<T>::ExternalField::volumeFractionBeginsAt);
    T *idPtr = c.getExternal(Descriptor<T>::ExternalField::particleIdBeginsAt);
    T *uPtr = c.getExternal(Descriptor<T>::ExternalField::boundaryVelocityBeginsAt);
    
    uPtr[0] = v[0];
    uPtr[1] = v[1];
    uPtr[2] = v[2];
    if(omega != 0){
      uPtr[0] += omega[1]*dz - omega[2]*dy;
      uPtr[1] += -omega[0]*dz + omega[2]*dx; 
      uPtr[2] += omega[0]*dy - omega[1]*dx;
    }
    *sfPtr = sf;
    *idPtr = (T) id;
    
    if(initVelFlag && sf > SOLFRAC_MAX){
      Array<T,3> u(uPtr[0],uPtr[1],uPtr[2]);
      T uSqr = uPtr[0]+uPtr[0] + uPtr[1]+uPtr[1] + uPtr[2]+uPtr[2];
      for(plint i=0;i<Descriptor<T>::q;i++)
        c[i] = c.computeEquilibrium(i,0.,u,uSqr);
    }
  }
  
  template<typename T, template<typename U> class Descriptor>
  void SetSingleSphere3D<T,Descriptor>::setToZero(Cell<T,Descriptor>& c)
  {
    T *sfPtr = c.getExternal(Descriptor<T>::ExternalField::volumeFractionBeginsAt);
    T *idPtr = c.getExternal(Descriptor<T>::ExternalField::particleIdBeginsAt);
    T *uPtr = c.getExternal(Descriptor<T>::ExternalField::boundaryVelocityBeginsAt);
    
    uPtr[0] = 0;
    uPtr[1] = 0;
    uPtr[2] = 0;
    *sfPtr = 0;
    *idPtr = (T) 0;//-1;
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
    static plint partId = Descriptor<T>::ExternalField::particleIdBeginsAt,
      fx = Descriptor<T>::ExternalField::hydrodynamicForceBeginsAt,
      fy = Descriptor<T>::ExternalField::hydrodynamicForceBeginsAt+1,
      fz = Descriptor<T>::ExternalField::hydrodynamicForceBeginsAt+2;

    Dot3D relativePosition = lattice.getLocation();
    
    plint nx = lattice.getNx(), ny = lattice.getNy(), nz = lattice.getNz();

    plint r=global::mpi().getRank();
    
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
      for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
        for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
          
          Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
          // LIGGGHTS indices start at 1
          plint id = (plint) round(*(cell.getExternal(partId)));

          if(id < 1) continue; // no particle here
          plint ind = wrapper.lmp->atom->map(id);
          
          T xGlobal = (T) (relativePosition.x + iX);
          T yGlobal = (T) (relativePosition.y + iY);
          T zGlobal = (T) (relativePosition.z + iZ);
          
          T forceX = (*(cell.getExternal(fx)));
          T forceY = (*(cell.getExternal(fy)));
          T forceZ = (*(cell.getExternal(fz)));
          
          addForce(ind,0,forceX);
          addForce(ind,1,forceY);
          addForce(ind,2,forceZ);
          
          //TODO: get torque evaluation right for periodic boundary conditions
          T dx = xGlobal - x[ind][0];
          T dy = yGlobal - x[ind][1];
          T dz = zGlobal - x[ind][2];
          pcerr << dx << " " << dy << " " << dz << std::endl;
          // minimum image convention
          // if(dx>nx/2) dx -= nx; if(dx<-nx/2) dx += nx;
          // if(dy>ny/2) dy -= ny; if(dy<-ny/2) dy += ny;
          // if(dz>nz/2) dz -= nz; if(dz<-nz/2) dz += nz;
          
          addTorque(ind,0,dy*forceZ - dz*forceY);
          addTorque(ind,1,-dx*forceZ + dz*forceX);
          addTorque(ind,2,dx*forceY - dy*forceX);
          
        }
      }
    }
  }
        

  template<typename T, template<typename U> class Descriptor>
  void SumForceTorque3D<T,Descriptor>::addForce(plint partId, plint coord, T value)
  {
    force[3*partId+coord] += value;
  }
  template<typename T, template<typename U> class Descriptor>
  void SumForceTorque3D<T,Descriptor>::addTorque(plint partId, plint coord, T value)
  {
    pcerr << 3*partId+coord << " " << torque[3*partId+coord] << " " << value;
    torque[3*partId+coord] += value;
    pcerr << " " << torque[3*partId+coord] << std::endl;
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
      T r; // id is of type T because external variables are of this type
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

    pcerr << "getforcesfromlattice ";
    for(plint iPart=0;iPart<nPart;iPart++){
      int tag = wrapper.lmp->atom->tag[iPart];
      int liggghts_ind = wrapper.lmp->atom->map(tag);
      for(plint i=0;i<3;i++){
        f_liggghts[liggghts_ind][i] = units.getPhysForce(force[3*iPart+i]);
        t_liggghts[liggghts_ind][i] = units.getPhysTorque(torque[3*iPart+i]);
        pcerr << units.getPhysTorque(torque[3*iPart+i]) << " ";
      }
      pcerr << std::endl;
    }
    couplingFix->comm_force_torque();

  }

  
};
