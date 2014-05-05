/*
 * ibProcessors3D.hh
 */

#include "ibDef.h"

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
                                                   T *force_, T *torque_)
    : x(x_), force(force_), torque(torque_)
  {
    
    /*
     * convention: for particle with id n, the forces are tracked in 
     * fx: n - fy: n+1 - fz: n+2
     * and the torques are tracked in
     * tx: n+3 - ty: n+4 - tz: n+5
     */
  }

  template<typename T, template<typename U> class Descriptor>
  SumForceTorque3D<T,Descriptor>::SumForceTorque3D(SumForceTorque3D<T,Descriptor> const &orig)
    : BoxProcessingFunctional3D_L<T,Descriptor>(orig), x(orig.x),
      force(orig.force), torque(orig.torque) {}
  
  template<typename T, template<typename U> class Descriptor>
  void SumForceTorque3D<T,Descriptor>::process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice)
  {
    static plint partId = Descriptor<T>::ExternalField::particleIdBeginsAt,
      fx = Descriptor<T>::ExternalField::hydrodynamicForceBeginsAt,
      fy = Descriptor<T>::ExternalField::hydrodynamicForceBeginsAt+1,
      fz = Descriptor<T>::ExternalField::hydrodynamicForceBeginsAt+2;

    Dot3D relativePosition = lattice.getLocation();
    
    plint nx = lattice.getNx(), ny = lattice.getNy(), nz = lattice.getNz();

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
      for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
        for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
          plint i=1;
          plint r=global::mpi().getRank();
          
          Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);

          // LIGGGHTS indices start at 1
          plint id = (plint) *(cell.getExternal(partId))-1;
          if(id < 0) continue; // no particle here
          // std::cout << r << " : " << iX << " " << iY << " " << iZ << std::endl;
          // std::cout << r << " | " << i++ << std::endl;

          T xGlobal = (T) (relativePosition.x + iX);
          T yGlobal = (T) (relativePosition.y + iY);
          T zGlobal = (T) (relativePosition.z + iZ);

          // std::cout << r << " | " << i++ << std::endl;
          
          T forceX = (*(cell.getExternal(fx)));
          T forceY = (*(cell.getExternal(fy)));
          T forceZ = (*(cell.getExternal(fz)));

          // std::cout << r << " | " << i++ << std::endl;
        
          addForce(id,0,forceX);
          addForce(id,1,forceY);
          addForce(id,2,forceZ);

          // std::cout << r << " | " << i++ << std::endl;
          // TODO: get torque evaluation right for periodic boundary conditions
          T dx = xGlobal - x[id][0];
          T dy = yGlobal - x[id][1];
          T dz = zGlobal - x[id][2];

          // minimum image convention
          if(dx>nx/2) dx -= nx; if(dx<-nx/2) dx += nx;
          if(dy>ny/2) dy -= ny; if(dy<-ny/2) dy += ny;
          if(dz>nz/2) dz -= nz; if(dz<-nz/2) dz += nz;

          addTorque(id,0,dy*forceZ - dz*forceY);
          addTorque(id,1,-dx*forceZ + dz*forceX);
          addTorque(id,2,dx*forceY - dy*forceX);

          // std::cout << r << " | " << i << std::endl;
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
    typename ParticleData<T>::ParticleDataArrayVector x_lb;
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
  
};
