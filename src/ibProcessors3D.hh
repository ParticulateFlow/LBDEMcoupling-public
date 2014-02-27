/*
 * ibProcessors3D.hh
 */

#include "ibDef.h"

namespace plb{

  /*
   * implementation of GetExternalScalarFunctional3D
   */

  template<typename T1, template<typename U> class Descriptor, typename T2>
  void GetExternalScalarFunctional3D<T1,Descriptor,T2>::process(Box3D domain, 
                                                                BlockLattice3D<T1,Descriptor>& lattice,
                                                                ScalarField3D<T2>& field){
    for(plint i=domain.x0;i<=domain.x1;i++)
      for(plint j=domain.y0;j<=domain.y1;j++)
        for(plint k=domain.z0;k<=domain.z1;k++){
          field.get(i,j,k) = *(lattice.get(i,j,k).getExternal(scalar));
        }
  }

  template<typename T1, template<typename U> class Descriptor, typename T2>
  GetExternalScalarFunctional3D<T1,Descriptor,T2>* GetExternalScalarFunctional3D<T1,Descriptor,T2>::clone() const{
    return new GetExternalScalarFunctional3D<T1,Descriptor,T2>(*this);
  }

  template<typename T1, template<typename U> class Descriptor, typename T2>
  void GetExternalScalarFunctional3D<T1,Descriptor,T2>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
  {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
  }

  /* ---------------------------------------------------- */

  /*
   * implementation of SetSingleSphere3D
   */

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
            T *uPtr = cell.getExternal(Descriptor<T>::ExternalField::boundaryVelocityBeginsAt);

          
            T dx = xGlobal - x[0];
            T dy = yGlobal - x[1];
            T dz = zGlobal - x[2];
            T sf = calcSolidFraction(dx,dy,dz,r);
            
            plint decFlag = (sf > SOLFRAC_MIN) + 2*(*sfPtr > SOLFRAC_MIN);


            switch(decFlag){
            case 0:
              // do nothing
              break;
            case 1:
              setValues(cell,sf,dx,dy,dz);
              break;
            case 2:
              if( (plint)(*idPtr) == id )
                setToZero(cell);
              // else do nothing
              break;
            case 3:
              if( sf >= *sfPtr || (plint)(*idPtr) == id )
                setValues(cell,sf,dx,dy,dz);
              // else do nothing
            }

          }
        }
      }
    
  }

  template<typename T, template<typename U> class Descriptor>
  T SetSingleSphere3D<T,Descriptor>::calcSolidFraction(T dx_, T dy_, T dz_, T r_)
  {
    plint const slicesPerDim = 6;
    
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

    T solFrac(0.);
    for(plint i=0;i<slicesPerDim;i++){
      T dx = dx_+i*sliceWidth;
      T dx_sq(dx*dx);
      for(plint j=0;j<slicesPerDim;j++){
        T dy = dy_+j*sliceWidth;
        T dy_sq(dy*dy);
        for(plint k=0;k<slicesPerDim;k++){
          T dz = dz_+k*sliceWidth;
          T dz_sq(dz*dz);
          if(dx_sq + dy_sq + dz_sq < r_sq)
            solFrac += fraction;
        }
      }
    }
    return solFrac;
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
    //bool initVelFlag = false;
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
    *idPtr = (T) -1;
  }

  /*
   * implementation of SumForceTorque3D
   */

  /* --------------------------------------------- */
  template<typename T, template<typename U> class Descriptor>
  SumForceTorque3D<T,Descriptor>::SumForceTorque3D(plint nPart_, T **x_)
    : x(x_)
  {
    
    /*
     * convention: for particle with id n, the forces are tracked in 
     * fx: n - fy: n+1 - fz: n+2
     * and the torques are tracked in
     * tx: n+3 - ty: n+4 - tz: n+5
     */
    for(plint i=0;i<6*nPart_;i++)
      sumId.push_back(this->getStatistics().subscribeSum());
  }

  template<typename T, template<typename U> class Descriptor>
  void SumForceTorque3D<T,Descriptor>::process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice)
  {
    static plint partId = Descriptor<T>::ExternalField::particleIdBeginsAt,
      fx = Descriptor<T>::ExternalField::hydrodynamicForceBeginsAt,
      fy = Descriptor<T>::ExternalField::hydrodynamicForceBeginsAt+1,
      fz = Descriptor<T>::ExternalField::hydrodynamicForceBeginsAt+2,
      //solfrac_id = Descriptor<T>::ExternalField::volumeFractionBeginsAt;
      solfrac_id = Descriptor<T>::ExternalField::bBeginsAt;

    Dot3D relativePosition = lattice.getLocation();
    // std::cout << "process " << global::mpi().getRank() << " | " 
    //           << domain.x0 << " " << domain.x1 << " "
    //           << domain.y0 << " " << domain.y1 << " "
    //           << domain.z0 << " " << domain.z1 << " | "
    //           << relativePosition.x << " " 
    //           << relativePosition.y << " " 
    //           << relativePosition.z << " | "
    //           << domain.x0 + relativePosition.x << " " << domain.x1 + relativePosition.x << " "
    //           << domain.y0 + relativePosition.y << " " << domain.y1 + relativePosition.y << " " 
    //           << domain.z0 + relativePosition.z << " " << domain.z1 + relativePosition.z << " "
    //           << std::endl;
    
    plint nx = lattice.getNx(), ny = lattice.getNy(), nz = lattice.getNz();

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
      for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
        for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
          Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
                    
          // LIGGGHTS indices start at 1
          plint id = (plint) *(cell.getExternal(partId))-1;
          if(id < 0) continue; // no particle here
          
          T xGlobal = (T) (relativePosition.x + iX);
          T yGlobal = (T) (relativePosition.y + iY);
          T zGlobal = (T) (relativePosition.z + iZ);
          
          T fs = *(cell.getExternal(solfrac_id));

          T forceX = fs*(*(cell.getExternal(fx)));
          T forceY = fs*(*(cell.getExternal(fy)));
          T forceZ = fs*(*(cell.getExternal(fz)));
        
          // pcout << id << " " << forceX << " " << forceY << " " << forceZ << std::endl;


          addForce(id,0,forceX);
          addForce(id,1,forceY);
          addForce(id,2,forceZ);

          // TODO: get torque evaluation right for periodic boundary conditions
          T dx = xGlobal - x[id][0];
          T dy = yGlobal - x[id][1];
          T dz = zGlobal - x[id][2];

          if(dx>nx/2) dx -= nx; if(dx<-nx/2) dx += nx;
          if(dy>ny/2) dy -= ny; if(dy<-ny/2) dy += ny;
          if(dz>nz/2) dz -= nz; if(dz<-nz/2) dz += nz;

          // minimum image convention
          addTorque(id,0,dy*forceZ - dz*forceY);
          addTorque(id,1,-dx*forceZ + dz*forceX);
          addTorque(id,2,dx*forceY - dy*forceX);

          // addTorque(id,0,0.);
          // addTorque(id,1,0.);
          // addTorque(id,2,0.);
        }
      }
    }
    // std::cout << global::mpi().getRank() << " | " << s << std::endl;
    // pcout << "###### ";
    // for(plint i=0;i<this->getStatistics().getSumVect().size();i++)
    //   pcout << this->getStatistics().getSumVect()[i] << " ";
    // pcout << std::endl;
  }

  template<typename T, template<typename U> class Descriptor>
  void SumForceTorque3D<T,Descriptor>::addForce(plint partId, plint coord, T value)
  {
    plint which = sumId[6*(partId)+coord];
    this->getStatistics().gatherSum(which,value);
  }
  template<typename T, template<typename U> class Descriptor>
  void SumForceTorque3D<T,Descriptor>::addTorque(plint partId, plint coord, T value)
  {
    plint which = sumId[6*(partId)+coord+3];
    this->getStatistics().gatherSum(which,value);
  }

  template<typename T, template<typename U> class Descriptor>
  double SumForceTorque3D<T,Descriptor>::getForce(plint partId, plint coord)
  {
    plint which = sumId[6*(partId)+coord];
    return this->getStatistics().getSum(which);
  }
  template<typename T, template<typename U> class Descriptor>
  double SumForceTorque3D<T,Descriptor>::getTorque(plint partId, plint coord)
  {
    plint which = sumId[6*(partId)+coord+3];
    return this->getStatistics().getSum(which);
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
    //    modified[1] = modif::staticVariables;
  }
};
