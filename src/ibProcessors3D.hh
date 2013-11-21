/*
 * ibProcessors3D.hh
 */

#include "ibDef.h"

namespace plb{

  /*
   * implementation of GetExternalScalarFunctional3D
   */

  template<typename T1, template<typename U> class Descriptor, typename T2>
  void GetExternalScalarFunctional3D<T1,Descriptor,T2>::process(Box3D domain, BlockLattice3D<T1,Descriptor>& lattice,
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
   * implementation of SetSpheres3D
   */
  template<typename T, template<typename U> class Descriptor>
  T SetSpheres3D<T,Descriptor>::calcDistSqr(T x0, T y0, T z0, T x1, T y1, T z1)
  {
    T dx = x0-x1; T dy = y0-y1; T dz = z0-z1; 
    return dx*dx + dy*dy + dz*dz;
  }

  template<typename T, template<typename U> class Descriptor>
  void SetSpheres3D<T,Descriptor>::process(Box3D domain, BlockLattice3D<T,Descriptor> &lattice)
  {
    
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
      for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
        for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
          Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
	
          Dot3D relativePosition = lattice.getLocation();
          T xGlobal = (T) (relativePosition.x + iX);
          T yGlobal = (T) (relativePosition.y + iY);
          T zGlobal = (T) (relativePosition.z + iZ);
          
          T *sfPtr = cell.getExternal(Descriptor<T>::ExternalField::volumeFractionBeginsAt);
          T *idPtr = cell.getExternal(Descriptor<T>::ExternalField::particleIdBeginsAt);
          T *uPtr = cell.getExternal(Descriptor<T>::ExternalField::boundaryVelocityBeginsAt);

          *idPtr = (T)-1;
          uPtr[0] = 0;
          uPtr[1] = 0;
          uPtr[2] = 0;
          *sfPtr = 0;

          for(plint iP=0;iP<nSpheres;iP++){
          
            T dx = xGlobal - x[iP][0];
            T dy = yGlobal - x[iP][1];
            T dz = zGlobal - x[iP][2];
	
           
            distSqr[iP] = dx*dx+dy*dy+dz*dz;
            if(dx*dx+dy*dy+dz*dz > (r[iP]+1)*(r[iP]+1))
              continue;
       
            T sf = calcSolidFraction(dx,dy,dz,r[iP]);

            if(sf > SOLFRAC_MIN){
              plint ind(iP);
              if(*sfPtr > SOLFRAC_MIN && distSqr[iP] > distSqr[(plint)*idPtr])
                ind = *idPtr;
              
              *idPtr = (T)id[ind];
              *sfPtr = sf;
              uPtr[0] = v[ind][0] + omega[ind][1]*dz - omega[ind][2]*dy;
              uPtr[1] = v[ind][1] - omega[ind][0]*dz + omega[ind][2]*dx; 
              uPtr[2] = v[ind][2] + omega[ind][0]*dy - omega[ind][1]*dx;
            } 
          }
        }
      }
    }
    
  }

  template<typename T, template<typename U> class Descriptor>
  T SetSpheres3D<T,Descriptor>::calcSolidFraction(T dx_, T dy_, T dz_, T r_)
  {
    plint const slicesPerDim = 5;
    
    T sliceWidth = 1./((T)slicesPerDim-1);
    T fraction = 1./((T)slicesPerDim*slicesPerDim*slicesPerDim);
    

    if (dx_*dx_ + dy_*dy_ + dz_*dz_ > (r_+1)*(r_+1))
      return 0;

    if (dx_*dx_ + dy_*dy_ + dz_*dz_ < (r_-1)*(r_-1))
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
  SetSpheres3D<T,Descriptor>* SetSpheres3D<T,Descriptor>::clone() const
  { 
    return new SetSpheres3D<T,Descriptor>(*this);
  }

  template<typename T, template<typename U> class Descriptor>
  void SetSpheres3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
  {
    modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
  }

  /*
   * implementation of SumForceTorque3D
   */

  /* --------------------------------------------- */
  
  template<typename T, template<typename U> class Descriptor>
  void SumForceTorque3D<T,Descriptor>::process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice)
  {
    static plint pid = Descriptor<T>::ExternalField::particleIdBeginsAt,
      fx = Descriptor<T>::ExternalField::hydrodynamicForceBeginsAt,
      fy = Descriptor<T>::ExternalField::hydrodynamicForceBeginsAt+1,
      fz = Descriptor<T>::ExternalField::hydrodynamicForceBeginsAt+2,
      solfrac_id = Descriptor<T>::ExternalField::volumeFractionBeginsAt;
    
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
      for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
        for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
          Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
          
          plint id = (plint) (*(cell.getExternal(pid)));
          if(id < 0) continue; // no particle here
          
          Dot3D relativePosition = lattice.getLocation();
          T xGlobal = (T) (relativePosition.x + iX);
          T yGlobal = (T) (relativePosition.y + iY);
          T zGlobal = (T) (relativePosition.z + iZ);
          
          T fs = *(cell.getExternal(solfrac_id));

          T forceX = fs*(*(cell.getExternal(fx)));
          T forceY = fs*(*(cell.getExternal(fy)));
          T forceZ = fs*(*(cell.getExternal(fz)));
          
          
          force[id][0] += forceX;
          force[id][1] += forceY;
          force[id][2] += forceZ;
          
          T dx = xGlobal - x[id][0];
          T dy = yGlobal - x[id][1];
          T dz = zGlobal - x[id][2];

          torque[id][0] += dy*forceZ - dz*forceY;
          torque[id][1] += -dx*forceZ + dz*forceX;
          torque[id][2] += dx*forceY - dy*forceX;
        }
      }
    }
    
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
};
