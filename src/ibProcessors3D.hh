/*
 * ibProcessors3D.hh
 */

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
   * implementation of SetSphere3D
   */

  template<typename T, template<typename U> class Descriptor>
  void SetSphere3D<T,Descriptor>::process(Box3D domain, BlockLattice3D<T,Descriptor> &lattice)
  {
    static plint
      partId_ind = Descriptor<T>::ExternalField::particleIdBeginsAt,
      ux_ind = Descriptor<T>::ExternalField::boundaryVelocityBeginsAt,
      uy_ind = Descriptor<T>::ExternalField::boundaryVelocityBeginsAt+1,
      uz_ind = Descriptor<T>::ExternalField::boundaryVelocityBeginsAt+2;
    
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
      for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
        for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
          Cell<T,Descriptor>& cell = lattice.get(iX,iY,iZ);
	
          Dot3D relativePosition = lattice.getLocation();
          T xGlobal = (T) (relativePosition.x + iX);
          T yGlobal = (T) (relativePosition.y + iY);
          T zGlobal = (T) (relativePosition.z + iZ);
          
          T dx = xGlobal - x[0];
          T dy = yGlobal - x[1];
          T dz = zGlobal - x[2];
	
          T solfrac =  calcSolidFraction(dx,dy,dz);
	
          *(cell.getExternal(Descriptor<T>::ExternalField::volumeFractionBeginsAt)) = solfrac;
	
          if(solfrac > 0.001){
            *(cell.getExternal(partId_ind)) = (T)id;
            *(cell.getExternal(ux_ind)) = v[0] + dy*omega[2] - dz*omega[1];
            *(cell.getExternal(uy_ind)) = v[1] - dx*omega[2] + dz*omega[0];
            *(cell.getExternal(uz_ind)) = v[2] + dx*omega[1] - dy*omega[0];
          }
          else {
            *(cell.getExternal(partId_ind)) = (T)-1;
            *(cell.getExternal(ux_ind)) = 0.;
            *(cell.getExternal(uy_ind)) = 0.;
            *(cell.getExternal(uz_ind)) = 0.;
          }
        }
      }
    }
    
  }

  template<typename T, template<typename U> class Descriptor>
  T SetSphere3D<T,Descriptor>::calcSolidFraction(T dx_, T dy_, T dz_)
  {
    plint const slicesPerDim = 5;
    
    T sliceWidth = 1./((T)slicesPerDim-1);
    T fraction = 1./((T)slicesPerDim*slicesPerDim*slicesPerDim);
    
    T r_sq = r*r;

    if (dx_*dx_ + dy_*dy_ + dz_*dz_ > (r+1)*(r+1))
      return 0;

    if (dx_*dx_ + dy_*dy_ + dz_*dz_ < (r-1)*(r-1))
      return 1;

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
  SetSphere3D<T,Descriptor>* SetSphere3D<T,Descriptor>::clone() const
  { 
    return new SetSphere3D<T,Descriptor>(*this);
  }

  template<typename T, template<typename U> class Descriptor>
  void SetSphere3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
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
