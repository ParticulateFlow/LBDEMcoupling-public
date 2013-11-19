/*
 * ibProcessors3D.h
 */

#ifndef IB_PROCESSORS_3D_H
#define IB_PROCESSORS_3D_H

namespace plb{
  template<typename T1, template<typename U> class Descriptor, typename T2>
  struct GetExternalScalarFunctional3D : public BoxProcessingFunctional3D_LS<T1,Descriptor,T2> {
    GetExternalScalarFunctional3D(int const scalar_)
    : BoxProcessingFunctional3D_LS<T1,Descriptor,T2>(), scalar(scalar_) {}

    virtual void process(Box3D domain, BlockLattice3D<T1,Descriptor>& lattice,
                                       ScalarField3D<T2>& field){
      for(plint i=domain.x0;i<=domain.x1;i++)
        for(plint j=domain.y0;j<=domain.y1;j++)
          for(plint k=domain.z0;k<=domain.z1;k++){
            field.get(i,j,k) = *(lattice.get(i,j,k).getExternal(scalar));
          }
    }
    GetExternalScalarFunctional3D* clone() const{
      return new GetExternalScalarFunctional3D(*this);
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const
    {
      modified[0] = modif::nothing;
      modified[1] = modif::staticVariables;
    }
  private:
    int scalar;
  }; 

  template<typename T, template<typename U> class Descriptor>
  struct SetSphere3D : public BoxProcessingFunctional3D_L<T,Descriptor> {
  public:
    SetSphere3D(Array<T,3> x_, Array<T,3> v_, Array<T,3> omega_, T r_, plint id_)
      : x(x_),v(v_),omega(omega_),r(r_),id(id_) {}

    SetSphere3D* clone() const
    {
      return new SetSphere3D(*this);
    }
    void process(Box3D domain, BlockLattice3D<T,Descriptor> &lattice);
    T calcSolidFraction(T dx_, T dy_, T dz_);

    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const
    {
      modified[0] = modif::nothing;
      modified[1] = modif::staticVariables;
    }

  private:
    Array<T,3> x,v,omega;
    T r;
    plint id;
    
  };

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
    pluint const slicesPerDim = 5;
    
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

};

#endif // IB_PROCESSORS_3D_H
