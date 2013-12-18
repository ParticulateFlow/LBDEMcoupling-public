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
                         ScalarField3D<T2>& field);
    GetExternalScalarFunctional3D<T1,Descriptor,T2>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
  private:
    int scalar;
  }; 


  template<typename T, template<typename U> class Descriptor>
  struct SetSpheres3D : public BoxProcessingFunctional3D_L<T,Descriptor> {
  public:
    SetSpheres3D(T **x_, T **v_, T *r_, int **id_, plint nSpheres_)
      : x(x_),v(v_),omega(0),r(r_),id(id_),nSpheres(nSpheres_) 
    {
      distSqr = new T[nSpheres];
    }
    SetSpheres3D(T **x_, T **v_, T **omega_, T *r_, int **id_, plint nSpheres_)
      : x(x_),v(v_),omega(omega_),r(r_),id(id_),nSpheres(nSpheres_) 
    {
      distSqr = new T[nSpheres];
    }
    SetSpheres3D(const SetSpheres3D &orig)
      : x(orig.x),v(orig.v),omega(orig.omega),r(orig.r),id(orig.id),nSpheres(orig.nSpheres)
    {
      distSqr = new T[nSpheres];
      for(plint i=0;i<nSpheres;i++)
        distSqr[i] = orig.distSqr[i];
    }

    ~SetSpheres3D()
    {
      delete[] distSqr;
    }
    void process(Box3D domain, BlockLattice3D<T,Descriptor> &lattice);

    SetSpheres3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
  private:
    T calcSolidFraction(T dx_, T dy_, T dz_, T r_);
    T calcDistSqr(T x0, T y0, T z0, T x1, T y1, T z1);
    T **x,**v,**omega;
    T *r;
    int **id;
    plint nSpheres;
    T *distSqr;
    
  };

  template<typename T, template<typename U> class Descriptor>
  struct SetSingleSphere3D : public BoxProcessingFunctional3D_L<T,Descriptor> {
  public:
    SetSingleSphere3D(T *x_, T *v_, T *omega_, T r_, int id_)
      : x(x_),v(v_),omega(omega_),r(r_),id(id_) {}
    SetSingleSphere3D(T *x_, T *v_, T r_, int id_)
      : x(x_),v(v_),omega(0),r(r_),id(id_) {}
    SetSingleSphere3D(const SetSingleSphere3D &orig)
      : x(orig.x),v(orig.v),omega(orig.omega),r(orig.r),id(orig.id) {}

    Box3D getBoundingBox() const { return Box3D(x[0]-r-2,x[0]+r+2,
                                                x[1]-r-2,x[1]+r+2,
                                                x[2]-r-2,x[2]+r+2);}

    void process(Box3D domain, BlockLattice3D<T,Descriptor> &lattice)
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
              if(((plint)*idPtr) == id)
                setToZero(cell);
              // else do nothing
              break;
            case 3:
              if(sf > *sfPtr || ((plint)*idPtr) == id)
                setValues(cell,sf,dx,dy,dz);
              // else do nothing
            }
          }
        }
      }
    }
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const
    {
      modified[0] = modif::nothing;
      modified[1] = modif::staticVariables;
    }
    SetSingleSphere3D<T,Descriptor>* clone() const { return new SetSingleSphere3D<T,Descriptor>(*this);}
  private:
    T *x,*v,*omega,r;
    int id;
    T calcSolidFraction(T dx_, T dy_, T dz_, T r_)
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
      return solFrac;//solFrac>0.5 ? 1. : 0.;
    }
    void setValues(Cell<T,Descriptor>& c, T sf, T dx, T dy, T dz)
    {
      T *sfPtr = c.getExternal(Descriptor<T>::ExternalField::volumeFractionBeginsAt);
      T *idPtr = c.getExternal(Descriptor<T>::ExternalField::particleIdBeginsAt);
      T *uPtr = c.getExternal(Descriptor<T>::ExternalField::boundaryVelocityBeginsAt);
      
      uPtr[0] = v[0];
      uPtr[1] = v[1];
      uPtr[2] = v[2];
      if(omega != 0){
        uPtr[0] += omega[1]*dz - omega[2]*dy;
        uPtr[1] += omega[0]*dz + omega[2]*dx; 
        uPtr[2] += omega[0]*dy - omega[1]*dx;
      }
      *sfPtr = sf;
      *idPtr = (T) id;
    }
    void setToZero(Cell<T,Descriptor>& c)
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
    
  };
  
  template<typename T, template<typename U> class Descriptor>
  void setSpheresOnLattice(MultiBlockLattice3D<T,Descriptor> &lattice,
                           T **x, T **v, T **omega, T *r, int **id, plint nSpheres)
  {
    for(plint iS=0;iS<nSpheres;iS++){
      SetSingleSphere3D<T,Descriptor> *sss 
        = new SetSingleSphere3D<T,Descriptor>(x[iS],v[iS],omega==0 ? 0 : omega[iS],r[iS],id[0][iS]);
      applyProcessingFunctional(sss,sss->getBoundingBox(),lattice);
      //delete sss;
    }
  }

  
  template<typename T, template<typename U> class Descriptor>
  struct SumForceTorque3D : public ReductiveBoxProcessingFunctional3D_L<T,Descriptor> {
  public:
    SumForceTorque3D(plint nPart_, T **x_);
    // SumForceTorque3D(plint const *partId_, T **x_, T **force_, T **torque_)
    //   : partId(partId_), x(x_), force(force_),torque(torque_) {}
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice);

    SumForceTorque3D<T,Descriptor>* clone() const;
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
    std::vector<double>& getForceTorque(){return this->getStatistics().getSumVect();}
  private:
    std::vector<plint> sumId;
    plint const *partId;
    T **x;
    T **force, **torque;
    void addForce(plint partId, plint coord, T value);
    void addTorque(plint partId, plint coord, T value);
  };

}; // plb

#include "ibProcessors3D.hh"

#endif // IB_PROCESSORS_3D_H
