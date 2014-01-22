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

    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor> &lattice);
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const
    {
      modified[0] = modif::nothing;
      modified[1] = modif::staticVariables;
    }
    SetSingleSphere3D<T,Descriptor>* clone() const { return new SetSingleSphere3D<T,Descriptor>(*this);}
  private:
    T *x,*v,*omega,r;
    int id;
    T calcSolidFraction(T dx_, T dy_, T dz_, T r_);
    void setValues(Cell<T,Descriptor>& c, T sf, T dx, T dy, T dz);
    void setToZero(Cell<T,Descriptor>& c);
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
