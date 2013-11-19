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
  struct SetSphere3D : public BoxProcessingFunctional3D_L<T,Descriptor> {
  public:
    SetSphere3D(Array<T,3> x_, Array<T,3> v_, Array<T,3> omega_, T r_, plint id_)
      : x(x_),v(v_),omega(omega_),r(r_),id(id_) {}

    void process(Box3D domain, BlockLattice3D<T,Descriptor> &lattice);
    T calcSolidFraction(T dx_, T dy_, T dz_);

    SetSphere3D<T,Descriptor>* clone() const;
    virtual void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
  private:
    Array<T,3> x,v,omega;
    T r;
    plint id;
    
  };

  template<typename T, template<typename U> class Descriptor>
  struct SumForceTorque3D : public BoxProcessingFunctional3D_L<T,Descriptor> {
  public:
    SumForceTorque3D(plint const *partId_, T **x_, T **force_, T **torque_)
      : partId(partId_), x(x_), force(force_),torque(torque_) {}
    virtual void process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice);

    SumForceTorque3D<T,Descriptor>* clone() const;
    void getTypeOfModification(std::vector<modif::ModifT>& modified) const;
  private:
    plint const *partId;
    T **x;
    T **force, **torque;
  };

}; // plb

#include "ibProcessors3D.hh"

#endif // IB_PROCESSORS_3D_H
