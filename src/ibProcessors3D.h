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
    SetSpheres3D(T **x_, T **v_, T **omega_, T *r_, plint *id_, plint nSpheres_)
      : x(x_),v(v_),omega(omega_),r(r_),id(id_),nSpheres(nSpheres_) 
    {
      distSqr = new T[nSpheres];
    }
    SetSpheres3D(const SetSpheres3D &orig)
      : x(orig.x),v(orig.v),omega(orig.omega),r(orig.r),id(orig.id),nSpheres(orig.nSpheres)
    {
      distSqr = new T[nSpheres];
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
    plint *id, nSpheres;
    T *distSqr;
    
  };
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
