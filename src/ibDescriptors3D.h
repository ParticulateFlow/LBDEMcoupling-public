/*
 * ibDescriptors3D.h
 *
 **/

#ifndef IB_DESCRIPTORS_3D_H
#define IB_DESCRIPTORS_3D_H

namespace plb {
  namespace descriptors {
    struct IBDynamicsDescriptor3D {
      static const int volumeFractionBeginsAt = 0;
      static const int sizeOfVolumeFraction = 1;

      static const int particleIdBeginsAt = volumeFractionBeginsAt + sizeOfVolumeFraction;
      static const int sizeOfParticleId = 1;

      static const int boundaryVelocityBeginsAt = particleIdBeginsAt + sizeOfParticleId;
      static const int sizeOfBoundaryVelocity = 3;

      static const int hydrodynamicForceBeginsAt = boundaryVelocityBeginsAt + sizeOfBoundaryVelocity;
      static const int sizeOfHydrodynamicForce = 3;

      static const int bBeginsAt = hydrodynamicForceBeginsAt + sizeOfHydrodynamicForce;
      static const int sizeOfB = 1;

      static const int forceBeginsAt = 0;
      static const int sizeOfForce   = 0;

      static const int numScalars = sizeOfVolumeFraction + sizeOfParticleId
        + sizeOfBoundaryVelocity + sizeOfHydrodynamicForce + sizeOfB + sizeOfForce;
      static const int numSpecies = 5;
    };

    struct IBDynamicsDescriptorBase3D{
      typedef IBDynamicsDescriptor3D ExternalField;
    };

    // should describe a cell/dynamics necessary for immersed boundary method
    template <typename T> struct ImmersedBoundaryD3Q19Descriptor
        : public D3Q19DescriptorBase<T>, public IBDynamicsDescriptorBase3D
    {
      static const char name [];
    };

    template<typename T>
    const char ImmersedBoundaryD3Q19Descriptor<T>::name[] = "ImmersedBoundaryD3Q19";


  }; // descriptors
}; // plb

#endif // IB_DESCRIPTORS_3D_H
