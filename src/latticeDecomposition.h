#ifndef LATTICE_DECOMPOSITION_H
#define LATTICE_DECOMPOSITION_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include "lammps.h"

namespace plb{

class LatticeDecomposition {
public:
  LatticeDecomposition(plb::plint nx_, plb::plint ny_, plb::plint nz_, LAMMPS_NS::LAMMPS *lmp_);
  
  ~LatticeDecomposition();

  plb::SparseBlockStructure3D getBlockDistribution();
  plb::ExplicitThreadAttribution* getThreadAttribution();
private:
  plb::plint nx,ny,nz;
  LAMMPS_NS::LAMMPS &lmp;
  plb::plint npx,npy,npz;
  std::vector<plb::plint> xVal, yVal, zVal;
  plb::SparseBlockStructure3D *blockStructure;
  plb::ExplicitThreadAttribution *threadAttribution;
};

}

#endif /* LATTICE_DECOMPOSITION_H */
