/*
 * This file is part of the LBDEMcoupling software.
 *
 * LBDEMcoupling is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, version 3.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Copyright 2014 Johannes Kepler University Linz
 *
 * Author: Philippe Seil (philippe.seil@jku.at)
 */

/*
 * helper class to create lattice decomposition from LIGGGHTS parallelization
 */

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
