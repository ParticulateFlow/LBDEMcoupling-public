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

#include "latticeDecomposition.h"
 
// LAMMPS / LIGGGHTS includes
#include "comm.h"

namespace plb {

LatticeDecomposition::LatticeDecomposition(plb::plint nx_, plb::plint ny_, plb::plint nz_, LAMMPS_NS::LAMMPS *lmp_)
  : nx(nx_),ny(ny_),nz(nz_),lmp(*lmp_),
    blockStructure(0),threadAttribution(0)
{
  npx = lmp.comm->procgrid[0];
  npy = lmp.comm->procgrid[1];
  npz = lmp.comm->procgrid[2];

  for(plint i=0;i<=npx;i++)
    xVal.push_back(round( lmp.comm->xsplit[i]*(double)nx ));
  for(plint i=0;i<=npy;i++)
    yVal.push_back(round( lmp.comm->ysplit[i]*(double)ny ));
  for(plint i=0;i<=npz;i++)
    zVal.push_back(round( lmp.comm->zsplit[i]*(double)nz ));

  blockStructure = new SparseBlockStructure3D(Box3D(xVal[0], xVal.back()-1, 
                                                    yVal[0], yVal.back()-1, 
                                                    zVal[0], zVal.back()-1) );
    
  threadAttribution = new ExplicitThreadAttribution();

  for (plint iX=0; iX<xVal.size()-1; ++iX) {
    for (plint iY=0; iY<yVal.size()-1; ++iY) {
      for (plint iZ=0; iZ<zVal.size()-1; ++iZ) {
        plint id = blockStructure->nextIncrementalId();
        blockStructure->addBlock (Box3D( xVal[iX], xVal[iX+1]-1, yVal[iY],
                                         yVal[iY+1]-1, zVal[iZ], zVal[iZ+1]-1 ),
                                  id);
        threadAttribution->addBlock(id,(plint)lmp.comm->grid2proc[iX][iY][iZ]);
      }
    }
  }
    
}
  
LatticeDecomposition::~LatticeDecomposition()
{
  if(blockStructure != 0) delete blockStructure;
  if(threadAttribution != 0) delete threadAttribution;
}


SparseBlockStructure3D LatticeDecomposition::getBlockDistribution()
{
  return SparseBlockStructure3D(*blockStructure);
}
ExplicitThreadAttribution* LatticeDecomposition::getThreadAttribution()
{
  return new ExplicitThreadAttribution(*threadAttribution);
}

}
