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

#ifndef LIGGGHTSCOUPLINGWRAPPER_H
#define LIGGGHTSCOUPLINGWRAPPER_H

#include "palabos3D.h"
#include "palabos3D.hh"

// necessary LAMMPS/LIGGGHTS includes

#include "lammps.h"
#include "input.h"
#include "library.h"
#include "library_cfd_coupling.h"

class LiggghtsCouplingWrapper {
public:
  LiggghtsCouplingWrapper(char **argv, MPI_Comm communicator);
  void execFile(char* const fname);
  void execCommand(std::stringstream const &cmd);
  void execCommand(char* const cmd);
  void run(plb::plint nSteps);
  void runUpto(plb::plint nSteps);
  int getNumParticles();
  void setVariable(char const *name, double value);
  void setVariable(char const *name, std::string &value);

  //private:
  LAMMPS_NS::LAMMPS *lmp;
};

#endif /* LIGGGHTSCOUPLINGWRAPPER_H */
