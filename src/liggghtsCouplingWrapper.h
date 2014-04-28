
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
  void allocateVariables();
  void dataFromLiggghts();
  void dataToLiggghts();
  void run(plb::plint nSteps);
  void runUpto(plb::plint nSteps);
  plb::plint getNumParticles();
  void setVariable(char const *name, double value);
  void setVariable(char const *name, std::string &value);

  //private:
  LAMMPS_NS::LAMMPS *lmp;
  double **x, **v, **f, **t, **r, **omega;
  int **id;

};

#endif /* LIGGGHTSCOUPLINGWRAPPER_H */
