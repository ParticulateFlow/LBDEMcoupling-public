class LiggghtsCouplingWrapper {
public:
  LiggghtsCouplingWrapper(char **argv, MPI_Comm communicator)
    : lmp(0), x(0),v(0),f(0),t(0),r(0),omega(0),id(0)
  {
    // todo: get LAMMPS to recognize command line options
    int argc_lmp = 1;
    char **argv_lmp = 0;
    argv_lmp = new char*[1];
    argv_lmp[0] = argv[0];

    lmp = new LAMMPS_NS::LAMMPS(argc_lmp,argv_lmp,communicator);

    //    delete[] argv_lmp[0];
    delete[] argv_lmp;
  }
  void execFile(char* const fname)
  {
    lmp->input->file(fname);
  }
  void execCommand(std::stringstream const &cmd)
  {
    lmp->input->one(cmd.str().c_str());
  }
  void execCommand(char* const cmd)
  {
    lmp->input->one(cmd);
  }
  void allocateVariables()
  {
    allocate_external_double(x,3,"nparticles",1.,lmp);
    allocate_external_double(v,3,"nparticles",0.,lmp);
    allocate_external_double(f,3,"nparticles",0.,lmp);
    allocate_external_double(t,3,"nparticles",0.,lmp);
    allocate_external_double(omega,3,"nparticles",0.,lmp);
    allocate_external_double(r,1,"nparticles",0.,lmp);
    allocate_external_int(id,1,"nparticles",0,lmp);
  }
  void dataFromLiggghts()
  {
    data_liggghts_to_of("x","vector-atom",lmp,(void*&)x,"double");
    data_liggghts_to_of("v","vector-atom",lmp,(void*&)v,"double");
    data_liggghts_to_of("omega","vector-atom",lmp,(void*&)omega,"double");
    data_liggghts_to_of("radius","scalar-atom",lmp,(void*&)r,"double");
    data_liggghts_to_of("id","scalar-atom",lmp,(void*&)id,"int");
  }
  void dataToLiggghts()
  {
    data_of_to_liggghts("dragforce","vector-atom",lmp,(void*)f,"double");
    data_of_to_liggghts("hdtorque","vector-atom",lmp,(void*)t,"double");
  }
  plb::plint getNumParticles()
  {
    return lammps_get_natoms(lmp);  
  }
  //private:
  LAMMPS_NS::LAMMPS *lmp;
  double **x, **v, **f, **t, **r, **omega;
  int **id;

};

