#ifndef IBDATAEXCHANGEWRAPPERS_HH_LBDEM
#define IBDATAEXCHANGEWRAPPERS_HH_LBDEM


namespace plb {
  template<typename T, template<typename U> class Descriptor>
  void setSpheresOnLattice(MultiBlockLattice3D<T,Descriptor> &lattice,
                           LiggghtsCouplingWrapper &wrapper,
                           PhysUnits3D<T> const &units,
                           std::vector<plint> &excludeType,
                           bool initVelFlag)
  {
    plint rank = global::mpi().getRank();

    // this relies on the fact that there is exactly one block on each lattice
    plint iBlock = lattice.getLocalInfo().getBlocks()[0];
    std::map<plint,Box3D> blockmap = lattice.getSparseBlockStructure().getBulks();
    Box3D localBB = blockmap[iBlock];

    plint nx=lattice.getNx(), ny=lattice.getNy(), nz=lattice.getNz();
    plint nPart = wrapper.lmp->atom->nlocal + wrapper.lmp->atom->nghost;


    
    FixMultisphere *fixMultisphere = 0;
    plint nms = wrapper.lmp->modify->n_fixes_style("multisphere");
    if(nms)
      fixMultisphere = static_cast<FixMultisphere*>(wrapper.lmp->modify->find_fix_style("multisphere",0));

    for(plint iS=0;iS<nPart;iS++){
      plint type = (plint)wrapper.lmp->atom->type[iS];
      bool excludeFlag(false);
      for(plint iT=0;iT<excludeType.size();iT++){
        //        pcout << iS << " " << type << " " << excludeType[iT] << std::endl;
        if(type == excludeType[iT]){
          excludeFlag = true;
          break;
        }
      }

      if(excludeFlag) continue;

      T x[3],v[3],omega[3];
      T r;
      plint id = wrapper.lmp->atom->tag[iS];

      SetSingleSphere3D<T,Descriptor> *sss = 0;

      for(plint i=0;i<3;i++){
        x[i] = units.getLbLength(wrapper.lmp->atom->x[iS][i]);
      }
      r = units.getLbLength(wrapper.lmp->atom->radius[iS]);

      plint msBody = fixMultisphere ? fixMultisphere->belongs_to(id) : -1;

      // std::cout << "proc " << rank << " | " << "id " << id << " multi-body " << msBody << std::endl;
      


      // if particle belongs to multisphere, don't use omega
      if(msBody > -1){
        T x_com[3],v_com[3];
	fixMultisphere->data().xcm(x_com,msBody);
	fixMultisphere->data().vcm(v_com,msBody);
        for(plint i=0;i<3;i++){
          v[i] = units.getLbVel(wrapper.lmp->atom->v[iS][i]);
	  omega[i] = 0;//units.getLbFreq(wrapper.lmp->atom->omega[iS][i]);
          x_com[i] = units.getLbLength(x_com[i]);
          v_com[i] = units.getLbVel(v_com[i]);
        }
	// std::cout << "proc " << rank << " | " 
	// 	  << "x " << x[0] << " " << x[1] << " " << x[2] << " | "
	// 	  << "v " << v[0] << " " << v[1] << " " << v[2] << std::endl;
	// std::cout << "proc " << rank << " | " 
	// 	  << "x_com " << x_com[0] << " " << x_com[1] << " " << x_com[2] << " | "
	// 	  << "v_com " << v_com[0] << " " << v_com[1] << " " << " " << v_com[2] << std::endl;

	sss = new SetSingleSphere3D<T,Descriptor>(x,v,omega,x,r,id,initVelFlag);

      } else {
        for(plint i=0;i<3;i++){
          v[i] = units.getLbVel(wrapper.lmp->atom->v[iS][i]);
	  omega[i] = units.getLbFreq(wrapper.lmp->atom->omega[iS][i]);
        }

        // use sphere center as center of mass for simple spheres
        sss  = new SetSingleSphere3D<T,Descriptor>(x,v,omega,x,r,id,initVelFlag);
      }

      Box3D sss_box = sss->getBoundingBox();
      
      // only go over part that lies on local processor
      // to avoid unnecessary communication overhead
      Box3D sss_box_intersect(0,0,0,0,0,0);
      bool boxes_intersect = intersect(sss_box,localBB,sss_box_intersect);
      // std::cout << r << " intersect " << boxes_intersect << " "
      //           << sss_box_intersect.x0 << " " << sss_box_intersect.x1 << " "
      //           << sss_box_intersect.y0 << " " << sss_box_intersect.y1 << " "
      //           << sss_box_intersect.z0 << " " << sss_box_intersect.z1 << " "
      //           << std::endl;
      if(boxes_intersect)
        applyProcessingFunctional(sss,sss_box_intersect,lattice);
      else
        delete sss;
      // applyProcessingFunctional(sss,sss_box,lattice);
    }

    // this one returns modif::staticVariables and forces an update of those along processor
    // boundaries
    applyProcessingFunctional(new AttributeFunctional<T,Descriptor>(),lattice.getBoundingBox(),lattice);

  }

  template<typename T, template<typename U> class Descriptor>
  void setSpheresOnLattice(MultiBlockLattice3D<T,Descriptor> &lattice,
                           LiggghtsCouplingWrapper &wrapper,
                           PhysUnits3D<T> const &units,
                           bool initVelFlag)
  {
    std::vector<plint> dummyExcludeType;
    setSpheresOnLattice(lattice,wrapper,units,dummyExcludeType,initVelFlag);
  }



  template<typename T, template<typename U> class Descriptor>
  void getForcesFromLattice(MultiBlockLattice3D<T,Descriptor> &lattice,
                            LiggghtsCouplingWrapper &wrapper,
                            PhysUnits3D<T> const &units)
  {
    // debug stuff
    plint r = global::mpi().getRank();

    static std::vector<T> force,torque;
    static typename ParticleData<T>::ParticleDataArrayVector x_lb;

    plint const nPart = wrapper.lmp->atom->nlocal + wrapper.lmp->atom->nghost;
    plint const n_force = nPart*3;

    if(nPart == 0) return; // no particles - no work

    // std::cout << "proc " << r << " | ";
    

    if(nPart > x_lb.size()){
      for(plint iPart=0;iPart<x_lb.size();iPart++){
        x_lb[iPart][0] = units.getLbLength(wrapper.lmp->atom->x[iPart][0]);
        x_lb[iPart][1] = units.getLbLength(wrapper.lmp->atom->x[iPart][1]);
        x_lb[iPart][2] = units.getLbLength(wrapper.lmp->atom->x[iPart][2]);
      }
      for(plint iPart = x_lb.size();iPart < nPart; iPart++)
        x_lb.push_back( Array<T,3>( units.getLbLength(wrapper.lmp->atom->x[iPart][0]),
                                    units.getLbLength(wrapper.lmp->atom->x[iPart][1]),
                                    units.getLbLength(wrapper.lmp->atom->x[iPart][2]) ) );
    } else{
      for(plint iPart=0;iPart<nPart;iPart++){
        x_lb[iPart][0] = units.getLbLength(wrapper.lmp->atom->x[iPart][0]);
        x_lb[iPart][1] = units.getLbLength(wrapper.lmp->atom->x[iPart][1]);
        x_lb[iPart][2] = units.getLbLength(wrapper.lmp->atom->x[iPart][2]);
      }
    }

    if(n_force > force.size()){
      for(plint i=0;i<force.size();i++){
        force[i] = 0;
        torque[i] = 0;
      }
      for(plint i=force.size();i<n_force;i++){
        force.push_back(0.);
        torque.push_back(0.);
      }
    } else {
      for(plint i=0;i<n_force;i++){
        force[i] = 0;
        torque[i] = 0;
      }
    }

    // std::cout << "proc " << r << " | " << "allocation done" << std::endl;;    

    SumForceTorque3D<T,Descriptor> *sft = new SumForceTorque3D<T,Descriptor>(x_lb,
                                                                             &force.front(),&torque.front(),
                                                                             wrapper);
    
    // this relies on the fact that there is exactly one block on each processor
    plint iBlock = lattice.getLocalInfo().getBlocks()[0];
    std::map<plint,Box3D> blockmap = lattice.getSparseBlockStructure().getBulks();
    Box3D localBB = blockmap[iBlock];
    applyProcessingFunctional(sft,localBB, lattice);

    // for(plint i=0;i<force.size();i++)
    //   pcerr << force[i] << " ";
    // pcerr << "|";
    // for(plint i=0;i<torque.size();i++)
    //   pcerr << torque[i] << " ";
    // pcerr << std::endl;

    // // experimental....
    // SparseBlockStructure3D sparseBlock = lattice.getSparseBlockStructure();
    // std::vector<Box3D> boxVec;
    // std::vector<plint> localBlocks = lattice.getLocalInfo().getBlocks();
    // sparseBlock.intersect(lattice.getBoundingBox(),localBlocks,boxVec);
    // std::cout << r << " boxvec " << boxVec.size() << std::endl;
    // for(plint i=0;i<boxVec.size();i++)
    //   applyProcessingFunctional(sft,boxVec[i],lattice);
    
    // // basic version with complete lattice
    // applyProcessingFunctional(sft,lattice.getBoundingBox(), lattice);

    LAMMPS_NS::FixLbCouplingOnetoone 
      *couplingFix 
      = dynamic_cast<LAMMPS_NS::FixLbCouplingOnetoone*>
      (wrapper.lmp->modify->find_fix_style("couple/lb/onetoone",0));

    double **f_liggghts = couplingFix->get_force_ptr();
    double **t_liggghts = couplingFix->get_torque_ptr();

    for(plint iPart=0;iPart<nPart;iPart++){
      int tag = wrapper.lmp->atom->tag[iPart];
      int liggghts_ind = wrapper.lmp->atom->map(tag);
      for(plint i=0;i<3;i++){
        f_liggghts[liggghts_ind][i] = units.getPhysForce(force[3*iPart+i]);
        t_liggghts[liggghts_ind][i] = units.getPhysTorque(torque[3*iPart+i]);
      }
    }
    couplingFix->comm_force_torque();

  }

}; /* namespace plb */

#endif /* IBDATAEXCHANGEWRAPPERS_HH_LBDEM */
