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
 * implementation of data exchange wrapper functions
 */

#ifndef IBDATAEXCHANGEWRAPPERS_HH_LBDEM
#define IBDATAEXCHANGEWRAPPERS_HH_LBDEM

#ifdef LBDEM_USE_MULTISPHERE
#include "fix_multisphere.h"
#endif

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


      // LIGGGHTS-public does not support multisphere (yet?)
      // thus need to explicitly exclude multisphere method
#ifdef LBDEM_USE_MULTISPHERE
    
    FixMultisphere *fixMultisphere = 0;
    plint nms = wrapper.lmp->modify->n_fixes_style("multisphere");
    if(nms)
      fixMultisphere = static_cast<FixMultisphere*>(wrapper.lmp->modify->find_fix_style("multisphere",0));

#endif /* LBDEM_USE_MULTISPHERE */

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
        x[i] = units.getLbPosition(wrapper.lmp->atom->x[iS][i]);
      }
      r = units.getLbLength(wrapper.lmp->atom->radius[iS]);

      // LIGGGHTS-public does not support multisphere (yet?)
      // thus need to explicitly exclude multisphere method
#ifdef LBDEM_USE_MULTISPHERE

      plint msBody = fixMultisphere ? fixMultisphere->belongs_to(id) : -1;


      // if particle belongs to multisphere, don't use omega
      if(msBody > -1){
        T x_com[3],v_com[3];
	fixMultisphere->data().xcm(x_com,msBody);
	fixMultisphere->data().vcm(v_com,msBody);
        for(plint i=0;i<3;i++){
          v[i] = units.getLbVel(wrapper.lmp->atom->v[iS][i]);
	  omega[i] = 0;//units.getLbFreq(wrapper.lmp->atom->omega[iS][i]);
          x_com[i] = units.getLbPosition(x_com[i]);
          v_com[i] = units.getLbVel(v_com[i]);
        }

	sss = new SetSingleSphere3D<T,Descriptor>(x,v,omega,x,r,id,initVelFlag);

      } else {
        for(plint i=0;i<3;i++){
          v[i] = units.getLbVel(wrapper.lmp->atom->v[iS][i]);
	  omega[i] = units.getLbFreq(wrapper.lmp->atom->omega[iS][i]);
        }

        // use sphere center as center of mass for simple spheres
        sss  = new SetSingleSphere3D<T,Descriptor>(x,v,omega,x,r,id,initVelFlag);
      }

#else
      for(plint i=0;i<3;i++){
        v[i] = units.getLbVel(wrapper.lmp->atom->v[iS][i]);
        omega[i] = units.getLbFreq(wrapper.lmp->atom->omega[iS][i]);
      }
      r = units.getLbLength(wrapper.lmp->atom->radius[iS]);
      sss = new SetSingleSphere3D<T,Descriptor>(x,v,omega,x,r,id,initVelFlag);

#endif /* LBDEM_USE_MULTISPHERE */

      Box3D sss_box = sss->getBoundingBox();

      // only go over part that lies on local processor
      // to avoid unnecessary communication overhead
      Box3D sss_box_intersect(0,0,0,0,0,0);

      bool boxes_intersect = intersect(sss_box,localBB,sss_box_intersect);

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

#ifdef PLB_DEBUG
    pcerr << "entering getForcesFromLattice" << std::endl;
#endif
    
    static std::vector<T> force,torque;
    static typename ParticleData<T>::ParticleDataArrayVector x_lb;

    plint const nPart = wrapper.lmp->atom->nlocal + wrapper.lmp->atom->nghost;
    plint const n_force = nPart*3;

#ifdef PLB_DEBUG
    pcerr << "number of particles: " << nPart << std::endl;
#endif

    if(nPart == 0) return; // no particles - no work

    if(nPart > x_lb.size()){
      for(plint iPart=0;iPart<x_lb.size();iPart++){
        x_lb[iPart][0] = units.getLbPosition(wrapper.lmp->atom->x[iPart][0]);
        x_lb[iPart][1] = units.getLbPosition(wrapper.lmp->atom->x[iPart][1]);
        x_lb[iPart][2] = units.getLbPosition(wrapper.lmp->atom->x[iPart][2]);
      }
      for(plint iPart = x_lb.size();iPart < nPart; iPart++)
        x_lb.push_back( Array<T,3>( units.getLbPosition(wrapper.lmp->atom->x[iPart][0]),
                                    units.getLbPosition(wrapper.lmp->atom->x[iPart][1]),
                                    units.getLbPosition(wrapper.lmp->atom->x[iPart][2]) ) );
    } else{
      for(plint iPart=0;iPart<nPart;iPart++){
        x_lb[iPart][0] = units.getLbPosition(wrapper.lmp->atom->x[iPart][0]);
        x_lb[iPart][1] = units.getLbPosition(wrapper.lmp->atom->x[iPart][1]);
        x_lb[iPart][2] = units.getLbPosition(wrapper.lmp->atom->x[iPart][2]);
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

#ifdef PLB_DEBUG
    pcerr << "starting SumForceTorque" << std::endl;
#endif
    
    SumForceTorque3D<T,Descriptor> *sft = new SumForceTorque3D<T,Descriptor>(x_lb,
                                                                             &force.front(),&torque.front(),
                                                                             wrapper);
    
    // this relies on the fact that there is exactly one block on each processor
    plint iBlock = lattice.getLocalInfo().getBlocks()[0];
    std::map<plint,Box3D> blockmap = lattice.getSparseBlockStructure().getBulks();
    Box3D localBB = blockmap[iBlock];
    applyProcessingFunctional(sft,localBB, lattice);


    LAMMPS_NS::FixLbCouplingOnetoone 
      *couplingFix 
      = dynamic_cast<LAMMPS_NS::FixLbCouplingOnetoone*>
      (wrapper.lmp->modify->find_fix_style("couple/lb/onetoone",0));

    double **f_liggghts = couplingFix->get_force_ptr();
    double **t_liggghts = couplingFix->get_torque_ptr();

    for(plint iPart=0;iPart<nPart;iPart++)
      for(plint j=0;j<3;j++){
        f_liggghts[iPart][j] = 0;
        t_liggghts[iPart][j] = 0;
      }


    for(plint iPart=0;iPart<nPart;iPart++){
      int tag = wrapper.lmp->atom->tag[iPart];
      int liggghts_ind = wrapper.lmp->atom->map(tag);

      for(plint j=0;j<3;j++){
        f_liggghts[liggghts_ind][j] += units.getPhysForce(force[3*iPart+j]);
        t_liggghts[liggghts_ind][j] += units.getPhysTorque(torque[3*iPart+j]);
      }
    }
    couplingFix->comm_force_torque();
  }

}; /* namespace plb */

#endif /* IBDATAEXCHANGEWRAPPERS_HH_LBDEM */
