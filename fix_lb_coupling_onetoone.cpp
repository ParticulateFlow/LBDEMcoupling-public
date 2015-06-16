/* ----------------------------------------------------------------------
   LIGGGHTS - LAMMPS Improved for General Granular and Granular Heat
   Transfer Simulations

   LIGGGHTS is part of the CFDEMproject
   www.liggghts.com | www.cfdem.com

   Christoph Kloss, christoph.kloss@cfdem.com
   Copyright 2009-2012 JKU Linz
   Copyright 2012-     DCS Computing GmbH, Linz

   LIGGGHTS is based on LAMMPS
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   This software is distributed under the GNU General Public License.

   See the README file in the top-level directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing authors:
   Philippe Seil (JKU Linz)
------------------------------------------------------------------------- */

#include "fix_lb_coupling_onetoone.h"

#include "lammps.h"
#include "modify.h"
#include "error.h"
#include "atom.h"
#include "comm.h"
#include "timer.h"

#include "fix_property_atom.h"

#include <iomanip>
#include <iostream>

namespace LAMMPS_NS {
  FixLbCouplingOnetoone::FixLbCouplingOnetoone(LAMMPS *lmp, int narg, char **arg)
    : Fix(lmp,narg,arg), fix_dragforce_(0), fix_hdtorque_(0), use_torque_(1)
  {

    if(narg < 5) return;

    if(strcmp(arg[3],"use_torque") == 0){
      if(strcmp(arg[4],"yes") == 0){
	use_torque_ = 1;
      } else if(strcmp(arg[4],"no") == 0){
	use_torque_ = 0;
	error->warning(FLERR,"Torque from LB simulation not used");
      } else{
	// complain
	error->all(FLERR,"Illegal fix lbcoupling command: use_torque must be \"yes\" or \"no\"");
      }
    }
  }

  FixLbCouplingOnetoone::~FixLbCouplingOnetoone()
  {

  }

  int FixLbCouplingOnetoone::setmask()
  {
    int mask = 0;
    mask |= FixConst::POST_FORCE;
    mask |= FixConst::POST_RUN;
    return mask;
  }

  void FixLbCouplingOnetoone::init()
  {
    // make sure there is only one fix of this style
    if(modify->n_fixes_style(style) != 1)
      error->fix_error(FLERR,this,"More than one fix of this style is not allowed");


  }

  void FixLbCouplingOnetoone::post_create()
  {
    // register dragforce
    if(!fix_dragforce_)
      {
        const char *fixarg[] = {
          "dragforce", // fix id
          "all",       // fix group
          "property/atom", // fix style: property/atom
          "dragforce",     // property name
          "vector", // 1 vector per particle
          "yes",    // restart
          "yes",     // communicate ghost
          "yes",    // communicate rev
          "0.","0.","0." // default values
        };
        fix_dragforce_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
      }

    // register hydrodynamic torque
    if(!fix_hdtorque_)
      {
        const char *fixarg[] = {
          "hdtorque", // fix id
          "all",      // fix group
          "property/atom", // fix style: property/atom
          "hdtorque",      // property name
          "vector", // 1 vector per particle
          "yes",    // restart
          "yes",     // communicate ghost
          "yes",    // communicate rev
          "0.","0.","0."
        };
        fix_hdtorque_ = modify->add_fix_property_atom(11,const_cast<char**>(fixarg),style);
      }
  }

  void FixLbCouplingOnetoone::pre_delete(bool unfixflag)
  {
    if(!unfixflag) return;

    if(fix_dragforce_)
      modify->delete_fix(fix_dragforce_->id);
    if(fix_hdtorque_)
      modify->delete_fix(fix_hdtorque_->id);
  }

  void FixLbCouplingOnetoone::post_force(int)
  {
    double **f_ext = fix_dragforce_->array_atom;
    //double **t_ext = fix_hdtorque_->array_atom;
    double **t_ext = fix_hdtorque_->array_atom;
    double **f = atom->f;
    //double **t = atom->torque;

    // for(int i=0;i<atom->nlocal;i++)
    //   std::cout << comm->me << " force_liggghts "
    //             << std::setprecision(12) << f_ext[i][2] << std::endl;

    // std::cout << "before" << std::endl;
    // fix_dragforce_->do_reverse_comm();
    // fix_hdtorque_->do_reverse_comm();
    // std::cout << "after" << std::endl;

    // for(int i=0;i<atom->nlocal;i++)
    //   std::cout << comm->me << " force_liggghts_after "
    //             << std::setprecision(12) << f_ext[i][2] << " x " << atom->x[0][2] << std::endl;
    double **t = atom->torque;

    for(int i=0;i<atom->nlocal;i++){
      f[i][0] += f_ext[i][0];
      f[i][1] += f_ext[i][1];
      f[i][2] += f_ext[i][2];
      if(use_torque_){
	t[i][0] += t_ext[i][0];
	t[i][1] += t_ext[i][1];
	t[i][2] += t_ext[i][2];
      }
    }

  }

  // LBDEM coupling relies on correct values for velocities in ghost atoms, thus 

  void FixLbCouplingOnetoone::post_run()
  {
    // need one very last forward_comm to make sure 
    // that velocities on owned and ghost particles match
    timer->stamp();
    comm->forward_comm();
    timer->stamp(TIME_COMM);

  }

  double **FixLbCouplingOnetoone::get_force_ptr()
  {
    return fix_dragforce_->array_atom;
  }
  double **FixLbCouplingOnetoone::get_torque_ptr()
  {
    return fix_hdtorque_->array_atom;
  }
  void FixLbCouplingOnetoone::comm_force_torque()
  {
    fix_dragforce_->do_reverse_comm();
    fix_hdtorque_->do_reverse_comm();
  }


}; /* LAMMPS_NS */
