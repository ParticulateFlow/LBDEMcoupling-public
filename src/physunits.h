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

#ifndef PHYSUNITS_H
#define PHYSUNITS_H

#include "palabos3D.h"
#include "palabos3D.hh"
#include <cmath>

namespace plb {
  template<typename T>
  class PhysUnits3D {
  public:
    PhysUnits3D(T const l_p_, T const u_p_, T const nu_p_, T const lx_, T const ly_, T const lz_,
                plint const N, T const uMax, T const rho_p_)
      : l_p(l_p_), u_p(u_p_), t_p(l_p_/u_p_), nu_p(nu_p_), lx(lx_), ly(ly_), lz(lz_), rho_p(rho_p_),
        Re(l_p_*u_p_/nu_p_), param(uMax,Re,N,lx_/l_p_,ly_/l_p_,lz/l_p_),
        forceFactor(pow(param.getDeltaX(),4)/pow(param.getDeltaT(),2) * rho_p * pow(l_p,4)/pow(t_p,2)),
        torqueFactor(pow(param.getDeltaX(),5)/pow(param.getDeltaT(),2) * rho_p * pow(l_p,5)/pow(t_p,2)),
        lbOffsetX(0.5), lbOffsetY(0.5), lbOffsetZ(0.5)
    {}

    T getRe() const {return Re; };

    void setLbOffset(T const lbOffsetX_,T const lbOffsetY_,T const lbOffsetZ_)
    {
      lbOffsetX = lbOffsetX_;
      lbOffsetY = lbOffsetY_;
      lbOffsetZ = lbOffsetZ_;
    }

    T getLbLength(T const physL) const { return physL / l_p / param.getDeltaX() ; };
    T getLbPosition(T const physL) const { return physL / l_p / param.getDeltaX() - 0.5; };
    T getLbPositionX(T const physL) const { return physL / l_p / param.getDeltaX() - lbOffsetX ; };
    T getLbPositionY(T const physL) const { return physL / l_p / param.getDeltaX() - lbOffsetY ; };
    T getLbPositionZ(T const physL) const { return physL / l_p / param.getDeltaX() - lbOffsetZ ; };
    T getLbVel(T const physVel) const { return physVel / u_p * param.getLatticeU(); };
    T getLbTime(T const physTime) const { return physTime / t_p / param.getDeltaT(); };
    T getLbDensity(T const physDensity) const { return physDensity / rho_p; };
    T getLbAccel(T const physAccel) const { return physAccel * (t_p*t_p/l_p)
        * (param.getDeltaT()*param.getDeltaT()/param.getDeltaX()); };
    T getLbFreq(T const physFreq) const { return physFreq * t_p * param.getDeltaT(); };
    plint getLbSteps(T const physTime) const { return (plint) (physTime/param.getDeltaT()/t_p + 0.5); };
    T getLbRho(T const physPress, T const pOffset = 0) const { return (physPress - pOffset) * 3. 
        / getPhysForce(1.) * pow(getPhysLength(1.),2); };
    
    T getPhysLength(T const lbLength) const { return lbLength*param.getDeltaX()*l_p; };
    T getPhysPosition(T const lbLength) const { return (lbLength+0.5)*param.getDeltaX()*l_p; };
    T getPhysPositionX(T const lbLength) const { return (lbLength+lbOffsetX)*param.getDeltaX()*l_p; };
    T getPhysPositionY(T const lbLength) const { return (lbLength+lbOffsetY)*param.getDeltaX()*l_p; };
    T getPhysPositionZ(T const lbLength) const { return (lbLength+lbOffsetZ)*param.getDeltaX()*l_p; };
    T getPhysVel(T const lbVel) const { return lbVel/param.getLatticeU()*u_p; };
    T getPhysTime(T const lbTime) const { return lbTime*param.getDeltaT()*t_p; };
    T getPhysDensity(T const lbDensity) const { return lbDensity * rho_p; };
    T getPhysAccel(T const lbAccel) const { return lbAccel / getLbAccel(1.); };    
    T getPhysForce(T const lbForce) const { return lbForce * forceFactor; };
    T getPhysTorque(T const lbTorque) const { return lbTorque * torqueFactor; };
    T getPhysPress(T const lbRho, T pOffset = 0) const { return (lbRho - 1.) / 3. * getPhysForce(1) 
        / pow(getPhysLength(1.),2) + pOffset; };


    const IncomprFlowParam<T>& getLbParam() const
    {
      return param;
    }

  private:
    T const l_p, u_p, t_p, nu_p;
    T const lx,ly,lz, rho_p;
    T const Re;
    IncomprFlowParam<T> param;
    T const forceFactor, torqueFactor;
    T lbOffsetX, lbOffsetY, lbOffsetZ;
  };

};

#endif // PHYSUNITS_H
