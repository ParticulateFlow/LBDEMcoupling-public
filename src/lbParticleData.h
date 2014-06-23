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
template<typename T>
class LbParticleData{
public:
  typedef std::vector<Array<T,3> > ParticleDataArrayVector;
  typedef std::vector<T> ParticleDataScalarVector;
  //  LbParticleData()
  updateParticleData(LiggghtsCouplingWrapper const &wrapper, PhysUnits3D const &units)
  {
    plint nPart = x.size();
    plint nPartNew = wrapper.getNumParticles();
    for(plint iPart=0;iPart<nPart;i++){
      r[iPart] = units.getLbLength(wrapper.r[iPart][0]);
      x[iPart].from_cArray(wrapper.x[iPart]);
      v[iPart].from_cArray(wrapper.v[iPart]);
      omega[iPart].from_cArray(wrapper.omega[iPart]);
      id[iPart] = (T)wrapper.id[iPart][0];
    }
    for(plint iPart=nPart;iPart<nPartNew;i++){
      r.push_back(wrapper.r[iPart][0]);
      Array<T,3> tmp;
      tmp.from_cArray(wrapper.x[iPart]);
      x.push_back(tmp);
      tmp.from_cArray(wrapper.v[iPart]);
      v.push_back(tmp);
      tmp.from_cArray(wrapper.omega[iPart]);
      omega.push_back(tmp);
      id.push_back( round( (T)wrapper.id[iPart][0] ) );
      // add zeros for force and torque
      f.push_back(Array<T,3>(0.,0.,0.));
      t.push_back(Array<T,3>(0.,0.,0.));
    }
  }

  ParticleDataArrayVector x,v,f,t,omega;
  ParticleDataScalarVector r,id;

};
