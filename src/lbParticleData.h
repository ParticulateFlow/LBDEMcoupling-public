
template<typename T>
class LbParticleData{
public:
  typedef std::vector<Array<T,3> > ParticleDataArrayVector;
  typedef std::vector<T> ParticleDataScalarVector;
  //  LbParticleData()
  updateParticledata(LiggghtsCouplingWrapper const &wrapper, PhysUnits3D const &units)
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
      id.push_back((T)wrapper.id[iPart][0]);
      // add zeros for force and torque
      f.push_back(Array<T,3>(0.,0.,0.));
      t.push_back(Array<T,3>(0.,0.,0.));
    }
  }

  ParticleDataArrayVector x,v,f,t,omega;
  ParticleDataScalarVector r,id;

};
