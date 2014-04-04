namespace plb {

  /*
   * implementation of PeriodicPressureFunctional3D
   */

  template<typename T, template<typename U> class Descriptor>
  PeriodicPressureFunctional3D<T,Descriptor>::PeriodicPressureFunctional3D(T const rhoTarget_, 
                                                                           T const rhoAvg_,
                                                                           plint const dimension_, 
                                                                           plint const direction_)
    : rhoTarget(rhoTarget_), rhoAvg(rhoAvg_), dimension(dimension_), direction(direction_)
  {

    // to do: perform error checking here
    for(plint iPop=0;iPop<Descriptor<T>::q;iPop++){
      if(Descriptor<T>::c[iPop][dimension] == direction){
        rescalePop.push_back(iPop);      
      }
    }
  }
  
  template<typename T, template<typename U> class Descriptor>
  T PeriodicPressureFunctional3D<T,Descriptor>::computeEquilibrium(Cell<T,Descriptor>& cell,
                                                                   plint iEq,
                                                                   T rho,
                                                                   Array<T,Descriptor<T>::d> u)
  {
    T jSqr = 0;
    for(plint i=0;i<Descriptor<T>::d;i++){
      u[i] *= rho;
      jSqr += u[i]*u[i]; // u now holds j (momentum) but no additional variable needs to be declared
    }
    return cell.computeEquilibrium(iEq,rho,u,jSqr);
  }

  template<typename T, template<typename U> class Descriptor>
  void PeriodicPressureFunctional3D<T,Descriptor>::process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice)
  {
    plint dx = dimension == 0 ? -direction : 0;
    plint dy = dimension == 1 ? -direction : 0;
    plint dz = dimension == 2 ? -direction : 0;

    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
      for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
        for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
          Cell<T,Descriptor> c = lattice.get(iX,iY,iZ);
          Cell<T,Descriptor> cNeigh = lattice.get(iX+dx,iY+dy,iZ+dz);

          T rhoNeigh = cNeigh.computeDensity();
          Array<T,Descriptor<T>::d> u;
          cNeigh.computeVelocity(u);

          T rhoTargetCell = rhoNeigh + (rhoTarget - rhoAvg);
          pcout << "rho " 
                << rhoNeigh << " " << rhoAvg << " " 
                << rhoTarget << " " << rhoTargetCell << std::endl;
          pcout << "pop ";
          for(IndexVec::iterator it = rescalePop.begin();
              it != rescalePop.end(); it++){
            plint ii = *it;
            pcout << ii << " " << c[ii] << " ";
            T fEq = computeEquilibrium(cNeigh,ii,rhoNeigh,u);
            c[ii] += (computeEquilibrium(c,ii,rhoTargetCell,u) - computeEquilibrium(cNeigh,ii,rhoNeigh,u));
            T fNeq = c[ii] - computeEquilibrium(cNeigh,ii,rhoNeigh,u);
            c[ii] = fNeq + computeEquilibrium(c,ii,rhoTargetCell,u);
            pcout << c[ii] << " | ";
          }
          pcout << std::endl;
        }
      }
    }
  }

  template<typename T, template<typename U> class Descriptor>
  void PeriodicPressureFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
  {
    modified[0] = modif::nothing; //modif::staticVariables;
  }

  template<typename T, template<typename U> class Descriptor>
  PeriodicPressureFunctional3D<T,Descriptor>* PeriodicPressureFunctional3D<T,Descriptor>::clone() const 
  { 
    return new PeriodicPressureFunctional3D<T,Descriptor>(*this);
  }

  
};
