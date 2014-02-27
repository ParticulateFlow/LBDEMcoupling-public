namespace plb {

  /*
   * implementation of PeriodicPressureFunctional3D
   */

  template<typename T, template<typename U> class Descriptor>
  PeriodicPressureFunctional3D<T,Descriptor>::PeriodicPressureFunctional3D(T const deltaRho, T const rhoAvg, 
                                                                           plint const dimension, plint const direction)
    : /*rescaleFactor((1.+3.*deltaRho)/rhoAvg)*/rescaleFactor((1.+deltaRho)/rhoAvg)
  {
    pcout << rescaleFactor - 1. << std::endl;
    // to do: perform error checking here
    for(plint iPop=0;iPop<Descriptor<T>::q;iPop++){
      if(Descriptor<T>::c[iPop][dimension] == direction){
        rescalePop.push_back(iPop);      
      }
    }
  }
  
template<typename T, template<typename U> class Descriptor>
  void PeriodicPressureFunctional3D<T,Descriptor>::processCell(Cell<T,Descriptor>& cell)
  {
    for(IndexVec::iterator it = rescalePop.begin();
        it != rescalePop.end(); it++){
      plint ii = *it;
      T fTmp = cell[ii] + Descriptor<T>::t[ii];
      fTmp *= rescaleFactor;
      cell[ii] = fTmp - Descriptor<T>::t[ii];
    }
  }

  template<typename T, template<typename U> class Descriptor>
  void PeriodicPressureFunctional3D<T,Descriptor>::process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice)
  {
    for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
      for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
        for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
          processCell(lattice.get(iX,iY,iZ));
        }
      }
    }
  }

  template<typename T, template<typename U> class Descriptor>
  void PeriodicPressureFunctional3D<T,Descriptor>::getTypeOfModification(std::vector<modif::ModifT>& modified) const
  {
    //modified[0] = modif::nothing;
    modified[1] = modif::staticVariables;
  }

  template<typename T, template<typename U> class Descriptor>
  PeriodicPressureFunctional3D<T,Descriptor>* PeriodicPressureFunctional3D<T,Descriptor>::clone() const 
  { 
    return new PeriodicPressureFunctional3D<T,Descriptor>(*this);
  }

  /* ---------------------------------------------------- */

  // template<typename T, template<typename U> class Descriptor>
  // StripeOffEquilibriumFunctional3D<T,Descriptor>::StripeOffEquilibriumFunctional3D(plint const dimension, plint const direction)
  //   : xoff(dimension==0?-direction:0),
  //     yoff(dimension==1?-direction:0),
  //     zoff(dimension==2?-direction:0)
  // {
  //   for(plint iPop=0;iPop<Descriptor<T>::q;iPop++){
  //     if(Descriptor<T>::c[iPop][dimension] == direction){
  //       rescalePop.push_back(iPop);      
  //     }
  //   }
    
  // }

  // template<typename T, template<typename U> class Descriptor>
  // void StripeOffEquilibriumFunctional3D<T,Descriptor>::process(Box3D domain, BlockLattice3D<T,Descriptor>& lattice)
  // {
  //   for (plint iX=domain.x0; iX<=domain.x1; ++iX) {
  //     for (plint iY=domain.y0; iY<=domain.y1; ++iY) {
  //       for (plint iZ=domain.z0; iZ<=domain.z1; ++iZ) {
  //         Cell<T,Descriptor> c = lattice.get(iX,iY,iZ);
  //         Cell<T,Descriptor> c2 = lattice.get(iX+xoff,iY+yoff,iZ+zoff);
  //         T rhoBar;
  //         Array<T,Descriptor<T>::d> j;
  //         momentTemplates<T,Descriptor>::get_rhoBar_j(c2, rhoBar, j);
  //         jSqr = j[0]*j[0]+j[1]*j[1]+j[2]*j[2];
  //         for(plint i=0;i<rescalePop.size();i++){
  //           c[rescalePop[i]] -= c.computeEquilibrium(rescalePop[i],rhoBar,j,jSqr);
  //         }
  //       }
  //     }
  //   }
  // }

};
