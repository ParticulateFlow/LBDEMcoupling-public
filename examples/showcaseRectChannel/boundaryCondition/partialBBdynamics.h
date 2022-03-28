#ifndef PARTIAL_BB_DYNAMICS_H
#define PARTIAL_BB_DYNAMICS_H

#include "basicDynamics/isoThermalDynamics.h"
#include "core/dynamics.h"
#include "core/globalDefs.h"

namespace plb {

/// Implementation of partially saturated bounce back dynamics (Partially Saturated Method)
template <typename T, template <typename U> class Descriptor>
class PartialBBdynamics : public IsoThermalBulkDynamics<T, Descriptor> {
public:
    /* *************** Construction / Destruction ************************ */
    PartialBBdynamics(T omega_);
    PartialBBdynamics(HierarchicUnserializer &unserializer);

    /// Clone the object on its dynamic type.
    virtual PartialBBdynamics<T, Descriptor> *clone() const;

    /// Return a unique ID for this class.
    virtual int getId() const;

    /* *************** Collision and Equilibrium ************************* */

    /// Implementation of the collision step
    virtual void collide(Cell<T, Descriptor> &cell, BlockStatistics &statistics_);
    /// Implementation of the collision step, with imposed macroscopic variables
    virtual void collideExternal(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j, T thetaBar,
        BlockStatistics &stat);
    /// Compute equilibrium distribution function
    virtual T computeEquilibrium(
        plint iPop, T rhoBar, Array<T, Descriptor<T>::d> const &j, T jSqr, T thetaBar = T()) const;

private:
    // handle the partially saturated method collision
    static T PSMCollision(
        Cell<T, Descriptor> &cell, T rhoBar, Array<T, Descriptor<T>::d> const &j,
        Array<T, Descriptor<T>::d> const &wallVelocity, T &solidFraction, T omega);
    virtual void decomposeOrder0(Cell<T, Descriptor> const &cell, std::vector<T> &rawData) const;
    virtual void recomposeOrder0(Cell<T, Descriptor> &cell, std::vector<T> const &rawData) const;

private:
    static int id;
};

}  // namespace plb

#endif  // PARTIAL_BB_DYNAMICS_H