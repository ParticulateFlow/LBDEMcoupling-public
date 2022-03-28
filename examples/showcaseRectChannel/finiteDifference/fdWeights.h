/* This file is part of the Palabos library.
 *
 * The Palabos softare is developed since 2011 by FlowKit-Numeca Group Sarl
 * (Switzerland) and the University of Geneva (Switzerland), which jointly
 * own the IP rights for most of the code base. Since October 2019, the
 * Palabos project is maintained by the University of Geneva and accepts
 * source code contributions from the community.
 *
 * Contact:
 * Jonas Latt
 * Computer Science Department
 * University of Geneva
 * 7 Route de Drize
 * 1227 Carouge, Switzerland
 * jonas.latt@unige.ch
 *
 * The most recent release of Palabos can be downloaded at
 * <https://palabos.unige.ch/>
 *
 * The library Palabos is free software: you can redistribute it and/or
 * modify it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * The library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef FD_WEIGHTS_H
#define FD_WEIGHTS_H

namespace plb {

// An FdWeights object is instantiated to be used in the following way. The weights of
// finite difference stencils will be calculated according to the user's requirements (see the
// "getWeights" functions). The order of the derivative to be approximated in all cases
// is the template parameter "order" and it must be greater than or equal to 0. If it is 0, then
// function interpolation coeficients on a uniform grid are computed. For a given object
// the derivative order is fixed, and the user can compute weights of stencils of
// different widths. The weights are hashed, so if the user asks twice for the same
// weights, they are not recomputed for efficiency reasons. To achieve hashing, the
// template parameter "maxWidth" is needed. The user is commited to never request the computation
// of weights for a stencil with width greater than this parameter (maxWidth >= 1).
// In all cases the derivative of the specified order is approximated as the sum of the
// products of the computed weights times the function values on the points of the stencil.
// The FdWeights class is singleton. It must be used only by calling the "fdWeights" function.
template <typename T, int order, int maxWidth>
class FdWeights {
public:
    // This function computes the weights of a finite difference stencil to calculate an
    // approximation of the derivative of order "order" for a grid-based function. It returns a
    // pointer to the beginning of "width" number of weights. The user of this function, must not
    // deallocate the memory returned by this pointer after its use. The stencil extends over
    // "width" equispaced grid points (with dx = 1, and width <= maxWidth). The "position" is the
    // index of the specific grid location in the stencil to compute the derivative on. So, if we
    // have order = 2, width = 3, and position = 1, then we get the weights to compute the second
    // derivative central finite difference approximation with 3 points. If position = 0, then we
    // get the forward difference, and if position = 2, we get the backward difference
    // approximation. This function is optimized, and when it is called with the same (width,
    // position) combination, it does not recompute the weights. So, it is important that for
    // approximations of the same derivative order with different stencils, only one instantiation
    // of FdWeights is used.
    T const *getWeights(int width, int position) const;

    // The weights can be computed also for any real position inside the interval spanned by the
    // stencil. This position does not need to be an integer corresponding to a grid location (like
    // in the function above). When the position is a floating-point variable "z", the weights are
    // computed, but they are not hashed as done in the integer case. The user needs to provide (and
    // be responsible) for the memory to store these weights. A pointer to this memory (of size at
    // least equal to "width"), is provided to this function as the parameter "w". This version of
    // the "getWeights" function, might be useful to someone who wants to perform interpolation
    // (order = 0) at an off-lattice position, using a wide stencil. Naturally, it can also be used
    // to evaluate derivatives of any order at an off-lattice position. Since, there is no hashing
    // involved when using this function, "width" does not need to be less equal to "maxWidth". "z"
    // must belong to the interval [0, width-1]. If one knows that he will be using only this
    // function (and never the one above), then no hashing will ever be performed. To save memory,
    // the user can instantiate an object with "maxStencilWidth" equal to 1, and then use the
    // following "getWeights" function with arbitrary "width".
    void getWeights(int width, T z, T *w) const;

private:
    // The construction and distruction is private to make this class singleton.
    FdWeights();
    FdWeights(FdWeights<T, order, maxWidth> const &rhs);
    void swap(FdWeights<T, order, maxWidth> &rhs);
    FdWeights<T, order, maxWidth> &operator=(FdWeights<T, order, maxWidth> const &rhs);
    FdWeights<T, order, maxWidth> *clone() const;
    ~FdWeights();

    void computeWeights(int width, T z, T *w) const;

private:
    mutable T *
        *weights;  // Finite difference weights for a many stencil widths, and stencil positions.

    template <typename T_, int order_, int maxWidth_>
    friend FdWeights<T_, order_, maxWidth_> const &fdWeights();
};

template <typename T, int order, int maxWidth>
inline FdWeights<T, order, maxWidth> const &fdWeights()
{
    static FdWeights<T, order, maxWidth> fdWeightsInstance;
    return fdWeightsInstance;
}

}  // namespace plb

#endif  // FD_WEIGHTS_H
