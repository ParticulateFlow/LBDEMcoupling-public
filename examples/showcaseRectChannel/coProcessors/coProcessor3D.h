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

/** \file
 * Helper functions for domain initialization -- header file.
 */
#ifndef CO_PROCESSOR_3D_H
#define CO_PROCESSOR_3D_H

#include <map>

#include "atomicBlock/blockLattice3D.h"
#include "core/geometry3D.h"
#include "core/globalDefs.h"
#include "latticeBoltzmann/nearestNeighborLattices3D.h"

namespace plb {

/// A co-processor provides access to a computational hardware unit, such as a
///   GPU or a FPGA.
/** An instance of the CoProcessor3D class is considered to represent a single
 *  hardware unit. In a multi-GPU machine for instance, a new CoProcessor3D is
 *  instantiated for each GPU.
 *
 *  A co-processor acts exclusively on rectangular domains, and can be
 *  responsible for more than one domain. The method addDomain is used to add
 *  new domains for which the co-processor is reponsible.
 *
 *  The memory is considered to be duplicated. It allocated once on the CPU by
 *  Palabos and once on the device by the co-processor. The send() and
 *  receive() methods are responsible for communication between the two memory
 *  spaces, while the collideAndStream() method works on device memory only.
 *
 *  At this stage, co-processors implement only BGK dynamics on a D3Q19
 *  lattice. Also, only the collide-and-stream operation is performed by the
 *  device at this point. Both these aspects will be generalized in the future.
 **/
template <typename T>
class CoProcessor3D {
public:
    virtual ~CoProcessor3D() { }

    /// Add a domain for which the co-processor will perform computations.
    /** All domains range from 0 to nx-1, from 0 to ny-1, and from 0 to nz-1 at
     *  the present interface representation, no matter where they are actually
     *  placed in the physical space.
     *
     *  The relaxation parameter omega is used to implement the BGK collision rule
     *  on the device.
     *
     *  A handle "domainHandle" is returned by the co-processor, and is subsequently
     *  used to identify the various domains during the calls to send(), receive(),
     *  and collideAndStream().
     *
     *  The method returns an error code: 1=success, 0=failure.
     **/
    virtual int addDomain(plint nx, plint ny, plint nz, T omega, int &domainHandle) = 0;

    /// Copy data from Palabos' CPU memory to the co-processors' device memory.
    /** The method returns an error code: 1=success, 0=failure.
     *  Please note that the memory of a std::vector is always contiguous, which means that
     *  you can get a c-array representation of the data through the syntax
     *  T const* carray = &data[0].
     *
     *  The memory layout must respect the following ordering:
     *  - The fastest running index is for the 19 populations, with an ordering specified
     *    in the structure "D3Q19Constants" in the file
     *    "latticeBoltzmann/nearestNeighborLattices3D.hh".
     *  - The space indices are ordered according to the C convention, meaning that, if you
     *    take the space matrix to be declared as matrix[nx][ny][nz], then the z-index is
     *    fastest running.
     **/
    virtual int send(int domainHandle, Box3D const &subDomain, std::vector<char> const &data) = 0;

    /// Copy data from the co-processors' device memory to Palabos' CPU memory.
    /** The method returns an error code: 1=success, 0=failure.
     *  Further information on the memory layout is available in the documentation of the
     *  method send().
     *
     *  Attention: it is the responsibility of the receive method to resize the data vector
     *  so it is big enough.
     **/
    virtual int receive(
        int domainHandle, Box3D const &subDomain, std::vector<char> &data) const = 0;

    /// Execute a collision step on each cell, and then a streaming on the full domain.
    /** Note that the result of the streaming step is undefined in a one-cell layer at the
     *  outer border of the domain. The method collideAndStream() is free to produce
     *  whatever result it wishes inside this layer.
     *
     *  It is also mentioned that the collideAndStream() operation is blocking: it does
     *  not terminated before the operation is fully completed. In order to overlay
     *  computations, you must use the MPI-based multi-thread mechanism in Palabos.
     **/
    virtual int collideAndStream(int domainHandle) = 0;
};

/// This place-holder co-processor does nothing else than implement all functionalities
///   of a co-processor on the CPU, using Palabos library calls. It is in principle
///   useless, but can be used to debug or illustrate the mechanism for calling
///   co-processors in Palabos.
template <typename T>
class D3Q19ExampleCoProcessor3D : public CoProcessor3D<T> {
public:
    virtual int addDomain(plint nx, plint ny, plint nz, T omega, int &domainHandle);
    virtual int send(int domainHandle, Box3D const &subDomain, std::vector<char> const &data);
    virtual int receive(int domainHandle, Box3D const &subDomain, std::vector<char> &data) const;
    virtual int collideAndStream(int domainHandle);

private:
    std::map<int, BlockLattice3D<T, descriptors::D3Q19Descriptor> > domains;
};

template <typename T>
class D3Q19CudaCoProcessor3D : public CoProcessor3D<T> {
public:
    D3Q19CudaCoProcessor3D();
    virtual ~D3Q19CudaCoProcessor3D();
    /// I'd suggest to allocate the domain in device memory at this point.
    virtual int addDomain(plint nx, plint ny, plint nz, T omega, int &domainHandle);
    //    virtual int addDomain(Box3D const& domain, T omega);
    /// Copy from CPU to device memory.
    virtual int send(int domainHandle, Box3D const &subDomain, std::vector<char> const &data);
    /// Copy from device memory to CPU.
    virtual int receive(int domainHandle, Box3D const &subDomain, std::vector<char> &data) const;
    /// Execute BGK collision and streaming in device memory.
    virtual int collideAndStream(int domainHandle);

private:
    /*
    PyObject *createCUDASimulation(void);
    void prepareSimulationToRun(void);

    PyObject *geometryClass;
    PyObject *simulationClass;
    PyObject *simulation;
    PyObject *pyCUDA;
    PyObject *copyFunctionDtoH;
    PyObject *copyFunctionHtoD;
    PyObject *buffer;
    float *data;
    int nx;
    int ny;
    int nz;
    int basis;
    int size;
    int floatSize;
    int strideXSailfish;
    int strideYSailfish;
    int strideZSailfish;
    int strideYPalabos;
    int strideZPalabos;
    */
};

namespace global {

template <typename T>
inline CoProcessor3D<T> &defaultCoProcessor3D()
{
    /// @Tomasz: At this stage, I'd suggest that you simply replace the
    //  D3Q19ExampleCoProcessor3D singleton by a
    //  D3Q19CudaCoProcessor3D singleton to get your code running.
    static D3Q19ExampleCoProcessor3D<T> singleton;
    // static D3Q19CudaCoProcessor3D<T> singleton;
    return singleton;
}

}  // namespace global

}  // namespace plb

#endif  // CO_PROCESSOR_3D_H
