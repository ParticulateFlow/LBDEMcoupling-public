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

#ifndef CO_PROCESSOR_3D_HH
#define CO_PROCESSOR_3D_HH

#include <limits>

#include "basicDynamics/isoThermalDynamics.h"
#include "coProcessors/coProcessor3D.h"

namespace plb {

template <typename T>
int D3Q19ExampleCoProcessor3D<T>::addDomain(
    plint nx, plint ny, plint nz, T omega, int &domainHandle)
{
    PLB_ASSERT((int)domains.size() < std::numeric_limits<int>::max());
    domainHandle = (int)domains.size();
    Dynamics<T, descriptors::D3Q19Descriptor> *dynamics =
        new BGKdynamics<T, descriptors::D3Q19Descriptor>(omega);
    domains.insert(std::pair<int, BlockLattice3D<T, descriptors::D3Q19Descriptor> >(
        domainHandle, BlockLattice3D<T, descriptors::D3Q19Descriptor>(nx, ny, nz, dynamics)));
    return 1;  // Success.
}

template <typename T>
int D3Q19ExampleCoProcessor3D<T>::send(
    int domainHandle, Box3D const &subDomain, std::vector<char> const &data)
{
    typename std::map<int, BlockLattice3D<T, descriptors::D3Q19Descriptor> >::iterator it =
        domains.find(domainHandle);
    PLB_ASSERT(it != domains.end());
    BlockLattice3D<T, descriptors::D3Q19Descriptor> &lattice = it->second;
    lattice.getDataTransfer().receive(subDomain, data, modif::staticVariables);
    return 1;  // Success.
}

template <typename T>
int D3Q19ExampleCoProcessor3D<T>::receive(
    int domainHandle, Box3D const &subDomain, std::vector<char> &data) const
{
    typename std::map<int, BlockLattice3D<T, descriptors::D3Q19Descriptor> >::const_iterator it =
        domains.find(domainHandle);
    PLB_ASSERT(it != domains.end());
    BlockLattice3D<T, descriptors::D3Q19Descriptor> const &lattice = it->second;
    lattice.getDataTransfer().send(subDomain, data, modif::staticVariables);
    return 1;  // Success.
}

template <typename T>
int D3Q19ExampleCoProcessor3D<T>::collideAndStream(int domainHandle)
{
    typename std::map<int, BlockLattice3D<T, descriptors::D3Q19Descriptor> >::iterator it =
        domains.find(domainHandle);
    PLB_ASSERT(it != domains.end());
    BlockLattice3D<T, descriptors::D3Q19Descriptor> &lattice = it->second;
    lattice.collideAndStream(lattice.getBoundingBox());
    return 1;  // Success.
}

template <typename T>
D3Q19CudaCoProcessor3D<T>::D3Q19CudaCoProcessor3D()
{
    /*
    PyObject *pName, *pModule, *pFunc, *driver, *value;

    char **sailfishArgs;
    sailfishArgs = new char*[1];
    sailfishArgs[0] = "cavity3D.py";
    Py_Initialize();
    PySys_SetArgvEx(1, sailfishArgs, 1);
    delete [] sailfishArgs;

    pName = PyString_FromString("cavity3D");
    pModule = PyImport_Import(pName);
    Py_CLEAR(pName);
    this->geometryClass = PyObject_GetAttrString(pModule, "Geometry");
    this->simulationClass = PyObject_GetAttrString(pModule, "Simulation");
    Py_CLEAR(pModule);

    pName = PyString_FromString("pycuda");
    this->pyCUDA = PyImport_Import(pName);
    driver = PyObject_GetAttrString(pyCUDA, "driver");
    this->copyFunctionHtoD = PyObject_GetAttrString(driver, "memcpy_htod");
    this->copyFunctionDtoH = PyObject_GetAttrString(driver, "memcpy_dtoh");
    Py_CLEAR(pName);
    */
}

template <typename T>
D3Q19CudaCoProcessor3D<T>::~D3Q19CudaCoProcessor3D()
{
    /*
    Py_CLEAR(this->buffer);
    delete [] this->data;

    Py_CLEAR(this->copyFunctionDtoH);
    Py_CLEAR(this->copyFunctionHtoD);
    Py_CLEAR(this->simulation);
    Py_CLEAR(this->pyCUDA);

    Py_CLEAR(this->simulationClass);
    Py_CLEAR(this->geometryClass);
    Py_Finalize();
    */
}

template <typename T>
int D3Q19CudaCoProcessor3D<T>::addDomain(plint nx, plint ny, plint nz, T omega, int &domainHandle)
// int D3Q19CudaCoProcessor3D<T>::addDomain(Box3D const& domain, T omega)
{
    /*

    this->nx = nx;
    this->ny = ny;
    this->nz = nz;
    this->basis = 19;
    this->simulation = createCUDASimulation();
    prepareSimulationToRun();

    PyObject *name, *grid, *s, *strides;
    name = PyString_FromString("grid");
    grid = PyObject_GetAttr(this->simulation, name);
    Py_CLEAR(name);

    name = PyString_FromString("get_dist_bytes");
    s = PyObject_CallMethodObjArgs(this->simulation, name, grid, NULL);
    Py_CLEAR(name);
    Py_CLEAR(grid);
    size = PyInt_AsLong(s);
    Py_CLEAR(s);

    name = PyString_FromString("float");
    grid = PyObject_GetAttr(this->simulation, name);
    Py_CLEAR(name);
    name = PyString_FromString("_get_strides");
    s = PyObject_CallMethodObjArgs(this->simulation, name, grid, NULL);
    Py_CLEAR(name);
    Py_CLEAR(grid);
    strides = PyTuple_GetItem(s, 0);
    this->floatSize = PyInt_AsLong(PyTuple_GetItem(strides, 2));
    size /= floatSize;
    this->strideZSailfish = PyInt_AsLong(PyTuple_GetItem(s, 1));
    this->strideYSailfish = PyInt_AsLong(PyTuple_GetItem(strides, 0))/floatSize;
    this->strideXSailfish = PyInt_AsLong(PyTuple_GetItem(strides, 1))/floatSize;
    Py_CLEAR(s);

    this->data = new float[size];
    this->buffer = PyBuffer_FromReadWriteMemory(data, size*floatSize);

    this->strideYPalabos = this->basis*this->nz*this->ny;
    this->strideZPalabos = this->basis*this->nz;
    */

    return 1;
}

#define SAILFISHD3Q19(x, y, z, b) \
    (b * strideZSailfish + z * strideYSailfish + y * strideXSailfish + x)
#define PALABOSD3Q19(x, y, z, b) (x * strideYPalabos + y * strideZPalabos + z * basis + b)

template <typename T>
int D3Q19CudaCoProcessor3D<T>::send(
    int domainHandle, Box3D const &subDomain, std::vector<char> const &data)
{
    /*
    PyObject *memories, *memory, *args;
    const T *palabosData = reinterpret_cast<const T *>(&data[0]);

    memories = PyObject_CallMethod(this->simulation, "curr_dists", NULL);
    memory = PyList_GetItem(memories, 0);
    Py_INCREF(memory);
    Py_CLEAR(memories);

    int x0, x1, y0, y1, z0, z1;
    x0 = subDomain.x0;
    x1 = subDomain.x1;
    y0 = subDomain.y0;
    y1 = subDomain.y1;
    z0 = subDomain.z0;
    z1 = subDomain.z1;
// Palabos population: src/latticeBoltzman/nearestNeighborLattices3D.hh:159
// Sailfish population: sailfish/sym.py:299
#if 0
    std::cout << "s " << palabosData[PALABOSD3Q19(1, 1, 1, 0)] << " ";
    std::cout << palabosData[PALABOSD3Q19(1, 1, 1, 10)] << " ";
    std::cout << palabosData[PALABOSD3Q19(10, 10, 10, 0)] << " ";
    std::cout << palabosData[PALABOSD3Q19(10, 10, 10, 17)] << std::endl;
for (int i = 0; i < size; i++) {
    this->data[i] = 0.0f;
}
#endif
    for (int i = z0; i <= z1; i++) {
        for (int j = y0; j <= y1; j++) {
            for (int k = x0; k <= x1; k++) {
                this->data[SAILFISHD3Q19(k, j, i, 0)] =
                    palabosData[PALABOSD3Q19(k, j, i, 0)] +descriptors::D3Q19Descriptor<T>::t[0];
                this->data[SAILFISHD3Q19(k, j, i, 1)] =
                    palabosData[PALABOSD3Q19(k, j, i, 10)]+descriptors::D3Q19Descriptor<T>::t[10];
                this->data[SAILFISHD3Q19(k, j, i, 2)] =
                    palabosData[PALABOSD3Q19(k, j, i, 1)] +descriptors::D3Q19Descriptor<T>::t[1];
                this->data[SAILFISHD3Q19(k, j, i, 3)] =
                    palabosData[PALABOSD3Q19(k, j, i, 11)]+descriptors::D3Q19Descriptor<T>::t[11];
                this->data[SAILFISHD3Q19(k, j, i, 4)] =
                    palabosData[PALABOSD3Q19(k, j, i, 2)] +descriptors::D3Q19Descriptor<T>::t[2];
                this->data[SAILFISHD3Q19(k, j, i, 5)] =
                    palabosData[PALABOSD3Q19(k, j, i, 12)]+descriptors::D3Q19Descriptor<T>::t[12];
                this->data[SAILFISHD3Q19(k, j, i, 6)] =
                    palabosData[PALABOSD3Q19(k, j, i, 3)] +descriptors::D3Q19Descriptor<T>::t[3];
                this->data[SAILFISHD3Q19(k, j, i, 7)] =
                    palabosData[PALABOSD3Q19(k, j, i, 13)]+descriptors::D3Q19Descriptor<T>::t[13];
                this->data[SAILFISHD3Q19(k, j, i, 8)] =
                    palabosData[PALABOSD3Q19(k, j, i, 5)] +descriptors::D3Q19Descriptor<T>::t[5];
                this->data[SAILFISHD3Q19(k, j, i, 9)] =
                    palabosData[PALABOSD3Q19(k, j, i, 14)]+descriptors::D3Q19Descriptor<T>::t[14];
                this->data[SAILFISHD3Q19(k, j, i, 10)] =
                    palabosData[PALABOSD3Q19(k, j, i, 4)] +descriptors::D3Q19Descriptor<T>::t[4];
                this->data[SAILFISHD3Q19(k, j, i, 11)] =
                    palabosData[PALABOSD3Q19(k, j, i, 17)]+descriptors::D3Q19Descriptor<T>::t[17];
                this->data[SAILFISHD3Q19(k, j, i, 12)] =
                    palabosData[PALABOSD3Q19(k, j, i, 9)] +descriptors::D3Q19Descriptor<T>::t[9];
                this->data[SAILFISHD3Q19(k, j, i, 13)] =
                    palabosData[PALABOSD3Q19(k, j, i, 18)]+descriptors::D3Q19Descriptor<T>::t[18];
                this->data[SAILFISHD3Q19(k, j, i, 14)] =
                    palabosData[PALABOSD3Q19(k, j, i, 8)] +descriptors::D3Q19Descriptor<T>::t[8];
                this->data[SAILFISHD3Q19(k, j, i, 15)] =
                    palabosData[PALABOSD3Q19(k, j, i, 15)]+descriptors::D3Q19Descriptor<T>::t[15];
                this->data[SAILFISHD3Q19(k, j, i, 16)] =
                    palabosData[PALABOSD3Q19(k, j, i, 7)] +descriptors::D3Q19Descriptor<T>::t[7];
                this->data[SAILFISHD3Q19(k, j, i, 17)] =
                    palabosData[PALABOSD3Q19(k, j, i, 16)]+descriptors::D3Q19Descriptor<T>::t[16];
                this->data[SAILFISHD3Q19(k, j, i, 18)] =
                    palabosData[PALABOSD3Q19(k, j, i, 6)] +descriptors::D3Q19Descriptor<T>::t[6];
            }
        }
    }
#if 0
    std::cout << "s " << this->data[SAILFISHD3Q19(10, 10, 10, 0)] << std::endl;
#endif
    args = PyTuple_New(2);
    PyTuple_SetItem(args, 0, memory);
    Py_INCREF(this->buffer);
    PyTuple_SetItem(args, 1, this->buffer);
    PyObject_Call(this->copyFunctionHtoD, args, NULL);
    Py_CLEAR(args);
#if 0
    std::cout << "s " << this->data[SAILFISHD3Q19(10, 10, 10, 0)] << std::endl;
#endif
    */
    return 1;
}
template <typename T>
int D3Q19CudaCoProcessor3D<T>::receive(
    int domainHandle, Box3D const &subDomain, std::vector<char> &data) const
{
    /*
    PyObject *memories, *memory, *args;

#if 0
for (int i = 0; i < this->size; i++) {
    this->data[i] = 0.0f;
}
#endif

    memories = PyObject_CallMethod(this->simulation, "curr_dists", NULL);
    memory = PyList_GetItem(memories, 0);
    Py_INCREF(memory);
    Py_CLEAR(memories);

    args = PyTuple_New(2);
    Py_INCREF(this->buffer);
    PyTuple_SetItem(args, 0, this->buffer);
    PyTuple_SetItem(args, 1, memory);
    PyObject_Call(this->copyFunctionDtoH, args, NULL);
    Py_CLEAR(args);

    int x0, x1, y0, y1, z0, z1, size;
    x0 = subDomain.x0;
    x1 = subDomain.x1;
    y0 = subDomain.y0;
    y1 = subDomain.y1;
    z0 = subDomain.z0;
    z1 = subDomain.z1;
    size = (x1-x0+1)*(y1-y0+1)*(z1-z0+1)*this->basis*sizeof(T);
    if (data.size() != size) {
        data.resize(size);
    }
    T *palabosData = reinterpret_cast<T *>(&data[0]);
#if 0
for (int i = 0; i < size/sizeof(T); i++) {
    palabosData[i] = 0.0;
}

    std::cout << "r " << this->data[SAILFISHD3Q19(1, 1, 1, 0)] << std::endl;
#endif
    for (int i = z0; i <= z1; i++) {
        for (int j = y0; j <= y1; j++) {
            for (int k = x0; k <= x1; k++) {
                palabosData[PALABOSD3Q19(k, j, i, 0)] =
                    this->data[SAILFISHD3Q19(k, j, i, 0)] -descriptors::D3Q19Descriptor<T>::t[0];
                palabosData[PALABOSD3Q19(k, j, i, 10)] =
                    this->data[SAILFISHD3Q19(k, j, i, 1)] -descriptors::D3Q19Descriptor<T>::t[10];
                palabosData[PALABOSD3Q19(k, j, i, 1)] =
                    this->data[SAILFISHD3Q19(k, j, i, 2)] -descriptors::D3Q19Descriptor<T>::t[1];
                palabosData[PALABOSD3Q19(k, j, i, 11)] =
                    this->data[SAILFISHD3Q19(k, j, i, 3)] -descriptors::D3Q19Descriptor<T>::t[11];
                palabosData[PALABOSD3Q19(k, j, i, 2)] =
                    this->data[SAILFISHD3Q19(k, j, i, 4)] -descriptors::D3Q19Descriptor<T>::t[2];
                palabosData[PALABOSD3Q19(k, j, i, 12)] =
                    this->data[SAILFISHD3Q19(k, j, i, 5)] -descriptors::D3Q19Descriptor<T>::t[12];
                palabosData[PALABOSD3Q19(k, j, i, 3)] =
                    this->data[SAILFISHD3Q19(k, j, i, 6)] -descriptors::D3Q19Descriptor<T>::t[3];
                palabosData[PALABOSD3Q19(k, j, i, 13)] =
                    this->data[SAILFISHD3Q19(k, j, i, 7)] -descriptors::D3Q19Descriptor<T>::t[13];
                palabosData[PALABOSD3Q19(k, j, i, 5)] =
                    this->data[SAILFISHD3Q19(k, j, i, 8)] -descriptors::D3Q19Descriptor<T>::t[5];
                palabosData[PALABOSD3Q19(k, j, i, 14)] =
                    this->data[SAILFISHD3Q19(k, j, i, 8)] -descriptors::D3Q19Descriptor<T>::t[14];
                palabosData[PALABOSD3Q19(k, j, i, 4)] =
                    this->data[SAILFISHD3Q19(k, j, i, 10)]-descriptors::D3Q19Descriptor<T>::t[4];
                palabosData[PALABOSD3Q19(k, j, i, 17)] =
                    this->data[SAILFISHD3Q19(k, j, i, 11)]-descriptors::D3Q19Descriptor<T>::t[17];
                palabosData[PALABOSD3Q19(k, j, i, 9)] =
                    this->data[SAILFISHD3Q19(k, j, i, 12)]-descriptors::D3Q19Descriptor<T>::t[9];
                palabosData[PALABOSD3Q19(k, j, i, 18)] =
                    this->data[SAILFISHD3Q19(k, j, i, 13)]-descriptors::D3Q19Descriptor<T>::t[18];
                palabosData[PALABOSD3Q19(k, j, i, 8)] =
                    this->data[SAILFISHD3Q19(k, j, i, 14)]-descriptors::D3Q19Descriptor<T>::t[8];
                palabosData[PALABOSD3Q19(k, j, i, 15)] =
                    this->data[SAILFISHD3Q19(k, j, i, 15)]-descriptors::D3Q19Descriptor<T>::t[15];
                palabosData[PALABOSD3Q19(k, j, i, 7)] =
                    this->data[SAILFISHD3Q19(k, j, i, 16)]-descriptors::D3Q19Descriptor<T>::t[7];
                palabosData[PALABOSD3Q19(k, j, i, 16)] =
                    this->data[SAILFISHD3Q19(k, j, i, 17)]-descriptors::D3Q19Descriptor<T>::t[16];
                palabosData[PALABOSD3Q19(k, j, i, 6)] =
                    this->data[SAILFISHD3Q19(k, j, i, 18)]-descriptors::D3Q19Descriptor<T>::t[6];
            }
        }
    }
    */
    return 1;
}

template <typename T>
int D3Q19CudaCoProcessor3D<T>::collideAndStream(int domainHandle)
{
    /*
    PyObject_CallMethod(this->simulation, "sim_step", NULL);
    */
    return 1;
}

/*
template<typename T>
PyObject * D3Q19CudaCoProcessor3D<T>::createCUDASimulation(void)
{
PyObject *result = NULL;
PyObject *pSailfishArgsCUDA, *args, *kw, *value;

pSailfishArgsCUDA = PyDict_New();
value = PyString_FromString("cuda");
PyDict_SetItemString(pSailfishArgsCUDA, "backend", value);
Py_CLEAR(value);
value = PyString_FromString("D3Q19");
PyDict_SetItemString(pSailfishArgsCUDA, "grid", value);
Py_CLEAR(value);
value = PyInt_FromLong(this->nx);
PyDict_SetItemString(pSailfishArgsCUDA, "lat_nx", value);
Py_CLEAR(value);
value = PyInt_FromLong(this->ny);
PyDict_SetItemString(pSailfishArgsCUDA, "lat_ny", value);
Py_CLEAR(value);
value = PyInt_FromLong(this->nz);
PyDict_SetItemString(pSailfishArgsCUDA, "lat_nz", value);
Py_CLEAR(value);
PyDict_SetItemString(pSailfishArgsCUDA, "verbose", Py_True);
PyDict_SetItemString(pSailfishArgsCUDA, "batch", Py_True);
value = PyInt_FromLong(64);
PyDict_SetItemString(pSailfishArgsCUDA, "max_iters", value);
Py_CLEAR(value);

args = PyTuple_New(1);
Py_INCREF(this->geometryClass);
PyTuple_SetItem(args, 0, this->geometryClass);

kw = PyDict_New();
PyDict_SetItemString(kw, "defaults", pSailfishArgsCUDA);
Py_CLEAR(pSailfishArgsCUDA);

result = PyObject_Call(this->simulationClass, args, kw);

Py_CLEAR(kw);
Py_CLEAR(args);

return result;
}
*/
/*
template<typename T>
void D3Q19CudaCoProcessor3D<T>:: prepareSimulationToRun(void)
{
PyObject_CallMethod(this->simulation, "_init_shape", NULL);
PyObject_CallMethod(this->simulation, "_init_vis", NULL);
PyObject_CallMethod(this->simulation, "_init_geo", NULL);
PyObject_CallMethod(this->simulation, "_init_post_geo", NULL);
PyObject_CallMethod(this->simulation, "_init_code", NULL);
PyObject_CallMethod(this->simulation, "_init_compute_fields", NULL);
PyObject_CallMethod(this->simulation, "_init_compute_kernels", NULL);
PyObject_CallMethod(this->simulation, "_init_compute_ic", NULL);
PyObject_CallMethod(this->simulation, "_init_output", NULL);
}
*/

}  // namespace plb

#endif  // CO_PROCESSOR_3D_HH
