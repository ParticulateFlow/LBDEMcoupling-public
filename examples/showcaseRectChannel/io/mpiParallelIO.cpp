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

#include "io/mpiParallelIO.h"

#include <cstdio>

#include "core/util.h"
#include "io/parallelIO.h"
#include "io/plbFiles.h"
#include "parallelism/mpiManager.h"

namespace plb {

namespace parallelIO {

void writeRawData_mpi(
    FileName fName, std::vector<plint> const &myBlockIds, std::vector<plint> const &offset,
    std::vector<std::vector<char> > &data, bool appendMode)
{
#ifdef PLB_MPI_PARALLEL
    char fNameBuf[1024];
    if (fName.get().size() < 1024) {
        strcpy(fNameBuf, fName.get().c_str());
    } else {
        plbIOError(std::string("File name is too long: ") + fName.get());
    }
    MPI_File fh;
    plint globalOffset = 0;
    if (appendMode) {
        int err = MPI_File_open(
            MPI_COMM_SELF, fNameBuf, MPI_MODE_APPEND | MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
        plbIOError(err != MPI_SUCCESS, "Could not open file " + fName.get());
        MPI_Offset offset;
        MPI_File_get_position(fh, &offset);
        globalOffset = offset;
        MPI_File_close(&fh);
    }
    int amode = MPI_MODE_CREATE | MPI_MODE_WRONLY;
    if (appendMode) {
        amode = MPI_MODE_WRONLY;
    }
    int err = MPI_File_open(MPI_COMM_SELF, fNameBuf, amode, MPI_INFO_NULL, &fh);
    plbIOError(err != MPI_SUCCESS, "Could not open file " + fName.get());
    bool ioError = false;
    for (plint iBlock = 0; iBlock < (plint)myBlockIds.size(); ++iBlock) {
        plint blockId = myBlockIds[iBlock];
        MPI_Offset nextOffset = 0;
        if (blockId == 0) {
            PLB_ASSERT(offset[blockId] == (plint)data[iBlock].size());
        } else {
            PLB_ASSERT(offset[blockId] - offset[blockId - 1] == (plint)data[iBlock].size());
            nextOffset = offset[blockId - 1];
        }
        err = MPI_File_seek(fh, nextOffset + globalOffset, MPI_SEEK_SET);
        if (err != MPI_SUCCESS) {
            ioError = true;
            break;
        }
        MPI_Status status;
        plint dataSize = (plint)data[iBlock].size();
        const plint maxDataSize = 1000000000;  // 1 GB.
        plint remainingDataSize = dataSize;
        plint numWritten = 0;
        while (remainingDataSize > 0) {
            plint nextDataSize = std::min(maxDataSize, remainingDataSize);
            err = MPI_File_write(fh, &data[iBlock][numWritten], nextDataSize, MPI_CHAR, &status);
            if (err != MPI_SUCCESS) {
                ioError = true;
                break;
            }
            remainingDataSize -= nextDataSize;
            numWritten += nextDataSize;
        }
        if (ioError)
            break;
    }
    err = MPI_File_close(&fh);
    if (err != MPI_SUCCESS) {
        ioError = true;
    }
    plbIOError(ioError, std::string("File access unsuccessful in file ") + fName.get());
#endif
}

void writeRawData_posix(
    FileName fName, std::vector<plint> const &myBlockIds, std::vector<plint> const &offset,
    std::vector<std::vector<char> > &data, bool appendMode)
{
    bool errorFlag = false;
    plint globalOffset = 0;
    if (appendMode) {
        if (global::mpi().getRank() == 0) {
            FILE *fp = fopen(fName.get().c_str(), "a");
            errorFlag = !fp;
            if (!errorFlag) {
                globalOffset = (plint)ftell(fp);
                fclose(fp);
            }
        }
        plbIOError(errorFlag, std::string("Unsuccessful writing into file.") + fName.get());
        global::mpi().bCast(&globalOffset, 1);
    }
    for (plint iProcess = 0; iProcess < global::mpi().getSize(); ++iProcess) {
        if (global::mpi().getRank() == iProcess) {
            FILE *fp = 0;
            if (iProcess == 0) {
                if (appendMode) {
#if defined PLB_MAC_OS_X || defined PLB_BSD || defined PLB_WINDOWS
                    fp = fopen(fName.get().c_str(), "r+b");
#else
                    fp = fopen64(fName.get().c_str(), "r+b");
#endif
                } else {
#if defined PLB_MAC_OS_X || defined PLB_BSD || defined PLB_WINDOWS
                    fp = fopen(fName.get().c_str(), "wb");
#else
                    fp = fopen64(fName.get().c_str(), "wb");
#endif
                }
            } else {
#if defined PLB_MAC_OS_X || defined PLB_BSD || defined PLB_WINDOWS
                fp = fopen(fName.get().c_str(), "r+b");
#else
                fp = fopen64(fName.get().c_str(), "r+b");
#endif
            }
            errorFlag = !fp;
            for (plint iBlock = 0; iBlock < (plint)myBlockIds.size() && !errorFlag; ++iBlock) {
                plint blockId = myBlockIds[iBlock];
                plint nextOffset = 0;
                if (blockId == 0) {
                    PLB_ASSERT(offset[blockId] == (plint)data[iBlock].size());
                } else {
                    PLB_ASSERT(offset[blockId] - offset[blockId - 1] == (plint)data[iBlock].size());
                    nextOffset = offset[blockId - 1];
                }
#if defined PLB_MAC_OS_X || defined PLB_BSD || defined PLB_WINDOWS
                int fSeekVal = fseek(fp, (long int)(nextOffset + globalOffset), SEEK_SET);
#else
                int fSeekVal = fseeko64(fp, nextOffset + globalOffset, SEEK_SET);
#endif
                errorFlag = fSeekVal != 0;
                if (!errorFlag) {
                    plint numWritten = (plint)fwrite(&data[iBlock][0], 1, data[iBlock].size(), fp);
                    errorFlag = numWritten != (plint)data[iBlock].size();
                }
            }
            fclose(fp);
        }
        // IMPORTANT: the following plbIOError implies a synchronization between
        //   MPI threads which is required for algorithmic reasons. If you
        //   remove this line, you should replace it by an mpi barrier.
        plbIOError(errorFlag, std::string("Unsuccessful writing into file.") + fName.get());
    }
}

void writeRawData(
    FileName fName, std::vector<plint> const &myBlockIds, std::vector<plint> const &offset,
    std::vector<std::vector<char> > &data, bool appendMode)
{
    PLB_ASSERT(myBlockIds.size() == data.size());
    if (!appendMode) {
        fName.defaultPath(global::directories().getOutputDir());
        fName.defaultExt("dat");
    }
    if (global::IOpolicy().useParallelIO() && global::mpi().getSize() > 1) {
        writeRawData_mpi(fName, myBlockIds, offset, data, appendMode);
    } else {
        // Works in parallel too, but has no parallel efficiency.
        writeRawData_posix(fName, myBlockIds, offset, data, appendMode);
    }
}

void loadRawData_mpi(
    FileName fName, std::vector<plint> const &myBlockIds, std::vector<plint> const &offset,
    std::vector<std::vector<char> > &data)
{
#ifdef PLB_MPI_PARALLEL
    char fNameBuf[1024];
    if (fName.get().size() < 1024) {
        strcpy(fNameBuf, fName.get().c_str());
    } else {
        strcpy(fNameBuf, "");
    }
    MPI_File fh;
    int err = MPI_File_open(MPI_COMM_SELF, fNameBuf, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
    plbIOError(err != MPI_SUCCESS, "Could not open file " + fName.get());
    bool ioError = false;
    for (plint iBlock = 0; iBlock < (plint)myBlockIds.size(); ++iBlock) {
        plint blockId = myBlockIds[iBlock];
        plint nextOffset = 0;
        plint nextSize = 0;
        if (blockId == 0) {
            nextSize = offset[0];
        } else {
            nextSize = offset[blockId] - offset[blockId - 1];
            nextOffset = offset[blockId - 1];
        }
        data[iBlock].resize(nextSize);
        MPI_File_seek(fh, nextOffset, MPI_SEEK_SET);
        if (err != MPI_SUCCESS) {
            ioError = true;
            break;
        }
        MPI_Status status;
        const plint maxDataSize = 1000000000;  // 1 GB.
        plint remainingDataSize = nextSize;
        plint numRead = 0;
        while (remainingDataSize > 0) {
            plint nextPiece = std::min(maxDataSize, remainingDataSize);
            err = MPI_File_read(fh, &data[iBlock][numRead], nextPiece, MPI_CHAR, &status);
            if (err != MPI_SUCCESS) {
                ioError = true;
                break;
            }
            remainingDataSize -= nextPiece;
            numRead += nextPiece;
        }
        if (ioError)
            break;
    }
    err = MPI_File_close(&fh);
    if (err != MPI_SUCCESS) {
        ioError = true;
    }
    plbIOError(ioError, std::string("File access unsuccessful in file ") + fName.get());
#endif
}

void loadRawData_posix(
    FileName fName, std::vector<plint> const &myBlockIds, std::vector<plint> const &offset,
    std::vector<std::vector<char> > &data)
{
    for (plint iProcess = 0; iProcess < global::mpi().getSize(); ++iProcess) {
        bool errorFlag = false;
        if (global::mpi().getRank() == iProcess) {
#if defined PLB_MAC_OS_X || defined PLB_BSD || defined PLB_WINDOWS
            FILE *fp = fopen(fName.get().c_str(), "rb");
#else
            FILE *fp = fopen64(fName.get().c_str(), "rb");
#endif
            errorFlag = !fp;
            for (plint iBlock = 0; iBlock < (plint)myBlockIds.size() && !errorFlag; ++iBlock) {
                plint blockId = myBlockIds[iBlock];
                plint nextOffset = 0;
                plint nextSize = 0;
                if (blockId == 0) {
                    nextSize = offset[0];
                } else {
                    nextSize = offset[blockId] - offset[blockId - 1];
                    nextOffset = offset[blockId - 1];
                }
                data[iBlock].resize(nextSize);
#if defined PLB_MAC_OS_X || defined PLB_BSD || defined PLB_WINDOWS
                int fSeekVal = fseek(fp, (long int)nextOffset, SEEK_SET);
#else
                int fSeekVal = fseeko64(fp, nextOffset, SEEK_SET);
#endif
                errorFlag = fSeekVal != 0;
                if (!errorFlag && nextSize > 0) {
                    plint numRead = (plint)fread(&data[iBlock][0], 1, nextSize, fp);
                    errorFlag = numRead != nextSize;
                }
            }
            if (!errorFlag) {
                fclose(fp);
            }
        }
        // IMPORTANT: the following plbIOError implies a synchronization between
        //   MPI threads which is required for algorithmic reasons. If you
        //   remove this line, you should replace it by an mpi barrier.
        plbIOError(errorFlag, std::string("Unsuccessful reading from file ") + fName.get());
    }
}

void loadRawData(
    FileName fName, std::vector<plint> const &myBlockIds, std::vector<plint> const &offset,
    std::vector<std::vector<char> > &data)
{
    PLB_ASSERT(myBlockIds.size() == data.size());
    fName.defaultPath(global::directories().getInputDir());
    fName.defaultExt("dat");
    if (global::IOpolicy().useParallelIO() && global::mpi().getSize() > 1) {
        loadRawData_mpi(fName, myBlockIds, offset, data);
    } else {
        // Works in parallel too, but has no parallel efficiency.
        loadRawData_posix(fName, myBlockIds, offset, data);
    }
}

}  // namespace parallelIO

}  // namespace plb
