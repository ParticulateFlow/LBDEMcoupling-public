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

#ifndef PARALLEL_IO_H
#define PARALLEL_IO_H

#include <cstdio>
#include <fstream>
#include <iostream>
#include <istream>
#include <ostream>
#include <streambuf>

#include "core/globalDefs.h"
#include "parallelism/mpiManager.h"

namespace plb {

/// A buffer which reads and writes nothing.
/** This buffer is designed for use in parallel programming.
 *  It is assigned to processors which have potentially no access
 *  to the file system.
 */
class DevNullBuffer : public std::streambuf {
protected:
    virtual int_type overflow(int_type c)
    {
        return EOF;
    }
    virtual int_type underflow()
    {
        return EOF;
    }
};

// Parallel Output Streams.

struct Parallel_ostream {
    virtual ~Parallel_ostream() { }
    virtual std::ostream &getOriginalStream() = 0;
};

class Parallel_referring_ostream : public Parallel_ostream {
public:
    Parallel_referring_ostream(std::ostream &original_ostream_) :
        devNullStream(&devNullBuffer), original_ostream(original_ostream_)
    { }
    virtual std::ostream &getOriginalStream()
    {
        if (global::mpi().isMainProcessor()) {
            return original_ostream;
        } else {
            return devNullStream;
        }
    }

private:
    DevNullBuffer devNullBuffer;
    std::ostream devNullStream;
    std::ostream &original_ostream;
};

template <typename Value>
Parallel_ostream &operator<<(Parallel_ostream &lhs, Value const &rhs)
{
    lhs.getOriginalStream() << rhs;
    return lhs;
}

inline Parallel_ostream &operator<<(Parallel_ostream &lhs, std::ostream &(*op)(std::ostream &))
{
    lhs.getOriginalStream() << op;
    return lhs;
}

class plb_ofstream : public Parallel_ostream {
public:
    plb_ofstream();
    explicit plb_ofstream(
        const char *filename,
        std::ostream::openmode mode = std::ostream::out | std::ostream::trunc);
    ~plb_ofstream();
    virtual std::ostream &getOriginalStream();

    bool is_open();
    plint tellp();
    void open(
        const char *filename,
        std::ostream::openmode mode = std::ostream::out | std::ostream::trunc);
    void close();

private:
    plb_ofstream(plb_ofstream const &rhs);
    plb_ofstream &operator=(plb_ofstream const &rhs);

private:
    DevNullBuffer devNullBuffer;
    std::ostream devNullStream;
    std::ofstream *original;
};

extern Parallel_referring_ostream pcout;
extern Parallel_referring_ostream pcerr;
extern Parallel_referring_ostream pclog;

// Parallel Input Streams.

struct Parallel_istream {
    virtual ~Parallel_istream() { }
    virtual std::istream &getOriginalStream() = 0;
};

class Parallel_referring_istream : public Parallel_istream {
public:
    Parallel_referring_istream(std::istream &original_istream_) :
        original_istream(original_istream_)
    { }
    virtual std::istream &getOriginalStream()
    {
        return original_istream;
    }

private:
    std::istream &original_istream;
};

template <typename Value>
Parallel_istream &operator>>(Parallel_istream &lhs, Value &rhs)
{
    lhs.getOriginalStream() >> rhs;
    return lhs;
}

inline Parallel_istream &operator>>(Parallel_istream &lhs, std::istream &(*op)(std::istream &))
{
    lhs.getOriginalStream() >> op;
    return lhs;
}

class plb_ifstream : public Parallel_istream {
public:
    plb_ifstream();
    explicit plb_ifstream(const char *filename, std::istream::openmode mode = std::ostream::in);
    ~plb_ifstream();
    virtual std::istream &getOriginalStream();

    bool is_open();
    void open(const char *filename, std::istream::openmode mode = std::ostream::in);
    void close();
    bool good();

private:
    plb_ifstream(plb_ifstream const &rhs);
    plb_ifstream &operator=(plb_ifstream const &rhs);

private:
    DevNullBuffer devNullBuffer;
    std::istream devNullStream;
    std::ifstream *original;
};

// General utility functions

void plbIOErrorIfCannotOpenFileForReading(std::string fileName);
void plbIOErrorIfCanOpenFileForReading(std::string fileName);
// Caution: The directory name "dirName" must include the separator (/ for Unix-like, or \ for
// Windows).
//          The file created is opened with the "w" mode, so it is erased if it already exists and
//          then removed.
void plbIOErrorIfCannotCreateFileInDir(std::string dirName, std::string fileName);
void plbIOErrorIfFileErrorOccurred(FILE *fp);

void plbMainProcIOErrorIfCannotOpenFileForReading(std::string fileName);
void plbMainProcIOErrorIfCanOpenFileForReading(std::string fileName);
// Caution: The directory name "dirName" must include the separator (/ for Unix-like, or \ for
// Windows).
//          The file created is opened with the "w" mode, so it is erased if it already exists and
//          then removed.
void plbMainProcIOErrorIfCannotCreateFileInDir(std::string dirName, std::string fileName);
void plbMainProcIOErrorIfFileErrorOccurred(FILE *fp);

void abortIfCannotOpenFileForReading(std::string fileName);
void abortIfCanOpenFileForReading(std::string fileName);
// Caution: The directory name "dirName" must include the separator (/ for Unix-like, or \ for
// Windows).
//          The file created is opened with the "w" mode, so it is erased if it already exists and
//          then removed.
void abortIfCannotCreateFileInDir(std::string dirName, std::string fileName);
void abortIfFileErrorOccurred(FILE *fp);

void abortIfCannotOpenFileForReadingAtMainProc(std::string fileName);
void abortIfCanOpenFileForReadingAtMainProc(std::string fileName);
// Caution: The directory name "dirName" must include the separator (/ for Unix-like, or \ for
// Windows).
//          The file created is opened with the "w" mode, so it is erased if it already exists and
//          then removed.
void abortIfCannotCreateFileInDirAtMainProc(std::string dirName, std::string fileName);
void abortIfFileErrorOccurredAtMainProc(FILE *fp);

void makeDirectory(std::string dirName, bool abortIfExists = true);

}  // namespace plb

#endif
