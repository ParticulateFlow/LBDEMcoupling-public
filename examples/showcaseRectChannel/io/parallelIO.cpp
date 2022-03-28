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

#include "io/parallelIO.h"

#include <cerrno>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

#include "core/globalDefs.h"
#include "core/runTimeDiagnostics.h"
#include "core/util.h"
#include "parallelism/mpiManager.h"

#ifdef PLB_USE_POSIX
#include <sys/stat.h>
#include <sys/types.h>
#else
#ifdef PLB_WINDOWS
#include <direct.h>
#endif
#endif

namespace plb {

Parallel_referring_ostream pcout(std::cout);
Parallel_referring_ostream pcerr(std::cerr);
Parallel_referring_ostream pclog(std::clog);

/* *************** Class plb_ofstream ******************************** */

plb_ofstream::plb_ofstream() :
    devNullStream(&devNullBuffer), original(global::mpi().isMainProcessor() ? new std::ofstream : 0)
{ }

plb_ofstream::plb_ofstream(const char *filename, std::ostream::openmode mode) :
    devNullStream(&devNullBuffer),
    original(global::mpi().isMainProcessor() ? new std::ofstream(filename, mode) : 0)
{ }

// QUESTION: Why is the copy const and equal doing nothing?
plb_ofstream::plb_ofstream(plb_ofstream const &rhs) : devNullStream(&devNullBuffer), original(0) { }

plb_ofstream &plb_ofstream::operator=(plb_ofstream const &rhs)
{
    return *this;
}

plb_ofstream::~plb_ofstream()
{
    delete original;
}

std::ostream &plb_ofstream::getOriginalStream()
{
    if (global::mpi().isMainProcessor()) {
        return *original;
    } else {
        return devNullStream;
    }
}

bool plb_ofstream::is_open()
{
#ifdef PLB_MPI_PARALLEL
    int open = false;
    if (global::mpi().isMainProcessor()) {
        open = original->is_open();
    }
    global::mpi().bCast(&open, 1);
    return open;
#else
    return original->is_open();
#endif
}

plint plb_ofstream::tellp()
{
#ifdef PLB_MPI_PARALLEL
    plint pos = -1;
    if (global::mpi().isMainProcessor()) {
        pos = original->tellp();
    }
    global::mpi().bCast(&pos, 1);
    return pos;
#else
    return original->tellp();
#endif
}

void plb_ofstream::open(const char *filename, std::ostream::openmode mode)
{
    if (global::mpi().isMainProcessor()) {
        original->open(filename, mode);
    }
}

void plb_ofstream::close()
{
    if (global::mpi().isMainProcessor()) {
        original->close();
    }
}

/* *************** Class plb_ifstream ******************************** */

plb_ifstream::plb_ifstream() :
    devNullStream(&devNullBuffer), original(global::mpi().isMainProcessor() ? new std::ifstream : 0)
{ }

plb_ifstream::plb_ifstream(const char *filename, std::istream::openmode mode) :
    devNullStream(&devNullBuffer),
    original(global::mpi().isMainProcessor() ? new std::ifstream(filename, mode) : 0)
{ }

// QUESTION: Why copy and equal are doing nothing?
plb_ifstream::plb_ifstream(plb_ifstream const &rhs) : devNullStream(&devNullBuffer), original(0) { }

plb_ifstream &plb_ifstream::operator=(plb_ifstream const &rhs)
{
    return *this;
}

plb_ifstream::~plb_ifstream()
{
    delete original;
}

std::istream &plb_ifstream::getOriginalStream()
{
    if (global::mpi().isMainProcessor()) {
        return *original;
    } else {
        return devNullStream;
    }
}

bool plb_ifstream::is_open()
{
#ifdef PLB_MPI_PARALLEL
    int open = false;
    if (global::mpi().isMainProcessor()) {
        open = original->is_open();
    }
    global::mpi().bCast(&open, 1);
    return open;
#else
    return original->is_open();
#endif
}

bool plb_ifstream::good()
{
#ifdef PLB_MPI_PARALLEL
    int open = false;
    if (global::mpi().isMainProcessor()) {
        open = original->good();
    }
    global::mpi().bCast(&open, 1);
    return open;
#else
    return original->good();
#endif
}

void plb_ifstream::open(const char *filename, std::istream::openmode mode)
{
    if (global::mpi().isMainProcessor()) {
        original->open(filename, mode);
    }
}

void plb_ifstream::close()
{
    if (global::mpi().isMainProcessor()) {
        original->close();
    }
}

/* *************** General utility functions ******************************** */

void plbIOErrorIfCannotOpenFileForReading(std::string fileName)
{
    bool issueError = false;
    std::string message = fileName + ": No error occurred in this process";
    FILE *fp = fopen(fileName.c_str(), "rb");
    if (fp == 0) {
        issueError = true;
        message = fileName + ": Cannot open file for reading; " + strerror(errno);
        global::plbErrors().registerIOError(message);
    } else {
        fclose(fp);
    }
    plbIOError(issueError, message);
}

void plbIOErrorIfCanOpenFileForReading(std::string fileName)
{
    bool issueError = false;
    std::string message = fileName + ": No error occurred in this process";
    FILE *fp = fopen(fileName.c_str(), "rb");
    if (fp != 0) {
        issueError = true;
        message = fileName + ": Can open file for reading";
        global::plbErrors().registerIOError(message);
        fclose(fp);
    }
    plbIOError(issueError, message);
}

void plbIOErrorIfCannotCreateFileInDir(std::string dirName, std::string dummyFileName)
{
    bool issueError = false;
    std::string message = dirName + ": No error occurred in this process";
    // Caution: The directory name "dirName" must include the separator (/ for Unix-like, or \ for
    // Windows).
    //          The file created is opened with the "w" mode, so it is erased if it already exists
    //          and then removed.
    std::string fileName = dirName + dummyFileName + "_" + util::val2str(global::mpi().getRank());
    FILE *fp = fopen(fileName.c_str(), "w");
    if (fp == 0) {
        issueError = true;
        message =
            dirName + ": Cannot create a file for writing in this directory; " + strerror(errno);
        global::plbErrors().registerIOError(message);
    } else {
        fclose(fp);
        remove(fileName.c_str());
    }
    plbIOError(issueError, message);
}

void plbIOErrorIfFileErrorOccurred(FILE *fp)
{
    bool issueError = false;
    std::string message = "No error occurred in this process";
    global::mpi().barrier();  // Normally we should lock the file, but there is no portable way of
                              // doing so...
    if (ferror(fp)) {
        issueError = true;
        message = "File stream error; " + std::string(strerror(errno));
        global::plbErrors().registerIOError(message);
    }
    plbIOError(issueError, message);
}

void plbMainProcIOErrorIfCannotOpenFileForReading(std::string fileName)
{
    bool issueError = false;
    std::string message = fileName + ": No error occurred in this process";
    if (global::mpi().isMainProcessor()) {
        FILE *fp = fopen(fileName.c_str(), "rb");
        if (fp == 0) {
            issueError = true;
            message = fileName + ": Cannot open file for reading; " + strerror(errno);
            global::plbErrors().registerIOError(message);
        } else {
            fclose(fp);
        }
    }
    plbMainProcIOError(issueError, message);
}

void plbMainProcIOErrorIfCanOpenFileForReading(std::string fileName)
{
    bool issueError = false;
    std::string message = fileName + ": No error occurred in this process";
    if (global::mpi().isMainProcessor()) {
        FILE *fp = fopen(fileName.c_str(), "rb");
        if (fp != 0) {
            issueError = true;
            message = fileName + ": Can open file for reading";
            global::plbErrors().registerIOError(message);
            fclose(fp);
        }
    }
    plbMainProcIOError(issueError, message);
}

void plbMainProcIOErrorIfCannotCreateFileInDir(std::string dirName, std::string dummyFileName)
{
    bool issueError = false;
    std::string message = dirName + ": No error occurred in this process";
    if (global::mpi().isMainProcessor()) {
        // Caution: The directory name "dirName" must include the separator (/ for Unix-like, or
        // \ for Windows).
        //          The file created is opened with the "w" mode, so it is erased if it already
        //          exists and then removed.
        std::string fileName = dirName + dummyFileName;
        FILE *fp = fopen(fileName.c_str(), "w");
        if (fp == 0) {
            issueError = true;
            message = dirName + ": Cannot create a file for writing in this directory; "
                      + strerror(errno);
            global::plbErrors().registerIOError(message);
        } else {
            fclose(fp);
            remove(fileName.c_str());
        }
    }
    plbMainProcIOError(issueError, message);
}

void plbMainProcIOErrorIfFileErrorOccurred(FILE *fp)
{
    bool issueError = false;
    std::string message = "No error occurred in this process";
    global::mpi().barrier();  // Normally we should lock the file, but there is no portable way of
                              // doing so...
    if (global::mpi().isMainProcessor()) {
        if (ferror(fp)) {
            issueError = true;
            message = "File stream error; " + std::string(strerror(errno));
            global::plbErrors().registerIOError(message);
        }
    }
    plbMainProcIOError(issueError, message);
}

void abortIfCannotOpenFileForReading(std::string fileName)
{
    try {
        plbIOErrorIfCannotOpenFileForReading(fileName);
    } catch (PlbException &exception) {
        std::string message = "Caught exception in process "
                              + util::val2str(global::mpi().getRank()) + ". " + exception.what()
                              + "\n";
        printSerially(stderr, message);
        exit(-1);
    }
}

void abortIfCanOpenFileForReading(std::string fileName)
{
    try {
        plbIOErrorIfCanOpenFileForReading(fileName);
    } catch (PlbException &exception) {
        std::string message = "Caught exception in process "
                              + util::val2str(global::mpi().getRank()) + ". " + exception.what()
                              + "\n";
        printSerially(stderr, message);
        exit(-1);
    }
}

void abortIfCannotCreateFileInDir(std::string dirName, std::string dummyFileName)
{
    try {
        plbIOErrorIfCannotCreateFileInDir(dirName, dummyFileName);
    } catch (PlbException &exception) {
        std::string message = "Caught exception in process "
                              + util::val2str(global::mpi().getRank()) + ". " + exception.what()
                              + "\n";
        printSerially(stderr, message);
        exit(-1);
    }
}

void abortIfFileErrorOccurred(FILE *fp)
{
    try {
        plbIOErrorIfFileErrorOccurred(fp);
    } catch (PlbException &exception) {
        std::string message = "Caught exception in process "
                              + util::val2str(global::mpi().getRank()) + ". " + exception.what()
                              + "\n";
        printSerially(stderr, message);
        exit(-1);
    }
}

void abortIfCannotOpenFileForReadingAtMainProc(std::string fileName)
{
    try {
        plbMainProcIOErrorIfCannotOpenFileForReading(fileName);
    } catch (PlbException &exception) {
        std::string message = "Caught exception in process "
                              + util::val2str(global::mpi().getRank()) + ". " + exception.what()
                              + "\n";
        printSerially(stderr, message);
        exit(-1);
    }
}

void abortIfCanOpenFileForReadingAtMainProc(std::string fileName)
{
    try {
        plbMainProcIOErrorIfCanOpenFileForReading(fileName);
    } catch (PlbException &exception) {
        std::string message = "Caught exception in process "
                              + util::val2str(global::mpi().getRank()) + ". " + exception.what()
                              + "\n";
        printSerially(stderr, message);
        exit(-1);
    }
}

void abortIfCannotCreateFileInDirAtMainProc(std::string dirName, std::string dummyFileName)
{
    try {
        plbMainProcIOErrorIfCannotCreateFileInDir(dirName, dummyFileName);
    } catch (PlbException &exception) {
        std::string message = "Caught exception in process "
                              + util::val2str(global::mpi().getRank()) + ". " + exception.what()
                              + "\n";
        printSerially(stderr, message);
        exit(-1);
    }
}

void abortIfFileErrorOccurredAtMainProc(FILE *fp)
{
    try {
        plbMainProcIOErrorIfFileErrorOccurred(fp);
    } catch (PlbException &exception) {
        std::string message = "Caught exception in process "
                              + util::val2str(global::mpi().getRank()) + ". " + exception.what()
                              + "\n";
        printSerially(stderr, message);
        exit(-1);
    }
}

void makeDirectory(std::string dirName, bool abortIfExists)
{
    bool issueError = false;
    std::string message = dirName + ": No error occurred in this process";
    if (global::mpi().isMainProcessor()) {
        bool useErrno = true;
        int err = 0;
#ifdef PLB_USE_POSIX
        err = mkdir(dirName.c_str(), 0777);
#else
#ifdef PLB_WINDOWS
        err = _mkdir(dirName.c_str());
#else
        issueError = true;
        useErrno = false;
        err = 0;
#endif
#endif
        if (err != 0) {
            issueError = true;
            if (errno == EEXIST && !abortIfExists) {
                issueError = false;
            }
        }
        if (issueError) {
            std::string explanation;
            if (useErrno) {
                explanation = strerror(errno);
            } else {
                explanation = " Operating system not supported";
            }
            message = dirName + ": Cannot make directory; " + explanation;
            global::plbErrors().registerIOError(message);
        }
    }
    plbMainProcIOError(issueError, message);
}

}  // namespace plb
