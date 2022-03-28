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
 * A timer class for benchmarking program parts -- header file.
 */
#ifndef PLB_FILES_H
#define PLB_FILES_H

#include <string>

namespace plb {

class FileName {
public:
    FileName() { }
    FileName(const char *file);
    FileName(std::string file);
    std::string get() const;
    FileName &setPath(std::string path_);
    FileName &setName(std::string name_);
    FileName &setExt(std::string ext_);
    FileName &defaultPath(std::string path_);
    FileName &defaultName(std::string name_);
    FileName &defaultExt(std::string ext_);
    std::string getPath() const
    {
        return path;
    }
    std::string getName() const
    {
        return name;
    }
    std::string getExt() const
    {
        return ext;
    }
    operator std::string() const
    {
        return get();
    }

private:
    void initialize(std::string file);

private:
    std::string path, name, ext;
};

}  // namespace plb

#endif  // PLB_FILES_H
