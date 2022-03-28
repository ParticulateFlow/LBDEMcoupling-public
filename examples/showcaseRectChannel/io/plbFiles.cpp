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

#include "io/plbFiles.h"

#include "core/plbDebug.h"

namespace plb {

FileName::FileName(const char *file)
{
    initialize(file);
}

FileName::FileName(std::string file)
{
    initialize(file);
}

void FileName::initialize(std::string file)
{
    std::string fullName;
    size_t sep = file.find_last_of("\\/");
    if (sep != std::string::npos) {
        fullName = file.substr(sep + 1, file.size() - sep - 1);
        path = file.substr(0, sep);
    } else {
        fullName = file;
    }

    size_t dot = fullName.find_last_of(".");
    if (dot != std::string::npos) {
        name = fullName.substr(0, dot);
        ext = fullName.substr(dot + 1, fullName.size() - dot - 1);
    } else {
        name = fullName;
        ext = "";
    }
}

std::string FileName::get() const
{
    std::string composite;
    if (path != "") {
        // First determine the separator.
        std::string separator;
        size_t sep = path.find_last_of("/");
        if (sep != std::string::npos) {
            // Unix separator.
            separator = "/";
        } else {
            sep = path.find_last_of("\\");
            if (sep != std::string::npos) {
                // Windows separator.
                separator = "\\";
            } else {
                // By default we choose the Unix separator.
                // This needs to change, and to find a better
                // way to detect when the separator is
                // "/" or "\".
                separator = "/";
            }
        }

        composite = path + separator;
    }
    composite += name;
    if (ext != "") {
        composite += "." + ext;
    }
    return composite;
}

FileName &FileName::setPath(std::string path_)
{
    path = path_;
    if (!path.empty() && (path[path.size() - 1] == '/' || path[path.size() - 1] == '\\')) {
        path.erase(path.end() - 1);
    }
    return *this;
}

FileName &FileName::setName(std::string name_)
{
    name = name_;
    return *this;
}

FileName &FileName::setExt(std::string ext_)
{
    ext = ext_;
    if (!ext.empty() && ext[0] == '.') {
        ext.erase(ext.begin());
    }
    return *this;
}

FileName &FileName::defaultPath(std::string path_)
{
    if (path == "") {
        return setPath(path_);
    } else {
        return *this;
    }
}

FileName &FileName::defaultName(std::string name_)
{
    if (name == "") {
        return setName(name_);
    } else {
        return *this;
    }
}

FileName &FileName::defaultExt(std::string ext_)
{
    if (ext == "") {
        return setExt(ext_);
    } else {
        return *this;
    }
}

}  // namespace plb
