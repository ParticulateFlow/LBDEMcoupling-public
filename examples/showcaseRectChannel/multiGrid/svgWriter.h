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
 * SVG writer for refined
 */
#ifndef SVG_WRITER_2D_H
#define SVG_WRITER_2D_H

#include <string>

#include "core/util.h"
#include "multiGrid/multiGridManagement2D.h"

namespace plb {

class SVGWriter2D {
public:
    SVGWriter2D(MultiGridManagement2D management_) : management(management_) { }

    void writePatterns()
    {
        plint numLevels = management.getNumLevels();
        out << "<defs>\n";
        for (plint iLevel = 0; iLevel < numLevels; ++iLevel) {
            double size = 20.0 * util::twoToThePower(-iLevel);
            out << "<pattern id=\"Pat" << size << "\" width=\"" << size << "\" height=\"" << size
                << "\" patternUnits=\"userSpaceOnUse\" >\n"
                << "\t<rect width=\"" << size << "\" height=\"" << size
                << "\" fill=\"none\" stroke=\"#000000\" stroke-width=\"2.0\" />\n</pattern>";
            out << std::endl;
        }
        out << "</defs>\n";
    }

    void writeDomainsWithDynamicsInfo(
        std::string fileName, int dynamicsNumber,
        std::vector<MultiScalarField2D<int> > dynamicsInfo,
        std::vector<std::map<int, std::string> > &idToName, std::map<std::string, int> &nameToColor)
    {
        out.open(fileName.c_str());
        Box2D finestBoundingBox = management.getBoundingBox(management.getNumLevels() - 1);
        Dot2D size(finestBoundingBox.getNx(), finestBoundingBox.getNy());
        writeHeader();
        writePatterns();

        for (plint iLevel = 0; iLevel < management.getNumLevels(); ++iLevel) {
            // write each block of the current level
            std::vector<Box2D> bulks = management.getBulks(iLevel);
            double size = util::twoToThePower(-iLevel);
            for (pluint iBlock = 0; iBlock < bulks.size(); ++iBlock) {
                Box2D currentBlock = bulks[iBlock];
                writeBG(currentBlock, size * 20.0);
                for (plint iX = currentBlock.x0; iX <= currentBlock.x1; ++iX) {
                    for (plint iY = currentBlock.y0; iY <= currentBlock.y1; ++iY) {
                        // std::string color = colors[dynamicsInfo[iLevel].get(iX,iY)];
                        int id = nameToColor[idToName[iLevel][dynamicsInfo[iLevel].get(iX, iY)]];
                        if (id != -1) {
                            Array<double, 2> location(
                                (double)iX * size * 20.0, (double)iY * size * 20.0);
                            if (id >= 128)
                                id = 127;
                            // writeRectangle(location, size, id);
                            //                                     writeCircle(location, size*10.0,
                            //                                     id);
                            // writeCircle(location, 7.0*size, id);
                        }
                    }
                }
            }
        }
        writeEnd();
    }

    void writeDomainsWithDynamicsInfo(
        std::string fileName, int dynamicsNumber, std::vector<std::vector<Box2D> > blocks,
        std::vector<MultiScalarField2D<int> > dynamicsInfo,
        std::vector<std::map<int, std::string> > &idToName, std::map<std::string, int> &nameToColor)
    {
        out.open(fileName.c_str());
        Box2D finestBoundingBox = management.getBoundingBox(management.getNumLevels() - 1);
        Dot2D size(finestBoundingBox.getNx(), finestBoundingBox.getNy());
        writeHeader();
        writePatterns();

        for (plint iLevel = 0; iLevel < management.getNumLevels(); ++iLevel) {
            // write each block of the current level
            std::vector<Box2D> bulks = blocks[iLevel];
            double size = util::twoToThePower(-iLevel);
            for (pluint iBlock = 0; iBlock < bulks.size(); ++iBlock) {
                Box2D currentBlock = bulks[iBlock];
                writeBG(currentBlock, size * 20.0);
                for (plint iX = currentBlock.x0; iX <= currentBlock.x1; ++iX) {
                    for (plint iY = currentBlock.y0; iY <= currentBlock.y1; ++iY) {
                        // std::string color = colors[dynamicsInfo[iLevel].get(iX,iY)];
                        int id = nameToColor[idToName[iLevel][dynamicsInfo[iLevel].get(iX, iY)]];
                        if (id != -1) {
                            Array<double, 2> location(
                                (double)iX * size * 20.0, (double)iY * size * 20.0);
                        }
                    }
                }
            }
        }
        writeEnd();
    }

    void writeBG(Box2D block, double size)
    {
        double height = block.y1 - block.y0 == 0 ? -size : (double)(block.y1 - block.y0) * size;
        double width = block.x1 - block.x0 == 0 ? -size : (double)(block.x1 - block.x0) * size;
        double x0 = (double)block.x0 * size;
        double y0 = (double)block.y0 * size;
        out << "<rect x=\"" << x0 << "\" y=\"" << y0 << "\" width=\"" << width << "\" height=\""
            << height << "\" "
            << "style=\"stroke:black;fill:url(#Pat" << size << ")\"/>\n";
    }

private:
    //         void writeHeader(Dot2D size){
    //             out << "<?xml version=\"1.0\"?>\n<svg xmlns=\"http://www.w3.org/2000/svg\"
    //             width=\""
    //                 << size.x << "cm\"" << "height=\""<< size.y << "cm\">\n";
    //             out << "\t<g style=\"fill-opacity:0.7; stroke:black; stroke-width:0.1cm;\">" <<
    //             std::endl;
    //         }

    void writeHeader()
    {
        out << "<?xml version=\"1.0\" standalone=\"no\"?>\n";
        out << "<!DOCTYPE svg PUBLIC"
            << " \"-//W3C//DTD SVG 1.1//EN\" "
               "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n";
        out << "\t<svg width=\"100%\" height=\"100%\" version=\"1.1\" "
               "xmlns=\"http://www.w3.org/2000/svg\">";
        out << std::endl;
    }

    /*void writeCircle(Array<double,2> location, double radius, int color){
        PLB_ASSERT( color>=0 && color <128 );
        out << "<circle cx=\"" << location[0] << "\" cy=\"" << location[1]
            << "\" r=\"" << radius << "\" style=\"stroke:none; " <<  "fill-opacity:0.5; "
            <<"fill:" << "rgb(" << colors[color][0] << ","
            << colors[color][1] << "," << colors[color][2]  << ")" << "\"/>" << std::endl;

    }*/

    void writeEnd()
    {
        out << "</svg>" << std::endl;
    }

    MultiGridManagement2D management;
    plb_ofstream out;
};

}  // namespace plb

#endif  // SVG_WRITER_2D_H
