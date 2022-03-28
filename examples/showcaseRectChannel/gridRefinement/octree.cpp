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

#include "gridRefinement/octree.h"

#include <vector>

#define LDB 0
#define LDF 1
#define LUB 2
#define LUF 3
#define RDB 4
#define RDF 5
#define RUB 6
#define RUF 7

#define LD 8
#define LU 9
#define LB 10
#define LF 11
#define RD 12
#define RU 13
#define RB 14
#define RF 15
#define DB 16
#define DF 17
#define UB 18
#define UF 19

#define L 20
#define R 21
#define D 22
#define U 23
#define B 24
#define F 25

#define UNDEF -1
#define O     UNDEF

#define BRDR -2
#define SMLR -3

namespace plb {

int OctreeTables::surface0N()
{
    return (L);
}
int OctreeTables::surface0P()
{
    return (R);
}
int OctreeTables::surface1N()
{
    return (D);
}
int OctreeTables::surface1P()
{
    return (U);
}
int OctreeTables::surface2N()
{
    return (B);
}
int OctreeTables::surface2P()
{
    return (F);
}

int OctreeTables::edge0NN()
{
    return (DB);
}
int OctreeTables::edge0NP()
{
    return (DF);
}
int OctreeTables::edge0PN()
{
    return (UB);
}
int OctreeTables::edge0PP()
{
    return (UF);
}
int OctreeTables::edge1NN()
{
    return (LB);
}
int OctreeTables::edge1NP()
{
    return (RB);
}
int OctreeTables::edge1PN()
{
    return (LF);
}
int OctreeTables::edge1PP()
{
    return (RF);
}
int OctreeTables::edge2NN()
{
    return (LD);
}
int OctreeTables::edge2NP()
{
    return (LU);
}
int OctreeTables::edge2PN()
{
    return (RD);
}
int OctreeTables::edge2PP()
{
    return (RU);
}

int OctreeTables::cornerNNN()
{
    return (LDB);
}
int OctreeTables::cornerNNP()
{
    return (LDF);
}
int OctreeTables::cornerNPN()
{
    return (LUB);
}
int OctreeTables::cornerNPP()
{
    return (LUF);
}
int OctreeTables::cornerPNN()
{
    return (RDB);
}
int OctreeTables::cornerPNP()
{
    return (RDF);
}
int OctreeTables::cornerPPN()
{
    return (RUB);
}
int OctreeTables::cornerPPP()
{
    return (RUF);
}

int OctreeTables::border()
{
    return (BRDR);
}
int OctreeTables::smaller()
{
    return (SMLR);
}

int OctreeTables::adj(int dir, int oct)
{
    return (adjTab[dir][oct]);
}
int OctreeTables::reflect(int dir, int oct)
{
    return (reflTab[dir][oct]);
}
int OctreeTables::commonFace(int dir, int oct)
{
    return (comFaceTab[dir][oct]);
}
int OctreeTables::commonEdge(int dir, int oct)
{
    return (comEdgeTab[dir][oct]);
}
int OctreeTables::undef()
{
    return (UNDEF);
}

const int OctreeTables::adjTab[26][8] = {
    {1, 0, 0, 0, 0, 0, 0, 0}, {0, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 0, 0, 0, 0, 0},
    {0, 0, 0, 1, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 0, 0, 0}, {0, 0, 0, 0, 0, 1, 0, 0},
    {0, 0, 0, 0, 0, 0, 1, 0}, {0, 0, 0, 0, 0, 0, 0, 1},

    {1, 1, 0, 0, 0, 0, 0, 0}, {0, 0, 1, 1, 0, 0, 0, 0}, {1, 0, 1, 0, 0, 0, 0, 0},
    {0, 1, 0, 1, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 1, 0, 0}, {0, 0, 0, 0, 0, 0, 1, 1},
    {0, 0, 0, 0, 1, 0, 1, 0}, {0, 0, 0, 0, 0, 1, 0, 1}, {1, 0, 0, 0, 1, 0, 0, 0},
    {0, 1, 0, 0, 0, 1, 0, 0}, {0, 0, 1, 0, 0, 0, 1, 0}, {0, 0, 0, 1, 0, 0, 0, 1},

    {1, 1, 1, 1, 0, 0, 0, 0}, {0, 0, 0, 0, 1, 1, 1, 1}, {1, 1, 0, 0, 1, 1, 0, 0},
    {0, 0, 1, 1, 0, 0, 1, 1}, {1, 0, 1, 0, 1, 0, 1, 0}, {0, 1, 0, 1, 0, 1, 0, 1}};

const int OctreeTables::reflTab[26][8] = {
    {RUF, RUB, RDF, RDB, LUF, LUB, LDF, LDB}, {RUF, RUB, RDF, RDB, LUF, LUB, LDF, LDB},
    {RUF, RUB, RDF, RDB, LUF, LUB, LDF, LDB}, {RUF, RUB, RDF, RDB, LUF, LUB, LDF, LDB},
    {RUF, RUB, RDF, RDB, LUF, LUB, LDF, LDB}, {RUF, RUB, RDF, RDB, LUF, LUB, LDF, LDB},
    {RUF, RUB, RDF, RDB, LUF, LUB, LDF, LDB}, {RUF, RUB, RDF, RDB, LUF, LUB, LDF, LDB},

    {RUB, RUF, RDB, RDF, LUB, LUF, LDB, LDF}, {RUB, RUF, RDB, RDF, LUB, LUF, LDB, LDF},
    {RDF, RDB, RUF, RUB, LDF, LDB, LUF, LUB}, {RDF, RDB, RUF, RUB, LDF, LDB, LUF, LUB},
    {RUB, RUF, RDB, RDF, LUB, LUF, LDB, LDF}, {RUB, RUF, RDB, RDF, LUB, LUF, LDB, LDF},
    {RDF, RDB, RUF, RUB, LDF, LDB, LUF, LUB}, {RDF, RDB, RUF, RUB, LDF, LDB, LUF, LUB},
    {LUF, LUB, LDF, LDB, RUF, RUB, RDF, RDB}, {LUF, LUB, LDF, LDB, RUF, RUB, RDF, RDB},
    {LUF, LUB, LDF, LDB, RUF, RUB, RDF, RDB}, {LUF, LUB, LDF, LDB, RUF, RUB, RDF, RDB},

    {RDB, RDF, RUB, RUF, LDB, LDF, LUB, LUF}, {RDB, RDF, RUB, RUF, LDB, LDF, LUB, LUF},
    {LUB, LUF, LDB, LDF, RUB, RUF, RDB, RDF}, {LUB, LUF, LDB, LDF, RUB, RUF, RDB, RDF},
    {LDF, LDB, LUF, LUB, RDF, RDB, RUF, RUB}, {LDF, LDB, LUF, LUB, RDF, RDB, RUF, RUB}};

const int OctreeTables::comFaceTab[26][8] = {
    {O, O, O, L, O, D, B, O}, {O, O, L, O, D, O, O, F}, {O, L, O, O, B, O, O, U},
    {L, O, O, O, O, F, U, O}, {O, D, B, O, O, O, O, R}, {D, O, O, F, O, O, R, O},
    {B, O, O, U, O, R, O, O}, {O, F, U, O, R, O, O, O},

    {O, O, L, L, D, D, O, O}, {L, L, O, O, O, O, U, U}, {O, L, O, L, B, O, B, O},
    {L, O, L, O, O, F, O, F}, {D, D, O, O, O, O, R, R}, {O, O, U, U, R, R, O, O},
    {B, O, B, O, O, R, O, R}, {O, F, O, F, R, O, R, O}, {O, D, B, O, O, D, B, O},
    {D, O, O, F, D, O, O, F}, {B, O, O, U, B, O, O, U}, {O, F, U, O, O, F, U, O},

    {O, O, O, O, O, O, O, O}, {O, O, O, O, O, O, O, O}, {O, O, O, O, O, O, O, O},
    {O, O, O, O, O, O, O, O}, {O, O, O, O, O, O, O, O}, {O, O, O, O, O, O, O, O}};

const int OctreeTables::comEdgeTab[26][8] = {
    {O, LD, LB, O, DB, O, O, O}, {LD, O, O, LF, O, DF, O, O}, {LB, O, O, LU, O, O, UB, O},
    {O, LF, LU, O, O, O, O, UF}, {DB, O, O, O, O, RD, RB, O}, {O, DF, O, O, RD, O, O, RF},
    {O, O, UB, O, RB, O, O, RU}, {O, O, O, UF, O, RF, RU, O},

    {O, O, O, O, O, O, O, O},    {O, O, O, O, O, O, O, O},    {O, O, O, O, O, O, O, O},
    {O, O, O, O, O, O, O, O},    {O, O, O, O, O, O, O, O},    {O, O, O, O, O, O, O, O},
    {O, O, O, O, O, O, O, O},    {O, O, O, O, O, O, O, O},    {O, O, O, O, O, O, O, O},
    {O, O, O, O, O, O, O, O},    {O, O, O, O, O, O, O, O},    {O, O, O, O, O, O, O, O},

    {O, O, O, O, O, O, O, O},    {O, O, O, O, O, O, O, O},    {O, O, O, O, O, O, O, O},
    {O, O, O, O, O, O, O, O},    {O, O, O, O, O, O, O, O},    {O, O, O, O, O, O, O, O}};

}  // namespace plb
