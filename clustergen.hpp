/* Copyright 2023 Daniel Grier and Luke Schaeffer

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License. */

/*
clustergen.hpp - Contains the logic for the divide-and-conquer measurement procedure on the grid
*/
#ifndef _CLUSTERGEN_HPP_
#define _CLUSTERGEN_HPP_

#include <vector>
#include "clifford.hpp"

typedef QUBIT PhysQubit;

class Frame;

class GridPoint {
  public:
    int x;
    int y;
    PhysQubit phys;
    bool isCornerOf(const Frame &f) const;

    GridPoint(int x, int y, PhysQubit phys) : x(x), y(y), phys(phys) { };

    int globalPosition(int len) const {
        return y*len + x;
    }
};

class Frame {
  public:
    int x;
    int y;
    int dx;
    int dy;

    Params range;
    GridPoint first;
    int firstIndexNW;

    Frame(int x, int y, int w, int h) : x(x), y(y), dx(w-1), dy(h-1), range(0,0), first(-1, -1, -1) {};

    void setDimensions(int x, int y, int w, int h) {
        this->x = x; this->y = y;
        dx = w-1; dy = h-1;
    }

    int indexNW(const GridPoint &q) const;
    int index(const GridPoint &q) const;

    void setStart(int x, int y, PhysQubit phys) {
        range.setParams(phys, phys+perimeter());
        first = GridPoint(x, y, phys);
        firstIndexNW = indexNW(first);
    }


    int x1() const {
      return x;
    }
    int y1() const {
      return y;
    }
    int x2() const {
      return x + dx;
    }
    int y2() const {
      return y + dy;
    }
    int perimeter() const {
      return 2*(dx + dy);
    }
    int width() const {
      return dx+1;
    }
    int height() const {
      return dy+1;
    }
    GridPoint begin() const {
        return first;
    }
    PhysQubit offset() const {
        return first.phys;
    }
    GridPoint operator[](int i) const;

    int wrap(int i) const {
        int p = perimeter();
        return ((i % p) + p) % p;
    }
};

// ClusterProg options
const int FIXEDBASES          = 0x100;
const int SPACE_OPT           = 0x200;
const int PRINT_OUTCOMES      = 0x400;
const int ROW_BY_ROW          = 0x800;
const int NAIVE               = 0x1000;
const int VERIFY              = 0x2000;
const int RANDOM_ORDER        = 0x4000;

class ClusterProg: public QProg {
public:
   int len;
   std::vector<Basis> bases;

   void printGrid() const;
   void generateCluster();
   void generateCluster(int len);
   void expand();
   void contract();
   void printOutcomes() const;

   // reads len many outcomes/bases from a file
   void readOutcomes(const std::string &fn);
   void readBases(const std::string &fn);

   ClusterProg(int len) : len(len) {
      setFlag(SPACE_OPT | PRINT_OUTCOMES);
      generateCluster();
   };
   ClusterProg() : len(0) {
      setFlag(SPACE_OPT | PRINT_OUTCOMES);
   };

private:
   int* qmap; // Mapping from qubits in the grid to their physical qubit
   std::vector<int> morder; // The order in which the qubits of the grid are measured

   void addSwap(PhysQubit a, PhysQubit b, const Params &p);
   void addCSIGN(PhysQubit a, PhysQubit b, const Params & p);
   void addHadamard(PhysQubit a); // adds a Hadamard for the purpose of changing qubit a from |0> to |+>
   void swapRange(PhysQubit a, PhysQubit m, PhysQubit b);
   void swapRange(PhysQubit a, PhysQubit m, PhysQubit b, const Params &p);
   void measurePoint(const GridPoint& i, const Params& p, bool reset=false);
   void measurePointAndReset(const GridPoint &q, const Params &p);
   void doFrame(const Frame& f);
   void doOuterFrame(const Frame& f);
};

// Outputs a basis (XBASIS, YBASIS, ZBASIS)
// based on their relative frequencies (xfreq, yfreq, zfreq)
Basis randombasis(double xfreq=XFREQ, double yfreq=YFREQ);

#endif
