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

#include <algorithm>
#include <iostream>
#include "clifford.hpp"
#include "clustergen.hpp"

using namespace std;

Basis randombasis(double xfreq, double yfreq) {
    double s = (double)(rand())/RAND_MAX;
    if(s <= xfreq) return XBASIS;
    if(s <= xfreq + yfreq) return YBASIS;
    return ZBASIS;
}

// If the gates on are the physical qubits, convert them to logical qubits
void ClusterProg::expand() {
    int m = 0;

    // keep track of which gates we've already converted to logical qubits
    int touched[gates.size()][2];
    for(int i=0; i < (int)gates.size(); i++) {
        touched[i][0] = 0; // g.a
        touched[i][1] = 0; // g.b
    }

    for(int i = 0; i < (int)gates.size(); i++) {
        Gate g = gates.at(i);
        if(g.op == MEASURE || g.op == MEASURE_RESET) { // update all previous gates matching this logical qubit
            int phys = g.a;
            int logi = morder.at(m);
            int j = i-1;
            while(j >= 0) {
                Gate& h = gates.at(j);
                if((h.op == MEASURE || h.op == MEASURE_RESET) && h.a == phys) {
                    break;
                } else if(h.op == SWAP) {
                    if(h.a == phys) phys = h.b;
                    else if(h.b == phys) phys = h.a;
                } else {
                    if(h.a == phys && !touched[j][0]) {
                        h.a = logi;
                        touched[j][0] = 1;
                    }
                    if(h.b == phys && !touched[j][1]) {
                        h.b = logi;
                        touched[j][1] = 1;
                    }
                }
                j--;
            }
            m++;
        }
    }
    // fix up the measure gates that weren't changed above
    m = 0;
    for(int i = 0; i < (int)gates.size(); i++) {
        Gate &g = gates.at(i);
        if(g.op == MEASURE || g.op == MEASURE_RESET) {
            g.a = morder.at(m);
            m++;
        } else if(g.op == RESET) { // this is a bit of a hack. convert things to SWAP that we want to kill
            g.op = SWAP; // no need to postselect when we have expanded the circuit and aren't reusing qubits
        }
    }
    // Params are now meaningless
    for(Gate &g : gates) {
        g.p.setParams(0,len*len);
    }
    this->n = len*len;

    // Remove all the SWAP gates
    gates.erase(std::remove(gates.begin(), gates.end(), SWAP), gates.end());
}

// If the gates on are the logical qubits, convert them to physical qubits
// Doesn't currently work if the physical implementation requires SWAP gates
void ClusterProg::contract() {
    for(Gate &g : gates) {
        g.a = qmap[g.a];
        g.b = qmap[g.b];
    }
    int max = 0;
    for(int i = 0; i < len*len; i++) {
        max = MAX(max, qmap[i]);
    }
    this->n = max;
}



/*
Prints the grid, used mainly for debugging
This function only works if space-optimization is off
*/
void ClusterProg::printGrid() const {
    cout << "len: " << len << " n: " << n << endl;
    int n = len*len;

    // this is unnecessarily wasteful with space
    char** A = new char*[n];
    for(int i = 0; i < n; i++) {
        A[i] = new char[n];
    }

    for(int i = 0; i< n; i++) {
        for(int j=0; j < n; j++)
        A[i][j]=0;
    }

    for (const Gate &g : gates)  {
        if (g.op==CZ) {
            A[g.a][g.b] += 1;
            A[g.b][g.a] += 1;
        } else if (g.op==MEASURE) {
            A[g.a][g.a] += 1;
        }
    }

    for(int i = 0; i<len; i++) {
        for(int j=0; j<len-1; j++) {
            int q = i*len + j;
            printf("%d", A[q][q]);
            if(A[q][q+1]) printf("-%d-", A[q][q+1]);
            else printf("   ");
        } printf("%d\n", A[i*len + len -1][i*len + len -1]);
        if(i!=len-1) {
            for(int j=0; j<len; j++) {
                int q = i*len + j;
                if(A[q][q+len]) printf("%d", A[q][q+len]);
                else printf(" ");
                printf("   ");
            }
            printf("\n");
        }
    }
}

ostream& operator<<(ostream& out, const GridPoint& q) {
    out << "(" << q.x << "," << q.y << ")" << "@" << q.phys;
    return out;
}

void ClusterProg::addCSIGN(PhysQubit a, PhysQubit b, const Params &p) {
    addGate(CZ, a, b, p);
}

void ClusterProg::addSwap(PhysQubit a, PhysQubit b, const Params &p) {
    addGate(SWAP, a, b, p);
}

void ClusterProg::addHadamard(PhysQubit a) {
    addGate(HADAMARD, a, Params(a));
}

// swap qubits [a, m)  [m, b), default params
void ClusterProg::swapRange(PhysQubit a, PhysQubit m, PhysQubit b) {
    Params p(a, b);
    swapRange(a, m, b, p);
}

// swap qubits [a, m)  [m, b)
void ClusterProg::swapRange(PhysQubit a, PhysQubit m, PhysQubit b, const Params &p) {
    int len1 = m - a;
    int len2 = b - m;

    if (len1 == 0 || len2 == 0) return;

    if (len1 >= len2) {
        int r = len1 % len2;
        // swap [a+r...m) with [m,b)
        for (int j = m-len2; j >= a+r; j -= len2) {
            for (int i = 0; i < len2; i++) {
                addSwap(i+j, i+j+len2, p);
            }
        }
        swapRange(a, a+r, a+r+len2, p);
    } else {
        int r = len2 % len1;
        // swap [a, m) with [m, b-r)
        for (int j = a; j < b-r-len1; j += len1) {
            for (int i = 0; i < len1; i++) {
                addSwap(i+j, i+j+len1, p);
            }
        }
        swapRange(b-r-len1, b-r, b, p);
    }
}

// If grid sufficiently large (i.e., greater that 3x3), then
// break the grid up into 4 subgrids before calling the more general doFrame() function
// saves 4 sqrt(n) many qubits
void ClusterProg::doOuterFrame(const Frame& f) {
    int w = f.width();
    int h = f.height();
    int w2 = w/2;
    int w1 = w - w2;
    int h2 = h/2;
    int h1 = h - h2;

    // set up subframes
    Frame sub11(0 ,  0, w1, h1); // Top left grid
    Frame sub12(w1,  0, w2, h1); // Top right grid
    Frame sub21(0 , h1, w1, h2); // Bottom left grid
    Frame sub22(w1, h1, w2, h2); // Bottom right grid

    int p11 = sub11.perimeter();
    int p12 = sub12.perimeter();
    int p21 = sub21.perimeter();
    int p22 = sub22.perimeter();

    // place first subframe starting with east edge, including northeast corner
    sub11.setStart(sub11.x2(), sub11.y1(), 0);

    doFrame(sub11);

    // Measure the qubits that are on the outside perimeter
    Params r11(0, p11);
    for (int i = p11-1; i >= h1 + w1 - 1; i--) {
        GridPoint q11 = sub11[i];
        measurePointAndReset(q11, r11);
        r11.ub--;
    }

    // place second subframe below, starting with northwest corner
    sub21.setStart(sub21.x1(), sub21.y1(), h1 + w1 - 1);

    doFrame(sub21);

    // Measure the qubits that are on the outside perimeter
    Params r21(h1 + w1 - 1, h1 + w1 - 1 + p21);
    for (int i = p21-1; i >= w1+h2-1; i--) {
        GridPoint q21 = sub21[i];
        measurePointAndReset(q21, r21);
        r21.ub--;
    }

    Params r(0, h + 2*w1-2);

    for(int i = 0; i < w1; i++) {
        GridPoint q1 = sub11[h1+i-1];
        GridPoint q2 = sub21[w1-i-1];
        addCSIGN(q1.phys, q2.phys, r);
    }

    // Swap the qubits on the east edge of subframe 21 with the
    // qubits on the south edge of subframe 11
    // and then measure all qubits not on the east edge of both frames
    if(!(h % 2)) { // logic is slightly different for odd and even sidelengths
        for (int i = 0; i < h2-1; i++) {
            GridPoint q1 = sub11[h1+i];
            GridPoint q2 = sub21[w1-1+i];
            addSwap(q1.phys, q2.phys, r);
        }

        addSwap(sub21[0].phys, sub21[w1+h2-2].phys, r);
        measurePointAndReset(GridPoint(sub21[0].x, sub21[0].y, r.ub-1), r);
        r.ub--;

        for(int i = h1+w1-2; i >= h1; i--) {
            GridPoint q = sub11[i];
            measurePointAndReset(GridPoint(q.x, q.y, r.ub-1), r);
            r.ub--;
        }
        for(int i = w1-2; i >= 1; i--) {
            GridPoint q = sub21[i];
            measurePointAndReset(GridPoint(q.x, q.y, r.ub-1), r);
            r.ub--;
        }
    } else {
        for (int i = 0; i < h2; i++) {
            GridPoint q1 = sub11[h1+i];
            GridPoint q2 = sub21[w1-1+i];
            addSwap(q1.phys, q2.phys, r);
        }

        for(int i = h1+w1-2; i >= h1; i--) {
            GridPoint q = sub11[i];
            measurePointAndReset(GridPoint(q.x, q.y, r.ub-1), r);
            r.ub--;
        }
        for(int i = w1-2; i >= 0; i--) {
            GridPoint q = sub21[i];
            measurePointAndReset(GridPoint(q.x, q.y, r.ub-1), r);
            r.ub--;
        }
    }

    // we now have the left frame is just a line
    Frame vert(0,0, w1, h);
    vert.setStart(vert.x2(), vert.y1(), 0);

    // place third subframe (22) starting with southwest corner
    sub22.setStart(sub22.x1(), sub22.y2(), h);

    doFrame(sub22);

    // Measure the qubits that are on the outside perimeter
    Params r22(h, h+p22);
    for (int i = p22-1; i >= h2 + w2 - 1; i--) {
        GridPoint q22 = sub22[i];
        measurePointAndReset(q22, r22);
        r22.ub--;
    }

    r.ub = h + w2 + h2 - 1;
    for(int i = 0; i < w2; i++) {
        GridPoint q1 = vert[h1+i];
        GridPoint q2 = sub22[h2-1-i];
        addCSIGN(q1.phys, q2.phys, r);
    }


    // Swap qubits on north edge of subframe 22
    // with qubits on the east edge of subframe 21
    for (int i = 0; i < h2; i++) {
        GridPoint q1 = vert[h1+i];
        GridPoint q2 = sub22[h2-1+i];
        addSwap(q1.phys, q2.phys, r);
    }

    for(int i = h-1; i >= h1; i--) {
        GridPoint q = vert[i];
        measurePointAndReset(GridPoint(q.x, q.y, r.ub-1), r);
        r.ub--;
    }
    for(int i = h2-2; i >= 0; i--) {
        GridPoint q = sub22[i];
        measurePointAndReset(GridPoint(q.x, q.y, r.ub-1), r);
        r.ub--;
    }

    // place fourth subframe (12) above, starting with southeast corner
    sub12.setStart(sub12.x2(), sub12.y2(), h1+w2);

    doFrame(sub12);

    // Measure the qubits that are on the outside perimeter
    Params r12(h1+w2, h1+w2+p12);
    for (int i = p12-1; i >= w2+h1-1; i--) {
        GridPoint q12 = sub12[i];
        measurePointAndReset(q12, r12);
        r12.ub--;
    }

    // Measure : the last physical qubit and the first physical qubit together
    // after applying CSIGN
    r.ub = 2*h1 + 2*w2 -1;
    for(int i = 0; i < h1-1; i++) {
        addCSIGN(r.lb, r.ub-1, r);
        measurePoint(GridPoint(w1-1, i, r.lb),r);
        r.lb++;
        measurePoint(GridPoint(w1, i, r.ub-1),r);
        r.ub--;
    }
    addCSIGN(r.lb, r.ub-1, r);
    measurePoint(GridPoint(w1-1, h1-1, r.lb),r);
    r.lb++;
    addCSIGN(r.lb, r.ub-1, r);
    measurePoint(GridPoint(w1, h1, r.lb),r);
    r.lb++;
    measurePoint(GridPoint(w1, h1-1, r.ub-1),r);
    r.ub--;

    for(int i = 0; i < w2 -1; i++) {
        addCSIGN(r.lb, r.ub-1, r);
        measurePoint(GridPoint(w1+1+i, h1, r.lb),r);
        r.lb++;
        measurePoint(GridPoint(w1+1+i, h1-1, r.ub-1),r);
        r.ub--;
    }
}

void ClusterProg::doFrame(const Frame& f) {
    int offset = f.offset(); // first physical qubit of the frame
    int w = f.width();
    int h = f.height();
    int p = f.perimeter();
    if (w < 4 && h < 4) { // base case
        //cout << "BASE CASE " << w << "x" << h << " OFFSET " << offset << endl;
        if (w == 3 && h == 3) { // the 3x3 case is the trickiest
            if (f[0].isCornerOf(f)) {
                // do initial hadamard gates
                for (int i = 0; i < 8; ++i) {
                    addHadamard(offset+i);
                }
                GridPoint mid(f.x + 1, f.y + 1, offset); // middle qubit goes in position 0
                for (int i = 1; i < 7; ++i) { // csigns 1-2, ..., 6-7
                    addCSIGN(offset+i, offset+i+1, f.range);
                    if (i % 2 == 1) {
                        addCSIGN(offset+i, offset, f.range); // every other qubit is an edge, connects to middle qubit
                    }
                }
                addCSIGN(offset, offset+7, f.range);
                measurePointAndReset(mid, f.range); // measure the middle uqbit
                addHadamard(mid.phys);
                addCSIGN(offset, offset+1, f.range); // csign 0-1
                addCSIGN(offset+7, offset, f.range); // csign 7-0
            } else {
                // do initial hadamard gates
                for (int i = 0; i < 8; ++i) {
                    addHadamard(offset+i);
                }
                GridPoint mid(f.x + 1, f.y + 1, offset+7); // middle qubit goes in position 7
                for (int i = 0; i < 6; i++) { // csigns 0-1, ..., 5-6
                    addCSIGN(offset+i, offset+i+1, f.range);
                    if (i % 2 == 0) {
                        addCSIGN(offset+i, offset+7, f.range);
                    }
                }
                addCSIGN(offset+6, offset+7, f.range);
                measurePointAndReset(mid, f.range);
                addHadamard(mid.phys);
                addCSIGN(offset+6, offset+7, f.range); // csign 6-7
                addCSIGN(offset+7, offset, f.range); // csign 7-0
            }
        } else {
            // do initial hadamard gates
            for (int i = 0; i < p; i++) {
                addHadamard(offset+i);
            }
            for (int i = 1; i < p; i++) {
                addCSIGN(offset+i-1, offset+i, f.range);
            }
            addCSIGN(offset, offset+p-1, f.range);
            if (w == 3) {
                GridPoint q1(f.x + 1, f.y, -1);
                GridPoint q2(f.x + 1, f.y + 1, -1);
                addCSIGN(offset + f.index(q1), offset + f.index(q2), f.range);
            } else if (h == 3) {
                GridPoint q1(f.x, f.y + 1, -1);
                GridPoint q2(f.x + 1, f.y + 1, -1);
                addCSIGN(offset + f.index(q1), offset + f.index(q2), f.range);
            }
        }
    } else {
        // recursive case
        //cout << "RECURSIVE CASE: " << w << "x" << h << endl;
        if (f.width() >= f.height()) {
            // set up subframes
            int w2 = w/2;
            int w1 = w - w2;
            Frame sub1(f.x1(),      f.y1(), w1, h);
            Frame sub2(f.x1() + w1, f.y1(), w2, h);
            int p1 = sub1.perimeter();
            int p2 = sub2.perimeter();

            // place first subframe at offset, starting with east edge, excluding northeast corner
            sub1.setStart(sub1.x2(), sub1.y1()+1, offset);
            // place second subframe immediately after, start from northwest corner
            sub2.setStart(sub2.x1(), sub2.y1(), offset + p1);

            // do each subframe recursively
            doFrame(sub1);
            doFrame(sub2);

            // cross-frame CSIGNS
            Params r(offset, offset+p1+p2);
            for (int i = 0; i < f.dy-1; i++) {
                GridPoint q1 = sub1[i];
                GridPoint q2 = sub2[-1-i];
                addCSIGN(q1.phys, q2.phys, r);
                measurePointAndReset(q2, r);
                r.ub--;
                measurePointAndReset(q1, r);
                r.lb++;
            }
            int beg = offset;
            int end = offset+p1+p2-1;
            // south pair of cross qubits
            addCSIGN(beg+f.dy-1, end-f.dy+1, r);
            // north pair of cross qubits
            addCSIGN(offset+p1-1, offset+p1, r);

            for (int i = h-3; i >= 0; i--) {
                addSwap(beg+i, beg+p+i, Params(beg, beg+p1+p2));
            }

            // need to rotate the qubits to correct orientation
            // identify the first qubit before rotation
            int index = f.index(sub1[f.dy-1])-h+2;
            // rotate by -index
            swapRange(offset, offset + f.wrap(-index), offset + f.perimeter(), f.range);
        } else {
            int h2 = h/2;
            int h1 = h - h2;
            Frame sub1(f.x1(), f.y1(),      w, h1);
            Frame sub2(f.x1(), f.y1() + h1, w, h2);
            int p1 = sub1.perimeter();
            int p2 = sub2.perimeter();
            // place first subframe at offset, starting on south edge, excluding southeast corner
            sub1.setStart(sub1.x2() - 1, sub1.y2(), offset);
            // place second subframe immediately after, start from northeast corner
            sub2.setStart(sub2.x2(), sub2.y1(), offset + p1);

            // do each subframe recursively
            doFrame(sub1);
            doFrame(sub2);

            // cross-frame CSIGNS
            Params r(offset, offset+p1+p2);
            for (int i = 0; i < f.dx-1; i++) {
                GridPoint q1 = sub1[i];
                GridPoint q2 = sub2[-1-i];
                addCSIGN(q1.phys, q2.phys, r);
                measurePointAndReset(q2, r);
                r.ub--;
                measurePointAndReset(q1, r);
                r.lb++;
            }
            int beg = offset;
            int end = offset+p1+p2-1;
            // west pair of cross qubits
            addCSIGN(beg+f.dx-1, end-f.dx+1, r);
            // east pair of cross qubits
            addCSIGN(offset+p1-1, offset+p1, r);

            for (int i = w-3; i >= 0; i--) {
                addSwap(beg+i, beg+p+i, Params(beg, beg+p1+p2));
            }
            // need to rotate the qubits to correct orientation
            // identify the first qubit before rotation and move w-2 things in front
            int index = f.index(sub1[f.dx-1])-w+2;
            // rotate by -index
            swapRange(offset, offset + f.wrap(-index), offset + f.perimeter(), f.range);
        }
    }
    if (!isFlagSet(PRINT_OUTCOMES)) { // if we're not printing the circuit, run it as we build it to save memory
        run();
    }
}

bool GridPoint::isCornerOf(const Frame &f) const {
    return (x == f.x1() || x == f.x2()) &&
    (y == f.y1() || y == f.y2());
}

int Frame::indexNW(const GridPoint &q) const {
    if (q.y == y1()) { // north edge
        return q.x - x1();
    } else if (q.x == x2()) { //east edge
        return q.y - y1() + dx;
    } else if (q.y == y2()) { // south edge
        return x2() - q.x + dx + dy;
    } else {
        return y2() - q.y + 2*dx + dy;
    }
}

int Frame::index(const GridPoint &q) const {
    return wrap(indexNW(q) - firstIndexNW);
}

GridPoint Frame::operator[](int ind) const {
    int i = wrap(ind + firstIndexNW);
    if (i < dx) {
        return GridPoint(x1() + i, y1(), first.phys + wrap(i - firstIndexNW));
    } else if (i < dx + dy) {
        return GridPoint(x2(), y1() + i - dx, first.phys + wrap(i - firstIndexNW));
    } else if (i < 2*dx + dy) {
        return GridPoint(x2() - i + dx + dy, y2(), first.phys + wrap(i - firstIndexNW));
    } else {
        return GridPoint(x1(), y2() - i + 2*dx + dy, first.phys + wrap(i - firstIndexNW));
    }
}

void ClusterProg::measurePoint(const GridPoint &q, const Params &p, bool reset) {
    Basis b;
    int logi = q.globalPosition(len);
    PhysQubit phys = q.phys;

    if(bases.empty()) {
        b = randombasis();
    } else {
        b = bases.at(logi);
    }

    if(reset) {
        addMeasureReset(phys, p, b);
    } else {
        addMeasurement(phys, p, b);
    }

    qmap[logi] = phys;
    morder.push_back(logi);
}

void ClusterProg::measurePointAndReset(const GridPoint &q, const Params &p) {
    measurePoint(q, p, true);
}

void ClusterProg::generateCluster() {
    generateCluster(len);
}

void ClusterProg::generateCluster(int len) {
    qmap = new int[len*len];

    if (isFlagSet(PRINT_OUTCOMES)) {
        // Prepare outcome vector
        outcomes = new vector<Outcome>();
        outcomes->reserve(len*len); // we should expect len*len measurement outcomes for the grid, clearly
    }

    if(isFlagSet(NAIVE)) {
        Params p(0, len*len);
        // hadamards first
        for (int i = 0; i < len*len; i++) {
            addHadamard(i);
        }

        for(int y = 0; y < len; y++) {
            for(int x = 0; x < len; x++) {
                if(x < len-1) addCSIGN(y*len+x,y*len+x+1,p);
                if(y < len-1) addCSIGN(y*len+x,(y+1)*len+x,p);
                measurePoint(GridPoint(x,y,y*len+x), p);
                p.lb++;
            }
        }
    } else if(isFlagSet(RANDOM_ORDER)) {
        Params p(0, len*len);
        // hadamards first
        for (int i = 0; i < len*len; i++) {
            addHadamard(i);
        }
        // create adjacency list for those edges not yet added
        vector<int> edges[len*len];
        for(int x = 0; x < len; x++) {
            for(int y = 0; y < len; y++) {
                int v = y*len + x;
                if(x > 0) edges[v].push_back(y*len + x-1);
                if(x < len-1) edges[v].push_back(y*len + x+1);
                if(y > 0) edges[v].push_back((y-1)*len + x);
                if(y < len-1) edges[v].push_back((y+1)*len + x);
            }
        }

        int verts[len*len];
        for(int i = 0; i < len*len; i++) {
            verts[i] = i;
        }
        int vr = len*len;
        for(int i = 0; i < len*len; i++) {
            int r = rand() % vr;
            int v = verts[r];

            // apply all CSIGNS that haven't yet been applied to v
            for(int u : edges[v]) {
                addCSIGN(v, u, p);
                // remove v from the adjacency list of u
                int idx = 0;
                for(int z = 0; z < (int)(edges[u].size()); z++) {
                    if(edges[u][z] == v) {
                        idx = z;
                        break;
                    }
                }
                edges[u][idx] = edges[u][edges[u].size()-1];
                edges[u].pop_back();
            }
            edges[v].clear();
            measurePoint(GridPoint(v % len, v/len, v), p);

            vr--;
            verts[r] = verts[vr];
        }

    } else if(isFlagSet(ROW_BY_ROW)) {
        /* INTERLEAVED MEASUREMENTS */
        Params p(0, 2*len);
        for (int i = 0; i < 2*len; i++) {
            addHadamard(i);
        }
        for(int y = 0; y < len; y++) {
            if(y % 2) {
                for(int x = len-1; x >= 0; x--) {
                    Params p(x, len+x+1);

                    if(x > 0) addCSIGN(x-1+len, x+len, p);
                    if(y < len-1) addCSIGN(x+len, x, p);
                    measurePointAndReset(GridPoint(x,y,len+x), p);
                    addHadamard(len+x);
                }
            } else {
                for(int x = 0; x < len; x++) {
                    Params p(x, len+x+1);

                    if(x < len-1) addCSIGN(x, x+1, p);
                    if(y < len-1) addCSIGN(x, len+x, p);
                    measurePointAndReset(GridPoint(x,y,x), p);
                    addHadamard(x);
                }
            }
        }
    } else {
        n = 6*len; // estimate 6n qubits

        Frame f(0, 0, len, len);

        if (len < 4) {
            // set first qubit of frame to be something convenient if vertically split
            f.setStart(f.x1() + f.width() - f.width()/2 + f.height() - 3, f.y2(), 0);

            doFrame(f);

            // measure remaining qubits
            int n = f.perimeter();
            for (int i = n-1; i >= 0; i--) {
                measurePoint(f[i], Params(0,i+1));
            }
        } else { // assumes we're not in a base case
            doOuterFrame(f);
        }

    }

    if(isFlagSet(VERIFY)) {
        gates.clear();  // reset gates for naive implementation
        vector<int> old_morder; // but keep the order of measurements used
        for(int i : morder) { old_morder.push_back(i); }
        morder.clear();
        // naive implementation of grid
        for (int i = 0; i < len*len; i++) { addHadamard(i); }
        Params p(0, len*len);
        for(int y = 0; y < len; y++) {
            for(int x = 0; x < len; x++) {
                if(x < len-1) addCSIGN(y*len+x,y*len+x+1,p);
                if(y < len-1) addCSIGN(y*len+x,(y+1)*len+x,p);
            }
        }
        for(int i : old_morder) {
            int x = i % len;
            int y = i / len;
            measurePoint(GridPoint(x,y,i), p);
        }
    }
}

void ClusterProg::printOutcomes() const {
    for(unsigned int i = 0; i < morder.size(); i++) {
        /* CHP-type output */
        // cout << "Outcome of measuring qubit " << morder[i] << ": " << (*outcomes)[i] << endl;
        /* Stim-type output */
        cout << outcomeToBinary((*outcomes)[i]);
    }
}
