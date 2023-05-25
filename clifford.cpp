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

#include <iostream>
#include <string>
#include "clifford.hpp"
#include <unordered_map>
#include <chrono>
#include <math.h>

using namespace std;

bool operator==(const Gate& g, const Opcode op) {
    return g.op == op;
}

ostream& operator<<(ostream& out, Outcome o) {
    switch (o) {
        case ZERO:
            out << 0;
            break;
        case ONE:
            out << 1;
            break;
        case ZERO_RAND:
            out << "0 (random)";
            break;
        case ONE_RAND:
            out << "1 (random)";
            break;
        case RANDOM:
            out << "unspecified (random)";
            break;
    }

    return out;
}
int outcomeToBinary(Outcome o) {
    switch (o) {
        case ZERO:
        case ZERO_RAND:
            return 0;
        case ONE:
        case ONE_RAND:
            return 1;
        default:
            return -1;
    }
}

void Params::setParams(QUBIT lb) {
    this->lb = lb;
    this->ub = lb+1;
}
void Params::setParams(QUBIT lb, QUBIT ub) {
    this->lb = lb;
    this->ub = ub;
}
void Params::setParams(const Params& p) {
    this->lb = p.lb;
    this->ub = p.ub;
}

ostream& operator<<(ostream& out, const Params& p) {
    out << "[" << p.lb << ":" << p.ub << "]";
    return out;
}

// assumes gate params and input params are both correct
// chooses smallest param interval
void Gate::restrictParams(const Params& p) {
    this->p.lb = max(this->p.lb, p.lb);
    this->p.ub = min(this->p.ub, p.ub);
}
void Gate::print() const {
    const GateDef & gd = GATE_DEFINITIONS[op];
    cout << gd.name << " ";
    if (op == MEASURE_RESET) {
        cout << a << " " << (Outcome)b;
    } else if (gd.args == 1) {
        cout << a;
    } else if (gd.args == 2) {
        cout << a << " " << b;
    }
    cout << " " << p.lb << " " << p.ub << endl;
}


// ************ ROW DEFINITIONS **************
WORD & Row::blockX(QUBIT i) {
    return x[i >> SHIFT].u64[(i >> 6) % INTS_PER_BLOCK];
}
WORD & Row::blockZ(QUBIT i) {
    return z[i >> SHIFT].u64[(i >> 6) % INTS_PER_BLOCK];
}
const WORD & Row::blockX(QUBIT i) const {
    return x[i >> SHIFT].u64[(i >> 6) % INTS_PER_BLOCK];
}
const WORD & Row::blockZ(QUBIT i) const {
    return z[i >> SHIFT].u64[(i >> 6) % INTS_PER_BLOCK];
}
WORD Row::bitX(QUBIT i) const {
    return blockX(i) & BIT(i);
}
WORD Row::bitZ(QUBIT i) const {
    return blockZ(i) & BIT(i);
}

SIGNBITS Row::pauli(QUBIT i) const {
    if (bitX(i)) {
        if (bitZ(i)) {
            return 'Y';
        } else {
            return 'X';
        }
    } else {
        if (bitZ(i)) {
            return 'Z';
        } else {
            return 'I';
        }
    }
}

// computes the overall phase of the row
SIGNBITS Row::phase() const { return phase(Params(0,n)); }
SIGNBITS Row::phase(const Params &p) const {
    TBLX c1{};
    TBLX c2{};
    for (QUBIT js = p.lbs(); js < p.ubs(); js++) {
        simd tmp = x[js] & z[js];
        c2 ^= c1 & tmp;
        c1 ^= tmp;
    }
    return (r + 3*(2*c2.popcount() + c1.popcount())) % 4;
}

ostream& operator<<(ostream& out, const Row& r) {

    SIGNBITS ph = r.phase();
    if (ph == 2) {
        out << "-";
    } else if (ph == 0) {
        out << " ";
    } else {
        out << "!";
    }
    for (QUBIT i = 0; i < r.n; i++) {
        out << r.pauli(i);
    }
    //assert(ph == 0 || ph == 2);
    return out;
}
int num_row_mults = 0;

void QState::rowMultBatch(Row *s, Row* t, QUBIT piv, QUBIT* idxs, int size, const Params &p) {
    num_row_mults+=size;

    TBLX accs[size];
    for(QUBIT i=0; i<size; i++) accs[i].zero();
    TBLX* piv_x = s[piv].x;
    TBLX* piv_z = s[piv].z;
    char piv_r = s[piv].r;
    QUBIT start = p.lbs();
    QUBIT stop = p.ubs();
    for (QUBIT js = start; js < stop; js++) {
        for(QUBIT i=0; i<size; i++) {
            accs[i] ^= piv_z[js] & t[idxs[i]].x[js];

            t[idxs[i]].x[js] ^= piv_x[js];
            t[idxs[i]].z[js] ^= piv_z[js];
        }
    }
    for(QUBIT i=0; i<size; i++) {
        t[idxs[i]].r = (piv_r + t[idxs[i]].r + 2*accs[i].parity()) % 4;
    }
}


void Row::leftMultBy(const Row& a, const Params& p, bool withPhase) {
    num_row_mults++;

    TBLX * ax = a.x;
    TBLX * az = a.z;
    QUBIT start = p.lbs();
    QUBIT stop = p.ubs();
    if (withPhase) {
        TBLX acc{};
        acc.zero(); // added for ARM
        for (QUBIT js = start; js < stop; js++) {
            acc ^= (az[js] & x[js]);
            x[js] ^= ax[js];
            z[js] ^= az[js];
        }
        r = (a.r + r + 2*acc.parity()) % 4;
    } else {
        for (QUBIT js = start; js < stop; js++) {
            x[js] ^= ax[js];
            z[js] ^= az[js];
        }
    }
}


/*
Initialize state q to have n qubits
*/
QState::QState(QUBIT n, bool destabilizers) : n(n), blocks((n & ((1ull << SHIFT) - 1)) == 0 ? (n >> SHIFT) : (n >> SHIFT) + 1) {
    try {

        int cache_fix_blocks = 1 - blocks % 2; // don't want the blocks to align in large array
        const QUBIT ROW_LENGTH = blocks + cache_fix_blocks;

        stab = new Row[n]; // pointers to the stabilizers
        rows.push_back(stab); // simulation always includes stabilizers

        if(destabilizers) {
            dest = new Row[n]; // pointers to the destabilizers
            rows.push_back(dest);

            master = new TBLX[4*n*ROW_LENGTH];
        } else {
            dest = NULL;
            master = new TBLX[2*n*ROW_LENGTH];
        }


        for (QUBIT i = 0; i < n; i++) {
            stab[i].x = master + (2*i + 0)*ROW_LENGTH;
            stab[i].z = master + (2*i + 1)*ROW_LENGTH;

            if(destabilizers) {
                dest[i].x = master + (2*(i+n) + 0)*ROW_LENGTH;
                dest[i].z = master + (2*(i+n) + 1)*ROW_LENGTH;
            }

            for (QUBIT j = 0; j < ROW_LENGTH; ++j) {
                stab[i].x[j].zero();
                stab[i].z[j].zero();
                if(destabilizers) {
                    dest[i].x[j].zero();
                    dest[i].z[j].zero();
                }
            }
        }

        // initialize to all |0> state
        for (QUBIT i = 0; i < n; i++) {
            stab[i].n = n;
            stab[i].r = 0;
            stab[i].blockZ(i) = BIT(i);
            if(destabilizers) {
                dest[i].n = n;
                dest[i].r = 0;
                dest[i].blockX(i) = BIT(i);
            }
        }
    } catch (const std::bad_alloc &e) {
        cout << "Allocation failed: " << e.what() << endl;
        exit(0);
    }
}


// Print the destabilizers and stabilizers for state q
void QState::print() const {
    if(dest != NULL) { // print destabilizers
        for (QUBIT i = 0; i < n; i++) {
            cout << dest[i] << endl;
        }
        cout << " ";
        for (QUBIT i = 0; i < n; i++) {
            cout << "-";
        }
        cout << endl;
    }
    // always print stabilizers
    for (QUBIT i = 0; i < n; i++) {
        cout << stab[i] << endl;
    }
}


/* PAULI_X */
void QState::PauliX(QUBIT a) { PauliX(a, Params(0,n)); }
void QState::PauliX(QUBIT a, const Params &p) {
    for(Row* row : rows) {
        for (QUBIT i = p.lb; i < p.ub; i++) {
            if (row[i].bitZ(a))
                row[i].r = (row[i].r + 2) % 4;
        }
    }
}

/* PAULI_Y */
void QState::PauliY(QUBIT a) { PauliY(a, Params(0,n)); }
void QState::PauliY(QUBIT a, const Params &p) {
    for(Row* row : rows) {
        for (QUBIT i = p.lb; i < p.ub; i++) {
            if (row[i].bitZ(a) ^ row[i].bitX(a))
                row[i].r = (row[i].r + 2) % 4;
        }
    }
}

/* PAULI_Z */
void QState::PauliZ(QUBIT a) { PauliZ(a, Params(0,n)); }
void QState::PauliZ(QUBIT a, const Params &p) {
    for(Row* row : rows) {
        for (QUBIT i = p.lb; i < p.ub; i++) {
            if (row[i].bitX(a))
                row[i].r = (row[i].r + 2) % 4;
        }
    }
}

/* HADAMARD_XY */
void QState::hadamard_xy(QUBIT a) { hadamard_xy(a,Params(0,n)); }
void QState::hadamard_xy(QUBIT a, const Params &p) {
    for(Row* row : rows) {
        for (QUBIT i = p.lb; i < p.ub; i++) {
            if (row[i].bitX(a))
            row[i].r = (row[i].r + 1) % 4;
            if (row[i].bitZ(a))
            row[i].r = (row[i].r + 2) % 4;

            row[i].blockZ(a) ^= row[i].bitX(a);
        }
    }
}

/* HADAMARD */
void QState::hadamard(QUBIT a) { hadamard(a,Params(0,n)); }
void QState::hadamard(QUBIT a, const Params &p) {
    for(Row* row : rows) {
        for (QUBIT i = p.lb; i < p.ub; i++) {
            WORD x_part = row[i].bitX(a);
            WORD z_part = row[i].bitZ(a);

            WORD flip = x_part ^ z_part;

            row[i].blockX(a) ^= flip;
            row[i].blockZ(a) ^= flip;

            if (x_part & z_part)
                row[i].r = (row[i].r + 2) % 4;
        }
    }
}

/* HADAMARD_YZ */
void QState::hadamard_yz(QUBIT a) { hadamard_yz(a,Params(0,n)); }
void QState::hadamard_yz(QUBIT a, const Params &p) {
    for(Row* row : rows) {
        for (QUBIT i = p.lb; i < p.ub; i++) {
            if (row[i].bitX(a))
            row[i].r = (row[i].r + 2) % 4;
            if (row[i].bitZ(a))
            row[i].r = (row[i].r + 1) % 4;

            row[i].blockX(a) ^= row[i].bitZ(a);
        }
    }
}

/* SQRT_X */
void QState::rotX90(QUBIT a) { rotX90(a, Params(0,n)); }
void QState::rotX90(QUBIT a, const Params &p) {
    for(Row* row : rows) {
        for (QUBIT i = p.lb; i < p.ub; i++) {
            row[i].blockX(a) ^= row[i].bitZ(a);

            if (row[i].bitZ(a))
            row[i].r = (row[i].r + 3) % 4;
        }
    }
}

/* SQRT_X_DAG */
void QState::rotX270(QUBIT a) { rotX270(a, Params(0,n)); }
void QState::rotX270(QUBIT a, const Params &p) {
    for(Row* row : rows) {
        for (QUBIT i = p.lb; i < p.ub; i++) {
            row[i].blockX(a) ^= row[i].bitZ(a);

            if (row[i].bitZ(a))
            row[i].r = (row[i].r + 1) % 4;
        }
    }
}

/* SQRT_Y */
void QState::rotY90(QUBIT a) { rotY90(a, Params(0,n)); }
void QState::rotY90(QUBIT a, const Params &p) {
    for(Row* row : rows) {
        for (QUBIT i = p.lb; i < p.ub; i++) {
            row[i].blockZ(a) ^= row[i].bitX(a);
            row[i].blockX(a) ^= row[i].bitZ(a);
            row[i].blockZ(a) ^= row[i].bitX(a);

            if (row[i].bitZ(a) && !row[i].bitX(a))
                row[i].r = (row[i].r + 2) % 4;
        }
    }
}

/* SQRT_Y_DAG */
void QState::rotY270(QUBIT a) { rotY270(a, Params(0,n)); }
void QState::rotY270(QUBIT a, const Params &p) {
    for(Row* row : rows) {
        for (QUBIT i = p.lb; i < p.ub; i++) {
            row[i].blockZ(a) ^= row[i].bitX(a);
            row[i].blockX(a) ^= row[i].bitZ(a);
            row[i].blockZ(a) ^= row[i].bitX(a);

            if (row[i].bitX(a) && !row[i].bitZ(a))
            row[i].r = (row[i].r + 2) % 4;
        }
    }
}

/* S - Phase */
void QState::phase(QUBIT a) { rotZ90(a, Params(0,n)); }
void QState::phase(QUBIT a, const Params &p) { rotZ90(a, p); }
void QState::rotZ90(QUBIT a) { rotZ90(a, Params(0,n)); }
void QState::rotZ90(QUBIT a, const Params &p) {
    for(Row* row : rows) {
        for (QUBIT i = p.lb; i < p.ub; i++) {
            row[i].blockZ(a) ^= row[i].bitX(a);

            if (row[i].bitX(a))
            row[i].r = (row[i].r + 1) % 4;
        }
    }
}

/* S_DAG */
void QState::rotZ270(QUBIT a) { rotZ270(a, Params(0,n)); }
void QState::rotZ270(QUBIT a, const Params &p) {
    for(Row* row : rows) {
        for (QUBIT i = p.lb; i < p.ub; i++) {
            row[i].blockZ(a) ^= row[i].bitX(a);

            if (row[i].bitX(a))
            row[i].r = (row[i].r + 3) % 4;
        }
    }
}

/* SWAP */
void QState::swapGate(QUBIT a, QUBIT b) { swapGate(a,b,Params(0,n)); }
void QState::swapGate(QUBIT a, QUBIT b, const Params &p) {
    QUBIT maskA = a & MASK;
    QUBIT maskB = b & MASK;
    if (maskA < maskB) {
        for(Row* row : rows) {
            for (QUBIT i = p.lb; i < p.ub; i++) {
                WORD x_xor = (row[i].bitX(a) << (maskB-maskA)) ^ (row[i].bitX(b));
                row[i].blockX(a) ^= x_xor >> (maskB-maskA);
                row[i].blockX(b) ^= x_xor;
                WORD z_xor = (row[i].bitZ(a) << (maskB-maskA)) ^ (row[i].bitZ(b));
                row[i].blockZ(a) ^= z_xor >> (maskB-maskA);
                row[i].blockZ(b) ^= z_xor;
            }
            swap(row[a], row[b]); // swap rows
        }
    } else {
        for(Row* row : rows) {
            for (QUBIT i = p.lb; i < p.ub; i++) {
                WORD x_xor = (row[i].bitX(a) >> (maskA-maskB)) ^ (row[i].bitX(b));
                row[i].blockX(a) ^= x_xor << (maskA-maskB);
                row[i].blockX(b) ^= x_xor;
                WORD z_xor = (row[i].bitZ(a) >> (maskA-maskB)) ^ (row[i].bitZ(b));
                row[i].blockZ(a) ^= z_xor << (maskA-maskB);
                row[i].blockZ(b) ^= z_xor;
            }
            swap(row[a], row[b]); // swap rows
        }
    }
}

/* iSWAP = SWAP * CZ * (S x S) */
void QState::iSWAP(QUBIT a, QUBIT b) { iSWAP(a,b,Params(0,n)); }
void QState::iSWAP(QUBIT a, QUBIT b, const Params &p) {
    QUBIT maskA = a & MASK;
    QUBIT maskB = b & MASK;
    if (maskA > maskB) {
        swap(a, b);
        maskA = a & MASK;
        maskB = b & MASK;
    }
    for(Row* row : rows) {
        for (QUBIT i = p.lb; i < p.ub; i++) {
            WORD z_xor = (row[i].bitZ(a) << (maskB-maskA)) ^ (row[i].bitZ(b));
            row[i].blockZ(a) ^= z_xor >> (maskB-maskA);
            row[i].blockZ(b) ^= z_xor;

            WORD x_a = (row[i].bitX(a) << (maskB-maskA));
            WORD x_b = (row[i].bitX(b));
            WORD x_xor = x_a ^ x_b;

            if(x_xor) {
                row[i].blockX(a) ^= BIT(a);
                row[i].blockX(b) ^= BIT(b);
                row[i].blockZ(a) ^= BIT(a);
                row[i].blockZ(b) ^= BIT(b);

                if(!(x_a & x_b))
                    row[i].r = (row[i].r + 1) % 4;
            }
        }
    }
}

/* iSWAP_DAG = SWAP * CZ * (S_DAG x S_DAG) */
void QState::iSWAP_dag(QUBIT a, QUBIT b) { iSWAP_dag(a,b,Params(0,n)); }
void QState::iSWAP_dag(QUBIT a, QUBIT b, const Params &p) {
    QUBIT maskA = a & MASK;
    QUBIT maskB = b & MASK;
    if (maskA > maskB) {
        swap(a, b);
        maskA = a & MASK;
        maskB = b & MASK;
    }
    for(Row* row : rows) {
        for (QUBIT i = p.lb; i < p.ub; i++) {
            WORD z_xor = (row[i].bitZ(a) << (maskB-maskA)) ^ (row[i].bitZ(b));
            row[i].blockZ(a) ^= z_xor >> (maskB-maskA);
            row[i].blockZ(b) ^= z_xor;

            WORD x_a = (row[i].bitX(a) << (maskB-maskA));
            WORD x_b = (row[i].bitX(b));
            WORD x_xor = x_a ^ x_b;

            if(x_xor) {
                row[i].blockX(a) ^= BIT(a);
                row[i].blockX(b) ^= BIT(b);
                row[i].blockZ(a) ^= BIT(a);
                row[i].blockZ(b) ^= BIT(b);

                if(!(x_a & x_b))
                    row[i].r = (row[i].r + 3) % 4;
            }
        }
    }
}
/* XCX - X-Controlled-X */
void QState::xcx(QUBIT a, QUBIT b) { xcx(a,b,Params(0,n)); }
void QState::xcx(QUBIT a, QUBIT b, const Params &p) {
    for(Row* row : rows) {
        for (QUBIT i = p.lb; i < p.ub; i++) {
            if (row[i].bitZ(a)) row[i].blockX(b) ^= BIT(b);
            if (row[i].bitZ(b)) row[i].blockX(a) ^= BIT(a);

            if (row[i].bitZ(a) && row[i].bitZ(b))
                row[i].r = (row[i].r + 2) % 4;
        }
    }
}

/* XCY - X-Controlled-Y */
void QState::xcy(QUBIT a, QUBIT b) { xcy(a,b,Params(0,n)); }
void QState::xcy(QUBIT a, QUBIT b, const Params &p) {
    for(Row* row : rows) {
        for (QUBIT i = p.lb; i < p.ub; i++) {
            if(row[i].bitX(b) ^ row[i].bitZ(b)) row[i].blockX(a) ^= BIT(a);
            if(row[i].bitZ(a)) {
                row[i].blockZ(b) ^= BIT(b);
                row[i].blockX(b) ^= BIT(b);

                if(row[i].bitZ(b)) {
                    row[i].r = (row[i].r + 1) % 4;
                } else {
                    row[i].r = (row[i].r + 3) % 4;
                }
            }
        }
    }
}

/* XCZ - X-Controlled-Z */
void QState::xcz(QUBIT a, QUBIT b) { cx(b,a,Params(0,n)); }
void QState::xcz(QUBIT a, QUBIT b, const Params &p) { cx(b,a,p); }

/* YCX - Y-Controlled-X */
void QState::ycx(QUBIT a, QUBIT b) { xcy(b,a,Params(0,n)); }
void QState::ycx(QUBIT a, QUBIT b, const Params &p) { xcy(b,a,p); }

/* YCY - Y-Controlled-Y */
void QState::ycy(QUBIT a, QUBIT b) { ycy(a,b,Params(0,n)); }
void QState::ycy(QUBIT a, QUBIT b, const Params &p) {
    for(Row* row : rows) {
        for (QUBIT i = p.lb; i < p.ub; i++) {
            WORD a_z = row[i].bitZ(a);
            WORD b_x = row[i].bitX(b);
            WORD b_xor = row[i].bitZ(b) ^ b_x;
            WORD a_xor = row[i].bitX(a) ^ a_z;

            if(b_xor) {
                row[i].blockX(a) ^= BIT(a);
                row[i].blockZ(a) ^= BIT(a);
                row[i].r += 1;
                if(a_z) row[i].r += 2;
            }
            if(a_xor) {
                row[i].blockX(b) ^= BIT(b);
                row[i].blockZ(b) ^= BIT(b);
                row[i].r += 1;
                if(b_x) row[i].r += 2;
            }
            row[i].r %= 4;
        }
    }
}

/* YCZ - Y-Controlled-Z */
void QState::ycz(QUBIT a, QUBIT b) { cy(b,a,Params(0,n)); }
void QState::ycz(QUBIT a, QUBIT b, const Params &p) { cy(b,a,p); }

/* CX - CNOT */
void QState::cx(QUBIT a, QUBIT b) { cnot(a,b,Params(0,n)); }
void QState::cx(QUBIT a, QUBIT b, const Params &p) { cnot(a,b,p); }
void QState::cnot(QUBIT a, QUBIT b) { cnot(a,b,Params(0,n)); }
void QState::cnot(QUBIT a, QUBIT b, const Params &p) {
    for(Row* row : rows) {
        for (QUBIT i = p.lb; i < p.ub; i++) {
            if (row[i].bitX(a)) row[i].blockX(b) ^= BIT(b);
            if (row[i].bitZ(b)) row[i].blockZ(a) ^= BIT(a);
        }
    }
}

/* CY - Controlled-Y */
void QState::cy(QUBIT a, QUBIT b) { cy(a,b,Params(0,n)); }
void QState::cy(QUBIT a, QUBIT b, const Params &p) {
    for(Row* row : rows) {
        for (QUBIT i = p.lb; i < p.ub; i++) {
            if(row[i].bitZ(b) ^ row[i].bitX(b)) row[i].blockZ(a) ^= BIT(a);
            if(row[i].bitX(a)) {
                row[i].blockX(b) ^= BIT(b);
                row[i].blockZ(b) ^= BIT(b);

                if(row[i].bitX(b)) {
                    row[i].r = (row[i].r + 1) % 4;
                } else {
                    row[i].r = (row[i].r + 3) % 4;
                }
            }
        }
    }
}

/* CZ - CSIGN */
void QState::cz(QUBIT a, QUBIT b) { csign(a,b,Params(0,n)); }
void QState::cz(QUBIT a, QUBIT b, const Params &p) { csign(a,b,p); }
void QState::csign(QUBIT a, QUBIT b) { csign(a,b,Params(0,n)); }
void QState::csign(QUBIT a, QUBIT b, const Params &p) {
    for(Row* row : rows) {
        for (QUBIT i = p.lb; i < p.ub; i++) {
            if (row[i].bitX(a)) row[i].blockZ(b) ^= BIT(b);
            if (row[i].bitX(b)) row[i].blockZ(a) ^= BIT(a);

            if (row[i].bitX(a) && row[i].bitX(b))
                row[i].r = (row[i].r + 2) % 4;
        }
    }
}


/*
MEASURE (in the Z-basis)

Returns:
ZERO if outcome would always be 0
ONE if outcome would always be 1
ZERO_RAND if outcome was random and 0 was chosen
ONE_RAND if outcome was random and 1 was chosen

QUBIT a - qubit to be measured
Outcome o - preferred outcome of the measurement
        i.e., if the measurement is random, then outcome "o" is returned;
        if the measurement is deterministic, then "o" is ignored
*/
Outcome QState::measure(QUBIT a, Outcome o, const Params &p) {
    Outcome outcome; // measurement outcome

    char random = 0; // 1: measurement outcome is random 0: measurement outcome is deterministic
    QUBIT piv;       // pivot row in stabilizer/destabilizer
    const int BATCH_NUM = 4;
    QUBIT idxs[BATCH_NUM];
    int unprocessed = 0;

    // loop over stabilizers to find X-Pauli
    // if found, measurement is random
    for (piv = p.lb; piv < p.ub; piv++) {
        if (stab[piv].bitX(a)) {
            random = 1;
            break;
        }
    }


    if (random) {
        // Multiply stabilizers so only generator "piv" has Pauli-X term
        for (QUBIT i = piv + 1; i < p.ub; i++) {
            if (i != piv && stab[i].bitX(a)) {
                idxs[unprocessed] = i;

                if(unprocessed < BATCH_NUM-1) {
                    unprocessed++;
                } else {
                    unprocessed=0;
                    rowMultBatch(stab, stab, piv, idxs, BATCH_NUM, p);
                }
            }
        }
        if(unprocessed) {
            rowMultBatch(stab, stab, piv, idxs, unprocessed, p);
        }

        // if there are destabilizers, we must simplify them
        if(dest != NULL) {
            unprocessed=0;
            // Multiply destabilizers so none has Pauli-X term
            for (QUBIT i = p.lb; i < p.ub; i++) {
                if (i != piv && dest[i].bitX(a)) {
                    idxs[unprocessed] = i;

                    if(unprocessed < BATCH_NUM-1) {
                        unprocessed++;
                    } else {
                        unprocessed=0;
                        rowMultBatch(stab, dest, piv, idxs, BATCH_NUM, p);
                    }

                }
            }
            if(unprocessed) {
                rowMultBatch(stab, dest, piv, idxs, unprocessed, p);
            }
        }

        // replace the stabilizer with I...IZI...I
        rowSetZ(stab[piv], a, p);

        if (o != RANDOM) {
            stab[piv].r = 2 * o; // fixed outcome for testing
        } else {
            stab[piv].r = 2 * (rand() % 2); // if not fixed, choose random outcome
        }
        if (stab[piv].r)
            outcome = ONE_RAND;
        else
            outcome = ZERO_RAND;
    } else { /* Deterministic measurement */
        if(dest != NULL) {
            // Search destabilizers for X-Pauli (must exist since measurement is not random)
            for (piv = p.lb; piv < p.ub; piv++) {
                if (dest[piv].bitX(a))
                    break;
            }

            // Multiply destabilizers so no other has Pauli-X term
            // Multiply stabilizer so that tableau is still symplectic
            for (QUBIT i = piv + 1; i < p.ub; i++) {
                if (dest[i].bitX(a)) {
                    rowMult(dest[piv], dest[i], p, false);
                    rowMult(stab[i], stab[piv], p, true);
                }
            }
        } else {
            piv = gaussianEliminate(a, p);
        }

        outcome = (Outcome)(stab[piv].phase(p) / 2);
    }

    /* Restore the tableau to direct-product structure */
    // Swap the a'th and piv'th generators
    if (a != piv) {
        rowSwap(stab[a], stab[piv], p);
        if(dest != NULL) rowSwap(dest[a], dest[piv], p);
    }

    // Remove all other Pauli Z
    for (QUBIT i = p.lb; i < p.ub; i++) {
        if(dest != NULL) dest[i].blockZ(a) &= ~(BIT(a));

        if (i != a && stab[i].bitZ(a)) {
            stab[i].blockZ(a) ^= BIT(a);
            stab[i].r = (stab[i].r + stab[a].r) % 4;
        }
    }

    // The destabilizers can be simplified
    if(dest != NULL) rowSetX(dest[a], a, p);

    return outcome;
}


// returns stabilizer row that will eventually have IIZII with the Z in position 'a'
// or -1 if no such stabilizer row exists
void QState::gaussianEliminate() { gaussianEliminate(0, Params(0,n)); }
QUBIT QState::gaussianEliminate(QUBIT a, const Params &p) {
    int rank = 0; // rank of the stabilizer Pauli-X matrix
    QUBIT piv = 0; // pivot row for Gaussian elimination
    QUBIT z_row = -1; // this is the row ID that will eventually have IIZII with the Z in position 'a'

    //reduced row echelon form
    for (QUBIT j = p.lb; j < p.ub; j++) {
        piv = p.lb+rank;
        // find pivot row with Pauli X in column j, starting with row j
        for(QUBIT i = p.lb+rank; i < p.ub; i++) {
            if(stab[i].bitX(j)) break;
            piv++;
        }
        //if no Pauli-X term then nothing to do for this column
        if(piv == p.ub) continue;
        // otherwise, reduce by the pivot row
        for(QUBIT i = p.lb; i < p.ub; i++) {
            if(i!=piv && stab[i].bitX(j)) {
                rowMult(stab[piv], stab[i], p, true);
                if(dest != NULL) rowMult(dest[i], dest[piv], p, false);
            }
        }
        // Swap the rank'th and piv'th generators
        if (p.lb+rank != piv) {
            rowSwap(stab[p.lb+rank], stab[piv], p);
            if(dest != NULL) rowSwap(dest[p.lb+rank], dest[piv], p);
        }
        rank++;
    }

    // We now perform Gaussian elimination on the Z-parts of the tableau for the zero rows above
    int rx = 0; // we will keep track of the rank of the Pauli-X matrix as we include more rows
    int rz = 0; // rank z
    for(QUBIT j = p.lb; j < p.ub; j++) {
        if(stab[p.lb+rx].bitX(j)) { // this was a pivot row
            rx++;
            continue;
        }
        // otherwise, search for Pauli Z (which must exist)
        for(QUBIT i = p.lb+rank+rz; i < p.ub; i++) {
            if(stab[i].bitZ(j)) {
                piv = i;
                if(j==a) z_row = p.lb+rank+rz; // this row will end up being IIIZIII
                break;
            }
        }
        // eliminate other Z
        for(QUBIT i = p.lb+rank; i < p.ub; i++) {
            if(i!=piv && stab[i].bitZ(j)) {
                rowMult(stab[piv], stab[i], p, true);
                if(dest != NULL) rowMult(dest[i], dest[piv], p, false);
            }
        }
        // Swap the (rank+rz)'th and piv'th generators
        if (p.lb+rank+rz != piv) {
            rowSwap(stab[p.lb+rank+rz], stab[piv], p);
            if(dest != NULL) rowSwap(dest[p.lb+rank+rz], dest[piv], p);
        }
        rz++;
    }
    return z_row;
}

int QState::postselect(QUBIT a, Outcome o, const Params &p) {
    Outcome outcome = measure(a, o, p);

    if (outcome == ZERO_RAND || outcome == ONE_RAND) { // measurement was random
        return 1; // measurement probability = 1/2
    } else if (outcome == o) { // measurement outcome was deterministic and correct
        return 2; // measurement probability = 2/2
    } else { // measurement outcome was deterministic and incorrect
        stab[a].r = o == ZERO ? 0 : 2; // fix sign bit
        return 0; //measurement probability = 0/2
    }
}

void QState::zeroOutPhase(QUBIT a) { stab[a].r = 0; }

void QState::reset(QUBIT a) { reset(a, Params(0,n)); }
void QState::reset(QUBIT a, const Params &p) {
    measure(a, RANDOM, p);
    zeroOutPhase(a);
}



// Copy row a onto row b
void QState::rowCopy(const Row & a, Row & b, const Params &cols) {
    for (QUBIT js = cols.lbs(); js < cols.ubs(); js++) {
        b.x[js] = a.x[js];
        b.z[js] = a.z[js];
    }
    b.r = a.r;
}

void QState::rowSwap(Row & a, Row & b, const Params &cols) {
    for (QUBIT js = cols.lbs(); js < cols.ubs(); js++) {
        swap(a.x[js], b.x[js]);
        swap(a.z[js], b.z[js]);
    }
    swap(a.r, b.r);
    swap(a.n, b.n);
}

// Set a row to zero, i.e., +III...III.
void QState::rowZero(Row & a, const Params &cols) {
    for (QUBIT js = cols.lbs(); js < cols.ubs(); js++) {
        a.x[js].zero();
        a.z[js].zero();
    }
    a.r = 0;
}

void QState::rowSetX(Row & a, QUBIT i, const Params &cols) {
    rowZero(a,cols);
    a.blockX(i) = BIT(i);
}
void QState::rowSetZ(Row & a, QUBIT i, const Params &cols) {
    rowZero(a,cols);
    a.blockZ(i) = BIT(i);
}

// Left-multiply row a into row b
void QState::rowMult(const Row & a, Row & b, const Params &cols, bool withPhase) {
    b.leftMultBy(a, cols, withPhase);
}

/*
Prints the quantum program in a format understood by CHP
*/
void QProg::printCHP() const {
    for (const Gate &g : gates) {
        const GateDef & gd = GATE_DEFINITIONS[g.op];
        if (gd.op == POSTSELECT) {
            cout << gd.name << " " << g.a << " " << g.o << endl;
        } else if (gd.args == 1) {
            cout << gd.name << " " << g.a << endl;
        } else if (gd.args == 2) {
            cout << gd.name << " " << g.a << " " << g.b << endl;
        } else {
            throw domain_error("attempting to print unimplemented gate");
        }
    }
}

void QProg::print() const {
    cout << "Circuit on " << n << " qubits" << endl;
    for (const Gate &g : gates) {
        g.print();
     }
}

bool QProg::isFlagSet(int o) {
    return flags & o;
}
void QProg::setFlag(int o) {
    flags |= o;
}
void QProg::clearFlag(int o) {
    flags &= ~o;
}

void QProg::addGate(const Gate &g) {
    if(g.op == CZ || g.op == CX) {
        n = MAX(n, MAX(g.a + 1, g.b + 1));
    } else {
        n = MAX(n, g.a + 1);
    }

    gates.push_back(g);
}
void QProg::addGate(Opcode op) {
    gates.push_back(Gate(op));
}
void QProg::addGate(Opcode op, QUBIT a) {
    n = MAX(n, a+1);
    gates.push_back(Gate(op, a, Params(0, n)));
}
void QProg::addGate(Opcode op, QUBIT a, QUBIT b) {
    n = MAX(n, MAX(a+1, b+1));
    gates.push_back(Gate(op, a, b, Params(0, n)));
}
void QProg::addGate(Opcode op, QUBIT a, Outcome o) {
    n = MAX(n, a+1);
    gates.push_back(Gate(op, a, o, Params(0, n)));
}

void QProg::addGate(Opcode op, QUBIT a, QUBIT b, const Params &p) {
    n = MAX(n, MAX(a+1, b+1));
    gates.push_back(Gate(op, a, b, p));
}
void QProg::addGate(Opcode op, QUBIT a, const Params &p) {
    n = MAX(n, a+1);
    gates.push_back(Gate(op, a, p));
}
void QProg::addGate(Opcode op, QUBIT a, Outcome o, const Params &p) {
    n = MAX(n, a+1);
    gates.push_back(Gate(op, a, o, p));
}

void QProg::addMeasurement(QUBIT a, const Params &p, Basis basis) {
    if (basis == XBASIS) {
        addGate(HADAMARD, a, p);
    } else if (basis == YBASIS) {
        addGate(SQRT_X, a, p);
    }

    addGate(MEASURE, a, p);
}

void QProg::addMeasureReset(QUBIT a, const Params &p, Basis basis) {
    if (basis == XBASIS) {
        addGate(HADAMARD, a, p);
    } else if (basis == YBASIS) {
        addGate(SQRT_X, a, p);
    }

    addGate(MEASURE_RESET, a, p);
}

void QProg::addPostselect(QUBIT a, Outcome o, const Params &p, Basis basis) {
    if (basis == XBASIS) {
        addGate(HADAMARD, a, p);
    } else if (basis == YBASIS) {
        addGate(SQRT_X, a, p);
    }

    addGate(POSTSELECT, a, o, p);

    if (basis == XBASIS) {
        addGate(HADAMARD, a, p);
    } else if (basis == YBASIS) {
        addGate(SQRT_X, a, p);
    }
}

/*
Simulate the quantum circuit
*/
void QProg::run() {
    Outcome m; // measurement result
    // int pr; // probability of getting postselected outcome

    if (q == NULL) { // if this is the first time running, we construct the state
        // WARNING: if we apply gates after this point on qubits >= n, it will fail
        if(isFlagSet(NO_DESTABILIZERS)) {
            q = new QState(n, false);
        } else {
            q = new QState(n, true);
        }
        auto_params = new AutoParams(n);

    }

    // Automatic parameterization (on, by default, for now)
    for (Gate &g : gates) {
        const GateDef & gd = GATE_DEFINITIONS[g.op];
        if(gd.args == 1) {
            g.restrictParams(auto_params->getRange(g.a));
        } else {
            g.restrictParams(auto_params->getRange(g.a, g.b));
        }
    }

    for (const Gate &g : gates) {
        switch (g.op) {
            case IDENTITY:
                break;
            case PAULI_X:
                q->PauliX(g.a, g.p);
                break;
            case PAULI_Y:
                q->PauliY(g.a, g.p);
                break;
            case PAULI_Z:
                q->PauliZ(g.a, g.p);
                break;
            case HADAMARD_XY:
                q->hadamard_xy(g.a, g.p);
                break;
            case HADAMARD:
                q->hadamard(g.a, g.p);
                break;
            case HADAMARD_YZ:
                q->hadamard_yz(g.a, g.p);
                break;
            case SQRT_X:
                q->rotX90(g.a, g.p);
                break;
            case SQRT_X_DAG:
                q->rotX270(g.a, g.p);
                break;
            case SQRT_Y:
                q->rotY90(g.a, g.p);
                break;
            case SQRT_Y_DAG:
                q->rotY270(g.a, g.p);
                break;
            case S:
                q->phase(g.a, g.p);
                break;
            case S_DAG:
                q->rotZ270(g.a, g.p);
                break;
            case SWAP:
                q->swapGate(g.a, g.b, g.p);
                break;
            case ISWAP:
                q->iSWAP(g.a, g.b, g.p);
                break;
            case ISWAP_DAG:
                q->iSWAP_dag(g.a, g.b, g.p);
                break;
            case XCX:
                q->xcx(g.a, g.b, g.p);
                break;
            case XCY:
                q->xcy(g.a, g.b, g.p);
                break;
            case XCZ:
                q->xcz(g.a, g.b, g.p);
                break;
            case YCX:
                q->ycx(g.a, g.b, g.p);
                break;
            case YCY:
                q->ycy(g.a, g.b, g.p);
                break;
            case YCZ:
                q->ycz(g.a, g.b, g.p);
                break;
            case CX:
                q->cnot(g.a, g.b, g.p);
                break;
            case CY:
                q->cy(g.a, g.b, g.p);
                break;
            case CZ:
                q->csign(g.a, g.b, g.p);
                break;
            case RESET:
                q->reset(g.a, g.p);
                break;
            case POSTSELECT:
            {
                q->postselect(g.a, g.o, g.p);

                // int pr = q->postselect(g.a, g.o, g.p);
                // if (!isFlagSet(SILENT)) {
                //     cout << "Probability of postselection on qubit " << g.a << ": " << pr/2.0 << endl;
                // }
            }
                break;

            case MEASURE_RESET:
            case MEASURE:

                if (random_bits.empty()) {
                    m = RANDOM;
                } else {
                    m = random_bits.back(); // whether we prefer a 0 or a 1
                    random_bits.pop_back();
                }
                m = q->measure(g.a, m, g.p);

                if (outcomes != NULL)
                    outcomes->push_back(m);

                if (!isFlagSet(SILENT)) {
                    /* CHP-type output */
                    // cout << "Outcome of measuring qubit " << g.a << ": " << m << endl;

                    /* Stim-type output */
                    cout << outcomeToBinary(m);
                }
                if(g.op == MEASURE_RESET) {
                    q->zeroOutPhase(g.a);
                }

                break;
            case DEBUG:
                q->print();
                break;
            default:
                throw domain_error("unimplemented gate");
                break;
        }
    }
    gates.clear(); // remove the gates we've just run
}

void Row::printPhase() const {
    switch (r) {
        case 0: cout << " + "; break;
        case 1: cout << " + i"; break;
        case 2: cout << " - "; break;
        case 3: cout << " - i"; break;
    }
}
void Row::printKet() const {
    cout << "|";
    for(QUBIT i = 0; i < n; i++) cout << (bitX(i) ? 1 : 0);
    cout << ">";
}

// WARNING: Printing requires a Gaussian elimination step, which may ruin Params
void QState::printZBasis() {

    gaussianEliminate();
    unsigned rank = 0; // rank of the Pauli-X matrix for the stabilizers
    for(QUBIT j = 0; j < n; j++) {
        if(stab[rank].bitX(j)) { // this was a pivot row
            rank++;
        }
    }

    Row scratch;  //define a scratch space row
    scratch.n = n;
    TBLX tmp_x[blocks], tmp_z[blocks];
    for(QUBIT i = 0; i < blocks; i++) {
        tmp_x[i].zero();
        tmp_z[i].zero();
    }
    scratch.r = 0;
    scratch.x = tmp_x;
    scratch.z = tmp_z;

    int rx = 0; // x-matrix rank as we traverse stabilizers
    int rz = 0; // z-matrix rank

    // We first determine the lexicographically least basis state
    for(QUBIT j = 0; j < n; j++) {
        if(stab[rx].bitX(j)) { // this was a pivot row
            rx++;
        } else {
            if(stab[rank+rz].phase(Params(0,n)) == 2) {
                scratch.blockX(j) ^= BIT(j);
            }
            rz++;
        }
    }
    if(rank > 31) {
        cout << "Too many states (2^" << rank << ") to print" << endl;
        return;
    }
    scratch.printKet();
    WORD p; // pattern of rows to select to print
    for(WORD r = 1; r < (1ull<<rank); r++) {
        for(WORD i = 0; i < rank; i++) {
            p = r^(r-1);
            if(p&(1ull<<i))
                rowMult(stab[rank-1-i], scratch, Params(0,n), true);
        }
        scratch.printPhase();
        scratch.printKet();
    }

    cout << endl;
}

AutoParams::AutoParams(QUBIT n) : n(n)  {
    heap_size = pow(2, ceil(log(n)/log(2)));
    heap = new Params[heap_size];
    initializeHeap(0, 0, n);
}
void AutoParams::initializeHeap(QUBIT i, QUBIT offset, QUBIT len) {
    if(len == 0) return;
    heap[i] = Params(offset + len/2);
    initializeHeap(leftChild(i), offset, len/2);
    initializeHeap(rightChild(i), offset+len/2+1, (len-1)/2);
}
void AutoParams::printHeap() const {
    for(QUBIT i=0; i<heap_size; i++) {
        cout << "[" << heap[i].lb << ":" << heap[i].ub << "] ";
    } cout << endl;
}

Params& AutoParams::getRange(QUBIT x) const {
    QUBIT i = 0;
    while(true) {
        if(heap[i].contains(x)) {
            return heap[i];
        } else if(x < heap[i].lb) {
            i = leftChild(i);
        } else {
            i = rightChild(i);
        }
    }
}

Params& AutoParams::getRange(QUBIT x, QUBIT y) {
    QUBIT i = 0;
    Params xp = getRange(x);
    Params yp = getRange(y);
    Params p = xp.join(yp);

    while(true) {
        if(heap[i].intersects(p)) {
            heap[i].setParams(p);
            return heap[i];
        } else if(p.ub <= heap[i].lb) {
            i = leftChild(i);
        } else {
            i = rightChild(i);
        }
    }
}
QUBIT AutoParams::parent(QUBIT i) const { return (i-1)/2; }
QUBIT AutoParams::leftChild(QUBIT i) const { return 2*i+1; }
QUBIT AutoParams::rightChild(QUBIT i) const { return 2*i+2; }
