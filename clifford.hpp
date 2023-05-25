#ifndef _CLIFFORD_HPP_
#define _CLIFFORD_HPP_

#include <vector>
#include <string>
#include <unordered_map>
#include <bit>

#ifdef __GNUC__
#define attribute(x) __attribute__((x));
#else
#define attribute(x)
#endif

#define MIN(a, b) (((a) < (b)) ? (a) : (b))
#define MAX(a, b) (((a) > (b)) ? (a) : (b))

typedef int QUBIT;
typedef char SIGNBITS;

#ifdef __x86_64__
#include <immintrin.h>
#define SHIFT 8 // should be 6 for 64 bit architecture, 8 for 256 bit SIMD

struct alignas(16) simd {
     union {
         __m256i val;
         __m128i u128[2];
         uint64_t u64[4];
         uint8_t u8[32];
     };

     inline simd() : val(__m256i{}) {}
     inline simd(__m256i val) : val(val) {}

     inline simd& operator^=(const simd &a) {
         val = _mm256_xor_si256(val, a.val);
         return *this;
     }
     inline simd& operator&=(const simd &a) {
         val = _mm256_and_si256(val, a.val);
         return *this;
     }
     inline simd& operator|=(const simd &a) {
         val = _mm256_or_si256(val, a.val);
         return *this;
     }
     inline simd operator^(const simd &a) const {
         return { _mm256_xor_si256(val, a.val) };
     }
     inline simd operator&(const simd &a) const {
         return { _mm256_and_si256(val, a.val) };
     }
     inline simd operator|(const simd &a) const {
         return { _mm256_or_si256(val, a.val) };
     }
     inline void zero() {
         u64[0] = u64[1] = u64[2] = u64[3] = 0;
     }
     inline uint64_t popcount() const {
         return _mm_popcnt_u64(u64[0]) +
                _mm_popcnt_u64(u64[1]) +
                _mm_popcnt_u64(u64[2]) +
                _mm_popcnt_u64(u64[3]);
     }

     inline uint64_t parity() const {
         return popcount() % 2;
     }
 };

#else

#define SHIFT 6
struct alignas(16) simd {
    union{
      uint64_t val;
      uint64_t u64[1];
    };

     inline simd() : val(0) {}
     inline simd(uint64_t val) : val(val) {}

     inline simd& operator^=(const simd &a) {
         val ^= a.val;
         return *this;
     }
     inline simd& operator&=(const simd &a) {
         val &= a.val;
         return *this;
     }
     inline simd& operator|=(const simd &a) {
         val |= a.val;
         return *this;
     }
     inline simd operator^(const simd &a) const {
         return { val^a.val };
     }
     inline simd operator&(const simd &a) const {
         return { val&a.val };
     }
     inline simd operator|(const simd &a) const {
         return { val | a.val };
     }
     inline void zero() {
         val = 0;
     }

     inline uint64_t popcount() const {
         return std::__popcount(val);
     }

     inline uint64_t parity() const {
         return popcount() % 2;
     }
 };
#endif
#define INTS_PER_BLOCK ((uint64_t)(pow(2, SHIFT-6)))




enum Basis { XBASIS = 0, YBASIS, ZBASIS, UNSPECIFIED } attribute(packed);
enum Opcode {
    NO_GATE = 0,
    MEASURE,
    MEASURE_RESET,
    RESET,
    POSTSELECT,
    IDENTITY,
    PAULI_X,
    PAULI_Y,
    PAULI_Z,
    HADAMARD_XY,
    HADAMARD,
    HADAMARD_YZ,
    SQRT_X,
    SQRT_X_DAG,
    SQRT_Y,
    SQRT_Y_DAG,
    S,
    S_DAG,
    SWAP,
    ISWAP,
    ISWAP_DAG,
    XCX,
    XCY,
    XCZ,
    YCX,
    YCY,
    YCZ,
    CX,
    CY,
    CZ,
    DEBUG } attribute(packed);

enum Outcome { ZERO = 0, ONE, ZERO_RAND, ONE_RAND, RANDOM } attribute(packed);
int outcomeToBinary(Outcome o);

std::ostream& operator<<(std::ostream& out, Outcome o);

const int DISPTIME            = 0x1;
const int SILENT              = 0x2;
const int FIXEDMEASUREMENTS   = 0x4;
const int LACONIC             = 0x8;
const int NO_DESTABILIZERS    = 0x10;


typedef uint64_t WORD;
typedef simd TBLX;

#define ENTRY_SIZE (1 << 6)
#define MASK (ENTRY_SIZE - 1)
#define BIT(i) (1ull << ((i)&MASK))

// Relative frequencies (should add to 1)
#define XFREQ 0.5
#define YFREQ 0.5
#define ZFREQ 0

#define NO_QUBIT 0 // For the second qubit in single qubit gates

// These parameters specific which elements of the tableau should be accessed for a given operation
// The operation will be applied within the two intervals specified by a lower bound (lb1, lb2) and an upper bound (ub1, ub2)
// The program is only guaranteed to be correct when all of the following conditions hold:
// 0 <= lb1 <= ub1 <= lb2 <= ub2 <= elements of the tableau
class Params {
  public:
    QUBIT lb;  // lower bound
    QUBIT ub;  // upper bound

    void setParams(QUBIT lb);
    void setParams(QUBIT lb, QUBIT ub);
    void setParams(const Params& p);

    Params() : lb(0), ub(0) {};
    Params(QUBIT lb, QUBIT ub) : lb(lb), ub(ub) {};
    Params(const Params &p) : lb(p.lb), ub(p.ub) {};
    Params(QUBIT a) : lb(a), ub(a+1) {};

    QUBIT lbs() const {
        return lb >> SHIFT;
    };
    QUBIT ubs() const {
        return 1 + ((ub-1) >> SHIFT);
    };
    bool contains(QUBIT i) const {
        return lb <= i && i < ub;
    };
    bool intersects(const Params& p) const {
        if(p.lb < lb) return lb < p.ub;
        else if(lb < p.lb) return p.lb < ub;
        else return true;
    };
    Params join(const Params& p) const {
        return Params(std::min(p.lb, lb), std::max(p.ub, ub));
    }
};
/* Binary heap used for automatic params */
class AutoParams {
public:
    QUBIT n;
    Params& getRange(QUBIT x) const;
    Params& getRange(QUBIT x, QUBIT y);

    AutoParams(QUBIT n);

private:
    Params* heap;
    QUBIT heap_size;

    QUBIT parent(QUBIT i) const;
    QUBIT leftChild(QUBIT i) const;
    QUBIT rightChild(QUBIT i) const;
    void initializeHeap(QUBIT i, QUBIT offset, QUBIT len);
    void printHeap() const;
};


/* The definition of a quantum gate for the purposes of parsing and display */
class GateDef {
  public:
    std::string name;
    Opcode op;
    int args;
    bool takesParams;

    GateDef(std::string name, Opcode op, int args, bool takesParams = false) :
        name(name),
        op(op),
        args(args),
        takesParams(takesParams){};
};


static const GateDef GATE_DEFINITIONS [] = {
        GateDef("NO_GATE",   NO_GATE,      0, false),
        GateDef("M",         MEASURE,      1, true),
        GateDef("MR",        MEASURE_RESET,1, true),
        GateDef("R",         RESET,        1, true),
        GateDef("F",         POSTSELECT,   1, true),
        GateDef("I",         IDENTITY,     1, true),
        GateDef("X",         PAULI_X,      1, true),
        GateDef("Y",         PAULI_Y,      1, true),
        GateDef("Z",         PAULI_Z,      1, true),
        GateDef("H_XY",      HADAMARD_XY,  1, true),
        GateDef("H",         HADAMARD,     1, true),
        GateDef("H_YZ",      HADAMARD_YZ,  1, true),
        GateDef("SQRT_X",    SQRT_X,       1, true),
        GateDef("SQRT_X_DAG",SQRT_X_DAG,   1, true),
        GateDef("SQRT_Y",    SQRT_Y,       1, true),
        GateDef("SQRT_Y_DAG",SQRT_Y_DAG,   1, true),
        GateDef("S",         S,            1, true),
        GateDef("S_DAG",     S_DAG,        1, true),
        GateDef("SWAP",      SWAP,         2, true),
        GateDef("ISWAP",     ISWAP,        2, true),
        GateDef("ISWAP_DAG", ISWAP_DAG,    2, true),
        GateDef("XCX",       XCX,          2, true),
        GateDef("XCY",       XCY,          2, true),
        GateDef("XCZ",       XCZ,          2, true),
        GateDef("YCX",       YCX,          2, true),
        GateDef("YCY",       YCY,          2, true),
        GateDef("YCZ",       YCZ,          2, true),
        GateDef("CX",        CX,           2, true),
        GateDef("CY",        CY,           2, true),
        GateDef("CZ",        CZ,           2, true),
        GateDef("DEBUG",     DEBUG,        0, true)
    };

/* An instance of a quantum gate in a circuit */
class Gate {
  public:
    QUBIT a, b;  // Gate acts on qubits a and b
    Opcode op;
    Outcome o;   // postselected outcome for postselection
    Params p;
    void print() const;
    void restrictParams(const Params& p);

    Gate(Opcode op) : a(NO_QUBIT), b(NO_QUBIT), op(op), o(RANDOM), p(){};
    Gate(Opcode op, QUBIT a) : a(a), b(NO_QUBIT), op(op), o(RANDOM), p(){};
    Gate(Opcode op, QUBIT a, QUBIT b) : a(a), b(b), op(op), o(RANDOM), p(){};
    Gate(Opcode op, QUBIT a, Outcome o) : a(a), b(NO_QUBIT), op(op), o(o), p(){};

    Gate(Opcode op, QUBIT a, QUBIT b, const Params &p) : a(a), b(b), op(op), o(RANDOM), p(p){};
    Gate(Opcode op, QUBIT a, const Params &p) : a(a), b(NO_QUBIT), op(op), o(RANDOM), p(p){};
    Gate(Opcode op, QUBIT a, Outcome o, const Params &p) : a(a), b(NO_QUBIT), op(op), o(o), p(p){};
};

bool operator==(const Gate& g, const Opcode op);

struct Row {
    QUBIT n; // unfortunate hack to support printing -- consider removing later
    TBLX *x;
    TBLX *z;
    char r;

    const WORD & blockX(QUBIT i) const;
    WORD & blockX(QUBIT i);
    const WORD & blockZ(QUBIT i) const;
    WORD & blockZ(QUBIT i);

    WORD bitX(QUBIT i) const;
    WORD bitZ(QUBIT i) const;

    char pauli(QUBIT i) const;
    SIGNBITS phase() const;
    SIGNBITS phase(const Params &p) const;

    void leftMultBy(const Row &other, const Params &p, bool withPhase);

    void printPhase() const;
    void printKet() const;
};

/* Quantum State */
class QState {
  public:
    const QUBIT n;      // n qubits
    const QUBIT blocks; // ceil(n/ENTRY_SIZE)

    void PauliX(QUBIT a);
    void PauliX(QUBIT a, const Params &p);
    void PauliY(QUBIT a);
    void PauliY(QUBIT a, const Params &p);
    void PauliZ(QUBIT a);
    void PauliZ(QUBIT a, const Params &p);

    void hadamard_xy(QUBIT a);
    void hadamard_xy(QUBIT a, const Params &p);
    void hadamard(QUBIT a);
    void hadamard(QUBIT a, const Params &p);
    void hadamard_yz(QUBIT a);
    void hadamard_yz(QUBIT a, const Params &p);

    void rotX90(QUBIT a);
    void rotX90(QUBIT a, const Params &p);
    void rotX270(QUBIT a);
    void rotX270(QUBIT a, const Params &p);
    void rotY90(QUBIT a);
    void rotY90(QUBIT a, const Params &p);
    void rotY270(QUBIT a);
    void rotY270(QUBIT a, const Params &p);
    void phase(QUBIT a);
    void phase(QUBIT a, const Params &p);
    void rotZ90(QUBIT a);
    void rotZ90(QUBIT a, const Params &p);
    void rotZ270(QUBIT a);
    void rotZ270(QUBIT a, const Params &p);

    void swapGate(QUBIT a, QUBIT b);
    void swapGate(QUBIT a, QUBIT b, const Params &p);
    void iSWAP(QUBIT a, QUBIT b);
    void iSWAP(QUBIT a, QUBIT b, const Params &p);
    void iSWAP_dag(QUBIT a, QUBIT b);
    void iSWAP_dag(QUBIT a, QUBIT b, const Params &p);

    void xcx(QUBIT a, QUBIT b);
    void xcx(QUBIT a, QUBIT b, const Params &p);
    void xcy(QUBIT a, QUBIT b);
    void xcy(QUBIT a, QUBIT b, const Params &p);
    void xcz(QUBIT a, QUBIT b);
    void xcz(QUBIT a, QUBIT b, const Params &p);

    void ycx(QUBIT a, QUBIT b);
    void ycx(QUBIT a, QUBIT b, const Params &p);
    void ycy(QUBIT a, QUBIT b);
    void ycy(QUBIT a, QUBIT b, const Params &p);
    void ycz(QUBIT a, QUBIT b);
    void ycz(QUBIT a, QUBIT b, const Params &p);

    void cx(QUBIT a, QUBIT b);
    void cx(QUBIT a, QUBIT b, const Params &p);
    void cnot(QUBIT a, QUBIT b);
    void cnot(QUBIT a, QUBIT b, const Params &p);
    void cy(QUBIT a, QUBIT b);
    void cy(QUBIT a, QUBIT b, const Params &p);
    void cz(QUBIT a, QUBIT b);
    void cz(QUBIT a, QUBIT b, const Params &p);
    void csign(QUBIT a, QUBIT b);
    void csign(QUBIT a, QUBIT b, const Params &p);

    void reset(QUBIT a);
    void reset(QUBIT a, const Params &p);
    Outcome measure(QUBIT a, Outcome o, const Params &p);
    int postselect(QUBIT a, Outcome o, const Params &p);
    void zeroOutPhase(QUBIT a);
    void gaussianEliminate();
    QUBIT gaussianEliminate(QUBIT a, const Params &p);

    void print() const;
    void printZBasis();

    QState(QUBIT n, bool destabilizers=true);

  private:
    Row * stab;
    Row * dest;
    std::vector<Row*> rows;
    TBLX *master;

    void rowCopy(const Row & a, Row & b, const Params &p);
    void rowSwap(Row & a, Row & b, const Params &p);
    void rowZero(Row & a, const Params &p);
    void rowSetX(Row & a, QUBIT i, const Params &p);
    void rowSetZ(Row & a, QUBIT i, const Params &p);

    void rowMult(const Row & a, Row & b, const Params &p, bool withPhase);
    // void rowMultBatch(Row *s, Row* t, QUBIT piv, QUBIT* idxs, int size, const Params &p, bool withPhase);
    void rowMultBatch(Row *s, Row* t, QUBIT piv, QUBIT* idxs, int size, const Params &p);
};

/* Quantum circuit */
class QProg {
  public:
    QUBIT n;
    std::vector<Gate> gates;
    std::vector<Outcome> random_bits;
    std::vector<Outcome> *outcomes;
    QState *q;
    AutoParams *auto_params;

  private:
    int flags;

  public:
    bool isFlagSet(int o);
    void setFlag(int o);
    void clearFlag(int o);

    void printCHP() const;
    void print() const;

    void addGate(const Gate &g);
    void addGate(Opcode op);
    void addGate(Opcode op, QUBIT a);
    void addGate(Opcode op, QUBIT a, QUBIT b);
    void addGate(Opcode op, QUBIT a, Outcome o);
    void addGate(Opcode op, QUBIT a, QUBIT b, const Params &p);
    void addGate(Opcode op, QUBIT a, const Params &p);
    void addGate(Opcode op, QUBIT a, Outcome o, const Params &p);


    void addMeasurement(QUBIT a, const Params &p, Basis basis);
    void addMeasurement(QUBIT a, Outcome o, const Params &p, Basis basis);
    void addMeasureReset(QUBIT a, const Params &p, Basis basis);
    void addPostselect(QUBIT a, Outcome o, const Params &p, Basis basis);
    void run();

    QProg() : n(0), outcomes(NULL), q(NULL), flags(0) {};
};



#endif
