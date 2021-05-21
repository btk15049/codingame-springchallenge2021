#ifndef VSCODE
// clang-format off
 #pragma GCC optimize("Ofast")
 #pragma GCC target("sse,sse2,sse3,ssse3,sse4,popcnt,abm,mmx,avx")
 #pragma GCC optimize("O3,omit-frame-pointer,inline")
 #pragma GCC optimize("unroll-loops")
// clang-format on
#endif

#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <sys/types.h>
#include <unistd.h>
#include <vector>
#include <queue>
#include <cstdint>
#include <ctime>
#include <fstream>

using namespace std;

constexpr bool isDebug() { return false; }
string debugInfo;

constexpr array<double, 4> ownedBias = {
#ifdef OWNED_BIAS
    OWNED_BIAS
#else
    0, 1, 2, 3
#endif
};


namespace util {
    using namespace std::chrono;

    class Timer {
      private:
        std::chrono::high_resolution_clock::time_point bg;
        int64_t cache = 0;

      public:
        Timer() { bg = high_resolution_clock::now(); }
        inline int64_t updateAndGet() {
            return cache = duration_cast<milliseconds>(
                               high_resolution_clock::now() - bg)
                               .count();
        }

        inline int64_t getCache() { return cache; }
    };

    namespace xorshift {
        constexpr uint64_t next(uint64_t p) {
            p = p ^ (p << 13);
            p = p ^ (p >> 7);
            return p ^ (p << 17);
        }
    } // namespace xorshift

    template <typename hash_t, size_t size>
    constexpr array<hash_t, size> genHashes(uint64_t seed) {
        array<hash_t, size> ret = {};
        for (size_t i = 0; i < size; i++) {
            ret[i] = seed = xorshift::next(seed);
        }
        return ret;
    }

    template <typename hash_t, size_t size1, size_t size2>
    constexpr array<array<hash_t, size2>, size1> genHashes(uint64_t seed) {
        array<array<hash_t, size2>, size1> ret = {};
        for (size_t i = 0; i < size1; i++) {
            for (size_t j = 0; j < size2; j++) {
                ret[i][j] = seed = xorshift::next(seed);
            }
        }
        return ret;
    }
} // namespace util

constexpr int TIMEOUT_MS = 95;

// constants for beam search
constexpr int MAX_DAY_DEPTH  = 6;
constexpr int MAX_TURN       = 8;
constexpr int STATE_SIZE_CAP = 350000;

// constants for game info
constexpr int MAX_DAY      = 24;
constexpr int CELL_NUM     = 37;
constexpr int NEIGHBOR_NUM = 6;

using hash_t                  = uint32_t;
constexpr uint64_t HASH_SPACE = 1u << 17;
constexpr hash_t HASH_MASK    = HASH_SPACE - 1;

inline double calcOwnedBonus(const array<int, 4>& owned) {
    double ret = 0;
    ret += ownedBias[0] * owned[0];
    ret += ownedBias[1] * owned[1];
    ret += ownedBias[2] * owned[2];
    ret += ownedBias[3] * owned[3];
    return ret;
}

constexpr auto activeHash =
    util::genHashes<hash_t, CELL_NUM>(15002898780217745505llu);
constexpr auto mineHash =
    util::genHashes<hash_t, CELL_NUM>(8300806424620601649llu);
constexpr auto enemiesHash =
    util::genHashes<hash_t, CELL_NUM>(9658931782774167505llu);
constexpr auto treeSizeHash =
    util::genHashes<hash_t, 4, CELL_NUM>(2763545743315817960llu);


namespace stat {
    util::Timer timer;
    int loopCount;
    int expandCount;
    int dedupeCount;

    void init() {
        timer       = util::Timer();
        loopCount   = 0;
        expandCount = 0;
        dedupeCount = 0;
    }

    void show() {
        string msg = "loop_count: " + to_string(loopCount)
                     + "\nexpand_count: " + to_string(expandCount)
                     + "\ndedupe_count: " + to_string(dedupeCount)
                     + "\nelapsed_time: " + to_string(timer.updateAndGet());

        cerr << msg << endl;
    }
} // namespace stat

namespace input {

    struct CellInfo {
        char cellIndex;
        char richness;
        array<char, NEIGHBOR_NUM> neighbors;
        CellInfo() {}
        inline int getRichness() { return richness; }
        inline int getNeighbor(int d) { return neighbors[d]; }
    };
    array<CellInfo, CELL_NUM> cells;

    array<array<int, CELL_NUM>, CELL_NUM> dist;
    array<array<vector<int>, 7>, CELL_NUM> dist2cells;
    array<bitset<CELL_NUM>, CELL_NUM> adj;
    array<bitset<CELL_NUM>, CELL_NUM> possibleShadow;

    void initDistance() {
        for (int i = 0; i < CELL_NUM; i++) {
            for (int j = 0; j < CELL_NUM; j++) {
                dist[i][j] = 100;
            }
        }

        for (int i = 0; i < CELL_NUM; i++) {
            dist[i][i] = 0;
            for (int j : cells[i].neighbors) {
                if (j == -1) continue;
                dist[i][j] = 1;
            }
        }
        for (auto& a : dist2cells) {
            for (auto& b : a) {
                b.clear();
            }
        }

        for (int k = 0; k < CELL_NUM; k++) {
            for (int i = 0; i < CELL_NUM; i++) {
                for (int j = 0; j < CELL_NUM; j++) {
                    dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j]);
                }
            }
        }
        for (int i = 0; i < CELL_NUM; i++) {
            adj[i]            = bitset<CELL_NUM>();
            possibleShadow[i] = bitset<CELL_NUM>();
            for (int j = 0; j < CELL_NUM; j++) {
                dist2cells[i][dist[i][j]].push_back(j);
                if (dist[i][j] == 1) {
                    adj[i][j] = true;
                }
            }
        }
        for (int i = 0; i < CELL_NUM; i++) {
            for (int dir = 0; dir < NEIGHBOR_NUM; dir++) {
                int j = i;
                for (int k = 0; k < 3; k++) {
                    j = cells[j].getNeighbor(dir);
                    if (j == -1) break;
                    possibleShadow[i].set(j);
                }
            }
        }
    }

    istream& operator>>(istream& is, CellInfo& c) {
        int richness; // 0 if the cell is unusable, 1-3 for usable cells
        int neigh0;   // the index of the neighbouring cell for each direction
        int neigh1;
        int neigh2;
        int neigh3;
        int neigh4;
        int neigh5;
        is >> richness >> neigh0 >> neigh1 >> neigh2 >> neigh3 >> neigh4
            >> neigh5;
        c.richness     = richness;
        c.neighbors[0] = neigh0;
        c.neighbors[1] = neigh1;
        c.neighbors[2] = neigh2;
        c.neighbors[3] = neigh3;
        c.neighbors[4] = neigh4;
        c.neighbors[5] = neigh5;
        return is;
    }

    struct TreeInfo {
        bool isMine;    // true if this is your tree
        bool isDormant; // true if this tree is dormant
        char cellIndex; // location of this tree
        char size;      // size of this tree: 0-3
        TreeInfo() {}
        inline int getCellIndex() const { return cellIndex; }
        inline int getSize() const { return size; }
    };
    istream& operator>>(istream& is, TreeInfo& t) {
        int cellIndex;  // location of this tree
        int size;       // size of this tree: 0-3
        bool isMine;    // 1 if this is your tree
        bool isDormant; // 1 if this tree is dormant
        is >> cellIndex >> size >> isMine >> isDormant;
        t.isMine    = isMine;
        t.isDormant = isDormant;
        t.size      = size;
        t.cellIndex = cellIndex;
        return is;
    }

    struct PersonalInfo {
        int sun;
        int score;
        bool isWaiting;
        vector<string> possibleMoves;
    };

    int day;
    int nutrients;
    PersonalInfo me, enemy;
    vector<TreeInfo> trees;

    void initGame(istream& in) {
        int noc;
        in >> noc;
        in.ignore();
        for (int i = 0; i < CELL_NUM; i++) {
            int cellIndex;
            in >> cellIndex;
            cells[cellIndex].cellIndex = cellIndex;
            in >> cells[cellIndex];
        }
    }


    void initTurn(istream& in) {
        in >> day;
        in.ignore();
        stat::init();
        in >> nutrients;
        in.ignore();
        in >> me.sun >> me.score;
        in.ignore();
        in >> enemy.sun >> enemy.score >> enemy.isWaiting;
        in.ignore();
        int numberOfTrees; // the current amount of trees
        in >> numberOfTrees;
        in.ignore();
        trees = vector<TreeInfo>(numberOfTrees);
        for (auto& tree : trees) {
            in >> tree;
            in.ignore();
        }
        int numberOfPossibleMoves;
        in >> numberOfPossibleMoves;
        in.ignore();
        me.possibleMoves = vector<string>(numberOfPossibleMoves);
        for (string& possibleMove : me.possibleMoves) {
            getline(in, possibleMove);
        }
    }
} // namespace input

const string ops[] = {"WAIT", "COMPLETE", "GROW", "SEED"};
struct Operation {
    static constexpr int WAIT     = 0;
    static constexpr int COMPLETE = 1;
    static constexpr int GROW     = 2;
    static constexpr int SEED     = 3;
    int op;
    int arg1;
    int arg2;

    Operation(int op = 0, int arg1 = 0, int arg2 = 0)
        : op(op), arg1(arg1), arg2(arg2) {}

    Operation(const uint8_t* src) {
        op   = static_cast<int>(src[0] >> 6);
        arg1 = static_cast<int>(src[0]) & ((1 << 6) - 1);
        arg2 = static_cast<int>(src[1]);
    }

    string toString() const {
        string ret = ops[op];
        if (op != WAIT) {
            ret += " " + to_string(arg1);
            if (op == SEED) {
                ret += " " + to_string(arg2);
            }
        }
        return ret;
    }

    static Operation wait() { return Operation(WAIT); }
    static Operation complete(int id) { return Operation(COMPLETE, id); }
    static Operation grow(int id) { return Operation(GROW, id); }
    static Operation seed(int from, int to) {
        return Operation(SEED, from, to);
    }
};


constexpr int heuristicBiasByDay[MAX_DAY + 1] = {
    /*0*/ 0,   0,  0,  0,  0,
    /*5*/ 0,   0,  0,  0,  1,
    /*10*/ 1,  2,  3,  4,  5,
    /*15*/ 6,  7,  8,  9,  10,
    /*20*/ 10, 10, 10, 10, 10};

constexpr int growCostOffset[] = {1, 3, 7};
array<hash_t, CELL_NUM> hashCache;

struct State {
    hash_t hash;
    int nutrients;
    int sunPoint;
    int nutrientScore;
    int bonusScore;
    int estimateNutrientScore;
    bitset<MAX_DAY_DEPTH> existCompleteLog;
    Operation firstOperation;
    array<int, 4> owned;
    bitset<CELL_NUM> isMine;
    bitset<CELL_NUM> isEnemies;
    array<uint8_t, MAX_DAY_DEPTH> completeLog;

    struct Tree {
        bool isActive;
        int size;
    };
    array<Tree, CELL_NUM> tree;
    array<bitset<CELL_NUM>, 4> sizeBits;

    inline void addCompleteLog(int day) {
        const int p = nutrients
                      - (heuristicBiasByDay[input::day + day]
                         - heuristicBiasByDay[input::day]);
        estimateNutrientScore += p > 0 ? p : 0;
        existCompleteLog.set(day);
        completeLog[day]++;
    }

    inline void subCompleteLog(int day) {
        const int p = nutrients
                      - (heuristicBiasByDay[input::day + day]
                         - heuristicBiasByDay[input::day]);
        estimateNutrientScore -= p > 0 ? p : 0;
        completeLog[day]--;
        existCompleteLog[day] = completeLog[day] != 0;
    }

    template <typename F>
    inline void iterateMyTree(F f) {
        for (int i = isMine._Find_first(), sz = isMine.size(); i < sz;
             i = isMine._Find_next(i)) {
            f(i);
        }
    }

    template <typename F>
    inline void iterateMyTree(F f) const {
        for (int i = isMine._Find_first(), sz = isMine.size(); i < sz;
             i = isMine._Find_next(i)) {
            f(i);
        }
    }


    inline void nextDay(int day) {
        int gainSunPoint       = 0;
        const int sunDirection = (day + 3) % 6;
        iterateMyTree([&](int id) {
            toggleHashOnlyMine(id);
            tree[id].isActive = true;
            // 日陰でなかったら gainSunPoint を更新
            int gain = tree[id].size;
            for (int i = 1, v = id; i <= 3; i++) {
                v = input::cells[v].getNeighbor(sunDirection);
                if (v == -1) break;
                if (tree[v].size >= i && tree[v].size >= tree[id].size) {
                    gain = 0;
                    break;
                }
            }
            gainSunPoint += gain;
            toggleHashOnlyMine(id);
        });

        // 必要であれば enemy も hash の管理を やる

        sunPoint += gainSunPoint;
    }


    inline void toggleHash(int id) {
        hash ^= tree[id].isActive * activeHash[id];
        hash ^= treeSizeHash[tree[id].size][id];
        hash ^= isMine[id] * mineHash[id];
        hash ^= isEnemies[id] * enemiesHash[id];
    }

    inline void toggleHashOnlyMineUsingCache(int id) { hash ^= hashCache[id]; }

    inline void toggleHashOnlyMineAndSetCache(int id) {
        hashCache[id] = (tree[id].isActive * activeHash[id])
                        ^ (treeSizeHash[tree[id].size][id])
                        ^ (isMine[id] * mineHash[id]);
        hash ^= hashCache[id];
    }

    inline void toggleHashOnlyMine(int id) {
        hash ^= tree[id].isActive * activeHash[id];
        hash ^= treeSizeHash[tree[id].size][id];
        hash ^= isMine[id] * mineHash[id];
    }

    inline void setHashToCacheOnlyMine(int id) {
        hashCache[id] = (tree[id].isActive * activeHash[id])
                        ^ (treeSizeHash[tree[id].size][id])
                        ^ (isMine[id] * mineHash[id]);
    }

    inline void applyGrow(int id) {
        toggleHashOnlyMineUsingCache(id);
        tree[id].isActive = false;
        sunPoint -= growCost(id);
        owned[tree[id].size]--;
        sizeBits[tree[id].size].reset(id);
        tree[id].size++;
        owned[tree[id].size]++;
        sizeBits[tree[id].size].set(id);
        toggleHashOnlyMineAndSetCache(id);
    }

    inline void applyComplete(int id) {
        toggleHashOnlyMineUsingCache(id);
        sunPoint -= completeCost();
        owned[3]--;
        sizeBits[3].reset(id);
        tree[id].isActive = false;
        tree[id].size     = 0;
        isMine[id]        = false;
        nutrientScore += nutrients > 0 ? nutrients : 0;
        bonusScore += (input::cells[id].getRichness() - 1) * 2;
        nutrients--;
    }

    inline void applySeed(int from, int to) {
        toggleHashOnlyMineUsingCache(from);
        tree[from].isActive = false;
        tree[to].isActive   = false;
        sunPoint -= seedCost();
        owned[0]++;
        sizeBits[0].set(to);
        isMine[to] = true;
        toggleHashOnlyMineAndSetCache(from);
        toggleHashOnlyMineAndSetCache(to);
    }

    inline void rollbackGrow(int id) {
        toggleHashOnlyMineUsingCache(id);
        tree[id].isActive = true;
        owned[tree[id].size]--;
        sizeBits[tree[id].size].reset(id);
        tree[id].size--;
        owned[tree[id].size]++;
        sizeBits[tree[id].size].set(id);
        sunPoint += growCost(id);
        toggleHashOnlyMineAndSetCache(id);
    }

    inline void rollbackComplete(int id) {
        owned[3]++;
        sizeBits[3].set(id);
        tree[id].size = 3;
        isMine[id]    = true;
        sunPoint += completeCost();
        tree[id].isActive = true;
        nutrients++;
        nutrientScore -= nutrients > 0 ? nutrients : 0;
        bonusScore -= (input::cells[id].getRichness() - 1) * 2;
        toggleHashOnlyMine(id);
    }

    inline void rollbackSeed(int from, int to) {
        toggleHashOnlyMineUsingCache(from);
        toggleHashOnlyMineUsingCache(to);
        tree[from].isActive = true;
        // tree[to].isActive   = true;
        owned[0]--;
        sizeBits[0].reset(to);
        isMine[to] = false;
        sunPoint += seedCost();
        toggleHashOnlyMineAndSetCache(from);
    }

    inline void apply(Operation op) {
        switch (op.op) {
            case Operation::GROW:
                applyGrow(op.arg1);
                break;
            case Operation::COMPLETE:
                applyComplete(op.arg1);
                break;
            case Operation::SEED:
                applySeed(op.arg1, op.arg2);
                break;
            case Operation::WAIT:
                break;
        }
    }

    inline void rollback(Operation op) {
        switch (op.op) {
            case Operation::GROW:
                rollbackGrow(op.arg1);
                break;
            case Operation::COMPLETE:
                rollbackComplete(op.arg1);
                break;
            case Operation::SEED:
                rollbackSeed(op.arg1, op.arg2);
                break;
            case Operation::WAIT:
                break;
        }
    }

    inline int growCost(int i) const {
        return growCostOffset[tree[i].size] + owned[tree[i].size + 1];
    }

    constexpr int completeCost() const { return 4; }

    inline int seedCost() const { return owned[0]; }

    vector<Operation> getPossibleMoves() const {
        vector<Operation> ret;
        ret.reserve(100);
        for (int i = isMine._Find_first(), sz = isMine.size(); i < sz;
             i = isMine._Find_next(i)) {
            if (tree[i].isActive) {
                if (sunPoint >= seedCost()) {
                    for (int d = 1; d <= tree[i].size; d++) {
                        for (int j : input::dist2cells[i][d]) {
                            if (!isUsed(j)) {
                                ret.emplace_back(Operation::seed(i, j));
                            }
                        }
                    }
                }
                if (tree[i].size == 3) {
                    if (sunPoint >= completeCost()) {
                        ret.emplace_back(Operation::complete(i));
                    }
                }
                else {
                    if (sunPoint >= growCost(i)) {
                        ret.emplace_back(Operation::grow(i));
                    }
                }
            }
        }
        ret.emplace_back(Operation::wait());
        return ret;
    }

    inline bool isUsed(int i) const {
        return (input::cells[i].richness == 0) || isMine[i] || isEnemies[i];
    }

    inline double calcScore(int day, [[maybe_unused]] int turn) const {
        int shadowPenalty = 0;
        // const auto treeBits   = isMine | isEnemies;
        bitset<CELL_NUM> treeBits[4];
        {
            treeBits[3] = isMine;
            auto b      = sizeBits[3];
            treeBits[2] = isMine | (isEnemies & b);
            b |= sizeBits[2];
            treeBits[1] = isMine | (isEnemies & b);
            b |= sizeBits[1];
            treeBits[0] = isMine | (isEnemies & b);
        }
        {
            uint64_t mb = isMine.to_ullong();
            for (int i = 0; i < CELL_NUM; i++, mb >>= 1) {
                if (mb & 1) {
                    shadowPenalty +=
                        (input::possibleShadow[i] & treeBits[tree[i].size])
                            .count();
                }
            }
        }
        constexpr int DAY_CAP = 20;
        const double rDay     = max(DAY_CAP - day, 0);
        // const double rTurn = max(8 - turn, 0);
        return input::me.score
               + day * (estimateNutrientScore + bonusScore) / 24.0
               + sunPoint / 3.0 + rDay * calcOwnedBonus(owned) * 1.0
               - rDay * shadowPenalty / double(DAY_CAP);
    }

    void assertHash() {
        hash_t oldHash = hash;
        hash           = hash_t(0);

        auto t = isMine | isEnemies;
        for (int i = 0; i < CELL_NUM; i++) {
            if (t[i]) {
                toggleHash(i);
            }
        }
        assert(oldHash == hash);
    }

    /* ##################
         Constructors
    ################## */
    State() {}

    State(int nutrients, const input::PersonalInfo& me,
          const vector<input::TreeInfo>& trees)
        : nutrients(nutrients),
          sunPoint(me.sun),
          nutrientScore(0),
          bonusScore(0),
          estimateNutrientScore(0) {
        fill(owned.begin(), owned.end(), 0);
        fill(completeLog.begin(), completeLog.end(), 0);
        for (int i = 0; i < CELL_NUM; i++) {
            tree[i].isActive = false;
            tree[i].size     = 0;
            isMine[i]        = false;
            isEnemies[i]     = false;
        }
        hash = 0;
        for (auto& it : trees) {
            const int id      = it.getCellIndex();
            isMine[id]        = it.isMine;
            isEnemies[id]     = !it.isMine;
            tree[id].isActive = !it.isDormant;
            tree[id].size     = it.getSize();
            sizeBits[it.size].set(id);
            if (it.isMine) {
                owned[it.size]++;
            }
            toggleHash(id);
        }
    }
};
static_assert(sizeof(State) <= 450);
void assertEqual(const State& actual, const State& expected,
                 bool checkOperation = false, bool checkCompleteLog = false) {
    auto as = [](string name, auto a, auto e) {
        if (a != e) {
            cerr << "assert failed!" << endl;
            cerr << "actual." << name << ": " << a << endl;
            cerr << "expected." << name << ": " << e << endl;
            cerr << "debugInfo: " << debugInfo << endl;
            exit(1);
        }
    };
    as("sunPoint", actual.sunPoint, expected.sunPoint);
    as("nutrientScore", actual.nutrientScore, expected.nutrientScore);
    as("nutrientScore", actual.bonusScore, expected.bonusScore);
    as("nutrients", actual.nutrients, expected.nutrients);
    for (int i = 0; i < CELL_NUM; i++) {
        auto id = to_string(i);
        as(id + ".size", actual.tree[i].size, expected.tree[i].size);
        as(id + ".isActive", actual.tree[i].isActive,
           expected.tree[i].isActive);
        as(id + ".isMine", actual.isMine[i], expected.isMine[i]);
        as(id + ".isEnemies", actual.isEnemies[i], expected.isEnemies[i]);
        as(id + ".sizeBits[0]", actual.sizeBits[0][i], actual.sizeBits[0][i]);
        as(id + ".sizeBits[1]", actual.sizeBits[1][i], actual.sizeBits[1][i]);
        as(id + ".sizeBits[2]", actual.sizeBits[2][i], actual.sizeBits[2][i]);
        as(id + ".sizeBits[3]", actual.sizeBits[3][i], actual.sizeBits[3][i]);
    }
    for (int i = 0; i < 4; i++) {
        auto id = to_string(i);
        as(id + ".owned", actual.owned[i], expected.owned[i]);
    }
    if (checkOperation) {
        as("operation.op", actual.firstOperation.op,
           expected.firstOperation.op);
        as("operation.arg1", actual.firstOperation.arg1,
           expected.firstOperation.arg1);
        as("operation.arg2", actual.firstOperation.arg2,
           expected.firstOperation.arg2);
    }
    if (checkCompleteLog) {
        for (int i = 0; i < MAX_DAY_DEPTH; i++) {
            auto id = to_string(i);
            as("completeLog." + id, int(actual.completeLog[i]),
               int(expected.completeLog[i]));
            as("existCompleteLog." + id, actual.existCompleteLog[i],
               expected.existCompleteLog[i]);
        }
    }
}

struct Result {
    double score;
    Operation op;
    Result(double score, Operation op) : score(score), op(op) {}
};

bool operator<(const Result& lhs, const Result& rhs) {
    return lhs.score < rhs.score;
}

void debugPossibleMove(const State& s) {
    vector<string> pm;
    for (auto it : s.getPossibleMoves()) {
        pm.push_back(it.toString());
    }

    sort(pm.begin(), pm.end());
    sort(input::me.possibleMoves.begin(), input::me.possibleMoves.end());

    cerr << "actual_possibleMoves:" << endl;
    for (auto& it : pm) {
        cerr << "- " << it << endl;
    }
    cerr << "expected_possibleMoves:" << endl;
    for (auto& it : input::me.possibleMoves) {
        cerr << "- " << it << endl;
    }
    assert(input::me.possibleMoves == pm);
}


namespace Queue {
    vector<State> st;
    struct StateId {
        hash_t hash;
        int id;
        double score;
        StateId(hash_t hash, int id = 0, double score = 0)
            : hash(hash), id(id), score(score) {}
    };
    bool operator<(StateId lhs, StateId rhs) { return lhs.score < rhs.score; }

    using queue_t = priority_queue<StateId, vector<StateId>>;
    array<queue_t, MAX_DAY_DEPTH * MAX_TURN> que;

    array<array<double, HASH_SPACE>, MAX_DAY_DEPTH * MAX_TURN> hashes;
    bitset<MAX_DAY_DEPTH * MAX_TURN> exist;

    void initTurn() {
        util::Timer timer;
        for (auto& it : que) {
            vector<StateId> c;
            c.reserve(5000);
            it = queue_t(less<StateId>(), move(c));
        }
        exist.reset();
        st.clear();
        st.reserve(STATE_SIZE_CAP);
        for (int i = 0; i < (int)hashes.size(); i++) {
            fill(hashes[i].begin(), hashes[i].end(), 0.0);
        }
        cerr << "clean: " << timer.updateAndGet() << endl;
    }

    bool push(int day, int turn, const State& state, double score) {
        const auto dayTurn      = (day - input::day) * MAX_TURN + turn;
        const hash_t maskedHash = state.hash & HASH_MASK;

        if (hashes[dayTurn][maskedHash] < score) {
            exist[dayTurn] = true;
            que[dayTurn].emplace(maskedHash, st.size(), score);
            hashes[dayTurn][maskedHash] = score;
            st.push_back(state);
        }
        else {
            stat::dedupeCount++;
            // cerr << "deduped: " << maskedHash << " " << score << " "
            //      << hashes[dayTurn][maskedHash] << endl;
        }
        return false;
    }

    inline State& pop(int day, int turn) {
        const int dayTurn = (day - input::day) * MAX_TURN + turn;
        const int id      = que[dayTurn].top().id;
        que[dayTurn].pop();
        exist[dayTurn] = !que[dayTurn].empty();
        return st[id];
    }

    inline pair<int, int> getNext(int day, int turn) {
        int ret = exist._Find_next((day - input::day) * MAX_TURN + turn);
        if (ret == (int)exist.size()) {
            ret = exist._Find_first();
        }
        if (ret == (int)exist.size()) {
            return {-1, -1};
        }
        day     = ret / MAX_TURN + input::day;
        turn    = ret % MAX_TURN;
        auto& q = que[ret];
        while (!q.empty()) {
            if (hashes[ret][q.top().hash] != q.top().score) {
                stat::dedupeCount++;
                q.pop();
                exist[ret] = !q.empty();
            }
            else {
                return {day, turn};
            }
        }
        return getNext(day, turn);
    }


} // namespace Queue

Result solve() {
    State today(input::nutrients, input::me, input::trees);
    if (isDebug()) {
        debugPossibleMove(today);
    }

    Queue::push(input::day, 0, today, 1);

    Result best(0, Operation::wait());
    const int dayCap = min(24, input::day + MAX_DAY_DEPTH);
    int lastDay      = 0;
    int lastTurn     = 0;
    while (stat::timer.updateAndGet() < TIMEOUT_MS) {
        auto [d, turn] = Queue::getNext(lastDay, lastTurn);
        lastDay        = d;
        lastTurn       = turn;
        if (d == -1) {
            cerr << "all route were searched" << endl;
            break;
        }
        stat::loopCount++;
        State& cur = Queue::pop(d, turn);

        // expand
        const auto ops = cur.getPossibleMoves();

        cur.iterateMyTree([&](int v) { cur.setHashToCacheOnlyMine(v); });

        // State copied = cur;

        //  bool hasNoCostSeed = false;
        if (turn + 1 < MAX_TURN) {
            for (auto& op : ops) {
                // if (stat::timer.updateAndGet() >= TIMEOUT_MS) break;

                if (op.op == Operation::WAIT) {
                    continue;
                }
                if (op.op == Operation::SEED) {
                    if ((cur.isMine & input::adj[op.arg2]).any()) {
                        continue;
                    }
                }
                if (op.op == Operation::COMPLETE) {
                    if (cur.completeLog[d - input::day] != turn) {
                        continue;
                    }
                }
                stat::expandCount++;


                if (d == input::day && turn == 0) {
                    cur.firstOperation = op;
                }
                if (op.op == Operation::COMPLETE) {
                    cur.addCompleteLog(d - input::day);
                }

                cur.apply(op);
                Queue::push(d, turn + 1, cur, cur.calcScore(d, turn + 1));
                cur.rollback(op);
                if (op.op == Operation::COMPLETE) {
                    cur.subCompleteLog(d - input::day);
                }
            }
        }
        if (d == input::day && turn == 0) {
            cur.firstOperation = Operation::wait();
        }
        if (d + 1 < dayCap) {
            cur.nextDay(d + 1);
            Queue::push(d + 1, 0, cur, cur.calcScore(d + 1, 0));
        }
        else {
            const Result x(cur.calcScore(d + 1, 0), cur.firstOperation);
            if (best < x) {
                best = x;
            }
        }
        // }
    }
    return best;
}

int getHourAsJST() {
    auto currentUnixTime = time(0);
    currentUnixTime /= 60 * 60; // 時刻に変換
    currentUnixTime += 9;       // to JST
    currentUnixTime %= 24;
    return currentUnixTime;
}

string getComment(int h) {
    string comment;
    switch (h) {
        case 4:
        case 5:
        case 6:
        case 7:
        case 8:
        case 9:
            comment = "Ohayochiwawa";
            break;
        case 10:
        case 11:
        case 12:
        case 13:
        case 14:
        case 15:
        case 16:
        case 17:
            comment = "Konnichiwawa";
            break;
        default:
            comment = "Konbanchiwawa";
            break;
    }
    return comment;
}

void run() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    int h = getHourAsJST();
    cerr << "currentHour: " << h << endl;
    const string comment = getComment(h);

    input::initGame(cin);
    input::initDistance();

    State lastState;


    while (1) {
        Queue::initTurn();
        input::initTurn(cin);
        cerr << "input: " << stat::timer.updateAndGet() << endl;
        auto result = solve();
        cout << result.op.toString() << " " << comment << endl;
        stat::show();
        cerr << "estimated_score: " << result.score << endl;

        lastState = State(input::nutrients, input::me, input::trees);
        lastState.firstOperation = result.op;
    }
}

#ifndef TEST
int main() { run(); }
#endif
