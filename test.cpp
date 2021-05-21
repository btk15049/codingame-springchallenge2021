#define TEST

#include "main.cpp"

void test1_1(State A, State B) {
    const string testName = "1_1";
    cerr << "-- Run test '" + testName + "'" << endl;

    for (int id = 1; id <= 4; id++) {
        A.iterateMyTree([&](int i) { A.setHashToCacheOnlyMine(i); });
        A.applyGrow(id);
    }

    for (int id = 4; id >= 1; id--) {
        B.iterateMyTree([&](int i) { B.setHashToCacheOnlyMine(i); });
        B.applyGrow(id);
    }

    assert(A.hash == B.hash);

    A.iterateMyTree([&](int i) { A.setHashToCacheOnlyMine(i); });
    A.nextDay(input::day + 1);
    B.iterateMyTree([&](int i) { B.setHashToCacheOnlyMine(i); });
    B.nextDay(input::day + 1);
    assert(A.hash == B.hash);

    cerr << "-- Finish test '" + testName + "'" << endl;
}

void test1_2(State A, State B) {
    const string testName = "1_2";
    cerr << "-- Run test '" + testName + "'" << endl;

    for (int id = 1; id <= 4; id++) {
        A.iterateMyTree([&](int i) { A.setHashToCacheOnlyMine(i); });
        A.applyGrow(id);
    }

    A.iterateMyTree([&](int i) { A.setHashToCacheOnlyMine(i); });
    A.nextDay(input::day + 1);

    B.iterateMyTree([&](int i) { B.setHashToCacheOnlyMine(i); });
    B.nextDay(input::day + 1);


    for (int id = 4; id >= 1; id--) {
        B.iterateMyTree([&](int i) { B.setHashToCacheOnlyMine(i); });
        B.applyGrow(id);
    }

    A.iterateMyTree([&](int i) { A.setHashToCacheOnlyMine(i); });
    A.nextDay(input::day + 1);
    B.iterateMyTree([&](int i) { B.setHashToCacheOnlyMine(i); });
    B.nextDay(input::day + 1);
    assert(A.hash == B.hash);

    cerr << "-- Finish test '" + testName + "'" << endl;
}

void test1_3(State A, State B) {
    const string testName = "1_3";
    cerr << "-- Run test '" + testName + "'" << endl;

    A.iterateMyTree([&](int i) { A.setHashToCacheOnlyMine(i); });
    for (int to : {10, 11, 12}) {
        A.applySeed(1, to);
        A.rollbackSeed(1, to);
    }
    A.applySeed(1, 0);

    B.iterateMyTree([&](int i) { B.setHashToCacheOnlyMine(i); });
    for (int to : {13, 11, 12}) {
        B.applySeed(1, to);
        B.rollbackSeed(1, to);
    }
    B.applySeed(1, 0);

    assertEqual(A, B);

    A.iterateMyTree([&](int i) { A.setHashToCacheOnlyMine(i); });
    A.nextDay(input::day + 1);
    B.iterateMyTree([&](int i) { B.setHashToCacheOnlyMine(i); });
    B.nextDay(input::day + 1);

    assertEqual(A, B);

    cerr << "-- Finish test '" + testName + "'" << endl;
}

void test1_4(State A, State B) {
    const string testName = "1_4";
    cerr << "-- Run test '" + testName + "'" << endl;

    A.iterateMyTree([&](int i) { A.setHashToCacheOnlyMine(i); });
    for (int id : {5, 6}) {
        A.applyComplete(id);
        A.rollbackComplete(id);
    }
    for (int to : {10, 11, 12}) {
        A.applySeed(1, to);
        A.rollbackSeed(1, to);
    }
    A.applyComplete(7);

    B.iterateMyTree([&](int i) { B.setHashToCacheOnlyMine(i); });
    for (int id : {6, 5}) {
        B.applyComplete(id);
        B.rollbackComplete(id);
    }
    for (int to : {13, 11, 12}) {
        B.applySeed(1, to);
        B.rollbackSeed(1, to);
    }
    B.applyComplete(7);

    assertEqual(A, B);

    A.iterateMyTree([&](int i) { A.setHashToCacheOnlyMine(i); });
    A.nextDay(input::day + 1);
    B.iterateMyTree([&](int i) { B.setHashToCacheOnlyMine(i); });
    B.nextDay(input::day + 1);

    assertEqual(A, B);

    cerr << "-- Finish test '" + testName + "'" << endl;
}

void test1() {
    cerr << "- Run test 1" << endl;
    {
        ifstream fin("game.in");
        input::initGame(fin);
        fin.close();
    }
    input::initDistance();
    {
        ifstream fin("turn1.in");
        input::initTurn(fin);
        fin.close();
    }
    State A(input::nutrients, input::me, input::trees);
    State B(input::nutrients, input::me, input::trees);

    test1_1(A, B);
    test1_2(A, B);
    test1_3(A, B);
    test1_4(A, B);

    cerr << "- Finish test 1" << endl;
}

int main() { test1(); }
