// tdma.h
#ifndef TDMA_H
#define TDMA_H

#include <vector>

// Solves a tridiagonal linear system using the TDMA (Thomas algorithm).
// a: sub-diagonal (a[0] unused or zero)
// b: main diagonal
// c: super-diagonal (c[n-1] unused or zero)
// d: right-hand side
std::vector<double> solveTridiagonal(const std::vector<double>& a,
    const std::vector<double>& b,
    const std::vector<double>& c,
    const std::vector<double>& d);

#endif // TDMA_H
