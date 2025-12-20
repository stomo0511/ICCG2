#ifndef CRS_HPP_
#define CRS_HPP_

#include <vector>
#include <string>

struct CRS {
    int n;                        // matrix size (n x n)
    std::vector<int> rowptr;      // size n+1
    std::vector<int> colind;      // size nnz
    std::vector<double> val;      // size nnz
};

// Matrix Market (*.mtx) を読み込んで CRS を返す（general/symmetric, real/integer/pattern）
CRS read_mm2crs(const std::string& filepath);

#endif // CRS_HPP_