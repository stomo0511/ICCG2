#pragma once

#include <vector>
#include <string>
#include <stdexcept>
#include <cmath>
#include "crs.hpp"

// ic0.hpp 
struct IC0
{
    // 係数の安全マージン（負のsqrt回避用）：A_ii - Σ L^2 + shift を使う
    explicit IC0(const CRS& A);

    // z = (LL^T)^{-1} r  を適用
    void apply(const std::vector<double>& r, std::vector<double>& z) const;

    std::string name() const { return "IC preconditioner"; }

private:
    int n = 0;
    // L（下三角＋対角）のCRS
    std::vector<int> L_rowptr, L_colind;
    std::vector<double> L_val;

    // L^T 用 CSC（列ポインタ配列）
    std::vector<int> LT_colptr, LT_rowind;
    std::vector<double> LT_val;

    // 対角位置（Lの各行の対角要素位置インデックス）
    std::vector<int> diag_pos;

    static CRS extract_lower_with_diag(const CRS& A);   // A の下三角＋対角のパターン・値を取り出し
    void build_ic0(const CRS& Alo);                     // Alo（下＋対角）から IC(0) を構築
    void build_csc_from_crs();                          // L(CRS) → L^T(CSC)
    void forward_solve(const std::vector<double>& b, std::vector<double>& y) const;  // Ly=b
    void backward_solve(const std::vector<double>& y, std::vector<double>& x) const; // L^T x=y
};
