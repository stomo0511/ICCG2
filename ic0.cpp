#include "ic0.hpp"

// ic0.cpp：本体
IC0::IC0(const CRS& A) {
    // 下＋対角を抽出（値も持ってくる）
    CRS Alo = extract_lower_with_diag(A);
    n = Alo.n;
    build_ic0(Alo);
    build_csc_from_crs();
}

CRS IC0::extract_lower_with_diag(const CRS& A) {
    CRS B; B.n = A.n;
    B.rowptr.assign(A.n+1, 0);
    std::vector<int> cols;
    std::vector<double> vals;
    cols.reserve(A.colind.size());
    vals.reserve(A.val.size());
    int nnz = 0;

    for (int i = 0; i < A.n; ++i) {
        B.rowptr[i] = nnz;
        for (int k = A.rowptr[i]; k < A.rowptr[i+1]; ++k) {
            int j = A.colind[k];
            if (j <= i) { // 下三角＋対角のみ
                cols.push_back(j);
                vals.push_back(A.val[k]);
                ++nnz;
            }
        }
    }
    B.rowptr[A.n] = nnz;
    B.colind.swap(cols);
    B.val.swap(vals);
    return B;
}

void IC0::build_ic0(const CRS& Alo) {
    // L のパターン＝Alo のパターン（同一）
    L_rowptr = Alo.rowptr;
    L_colind = Alo.colind;
    L_val    = Alo.val;   // いったん A の値をコピー。上書きで L にしていく。
    n = Alo.n;

    diag_pos.assign(n, -1);
    for (int i = 0; i < n; ++i) {
        for (int k = L_rowptr[i]; k < L_rowptr[i+1]; ++k) {
            if (L_colind[k] == i) { diag_pos[i] = k; break; }
        }
        if (diag_pos[i] < 0)
            throw std::runtime_error("IC(0): diagonal missing at row " + std::to_string(i));
    }

    // 行ごとに L を求める
    for (int i = 0; i < n; ++i) {
        // 1) 行iの各 j<i について L(i,j) を更新
        for (int kij = L_rowptr[i]; kij < L_rowptr[i+1]; ++kij) {
            int j = L_colind[kij];
            if (j >= i)
                break; // 行は昇順なのでここで止められる

            // L(i,0..j-1) と L(j,0..j-1) のドットを「二本ポインタ」で計算
            double s = 0.0;
            int pi = L_rowptr[i], pj = L_rowptr[j];
            const int ei = L_rowptr[i+1], ej = L_rowptr[j+1];

            while (pi < ei && pj < ej) {
                int ci = L_colind[pi];
                int cj = L_colind[pj];
                if (ci >= j || cj >= j)
                    break;      // どちらかが j に到達したら終了（<j のみ対象）
                if (ci == cj) {
                    s += L_val[pi] * L_val[pj]; ++pi; ++pj;
                }
                else if (ci < cj)
                    ++pi;
                else
                    ++pj;
            }

            const int dj = diag_pos[j];
            const double Ljj = L_val[dj];

            if (std::abs(Ljj) < 1e-30)
                throw std::runtime_error("IC(0): zero/near-zero pivot at row " + std::to_string(j));
            L_val[kij] = (L_val[kij] - s) / Ljj; // overwrite with L(i,j)
        }

        // 2) 対角 L(i,i)
        double s = 0.0;
        for (int k = L_rowptr[i]; k < L_rowptr[i+1]; ++k) {
            int j = L_colind[k];
            if (j >= i)
                break;
            double lij = L_val[k];
            s += lij * lij;
        }
        int di = diag_pos[i];
        double diag = L_val[di] - s;
        if (diag <= 0.0)
            throw std::runtime_error("IC(0): non-positive pivot at row " + std::to_string(i));
        L_val[di] = std::sqrt(diag);
    }
}

void IC0::build_csc_from_crs() {
    // L(CRS) → L^T の CSC を構築（同じ値を列構造に並べる）
    LT_colptr.assign(n+1, 0);
    const int nnz = (int)L_colind.size();
    LT_rowind.resize(nnz);
    LT_val.resize(nnz);

    // 列ごとの数をカウント
    for (int k = 0; k < nnz; ++k) {
        int c = L_colind[k];
        ++LT_colptr[c+1];
    }
    // 累積和
    for (int c = 0; c < n; ++c)
        LT_colptr[c+1] += LT_colptr[c];

    // 一時カウンタ
    std::vector<int> ctr = LT_colptr;
    for (int i = 0; i < n; ++i)
    {
        for (int k = L_rowptr[i]; k < L_rowptr[i+1]; ++k)
        {
            int c = L_colind[k];
            int pos = ctr[c]++;
            LT_rowind[pos] = i;        // 列c に行i が来る（= L^T の (c,i)）
            LT_val[pos]    = L_val[k];
        }
    }
}

void IC0::forward_solve(const std::vector<double>& b, std::vector<double>& y) const
{
    // y.assign(n, 0.0);
    if ((int)y.size() != n)
        y.resize(n);    // y.assign(n, 0.0); をやめる

    for (int i = 0; i < n; ++i)
    {
        double s = b[i];
        double di = 1.0;

        for (int k = L_rowptr[i]; k < L_rowptr[i+1]; ++k)
        {
            int j = L_colind[k];
            if (j < i)
                s -= L_val[k] * y[j];
            else if (j == i)
                di = L_val[k];
            else
                break;
        }
        y[i] = s / di;
    }
}

void IC0::backward_solve(const std::vector<double>& y, std::vector<double>& x) const
{
    // x.assign(n, 0.0);
    if ((int)x.size() != n)
        x.resize(n);    // x.assign(n, 0.0); をやめる

    for (int i = n-1; i >= 0; --i)
    {
        double s = y[i];
        double di = 1.0;

        // L^T の行 i は、L の列 i。CSC なら列 i がそのままアクセス可能。
        for (int k = LT_colptr[i]; k < LT_colptr[i+1]; ++k)
        {
            int j = LT_rowind[k];      // これは L^T の列 index に当たる
            if (j > i)
                s -= LT_val[k] * x[j];  // 上三角項
            else if (j == i)
                di = LT_val[k];   // 対角
            // j < i は関与しない
        }
        x[i] = s / di;
    }
}

void IC0::apply(const std::vector<double>& r, std::vector<double>& z) const
{
    std::vector<double> y;
    forward_solve(r, y);          // L y = r
    backward_solve(y, z);         // L^T z = y
}
