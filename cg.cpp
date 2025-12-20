//
// CRS形式の疎行列に対するJacobi前処理付き共役勾配(CG)法
//
#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>
#include <map>
#include <set>
#include <queue>
#include <cmath>
#include <limits>

#include <omp.h>
#include "crs.hpp"

using namespace std;

// 疎行列ベクトル積
// y = A*x
static void spmv(const CRS& A, const vector<double>& x, vector<double>& y) {
    int n = A.n;
    y.assign(n, 0.0);
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int k = A.rowptr[i]; k < A.rowptr[i+1]; ++k) {
            sum += A.val[k] * x[A.colind[k]];
        }
        y[i] = sum;
    }
}

// Jacobi preconditioner: z = M^{-1} r with M = diag(A)
struct Jacobi {
    vector<double> inv_diag;
    explicit Jacobi(const CRS& A) {
        inv_diag.assign(A.n, 1.0);
        for (int i = 0; i < A.n; ++i) {
            double d = 0.0;
            for (int k = A.rowptr[i]; k < A.rowptr[i+1]; ++k) {
                if (A.colind[k] == i) { d = A.val[k]; break; }
            }
            // 念のための保護（ゼロ対角はCGの前提を満たさない）
            inv_diag[i] = (fabs(d) > 0.0) ? 1.0 / d : 1.0;
        }
    }
    void apply(const vector<double>& r, vector<double>& z) const {
        z.resize(r.size());
        for (size_t i = 0; i < r.size(); ++i) z[i] = inv_diag[i] * r[i];
    }
};

struct CGResult {
    int iters = 0;
    double rel_resid = NAN;
    bool converged = false;
};

// CG法
CGResult conjugate_gradient(
    const CRS& A,
    const vector<double>& b,
    vector<double>& x,           // in: initial guess, out: solution
    int max_iter = 1000,
    double tol = 1e-8,
    bool use_jacobi = true
) {
    const int n = A.n;
    vector<double> r(n), p(n), Ap(n), z(n);

    // r0 = b - A x0
    spmv(A, x, r);
    for (int i = 0; i < n; ++i)
        r[i] = b[i] - r[i];

    double normb = 0.0;
    for (int i = 0; i < n; ++i)
        normb += b[i]*b[i];
    normb = sqrt(normb);
    if (normb == 0.0)
        normb = 1.0; // b=0 でも動くように

    Jacobi M(A);
    if (use_jacobi)
        M.apply(r, z); else z = r;

    // p0 = z0
    p = z;

    // rsold = r^T z
    auto dot = [](const vector<double>& u, const vector<double>& v) {
        long double s = 0.0L;
        for (size_t i = 0; i < u.size(); ++i)
            s += (long double)u[i]*v[i];
        return (double)s;
    };
    double rsold = dot(r, z);

    CGResult res;
    for (int k = 0; k < max_iter; ++k) {
        spmv(A, p, Ap);
        double pAp = dot(p, Ap);
        // SPD なら pAp > 0 のはず（数値誤差で 0 に近いと破綻）
        if (fabs(pAp) < 1e-30) {
            res.iters = k; break;
        }

        double alpha = rsold / pAp;

        // x_{k+1} = x_k + alpha p_k
        for (int i = 0; i < n; ++i)
            x[i] += alpha * p[i];
        // r_{k+1} = r_k - alpha A p_k
        for (int i = 0; i < n; ++i)
            r[i] -= alpha * Ap[i];

        // 収束判定
        double nr = 0.0;
        for (int i = 0; i < n; ++i)
            nr += r[i]*r[i];
        nr = sqrt(nr);
        double rel = nr / normb;
        if (rel < tol) {
            res.iters = k+1;
            res.rel_resid = rel;
            res.converged = true;
            return res;
        }

        // 前処理
        if (use_jacobi)
            M.apply(r, z);
        else
            z = r;

        double rsnew = dot(r, z);
        double beta = rsnew / rsold;

        // p_{k+1} = z_{k+1} + beta p_k
        for (int i = 0; i < n; ++i)
            p[i] = z[i] + beta * p[i];

        rsold = rsnew;
        res.iters = k+1;
        res.rel_resid = rel;
    }
    return res;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " matrix.mtx [tol=1e-10] [max_iter=10000]\n";
        return 1;
    }
    std::string path = argv[1];
    double tol = (argc >= 3) ? std::atof(argv[2]) : 1e-10;    // 収束判定
    int max_iter = (argc >= 4) ? std::atoi(argv[3]) : 10000;  // 最大反復回数

    // 行列データの読み込み
    CRS A = read_mm2crs(path);

    if (A.n <= 0) {
        std::cerr << "Invalid matrix size.\n"; return 1;
    }

    std::vector<double> x(A.n, 0.0);  // 解ベクトル
    std::vector<double> b(A.n, 1.0);  // RHSベクトル（すべて 1）

    double t0 = omp_get_wtime();
    auto out = conjugate_gradient(A, b, x, max_iter, tol, /*Jacobi*/true);
    double t1 = omp_get_wtime();

    std::cout << "cg_time [ms]=" << (t1 - t0) * 1000.0
              << ", converged=" << out.converged
              << ", iters=" << out.iters
              << ", rel_resid=" << out.rel_resid << "\n";

    // 仕上げの残差チェック
    std::vector<double> Ax;
    spmv(A, x, Ax);
    double nr=0.0, nb=0.0;
    for (int i = 0; i < A.n; ++i) {
        double d = Ax[i]-b[i];
        nr += d*d; nb += b[i]*b[i];
    }
    std::cout << "final rel ||Ax-b||/||b|| = " << std::sqrt(nr)/std::sqrt(nb) << "\n";

    return 0;
}
