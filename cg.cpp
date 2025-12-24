//
// CRS形式の疎行列に対するJacobi前処理付き共役勾配(CG)法
//
#include <iostream>
#include <vector>
#include <numeric>
#include <cmath>
#include <omp.h>
#include "crs.hpp"
#include "precond.hpp"

using namespace std;

// 内積
static inline double dot(const vector<double>& a, const vector<double>& b)
{
    long double s = 0.0L;
    const size_t n = a.size();
    #pragma omp parallel for reduction(+:s) schedule(static)
    for (size_t i = 0; i < n; ++i)
        s += (long double)a[i] * b[i];
    return (double)s;
}

// 2ノルム
static inline double nrm2(const vector<double>& a)
{
    long double s = 0.0L;
    const size_t n = a.size();
    #pragma omp parallel for reduction(+:s) schedule(static)
    for (size_t i = 0; i < n; ++i)
        s += (long double)a[i] * a[i];
    return (double)sqrt((double)s);
}

// 疎行列ベクトル積
// y = A*x
static void spmv(const CRS& A, const vector<double>& x, vector<double>& y)
{
    int n = A.n;
    // y.assign(n, 0.0);
    if ((int)y.size() != n) y.resize(n);  // 再割当て/ゼロ埋めを避ける

#ifdef USEMKL
    const double alpha = 1.0, beta = 0.0;
    mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, alpha, A.handle, A.descr, x.data(), beta, y.data());
#else
    const int*    rowptr = A.rowptr.data();
    const int*    colind = A.colind.data();
    const double* aval   = A.val.data();
    const double* xp     = x.data();
    double*       yp     = y.data();
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
    {
        double sum = 0.0;
        const int rs = rowptr[i];
        const int re = rowptr[i+1];
        for (int k = rs; k < re; ++k)
            sum += aval[k] * xp[colind[k]];

        yp[i] = sum;
    }
#endif
}

// #ifdef JAC
// // Jacobi preconditioner: z = M^{-1} r with M = diag(A)
// struct Jacobi
// {
//     vector<double> inv_diag;  // 逆対角成分

//     explicit Jacobi(const CRS& A)
//     {
//         inv_diag.assign(A.n, 1.0);
//         for (int i = 0; i < A.n; ++i)
//         {
//             double d = 0.0;
//             for (int k = A.rowptr[i]; k < A.rowptr[i+1]; ++k)
//             {
//                 if (A.colind[k] == i)
//                 {
//                     d = A.val[k];
//                     break;
//                 }
//             }
//             // 念のための保護（ゼロ対角はCGの前提を満たさない）
//             inv_diag[i] = (fabs(d) > 0.0) ? 1.0 / d : 1.0;
//         }
//     }

//     void apply(const vector<double>& r, vector<double>& z) const
//     {
//         z.resize(r.size());
//         for (size_t i = 0; i < r.size(); ++i)
//             z[i] = inv_diag[i] * r[i];
//     }
// };
// #endif

struct CGResult
{
    int iters = 0;
    double rel_resid = NAN;
    bool converged = false;
};

// CG法
template <class Preconditioner>
CGResult conjugate_gradient(
    const CRS& A,
    const vector<double>& b,
    vector<double>& x,           // in: initial guess, out: solution
    const Preconditioner& M,
    int max_iter = 10000,
    double tol = 1e-8)
{
    const int n = A.n;
    vector<double> r(n), p(n), Ap(n), z(n);

    // r0 = b - A x0
    spmv(A, x, r);

    #pragma omp parallel for schedule(static)
    for (int i = 0; i < n; ++i)
        r[i] = b[i] - r[i];

    double normb = nrm2(b);
    if (normb == 0.0)
        normb = 1.0; // b=0 でも動くように

    M.apply(r, z);

    // p0 = z0
    p = z;

    // rsold = r^T z
    double rsold = dot(r, z);

    CGResult res;

    // Main loop
    for (int k = 0; k < max_iter; ++k) {
        spmv(A, p, Ap);
        double pAp = dot(p, Ap);
        // SPD なら pAp > 0 のはず（数値誤差で 0 に近いと破綻）
        if (fabs(pAp) < 1e-30)
        {
            res.iters = k;
            break;
        }

        double alpha = rsold / pAp;

        #pragma omp parallel for schedule(static)
        for (int i = 0; i < n; ++i)
        {
            x[i] += alpha * p[i];
            r[i] -= alpha * Ap[i];
        }

        // 収束判定
        double rel = nrm2(r) / normb;
        if (rel < tol)
        {
            res.iters = k+1;
            res.rel_resid = rel;
            res.converged = true;
            return res;
        }

        M.apply(r, z);

        double rsnew = dot(r, z);
        double beta = rsnew / rsold;

        #pragma omp parallel for schedule(static)
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
        cerr << "Usage: " << argv[0] << " matrix.mtx [max_iter=10000] [tol=1e-08]\n";
        return 1;
    }
    string path = argv[1];
    int max_iter = (argc >= 3) ? atoi(argv[2]) : 10000;  // 最大反復回数
    double tol = (argc >= 4) ? atof(argv[3]) : 1e-08;    // 収束判定

    // 行列データの読み込み
    CRS A = read_mm2crs(path);

    if (A.n <= 0) {
        cerr << "Invalid matrix size.\n"; return 1;
    }

    vector<double> x(A.n, 0.0);  // 解ベクトル
    vector<double> b(A.n, 1.0);  // RHSベクトル（すべて 1）

    ///////////////////////////////////////////////////////////////
    double t0 = omp_get_wtime(); // timer start

    #ifdef JAC
        Jacobi M(A);
    #else
        Identity M(A);
    #endif
    auto out = conjugate_gradient(A, b, x, M, max_iter, tol);

    double t1 = omp_get_wtime(); // timer stop
   ///////////////////////////////////////////////////////////////

    cout << M.name() << " applied.\n";
    cout << "cg_time [ms]=" << (t1 - t0) * 1000.0
         << ", converged=" << out.converged
         << ", iters=" << out.iters
         << ", rel_resid=" << out.rel_resid;

    // 仕上げの残差チェック
    vector<double> Ax;
    spmv(A, x, Ax);
    double nr=0.0, nb=0.0;
    for (int i = 0; i < A.n; ++i) {
        double d = Ax[i]-b[i];
        nr += d*d;
        nb += b[i]*b[i];
    }
    cout << ", ||Ax-b||/||b|| = " << sqrt(nr)/sqrt(nb) << "\n";

    return 0;
}
