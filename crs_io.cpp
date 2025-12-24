//
// ファイル読み込み
//

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <numeric>
#include <string>
#include <algorithm>
#include <cctype>
#include "crs.hpp"

// ---- 内部ユーティリティ ----
static std::string to_lower(std::string s) {
    for (auto& ch : s) ch = std::tolower(static_cast<unsigned char>(ch));
    return s;
}

static CRS coo_to_crs(int nrows, int ncols,
                      std::vector<int>& I, std::vector<int>& J, std::vector<double>& V)
{
    (void)ncols; // 未使用だが将来チェック用途に残す
    const size_t nnz = I.size();
    std::vector<size_t> ord(nnz);
    std::iota(ord.begin(), ord.end(), 0);

    std::sort(ord.begin(), ord.end(), [&](size_t a, size_t b) {
        if (I[a] != I[b]) return I[a] < I[b];
        return J[a] < J[b];
    });

    std::vector<int> rows_u; rows_u.reserve(nnz);
    std::vector<int> cols_u; cols_u.reserve(nnz);
    std::vector<double> vals_u; vals_u.reserve(nnz);

    for (size_t idx = 0; idx < nnz; ) {
        int r = I[ord[idx]];
        int c = J[ord[idx]];
        double s = 0.0;

        do {
            s += V[ord[idx]];
            ++idx;
        } while (idx < nnz && I[ord[idx]] == r && J[ord[idx]] == c);

        rows_u.push_back(r);
        cols_u.push_back(c);
        vals_u.push_back(s);
    }

    CRS A;
    A.n = nrows;
    A.rowptr.assign(nrows+1, 0);
    A.colind.resize(cols_u.size());
    A.val.resize(vals_u.size());

    for (size_t k = 0; k < rows_u.size(); ++k)
        ++A.rowptr[rows_u[k]+1];
    for (int i = 0; i < nrows; ++i)
        A.rowptr[i+1] += A.rowptr[i];

    std::vector<int> ctr = A.rowptr;
    for (size_t k = 0; k < rows_u.size(); ++k) {
        int r = rows_u[k];
        int pos = ctr[r]++;
        A.colind[pos] = cols_u[k];
        A.val[pos] = vals_u[k];
    }
    return A;
}

// ---- 公開関数 ----
CRS read_mm2crs(const std::string& filepath) {
    std::ifstream fin(filepath);
    if (!fin)
        throw std::runtime_error("Failed to open: " + filepath);

    std::string banner, mtx, storage, field, symmetry;
    {
        std::string line;
        if (!std::getline(fin, line))
            throw std::runtime_error("Empty file");
        std::istringstream iss(line);
        if (!(iss >> banner >> mtx >> storage >> field >> symmetry))
            throw std::runtime_error("Invalid MatrixMarket header");
        if (banner != "%%MatrixMarket")
            throw std::runtime_error("Not a MatrixMarket file");
        mtx = to_lower(mtx);
        storage = to_lower(storage);
        field = to_lower(field);
        symmetry = to_lower(symmetry);
        if (mtx != "matrix" || storage != "coordinate")
            throw std::runtime_error("Only 'matrix coordinate' is supported");
        if (!(field == "real" || field == "integer" || field == "pattern"))
            throw std::runtime_error("Unsupported field: " + field);
        if (!(symmetry == "general" || symmetry == "symmetric"))
            throw std::runtime_error("Unsupported symmetry: " + symmetry);
    }

    auto skip_comments = [&](std::istream& in) {
        std::streampos pos;
        std::string line;
        while (true) {
            pos = in.tellg();
            if (!std::getline(in, line)) break;
            if (!line.empty() && line[0] == '%') continue;
            in.seekg(pos);
            break;
        }
    };
    skip_comments(fin);

    int nrows=0, ncols=0, nnz=0;
    {
        std::string line;
        if (!std::getline(fin, line))
            throw std::runtime_error("Missing size line");
        std::istringstream iss(line);
        if (!(iss >> nrows >> ncols >> nnz))
            throw std::runtime_error("Invalid size line");
    }

    std::vector<int> I; I.reserve(nnz * (symmetry=="symmetric" ? 2 : 1));
    std::vector<int> J; J.reserve(I.capacity());
    std::vector<double> V; V.reserve(I.capacity());

    auto push_entry = [&](int i, int j, double v) {
        if (i < 0 || i >= nrows || j < 0 || j >= ncols)
            throw std::runtime_error("Index out of range");
        I.push_back(i); J.push_back(j); V.push_back(v);
    };

    for (int k = 0; k < nnz; ++k) {
        std::string line;
        if (!std::getline(fin, line))
            throw std::runtime_error("Unexpected EOF");
        if (line.empty() || line[0] == '%') {
            --k;
            continue;
        }
        std::istringstream iss(line);
        long long ii, jj;
        if (!(iss >> ii >> jj))
            throw std::runtime_error("Invalid entry line (indices)");
        double vv = 1.0;
        if (field == "real" || field == "integer") {
            if (!(iss >> vv))
                throw std::runtime_error("Invalid entry line (value)");
        } else {
            vv = 1.0;
        }
        int i = static_cast<int>(ii - 1);
        int j = static_cast<int>(jj - 1);
        push_entry(i, j, vv);
        if (symmetry == "symmetric" && i != j) push_entry(j, i, vv);
    }

    // 正方行列でなければエラー
    if (nrows != ncols)
        throw std::runtime_error("CG requires a square matrix (nrows == ncols).");

    return coo_to_crs(nrows, ncols, I, J, V);
}
