// SPDX-License-Identifier: MIT OR GPL-3.0-or-later
// Copyright (c) 2026 Andrew Huchala
// Adapted with permission from MathAgents/andrew-huchala,
// commit 708c1efbf00b9675ca3857f92382aaf200063551. See PROVENANCE.md.
#include "harvey_window.hpp"
#include <algorithm>
#include <exception>
#include <array>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <vector>
#ifndef NV
#define NV 4
#endif
using Exp = std::array<int, NV>;
using u32 = uint32_t;
using u64 = uint64_t;
using PowerWord = uint16_t; // p=11,13 only: power_mod <= 13^4 < 2^15
struct Term { Exp e; u32 c; };

static int P, N, D;
static u64 MOD;
static std::vector<Term> terms;

// Number of products of canonical residues that can be added after a reduction.
static u64 product_budget(u64 modulus) {
    const u64 high = modulus - 1;
    return high ? (UINT64_MAX - high) / (high * high) : UINT64_MAX;
}
static std::vector<std::array<u64, NV + 1>> choose;

static u64 binom(int a, int b) { return (a < 0 || b < 0 || b > a) ? 0 : choose[a][b]; }
static u64 powmod(u64 a, u64 e, u64 m) { u64 r = 1 % m; a %= m; while (e) { if (e & 1) r = r * a % m; a = a * a % m; e >>= 1; } return r; }
static u64 invmod(u64 a, u64 m) {           // a unit modulo m = p^N
    int64_t t = 0, nt = 1, r = (int64_t)m, nr = (int64_t)(a % m);
    while (nr) { int64_t q = r / nr; int64_t x = t - q * nt; t = nt; nt = x; x = r - q * nr; r = nr; nr = x; }
    return (u64)((t % (int64_t)m + (int64_t)m) % (int64_t)m);
}
static double now() { return std::chrono::duration<double>(std::chrono::steady_clock::now().time_since_epoch()).count(); }

// ---- monomials: compositions in lexicographic order, as projective_dwork.compositions ----------
static void comp_rec(int total, int parts, int i, Exp &c, std::vector<Exp> &out) {
    if (i == parts - 1) { c[i] = total; out.push_back(c); return; }
    for (int a = 0; a <= total; ++a) { c[i] = a; comp_rec(total - a, parts, i + 1, c, out); }
}
static void compositions(int total, int parts, std::vector<Exp> &out) {
    if (total < 0) return;
    Exp c{};
    comp_rec(total, parts, 0, c, out);
}
// rank of a composition c of `total` into NV parts in that order
static u64 lexrank(const Exp &c, int total) {
    u64 r = 0;
    for (int i = 0; i < NV - 1; ++i) {
        int rest = NV - 1 - i;                       // parts after position i
        r += binom(total + rest, rest) - binom(total - c[i] + rest, rest);
        total -= c[i];
    }
    return r;
}

struct Level {
    int k = 0;
    std::vector<Exp> cols;                // nu (all >= 1), |nu| = k d
    std::vector<int> basis;               // non-pivot columns, increasing
    std::vector<u32> Bas, Low;            // Bas: basis.size() x cols.size();  Low: (cols of level k-1) x cols.size()
    int index(const Exp &nu) const { Exp c; for (int i = 0; i < NV; ++i) c[i] = nu[i] - 1; return (int)lexrank(c, k * D - NV); }
};
static std::vector<Level> levels;         // levels[k], k = 0..K

static long fermat_hilb(int k) {          // #{nu in [1,d-1]^NV : |nu| = k d}
    int t = k * D - NV;
    if (t < 0) return 0;
    std::vector<long> ways(t + 1, 0); ways[0] = 1;
    for (int v = 0; v < NV; ++v) {
        std::vector<long> nw(t + 1, 0);
        for (int s = 0; s <= t; ++s) if (ways[s]) for (int a = 0; a <= D - 2 && s + a <= t; ++a) nw[s + a] += ways[s];
        ways.swap(nw);
    }
    return ways[t];
}

// inverse of a square matrix modulo MOD (invertible modulo p), Gauss-Jordan with unit pivots.
// Reduce all accumulated entries only when the number of pivot updates exhausts a word.
static bool invert(const std::vector<u32> &A32, int n, std::vector<u32> &Xout) {
    if (n >= (1 << 14)) return false;
    std::vector<u64> A((size_t)n * n), X((size_t)n * n, 0);
    for (size_t i = 0; i < A.size(); ++i) A[i] = A32[i];
    for (int i = 0; i < n; ++i) X[(size_t)i * n + i] = 1;
    std::vector<u32> ac(n), xc(n);
    const u64 budget = product_budget(MOD);
    for (int c = 0; c < n; ++c) {
        if (c && (u64)c % budget == 0) {
            for (auto &v : A) v %= MOD;
            for (auto &v : X) v %= MOD;
        }
        int piv = -1;
        for (int r = c; r < n; ++r) if (A[(size_t)r * n + c] % P) { piv = r; break; }
        if (piv < 0) return false;
        if (piv != c) {
            std::swap_ranges(A.begin() + (size_t)piv * n, A.begin() + (size_t)(piv + 1) * n, A.begin() + (size_t)c * n);
            std::swap_ranges(X.begin() + (size_t)piv * n, X.begin() + (size_t)(piv + 1) * n, X.begin() + (size_t)c * n);
        }
        u64 *arow = &A[(size_t)c * n], *xrow = &X[(size_t)c * n];
        u64 inv = invmod(arow[c] % MOD, MOD);
        for (int j = 0; j < n; ++j) { ac[j] = (u32)(arow[j] % MOD * inv % MOD); arow[j] = ac[j]; xc[j] = (u32)(xrow[j] % MOD * inv % MOD); xrow[j] = xc[j]; }
        const u32 *acp = ac.data(), *xcp = xc.data();
        for (int r = 0; r < n; ++r) {
            if (r == c) continue;
            u64 *ar = &A[(size_t)r * n], *xr = &X[(size_t)r * n];
            u64 f = ar[c] % MOD;
            if (!f) { ar[c] = 0; continue; }
            f = MOD - f;
            for (int j = c + 1; j < n; ++j) ar[j] += f * acp[j];
            ar[c] = 0;
            for (int j = 0; j < n; ++j) xr[j] += f * xcp[j];
        }
    }
    Xout.resize((size_t)n * n);
    for (size_t i = 0; i < X.size(); ++i) Xout[i] = (u32)(X[i] % MOD);
    return true;
}

static bool build_level(int k, bool need_relations) {
    Level &L = levels[k];
    L.k = k;
    std::vector<Exp> cs;
    compositions(k * D - NV, NV, cs);
    for (auto &c : cs) { Exp nu; for (int i = 0; i < NV; ++i) nu[i] = c[i] + 1; L.cols.push_back(nu); }
    int nc = (int)L.cols.size();
    const Level &Lw = levels[k - 1];
    int nlow = (int)Lw.cols.size();
    // relations (m, i): m in L^{S\i} of degree k-1
    struct Row { Exp m; int i; };
    std::vector<Row> rows;
    std::vector<Exp> ms;
    compositions((k - 1) * D - (NV - 1), NV, ms);
    for (int i = 0; i < NV; ++i)
        for (auto &c : ms) { Row r; r.i = i; for (int j = 0; j < NV; ++j) r.m[j] = c[j] + (j != i); rows.push_back(r); }
    int nr = (int)rows.size();
    // sparse rows of J: (column, value mod p^N)
    std::vector<std::vector<std::pair<int, u32>>> J(nr);
    for (int r = 0; r < nr; ++r)
        for (const Term &t : terms) {
            int ui = t.e[rows[r].i];
            if (!ui) continue;
            Exp nu; for (int j = 0; j < NV; ++j) nu[j] = rows[r].m[j] + t.e[j];
            J[r].push_back({L.index(nu), (u32)((u64)t.c * ui % MOD)});
        }
    // greedy pivot columns modulo p (forward elimination on the unused rows)
    std::vector<u32> M((size_t)nr * nc, 0);
    for (int r = 0; r < nr; ++r) for (auto &e : J[r]) M[(size_t)r * nc + e.first] = ((u64)M[(size_t)r * nc + e.first] + e.second) % P;
    std::vector<char> used(nr, 0);
    std::vector<int> prow, pcol;
    // At one digit p may itself nearly fill a word. Compute just the inverses
    // needed by this block when a p-entry lookup would exceed its dimensions.
    std::vector<u32> invp(P <= nc ? P : 0), mp(nc);
    for (size_t a = 1; a < invp.size(); ++a) invp[a] = (u32)powmod(a, P - 2, P);
    // entries stay below 2^32: each elimination adds less than p^2, and all rows are reduced before that could overflow
    long budget = (long)((0xFFFFFFFFu - P) / ((u64)P * P)), since = 0;
    for (int c = 0; c < nc; ++c) {
        int piv = -1;
        for (int r = 0; r < nr; ++r) if (!used[r] && M[(size_t)r * nc + c] % P) { piv = r; break; }
        if (piv < 0) { L.basis.push_back(c); continue; }
        used[piv] = 1; prow.push_back(piv); pcol.push_back(c);
        if (budget && ++since >= budget) {
            for (int r = 0; r < nr; ++r) if (!used[r]) for (int j = c; j < nc; ++j) M[(size_t)r * nc + j] %= P;
            since = 1;
        }
        u32 *mrow = &M[(size_t)piv * nc];
        u32 pivot = mrow[c] % P;
        u32 inv = invp.empty() ? (u32)powmod(pivot, P - 2, P) : invp[pivot];
        for (int j = c; j < nc; ++j) mp[j] = (u64)(mrow[j] % P) * inv % P;
        const u32 *mpp = mp.data();
        for (int r = 0; r < nr; ++r) {
            if (used[r]) continue;
            u32 *mr = &M[(size_t)r * nc];
            u32 f = mr[c] % P;
            if (!f) continue;
            f = P - f;
            if (budget) {
                for (int j = c; j < nc; ++j) mr[j] += f * mpp[j];
            } else {
                for (int j = c; j < nc; ++j) mr[j] = (mr[j] + (u64)f * mpp[j]) % P;
            }
        }
    }
    std::vector<u32>().swap(M);
    if ((long)L.basis.size() != fermat_hilb(k)) return false;
    if (!need_relations) return true;
    int np = (int)pcol.size(), nb = (int)L.basis.size();
    std::vector<int> colpos(nc, -1), baspos(nc, -1);
    for (int a = 0; a < np; ++a) colpos[pcol[a]] = a;
    for (int b = 0; b < nb; ++b) baspos[L.basis[b]] = b;
    // S = J[prow][pcol], inverted modulo p^N
    std::vector<u32> S((size_t)np * np, 0), Sinv;
    for (int a = 0; a < np; ++a) for (auto &e : J[prow[a]]) if (colpos[e.first] >= 0) {
        u32 &x = S[(size_t)a * np + colpos[e.first]]; x = (u32)(((u64)x + e.second) % MOD); }
    if (np && !invert(S, np, Sinv)) return false;
    std::vector<u32>().swap(S);
    // x_P = Sinv (ell - J[prow][basis] x_V),  ell_a = -(m_a)_{i_a} x^{m_a}
    L.Bas.assign((size_t)nb * nc, 0);
    L.Low.assign((size_t)nlow * nc, 0);
    for (int b = 0; b < nb; ++b) L.Bas[(size_t)b * nc + L.basis[b]] = 1;
    std::vector<int> lowidx(np, -1);
    std::vector<u32> lowval(np, 0);
    for (int a = 0; a < np; ++a) {
        const Row &r = rows[prow[a]];
        int mi = r.m[r.i];
        if (mi) { lowidx[a] = Lw.index(r.m); lowval[a] = (u32)((MOD - (u64)mi % MOD) % MOD); }
    }
    const bool reduce_basis = (u64)np * terms.size() > product_budget(MOD);
    for (int q = 0; q < np; ++q) {                 // column pcol[q] of Bas and Low
        const u32 *si = &Sinv[(size_t)q * np];
        int c = pcol[q];
        std::vector<u64> accb(nb, 0), accl(nlow, 0);
        for (int a = 0; a < np; ++a) {
            u64 s = si[a];
            if (!s) continue;
            for (auto &e : J[prow[a]]) if (baspos[e.first] >= 0) {
                auto &sum = accb[baspos[e.first]];
                sum += s * e.second;
                if (reduce_basis) sum %= MOD;
            }
            if (lowidx[a] >= 0) accl[lowidx[a]] = (accl[lowidx[a]] + s * lowval[a]) % MOD;
        }
        for (int b = 0; b < nb; ++b) L.Bas[(size_t)b * nc + c] = (u32)((MOD - accb[b] % MOD) % MOD);
        for (int m = 0; m < nlow; ++m) L.Low[(size_t)m * nc + c] = (u32)accl[m];
    }
    return true;
}

// At p=11,13 only: half powers have degree <=104 and at most 198485
// coefficients, each <13^4. Every complete dot product fits in uint64.
static inline u64 row_start(int degree, const int *e) {
    u64 r = 0;
    for (int i = 0; i < NV - 2; ++i) { degree -= e[i]; r += choose[degree + NV - 2 - i][NV - 1 - i]; }
    return r;
}

struct Pull {                               // one coefficient row of f^h from f^(h-1), race-free
    const std::vector<PowerWord> *src; std::vector<PowerWord> *dst; int degree; u64 power_mod;      // degree of the destination
    const std::vector<Term> &power_terms;
    void row(int *e, int rem) {             // e[0..NV-3] fixed, rem = degree left for the last two exponents
        static thread_local std::vector<u64> buf;
        buf.assign(rem + 1, 0);
        for (const Term &t : power_terms) {
            int f[NV]; bool ok = true; int left = degree - D;
            for (int l = 0; l < NV - 2; ++l) { f[l] = e[l] - t.e[l]; if (f[l] < 0) { ok = false; break; } left -= f[l]; }
            if (!ok || left < 0) continue;
            const PowerWord *in = src->data() + row_start(degree - D, f);
            u64 *out = buf.data() + t.e[NV - 2];
            u64 c = t.c;
            // source index j has last exponent left - j; the destination index is j + t.e[NV-2]
            for (int j = 0; j <= left; ++j) out[j] += in[j] * c;
        }
        PowerWord *out = dst->data() + row_start(degree, e);
        for (int j = 0; j <= rem; ++j) out[j] = static_cast<PowerWord>(buf[j] % power_mod);
    }
    void rec(int *e, int i, int rem) {
        if (i == NV - 2) { row(e, rem); return; }
        for (e[i] = 0; e[i] <= rem; ++e[i]) rec(e, i + 1, rem - e[i]);
    }
};

struct Box {
    const PowerWord *A, *B; const int *e; int suffix[NV + 1]; u64 total = 0;
    // In a square, the slices x_0=u and x_0=e[0]-u have equal sums.
    // Count the central slice once. The other slices can be doubled together.
    void square(int degree) {
        const int lo = std::max(0, e[0] - degree);
        auto slice = [&](int u) {
            const int ra = degree - u, rb = degree - e[0] + u;
            rec<1>(ra, choose[ra + NV - 2][NV - 1],
                                       choose[rb + NV - 2][NV - 1]);
        };
        for (int u = lo; 2*u < e[0]; ++u) slice(u);
        total *= 2;
        if ((e[0] & 1) == 0) slice(e[0] / 2);
    }
    void coefficient(int degree, bool is_square) {
        if constexpr (NV >= 3) {
            if (is_square) {
                square(degree);
                return;
            }
        }
        rec<0>(degree, 0, 0);
    }
    // A leaf is a contiguous row of each half power, read in opposite directions.
    // Reduce once per coefficient when the entire dot product fits in a word.
    template <int i>
    void rec(int remA, u64 rowA, u64 rowB) {
        const int remB = suffix[i] - remA;
        const int lo = std::max(0, e[i] - remB);
        const int hi = std::min(e[i], remA);
        if constexpr (i == NV - 2) {
            if (lo > hi) return;
            const PowerWord *a = A + rowA + lo;
            const PowerWord *b = B + rowB + (e[i] - lo);
            for (int j = 0; j <= hi - lo; ++j) total += (u64)a[j] * b[-j];
        } else {
            for (int u = lo; u <= hi; ++u) {
                const int ra = remA - u, rb = remB - e[i] + u;
                rec<i + 1>(ra, rowA + choose[ra + NV - 2 - i][NV - 1 - i],
                                               rowB + choose[rb + NV - 2 - i][NV - 1 - i]);
            }
        }
    }
};

int compute(u64 power_mod) {
    int K0 = NV * (D - 1) / D + 1, top = N - 1;
    if (top > K0) { fprintf(stderr, "needs N - 1 <= K0 = %d\n", K0); return 2; }
    int K = std::max(K0 - 1, top);
    int mmax = P * top;
    const bool use_window = P > 13;
    int cmax = use_window ? K0 * D + NV :
        std::max(D * ((mmax + 1) / 2 + 1), K * D) + NV + 2;
    choose.resize(cmax + 1);
    for (int i = 0; i <= cmax; ++i) { choose[i][0] = 1; for (int k = 1; k <= NV && k <= i; ++k) choose[i][k] = choose[i - 1][k - 1] + choose[i - 1][k]; }

    double t0 = now();
    levels.resize(K0 + 1);
    levels[0].k = 0;
    // level K0 has no basis exactly when X is smooth; it is built for that check even if no slab reaches it
    for (int k = 1; k <= K0; ++k) if (!build_level(k, k <= top)) { fprintf(stderr, "Jacobian rank in W-degree %d is not that of a smooth hypersurface\n", k); return 3; }
    double t1 = now();

    // Integral graded complements have dimensions 1, 19, 1.
    // basis, in the order of ProjectiveDwork.V (sorted by (k, nu))
    struct BasisEl { int k, col; };
    std::vector<BasisEl> V;
    std::vector<std::vector<int>> vpos(K + 1);
    for (int k = 1; k <= K; ++k) for (int c : levels[k].basis) { vpos[k].push_back((int)V.size()); V.push_back({k, c}); }
    int b = (int)V.size();

    // requests: for input j and level k <= top, m = p k - d_j, e = p nu - mu_j for all nu of level k
    std::vector<int> need(mmax + 1, 0);
    for (auto &v : V) for (int k = 1; k <= top; ++k) if (P * k - v.k >= 1) need[P * k - v.k] = 1;
    int H = 0;
    std::vector<int> keep(mmax + 2, 0);
    for (int m = 1; m <= mmax; ++m) if (need[m]) { keep[m / 2] = keep[m - m / 2] = 1; H = std::max(H, m - m / 2); }
    std::vector<std::vector<PowerWord>> layers(H + 1);
    if (!use_window) {
        auto power_terms = terms;
        for (auto &t : power_terms) t.c %= power_mod;
        layers[0] = {1};
        std::vector<PowerWord> tmp;
        const std::vector<PowerWord> *prev = &layers[0];
        for (int h = 1; h <= H; ++h) {
            std::vector<PowerWord> next(choose[D * h + NV - 1][NV - 1], 0);
            int degree = D * h;
            for (int e0 = 0; e0 <= degree; ++e0) {
                int e[NV]; e[0] = e0;
                Pull pl{prev, &next, degree, power_mod, power_terms};
                pl.rec(e, 1, degree - e0);
            }
            if (keep[h]) { layers[h].swap(next); prev = &layers[h]; std::vector<PowerWord>().swap(tmp); }
            else { tmp.swap(next); prev = &tmp; }
        }
    }
    double t2 = now();

    // slabs[k]: |cols_k| x b (row-major), accumulated over (k, b) terms
    std::vector<std::vector<u32>> slab(top + 1);
    for (int k = 1; k <= top; ++k) slab[k].assign(levels[k].cols.size() * (size_t)b, 0);
    // shift maps: index in level k+1 of nu + u, for each term
    std::vector<std::vector<std::vector<int>>> shift(top + 1);
    for (int k = 1; k < top; ++k) {
        shift[k].resize(terms.size());
        for (size_t t = 0; t < terms.size(); ++t) {
            shift[k][t].resize(levels[k].cols.size());
            for (size_t c = 0; c < levels[k].cols.size(); ++c) {
                Exp nu; for (int j = 0; j < NV; ++j) nu[j] = levels[k].cols[c][j] + terms[t].e[j];
                shift[k][t][c] = levels[k + 1].index(nu);
            }
        }
    }
    // For fixed k, nu and input pole order dj all targets are p*nu - mu.
    // Use readouts mu+pad at the common center p*nu+pad.
    std::vector<std::vector<u32>> window_values(top+1);
    if (use_window) {
        size_t groups = 0;
        std::vector<Exp> support;
        std::vector<u64> coeffs;
        for (const auto &t : terms) { support.push_back(t.e); coeffs.push_back(t.c); }
        try {
            for (int k=1; k<=top; ++k) {
                const auto &cols=levels[k].cols;
                window_values[k].resize(cols.size()*b);
                for (int dj=1; dj<=K; ++dj) {
                    int m=P*k-dj;
                    if (m<1 || levels[dj].basis.empty() || cols.empty()) continue;
                    int pad=std::max(0,(D-1)*NV+1-dj*D);
                    std::vector<Exp> readouts;
                    for (int col : levels[dj].basis) {
                        Exp slot=levels[dj].cols[col];
                        slot.back()+=pad;
                        readouts.push_back(slot);
                    }
                    harvey::Window<NV> window(D,m,P,power_mod,support,coeffs,readouts);
                    groups += cols.size();
                    std::exception_ptr error;
                    for (size_t c=0; c<cols.size(); ++c) {
                        try {
                            Exp center=cols[c];
                            for (int &x : center) x*=P;
                            center.back()+=pad;
                            auto values=window.extract(center);
                            for (size_t r=0; r<values.size(); ++r)
                                window_values[k][c*b+vpos[dj][r]]=values[r];
                        } catch (...) {
                            error=std::current_exception();
                        }
                    }
                    if (error) std::rethrow_exception(error);
                }
            }
            fprintf(stderr, "window_groups %zu\n", groups);
        } catch (const std::exception &e) {
            fprintf(stderr, "window: %s\n",e.what());
            return 2;
        }
    }
    u64 nreq = 0;
    for (int j = 0; j < b; ++j) {
        const Exp &mu = levels[V[j].k].cols[V[j].col];
        int dj = V[j].k;
        for (int k = 1; k <= top; ++k) {
            int m = P * k - dj;
            if (m < 1 || levels[k].cols.empty()) continue;
            int order = 0; u64 unit = 1;
            for (int x = 2; x <= m; ++x) { int y = x; while (y % P == 0) { y /= P; ++order; } unit = unit * y % MOD; }
            u64 factor = invmod(unit, MOD);
            for (int i = 0; i < k - order; ++i) factor = factor * P % MOD;       // k - order >= 0 for m < p^2
            if (k & 1) factor = (MOD - factor) % MOD;
            size_t nc = levels[k].cols.size();
            std::vector<u64> cur(nc, 0);
            for (size_t c = 0; c < nc; ++c) {
                Exp e;
                for (int l = 0; l < NV; ++l) e[l] = P * levels[k].cols[c][l] - mu[l];
                ++nreq;
                if (*std::min_element(e.begin(), e.end()) < 0) continue;
                if (use_window) {
                    cur[c] = factor * window_values[k][c*b+j] % MOD;
                } else {
                    int a = m/2, bb = m-a;
                    Box box;
                    box.A = layers[a].data(); box.B = layers[bb].data(); box.e = e.data();
                    box.suffix[NV] = 0;
                    for (int l=NV-1; l>=0; --l) box.suffix[l] = box.suffix[l+1]+e[l];
                    box.coefficient(D*a, a==bb);
                    cur[c] = (box.total % MOD) * factor % MOD;
                }
            }
            // multiply by exp(-W f): term bb has W-degree k + bb
            for (int bb = 0; k + bb <= top; ++bb) {
                std::vector<u32> &S = slab[k + bb];
                for (size_t c = 0; c < cur.size(); ++c) if (cur[c]) S[c * b + j] = (u32)((S[c * b + j] + cur[c]) % MOD);
                if (k + bb == top) break;
                std::vector<u64> nxt(levels[k + bb + 1].cols.size(), 0);
                for (size_t t = 0; t < terms.size(); ++t) {
                    const std::vector<int> &sh = shift[k + bb][t];
                    u64 a2 = terms[t].c;
                    for (size_t c = 0; c < cur.size(); ++c) if (cur[c]) nxt[sh[c]] = (nxt[sh[c]] + cur[c] * a2) % MOD;
                }
                u64 scale = (MOD - invmod(bb + 1, MOD)) % MOD;
                for (auto &x : nxt) x = x * scale % MOD;
                cur.swap(nxt);
            }
        }
    }
    double t3 = now();

    // reduction: from the top slab down
    std::vector<u32> A((size_t)b * b, 0);                 // A[i * b + j]: coefficient of V[j] in alpha(V[i])
    std::vector<u64> carry;
    for (int k = top; k >= 1; --k) {
        const Level &L = levels[k];
        size_t nc = L.cols.size();
        if (!nc) { carry.clear(); continue; }
        std::vector<u64> x(nc * (size_t)b);
        for (size_t i = 0; i < x.size(); ++i) x[i] = (slab[k][i] + (carry.empty() ? 0 : carry[i])) % MOD;
        size_t nb = L.basis.size(), nlow = levels[k - 1].cols.size();
        for (int j = 0; j < b; ++j)
            for (size_t r = 0; r < nb; ++r) {
                u64 s = 0; const u32 *row = &L.Bas[r * nc];
                for (size_t c = 0; c < nc; ++c) s = (s + row[c] * x[c * b + j]) % MOD;
                A[(size_t)j * b + vpos[k][r]] = (u32)s;
            }
        carry.assign(nlow * (size_t)b, 0);
        for (int j = 0; j < b; ++j)
            for (size_t m = 0; m < nlow; ++m) {
                u64 s = 0; const u32 *row = &L.Low[m * nc];
                for (size_t c = 0; c < nc; ++c) s = (s + row[c] * x[c * b + j]) % MOD;
                carry[m * b + j] = s;
            }
    }
    double t4 = now();
    fprintf(stderr, "backend %s\n", use_window ? "window" : "half_power");

    printf("%d\n", b);
    for (auto &v : V) { printf("%d", v.k); for (int l = 0; l < NV; ++l) printf(" %d", levels[v.k].cols[v.col][l]); printf("\n"); }
    for (int i = 0; i < b; ++i) { for (int j = 0; j < b; ++j) printf(j ? " %u" : "%u", A[(size_t)i * b + j]); printf("\n"); }
    fprintf(stderr, "rank %d levels %.6fs powers %.6fs coefficients+grouping %.6fs reduction %.6fs requests %llu largest_power %d\n",
            b, t1 - t0, t2 - t1, t3 - t2, t4 - t3, (unsigned long long)nreq, mmax);
    return 0;
}

int main() {
    int nt;
    if (!(std::cin >> P >> N >> nt)) return 2;
    // Protocol v1: full dense quartic, bounded prime and absolute precision.
    if (P < 11 || P > 101 || N < 3 || N > (P <= 47 ? 5 : 4) || nt != 35) return 2;
    for (int a = 2; a*a <= P; ++a) if (P % a == 0) return 2;
    MOD = 1;
    for (int i = 0; i < N; ++i) {
        if (MOD > UINT32_MAX / (u64)P) {
            fprintf(stderr, "needs p^N < 2^32\n");
            return 2;
        }
        MOD *= P;
    }
    terms.resize(nt);
    D = 0;
    std::map<Exp, bool> seen;
    for (auto &t : terms) {
        int s = 0; u64 c;
        for (int &x : t.e) { if (!(std::cin >> x) || x < 0 || x > 4) return 2; s += x; }
        if (!(std::cin >> c) || s != 4 || c == 0 || c >= (u64)P || !seen.emplace(t.e, true).second) return 2;
        D = s;
        t.c = (u32)powmod(c % P, MOD / P, MOD);          // Teichmuller lift
    }
    if (D % P == 0) return 2;
    // Every requested power has m = pk - d_j with d_j >= 1 and 0 < m < p^2.
    // Thus v_p(m!) = floor(m/p) <= k-1. The multiplier (-p)^k / m!
    // supplies at least one digit, so a coefficient modulo p^(N-1) suffices.
    const u64 power_mod = MOD / P;
    return compute(power_mod);
}
