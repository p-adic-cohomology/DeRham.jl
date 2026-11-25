// SPDX-License-Identifier: MIT OR GPL-3.0-or-later
// Copyright (c) 2026 Andrew Huchala
// Adapted with permission from MathAgents/andrew-huchala,
// commit 708c1efbf00b9675ca3857f92382aaf200063551. See PROVENANCE.md.
// Harvey (2015), Theorem 4.1. Slots t read the coefficient at center - t.
// Every slot has the same degree H > (degree-1)*n. Several output slots
// share the same recurrence, including its diagonal starting polynomial.
#pragma once
#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <map>
#include <stdexcept>
#include <utility>
#include <vector>

namespace harvey {
using Word = uint64_t;
using Wide = __uint128_t;
inline Word multiply(Word a, Word b, Word q) { return Wide(a)*b % q; }
inline Word power(Word a, Word n, Word q) {
    Word r = 1;
    for (; n; n >>= 1, a = multiply(a,a,q))
        if (n & 1) r = multiply(r,a,q);
    return r;
}

template<size_t Dim>
class Window {
public:
    using Exp = std::array<int,Dim>;
private:
    int degree, m;
    Word p, q, scale, inverse_scale_unit;
    std::vector<Exp> slots;
    std::vector<int> next_slot, group, outputs, order;
    std::vector<std::pair<int,int>> row_types;
    std::vector<Word> coeff, inv_factorial;
    std::vector<int> powers_by_type;
    Word factorial;

    template<class Sum>
    void advance(const Exp &center, std::vector<Word> &values) const {
        const size_t width = coeff.size();
        std::vector<Word> next(slots.size());
        std::vector<int64_t> weights(row_types.size()*width);
        for (int j=1; j<=m; ++j) {
            for (size_t g=0; g<row_types.size(); ++g) {
                auto [axis, base] = row_types[g];
                for (size_t a=0; a<width; ++a) {
                    int64_t factor = int64_t(center[axis])-base-int64_t(j)*powers_by_type[g*width+a];
                    weights[g*width+a] = factor*int64_t(coeff[a]);
                }
            }
            for (size_t r=0; r<slots.size(); ++r) {
                const int64_t *w = &weights[group[r]*width];
                const int *src = &next_slot[r*width];
                Sum sum = 0;
                for (size_t a=0; a<width; ++a)
                    sum += Sum(w[a])*int64_t(values[src[a]]);
                int64_t residue = sum % int64_t(q);
                next[r] = residue < 0 ? residue+q : residue;
            }
            values.swap(next);
        }
    }

public:
    Window(int degree_, int exponent, Word prime, Word outmod, const std::vector<Exp> &support,
           const std::vector<Word> &coefficients, const std::vector<Exp> &readouts)
        : degree(degree_), m(exponent), p(prime), q(outmod), scale(1), coeff(coefficients) {
        // Bounded word backend. A caller can fall back to arbitrary precision.
        constexpr Word limit = Word(1)<<40;
        if (!q || q>limit || p<2 || m<0 || !degree || degree%p==0 ||
            coeff.empty() || coeff.size()!=support.size() || readouts.empty())
            throw std::invalid_argument("unsupported Harvey window");
        for (Word a=m/p; a; a/=p) {
            for (Word i=0; i<a; ++i) {
                if (q>limit/p) throw std::overflow_error("Harvey window guard precision exceeds 40 bits");
                q*=p;
                scale*=p;
            }
        }
        for (auto &c : coeff) c%=q;
        // Factorials stripped of their p-part. One inversion initializes
        // all inverse factorials, including multinomials at the edges.
        order.resize(m+1);
        inv_factorial.resize(m+1);
        factorial=1;
        for (int k=1; k<=m; ++k) {
            Word u=k; order[k]=order[k-1];
            while (u%p==0) { u/=p; ++order[k]; }
            factorial=multiply(factorial,u,q);
        }
        inv_factorial[m]=power(factorial,q-q/p-1,q);
        for (int k=m; k>0; --k) {
            Word u=k;
            while (u%p==0) u/=p;
            inv_factorial[k-1]=multiply(inv_factorial[k],u,q);
        }
        inverse_scale_unit=power(multiply(power(degree,m,outmod),factorial,outmod),
                                 outmod-outmod/p-1,outmod);
        std::map<Exp,int> index;
        auto intern = [&](const Exp &t) {
            auto [it, fresh]=index.emplace(t, int(slots.size()));
            if (fresh) slots.push_back(t);
            return it->second;
        };
        for (const auto &t : readouts) outputs.push_back(intern(t));
        std::map<std::pair<int,int>,int> types;
        for (size_t r=0; r<slots.size(); ++r) {
            Exp base=slots[r];
            size_t axis=0;
            while (axis<Dim && base[axis]<degree) ++axis;
            if (axis==Dim) throw std::invalid_argument("window degree too small");
            base[axis]-=degree;
            auto key=std::make_pair(int(axis),base[axis]);
            auto [it,fresh]=types.emplace(key,int(types.size()));
            if (fresh) {
                row_types.push_back(key);
                for (const auto &y : support) powers_by_type.push_back(y[axis]);
            }
            group.push_back(it->second);
            for (const auto &y : support) {
                Exp t=base;
                for (size_t i=0; i<Dim; ++i) t[i]+=y[i];
                next_slot.push_back(intern(t));
            }
        }
    }

    std::vector<Word> extract(const Exp &center) const {
        std::vector<Word> values(slots.size());
        for (size_t r=0; r<slots.size(); ++r) {
            Word u=factorial;
            int v=order[m];
            for (size_t i=0; i<Dim; ++i) {
                int a=center[i]-slots[r][i];
                if (a<0 || a%degree || a/degree>m) { u=0; break; }
                u=multiply(u,inv_factorial[a/degree],q);
                v-=order[a/degree];
            }
            if (u) values[r]=multiply(u,power(p,v,q),q);
        }
        // Bound the signed edge weights and both dot-product accumulators
        // before entering the hot loop. All operands are canonical residues.
        Word span=*std::max_element(center.begin(),center.end());
        for (auto [axis,base] : row_types) span=std::max(span,Word(base)+Word(m)*degree);
        Wide w=Wide(span)*(*std::max_element(coeff.begin(),coeff.end()));
        if (w>INT64_MAX || Wide(coeff.size())*(q-1)>(Wide(-1)>>1)/std::max<Wide>(1,w))
            throw std::overflow_error("Harvey window dot product exceeds native range");
        if (Wide(coeff.size())*(q-1)*w<=INT64_MAX) advance<int64_t>(center,values);
        else advance<__int128_t>(center,values);
        std::vector<Word> result;
        for (int r : outputs) {
            if (values[r]%scale) throw std::logic_error("inexact Harvey window division");
            result.push_back(multiply(values[r]/scale,inverse_scale_unit,q/scale));
        }
        return result;
    }
};
}
