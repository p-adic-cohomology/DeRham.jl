// SPDX-License-Identifier: MIT OR GPL-3.0-or-later
// Direct enumeration, independent of Dwork matrices; small-prime test utility.
#include <cstdint>
#include <iostream>
#include <vector>
int main() {
    int p,n,d,m;
    if (!(std::cin>>p>>n>>d>>m) || p<2 || p>23 || n<4 || n>6 || d<3 || d>4 || m<1 || m>56) return 2;
    std::vector<std::vector<int>> e(m,std::vector<int>(n));
    std::vector<int> c(m);
    for (int j=0;j<m;++j) {
        for (int &a:e[j]) if (!(std::cin>>a) || a<0 || a>d) return 2;
        if (!(std::cin>>c[j]) || c[j]<0 || c[j]>=p) return 2;
    }
    int powers[23][5];
    for (int x=0;x<p;++x) {
        powers[x][0]=1;
        for (int k=1;k<=d;++k) powers[x][k]=powers[x][k-1]*x%p;
    }
    uint64_t count=0;
    for (int first=0;first<n;++first) {
        uint64_t size=1;
        for (int i=first+1;i<n;++i) size*=p;
        for (uint64_t index=0;index<size;++index) {
            int x[6]={}; x[first]=1;
            uint64_t tail=index;
            for (int i=first+1;i<n;++i) {x[i]=tail%p;tail/=p;}
            int value=0;
            for (int j=0;j<m;++j) {
                int term=c[j];
                for (int i=0;i<n;++i) term=term*powers[x[i]][e[j][i]]%p;
                value+=term;
            }
            count+=value%p==0;
        }
    }
    std::cout<<count<<'\n';
}
