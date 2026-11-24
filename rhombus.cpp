// Exact linear rhombus filter, same rules as rhombus.py, written independently (optional, faster).
// Input: one structure per line "k m a b c a b c ..."; output: 1 (no contradiction) or 0 (unrealisable).
// Exact rationals on 64-bit integers with 128-bit intermediates; exits with status 3 on overflow.
// Build: g++ -O2 -o rhombus rhombus.cpp
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <set>
#include <vector>
using namespace std;
typedef long long ll;
typedef __int128 LL;

struct Q { ll n, d; };
static ll chk(LL x) {
    if (x > (LL)4000000000000000000LL || x < -(LL)4000000000000000000LL) { fprintf(stderr, "overflow\n"); exit(3); }
    return (ll)x;
}
static Q mk(LL n, LL d) {
    if (d < 0) { n = -n; d = -d; }
    LL a = n < 0 ? -n : n, b = d;
    while (b) { LL t = a % b; a = b; b = t; }
    if (a > 1) { n /= a; d /= a; }
    return {chk(n), chk(d)};
}
static Q sub(Q a, Q b) { return mk((LL)a.n * b.d - (LL)b.n * a.d, (LL)a.d * b.d); }
static Q add(Q a, Q b) { return mk((LL)a.n * b.d + (LL)b.n * a.d, (LL)a.d * b.d); }
static Q mul(Q a, Q b) { return mk((LL)a.n * b.n, (LL)a.d * b.d); }
static Q dv(Q a, Q b) { return mk((LL)a.n * b.d, (LL)a.d * b.n); }
static bool operator<(const Q& a, const Q& b) { return a.n != b.n ? a.n < b.n : a.d < b.d; }
static bool operator==(const Q& a, const Q& b) { return a.n == b.n && a.d == b.d; }
typedef vector<Q> V;

static int N, K;
static vector<V> rows;
static vector<int> piv;
static vector<char> isPiv;
static vector<int> rowOf;

static void reduce(V& v) {
    for (size_t r = 0; r < rows.size(); r++) {
        Q c = v[piv[r]];
        if (c.n) for (int i = 0; i < N; i++) v[i] = sub(v[i], mul(c, rows[r][i]));
    }
}
static bool addRel(V v) {
    reduce(v);
    for (int i = 0; i < N; i++) if (v[i].n) {
        Q x = v[i];
        for (int t = 0; t < N; t++) v[t] = dv(v[t], x);
        for (size_t r = 0; r < rows.size(); r++) {
            Q c = rows[r][i];
            if (c.n) for (int t = 0; t < N; t++) rows[r][t] = sub(rows[r][t], mul(c, v[t]));
        }
        rowOf[i] = rows.size(); isPiv[i] = 1;
        rows.push_back(v); piv.push_back(i);
        return true;
    }
    return false;
}
static V nf(int x) {
    V v(N, Q{0, 1});
    if (isPiv[x]) { const V& r = rows[rowOf[x]]; for (int i = 0; i < N; i++) if (i != x) v[i] = Q{-r[i].n, r[i].d}; }
    else v[x] = Q{1, 1};
    return v;
}
static bool sameType(int x, int y) { return (x < K) == (y < K); }

static bool run(const vector<int>& tr) {
    int m = tr.size() / 3;
    N = K + m;
    rows.clear(); piv.clear(); isPiv.assign(N, 0); rowOf.assign(N, -1);
    vector<vector<char>> U(N, vector<char>(N, 0));
    for (int j = 0; j < m; j++) for (int t = 0; t < 3; t++) { int p = tr[3 * j + t]; U[p][K + j] = U[K + j][p] = 1; }
    while (true) {
        vector<V> nfs(N);
        for (int x = 0; x < N; x++) nfs[x] = nf(x);
        for (int x = 0; x < N; x++) for (int y = x + 1; y < N; y++)
            if (sameType(x, y) && nfs[x] == nfs[y]) return false;
        set<V> ev;
        map<V, Q> dirs;
        for (int a = 0; a < N; a++) for (int b = a + 1; b < N; b++) if (U[a][b]) {
            V d(N), nd(N), pd(N);
            int f = -1;
            for (int i = 0; i < N; i++) { d[i] = sub(nfs[a][i], nfs[b][i]); nd[i] = Q{-d[i].n, d[i].d}; if (f < 0 && d[i].n) f = i; }
            if (f < 0) return false;
            for (int i = 0; i < N; i++) pd[i] = dv(d[i], d[f]);
            Q q = d[f]; if (q.n < 0) q.n = -q.n;
            auto it = dirs.find(pd);
            if (it == dirs.end()) dirs[pd] = q; else if (!(it->second == q)) return false;
            ev.insert(d); ev.insert(nd);
        }
        bool changed = false;
        vector<vector<char>> dist(N, vector<char>(N, 0));
        for (int x = 0; x < N; x++) for (int y = x + 1; y < N; y++) {
            V d(N), pd(N, Q{0, 1});
            int f = -1;
            for (int i = 0; i < N; i++) { d[i] = sub(nfs[x][i], nfs[y][i]); if (f < 0 && d[i].n) f = i; }
            if (!U[x][y] && ev.count(d)) { U[x][y] = U[y][x] = 1; changed = true; }
            bool dd = sameType(x, y);
            if (!dd && f >= 0) { for (int i = 0; i < N; i++) pd[i] = dv(d[i], d[f]); dd = dirs.count(pd) > 0; }
            dist[x][y] = dist[y][x] = dd;
        }
        for (int x2 = 0; x2 < N; x2++) for (int x4 = x2 + 1; x4 < N; x4++) {
            if (!dist[x2][x4]) continue;
            vector<int> com;
            for (int c = 0; c < N; c++) if (U[x2][c] && U[x4][c]) com.push_back(c);
            for (size_t i = 0; i < com.size(); i++) for (size_t j = i + 1; j < com.size(); j++) {
                int x1 = com[i], x3 = com[j];
                if (!dist[x1][x3]) continue;
                V g(N, Q{0, 1});
                g[x1] = add(g[x1], Q{1, 1}); g[x3] = add(g[x3], Q{1, 1});
                g[x2] = sub(g[x2], Q{1, 1}); g[x4] = sub(g[x4], Q{1, 1});
                if (addRel(g)) changed = true;
            }
        }
        if (!changed) return true;
    }
}

int main() {
    int m;
    while (scanf("%d %d", &K, &m) == 2) {
        vector<int> tr(3 * m);
        for (int& x : tr) if (scanf("%d", &x) != 1) return 2;
        puts(run(tr) ? "1" : "0");
    }
    return 0;
}
