// beam_hash_patch.cpp
#include <bits/stdc++.h>
#include <omp.h>
using namespace std;

// ====== cấu hình Zobrist (giữ kiểu cũ) ======
static const int MAX_N = 24;   // sửa theo đề nếu cần
static const int MAX_V = 287;  // giá trị ô phải <= MAX_V

uint64_t ZR[MAX_N][MAX_N][MAX_V + 1];

void init_zobrist(int n, int maxv) {
    mt19937_64 rng(712367821);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int v = 0; v <= maxv; ++v)
                ZR[i][j][v] = rng();
}

uint64_t compute_hash(const vector<vector<int>>& a) {
    int n = (int)a.size();
    int m = (int)a[0].size();
    uint64_t h = 0;
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < m; ++j)
            h ^= ZR[i][j][a[i][j]];
    return h;
}

// ====== Action, apply_action (in-place rotate using tmp) ======
struct Action {
    int y, x, size;
    Action(int y=0, int x=0, int size=0) : y(y), x(x), size(size) {}
};

void apply_action(vector<vector<int>>& a, const Action& ac) {
    int y = ac.y, x = ac.x, k = ac.size;
    // tmp kxk
    vector<int> tmp;
    tmp.reserve(k*k);
    for (int i = 0; i < k; ++i)
        for (int j = 0; j < k; ++j)
            tmp.push_back(a[y + i][x + j]);
    // write rotated
    for (int i = 0; i < k; ++i)
        for (int j = 0; j < k; ++j) {
            int val = tmp[i*k + j];
            int ny = y + j;
            int nx = x + (k - 1 - i);
            a[ny][nx] = val;
        }
}

// ====== Heuristic cho block 2 x X ======
struct BlockEval {
    int pair_cnt = 0;
    int row_runs = 0;
    int left_bias = 0;
    int outer_pairs = 0;
};

// helper: trả về giá trị tại (i,j) *sau* khi áp dụng rotation ac lên cur**,
// nhưng không tạo ma trận mới.
// Nếu (i,j) ngoài vùng xoay, trả cur[i][j].
// Chú ý: ac có origin (y,x) và size k.
inline int value_after_rotate(const vector<vector<int>>& cur, const Action& ac, int i, int j) {
    int y = ac.y, x = ac.x, k = ac.size;
    // kiểm tra (i,j) có trong vùng đích của submatrix xoay không
    if (i < y || i >= y + k || j < x || j >= x + k) {
        return cur[i][j];
    }
    // tính vị trí (oy,ox) trong ma trận *gốc* mà sau khi xoay sẽ đến (i,j)
    // Dựa trên mapping: original (oi,oj) -> new (y+oj, x+k-1-oi)
    // => cho (i,j): oj = i - y; oi = k-1 - (j - x)
    int oj = i - y;
    int oi = k - 1 - (j - x);
    int oy = y + oi;
    int ox = x + oj;
    return cur[oy][ox];
}

BlockEval eval_block_after_rotate(const vector<vector<int>>& cur, const Action& ac, int sy, int sx, int X) {
    // compute BlockEval on the matrix that would result from applying ac to cur
    BlockEval e;
    int n = (int)cur.size();
    int m = (int)cur[0].size();
    int y0 = sy, y1 = sy + 1;

    // 1. pair_cnt trong block 2 x X (sử dụng value_after_rotate)
    for (int y = y0; y <= y1; ++y) {
        for (int x = sx; x < sx + X; ++x) {
            if (x + 1 < sx + X) {
                int v1 = value_after_rotate(cur, ac, y, x);
                int v2 = value_after_rotate(cur, ac, y, x+1);
                if (v1 == v2) e.pair_cnt++;
            }
            if (y + 1 <= y1) {
                int v1 = value_after_rotate(cur, ac, y, x);
                int v2 = value_after_rotate(cur, ac, y+1, x);
                if (v1 == v2) e.pair_cnt++;
            }
        }
    }

    // 2. row_runs trong block
    for (int row = y0; row <= y1; ++row) {
        int curv = value_after_rotate(cur, ac, row, sx);
        int len = 1;
        for (int x = sx + 1; x < sx + X; ++x) {
            int v = value_after_rotate(cur, ac, row, x);
            if (v == curv) len++;
            else {
                if (len >= 2) e.row_runs += (len - 1) * (len - 1);
                curv = v; len = 1;
            }
        }
        if (len >= 2) e.row_runs += (len - 1) * (len - 1);
    }

    // 3. left_bias trong block (sử dụng value_after_rotate)
    for (int x = sx; x < sx + X; ++x) {
        int weight = (sx + X - x);
        int v0 = value_after_rotate(cur, ac, y0, x);
        int v1 = value_after_rotate(cur, ac, y1, x);
        if (v0 == v1) e.left_bias += weight;
    }

    // 4. outer_pairs: cặp kề nhau nằm ngoài block 2xX
    // Need to compare neighbors across entire matrix, but with value_after_rotate
    for (int y = 0; y < n; ++y) {
        for (int x = 0; x < m; ++x) {
            bool in_block =
                (y >= sy && y < sy + 2 &&
                 x >= sx && x < sx + X);

            if (x + 1 < m) {
                bool in_block_r =
                    (y >= sy && y < sy + 2 &&
                     x + 1 >= sx && x + 1 < sx + X);
                if (!(in_block && in_block_r)) {
                    int v0 = value_after_rotate(cur, ac, y, x);
                    int v1 = value_after_rotate(cur, ac, y, x+1);
                    if (v0 == v1) e.outer_pairs++;
                }
            }
            if (y + 1 < n) {
                bool in_block_d =
                    (y + 1 >= sy && y + 1 < sy + 2 &&
                     x >= sx && x < sx + X);
                if (!(in_block && in_block_d)) {
                    int v0 = value_after_rotate(cur, ac, y, x);
                    int v1 = value_after_rotate(cur, ac, y+1, x);
                    if (v0 == v1) e.outer_pairs++;
                }
            }
        }
    }

    return e;
}

// combine
long long combine_score(const BlockEval& e) {
    const long long W_PAIR  = 10000000LL;
    const long long W_RUN   = 100LL;
    const long long W_LEFT  = 200LL;
    const long long W_OUTER = 2000LL;
    return e.pair_cnt * W_PAIR + e.row_runs * W_RUN + e.left_bias * W_LEFT + (long long)e.outer_pairs * W_OUTER;
}

// ====== State & Node structures ======
struct FieldState {
    vector<vector<int>> a; // lưu ma trận hiện tại (giá trị thật, <= MAX_V)
};

unordered_map<uint64_t, shared_ptr<FieldState>> state_pool;

// Node in beam (holds only hash and score; parent map lưu bên ngoài)
struct BeamNode {
    uint64_t h;
    long long score;
    BlockEval be;
};

struct Candidate {
    uint64_t h;
    long long score;
    BlockEval be;
    uint64_t parent_h;
    Action ac;
};

// parent map: child_hash -> (parent_hash, Action)
unordered_map<uint64_t, pair<uint64_t, Action>> parent_map;

// ====== incremental hash after rotate (không thay đổi cur) ======
// Tính hash mới bằng cách xor ra vị trí cũ và xor vào vị trí mới,
// sử dụng giá trị trong cur (không thay đổi cur).
inline uint64_t incremental_hash_after_rotate(const vector<vector<int>>& cur, uint64_t old_h, const Action& ac) {
    int y = ac.y, x = ac.x, k = ac.size;
    // copy submatrix values vào tmp (k*k) để biết old values
    // dùng vector để an toàn
    vector<int> tmp; tmp.reserve(k*k);
    for (int i = 0; i < k; ++i)
        for (int j = 0; j < k; ++j)
            tmp.push_back(cur[y + i][x + j]);

    uint64_t h = old_h;
    // cho mỗi ô (oi,oj) gốc -> nó đi đến (y+oj, x+k-1-oi)
    for (int oi = 0; oi < k; ++oi) {
        for (int oj = 0; oj < k; ++oj) {
            int oldy = y + oi, oldx = x + oj;
            int val = tmp[oi*k + oj];
            int newy = y + oj, newx = x + (k - 1 - oi);

            // bỏ contribution cũ tại (oldy,oldx,val) và thêm contribution tại (newy,newx,val)
            h ^= ZR[oldy][oldx][val];
            h ^= ZR[newy][newx][val];
        }
    }
    return h;
}

// ====== BEAM SEARCH with incremental hash & reduced copy ======
vector<Action> beam_search_build_block(
    const vector<vector<int>>& field,
    int sy, int sx, int X,
    int beam_width = 3000,
    int depth_limit = 50,
    int limit_top = 0,
    int limit_bottom = -1,
    int limit_left = 0,
    int limit_right = -1
) {
    int n = (int)field.size();
    int m = (int)field[0].size();
    if (limit_bottom < 0) limit_bottom = n;
    if (limit_right  < 0) limit_right  = m;

    state_pool.clear();
    parent_map.clear();

    // init Zobrist
    init_zobrist(n, MAX_V);

    // tạo state root
    auto root = make_shared<FieldState>();
    root->a = field;
    uint64_t root_h = compute_hash(field);
    state_pool[root_h] = root;

    BlockEval be0 = eval_block_after_rotate(field, Action(0,0,1), sy, sx, X); // rotate dummy no-op
    long long sc0 = combine_score(be0);
    vector<BeamNode> beam;
    beam.push_back({root_h, sc0, be0});

    int base_pairs = be0.pair_cnt;

    for (int depth = 0; depth < depth_limit; ++depth) {
        vector<Candidate> all_cands;
        all_cands.reserve((size_t)beam.size() * 100); // heuristic reserve

        // parallel expand
        #pragma omp parallel
        {
            vector<Candidate> local;
            local.reserve(1024);

            #pragma omp for schedule(dynamic,1)
            for (int bi = 0; bi < (int)beam.size(); ++bi) {
                uint64_t ph = beam[bi].h;
                auto itst = state_pool.find(ph);
                if (itst == state_pool.end()) continue;
                const vector<vector<int>>& cur = itst->second->a;

                for (int k = 2; k <= n; ++k) {
                    for (int y = limit_top; y + k <= limit_bottom; ++y) {
                        for (int x = limit_left; x + k <= limit_right; ++x) {
                            Action ac(y,x,k);
                            uint64_t nh = incremental_hash_after_rotate(cur, ph, ac);
                            BlockEval be = eval_block_after_rotate(cur, ac, sy, sx, X);
                            long long sc = combine_score(be);

                            local.push_back({nh, sc, be, ph, ac});
                        }
                    }
                }
            }

            #pragma omp critical
            all_cands.insert(all_cands.end(), local.begin(), local.end());
        } // end parallel

        if (all_cands.empty()) break;

        // keep best candidate per hash (dedup)
        sort(all_cands.begin(), all_cands.end(), [](const Candidate& a, const Candidate& b){
            if (a.h != b.h) return a.h < b.h;
            return a.score > b.score;
        });

        vector<Candidate> uniq;
        uniq.reserve(min((size_t)beam_width * 4, all_cands.size()));
        for (size_t i = 0; i < all_cands.size(); ++i) {
            if (i > 0 && all_cands[i].h == all_cands[i-1].h) continue;
            uniq.push_back(all_cands[i]);
        }
        all_cands.clear();

        // sort by score desc and keep top beam_width
        sort(uniq.begin(), uniq.end(), [](const Candidate& a, const Candidate& b){
            return a.score > b.score;
        });
        if ((int)uniq.size() > beam_width) uniq.resize(beam_width);

        // Insert selected candidates into state_pool: create child matrices by copying parent's matrix & apply_action
        vector<BeamNode> next_beam;
        next_beam.reserve(uniq.size());

        for (auto &c : uniq) {
            // If not in pool yet, construct child state from parent
            if (state_pool.find(c.h) == state_pool.end()) {
                auto itp = state_pool.find(c.parent_h);
                if (itp == state_pool.end()) continue; // parent missing (shouldn't)
                auto child = make_shared<FieldState>();
                child->a = itp->second->a; // full copy (only for stored nodes, beam_width count)
                apply_action(child->a, c.ac);
                // sanity: recompute hash to be sure (optional)
                uint64_t check = compute_hash(child->a);
                if (check != c.h) {
                    // Hash mismatch: shouldn't happen unless ZR or mapping wrong
                    // but we still insert using computed check to keep consistency
                    c.h = check;
                }
                state_pool[c.h] = child;
                parent_map[c.h] = {c.parent_h, c.ac};
            } else {
                // already exist: but ensure parent_map set (maybe from previous depth)
                if (!parent_map.count(c.h)) parent_map[c.h] = {c.parent_h, c.ac};
            }
            next_beam.push_back({c.h, c.score, c.be});
        }

        if (next_beam.empty()) break;

        // debug
        cerr << "Depth " << depth+1 << ": best_pairs = " << next_beam[0].be.pair_cnt
             << " / " << X << ", score = " << next_beam[0].score
             << " (beam candidates: " << next_beam.size() << ")\n";

        if (next_beam[0].be.pair_cnt >= X) {
            // reconstruct path from parent_map
            vector<Action> path;
            uint64_t curh = next_beam[0].h;
            while (curh != root_h) {
                auto it = parent_map.find(curh);
                if (it == parent_map.end()) break;
                path.push_back(it->second.second);
                curh = it->second.first;
            }
            reverse(path.begin(), path.end());
            return path;
        }

        // move next_beam to beam
        beam.swap(next_beam);

        // prune state_pool: keep only hashes present in beam and their ancestors (parents chain up to root)
        // keep map of needed hashes: beam hashes + their parent hashes recursively
        unordered_set<uint64_t> keep;
        keep.reserve(beam.size()*2 + 10);
        for (auto &bn : beam) keep.insert(bn.h);

        // find parents up to root (optional depth-limited)
        const int PARENT_DEPTH_KEEP = 3; // keep a few ancestors to avoid rebuilding often
        queue<uint64_t> q;
        for (auto h : keep) q.push(h);
        int depth_iter = 0;
        while (!q.empty() && depth_iter < PARENT_DEPTH_KEEP) {
            int sz = q.size();
            for (int i = 0; i < sz; ++i) {
                uint64_t ch = q.front(); q.pop();
                auto it = parent_map.find(ch);
                if (it == parent_map.end()) continue;
                uint64_t ph = it->second.first;
                if (!keep.count(ph)) {
                    keep.insert(ph);
                    q.push(ph);
                }
            }
            depth_iter++;
        }

        if (state_pool.size() > 500000) {
            unordered_map<uint64_t, shared_ptr<FieldState>> new_pool;
            new_pool.reserve(keep.size()*2 + 10);
            for (auto &k : keep) {
                auto it = state_pool.find(k);
                if (it != state_pool.end()) new_pool.emplace(it->first, it->second);
            }
            state_pool.swap(new_pool);
        }
    } // end depth loop

    // no perfect solution, return best if improved by pair_cnt over base
    if (!state_pool.empty()) {
        // find best beam entry
        // reconstruct path for best beam[0]
        // best is beam[0]
        // but need compare pair_cnt to base
        if (!beam.empty() && beam[0].be.pair_cnt > base_pairs) {
            vector<Action> path;
            uint64_t curh = beam[0].h;
            while (curh != 0 && state_pool.count(curh)) {
                auto it = parent_map.find(curh);
                if (it == parent_map.end()) break;
                path.push_back(it->second.second);
                curh = it->second.first;
            }
            reverse(path.begin(), path.end());
            cerr << "Không đạt 100%, trả về best: " << beam[0].be.pair_cnt << "/" << X << "\n";
            return path;
        }
    }

    return {};
}

// ====== Main ======
int main() {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);
    omp_set_num_threads(omp_get_max_threads());

    int n;
    if (!(cin >> n)) return 0;
    vector<vector<int>> F(n, vector<int>(n));
    int maxv = 0;
    for (int i=0;i<n;i++) for (int j=0;j<n;j++) { cin >> F[i][j]; maxv = max(maxv, F[i][j]); }

    if (maxv > MAX_V) {
        cerr << "Error: value > MAX_V = " << MAX_V << ". Please compress values or increase MAX_V.\n";
        return 0;
    }

    int sy = n/2 - 1, sx = 0;
    int X = n/2;

    int margin = 12;
    int top    = max(0, sy - margin);
    int bottom = min(n, sy + 2 + margin);
    int left   = max(0, sx - margin);
    int right  = min(n, sx + X + margin/2);

    auto result = beam_search_build_block(F, sy, sx, X,
                                          2300, 60,
                                          top, bottom, left, right);

    cout << "Found sequence of actions = " << result.size() << "\n";
    for (auto &ac : result) cout << "(" << ac.y << "," << ac.x << "," << ac.size << ")\n";

    // print final matrix after applying actions (for debug)
    vector<vector<int>> finalF = F;
    for (auto &ac : result) apply_action(finalF, ac);

    cerr << "\nFinal matrix:\n";
    for (auto &row : finalF) {
        for (int v : row) cerr << setw(3) << v << " ";
        cerr << "\n";
    }

    BlockEval bef = eval_block_after_rotate(finalF, Action(0,0,1), sy, sx, X);
    cerr << "Final pairs = " << bef.pair_cnt << " / " << X
         << ", row_runs = " << bef.row_runs
         << ", left_bias = " << bef.left_bias << "\n";

    return 0;
}
