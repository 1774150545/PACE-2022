/************************************************
** This is a local search solver for Minimum Vertex Cover.
************************************************/


#ifndef NuMVC_INCLUDED
#define NuMVC_INCLUDED

#include <algorithm>
#include <chrono>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <utility>
#include <vector>

//#include "internal/indexed_heap.hpp"

#include <functional>
#include <limits>
#include <type_traits>
#include <vector>


namespace internal {

    template <typename T, class Compare = std::less<T>>
    class Indexed_Heap {
        static_assert(std::is_integral<T>::value, "Type has to be integral type");

    private:
        using Container = std::vector<T>;
        using Index = std::vector<T>;

    public:
        using container_type = Container;
        using value_compare = Compare;
        using value_type = typename Container::value_type;
        using size_type = typename Container::size_type;
        using reference = typename Container::reference;
        using const_reference = typename Container::const_reference;

    private:
        Container data;
        Index index;
        Compare comp;

    private:
        void element_swap(const size_type& a, const size_type& b) {
            std::swap(data[a], data[b]);
            index[data[a]] = a;
            index[data[b]] = b;
        }

        bool is_leaf(const size_type& pos) const {
            if ((pos >= data.size() / 2) && (pos < data.size())) return true;
            return false;
        }

        size_type left_child(const size_type& pos) const { return 2 * pos + 1; }

        size_type right_child(const size_type& pos) const { return 2 * pos + 2; }

        size_type parent(const size_type& pos) const { return (pos - 1) / 2; }

        void shift_down(size_type pos) {
            while (!is_leaf(pos)) {
                // std::cout << "shift _ down pos :" << pos << std::endl;
                auto min = pos;
                auto lc = left_child(pos);
                auto rc = right_child(pos);
                if ((lc < data.size()) && (!comp(data[lc], data[min]))) min = lc;
                if ((rc < data.size()) && (!comp(data[rc], data[min]))) min = rc;
                if (min == pos) return;
                element_swap(pos, min);
                pos = min;
            }
        }

    public:
        explicit Indexed_Heap(const Compare& compare = Compare()) : comp(compare) {}

        explicit Indexed_Heap(const size_type& size,
            const Compare& compare = Compare())
            : index(size, std::numeric_limits<T>::max()), comp(compare) {
            data.reserve(size);
        }

    public:
        void resize(const size_type& size) {
            data.reserve(size);
            index.resize(size);
        }

        void clear() {
            data.clear();
            std::fill(index.begin(), index.end(), std::numeric_limits<T>::max());
        }

        const_reference top() const { return data[0]; }

        bool empty() const { return data.size() == 0; }

        size_type size() const { return data.size(); }

        void push(const value_type& value) {
            int curr = data.size();
            index[value] = curr;
            data.push_back(value);

            while (curr != 0 && !comp(data[curr], data[parent(curr)])) {
                element_swap(curr, parent(curr));
                curr = parent(curr);
            }
        }

        void pop() {
            index[0] = std::numeric_limits<T>::max();
            std::swap(data[0], data.back());
            data.pop_back();
            index[data[0]] = 0;
            if (data.size() != 0) shift_down(0);
        }

        /**
         * returns 1 if element is in heap, 0 otherwise
         */
        size_type count(const value_type& value) const {
            return index[value] != std::numeric_limits<T>::max();
        }

        /**
         * remove element at given position and restore heap
         */
        void erase(const size_type& pos) {
            if (pos == (data.size() - 1)) {
                index[data.back()] = std::numeric_limits<T>::max();
                data.pop_back();
            }
            else {
                element_swap(pos, data.size() - 1);
                index[data.back()] = std::numeric_limits<T>::max();
                data.pop_back();
                auto idx = pos;
                while ((idx != 0) && (!comp(data[idx], data[parent(idx)]))) {
                    element_swap(idx, parent(idx));
                    idx = parent(idx);
                }
                if (data.size() != 0) shift_down(idx);
            }
        }

        /**
         * return position in heap of given value
         */
        size_type operator[](const value_type& value) const { return index[value]; }
    };

} // namespace internal


class NuMVC {
public:
    // using Edge = std::pair<int, int>;
    struct Edge {
        int first;
        int second;
    };
private:
    using timepoint_t = std::chrono::time_point<std::chrono::system_clock>;
    using duration_ms = std::chrono::milliseconds;
    template<typename T>using Vec = std::vector<T>;
    using compare_t = std::function<bool(const int&, const int&)>;
    using heap_t = internal::Indexed_Heap<int, compare_t>;
    template<typename Value1, typename Value2>void Cout(Value1& value1, Value2& value2) {
#if COUTDATE
        std::cout << value1 << value2;
#endif//COUTDATE
    }
    template<typename T>void Cout(T& value1) {
#if COUTDATE
        std::cout << value1;
#endif //COUTDATE
    }
    void Cout() {
#if COUTDATE
        std::cout << std::endl;
#endif //COUTDATE
    }
private:
    static constexpr int try_step = 100000;

    bool verbose = false;

    timepoint_t start, finish;

    /*parameters of algorithm*/
    duration_ms cutoff_time;  // time limit
    long long step;
    int optimal_size;  // terminate the algorithm before step limit if it finds a
                       // vertex cover of optimal_size

    /*parameters of the instance*/
    int v_num;  //|V|: 0...v-1
    int e_num;  //|E|: 0...e-1

    /*structures about edge*/
    Vec<Edge> edge;
    Vec<int> edge_weight;
    int default_edge_weight = 1;

    /*structures about vertex*/
    Vec<int> dscore;  // dscore of v
    Vec<long long> time_stamp;
    int best_cov_v;

    // from vertex to it's edges and neighbors
    Vec<int> v_beg_idx;  // v_edges and v_adj is flattened 2-d array,
                                 // hence store indices
    Vec<int> v_edges;    // edges related to v, v_edges[i][k] means vertex
                                 // v_i's k_th edge
    Vec<int>
        v_adj;  // v_adj[v_i][k] = v_j (actually, that is v_i's k_th neighbor)

    /* structures about solution */
    // current candidate solution
    int c_size;                    // cardinality of C
    Vec<char> v_in_c;      // a flag indicates whether a vertex is in C
    Vec<int> remove_cand;  // remove candidates, an array consists of only
                                   // vertices in C, not including tabu_remove
    Vec<int> index_in_remove_cand;
    int remove_cand_size;

    // best solution found
    int best_c_size;
    Vec<char>
        best_v_in_c;  // a flag indicates whether a vertex is in best solution

    duration_ms best_comp_time;
    long best_step;

    // uncovered edge stack
    Vec<int> uncov_stack;  // store the uncov edge number
    int uncov_stack_fill_pointer;
    Vec<int>
        index_in_uncov_stack;  // which position is an edge in the uncov_stack

    // CC and taboo, pack densly to avoid cache misses
    Vec<bool> conf_change;
    int tabu_remove = 0;

    // smooth
    static constexpr float p_scale = 0.3;  // w=w*p
    int delta_total_weight = 0;
    int ave_weight = 1;

    // random
    std::mt19937 mt_rand;

public:
    template <typename Duration>
    NuMVC(
        /// Graph in edge-list format
        //const Vec<std::pair<int, int>> &edges,
        /// Number of vertices in graph
        //const int &num_vertices,
        /// Size of optimal vertex cover (set to 0 if not known)
        int optimal_size,
        /// Stop calculation after this duration (chrono duration)
        Duration cutoff_time,
        /// Print messages during calculation
        bool verbose = false,
        /// seed for random number generator
        unsigned int rnd_seed =
        std::chrono::high_resolution_clock::now().time_since_epoch().count())
        : verbose(verbose),
        cutoff_time(std::chrono::duration_cast<duration_ms>(cutoff_time)),
        optimal_size(optimal_size) {
        set_random_seed(rnd_seed);
        //build_instance(num_vertices, edges);
        //init_sol();
    }

private:
    void init_internal(int num_vertices, int num_edges);

    inline int degree(int v) const {
        return v_beg_idx[v + 1] - v_beg_idx[v];
    }

    inline void update_best_sol() {
        best_v_in_c = v_in_c;
        best_c_size = c_size;
        finish = std::chrono::system_clock::now();
        best_comp_time = std::chrono::duration_cast<duration_ms>(finish - start);
        best_step = step;
    }

    /** choose a vertex u in C with the highest dscore,
     * breaking ties in favor of the oldest one;
     */
    inline void update_best_cov_v() {
        best_cov_v = remove_cand[0];
        auto dscore_best = dscore[best_cov_v];

        for (int i = 0; i < remove_cand_size - 1; ++i) {
            int v = remove_cand[i];
            if (dscore[v] < dscore_best) continue;
            if (dscore[v] > dscore_best) {
                if (v != tabu_remove) {
                    best_cov_v = v;
                    dscore_best = dscore[v];
                }
            }
            else if (time_stamp[v] < time_stamp[best_cov_v]) {
                if (v != tabu_remove) {
                    best_cov_v = v;
                    dscore_best = dscore[v];
                }
            }
        }
    }

    template <typename Is>
    void build_instance(Is& str) {
        char line[1024];
        char tempstr1[10];
        char tempstr2[10];
        int e;

        char tmp;
        int v1, v2;

        /*** build problem data structures of the instance ***/
        str.getline(line, 1024);
        while (line[0] != 'p') str.getline(line, 1024);
        sscanf(line, "%s %s %d %d", tempstr1, tempstr2, &v_num, &e_num);
        init_internal(v_num, e_num);

        Vec<int> v_degree(v_num);
        for (e = 0; e < e_num; ++e) {
            str >> tmp >> v1 >> v2;
            // start counting at zero
            --v1;
            --v2;
            ++(v_degree[v1]);
            ++(v_degree[v2]);

            edge[e] = { v1, v2 };
        }
        update_instance_internal(v_degree);
    }

    void build_instance(const int& num_vertices,
        const Vec<Edge>& edges) {
        v_num = num_vertices;
        e_num = edges.size();

        init_internal(v_num, e_num);

        Vec<int> v_degree(v_num);
        for (unsigned int e = 0; e < edges.size(); ++e) {
            ++(v_degree[edges[e].first]);
            ++(v_degree[edges[e].second]);
        }
        edge = edges;
        update_instance_internal(v_degree);
    }

    /**
     * builds internal data structures for processing this instance.
     * has to be called after loading a new instance using
     * build_instance
     */
    void update_instance_internal(const Vec<int>& v_degree);

    inline void reset_remove_cand() {
        int j = 0;
        for (int v = 0; v < v_num; ++v) {
            if (v_in_c[v]) {  // && v!=tabu_remove)
                remove_cand[j] = v;
                index_in_remove_cand[v] = j;
                ++j;
            }
            else
                index_in_remove_cand[v] = 0;
        }
        remove_cand_size = j;
    }

    // kick out the worst vertex in current cover
    void update_target_size() {
        --c_size;

        int max_improvement = std::numeric_limits<int>::min();
        int max_vertex = 0;  // vertex with the highest improvement in C

        for (int v = 0; v < v_num; ++v) {
            if (!v_in_c[v]) continue;
            if (dscore[v] > max_improvement) {
                max_improvement = dscore[v];
                max_vertex = v;
            }
        }
        remove(max_vertex);

        reset_remove_cand();
    }

    inline void uncover(int e) {
        index_in_uncov_stack[e] = uncov_stack_fill_pointer;
        uncov_stack[uncov_stack_fill_pointer++] = e;
    }

    inline void cover(int e) {
        // since the edge is satisfied, its position can be reused to store the
        // last_uncov_edge
        int last_uncov_edge = uncov_stack[--uncov_stack_fill_pointer];
        int index = index_in_uncov_stack[e];
        uncov_stack[index] = last_uncov_edge;
        index_in_uncov_stack[last_uncov_edge] = index;
    }

    void init_sol();

    // add a vertex to current cover
    void add(int v);

    void add_init(int v, heap_t& v_heap);

    void remove(int v);

    void forget_edge_weights();

    void update_edge_weight();

    /**
     * calculate minimum vertex cover
     */
    void cover_LS() { cover_LS(nullptr); }

    /**
     * calculate minimum vertex cover and call callback after
     * each iteration. If callback returns true, stop calculation.
     */
    void cover_LS(
        const std::function<bool(const NuMVC&, bool)>& callback_on_update);
    /**
     * Check if the calculated solution is a valid vertex cover
     */
    bool check_solution() const {
        int e;
        for (e = 0; e < e_num; ++e) {
            if (!(best_v_in_c[edge[e].first] || best_v_in_c[edge[e].second])) {
                //std::cout << "uncovered edge " << e << std::endl;
                /*Cout("uncovered edge ", e);
                Cout();*/
                return false;
            }
        }

        return true;
    }

    /**
     * set the maximum duration after which the solver terminates
     *
     * \param d any \c std::chrono::duration
     */
    template <typename Duration>
    void set_cutoff_time(Duration d) {
        cutoff_time = std::chrono::duration_cast<duration_ms>(d);
    }

    /**
     * Set the size of the vertex cover at which the algorithm should stop
     */
    void set_optimal_size(int size) { optimal_size = size; }

    /**
     * set the seed used for the pseudo random number generator
     */
    void set_random_seed(unsigned int seed) { mt_rand.seed(seed); }

    /**
     * returns the current instance as a pair consisting of the
     * number of vertices and a vector of edges
     */
    std::pair<int, Vec<Edge>> get_instance_as_edgelist() const {
        return std::make_pair(v_num, edge);
    }

    /**
     * return vertex indices of current best vertex cover
     */
    Vec<int> get_cover(bool bias_by_one = true) const {
        Vec<int> cover;
        for (int i = 0; i < v_num; i++) {
            if (best_v_in_c[i]) {
                cover.push_back(i + bias_by_one);
            }
        }
        return cover;
    }

    /**
     * returns a vector of flags, where a true-flag at position i denots
     * that vertex i is covered
     */
    Vec<char> get_cover_as_flaglist() const { return best_v_in_c; }

    /**
     * return vertex indices of current best independent set
     */
    Vec<int> get_independent_set(bool bias_by_one = true) const {
        Vec<int> iset;
        for (int i = 0; i < v_num; i++) {
            if (!best_v_in_c[i]) {
                iset.push_back(i + bias_by_one);
            }
        }
        return iset;
    }

    /**
     * Number of vertices
     */
    inline int get_vertex_count() const { return v_num; }

    /**
     * Number of edges
     */
    inline int get_edge_count() const { return e_num; }

    /**
     * Size of current best vertex cover
     */
    inline int get_best_cover_size() const { return best_c_size; }

    /**
     * Tries necessary for current best vertex cover
     */
    inline long get_best_step() const { return best_step; }

    /**
     * duration for calculating current best vertex cover
     */
    inline std::chrono::milliseconds get_best_duration() const {
        return std::chrono::duration_cast<std::chrono::milliseconds>(
            best_comp_time);
    }

    /**
     * total duration since start of calculation
     */
    inline std::chrono::milliseconds get_total_duration() const {
        return std::chrono::duration_cast<std::chrono::milliseconds>(
            (std::chrono::system_clock::now() - start));
    }

    /**
     * Print statistics during calculation
     */
    static bool default_stats_printer(const NuMVC& solver,
        bool better_cover_found) {
        if (better_cover_found) {
            auto time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                solver.get_best_duration());
            /*std::cout << "Better MVC found.\tSize: " << solver.get_best_cover_size()
                << "\tTime: " << std::fixed << std::setw(4)
                << std::setprecision(4) << time_ms.count() << "ms" << std::endl;*/
        }
        return false;
    }
public:
    bool NuMVC_minVertexCovSolver(Vec<Edge>& edges, int NodeNum, Vec<int>& solution) {
        build_instance(NodeNum, edges);
        init_sol();
        cover_LS();
        int Num = 0;
        if (check_solution() == 1) {
            for (int i = 0; i < v_num; ++i) {
                if (best_v_in_c[i]) { solution.push_back(i); Num++; }
            }
            Cout("====*********",Num );
            Cout("c Best found vertex cover size = " ,best_c_size);
            Cout();
            Cout("c SearchSteps for best found vertex cover = ", best_step);
            Cout();
            auto tt = best_comp_time / 1000;
            Cout("c SearchTime for best found vertex cover = ",tt );
            Cout();
            return true;
        }
        return false;
    }
};

#endif
