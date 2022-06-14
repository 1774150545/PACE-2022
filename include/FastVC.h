#pragma once
#ifndef FASTVC_INCLUDED
#define FASTVC_INCLUDED
/************************************************
** This is a local search solver for Minimum Vertex Cover.
************************************************/

#include <algorithm>
#include <chrono>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>


/**
 * Local search solver for Minimum Vertex Cover
 */
class FastVC {
public:
    //using Edge = std::pair<int, int>;
    struct Edge {
        int first;
        int second;
    };

private:
    template<typename T>using Vec = std::vector<T>;
    using timepoint_t = std::chrono::time_point<std::chrono::system_clock>;
    using duration_ms = std::chrono::milliseconds;

    static constexpr int try_step = 10;  // check time each number of steps

    /* print messages while solving */
    bool verbose = false;

    timepoint_t start, finish;

    /*parameters of algorithm*/
    duration_ms cutoff_time;  // time limit
    long long step;
    int optimal_size;  // terminate the algorithm before step limit if it finds a
                       // vertex cover of optimal_size
    /*parameters of the instance*/
    int v_num;  //|V|: 1...v
    int e_num;  //|E|: 0...e-1
    /*structures about edge*/
    Vec<Edge> newEdge;
    /*structures about vertex*/
    Vec<int> dscore;  // dscore of v
    Vec<long long> time_stamp;
    Vec<int>edgeWeight;
    // from vertex to it's edges and neighbors
    Vec<int> v_beg_idx;  // v_edges and v_adj is flattened 2-d array,
                                 // hence store indices
    Vec<int> v_edges;    // edges related to v, v_edges[i][k] means vertex
                                 // v_i's k_th edge
    Vec<int>
        v_adj;  // v_adj[v_i][k] = v_j(actually, that is v_i's k_th neighbor)
    Vec<int> v_degree;  // amount of edges (neighbors) related to v

    /* structures about solution */
    // current candidate solution
    int c_size;                // cardinality of C
    Vec<char> v_in_c;  // a flag indicates whether a vertex is in C
    // remove candidates, an array consists of only vertices in C, not including
    // tabu_remove  use custom stack for perfomance reasons
    Vec<int> remove_cand;
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

    // random
    std::mt19937 mt_rand;
    template<typename Value1, typename Value2>void Cout(Value1& value1, Value2& value2) {
#if COUTDATE
       // std::cout << value1 << value2;
#endif//COUTDATE
    }
    template<typename T>void Cout(T& value1) {
#if COUTDATE
       // std::cout << value1;
#endif //COUTDATE
    }
    void Cout() {
#if COUTDATE
       // std::cout << std::endl;
#endif //COUTDATE
    }
public:
    template <typename Duration>
    FastVC(
        /// Graph in edge-list format
        //const Vec<std::pair<int, int>>& edges,
        /// Number of vertices in graph
        //const int& num_vertices,
        /// Size of optimal vertex cover (set to 0 if not known)
        int optimal_size,
        /// Stop calculation after this duration (chrono duration)
        Duration cutoff_time,
        /// Print messages during calculation
        bool verbose = true,
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
    void init_internal(int nodeNum, int edgeNum);

    void update_best_sol();

    template <typename Is>
    void build_instance(Is& str);

    void build_instance(const int& nodeNum, const Vec<Edge>& edges);

    /**
     * builds internal data structures for processing this instance.
     * has to be called after loading a new instance using
     * build_instance
     */
    void update_instance_internal();

    void reset_remove_cand();

    void update_target_size();

    // update the best vertex in C
    int choose_remove_v();

    inline void uncover(int e) {
        index_in_uncov_stack[e] = uncov_stack_fill_pointer;
        uncov_stack[uncov_stack_fill_pointer++] = e;
    }

    inline void cover(int e) {
        int index, last_uncov_edge;

        // since the edge is satisfied, its position can be reused to store the
        // last_uncov_edge
        last_uncov_edge = uncov_stack[--uncov_stack_fill_pointer];
        index = index_in_uncov_stack[e];
        uncov_stack[index] = last_uncov_edge;
        index_in_uncov_stack[last_uncov_edge] = index;
    }

    void init_sol();

    void add(int v);

    void remove(int v);

    /**
     * Check if the solution is valid
     */
    bool check_solution() const;

    /**
     * calculate minimum vertex cover
     */
    inline void cover_LS() { cover_LS(nullptr); }

    /**
     * calculate minimum vertex cover and call callback after
     * each iteration. If callback returns true, stop calculation.
     */
    void cover_LS(const std::function<bool(const FastVC&, bool)>& callback_on_update);

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
        return std::make_pair(v_num, newEdge);
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
    static bool default_stats_printer(const FastVC& solver,
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
    bool FastVC_minVertexCovSolver(Vec<Edge>&edges,int NodeNum,Vec<int>&solution);
};

#endif
