#pragma once
#ifndef MAXVERVC_INCLUDED
#define MAXVERVC_INCLUDED
#include <algorithm>
#include <chrono>
#include <cstring>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

class maxVerVC
{
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

    static constexpr int try_step = 10; 
    bool verbose = false;

    timepoint_t start, finish;

    duration_ms cutoff_time;  // time limit
    long long step;
    int optimal_size;  
    int v_num;  //|V|: 1...v
    int e_num;  
    Vec<Edge> newEdge;
    Vec<int> dscore;  // dscore of v
    Vec<long long> time_stamp;
    Vec<int>edgeWeight;
    Vec<int> v_beg_idx;  
    Vec<int> v_edges;    
    Vec<int>  v_adj;  
    Vec<int> v_degree;  
    int c_size;
    Vec<char> v_in_c;  
    Vec<int> remove_cand;
    Vec<int> index_in_remove_cand;
    int remove_cand_size;

    int best_c_size;
    Vec<char>
        best_v_in_c;  
    duration_ms best_comp_time;
    long best_step;

   
    Vec<int> uncov_stack;  
    int uncov_stack_fill_pointer;
    Vec<int>
        index_in_uncov_stack;  

    // random
    std::mt19937 mt_rand;
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
public:

    template <typename Duration>
    maxVerVC(
        int optimal_size,
        /// Stop calculation after this duration (chrono duration)
        Duration cutoff_time)
        //bool verbose ,
        //unsigned int rnd_seed)
        :cutoff_time(std::chrono::duration_cast<duration_ms>(cutoff_time)),
        optimal_size(optimal_size) {
        verbose = true;
        set_random_seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
        //build_instance(num_vertices, edges);
        //init_sol();
    }

private:
    void init_internal(int nodeNum, int edgeNum);
    void update_best_sol();

    template <typename Is>
    void build_instance(Is& str);
    void build_instance(const int& nodeNum,
        const Vec<Edge>& edges);
    void update_instance_internal();

    void reset_remove_cand();

    void update_target_size();

    // update the best vertex in C
    int choose_remove_v();

    inline void uncover(int e);

    inline void cover(int e);
    void init_sol();

    void add(int v);

    void remove(int v);
    /**
     * Check if the solution is valid
     */
    bool check_solution();
    /**
     * calculate minimum vertex cover
     */
    inline void cover_LS() { cover_LS(nullptr); }

   
    void cover_LS(const std::function<bool(const maxVerVC&, bool)>& callback_on_update);

   
    template <typename Duration>
    void set_cutoff_time(Duration d) {
        cutoff_time = std::chrono::duration_cast<duration_ms>(d);
    }

    void set_optimal_size(int size) { optimal_size = size; }

   
    void set_random_seed(unsigned int seed) { mt_rand.seed(seed); }

    std::pair<int, Vec<Edge>> get_instance_as_edgelist() const {
        return std::make_pair(v_num, newEdge);
    }

    Vec<int> get_cover(bool bias_by_one) const;
    Vec<char> get_cover_as_flaglist() const; 
    Vec<int> get_independent_set(bool bias_by_one) const;
   
    inline int get_vertex_count() const { return v_num; }
    inline int get_edge_count() const { return e_num; }    
    inline int get_best_cover_size() const { return best_c_size; }
    inline long get_best_step() const { return best_step; }

    inline std::chrono::milliseconds get_best_duration() const;
    /**
     * total duration since start of calculation
     */
    inline std::chrono::milliseconds get_total_duration() const;

    bool default_stats_printer(const maxVerVC& solver,
        bool better_cover_found);
public:
    bool FastVC_minVertexCovSolver(Vec<Edge>& edges, int NodeNum, Vec<int>& solution);
};
#endif
