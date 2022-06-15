#pragma once
#ifndef MINVERVC_INCLUDED
#define MINVERVC_INCLUDED

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
        size_type count(const value_type& value) const {
            return index[value] != std::numeric_limits<T>::max();
        }

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

        size_type operator[](const value_type& value) const { return index[value]; }
    };

} // namespace internal


class minVerVC {
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

    duration_ms cutoff_time;
    long long step;
    int optimal_size;
    int v_num;
    int e_num;

    /*structures about edge*/
    Vec<Edge> edge;
    Vec<int> edge_weight;
    int default_edge_weight = 1;

    /*structures about vertex*/
    Vec<int> dscore;  // dscore of v
    Vec<long long> time_stamp;
    int best_cov_v;
    Vec<int> v_beg_idx;
    Vec<int> v_edges;
    Vec<int>
        v_adj;
    int c_size;
    Vec<char> v_in_c;
    Vec<int> remove_cand;
    Vec<int> index_in_remove_cand;
    int remove_cand_size;

    // best solution found
    int best_c_size;
    Vec<char>
        best_v_in_c;
    duration_ms best_comp_time;
    long best_step;

    Vec<int> uncov_stack;  // store the uncov edge number
    int uncov_stack_fill_pointer;
    Vec<int>
        index_in_uncov_stack;
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
    minVerVC(
        int optimal_size,
        Duration cutoff_time) :
        cutoff_time(std::chrono::duration_cast<duration_ms>(cutoff_time)),
        optimal_size(optimal_size) {
        set_random_seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    }
private:
    void init_internal(int num_vertices, int num_edges);

    inline int degree(int v) const {
        return v_beg_idx[v + 1] - v_beg_idx[v];
    }

    inline void update_best_sol();

    inline void update_best_cov_v();

    template <typename Is>
    void build_instance(Is& str);

    void build_instance(const int& num_vertices,
        const Vec<Edge>& edges);
    void update_instance_internal(const Vec<int>& v_degree);

    inline void reset_remove_cand();
    void update_target_size();

    inline void uncover(int e);
    inline void cover(int e);
    void init_sol();
    void add(int v);
    void add_init(int v, heap_t& v_heap);
    void remove(int v);
    void forget_edge_weights();
    void update_edge_weight();
    void cover_LS();
    void cover_LS(const std::function<bool(const minVerVC&, bool)>& callback_on_update);
    bool check_solution() const;

    template <typename Duration>
    void set_cutoff_time(Duration d) {
        cutoff_time = std::chrono::duration_cast<duration_ms>(d);
    }
    void set_optimal_size(int size) { optimal_size = size; }

    void set_random_seed(unsigned int seed) { mt_rand.seed(seed); }

    std::pair<int, Vec<Edge>> get_instance_as_edgelist() const {
        return std::make_pair(v_num, edge);
    }
    Vec<int> get_cover(bool bias_by_one) const;
    Vec<char> get_cover_as_flaglist() const { return best_v_in_c; }
    Vec<int> get_independent_set(bool bias_by_one) const;

    inline int get_vertex_count() const { return v_num; }
    inline int get_edge_count() const { return e_num; }
    inline int get_best_cover_size() const { return best_c_size; }

    inline long get_best_step() const { return best_step; }

    inline std::chrono::milliseconds get_best_duration() const;

    inline std::chrono::milliseconds get_total_duration() const;

    static bool default_stats_printer(const minVerVC& solver,
        bool better_cover_found);
public:
    bool NuMVC_minVertexCovSolver(Vec<Edge>& edges, int NodeNum, Vec<int>& solution);
};

#endif


