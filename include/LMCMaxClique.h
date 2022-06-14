//Copyright <2016> <Chu-Min Li & Hua Jiang>
//
//Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
//
//The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <sys/types.h>
#include <limits.h>
#include <math.h>
#include <assert.h>
#include <iostream>
#include <vector>
class LMCMaxClique {
#define WORD_LENGTH 100
#define TRUE 1
#define FALSE 0
#define NONE -1
#define DELIMITER 0
#define PASSIVE 0
#define ACTIVE 1
#define P_TRUE 2
#define P_FALSE 0
#define NO_REASON -3
#define CONFLICT -1978
#define MAX_NODE 10000
#define max_expand_depth 100000
#define PUSH_OPT(item, stack) stack[stack ## _fill_pointer++] = item
#define ptr(stack) stack ## _fill_pointer
#define is_neibor(i,j) matrice[i][j]

#define CUR_CLQ_SIZE Clique_Stack_fill_pointer
#define CURSOR Cursor_Stack[Cursor_Stack_fill_pointer-1]
#define MIN(a,b) a<=b?a:b
#define BIT_MAP_SIZE 4097

#define SET_EDGE(row,col) ((*(Adj_Matrix + (row)* MATRIX_ROW_WIDTH + ((col) >> 3))) |= (1 << ((col) & 7)))
#define GET_EDGE(row,col) ((*(Adj_Matrix + (row)* MATRIX_ROW_WIDTH + ((col) >> 3))) & (1 << ((col) & 7)))

#define iMatrix(i) (Adj_Matrix+(i)*MATRIX_ROW_WIDTH)
#define Matrix(i,j) ((*((i) + ((j) >> 3))) & (1 << ((j) & 7)))

    int FORMAT = 1, nodeNum, NB_NODE_O=0,  edgeNum, NB_EDGE_O=0,
        MAX_CLQ_SIZE=0, MAX_ISET_SIZE=0, INIT_CLQ_SIZE=0,MATRIX_ROW_WIDTH=0, 
        MAX_VERTEX_NO=0, K_CORE_G = 0;

#define CORE_NO Vertex_UB
    int Max_Degree = 0;
    int Node_Degree[MAX_NODE];
    char Node_State[MAX_NODE];
    int** Node_Neibors;

    int Candidate_Stack_fill_pointer = 0;
    int Candidate_Stack[MAX_NODE * 2];
    int Vertex_UB[MAX_NODE * 2];
    int Clique_Stack_fill_pointer;
    int* Clique_Stack, * MaxCLQ_Stack;
    int Cursor_Stack[max_expand_depth];
    int Cursor_Stack_fill_pointer = 0;

    int* Node_Reason;
    unsigned char* Adj_Matrix;

    int iSET_COUNT = 0;
    int* iSET_Size;
    char* iSET_State;
    char* iSET_Used;
    char* iSET_Tested;
    int* iSET_Index;
    char* iSET_Involved;
    char* Is_Tested;
    int** iSET;

    int* REASON_STACK;
    int REASON_STACK_fill_pointer = 0;
    int* CONFLICT_ISET_STACK;
    int CONFLICT_ISET_STACK_fill_pointer;
    int* ADDED_NODE_iSET;
    int* REDUCED_iSET_STACK = Node_Degree;
    int REDUCED_iSET_STACK_fill_pointer = 0;
    int* PASSIVE_iSET_STACK;
    int PASSIVE_iSET_STACK_fill_pointer = 0;
    int* FIXED_NODE_STACK;
    int FIXED_NODE_STACK_fill_pointer = 0;
    int* UNIT_STACK;
    int UNIT_STACK_fill_pointer = 0;
    int* NEW_UNIT_STACK;
    int NEW_UNIT_STACK_fill_pointer = 0;

    int Rollback_Point;
    int Branching_Point;
    //int Matrix_Reallocation = 0;
    int* Old_Name;
    int* Second_Name;
    int NB_CANDIDATE = 0, FIRST_INDEX;
    int START_MAXSAT_THD = 15;

    int Extra_Node_Stack[100000];

    int Last_Idx = 0;
    int cut_ver = 0, total_cut_ver = 0;
    int cut_inc = 0, total_cut_inc = 0;
    int cut_iset = 0, total_cut_iset = 0;
    int cut_satz = 0, total_cut_satz = 0;
    long long Branches_Nodes[6];
    int LAST_IN;
    float Dynamic_Radio = 0.70;
    int REBUILD_MATRIX = FALSE;
    int CUR_MAX_NODE;
    int Branches[1200];
    int* Init_Adj_List;
    int BLOCK_COUNT = 0;
    int* BLOCK_LIST[100];
    double READ_TIME, INIT_TIME, SEARCH_TIME;
    long long N0_0 = 0, N0_1 = 0, N1_0 = 0, N1_1 = 0, L1 = 0;
    double D0 = 0, D1 = 0, D2 = 0, D3 = 0, Dt = 0, R1 = 0;
    long long SubNB = 0;
    template<typename T>using Vec = std::vector<T>;
public:
    struct NewEdge {
        int src;
        int dst;
    };
protected:
    static int int_cmp(const void* a, const void* b) {
        return *((int*)b) - *((int*)a);
    }
    int is_adjacent(int node1, int node2) {
        int neibor, * neibors;
        neibors = Node_Neibors[node1];
        for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
            if (neibor == node2) {
                return TRUE;
            }
        }
        return FALSE;
    }

    void check_clique() {
        int i, j;
        for (i = 0; i < MAX_CLQ_SIZE; i++) {
            for (j = i + 1; j < MAX_CLQ_SIZE; j++) {
                if (MaxCLQ_Stack[i] > MaxCLQ_Stack[j]) {
                    assert(is_adjacent(MaxCLQ_Stack[i], MaxCLQ_Stack[j]));
                }
                else {
                    assert(is_adjacent(MaxCLQ_Stack[j], MaxCLQ_Stack[i]));
                }
            }
        }
    }

    void check_clique_in_result_file(int* solvers) {
        int i = 0, j, node, clique_size = 0;
        Vec<int> clique;

        clique_size = sizeof(solvers) / sizeof(solvers[0]);
        if (clique_size > 0) {
            for (i = 0; i < clique_size; i++) {
                for (j = i + 1; j < clique_size; j++) {
                    if (is_adjacent(clique[i], clique[j]) == FALSE) {
                        //printf("find non-adjacent vertices: %d %d\n", clique[i],clique[j]);
                        //printf("#FALSE:%s\n", solvers);
                        return;
                    }

                }
            }
            //printf("#TRUE:%s\n", solvers);
        }
        else {
            //printf("#NONE:%s\n", solvers);
        }
    }


    void allcoate_memory_for_adjacency_list(int nb_node, int nb_edge,
        int offset) {
        int i, block_size = 40960000, free_size = 0;
        Init_Adj_List = (int*)malloc((2 * nb_edge + nb_node) * sizeof(int));
        if (Init_Adj_List == NULL) {
            for (i = 1; i <= nodeNum; i++) {
                if (Node_Degree[i - offset] + 1 > free_size) {
                    Node_Neibors[i] = (int*)malloc(block_size * sizeof(int));
                    BLOCK_LIST[BLOCK_COUNT++] = Node_Neibors[i];
                    free_size = block_size - (Node_Degree[i - offset] + 1);
                }
                else {
                    Node_Neibors[i] = Node_Neibors[i - 1]
                        + Node_Degree[i - 1 - offset] + 1;
                    free_size = free_size - (Node_Degree[i - offset] + 1);
                }
            }
        }
        else {
            BLOCK_COUNT = 1;
            BLOCK_LIST[BLOCK_COUNT - 1] = Init_Adj_List;
            Node_Neibors[1] = Init_Adj_List;
            for (i = 2; i <= nodeNum; i++) {
                Node_Neibors[i] = Node_Neibors[i - 1] + Node_Degree[i - 1 - offset]
                    + 1;
            }
        }
    }


    int sort_by_degeneracy_ordering();

    int re_number_adj(int node);

    int re_number(int node);

    int addIntoIsetTomitaBis_adj(int node);

    int addIntoIsetTomitaBis(int node);

    int cut_by_iset_less_vertices();

    int cut_by_iset_last_renumber();

#define assign_node(node, value, reason) {\
	Node_State[node] = value;\
	Node_Reason[node] = reason;\
	PUSH_OPT(node, FIXED_NODE_STACK);\
}
    int fix_newNode_for_iset(int fix_node, int fix_iset);

    int fix_oldNode_for_iset(int fix_node, int fix_iset);

#define fix_node(node,iset) ((node > nodeNum)? fix_newNode_for_iset(node, iset):fix_oldNode_for_iset(node, iset))

    int fix_node_iset(int fix_iset);

    int unit_iset_process();

    int unit_iset_process_used_first();

    void identify_conflict_sets(int iset_idx);

    void enlarge_conflict_sets(int & ADDED_NODE);

    void rollback_context_for_maxsatz(int start_fixed, int start_passive,
        int start_reduced);

    void reset_context_for_maxsatz();

    int further_test_reduced_iset(int start);

    int fix_anyNode_for_iset(int fix_node, int fix_iset) {
        if (fix_node > MAX_VERTEX_NO)
            return fix_newNode_for_iset(fix_node, fix_iset);
        else
            return fix_oldNode_for_iset(fix_node, fix_iset);
    }

    int inc_maxsatz_lookahead_by_fl2();

    int inc_maxsatz_on_last_iset(int &ADDED_NODE);

    int open_new_iset_old(int i);

    int simple_further_test_node(int start);

    int test_node_for_failed_nodes(int node, int iset);

    int test_by_eliminate_failed_nodes();

    int cut_by_inc_maxsat_eliminate_first();

    float compute_subgraph_density(int start);

    int compute_subgraph_degree(int start);

    int BRANCHING_COUNT = 0;
    void allocate_memory_for_maxsat() {
        Node_Reason = (int*)malloc((MAX_VERTEX_NO + 1) * 10 * sizeof(int));
        ADDED_NODE_iSET = (int*)malloc((MAX_VERTEX_NO + 1) * 2 * sizeof(int));
        FIXED_NODE_STACK = (int*)malloc((MAX_VERTEX_NO + 1) * 2 * sizeof(int));

        iSET_State = (char*)malloc((MAX_VERTEX_NO + 1) * sizeof(char));
        iSET_Used = (char*)malloc((MAX_VERTEX_NO + 1) * sizeof(char));
        iSET_Tested = (char*)malloc((MAX_VERTEX_NO + 1) * sizeof(char));
        UNIT_STACK = (int*)malloc((MAX_VERTEX_NO + 1) * sizeof(int));
        NEW_UNIT_STACK = (int*)malloc((MAX_VERTEX_NO + 1) * sizeof(int));
        PASSIVE_iSET_STACK = (int*)malloc((MAX_VERTEX_NO + 1) * sizeof(int));
        iSET_Involved = (char*)malloc((MAX_VERTEX_NO + 1) * sizeof(char));
        CONFLICT_ISET_STACK = (int*)malloc(
            (MAX_VERTEX_NO + 1) * 10 * sizeof(int));
        REASON_STACK = (int*)malloc((MAX_VERTEX_NO + 1) * 10 * sizeof(int));
        Is_Tested = (char*)malloc((MAX_VERTEX_NO + 1) * sizeof(char));
    }
    void store_maximum_clique(int node, int print_info);
    void store_maximum_clique2(int print_info);

    int reduce_subgraph(int start);

    void rebuild_matrix(int start);

    int cut_by_inc_ub();

    int find_3_clique(int node);

    void init_for_search(int using_init_clique);

    void allocate_memory() {
        int i;
        Second_Name = (int*)malloc((MAX_VERTEX_NO + 1) * sizeof(int));
        iSET = (int**)malloc((MAX_VERTEX_NO + 1) * sizeof(int*));
        iSET[0] = (int*)malloc(
            (MAX_VERTEX_NO + 1) * (MAX_VERTEX_NO + 1) * sizeof(int));
        for (i = 1; i < MAX_VERTEX_NO; i++) {
            iSET[i] = iSET[i - 1] + MAX_VERTEX_NO + 1;
        }
        iSET_Size = (int*)malloc((MAX_VERTEX_NO + 1) * sizeof(int));
        iSET_Index = (int*)malloc((nodeNum + 1) * sizeof(int));

        if (INIT_CLQ_SIZE >= START_MAXSAT_THD)
            allocate_memory_for_maxsat();
    }

    bool search_maxclique(int cutoff, int using_init_clique, double MaxTime);

    void printMaxClique(Vec<int>& solution) {
        int i;
        //printf("M ");
        if (INIT_CLQ_SIZE < MAX_CLQ_SIZE) {
            for (i = 0; i < MAX_CLQ_SIZE; i++) {
                //printf("%d ", Old_Name[MaxCLQ_Stack[i]]);
                solution.push_back(Old_Name[MaxCLQ_Stack[i]]);
            }

        }
        else {
            for (i = 0; i < MAX_CLQ_SIZE; i++) {
                //printf("%d ", MaxCLQ_Stack[i]);
                solution.push_back(MaxCLQ_Stack[i]);
            }

        }
        //printf("\n");
    }
    void build_init_matrix() {
        int node, neibor, * neibors;
        MATRIX_ROW_WIDTH = MAX_VERTEX_NO / 8 + 1;
        Adj_Matrix = (unsigned char*)malloc(
            (MAX_VERTEX_NO + 1) * MATRIX_ROW_WIDTH * sizeof(char));

        memset(Adj_Matrix, 0,
            (MAX_VERTEX_NO + 1) * MATRIX_ROW_WIDTH * sizeof(char));

        for (node = 1; node <= MAX_VERTEX_NO; node++) {
            Second_Name[node] = node;
            neibors = Node_Neibors[node];
            for (neibor = *neibors; neibor != NONE; neibor = *(++neibors)) {
                SET_EDGE(node, neibor);
                SET_EDGE(neibor, node);
            }
        }
    }
    int* Adj_List;
#define New_Name Node_Degree
    int search_in_2_k_core_graph();

    void free_block() {
        int i = 0;
        for (i = 0; i < BLOCK_COUNT; i++)
            free(BLOCK_LIST[i]);
    }

    void reduce_instance();

    int initialize(int ordering) {
        clock_t state = clock();
        int r = sort_by_degeneracy_ordering();
        if (r == FALSE) {
            //printf("*****************\n");
            reduce_instance();
            
            if (K_CORE_G <= 2) {
                r = search_in_2_k_core_graph();
             
            }
        }
        INIT_TIME = (clock() - state) * 1.0 / CLOCKS_PER_SEC;

        //printf("I the initial time is %4.2lf \n", INIT_TIME);
        // if(r==0) //printf("init!");
        return r;
    }

    void print_branches() {
        int i = 0;
        for (i = 0; i <= nodeNum; i++) {
            //printf("%3d", Branches[i]);
        }
        //printf("\n");
    }
    void print_version() {
        //printf("# Hello! I am LMC with Incremental MaxSAT Reasoning and 2 Level Cores Decomposition. build at %s %s.\n", __TIME__, __DATE__);
        return;
    }

    char* getInstanceName(char* s) {
        if (strrchr(s, '/') == NULL)
            return s;
        else
            return strrchr(s, '/') + 1;
    }
    bool InstanceLoad(Vec<NewEdge>& newEdges, int newNodeNum);

    bool loadOriginresult(Vec<NewEdge>& newEdges, int format);

    void check_result(Vec<NewEdge>& newEdges, int* solvers) {

        if (FORMAT == 1) {
            loadOriginresult(newEdges, 1);
        }
        else if (FORMAT == 2) {
            loadOriginresult(newEdges, 2);
        }
        else {
            loadOriginresult(newEdges, 1);
        }
        /*//printf("R Instance Information: #node=%d #edge=%d density=%9.8f\n", nodeNum,
            edgeNum, ((float)edgeNum * 2 / nodeNum / (nodeNum - 1)));*/
        check_clique_in_result_file(solvers);

    }

    void minVertexCovering(Vec<bool>& outputSolver, Vec<int>& maxCliqueSolution) {
       // printf("%d:\n", MAX_CLQ_SIZE);
        if (INIT_CLQ_SIZE < MAX_CLQ_SIZE) {
            for (int i = 0; i < MAX_CLQ_SIZE; i++) {
                outputSolver[Old_Name[MaxCLQ_Stack[i]]] = 1;
                //printf("%d ", Old_Name[MaxCLQ_Stack[i]]);
                maxCliqueSolution.push_back(Old_Name[MaxCLQ_Stack[i]]);
            }
        }
        else {
            for (int i = 0; i < MAX_CLQ_SIZE; i++) {
                outputSolver[MaxCLQ_Stack[i]] = 1;
                //printf("%d ", MaxCLQ_Stack[i]);
                maxCliqueSolution.push_back(MaxCLQ_Stack[i]);
            }
        }
    }

    bool Checker(Vec<int>& solution, Vec<NewEdge>& newEdges, int newNodeNum) {
        InstanceLoad(newEdges, newNodeNum);
        int clique_size = solution.size();
        int i = 0;
        if (clique_size > 0) {
            for (i = 0; i < clique_size; i++) {
                for (int j = i + 1; j < clique_size; j++) {
                    if (is_adjacent(solution[i], solution[j]) == false) {
                        return false;
                    }
                }
            }
            ////printf("#true\n");
        }
        else {
            ////printf("#NONE\n");
        }
        return true;
    }

public:
    bool maxCliqueSolver(Vec<NewEdge>& newEdges, int newNodeNum, Vec<bool>& Solver, double maxTime);
};