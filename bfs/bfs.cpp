#include "bfs.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstddef>
#include <omp.h>
#include <vector>

#include "../common/CycleTimer.h"
#include "../common/graph.h"

#define ROOT_NODE_ID 0
#define NOT_VISITED_MARKER -1

void vertex_set_clear(vertex_set* list) {
    list->count = 0;
}

void vertex_set_init(vertex_set* list, int count) {
    list->max_vertices = count;
    list->vertices = (int*)malloc(sizeof(int) * list->max_vertices);
    vertex_set_clear(list);
}

// Take one step of "top-down" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighboring vertices to the
// new_frontier.
void top_down_step(
    Graph g,
    vertex_set* frontier,
    vertex_set* new_frontier,
    int* distances)
{
    #pragma omp parallel
    {
        vertex_set list;
        vertex_set_init(&list, g->num_nodes);
        vertex_set* new_frontier_private = &list;

        // std::vector<int> nf_private;

        #pragma omp for
        for (int i=0; i<frontier->count; i++) {
            // printf("Hello %d\n", i);
            int node = frontier->vertices[i];

            int start_edge = g->outgoing_starts[node];
            int end_edge = (node == g->num_nodes - 1)
                            ? g->num_edges
                            : g->outgoing_starts[node + 1];

            // attempt to add all neighbors to the new frontier
            for (int neighbor=start_edge; neighbor<end_edge; neighbor++) {
                int outgoing = g->outgoing_edges[neighbor];

                if (distances[outgoing] == NOT_VISITED_MARKER) {
                    distances[outgoing] = distances[node] + 1;
                    
                    int index;
                    
                    // #pragma omp critical
                    index = new_frontier_private->count++;

                    new_frontier_private->vertices[index] = outgoing;
                    // nf_private.push_back(outgoing);
                }
            }
            // printf("End %d\n", i);
        }

        int start;
        #pragma omp critical
        {
            start = new_frontier->count;
            new_frontier->count += new_frontier_private->count;
            // new_frontier->count += nf_private.size();
        }

        for(int i=0; i<new_frontier_private->count; i++) {
            new_frontier->vertices[start++] = new_frontier_private->vertices[i];
        }

        // for(int i=0; i<nf_private.size(); i++) {
        //     new_frontier->vertices[start++] = nf_private[i];
        // }

        free(list.vertices);
    }
}

// Implements top-down BFS.
//
// Result of execution is that, for each node in the graph, the
// distance to the root is stored in sol.distances.
void bfs_top_down(Graph graph, solution* sol) {

    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set* frontier = &list1;
    vertex_set* new_frontier = &list2;

    // initialize all nodes to NOT_VISITED

    for (int i=0; i<graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    // setup frontier with the root node
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;

    while (frontier->count != 0) {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(new_frontier);

        top_down_step(graph, frontier, new_frontier, sol->distances);

#ifdef VERBOSE
    double end_time = CycleTimer::currentSeconds();
    printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif

        // swap pointers
        vertex_set* tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
    }
}



// Take one step of "bottom-up" BFS.
void bottom_up_step(
    Graph g,
    vertex_set* frontier,
    vertex_set* new_frontier,
    int* distances)
{
    bool fmap[g->num_nodes] = { false };

    for(int i=0; i<frontier->count; i++) {
        fmap[frontier->vertices[i]] = true;
    }

    #pragma omp parallel
    {
        vertex_set list;
        vertex_set_init(&list, g->num_nodes);
        vertex_set* new_frontier_private = &list;

        // std::vector<int> nf_private;

        #pragma omp for
        for (Vertex v=0; v<g->num_nodes; v++) {
            if(distances[v] == NOT_VISITED_MARKER) {
                int start_edge = g->incoming_starts[v];
                int end_edge = (v == g->num_nodes - 1)
                            ? g->num_edges
                            : g->incoming_starts[v + 1];
                            
                for(int j = start_edge; j<end_edge; j++) {
                    Vertex u = g->incoming_edges[j];
                    if(fmap[u]) {
                        distances[v] = distances[u] + 1;
                        
                        int index;
                    
                        // #pragma omp critical
                        index = new_frontier->count++;
                    
                        new_frontier_private->vertices[new_frontier_private->count++] = v;
                        // nf_private.push_back(v);
                        break;
                    }
                }
            }
        }

        int start;
        #pragma omp critical
        {
            start = new_frontier->count;
            new_frontier->count += new_frontier_private->count;
            // new_frontier->count += nf_private.size();
        }

        for(int i=0; i<new_frontier_private->count; i++) {
            new_frontier->vertices[start++] = new_frontier_private->vertices[i];
        }

        // for(int i=0; i<nf_private.size(); i++) {
        //     new_frontier->vertices[start++] = nf_private[i];
        // }

        free(list.vertices);
    }
}


void bfs_bottom_up(Graph graph, solution* sol)
{
    // You will need to implement the "bottom up" BFS here as
    // described in the handout.
    //
    // As a result of your code's execution, sol.distances should be
    // correctly populated for all nodes in the graph.
    //
    // As was done in the top-down case, you may wish to organize your
    // code by creating subroutine bottom_up_step() that is called in
    // each step of the BFS process.

    
    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set* frontier = &list1;
    vertex_set* new_frontier = &list2;

    // initialize all nodes to NOT_VISITED
    for (int i=0; i<graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    // setup frontier with the root node
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;

    while (frontier->count != 0) {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(new_frontier);

        bottom_up_step(graph, frontier, new_frontier, sol->distances);

#ifdef VERBOSE
    double end_time = CycleTimer::currentSeconds();
    printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif
        // swap pointers
        vertex_set* tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
    }
}

void bfs_hybrid(Graph graph, solution* sol)
{
    // You will need to implement the "hybrid" BFS here as
    // described in the handout.
    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set* frontier = &list1;
    vertex_set* new_frontier = &list2;

    int stage = 0;

    // initialize all nodes to NOT_VISITED
    #pragma omp parallel for
    for (int i=0; i<graph->num_nodes; i++)
        sol->distances[i] = NOT_VISITED_MARKER;

    // setup frontier with the root node
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;

    while (frontier->count != 0) {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(new_frontier);
        
        if(stage == 1) {
            if(frontier->count < (graph->num_nodes/24)) {
                stage = 3;
            }
        }

        if(stage == 0) {
            int mf = 0;
            int mu = 0;

            for(int i=0; i<frontier->count; i++) {
                mf += outgoing_size(graph, frontier->vertices[i]);
            }

            for(int i=0; i<graph->num_nodes; i++) {
                if(sol->distances[i] == NOT_VISITED_MARKER) {
                    mu += incoming_size(graph, i);
                }
            }

            if(mf > (mu/14)) {
                stage = 1;
            }
        }
        
        // alpha = 14, beta = 24.
        if(stage == 1) {
            bottom_up_step(graph, frontier, new_frontier, sol->distances);
        }
        else {
            top_down_step(graph, frontier, new_frontier, sol->distances);
        }

#ifdef VERBOSE
    double end_time = CycleTimer::currentSeconds();
    printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif
        // swap pointers
        vertex_set* tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
    }

}
