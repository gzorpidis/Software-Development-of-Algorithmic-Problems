#include "graph_search.hpp"

int main(int argc, char** argv) {
    
    // LSH::lsh_input input(argc, argv);
    // LSH lsh(input, euclideanDistance,500);
    Graph_Search_Ann::graph_search_input input(argc, argv);
    Graph_Search_Ann g(input, euclideanDistance);
    g.search();
    return 0;
}