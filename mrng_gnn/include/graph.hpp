#pragma once
#include <iostream>
#include <vector>
#include <map>

// Template class for a Graph using adjacency list representation
template <typename T>
class Graph {
public:
    Graph() {}

    // Add a vertex to the graph
    void add_vertex(T v) {
        if (vertices.find(v) == vertices.end()) {
            vertices[v] = std::vector<T>();
        }
    }

    // Add an edge between two vertices
    void add_edge(T v1, T v2) {
        vertices[v1].push_back(v2);
        // vertices[v2].push_back(v1); // For an undirected graph
    }

    // Display the graph
    void display_graph() {
        for (const auto& pair : vertices) {
            std::cout << "Vertex " << pair.first << " is connected to:";
            for (const T& neighbor : pair.second) {
                std::cout << " " << neighbor;
            }
            std::cout << std::endl;
        }
    }

    size_t size() { return vertices.size(); }
    
    
    std::vector<T> get_connected_to(T v) {
        if (vertices.find(v) != vertices.end()) {
            return vertices[v];
        } else {
            // Handle the case where the vertex is not in the graph
            // std::cerr << "Vertex " << v << " not found in the graph." << std::endl;
            return std::vector<T>();
        }
    }

private:
    std::map<T, std::vector<T>> vertices;
};