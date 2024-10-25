#pragma once

#include <functional>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <random>
#include <chrono>
#include <iostream>
#include <stdexcept>
#include <limits>
#include <map>
#include <set>
#include <assert.h>
#include <list>
#include <cmath>
#include <queue>
#include <algorithm>
#include <string>
#include <functional>
#include <unordered_set>
#include <algorithm>
#include "lsh.hpp"
#include "graph.hpp"
#include "common.hpp"
#include "hypercube.hpp"
#include <fstream>
#include <sstream>
#include <chrono>
#include <thread>
#include <mutex>

class Results {
    public:
        int image_number;
        std::vector<std::pair<int, double>> nn;
        std::vector<double> nn_true;
        std::vector<double> t_aproximates;
        std::vector<double> t_true;
        double maf;
        std::vector<double> aaf;
        std::ofstream outputFile;

        Results(string output_file, int m) : maf(0),outputFile(output_file, std::ios::trunc) {

            if (outputFile.is_open()) {
                assert(true);
            } else assert(false);

            if(m == 1) {
                outputFile << "GNNS Results" << endl;
            } else {
                outputFile << "MRNG Results" << endl;
            }
        };

        ~Results() {
            if(outputFile.is_open()) {
                outputFile.close();
            }
        }

        bool change_maf(double new_maf) {
            if (new_maf > maf) {
                maf = new_maf;
                return true;
            }
            return false;
        }


    void print_query(int N) {

        outputFile << "Query: " << this->image_number << endl;

        for(int i = 0 ; i < N; i++) {
            outputFile << "Nearest neighbor-" << i + 1 <<": "<< nn[i].first << endl;
            outputFile << "distanceApproximate: " << nn[i].second << endl;
            outputFile << "distanceTrue: " << nn_true[i] << endl;
        }

        outputFile << endl;
    } 
        
    void print_avarage_results() {
        outputFile << "tAvarageApproximate: " << accumulate(t_aproximates.begin(), t_aproximates.end(), 0.0) / t_aproximates.size() << endl;
        outputFile << "tAvarageTrue: " << accumulate(t_true.begin(), t_true.end(), 0.0) / t_true.size() << endl;
        outputFile << "Maf: " << maf << endl;
        assert(aaf.size() != 0);
        cout << "Aaf " << accumulate(aaf.begin(),aaf.end(),0.0) / aaf.size() << endl;
        outputFile << endl;
    
    }
};

class Graph_Search_Ann {
    public:
        // enum method {G_LSH, G_Hypercube};

    private:
        // enum method internal_method;
        
        // pathes to files
        std::string output_file;
        std::string data_file_path;
        std::string query_file_path;
        
        std::string input_file;
        std::string configuration_file;
        std::string graph_file_path;

        bool optional;
        int number_of_clusters;
        int number_of_vector_hash_tables;
        int number_of_vector_hash_functions;
        int max_number_M_hypercube;
        int number_of_hypercube_dimensions;
        int number_of_probes;

        Results* to_print;

        // number of local steps
        int T; 

        // Number of extensions in the search
        size_t E;
        
        // Number of closest neighbours in the graph creation
        int k;

        // Number of random Restarts
        int R;

        // Number of closest neighbours in the graph search
        size_t N;

        int m;

        int l;

        LSH* lsh;
        Hypercube* hc;

        // Dataset
        std::vector<std::vector<uint8_t>> training_points;

        std::vector<t_point> points_inside;
        // Graph
        Graph<int> search_graph;
        std::map<int, std::vector<int>> graph;

        // Metric function to be used
        std::function<double(const std::vector<uint8_t>&, const std::vector<uint8_t>&)> metric;


        struct CompareDistance {
            std::function<double(const std::vector<uint8_t>&, const std::vector<uint8_t>&)> metric;
            const std::vector<uint8_t>& query_point;  // Add a reference to the query point

            // Constructor
            CompareDistance(std::function<double(const std::vector<uint8_t>&, const std::vector<uint8_t>&)> m, const std::vector<uint8_t>& qp)
                : metric(m), query_point(qp) {}

            bool operator()(const t_point& point1, const t_point& point2) const {
                // Use the custom metric function to compute distances
                double distance1 = metric(point1.get_point(), query_point);
                double distance2 = metric(point2.get_point(), query_point);

                // Use the distances to the query point for comparison
                return distance1 > distance2;
            }
        };

    
    public:

        class graph_search_input {
            public:
            std::map<std::string, std::string> mapped_values;
            
            graph_search_input(int argc, char** argv) {
                std::vector<std::string> options = {
                    "d", "o", "q", "k", "E", "R", "N", "l", "m", "g"
                };

                Get_long_opts optobj (options, Get_long_opts::options::SINGLE);
                optobj.parse_values(argc,(const char**) argv);

                for(const std::string& opt : options) {
                    if (!optobj.get_vector_from_option(opt.c_str()).empty()) {
                        mapped_values[opt] = optobj.get_vector_from_option(opt.c_str()).at(0);
                    } else {
                        mapped_values[opt] = "";
                    }
                }
                // optobj.print();
            };
        };
        
        void set_inputs(graph_search_input& inputs);

        // Function to process a line and add vertices and edges to the graph
        void processLine(const std::string& line, std::map<int, std::vector<int>>& graph);

        Graph_Search_Ann(graph_search_input& input_file_, std::function<double(const std::vector<uint8_t>&, const std::vector<uint8_t>&)> metric_);

        void search();

        void gnns_search(const std::vector<uint8_t>& query_point);

        void graph_search(const std::vector<uint8_t>& query_point);
        
        void create_knn_graph(const std::string);

        void mrng_construction(const std::string);

        void get_real_points(const std::string, const std:: string, std::vector<std::pair<int, double>>&);

        vector<double> new_bf_nearest_neighbor(const vector<uint8_t>&,const vector<vector<uint8_t>>&,const int);

        void print();

        ~Graph_Search_Ann() {
            delete to_print;
            delete lsh;
        }
};