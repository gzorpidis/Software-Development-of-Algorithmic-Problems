// #include "headers.hpp"
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
#include "common.hpp"
#include "structures.hpp"
#include <functional>
#include <unordered_set>
#include <algorithm>

using namespace std;

    class Hypercube {
        private:

            // Some pre-defined metrics to be chosen by Metric_functions
            static double l1(const std::vector<uint8_t>& vector1, const std::vector<uint8_t>& vector2);
            static double l2(const vector<uint8_t>& vector1, const vector<uint8_t>& vector2);
            static double empty([[maybe_unused]] const vector<uint8_t>& vector1, [[maybe_unused]] const vector<uint8_t>& vector2) { 
                return 0;
            };

            int W;
            
            /* The dimension k(=d') given by the user */
            int k;

            /* The maximum number of points that will be checked */
            uint M;

            // Number of points to be checked during the search
            int N;

            // Number of vertices of the hypercube on which we will search on
            int R;

            // Number of vertices of the hypercube on which we will search on
            size_t probes;

            // The metric function to be used in the hypercube
            std::function<double(const std::vector<uint8_t>&, const std::vector<uint8_t>&)> metric;

            /* The training set copy of the hypercube */
            vector<vector<uint8_t>> training_set;

            /* A vector of h functions h_i, dimension of k (=d') */
            vector<t_hash> h_functions;

            /*  The actual "hashtable" of the hypercube
                    Each index corresponds to one of the 2^k different integers
                    that can be created.
                    At each index there is a vector saving the points that hash to
                    that index.
             */
            vector<t_point>* _hashtable;

            // pathes to files
            string output_file;
            string data_file_path;
            string query_file_path;


        
            /* The Fs map, mapping each f_i(<int>) to another map structure(map<int,int>)
                which then maps hash values to their correspodning binary values
                e.g. fs[0] given as access to the mappings of h_1
                and then h_1(value) will return weather its value is 0 or 1
                according to the coin toss we do after each evaluation
            */
            map<int, map<int,int>> fs;

            /*
            *   @brief 
            *       Given a point, and a temp_evaluations array
            *       FILL the evaluations array with the evaluations
            *       from each h_i, and then check if the evaluation
            *       is mapped to a binary representation (f(h(p))).
            *       If the value has not been evaluated before, and therefore
            *       does not have a binary number to it, flip a coin, and save the result.
            *       If value has been evaluted before, do nothing.
            *   
            *   @param p The given point to be hashed
            *   @param temp_evaluations An integer array of k (=d') dimension, for each evalution to be saved on
            *   @return void
            */
            void _map_point_to_f_values(const vector<uint8_t>& p, int*& temp_evaluations );

            /*
            *   @brief 
            *       Evaluate the hash value (= hypercube vertex)
            *       for a given hash_evaluations array
            *       computed by the seperate evaluations for a point
            *   
            *   @param hash_evaluations A integer array of length k
            *   with the corresponding hash_value from evaluation saved
            *   in the index position
            *   e.g. at h[0] there must be saved h[0].evaluate(p) and so on...
            *   @return Integer representation of the accumulated hash from
            *   the f_i functions as: p -> [f_1(h_1(p)) f_2(h_2(p)) ... f_k(h_k(p))]
            */
            int _evaluate(int*& hash_evaluations);

            /*
            *   @brief 
            *       Complete training for a single point of ID=id.
            *       It computes the hash values of h, saving the f(h(p))
            *       evaluations if needed, and after evaluating the integer
            *       representation of that point, it later add it to the hypercube
            *       vertex corresponding to that integer representation
            *   
            *   @param p The point to be "trained"
            *   @param id The id of that point
            *   @return void
            */
            void _train_point(const vector<uint8_t>& p, int& id) ;

            void _set_metric(std::function<double(const std::vector<uint8_t>&, const std::vector<uint8_t>&)> metric) {
                this->metric = metric;
            }

            /*
            *   @brief 
            *       Given a starting vertex position represented by an integer, and a distance d,
            *       return a vector of integers, representing all neighbouring
            *       vertices that are d-hamming distance away from the starting vertex
            *   
            *   @param p The point to be "trained"
            *   @param id The id of that point
            *   @return vector of ints representing integers of vertex positions
            */
            vector<uint> _find_vertices_at_hamming_distance(int starting_vertex, int d);

            vector<vector<t_point>> find_neighboring_vertices_buckets(int starting_vertex);
       
    public:
           class hypercube_input {
            public:
            std::map<std::string, std::string> mapped_values;
            
            hypercube_input(int argc, char** argv) {
                std::vector<std::string> options = {
                    "d", "q", "k", "M", "probes", "o", "N", "R"
                };

                Get_long_opts optobj (options, Get_long_opts::options::SINGLE);
                optobj.parse_values(argc,(const char**) argv);

                for(const string& opt : options) {
                    if (!optobj.get_vector_from_option(opt.c_str()).empty()) {
                        mapped_values[opt] = optobj.get_vector_from_option(opt.c_str()).at(0);
                    } else {
                        mapped_values[opt] = "";
                    }
                }
            };
        };
    
            /*
            *   @brief 
            *       A enumeration for pre-defined metrics
            *       Select which option you want to use, and pass it
            *       as a parameters in the Hypercube() constructor
            *       as hc::Hypercube::Metric_functions::<OPTION>
            *  
            */
            enum Metric_functions {
                Manhattan,
                Euclidean
            };
            Hypercube(){};

            /*
            *   @brief 
            *       Initialize a Hypercube projection object, using 
            *       <metric> as the metric to be used for comparisons.
            *   
            *   @param input A hypercube_input object, containing all the parameters needed for the hypercube
            *   @param metric The metric to be used, either user-defined or pre-defined in the Metric_functions enum              
            *   @param k The dimension in which the points will be projected onto
            *   @param training_set A 2D vector of images(vectors of <uint8_t>), used to initialize the structure
            *   @param M The maximum number of points that will be checked
            *   @param probes Number of vertices of the hypercube on which we will search on
            *   @param N Number of nearest neighbours to be searched for
            *   @param R Number of points to be checked during the search  
            * @return
            *   
            *   
            */
            Hypercube(hypercube_input& input,std::function<double(const std::vector<uint8_t>&, const std::vector<uint8_t>&)> metric,const int window = 500);
            
            Hypercube(hypercube_input& input,Metric_functions mf) : Hypercube(input,empty) {
                switch (mf) {
                case Manhattan:
                    _set_metric(l1);
                    break;
                case Euclidean:
                    _set_metric(l2);
                    break;
                default:
                    throw std::invalid_argument("No metric found for that option");
                    break;
                }
            }

            ~Hypercube() {
                delete[] _hashtable;
            };

            /*
            *   @brief 
            *       Function to train a given dataset, aka initialize all structures used for later searching
            *   
            *   @param training_set A 2D vector of images(vectors of <uint8_t>), used to initialize the structure
            */
            void train_dataset(vector<vector<uint8_t>> &training_set);
            
            /*
            *   @brief 
            *       Given a query point, return a vector of pairs containing the N-nearest neighbours to that point
            * 
            *       @param query_point A 1D vector representing an image (vector of <uint8_t>)
            */
            vector<pair<int, double>> find_knn(const vector<uint8_t>& query_point);

            void search();
 
            void set_inputs(hypercube_input& inputs);
        
        vector<double> bf_nearest_neighbor(const vector<uint8_t>& ,const int);

        set<int> range_search(const vector<uint8_t>&,const int);
       
        void print_output(vector<pair<int,double>>&,vector<double>&,set<int>&,const string,const double,const double,const int) const;
        
        private:
            Metric_functions metric_selection;
};
