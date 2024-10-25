#include "headers.hpp"
#pragma once

class LSH {
    private:
        
        /* Configuration for LSH structure */
        
        // pathes to files
        string output_file;
        string data_file_path;
        string query_file_path;

        string init_data_file_path;
        string init_query_file_path;
        bool latent;
        // Number of g functions
        int L;
        
        // Number of h functions per g
        int k;

        // Number of return neighbours
        int N;

        // Number of radius for range search
        int R;

        // Metric function to be used
        std::function<double(const std::vector<uint8_t>&, const std::vector<uint8_t>&)> metric;
        
        // Number of h functions it's going to be created
        int family_size;

        // Window size for the hash functions
        int W;
        
        // Training set
        vector<vector<uint8_t>> training_points;

        // H hash family
        vector<t_hash> hash_family;

        // G family
        vector<g_hash> g_hash_functions;

        // G family
        int get_rand_h_index() {
            return get_int_uniform_distribution(0,family_size);
        }

        // train the LSH
        void _training();

    public:

        // input prosessing
        class lsh_input {
            public:
            std::map<std::string, std::string> mapped_values;

            lsh_input(int argc, char** argv) {
                // std::vector<std::string> options = { "L" };
                std::vector<std::string> options = {
                    "d", "q", "o", "N", "R", "k", "L", "latent", "id", "iq"
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
        *   @brief Given a point p, return its closest neighbour
        *   
        *   @param query_point A 1D vector of binary values, representing a flattened image as a query image
        *   @return A std::pair<> of the actual image id of the point in the dataset, nearest to the query point (.first)
        *   and its calculated distance as defined by the metric function provided (.second)
        */

        pair<int,double> find_nearest_neighbour(vector<uint8_t>&);

        /*
        *   @brief Given a point p, return its K closest neighbours
        *   
        *   @param query_point A 1D vector of binary values, representing a flattened image as a query image
        *   @param k Number of closest neighbours to search for
        *   @return A vector of k-std::pair<> of the actual point neirest to the query point (.first)
        *   and its calculated distance as defined by the metric function provided (.second)
        */
       
        vector<pair<int,double>> find_k_nearest_neighbour(vector<uint8_t>&,const int);

        LSH(lsh_input& input_file,function<double(const std::vector<uint8_t>&, const std::vector<uint8_t>&)>,int family_size = 100,int window_size = 1000);
        
        /* @brief Performs the searching, KNN and Range search */
        void search();

        /* @brief Input processing*/
        void set_inputs(lsh_input& inputs);

        /*
        *   @brief Given a point query point, and a integer N
        *    return the N true nearest neighbours
        *   
        *   @param query_point A 1D vector of binary values, representing a flattened image as a query image
        *   @param k Number of closest neighbours to search for
        *   @return A vector of n-std::pair<> of the actual point neirest to the query point (.first)
        *   and its calculated distance as defined by the metric function provided (.second)
        */

        vector<double> bf_nearest_neighbor(const vector<uint8_t>&,const int);

        set<int> range_search(const vector<uint8_t>& ,const int);
       
        void print_output(vector<pair<int,double>>&,vector<double>&,set<int>&,const string,const double,const double,const int, const double=0) const;

        vector<double> new_bf_nearest_neighbor(const vector<uint8_t>& query,const vector<vector<uint8_t>>& real_points ,const int n);
};