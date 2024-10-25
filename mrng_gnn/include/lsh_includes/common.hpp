// #include "headers.hpp"
#pragma once
#include "get_long_opts_single.hpp"
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
#include "structures.hpp"
#include <functional>
#include <unordered_set>
#include <algorithm>

int get_int_uniform_distribution(int , int);
double get_double_uniform_distribution(double, double);
double gen_number_normal_distribution(void);
std::vector<double> generate_random_normal_vector(const int);
void fill_random_normal_vector(std::vector<double>&, int);
double euclideanDistance(const std::vector<uint8_t>&, const std::vector<uint8_t>&);

/*
*   @brief 
*       This is a MNIST dataset specific function. Provided a path to a dataset
*       , it is read in BINARY mode, and a 2D binary std::vector of pixels is returned.
*       Each row represents a 1x784 (28*28) black and white image.
*   
*   @param path The path to the dataset file, opened in binary
*   @return A 2D array of binary data, representing the dataset, each row constitutes a image
*   If the file cannot be opened a exception is thrown.
*/
std::vector<std::vector<uint8_t>> get_pointset(std::string path);

/*
*   @brief 
*       This is a MNIST dataset specific function. Provided a path to a test dataset
*       , it is read in BINARY mode, and a random image of the dataset is returned,
*       used to later query the dataset.
*   
*   @param path The path to the test dataset file, opened in binary
*   @return A 1D binary std::vector representing a black and white image, randomly selected from
*   the train set. \n
*   If the file cannot be opened a exception is thrown.
*   
*/
std::pair<std::vector<uint8_t>,int> get_query_point(std::string path, std::string labels_path="");

// Function to calculate the Hamming distance between two integers
int hamming_distance(int x, int y);

// class input_interface{
//     int _d;
//     int _q;
//     int _k;
//     int _o;
//     int _n;
//     int _r;
//     std::map<std::string, std::string> mapped_values;
    
//     public:
//     input_interface(int argc, char** argv) {
//         // std::std::vector<std::string> options = {
//         //     "d", "q", "o", "N", "R", "k"
//         // };

//         // Get_long_opts optobj (options, Get_long_opts::options::SINGLE);
//         // optobj.parse_values(argc,(const char**) argv);

//         // for(const string& opt : options) {
//         //     if (!optobj.get_std::vector_from_option(opt.c_str()).empty()) {
//         //         mapped_values[opt] = optobj.get_std::vector_from_option(opt.c_str()).at(0);
//         //     } else {
//         //         mapped_values[opt] = "";
//         //     }
//         // }
//     }
// };

