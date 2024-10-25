#pragma once
#include "headers.hpp"
#include "hypercube.hpp"
#include "lsh.hpp"

class cluster_item {
    private:
        vector<uint8_t> _vector;
        size_t _id;
        bool _marked;
        /* is queued to be inserted to a cluster, used in reverse assignment to resolve merge conflicts */
        bool _queued;
    public:
        cluster_item(const vector<uint8_t> data, size_t id, bool marked = false) : _vector(data), _id(id), _marked(marked), _queued(false)  {   
        }

        void set_marked() { _marked = true; }
        void unset_marked() { _marked = false; }
        
        void set_queued() { _queued = true; }
        void unset_queued() { _queued = false; }
        
        bool is_queued() const { return _queued; }
        bool is_marked() const { return _marked; } 
        vector<uint8_t> get_vector() const { return _vector; }
        size_t get_id() const { return _id; }
};

class Cluster_center {
    private:
        /* The id and the coordinates(vector) of the centroid*/
        int id;
        vector<uint8_t> point;
        vector<double> centroid_point;
        double silhouette_value;
        double mean_value;
        /* A vector containing all cluster_items(vectors) assigned to this cluster center*/
        vector<cluster_item> vectors_in_cluster;
        void _compute_silhouette();

    public:

        double evaluate_silhouette(double* a, double* b);
        void _add_to_average(const cluster_item& item);
        void _remove_from_average(const cluster_item& item );
        Cluster_center(const vector<uint8_t> p, int id);

        double get_silhouette() const { return silhouette_value; }
        void _print_coordinates (std::ostream& output_stream) const;

        void fill_average_distance_of_vectors_to_vectors_in_same_cluster(double*& distances, std::function<double(const vector<uint8_t> &, const vector <uint8_t>&)> metric);
        

        double average_distance_of_point_to_cluster(const cluster_item& point, std::function<double(const vector<uint8_t> &, const vector <uint8_t>&)> metric);
        /* @brief Returns a vector of <cluster_items> that are inside a Cluster_center class*/
        vector<cluster_item> get_cluster_items() const;
        
        /* @brief Returns a 2D vector of the vectors in the Cluster_center */
        vector<vector<uint8_t>> get_cluster_vectors() const;

        /* @brief Returns the size of the cluster, thus how many elements it encompasses */
        int get_cluster_size() const;

        /* @brief Return the coordinates of the centroid of the cluster */
        vector<uint8_t> get_vector() const;

        vector<uint8_t> get_vector_() const {return point;}


        int get_id() const;

        /* @brief Append a given cluster_item to the vector */
        void append_in_cluster( const cluster_item& item );
        
        /* @brief  Remove a given cluster_item from the vector */
        bool remove_from_cluster(const cluster_item& item) ;
        
        /* @brief Clear the vector containing assigned vectors to the cluster */
        void flush_cluster();

        /* @brief Compute a new centroid vector from the ones present in the cluster */
        /* @return bool True if center changed, false if center was not changed */
        bool evaluate_new_center( std::function<double(const vector<uint8_t> &, const vector <uint8_t>&)> metric, int threshold);

        void latent_fill_average_distance_of_vectors_to_vectors_in_same_cluster(double*&, std::function<double(const vector<uint8_t> &, const vector <uint8_t>&)>,const vector<vector<uint8_t>>);

        double latent_average_distance_of_point_to_cluster(const cluster_item&, std::function<double(const vector<uint8_t> &, const vector <uint8_t>&)>, const vector<vector<uint8_t>>);
};



class Cluster {
    public:
        enum method {LLOYDS, RA_LSH, RA_Hypercube};
    private:
        
        enum method internal_method;

        /* @brief returns the distance from a point to the nearest center, according to a metric*/
        double find_near_center(
            const vector<uint8_t>& point, 
            const vector<Cluster_center>& centers,
            std::function<double(const vector<uint8_t> &, const vector <uint8_t>&)> metric
            );

        double evaluate_probability(const double*& d,const int r);
        int bin_search(const vector<double>& array, const double index) ;
        
        string input_file;
        string output_file;
        string configuration_file;
        bool optional;
        int number_of_clusters;
        int number_of_vector_hash_tables;
        int number_of_vector_hash_functions;
        int max_number_M_hypercube;
        int number_of_hypercube_dimensions;
        int number_of_probes;

        int k;

        /* */
        vector<vector<uint8_t>> pointset;
        /* A vector of class cluster_center, containing all information about a cluster */
        vector<Cluster_center> centers;
        /* The overall dataset, a vector of cluster_item type*/
        vector<cluster_item> items;

        /* Silhouette values saved for each cluster, and the overall average silhouette value */
        vector<double> silhouettes;
        double overall_silhouette;

        /* LSH or Hypercube used for reverse assignment*/
        LSH* lsh;
        Hypercube* hc;

        /* The metric function*/
        std::function<double(const vector<uint8_t> &, const vector <uint8_t>&)> metric;

        /* Saves a mapping between image (vector) id to -> cluster id */
        map<int,int> point_to_cluster;

        /*
            @brief Empties each cluster center, which holds the points assigned to it, resetting it
        */
        void _empty_centroid_clusters();
        
        /*
            @brief Reset the markings of all the cluster_items
        */
        void _reset_item_marking();
        /*
            @brief Execution code for lloyd's
        */
        void _lloyds(bool &);

        /*
            @brief Reverse assignment of the points, depending on mode 
        */

       void _reverse_assignment();


        /*
            @brief Given a cluster item, return a pair<> holding the index of the closest/nearest centroid item
            and its distance
            @param point A cluster_item item, representing a point in a cluster
            @return pair<int,double> The index of the nearest centroid, and the distance from it
        */
        pair<int,double> _find_nearest_center(const cluster_item& point);

    public:

        // input processing    
        class cluster_input {
            public:
            std::map<std::string, std::string> mapped_values;
            
            cluster_input(int argc, char** argv) {
                std::vector<std::string> options = {
                    "i", "c", "o", "complete", "m"
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

    int minimum_dist_between_centers();

    pair<int,int> find_two_closest_centers(const vector<uint8_t>& p);


    Cluster(cluster_input& input, std::function<double(const vector<uint8_t> &, const vector <uint8_t>&)> metric) ;
    ~Cluster() ;
    
    /* @brief Apply kmeans++ initialization on the dataset */
    void kmeanspp(const vector<vector<uint8_t>> &points, std::function<double(const vector<uint8_t> &, const vector <uint8_t>&)> metric);

    /* @brief Apply clustering algorithm */
    void apply_clustering();

    /* @brief Sets cluster object inputs from the cluster_input objects, used for initilization */
    void set_inputs(cluster_input& inputs);
    void get_from_file();
    
    void print_results(const double) const;
    
    /* Computes silhouettes for each cluster */
    void _compute_silhouettes();

    void _compute_latent_silhouettes(string dataset);  
};
