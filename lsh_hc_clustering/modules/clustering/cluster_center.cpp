#include "headers.hpp"

void Cluster_center::_add_to_average(const cluster_item& item)  {
    vector<uint8_t> centroid = item.get_vector();

    for(size_t i = 0 ; i < centroid.size() ; i++) {
        centroid_point[i] = ( centroid_point[i] * (vectors_in_cluster.size()-1) + static_cast<double>(centroid[i]) ) / vectors_in_cluster.size();
        
        if (!(centroid_point[i] >= 0 && centroid_point[i] <= 255)) {
            cout << "Centroid value out of bounds: " << centroid_point[i] << endl;
        }
        assert(centroid_point[i] >= 0 && centroid_point[i] <= 255);
    }
}

void Cluster_center::_remove_from_average(const cluster_item& item) {

    vector<uint8_t> centroid = item.get_vector();

    for(size_t i = 0 ; i < centroid.size() ; i++) {
        if (vectors_in_cluster.size() - 1 == 0) {
            cout << "RE MALAKA XAZE DIAIREIS ME 0" << endl;
            cout << centroid_point[i] << endl;
        } else {
            // cout << "before evaluation " << centroid_point[i] <<endl;
            // cout << "before evaluation centroid[i] is " << static_cast<int>(centroid[i]) << endl;
            centroid_point[i] = ( centroid_point[i] * (vectors_in_cluster.size() ) - static_cast<double>(centroid[i]) ) / (vectors_in_cluster.size() - 1);
            // cout << i << endl;
            if (centroid_point[i] < 0) { centroid_point[i] = 0.0; }
            if (!(centroid_point[i] >= 0 && centroid_point[i] <= 255)) {
                cout << "Centroid value out of bounds: " << centroid_point[i] << endl;
                cout << "Removed value was " << static_cast<int>(centroid[i]) << endl;
                cout << "id " << item.get_id() << endl;
            }
            assert(centroid_point[i] >= 0 && centroid_point[i] <= 255);
        }
        
    }
}

/* @brief Returns a vector of <cluster_items> that are inside a Cluster_center class*/
vector<cluster_item> Cluster_center::get_cluster_items() const { return vectors_in_cluster; }

/* @brief Returns a 2D vector of the vectors in the Cluster_center */
vector<vector<uint8_t>> Cluster_center::get_cluster_vectors() const { 
    vector<vector<uint8_t>> result;
    for(const cluster_item& item : vectors_in_cluster) {
        result.push_back(item.get_vector());
    }
    return result;
}

/* @brief Returns the size of the cluster, thus how many elements it encompasses */
int Cluster_center::get_cluster_size() const { return vectors_in_cluster.size(); }

/* @brief Return the coordinates of the centroid of the cluster */
vector<uint8_t> Cluster_center::get_vector() const { 
    std::vector<uint8_t> to_return;
    to_return.reserve(point.size());

    for (const double& value : centroid_point) {
        // Convert the double value to uint8_t. You may want to add error handling
        // to handle cases where the double value is out of range for uint8_t.
        uint8_t convertedValue = static_cast<uint8_t>(value);
        to_return.push_back(convertedValue);
    }

    return to_return;
};

Cluster_center::Cluster_center(const vector<uint8_t> p, int id):id(id), point(p), mean_value(-1) { 

    centroid_point.reserve(point.size());
    for (const uint8_t& value : point) {
        double convertedValue = static_cast<double>(value);
        centroid_point.push_back(convertedValue);
    }
    
    // assert(centroid_point.size() == (28*28));

};

int Cluster_center::get_id() const { return id; };

/* @brief Append a given cluster_item to the vector */
void Cluster_center::append_in_cluster( const cluster_item& item ) { 
    vectors_in_cluster.push_back(item); 
    /* MacQueen find updating */
    // _update_midean(item);
}

/* @brief Clear the vector containing assigned vectors to the cluster */
void Cluster_center::flush_cluster() { vectors_in_cluster.clear(); };

/* @brief Compute a new centroid vector from the ones present in the cluster */
/* @return bool True if center changed, false if center was not changed */
bool Cluster_center::evaluate_new_center( std::function<double(const vector<uint8_t> &, const vector <uint8_t>&)> metric, int threshold) {
    // size_t data_dimension = 28*28;
    size_t data_dimension = point.size();
    assert(point.size() == (data_dimension));
    vector<uint8_t> new_center = vector<uint8_t>((data_dimension),0);
    
    // Construct the new centroid
    // for each coordinate 
    for(size_t i = 0; i < data_dimension; i++) {
        vector<uint8_t> temp;

        // Go through each item in the cluster
        // and add it to a temp vector to get compuate the median
        for(cluster_item& item: vectors_in_cluster) {
            vector<uint8_t> point = item.get_vector();
            temp.push_back(point[i]);
        }
        int median_position = temp.size() / 2;

        nth_element(temp.begin(), temp.begin()+median_position, temp.end());
        new_center[i] = temp[median_position];
    }

    assert(point.size() == new_center.size());

    if(metric(point, new_center) < threshold ) { return false; };

    return true;
}

void Cluster_center::_print_coordinates(std::ostream& output_stream) const {
    output_stream << "[ ";
    // assert(centroid_point.size() == (28*28));
    for(size_t i = 0; i < centroid_point.size() ; i++) {
        output_stream << centroid_point[i] << ", ";
    }
    output_stream << " ]";
}

double Cluster_center::evaluate_silhouette(double* a, double* b) {
    double average_si = 0;
    for(const cluster_item& item : vectors_in_cluster) {
        int index = item.get_id();
        average_si += (b[index]-a[index]) / std::max(b[index], a[index]);
    }

    assert(vectors_in_cluster.size() != 0);
    average_si /= vectors_in_cluster.size();
    return average_si;
}

void Cluster_center::fill_average_distance_of_vectors_to_vectors_in_same_cluster(double*& distances, std::function<double(const vector<uint8_t> &, const vector <uint8_t>&)> metric) {

    for(const cluster_item& item : vectors_in_cluster) {
        uint index = item.get_id();
        double total_distance = 0;
        
        for(const cluster_item& other_vectors : vectors_in_cluster) {
            double distance = metric(item.get_vector(), other_vectors.get_vector());
            total_distance += distance;
        }

        // -1 because we don't take into account this point,
        // distance from itself is already 0, so average is not affected
        distances[index] = total_distance / (vectors_in_cluster.size()-1);
    }            
}

void Cluster_center::latent_fill_average_distance_of_vectors_to_vectors_in_same_cluster(double*& distances, std::function<double(const vector<uint8_t> &, const vector <uint8_t>&)> metric,const vector<vector<uint8_t>> starting_space) {
    
    for(const cluster_item& item : vectors_in_cluster) {
        uint index = item.get_id();
        double total_distance = 0;
        
        vector<uint8_t> starting_point = starting_space[index];

        for(const cluster_item& other_vectors : vectors_in_cluster) {
            vector<uint8_t> other_point = starting_space[other_vectors.get_id()];
            double distance = metric(starting_point, other_point);
            total_distance += distance;
        }

        // -1 because we don't take into account this point,
        // distance from itself is already 0, so average is not affected
        distances[index] = total_distance / (vectors_in_cluster.size()-1);
    }            
}

double Cluster_center::average_distance_of_point_to_cluster(const cluster_item& point, std::function<double(const vector<uint8_t> &, const vector <uint8_t>&)> metric) {
    double total_distance = 0;
    for(const cluster_item& item : vectors_in_cluster) {
        total_distance += metric(item.get_vector(), point.get_vector());
    }

    return (double) (total_distance / vectors_in_cluster.size()) ; 
}

double Cluster_center::latent_average_distance_of_point_to_cluster(const cluster_item& point, std::function<double(const vector<uint8_t> &, const vector <uint8_t>&)> metric, const vector<vector<uint8_t>> starting_space) {
    double total_distance = 0;
    vector<uint8_t> starting_point = starting_space[point.get_id()];
    for(const cluster_item& item : vectors_in_cluster) {
        vector<uint8_t> other_point = starting_space[item.get_id()];
        total_distance += metric(starting_point, other_point);
    }
    
    return (double) (total_distance / vectors_in_cluster.size()) ; 
}

bool Cluster_center::remove_from_cluster(const cluster_item& item) {
    for (auto it = vectors_in_cluster.begin(); it != vectors_in_cluster.end();) {
        if (it->get_id() == item.get_id()) {
            it = vectors_in_cluster.erase(it); return true;
        } else {
            ++it;
        }
    }

    return false;
}

void Cluster::get_from_file() {
    std::ifstream configFile(configuration_file); // Replace "config.txt" with your file name

    number_of_clusters = 10;   
    number_of_vector_hash_tables = 3;
    number_of_vector_hash_functions = 4;
    max_number_M_hypercube = 10;
    number_of_hypercube_dimensions = 3;
    number_of_probes = 2;

    std::string line;
    while (std::getline(configFile, line)) {
        size_t found = line.find(':');
        if (found != std::string::npos) {
            std::string param = line.substr(0, found);
            std::string valueStr = line.substr(found + 1);
            valueStr.erase(0, valueStr.find_first_not_of(" \t\n\r\f\v")); // Trim leading whitespace

            if (param == "number_of_clusters") {
                try {
                    number_of_clusters = std::stoi(valueStr);
                } catch( invalid_argument& e ) {
                    cerr << "Please provide number" << endl;
                    exit(EXIT_FAILURE);
                }
            } else if (param == "number_of_vector_hash_tables") {
                try {
                    number_of_vector_hash_tables = std::stoi(valueStr);
                } catch( invalid_argument& e ) {
                    cerr << "Please provide number" << endl;
                    exit(EXIT_FAILURE);
                }
            } else if (param == "number_of_vector_hash_functions") {
                try { number_of_vector_hash_functions = std::stoi(valueStr); }
                catch( invalid_argument& e ) {
                    cerr << "Please provide number" << endl;
                    exit(EXIT_FAILURE);
                }
            } else if (param == "max_number_M_hypercube") {
                try{
                max_number_M_hypercube = std::stoi(valueStr);
                } catch( invalid_argument& e ) {
                    cerr << "Please provide number" << endl;
                    exit(EXIT_FAILURE);
                }
            } else if (param == "number_of_hypercube_dimensions") {
                try {
                    number_of_hypercube_dimensions = std::stoi(valueStr);
                }  catch( invalid_argument& e ) {
                    cerr << "Please provide number" << endl;
                    exit(EXIT_FAILURE);
                }
                k = number_of_hypercube_dimensions;
            } else if (param == "number_of_probes") {
                try {
                    number_of_probes = std::stoi(valueStr);
                }  catch( invalid_argument& e ) {
                    cerr << "Please provide number" << endl;
                    exit(EXIT_FAILURE);
                }
            }
        }
    }
    
    std::cout << "number_of_clusters: " << number_of_clusters << std::endl;
    std::cout << "number_of_vector_hash_tables: " << number_of_vector_hash_tables << std::endl;
    std::cout << "number_of_vector_hash_functions: " << number_of_vector_hash_functions << std::endl;
    std::cout << "max_number_M_hypercube: " << max_number_M_hypercube << std::endl;
    std::cout << "number_of_hypercube_dimensions: " << number_of_hypercube_dimensions << std::endl;
    std::cout << "number_of_probes: " << number_of_probes << std::endl;
    configFile.close();
}