#include "headers.hpp"
#include "get_long_opts_single.hpp"

pair<int,double> Cluster::_find_nearest_center(const cluster_item& point) {
        
    double min_distance = numeric_limits<double>::infinity();
    int index = -1;
    
    for (int i = 0; i < k ; i++) {
        // cout << centers[i].get_vector().size() <<  " " << point.get_vector().size() << endl;
        if (metric(centers[i].get_vector(), point.get_vector()) < min_distance) {
            min_distance = metric(centers[i].get_vector(), point.get_vector());
            index = i;
        }
    }
    return make_pair(index,min_distance);
}   

void Cluster::print_results(const double time) const {
    std::ofstream outputFile(output_file);
    if(internal_method == RA_LSH) {
        outputFile << "Algorithm: Range Search LSH" << endl;
    }
    if(internal_method == RA_Hypercube) {
        outputFile  << "Range Search Hypercube" << endl;
    }
    if(internal_method == LLOYDS) {
        outputFile  << "Algorithm: Lloyds" << endl;
    }
    outputFile << endl;
    for(int i = 0; i < number_of_clusters ; i++) {
        outputFile << "Cluster-"<<i+1<< " {size: "<<centers[i].get_cluster_size() <<", centroid:  ";
        centers[i]._print_coordinates(outputFile);
        outputFile << " }" << endl << endl;;
    }

    outputFile << endl;
    outputFile <<"clustering_time: " << time  << " seconds" << endl << endl;
    outputFile <<"Silhouette: [";
    for(int i =0; i < k; i++) {
        outputFile << silhouettes[i] << ", ";
    }  
    outputFile << " " << overall_silhouette ; 
    outputFile << "]" << endl;
    
    if (optional) {
        for(int i = 0; i < number_of_clusters ; i++) {
            outputFile << "CLUSTER-" <<i+1<< " {centroid, ";
            for(int j = 0; j < centers[i].get_cluster_size(); j++) {
                outputFile << centers[i].get_cluster_items()[j].get_id() << ", ";
            }
            outputFile << "}" << endl;
        }
    }
}



/* 

//// CLUSTER


*/

double Cluster::find_near_center(
    const vector<uint8_t>& point, 
    const vector<Cluster_center>& centers,
    std::function<double(const vector<uint8_t> &, const vector <uint8_t>&)> metric) {
    
    double min_distance = numeric_limits<double>::infinity();
    
    for (const Cluster_center& center : centers) {
        if (center.get_vector_().size() != point.size()) { 
            // cout << center.get_vector().size() << " " << point.size() << endl;
        }
        if (metric(center.get_vector_(), point) < min_distance) {
            min_distance = metric(center.get_vector(), point);
        }
    }
    
    return min_distance;
}

double Cluster::evaluate_probability(const double*& d,const int r) {
    double sum = 0;
    double max_d = 0;
    for(int i = 0; i < r; i++) {
        sum += pow(d[i],2);
        if(d[i] > max_d) {
            max_d = d[i];
        }
    }
    return sum / max_d;
}

int Cluster::bin_search(const vector<double>& array, const double index) {
    size_t left = 0;
    size_t right = array.size() - 1;
    
    while (left < right) {
        int mid = left + (right - left) / 2;
        
        if (index <= array[mid]) {
            right = mid;
        } else {
            left = mid + 1;
        }
    }
    
    return left;
}

/*
    @brief Empties each cluster center, which holds the points assigned to it, resetting it
*/
void Cluster::_empty_centroid_clusters() { for(Cluster_center& cc : centers) { cc.flush_cluster(); } };

/*
    @brief Reset the markings of all the cluster_items
*/
void Cluster::_reset_item_marking() { for(cluster_item& item: items) { item.unset_marked(); }; }


/*
*   @brief Lloyds algorithm implementation: assign each datapoint to the closest center
*/
void Cluster::_lloyds(bool &d) {

    // cout << "Gets in lloyds" << endl;
    vector<bool> changed_centers = vector<bool>(k,false);
    
    for(cluster_item& item : items) {
        // If the item has not been added to a cluster (= not marked)
        if (!item.is_marked()) {

            // Find its nearest center
            pair<int,double> closest = _find_nearest_center(item);
            // cout << "Out of nearest_center" << endl;
            // cout << "daino" << endl;
            assert(closest.first != -1);
            // Add it to that cluster
            centers[closest.first].append_in_cluster(item);

            point_to_cluster[item.get_id()] = closest.first;

            centers[closest.first]._add_to_average(item);
            // Now mark it as it has been added to a cluster
            item.set_marked();
            changed_centers[closest.first] = true;
        }
        for(int i = 0; i < k; i++) {
            if (changed_centers[i] == true) d = false;
        }
    }
}

void Cluster::kmeanspp(const vector<vector<uint8_t>> &points, std::function<double(const vector<uint8_t> &, const vector <uint8_t>&)> metric) {

    // Initialize list of centers
    int n = points.size();
    int t = number_of_clusters;
    set<int> centroid_point_ids;
    // double *dist = new double[n];

    // Get a random starting point from the dataset
    int result = get_int_uniform_distribution(0, n - 1);
    centroid_point_ids.insert(result);
    Cluster_center center(points[result], result);
    centers.push_back(center);
    // Very first point is assigned to the 0th index
    point_to_cluster[result] = 0;
    // cout << "Chosen center (" << 0 <<  ") is: " << result << endl;

    // Iterate through the dataset
    // Assign each point to the closest center


    // While k (=t) clusters have not beed selected
    for(int cluster = 1; cluster < t; cluster++) {
        vector<double> dist(n, 0.0);
        vector<double> p(n,0.0);
        double sum = 0;

        for(int i = 0; i < n; i++) {
            dist[i] = find_near_center(points[i], centers, metric);         
        }
        for(int i = 0; i < n; i++) {
            sum += dist[i];
            p[i] = sum;
        }

        float random = get_double_uniform_distribution(0, sum);
        int selected_next_centroid = bin_search(p,random);

        centers.emplace_back(Cluster_center(points[selected_next_centroid], selected_next_centroid));
        point_to_cluster[selected_next_centroid] = cluster ;
        // cout << "Chosen center (" << cluster <<  ") is: " << selected_next_centroid << endl;
    }


};


void Cluster::_reverse_assignment(void) {

    int iterations = 0;
    int max_iterations = 10;
    /* Image id to -> cluster id */
    int total_assigned_images = 0;
    map<int, int> assigned_clusters;
    int radius = minimum_dist_between_centers() / 2;
    int centroids_expanded = 0;

    do {
        centroids_expanded = 0;
        set<int> range_search_res;
        
        for(int i = 0; i < k; i++) {

            if (internal_method == RA_LSH) {
                const Cluster_center cent = centers[i];
                range_search_res = lsh->range_search(cent.get_vector(), radius);
            } else if (internal_method == RA_Hypercube) {
                const Cluster_center cent = centers[i];
                range_search_res = hc->range_search(cent.get_vector(), radius);
            }

            if (range_search_res.size()>0) centroids_expanded++;
            
            for(const int item_in_range : range_search_res) {
    
                assert(!(item_in_range<0));

                assert(item_in_range >= 0 && item_in_range <= 59999);

                cluster_item& cl_item = items[item_in_range];

                if (cl_item.is_marked()) continue;

                if (cl_item.is_queued()) {
                    int assigned_cluster_index = assigned_clusters[item_in_range];
                    if (assigned_clusters[item_in_range] == i) { continue; }

                    double distance_to_assigned = metric(cl_item.get_vector(), centers[assigned_cluster_index].get_vector());
                    double distance_to_candidate = metric(cl_item.get_vector(), centers[i].get_vector());
                    
                    if (distance_to_candidate < distance_to_assigned) {
                        // cout << "Change incoming: " << item_in_range << " -> " << assigned_clusters[item_in_range] << endl;
                        int previous_cluster = assigned_clusters[item_in_range];
                        assert(previous_cluster == assigned_cluster_index);
                        // cout << "Vector should change cluster, from " << previous_cluster << " and go to " << i << endl;
                        assigned_clusters[item_in_range] = i;
                    }

                } else {
                    cl_item.set_queued();

                    /* Add mapping of image_id (vector_id) -> cluster id */
                    // assigned_clusters[item_in_range] = i;
                    auto it = assigned_clusters.find(item_in_range);

                    if (it != assigned_clusters.end()) {
                        // Key exists, update the value
                        it->second = i;
                    } else {
                        // Key doesn't exist, handle the case as needed
                        // For example, you can print an error message or skip it
                        // or you can choose to insert a new key-value pair:
                        assigned_clusters.insert(make_pair(item_in_range, i));
                    }

                }
            }    
        }

        for(auto x : assigned_clusters) {
            // x.first -> image id
            // x.second -> cluster id
            cluster_item& cl_item = items[x.first];
            cl_item.set_marked();
            cl_item.unset_queued();
            centers[x.second].append_in_cluster(cl_item);
            centers[x.second]._add_to_average(cl_item);
            total_assigned_images++;
        }

        // cout <<  "Vectors in map: " << assigned_clusters.size() << endl;
        // for (int i = 0; i < k; i++) { cout << "Cluster " << i << " assigned: " << centers[i].get_cluster_size() << " with reverse assignment" << endl; }
        
        radius *= 2;
        iterations++;
        assigned_clusters.clear();
        // cout << "Percentage " << static_cast<double>(static_cast<double>(total_assigned_images) / static_cast<double>(items.size())) << "%" << endl;
    } while( total_assigned_images <= 0.7 * items.size() && iterations < max_iterations );

    // cout << "Assigned " << total_assigned_images << " out of " << items.size() << endl; 
    bool f = true;
    
    _lloyds(f);
}

pair<int,int> Cluster::find_two_closest_centers(const vector<uint8_t>& p) {
    // Cluster_center closest();
    int closest = 0;
    int second_closest = 1;
    double closest_dist = metric(p, centers[0].get_vector());
    double second_closest_dist = metric(p, centers[1].get_vector());

    for (size_t i = 2; i < centers.size(); ++i) {
        double dist = metric(p, centers[i].get_vector());
        if (dist < closest_dist) {
            second_closest = closest;
            second_closest_dist = closest_dist;

            closest = i;
            closest_dist = dist;
        } else if (dist < second_closest_dist) {
            second_closest = i;
            second_closest_dist = dist;
        }
    }

    return make_pair(closest,second_closest);
}



void Cluster::apply_clustering() {

    // Begin clustering with the flag set to true
    bool changed_centers = true;
    int max_iter = 5;
    int iterations = 0;

    while(changed_centers && iterations < max_iter) {
        
        changed_centers = false;

        _reset_item_marking();

        vector<size_t> old_sizes = vector<size_t>(k);
        vector<size_t> new_sizes = vector<size_t>(k);

        for( size_t i = 0; i < centers.size();i++) { old_sizes[i] = centers[i].get_cluster_size(); }

        _empty_centroid_clusters();

        for(size_t i = 0; i < centers.size();i++) {
            assert(centers[i].get_cluster_size() == 0);
        }
        
        for(const cluster_item& item : items ) {
            assert(!item.is_marked());
            assert(!item.is_queued());
        }

        if (internal_method == LLOYDS) {
            _lloyds(changed_centers);
        }
        else _reverse_assignment();

        for(cluster_item& item : items) { item.unset_queued(); }
        for (size_t i = 0; i < centers.size();i++) { new_sizes[i] = centers[i].get_cluster_size();}
        for (size_t i = 0; i < centers.size(); i++ ) { if (new_sizes[i] != old_sizes[i]) changed_centers = true; }
        // for (size_t i = 0; i < centers.size(); i++) {
        //     cout << "Center " << i << " has: " << centers[i].get_cluster_size() << endl;
        // }
        // cout << endl;
        iterations++;
    }

}

void Cluster::_compute_latent_silhouettes(string dataset) {
    
    vector<vector<uint8_t>> output = get_pointset(dataset);
    double* a = new double[items.size()];
    double* b = new double[items.size()];
    double* s = new double[items.size()];

    for(int cluster = 0 ; cluster < k ; cluster++) {
        centers[cluster].latent_fill_average_distance_of_vectors_to_vectors_in_same_cluster(a, metric, output);
    }

    for(const cluster_item& item : items) {
        pair<int,int> closest_clusters = find_two_closest_centers(item.get_vector());
        assert(closest_clusters.second != -1);
        assert(closest_clusters.first != -1);
        int i = item.get_id();
        assert(closest_clusters.first >= 0 && closest_clusters.first <= k-1);
        assert(closest_clusters.second >= 0 && closest_clusters.second <= k-1);
        b[i] = centers[closest_clusters.second].latent_average_distance_of_point_to_cluster(item, metric, output);
    }

    for(const cluster_item& item : items) {
        uint i = item.get_id();
        s[i] = (b[i] - a[i]) / std::max(b[i],a[i]);
    }

    double average = 0;
    for(size_t i = 0; i < centers.size(); i++) {
        silhouettes[i] = centers[i].evaluate_silhouette(a,b);
        average += silhouettes[i];
    }

    assert(centers.size() != 0);
    average /= centers.size();
    overall_silhouette = average;
    delete[] a; delete[] b; delete[] s;
}
void Cluster::_compute_silhouettes() {

    double* a = new double[items.size()];
    double* b = new double[items.size()];
    double* s = new double[items.size()];

    for(int cluster = 0 ; cluster < k ; cluster++) {
        centers[cluster].fill_average_distance_of_vectors_to_vectors_in_same_cluster(a, metric);
    }

    for(const cluster_item& item : items) {
        pair<int,int> closest_clusters = find_two_closest_centers(item.get_vector());
        assert(closest_clusters.second != -1);
        assert(closest_clusters.first != -1);
        int i = item.get_id();
        assert(closest_clusters.first >= 0 && closest_clusters.first <= k-1);
        assert(closest_clusters.second >= 0 && closest_clusters.second <= k-1);
        b[i] = centers[closest_clusters.second].average_distance_of_point_to_cluster(item, metric);
    }

    for(const cluster_item& item : items) {
        uint i = item.get_id();
        s[i] = (b[i] - a[i]) / std::max(b[i],a[i]);
    }

    double average = 0;
    for(size_t i = 0; i < centers.size(); i++) {
        silhouettes[i] = centers[i].evaluate_silhouette(a,b);
        // cout << "s[" << i << "] = " << silhouettes[i] << endl;
        average += silhouettes[i];
    }

    assert(centers.size() != 0);
    average /= centers.size();
    overall_silhouette = average;
    // for(int i = 0; i < k; i++) {
    //     // cout << "s[" << i << "] = " << silhouettes[i] << endl;
    // }
    // cout << "sTotal=" << overall_silhouette << endl;
    
    delete[] a; delete[] b; delete[] s;
}

int Cluster::minimum_dist_between_centers() {
    double min_distance = numeric_limits<double>::infinity();
    for(const Cluster_center& cc: centers) {
        for(const Cluster_center& rest_centers: centers) {
            double distance = metric(cc.get_vector(), rest_centers.get_vector());
            if (distance != 0) {
                if (distance < min_distance) {
                    min_distance = distance;
                }
            }    
        }
    }

    assert(min_distance != numeric_limits<double>::infinity());
    return static_cast<int>(min_distance);
}

Cluster::Cluster(cluster_input& input, std::function<double(const vector<uint8_t> &, const vector <uint8_t>&)> metric) : 
    silhouettes(std::vector<double>(k, 0)), overall_silhouette(0),lsh(nullptr), hc(nullptr), metric(metric) {
    set_inputs(input);
    get_from_file();

    if(internal_method == RA_LSH) {
        // Define the command-line arguments as separate strings
        std::vector<std::string> args = {
            "-d", input_file,
            "-q", "",
            "-k", std::to_string(number_of_vector_hash_functions)
        };

        // Create a vector of char* pointers
        std::vector<char*> argv;
        for (auto& arg : args) {
            argv.push_back(&arg[0]);
        }
        argv.push_back(nullptr); // Null-terminate the array
        
        // Create LSH::lsh_input object using the argv data
        LSH::lsh_input ip(static_cast<int>(argv.size()) - 1, argv.data());

        lsh = new LSH(ip, metric);
    } if(internal_method == RA_Hypercube) {
        // Define the command-line arguments as separate strings
        std::vector<std::string> args = {
            "-d", input_file,
            "-q", "",
            "-M", to_string(max_number_M_hypercube),
            "-probes" , to_string(number_of_probes),
            "-o" , output_file,
            "-k", to_string(k)
        };

        // Create a vector of char* pointers
        std::vector<char*> argv;
        for (auto& arg : args) {
            argv.push_back(&arg[0]);
        }
        argv.push_back(nullptr); // Null-terminate the array
        
        // Create LSH::lsh_input object using the argv data
        Hypercube::hypercube_input ip(static_cast<int>(argv.size()) - 1, argv.data());

        hc = new Hypercube(ip, metric);
    } 

    // cout << "Method chosen is " << internal_method << endl;
    k = 10;

    try {
        pointset = get_pointset(input_file);   
    } catch( std::runtime_error& e ) 
    {
        std::cout << e.what() << std::endl;
    }

    for(size_t i = 0 ; i < pointset.size(); i++) {
        items.emplace_back(cluster_item(pointset[i], i));
    };

    // TODO:: LSH OR HC Initilization

    kmeanspp(pointset,metric);


    };

Cluster::~Cluster() {
    if(internal_method == RA_LSH) delete lsh;
    if(internal_method == RA_Hypercube) delete hc; 
}

void Cluster::set_inputs(cluster_input& inputs) {
        
        std::string dataset_path;
        this->input_file = inputs.mapped_values["i"];
        this->output_file = inputs.mapped_values["o"];
        this->configuration_file = inputs.mapped_values["c"];
        if(inputs.mapped_values["complete"] == "Yes" || inputs.mapped_values["complete"] == "YES" || inputs.mapped_values["complete"] == "yes" || inputs.mapped_values["complete"] == "Y" || inputs.mapped_values["complete"] == "y") 
        { this->optional = true;} else { this->optional = false;} 

        // cout << inputs.mapped_values["m"];

        if (inputs.mapped_values["m"] == "LSH") {
            this->internal_method = RA_LSH;
        } else if (inputs.mapped_values["m"] == "Classic") {
            this->internal_method = LLOYDS;
        } else if (inputs.mapped_values["m"] == "Hypercube") {
            this->internal_method = RA_Hypercube;
        } else {
            cout << "Invalid method provided. Please select: LSH or Classic or Hypercube" << endl;
            exit(EXIT_FAILURE);
        }
    }