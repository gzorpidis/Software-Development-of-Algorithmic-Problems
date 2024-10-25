#include "headers.hpp"

Hypercube::Hypercube(hypercube_input& input,std::function<double(const std::vector<uint8_t>&, const std::vector<uint8_t>&)> m, 
    const int window_size) :  W(window_size), metric(m)
{

    set_inputs(input);

    int hashtable_size = pow(2, k);

    _hashtable = new vector<t_point>[hashtable_size];

     try {
        training_set = get_pointset(data_file_path);   
    } catch( std::runtime_error& e ) {
        std::cout << e.what() << std::endl;
    }

    // int n = training_set.size();
    int size = training_set[0].size();
    for(int i = 0; i < k; i++) {
        h_functions.emplace_back(t_hash(size, W, i));
    }

    train_dataset(training_set); 
    
}

void Hypercube::train_dataset(vector<vector<uint8_t>>& training_set) {
    int id = 0;
    for(const vector<uint8_t>& p : training_set) {
        _train_point(p, id);
        id++;
    }

    int total = 0;
    int hashtable_size = pow(2, k);
    for(int i = 0; i < hashtable_size; i++) {
        total += _hashtable[i].size();
    }
    assert(total == static_cast<int>(training_set.size()));
}

void Hypercube::_train_point(const vector<uint8_t>& p, int& id) {
   
    int* temp_evaluations = new int[k];

    _map_point_to_f_values(p, temp_evaluations);

    int hash_value = _evaluate(temp_evaluations);

    _hashtable[hash_value].push_back(t_point(hash_value,p,id));

    delete[] temp_evaluations;
}

vector<vector<t_point>> Hypercube::find_neighboring_vertices_buckets(int starting_hc_vertex) {
    vector<vector<t_point>> neighboring_vertices_buckets;
    set<int> vertices_inside;
    // cout << "Starting vertex: " << starting_hc_vertex << endl;
    // Go from distance 0 to k
    for(int i = 0; i < k + 1; i++) {
        // Search incrementally the distance from 0 up to k bits, and find the neighbouring vertices
        vector<uint> neighboring_vertices = _find_vertices_at_hamming_distance(starting_hc_vertex, i);
        if (vertices_inside.size() == probes) break;

        for(const int& vertex : neighboring_vertices) {
            auto it = vertices_inside.find(vertex);
            if (it == vertices_inside.end()) {
                vertices_inside.insert(vertex);
                neighboring_vertices_buckets.push_back(_hashtable[vertex]);
                if (vertices_inside.size() == probes) break;
            }
            if (vertices_inside.size() == probes) break;
        }
    }

    assert(vertices_inside.size() == probes);
    return neighboring_vertices_buckets;
}

vector<pair<int, double>> Hypercube::find_knn(const vector<uint8_t>& query_point) {
    
    // Create the vector, and initialize it with empty buckets, and infinite distance
    vector<pair<int,double>> nn;
    
    // Using for dublicates
    set<int> set_nn;
    for(int i = 0; i < N; i++) {
        nn.push_back(make_pair(0,numeric_limits<int>::max()));
    }
    
    int* temp_evaluations = new int[k];
    _map_point_to_f_values(query_point, temp_evaluations);
    int starting_hc_vertex = _evaluate(temp_evaluations);

    vector<vector<t_point>> neighboring_vertices_buckets = find_neighboring_vertices_buckets(starting_hc_vertex);

    class compare_pair_by_distance {
        public:
        bool operator()(const std::pair<int, double>& a, const std::pair<int, double>& b) const {
            return a.second < b.second;
        }
    };

    uint compared_points = 0; 

    // For each vertex (which contains points), and for each of their points
    // compare using the metric function and update the vector accordingly
    for(const vector<t_point>& hc_vertex : neighboring_vertices_buckets ) {
        for(const t_point& candidate_neighbour : hc_vertex) {

            double distance = metric(query_point, candidate_neighbour.get_point());

            if (distance == 0) {
                cout << "\t\tFound distance 0!!!" << candidate_neighbour.get_img_id() << endl;
            }
            // Compare this distance to the distance in the last place
            // and then sort, so at the last position the point that 
            // is most distant is replaced with some point that is closer
            if (distance < nn[N-1].second && set_nn.find(candidate_neighbour.get_img_id()) == set_nn.end()) {
                nn[N-1].first = candidate_neighbour.get_img_id();
                nn[N-1].second = distance;
                sort(nn.begin(), nn.end(), compare_pair_by_distance{});
            }
            
            // Check weather we have surpassed the threshold for points to be searched
            if (M > 0 && ++compared_points >= M) break;
        }
    }

    return nn;
}

// Pre-defined metrics that can be used for quick initialization
double Hypercube::l2(const vector<uint8_t>& vector1, const vector<uint8_t>& vector2) {
    // Check if the vectors have the same dimension
    if (vector1.size() != vector2.size()) {
        throw std::invalid_argument("Vectors must have the same dimension");
    }
    
    double distance = 0.0;

    for (size_t i = 0; i < vector1.size(); ++i) {
        // Calculate the squared difference between corresponding elements
        double diff = static_cast<double>(vector1[i]) - static_cast<double>(vector2[i]);
        distance += diff * diff;
    }

    // Take the square root of the sum of squared differences
    return std::sqrt(distance);
}

double Hypercube::l1(const std::vector<uint8_t>& vector1, const std::vector<uint8_t>& vector2) {
    // Check if the vectors have the same size
    if (vector1.size() != vector2.size()) {
        throw std::invalid_argument("Vectors must have the same size.");
    }

    double distance = 0.0; // Initialize the distance
    for (size_t i = 0; i < vector1.size(); ++i) {
        distance += std::abs(static_cast<double>(vector1[i]) - static_cast<double>(vector2[i]));
    }

    return distance;
};

void Hypercube::set_inputs(hypercube_input& inputs) {
    std::string dataset_path;

    /* if -d is not given in options, ask for it*/
    if (inputs.mapped_values["d"] == "") {
        std::cout << "Give path to dataset " << endl;
        std::cin >> dataset_path;
        inputs.mapped_values["d"] = dataset_path;
    }

    /* Check each initilization option, if not given given default values */
    if (inputs.mapped_values["k"] == "") { inputs.mapped_values["k"] = "14"; }
    if (inputs.mapped_values["M"] == "") { inputs.mapped_values["M"] = "10"; }
    if (inputs.mapped_values["probes"] == "") { inputs.mapped_values["probes"] = "2"; }
    if (inputs.mapped_values["N"] == "") { inputs.mapped_values["N"] = "1"; }
    if (inputs.mapped_values["R"] == "") { inputs.mapped_values["R"] = "10000"; }
    
    if (inputs.mapped_values["q"] == "") { this->query_file_path = ""; } else { this->query_file_path = inputs.mapped_values["q"]; }
    
    /* if -o is not given in options, ask for it*/
    if (inputs.mapped_values["o"] == "") {
        cout << "Give output file path/name " << endl;
        cin >> dataset_path;
        inputs.mapped_values["o"] = dataset_path;
    }

    this->output_file = inputs.mapped_values["o"];
    this->data_file_path = inputs.mapped_values["d"];
    
    this->k = std::stoi(inputs.mapped_values["k"]);
    this->N = std::stoi(inputs.mapped_values["N"]); 
    this->R = std::stoi(inputs.mapped_values["R"]);
    this->M = std::stoi(inputs.mapped_values["M"]);
    this->probes = std::stoi(inputs.mapped_values["probes"]);

    this->init_data_file_path = inputs.mapped_values["id"];
    this->init_query_file_path = inputs.mapped_values["iq"];


    if (inputs.mapped_values["latent"] == "" ||
        inputs.mapped_values["latent"] == "No" ||
        inputs.mapped_values["latent"] == "N" ||
        inputs.mapped_values["latent"] == "no" || 
        inputs.mapped_values["latent"] == "n" ||
        inputs.mapped_values["latent"] == "NO"
    ) {
        this->latent = false;
    } else if (
        inputs.mapped_values["latent"] == "Yes" ||
        inputs.mapped_values["latent"] == "Y" ||
        inputs.mapped_values["latent"] == "yes" || 
        inputs.mapped_values["latent"] == "y" ||
        inputs.mapped_values["latent"] == "YES"
    ){
        this->latent = true;
        if (init_data_file_path == "" || init_query_file_path == "") {
            cerr << "Latent is true, but initial dataset or query file is empty" << endl;
            throw std::runtime_error("Latent is true, but initial dataset or query file is empty\n");
        }
    } else {
        cerr << "Latent must be either empty or No or Yes" << endl;
        throw std::runtime_error("Latent must be either empty or No or Yes\n");
    }
}

int Hypercube::_evaluate(int*& hash_evaluations) {

    uint hash_value = 0;

    for(int i = 0; i < k; i++) {

        int f_result = fs[i][hash_evaluations[i]];

        // XOR with 1 to assign 1 in the current bit
        // f_result is 1, else if 0, not need to do anything
        if (f_result == 1) { hash_value |= 1; }

        // Shift to the left fot the next value
        // but do not shift left if we are not the 
        // last iteration
        if (i < k - 1) { hash_value <<= 1; }

    }
    return hash_value;
}

vector<uint> Hypercube::_find_vertices_at_hamming_distance(int starting_vertex, int d) {
    vector<uint> result;
    // Consider all possible combinations of 'd' bits
    for (uint i = 0; i < pow(2,k); i++) {
        // Create a candidate by flipping 'd' bits in the given integer
        if (hamming_distance(starting_vertex, i) == d)
            result.push_back(i);
    }

    return result;
}

void Hypercube::search() {
    pair<vector<uint8_t>,int> return_tuple;
    vector<uint8_t> query_point;
    int query_number;

    if (query_file_path == "") {
        cout << "Give path to dataset" << endl;
        cin >> this->query_file_path;
    }
    try {
        return_tuple = get_query_point(query_file_path);
        query_point = return_tuple.first;
        query_number = return_tuple.second;
        } catch(const std::runtime_error& e) {
            cout << e.what( ) << endl;
            exit(EXIT_FAILURE);
        
        }

        // Start the timer
        auto start1 = std::chrono::high_resolution_clock::now();
        vector<pair<int,double>> result = find_knn(query_point);

         // Stop the timer
        auto stop1 = std::chrono::high_resolution_clock::now();
        
        auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(stop1 - start1);
        double k_nearest_time = duration1.count();

        // auto start2 = std::chrono::high_resolution_clock::now();
        // vector<double> bf_result = bf_nearest_neighbor(query_point,N);
        // auto stop2 = std::chrono::high_resolution_clock::now();
        // auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(stop2 - start2);
        // double bf_time = duration2.count();
         
        set<int> range_result = range_search(query_point,R);    

        if (latent) {
            vector<vector<uint8_t>> real_pointset =  get_pointset("data/train-images.idx3-ubyte");
            vector<uint8_t> real_query_point = get_pointset("data/t10k-images-idx3-ubyte")[query_number];

            vector<pair<int,double>> real_nn;
            for(const pair<int,double>& neighbour : result) {
                real_nn.push_back(pair<int,double>(neighbour.first,metric(real_query_point,real_pointset[neighbour.first])));
            }

            auto start2 = std::chrono::high_resolution_clock::now();
            double af = 0;
            assert(real_nn.size()>=1);
            af = real_nn[0].second;
            vector<double> bf_result = new_bf_nearest_neighbor(real_query_point,real_pointset,N);
            // vector<double> bf_result = bf_nearest_neighbor(query_point,N);
            auto stop2 = std::chrono::high_resolution_clock::now();
            auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(stop2 - start2);
            double bf_time = duration2.count();
            af /= bf_result[0];
            
            print_output(real_nn,bf_result,range_result,output_file,k_nearest_time,bf_time,query_number,af);
            query_file_path = "";
        } else {
            auto start2 = std::chrono::high_resolution_clock::now();
            vector<double> bf_result = bf_nearest_neighbor(query_point,N);
            auto stop2 = std::chrono::high_resolution_clock::now();
            auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(stop2 - start2);
            double bf_time = duration2.count();
            print_output(result,bf_result,range_result,output_file,k_nearest_time,bf_time,query_number);
            query_file_path = "";
        }
}

vector<double> Hypercube::new_bf_nearest_neighbor(const vector<uint8_t>& query,const vector<vector<uint8_t>>& real_points ,const int n) {

    vector<double> return_vector;

    for(int i = 0; i < n; i++) {
        return_vector.push_back(numeric_limits<int>::max());
    }

    class compare_pair_by_distance {
        public:
        bool operator()(const double& a, const double& b) const {
            return a < b;
        }
    };
    int count = 0;
    for(const vector<uint8_t>& point : real_points) {
        assert(point.size() == 784);
        assert(real_points.size() == 60000);

        double distance = metric(query,point);

        // Compare this distance to the distance in the last place
        // and then sort, so at the last position the point that 
        // is most distant is replaced with some point that is closer
        if (distance != 0 && distance < return_vector[n-1]) {
            return_vector[n-1] = distance;
            sort(return_vector.begin(), return_vector.end(), compare_pair_by_distance{});
        }
        count++;
    }
    assert(count == 60000);

    assert(real_points.size() == 60000);
    assert(query.size() == 784);
    assert(real_points[0].size() == 784);   

    return return_vector;
}

set<int> Hypercube::range_search(const vector<uint8_t>& query_point,const int r) {

    set<int> selected_images;

    vector<t_point> temp_close_points;
    uint compared_points = 0;
    int* temp_evaluations = new int[k];
    _map_point_to_f_values(query_point, temp_evaluations);
    int starting_hc_vertex = _evaluate(temp_evaluations);

    vector<vector<t_point>> neighboring_vertices_buckets = find_neighboring_vertices_buckets(starting_hc_vertex);

    for(const vector<t_point>& hc_vertex : neighboring_vertices_buckets ) {
        for(const t_point& candidate_neighbour : hc_vertex) {

            double distance = metric(query_point, candidate_neighbour.get_point());
            if (distance == 0 ) continue;
            // Compare this distance to the distance in the last place
            // and then sort, so at the last position the point that 
            // is most distant is replaced with some point that is closer
            if (distance <= r) selected_images.insert(candidate_neighbour.get_img_id());
            
            // Check weather we have surpassed the threshold for points to be searched
            if (M > 0 && ++compared_points >= M) break;
        }
    }

    delete temp_evaluations;
    return selected_images;
}



vector<double> Hypercube::bf_nearest_neighbor(const vector<uint8_t>& query ,const int n) {

    vector<double> return_vector;

    for(int i = 0; i < n; i++) {
        return_vector.push_back(numeric_limits<int>::max());
    }

    class compare_pair_by_distance {
        public:
        bool operator()(const double& a, const double& b) const {
            return a < b;
        }
    };

    for(const vector<uint8_t>& point : training_set) {

        double distance = metric(query,point);

        // Compare this distance to the distance in the last place
        // and then sort, so at the last position the point that 
        // is most distant is replaced with some point that is closer
        if (distance != 0 && distance < return_vector[n-1]) {
            return_vector[n-1] = distance;
            sort(return_vector.begin(), return_vector.end(), compare_pair_by_distance{});
        }
    }

    return return_vector;
}

void Hypercube::print_output(vector<pair<int,double>>& kn_result,vector<double>& bf_result,set<int>& range_results,const string output_file,
const double kn_time ,const double bf_time,const int query_number,const double af) const {
    std::ofstream outputFile(output_file);
    if(!outputFile.is_open()) {
        throw std::runtime_error("Error opening output file\n");
    }
    outputFile << "Query: " << query_number << endl;
    for(int i = 0; i < N; i++) {
        outputFile << "Nearest neighbor-" << i+1 << ": " << kn_result[i].first << endl;
        outputFile << "distanceHypercube: " << kn_result[i].second << endl;
        outputFile << "distanceTrue: " << bf_result[i] << endl;
        if (i == 0) {
            cout << "[";
            long size = training_set[kn_result[i].first].size();
            for(long j = 0; j < size; j++) {
                int numeric_value = static_cast<int>(training_set[kn_result[i].first][j]);
                cout << numeric_value;
                // cout << training_points[kn_result[i].first][j];
                if (j != size-1) {
                    cout << ",";
                }
            }
            cout << "]" << endl;
        }
    }  
    outputFile << "tHypercube: " << kn_time << "ms" << endl;
    outputFile << "tTrue: " << bf_time << "ms" << endl;
    if (latent) outputFile << "AF: " << af << endl;
    outputFile << "R-near neighbors: " << range_results.size() << endl;
    for(int i : range_results) {
        outputFile << i << endl;
    }
}

void Hypercube::_map_point_to_f_values(const vector<uint8_t>& p, int*& temp_evaluations ) {
                
    // For each of the k (=d') h functions (and thus f functions)
    // evaluate the hash and then toss coin to determine
    // binary values, if not evaluated before
    for(int i = 0; i < k; i++) {
        int evaluation = h_functions[i].evaluate(p);
        temp_evaluations[i] = evaluation;
        // Search for the value in the corresponding map
        // if not found, a value for that hash value
        // has not been determined yet, so we have to toss a coin
        std::map<int, int>::iterator it = fs[i].find(evaluation);
        if (it == fs[i].end()) {
            fs[i][evaluation] = get_int_uniform_distribution(0,1);
        }
    }
}