#include "headers.hpp"
#include <algorithm>

LSH::LSH(lsh_input& input_file,std::function<double(const std::vector<uint8_t>&, const std::vector<uint8_t>&)> metric,
    int family_size,int window_size) :   
        metric(metric),family_size(family_size), W(window_size)
        {

            set_inputs(input_file);
            
            try {
               training_points = get_pointset(data_file_path);   
            } catch( std::runtime_error& e ) {
                std::cout << e.what() << std::endl;
            }

            int n = training_points.size();
            int vector_size = training_points[0].size();
            for(int i = 0; i < family_size; i++) {
                hash_family.push_back(t_hash(vector_size, W, i));
            }

            // cout << "Initialized L to : " << L << endl;
            for(int i = 0; i < L ;i++) {
                // Make L gs
                g_hash_functions.push_back(g_hash(hash_family, k, family_size, n/4));
                // g has collected all the indexes of h
            }

            // cout << "G hash functions are: " <<  g_hash_functions.size() << endl;
            
            _training();
            cout << "Training done" << endl;
        }


void LSH::_training() {
    int id = 0;
    // cout << "Training starts" << endl;
    // assert(training_points.size() == 60000);
    
    for(long unsigned int i = 0; i < training_points.size();i++) {
        for(g_hash& g : g_hash_functions) {
            g.evaluate(training_points[i],id);
        }
        id++;
    }

    for(const g_hash& g: g_hash_functions) {
        // int id = 0;
        for (const auto& pair : g.buckets) {
            vector<t_point> points = pair.second;
            points.size();
        }
        // assert(id == 60000);
    }

}

std::pair<int,double> LSH::find_nearest_neighbour(vector<uint8_t>& query_point) {

    double min_dist = std::numeric_limits<double>::infinity();
    int id;
    double dist;
    t_point closest_point;
    
    vector<t_point> temp_close_points;
    set<t_point> set_temp_close_points;
    for(g_hash& g : g_hash_functions) {
        t_point point_struct = g.evaluate(query_point,-1);
        temp_close_points = g.find_in_bucket(point_struct, point_struct.get_bucket() );    
        for(const t_point& p : temp_close_points) {
            if(set_temp_close_points.find(p) == set_temp_close_points.end()) {
                set_temp_close_points.insert(p);
            } else {
                continue;
            }
            if ( (dist = metric(query_point, p.get_point() ) )  < min_dist ) {
                if(dist != 0) {
                    min_dist = dist;
                    id = p.get_img_id();
                }
            }
        }   
    }
    return std::pair<int, double>(id, min_dist);
}

vector<pair<int,double>> LSH::find_k_nearest_neighbour(vector<uint8_t>& query_point,const int n) {
    
    if (k == 1) {
        return vector<pair<int,double>>(1,find_nearest_neighbour(query_point));
    }

    set<int> selected_images;

    vector<pair<int,double>> return_vector;
    
    for(int i = 0; i < n; i++) {
        pair<int,double> p = make_pair(0,numeric_limits<double>::max());
        return_vector.push_back(p);
    }

    class compare_pair_by_distance {
        public:
        bool operator()(const std::pair<int, double>& a, const std::pair<int, double>& b) const {
            return a.second < b.second;
        }
    };

    set<t_point> set_temp_close_points;
    vector<t_point> temp_close_points;

    for(g_hash& g : g_hash_functions) {
        t_point point_struct = g.evaluate(query_point,-1);
        temp_close_points = g.find_in_bucket(point_struct, point_struct.get_bucket());    
        for(const t_point& p : temp_close_points) {
            if(set_temp_close_points.find(p) == set_temp_close_points.end()) {
                set_temp_close_points.insert(p);
            } else {
                continue;
            }
            double dist = metric(query_point, p.get_point() );
            if (dist != 0 && dist < return_vector[n-1].second) {
                return_vector.push_back(pair<int,double>(p.get_img_id(),dist));
                int size = p.get_point().size();
                std::vector<uint8_t> t = p.get_point();
                for(int i = 0; i < size; i++) {
                    assert(training_points[p.get_img_id()][i] == t[i]);
                }
                sort(return_vector.begin(), return_vector.end(),compare_pair_by_distance{});
            }
        }   
    }

    return return_vector;
} 

vector<double> LSH::new_bf_nearest_neighbor(const vector<uint8_t>& query,const vector<vector<uint8_t>>& real_points ,const int n) {

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

void LSH::search() {
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
    vector<pair<int,double>> result = find_k_nearest_neighbour(query_point,N);
    // Stop the timer
    auto stop1 = std::chrono::high_resolution_clock::now();
    
    auto duration1 = std::chrono::duration_cast<std::chrono::milliseconds>(stop1 - start1);
    double k_nearest_time = duration1.count();
    
    set<int> range_result = range_search(query_point,R);    
    
    if (latent) {

        vector<vector<uint8_t>> real_pointset =  get_pointset(init_data_file_path);
        vector<uint8_t> real_query_point = get_pointset(init_query_file_path)[query_number];

        vector<pair<int,double>> real_nn;
        int counter = 0;
        double af = 0;
        for(const pair<int,double>& neighbour : result) {
            real_nn.push_back(pair<int,double>(neighbour.first,metric(real_query_point,real_pointset[neighbour.first])));
            if (counter == 0) {
                assert(real_nn.size()>=1);
                af = real_nn[0].second;
            }
            counter ++;
        }

        auto start2 = std::chrono::high_resolution_clock::now();

        vector<double> bf_result = new_bf_nearest_neighbor(real_query_point,real_pointset,N);
        // vector<double> bf_result = bf_nearest_neighbor(query_point,N);
        auto stop2 = std::chrono::high_resolution_clock::now();
        af /= bf_result[0];
        auto duration2 = std::chrono::duration_cast<std::chrono::milliseconds>(stop2 - start2);
        double bf_time = duration2.count();
        print_output(real_nn,bf_result,range_result,output_file,k_nearest_time,bf_time,query_number, af);
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

void LSH::set_inputs(lsh_input& inputs) {
    std::string dataset_path;


    /* if -d is not given in options, ask for it*/
    if (inputs.mapped_values["d"] == "") {
        std::cout << "Give path to dataset " << endl;
        std::cin >> dataset_path;
        inputs.mapped_values["d"] = dataset_path;
    }

    /* Check each initilization option, if not given given default values */
    if (inputs.mapped_values["L"] == "") { inputs.mapped_values["L"] = "5"; }
    if (inputs.mapped_values["k"] == "") { inputs.mapped_values["k"] = "4"; }
    if (inputs.mapped_values["N"] == "") { inputs.mapped_values["N"] = "1"; }
    if (inputs.mapped_values["R"] == "") { inputs.mapped_values["R"] = "10000"; }
    if (inputs.mapped_values["q"] == "") { this->query_file_path = ""; } else { this->query_file_path = inputs.mapped_values["q"]; }

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
    this->output_file = inputs.mapped_values["o"];
    this->data_file_path = inputs.mapped_values["d"];
    
    this->L = std::stoi(inputs.mapped_values["L"]);
    this->k = std::stoi(inputs.mapped_values["k"]);
    this->N = std::stoi(inputs.mapped_values["N"]); 
    this->R = std::stoi(inputs.mapped_values["R"]);
    
}

set<int> LSH::range_search(const vector<uint8_t>& query_point,const int r) {

    if ( r == 0 ) return set<int>();
    set<int>   selected_images;

    vector<t_point> temp_close_points;
    for(g_hash& g : g_hash_functions) {
        t_point point_struct = g.evaluate(query_point,-1);
        temp_close_points = g.find_in_bucket(point_struct, point_struct.get_bucket());    
        for(const t_point& p : temp_close_points) {

            double dist = metric(query_point, p.get_point());
            // cout << dist<< " " << r << endl;
            if (dist != 0 && dist < r) {
                selected_images.insert(p.get_img_id());
            }
        }   
    }

    return selected_images;
}

vector<double> LSH::bf_nearest_neighbor(const vector<uint8_t>& query ,const int n) {

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

    for(const vector<uint8_t>& point : training_points) {

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

void LSH::print_output(vector<pair<int,double>>& kn_result,vector<double>& bf_result,set<int>& range_results,const string output_file,
const double kn_time ,const double bf_time,const int query_number, const double af) const {
    std::ofstream outputFile(output_file);
    if(!outputFile.is_open()) {
        throw std::runtime_error("Error opening output file\n");
    }
    outputFile << "Query: " << query_number << endl;
    for(int i = 0; i < N; i++) {
        outputFile << "Nearest neighbor-" << i+1 << ": " << kn_result[i].first << endl;
        outputFile << "distanceLSH: " << kn_result[i].second  << endl;
        outputFile << "distanceTrue: " << bf_result[i] << endl;
        if (i == 0) {
            cout << "[";
            long size = training_points[kn_result[i].first].size();
            for(long j = 0; j < size; j++) {
                int numeric_value = static_cast<int>(training_points[kn_result[i].first][j]);
                cout << numeric_value;
                if (j != size-1) {
                    cout << ",";
                }
            }
            cout << "]" << endl;
        }
    }  
    outputFile << "tLSH: " << kn_time << "ms" << endl;
    outputFile << "tTrue: " << bf_time << "ms" << endl;
    if (latent)
        outputFile << "AF: " << af << endl;
    outputFile << "R-near neighbors: " << range_results.size() << endl;
    for(int i : range_results) {
        outputFile << i << endl;
    }
}