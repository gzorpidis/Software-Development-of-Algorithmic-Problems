
#include "graph_search.hpp"
#include <string>

void Graph_Search_Ann::set_inputs(graph_search_input& inputs) {
            
    std::string dataset_path;

    /* if -d is not given in options, ask for it*/
    if (inputs.mapped_values["d"] == "") {
        std::cout << "Give path to dataset from graph_search " << std::endl;
        std::cin >> dataset_path;
        inputs.mapped_values["d"] = dataset_path;
    }
   
    this->query_file_path = inputs.mapped_values["q"];
    this->output_file = inputs.mapped_values["o"];
    this->input_file = inputs.mapped_values["d"];
    this->graph_file_path = inputs.mapped_values["g"];

    if (inputs.mapped_values["k"] == "") { inputs.mapped_values["k"] = "50"; }
    if (inputs.mapped_values["E"] == "") { inputs.mapped_values["E"] = "30"; }
    if (inputs.mapped_values["R"] == "") { inputs.mapped_values["R"] = "1"; }
    if (inputs.mapped_values["N"] == "") { inputs.mapped_values["N"] = "1"; }
    if (inputs.mapped_values["l"] == "") { inputs.mapped_values["l"] = "20";}
    if (inputs.mapped_values["m"] == "") { inputs.mapped_values["m"] = "1"; }


    this->k =  stoi(inputs.mapped_values["k"]); 
    this->E =  stoi(inputs.mapped_values["E"]); 
    this->R =  stoi(inputs.mapped_values["R"]);
    this->N =  stoi(inputs.mapped_values["N"]);
    this->l =  stoi(inputs.mapped_values["l"]);
    this->m =  stoi(inputs.mapped_values["m"]);
    if (m > 2 || m <= 0) {
        std::cout << "Invalid algorithm, please provide 1 for GNN or 2 for MRNG" << std::endl;
        exit(EXIT_FAILURE);
    }

};


Graph_Search_Ann::Graph_Search_Ann(graph_search_input& input_file_, std::function<double(const std::vector<uint8_t>&, const std::vector<uint8_t>&)> metric_)
    : number_of_vector_hash_functions(4), to_print(nullptr), T(50) ,lsh(nullptr), hc(nullptr),  metric(metric_){
    
    set_inputs(input_file_);
    to_print = new Results(output_file,m);
    cout << "Input file: " << input_file << endl;
    {
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
        
        lsh = new LSH(ip, metric,100, 4);
        cout << "LSH created" << endl;
    }

    training_points = get_pointset(input_file);
    
    int id = 0;
    for(std::vector<uint8_t>& p : training_points) {
        points_inside.emplace_back(t_point(-1, p, id++));
    }

    if (graph_file_path != "") {
        cout << "Reading from graph file..." << endl;

        ifstream inputFile(graph_file_path);

        if (!inputFile.is_open()) {
            cerr << "Error opening the file: " << graph_file_path << endl;
            return;
        }

        std::string line;
        while (std::getline(inputFile, line)) {
            std::istringstream iss(line);

            int node;
            char colon;
            iss >> node >> colon;

            if (colon != ':') {
                std::cerr << "Error parsing line: " << line << " (colon expected)" << std::endl;
                continue;
            }
            search_graph.add_vertex(node);
            int neighbor;   
            size_t counter = 0;
            while (iss >> neighbor) {
                counter++;
                // graph[node].push_back(neighbor);
                search_graph.add_edge(node, neighbor);
                // Check for the comma after each neighbor
                char comma;
                if ( m == 1 && ( !(iss >> comma) || (comma != ',' && !isspace(comma)) )) {
                    break;
                }
            }

            assert(counter == search_graph.get_connected_to(node).size());
        }
        assert(search_graph.size() == 60000);
    }

    // If no graph_file_path was given, construct graph from start
    if (graph_file_path == "" && m == 1) {
        create_knn_graph("outputs/gnn_graph.txt");
    } else if (graph_file_path == "" && m == 2) {
        mrng_construction("outputs/mrng_graph.txt");
    }
};

class checked {
    int id;
    bool flagged;
    double distance;

    public:
        checked(int id_, double distance_, bool flagged_ = false) {
            id = id_;
            flagged = flagged_;
            distance = distance_;
        };

        void flag() { flagged = true; };
        void unflag() { flagged = false; };
        bool get_flag() const { return flagged; };
        int get_id() const { return id; };
        double get_distance() const { return distance; }
};

bool compareCheckedByDistance(const checked& a, const checked& b) {
    // Compare by distance in ascending order.
    return a.get_distance() < b.get_distance();
}

void Graph_Search_Ann::graph_search(const std::vector<uint8_t>& query_point) {
    
    auto start_time = std::chrono::high_resolution_clock::now();
    std::vector<checked> R;

    int selected = get_int_uniform_distribution(0, training_points.size()-1);
    const t_point& p = points_inside[selected];
    checked t(selected, metric(query_point, p.get_point()));
    R.push_back(t);

    int i = 1;
    while(i < l) {
        // cout << i << endl;
        int id = -1;
        
        for(checked& node : R) {
            if (!node.get_flag()) {
                node.flag();
                id = node.get_id();
                break;
            }
        }

        assert(id != -1);
        std::vector<int> neighbours = search_graph.get_connected_to(id);

        for(const int& neighbor : neighbours) {
            bool found = false;
            for(const checked& point : R) {
                if (point.get_id() == neighbor ) {
                    found = true;
                }
            }
            
            if (!found) {
                const t_point& N = points_inside[neighbor];
                checked t(neighbor, metric(query_point, N.get_point()));
                R.push_back(t);
                i++;
            }
        }
        
        std::sort(R.begin(), R.end(), compareCheckedByDistance);
        
    }

    std::vector< std::pair<int,double>> nearest_neighbors;

    int limit = -1;
    if (R.size() < N) {
        limit = R.size();
    } else {
        limit = N;
    }
    assert(limit >= 0);

    for(int i = 0; i < limit; i++) {
        nearest_neighbors.push_back(make_pair(R[i].get_id(), R[i].get_distance()));
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    double tApproximate = duration.count();

    
    // auto start_time_ = std::chrono::high_resolution_clock::now();
    // std::vector<double> true_nearest = lsh->bf_nearest_neighbor(query_point, N);
    // auto end_time_ = std::chrono::high_resolution_clock::now();
    // auto duration_ = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_ - start_time_);
    // double tTrue = duration_.count();

    get_real_points("data/train-images.idx3-ubyte", "data/t10k-images-idx3-ubyte", nearest_neighbors);
    
    // double af = nearest_neighbors[0].second/true_nearest[0];
    // to_print->change_maf(af);
    to_print->t_aproximates.push_back(tApproximate);
    // to_print->aaf.push_back(af);
    // to_print->t_true.push_back(tTrue);

}

void Graph_Search_Ann::search() {
    
    do{
        if(query_file_path == "") {
            std::cout << "Give path to query file from graph_search" << std::endl;
            std::cin >> query_file_path;
        }
            std::pair<std::vector<uint8_t>,int> query = get_query_point(query_file_path);
            to_print->image_number = query.second;
            std::vector<uint8_t> query_point = query.first;
        if (m == 1) {
            gnns_search(query_point);
        } else {
            graph_search(query_point);
        }
        
        // Print query results
        to_print->print_query(N);
        std::cout << "Type y/Y to terminate or give path to new query file" << endl;
        std::cin >> query_file_path;
   } while (query_file_path != "y" && query_file_path != "Y" && query_file_path != "yes" && query_file_path != "Yes");
   to_print->print_avarage_results();

}

void Graph_Search_Ann::gnns_search(const std::vector<uint8_t>& query_point) {
    
    auto start_time = std::chrono::high_resolution_clock::now();
    
    assert(search_graph.size() != 0);
    if (search_graph.size() == 0) {
        throw std::runtime_error("Graph has not been initialized. Use -g to load a graph from a text file, as a csv file.");
    }

    // Initialize S
    std::set<int> neighbours_extended_id;
    std::priority_queue<t_point, std::vector<t_point>, CompareDistance> pq(CompareDistance(metric, query_point));

    // After the graph is created we can search in it:
    
    // for r = 1, ... R do
    for(int r = 0; r < R; r++) {
        // Y_0: a random point uniformly over D
        int selected = get_int_uniform_distribution(0, training_points.size()-1);
        
        // for t = 1, ... T, do:
        for(int t = 0; t < T;t++) {
            
            neighbours_extended_id.insert(selected);
            
            t_point &selected_point = points_inside[selected];
            double distance_to_query = metric(query_point, selected_point.get_point());
            
            // std::vector<int> neighbours = graph[selected];
            std::vector<int> neighbours = search_graph.get_connected_to(selected);
            int limit = E;
            if(E > neighbours.size()) limit = neighbours.size();
            std::vector<int> e_neigbours(neighbours.begin(), neighbours.begin() + limit);

            double min_distance = std::numeric_limits<double>::infinity();
            int closest = -1;
            int count = 0;
            bool selected_bigger = false;
            for(const int& neighbour : e_neigbours) {
                count++;
                neighbours_extended_id.insert(neighbour);
                t_point& neigh = points_inside[neighbour];
                double distance = metric(query_point, neigh.get_point());
                if (distance < min_distance) {
                    closest = neighbour;
                    min_distance = distance;
                }

                if(distance < distance_to_query) {
                    selected_bigger = true;
                }
            }

            if (!selected_bigger) {
                break;
            }

            selected = closest;
            if (selected == -1) { assert(false); }

        }
    }
    

    for(int neighbour : neighbours_extended_id) {
        pq.push(points_inside[neighbour]);
    }

    std::vector<int> results;
    std::vector<std::pair<int, double>> nn;
    int limit = -1;

    if ( N >= pq.size()) { limit = pq.size(); }
    else {
        limit = N;
    }

    for(int i = 0 ; i < limit; i++) {
        t_point popped = pq.top();
        results.push_back(popped.get_img_id());
        double distance =  metric(query_point, popped.get_point());
        nn.push_back(make_pair(popped.get_img_id(), distance));
        pq.pop();
    }

    auto end_time = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    double tApproximate = duration.count();

    
    // auto start_time_ = std::chrono::high_resolution_clock::now();
    // std::vector<double> true_nearest = lsh->bf_nearest_neighbor(query_point, N);
    // auto end_time_ = std::chrono::high_resolution_clock::now();
    // auto duration_ = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_ - start_time_);
    // double tTrue = duration_.count();

    get_real_points("data/train-images.idx3-ubyte", "data/t10k-images-idx3-ubyte", nn);

    to_print->t_aproximates.push_back(tApproximate);
}

void Graph_Search_Ann::create_knn_graph(const std::string filename) {
    std::ofstream outputFile(filename);
    // Check if the file is successfully opened
    if (outputFile.is_open()) {
        assert(true);
    } else assert(false);

    
    for(t_point& t_p : points_inside) {
        // t_point graph_vertex_point(-1, p, id);
        // id++;
        // Add to the graph
        // cout << "at point: " << t_p.get_img_id() << endl;

        if (filename != "")
            outputFile << t_p.get_img_id() << " : ";
        
        search_graph.add_vertex(t_p.get_img_id());
        std::vector<uint8_t> p = t_p.get_point();
        std::vector<std::pair<int,double>> results = lsh->find_k_nearest_neighbour(p, k);

        for(size_t i = 0; i < results.size(); i++) {
            
            if (results[i].second == numeric_limits<int>::max() ) {
                continue;
            }
            const t_point& neighbour = points_inside[results[i].first];
            search_graph.add_vertex(neighbour.get_img_id());
            search_graph.add_edge(t_p.get_img_id(), neighbour.get_img_id());
            if (filename != "")
                outputFile << " " << neighbour.get_img_id() << "," ;
            
            // results[i].first
            // t_point neigbour(results[i].first, training_points[results[i].first], -1);
            // search_graph.add_edge(graph_vertex_point, neigbour);
            assert(search_graph.size() <= training_points.size());

        }
        if (filename != "")
            outputFile << endl;
    }
    outputFile.close();
}

void Graph_Search_Ann::mrng_construction([[maybe_unused]] const std::string filename) {
    
    cout << "Creating MRNG graph" << endl;
    
    for(const t_point& t_p : points_inside) {
        search_graph.add_vertex(t_p.get_img_id());
    }
    
    std::set<int> S;
    for(size_t i = 0; i < training_points.size();i++) { S.insert(i); }

    const int num_threads = 4;

    std::vector<std::thread> threads(num_threads);
    std::vector<std::thread::id> thread_ids;

    // Create a mutex to protect access to shared resources
    std::mutex mutex;
    std::mutex threadIdsMutex;  // Mutex to protect access to threadIds

    // Function to be executed by each thread
    auto process_range = [&](int start, int end) {
        size_t threadIdHash = std::hash<std::thread::id>{}(std::this_thread::get_id());
        std::string threadIdStr = std::to_string(threadIdHash);
        std::ofstream thread_output_file("outputs/mrng_graph_" + threadIdStr + ".txt");

        {
            std::lock_guard<std::mutex> lock(threadIdsMutex); // Lock to protect threadIds
            cout << "Created thread with ID " << std::this_thread::get_id() << " using range(" << start << "," << end << ")" << endl;   
        }

        for (int i = start; i <= end; ++i) {
            const t_point& p = points_inside[i];
            thread_output_file << p.get_img_id() << " : ";

            std::set<int> R_p = S;
            R_p.erase(p.get_img_id());

            std::set<int> L;
            vector<bool> in_lp = vector<bool>(training_points.size(), false);
            std::vector<uint8_t> vector = p.get_point();

            int closest_first = lsh->find_nearest_neighbour(vector).first;
            L.insert(closest_first);
            in_lp[closest_first] = true;
            R_p.erase(closest_first);
            

            for (const int& r : R_p) {
                bool condition = true;

                for (const int& t : L) {

                    // If r is not in L
                    if (!in_lp[r]) {
                        t_point& r_vertex = points_inside[r];
                        t_point& t_vertex = points_inside[t];
                        double edge_p_r = metric(vector, r_vertex.get_point());
                        double edge_p_t = metric(vector, t_vertex.get_point());

                        if (edge_p_r > edge_p_t) {
                            double edge_r_t = metric(r_vertex.get_point(), t_vertex.get_point());
                            if (edge_p_r > edge_r_t) {
                                condition = false;
                                break;
                            }
                        }

                    }
                }

                if (condition) {
                    L.insert(r);
                    in_lp[r] = true;
                    thread_output_file << r << " ";

                    {
                        std::lock_guard<std::mutex> lock(mutex); // Lock to protect shared resources
                        search_graph.add_edge(p.get_img_id(), r);
                    }
                }

            
            }
            thread_output_file << std::endl;
        }
        {
            std::lock_guard<std::mutex> lock(threadIdsMutex); // Lock to protect threadIds
            thread_ids.push_back(std::this_thread::get_id());
        }
        thread_output_file.close();
    };

    int points_per_thread = points_inside.size() / num_threads;
    size_t total_items = 0;

    for (size_t i = 0; i < num_threads; ++i) {
        int start = i * points_per_thread;
        int end = (i == num_threads - 1) ? points_inside.size() - 1 : start + points_per_thread - 1;
        total_items += (end - start + 1);
        cout << "Thread : " << i << " starts at: " << start << " ends at << " << end << endl;
    }

    if (total_items != points_inside.size()) { std::cerr << "Total points to be calculated by threads does not equal the whole dataset size."; return; }
    assert(total_items == points_inside.size());
    
    for(size_t i = 0; i < num_threads; i++) {
        int start = i * points_per_thread;
        int end = (i == num_threads - 1) ? points_inside.size() - 1 : start + points_per_thread - 1;
        threads[i] = std::thread(process_range, start, end);
    }


    cout << "Waiting for threads..." << endl;
    // Wait for all threads to finish
    for (int i = 0; i < num_threads; ++i) {
        threads[i].join();
    }

    for (const auto& id : thread_ids) {
        std::cout << "Thread ID: " << id << std::endl;
    }
}

vector<double> Graph_Search_Ann::new_bf_nearest_neighbor(const vector<uint8_t>& query,const vector<vector<uint8_t>>& real_points ,const int n) {

    vector<double> return_vector;

    for(int i = 0; i < n; i++) {
        return_vector.push_back(numeric_limits<int>::max());
    }
    assert(n == 1);

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

void Graph_Search_Ann::get_real_points(const std::string dataset, const std:: string query_set, std::vector<std::pair<int, double>>& nn_10) {

    vector<vector<uint8_t>> real_pointset = get_pointset(dataset);
    assert(real_pointset.size() == 60000);
    vector<vector<uint8_t>> real_query_pointset = get_pointset(query_set);
    assert(real_query_pointset.size() == 10000);
    vector<uint8_t> real_query_point = real_query_pointset[to_print->image_number];
    assert(real_query_point.size() == 784);
    vector<pair<int,double>> r_nn;

    for(const pair<int,double>& p : nn_10) {
        double real_distance = metric(real_query_point, real_pointset[p.first]);
        r_nn.push_back(pair<int,double>(p.first, real_distance));
    }

    auto start_time_ = std::chrono::high_resolution_clock::now();
    to_print->nn_true = new_bf_nearest_neighbor(real_query_point,real_pointset,N);
    auto end_time_ = std::chrono::high_resolution_clock::now();
    auto duration_ = std::chrono::duration_cast<std::chrono::milliseconds>(end_time_ - start_time_);
    double tTrue = duration_.count();


    to_print->t_true.push_back(tTrue);
    to_print->nn = r_nn;
    to_print->aaf.push_back(r_nn[0].second/to_print->nn_true[0]);
    to_print->change_maf(r_nn[0].second/to_print->nn_true[0]);
}