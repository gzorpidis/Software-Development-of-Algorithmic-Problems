#include "cluster.hpp"

int main(int argc, char** argv) {

    auto start1 = std::chrono::high_resolution_clock::now();
    Cluster::cluster_input input(argc,argv);
    Cluster cluster(input, euclideanDistance);

    cluster.apply_clustering();
    // Stop the timer
    auto stop1 = std::chrono::high_resolution_clock::now();
    auto duration1 = std::chrono::duration_cast<std::chrono::seconds>(stop1 - start1);
    double time = duration1.count();
    cout << "computing silhouettes" << endl;
    // cluster._compute_silhouettes();
    cluster._compute_latent_silhouettes("data/train-images.idx3-ubyte");
    cluster.print_results(time);
    return 0;
}