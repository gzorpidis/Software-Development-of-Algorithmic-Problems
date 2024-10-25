#include "headers.hpp"
#include "get_long_opts_single.hpp"

int main(int argc, char** argv) {

    Hypercube::hypercube_input input(argc, argv);
    Hypercube hypercube(input, euclideanDistance, 4);
    string response = "";


    while(1) {
        hypercube.search();
        cout << "Type (YES) to terminate, or anything else to continue with next query" << endl;
        cin >> response;
        if (response == "YES" || response == "y" || response == "Yes" || response == "Y" || response == "yes") { break; }
    };
    
}