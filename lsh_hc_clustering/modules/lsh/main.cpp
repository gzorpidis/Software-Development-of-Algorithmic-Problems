#include "headers.hpp"
#include "get_long_opts_single.hpp"


int main(int argc, char** argv) {

    LSH::lsh_input input(argc, argv);
    LSH lsh(input, euclideanDistance,500, 1000);

    string response = "";

    while(1) {
        lsh.search();

        cout << "Type (YES) to terminate, or anything else to continue with next query" << endl;
        cin >> response;
        if ((response == "YES" || response == "y" || response == "Yes" || response == "Y" || response == "yes")) { break; }
    };

    return 0;

}