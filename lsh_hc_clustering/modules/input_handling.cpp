#include "headers.hpp"

static int reverseInt (int i) {
    unsigned char c1, c2, c3, c4;

    c1 = i & 255;
    c2 = (i >> 8) & 255;
    c3 = (i >> 16) & 255;
    c4 = (i >> 24) & 255;

    return ((int)c1 << 24) + ((int)c2 << 16) + ((int)c3 << 8) + c4;
}

pair<vector<uint8_t>,int> get_query_point(string path, string labels_path) {

    std::ifstream file(path, std::ios::binary);
    std::ifstream labels;

    if (labels_path != "") {
        labels = std::ifstream(labels_path, std::ios::binary);
    }

    // std::ifstream labels(labels_path, std::ios::binary);
    int random_point = get_int_uniform_distribution(0,10);
    cout << "Selected: " << random_point << endl;
    if (!file.is_open()) {
        throw runtime_error("Error opening file\n");
    }

    // if (!labels.is_open()) {
    //     throw runtime_error("Error opening labels file\n");
    // }

    int magic_number = 0, test_set_size = 0;
    int mn_labels = 0;

    file.read((char *)&magic_number, sizeof(magic_number));
    // labels.read((char*)&mn_labels, sizeof(mn_labels));

    magic_number = reverseInt(magic_number);
    mn_labels = reverseInt(mn_labels);


    // if(magic_number != 32) throw std::runtime_error("Invalid MNIST image file!");
    // if(mn_labels != 2049)  throw std::runtime_error("Invalid MNIST labels file!");
    file.read((char *)&test_set_size, sizeof(test_set_size));
    test_set_size = reverseInt(test_set_size);

    int n_rows = 0, n_cols = 0;
    file.read((char *)&n_rows, sizeof(n_rows)), n_rows = reverseInt(n_rows);
    file.read((char *)&n_cols, sizeof(n_cols)), n_cols = reverseInt(n_cols);
    // assert(n_rows == 28);
    // assert(n_cols == 28);
    // if (test_set_size != 10000) throw std::runtime_error("Invalid MNIST test set size\n");
    cout << n_rows << " " << n_cols << endl;
    // labels.read((char*)&test_set_size, sizeof(test_set_size));
    test_set_size = reverseInt(test_set_size);
    // if (test_set_size != 10000) throw std::runtime_error("Invalid MNIST test set size\n");

    // cout << "Selected: " << random_point << endl;

    int image_size = n_rows * n_cols;
    std::vector<uint8_t> image;

    if (random_point != 0) {
        uint8_t* temp_buffer = new uint8_t[image_size*random_point];
        file.read((char *)temp_buffer, random_point * image_size);
        // cout << "Read: " << random_point * image_size << " before" << endl;
        delete[] temp_buffer;

        // uint8_t* temp_buffer_t = new uint8_t[random_point];
        // labels.read((char *)temp_buffer_t, sizeof(uint8_t) * random_point);
        // delete[] temp_buffer_t;
    }

    uint8_t* buffer = new uint8_t[image_size];
    // cout << "Now reading next " << image_size << " bytes";
    file.read((char*)buffer, image_size);
    cout << "Query: " << endl;
    for(int j = 0; j < image_size; j++) {
        image.push_back(buffer[j]);
        cout << static_cast<int>(buffer[j]) << " ";
    }
    cout << endl;

    // uint8_t label_assigned;
    // labels.read((char*)&label_assigned, sizeof(label_assigned));
    // uint l = label_assigned;
    // cout << "Selected label is: " << l << endl;
    delete[] buffer;

    // assert(image.size() == (28*28));
    return make_pair(image, random_point);

}

vector<vector<uint8_t>> get_pointset(string path) {

    cout << path << endl;
    std::ifstream file(path, std::ios::binary);

    if (!file.is_open()) {
        std::cerr << "Error opening file.\n";
        throw runtime_error("Error opening file\n");
    }

    int magic_number = 0, n_rows = 0, n_cols = 0;
    int number_of_images = 0;

    file.read((char *)&magic_number, sizeof(magic_number));
    
    magic_number = reverseInt(magic_number);
    
    // if(magic_number != 2051) throw std::runtime_error("Invalid MNIST image file!");
    
    file.read((char *)&number_of_images, sizeof(number_of_images)), number_of_images = reverseInt(number_of_images);
    file.read((char *)&n_rows, sizeof(n_rows)), n_rows = reverseInt(n_rows);
    file.read((char *)&n_cols, sizeof(n_cols)), n_cols = reverseInt(n_cols);
    
    int image_size = n_rows * n_cols;

    std::vector<std::vector<uint8_t>> images;

    for(int i = 0; i < number_of_images; i++) {
        images.push_back(std::vector<uint8_t>());
    }

    uint8_t* buffer = new uint8_t[image_size];

    for(int i = 0; i < number_of_images; i++) {    
        file.read((char *)buffer, image_size);
        for(int j = 0; j < image_size; j++) {
            images[i].push_back(buffer[j]);
        }
    }

    delete[] buffer;
    return images;
}