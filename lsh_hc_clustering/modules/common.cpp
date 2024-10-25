#include "headers.hpp"

int get_int_uniform_distribution(int min = 0, int max = std::numeric_limits<int>::max() ) {

    std::random_device rd;
    std::mt19937 gen(rd());
    
    std::uniform_int_distribution<int> distribution(min, max);
    
    return distribution(gen);
}

double get_double_uniform_distribution(double min = 0, double max = 1.0) {
    std::random_device rd;
    std::mt19937 gen(rd());
    
    std::uniform_int_distribution<int> distribution(min, max);

    return distribution(gen);
}

double gen_number_normal_distribution() {

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<double> distribution(0.0, 1.0);
    return distribution(gen);

}

void fill_random_normal_vector(std::vector<double>& vec, const int length) {
    
    if (!vec.empty()) {
        throw std::invalid_argument("Argugment vector should be empty");
    }

    if (length <= 0) {
        throw std::invalid_argument("Length should be positive");
    }
    
    
    for(int i = 0 ; i < length; i++) {
        vec.push_back(gen_number_normal_distribution());
    }

}

std::vector<double> generate_random_normal_vector(const int length) {
    std::vector<double> temp;
    fill_random_normal_vector(temp, length);
    return temp;
}

double euclideanDistance(const vector<uint8_t>& vector1, const vector<uint8_t>& vector2) {
    // Check if the vectors have the same dimension
    if (vector1.size() != vector2.size()) {
        throw std::invalid_argument("Vectors must have the same dimension");
    }

    double distance = 0.0;
    
    for (size_t i = 0; i < vector1.size(); i++) {
        // Calculate the squared difference between corresponding elements
        double diff = static_cast<double>(abs(vector1[i] - vector2[i]));
        distance += diff * diff;
    }
    
    // Take the square root of the sum of squared differences
    return std::sqrt(distance);
}

// Function to calculate the Hamming distance between two integers
int hamming_distance(int x, int y) {
    // Calculate the bitwise XOR of the two integers to find differing bits
    int xorResult = x ^ y;
    int distance = 0; // Initialize the Hamming distance to 0

    // Count the set bits in the XOR result to determine the Hamming distance
    while (xorResult > 0) {
        distance += xorResult & 1; // Increment distance if the lowest bit is set (1)
        xorResult >>= 1; // Right-shift to check the next bit
    }
    return distance; // Return the calculated Hamming distance
}