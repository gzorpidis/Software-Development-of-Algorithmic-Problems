#include "headers.hpp"

double t_hash::inner_product(const std::vector<uint8_t>& vector1, const std::vector<double>& vector2) {
    // cout << vector1.size() << " " << vector2.size() << endl;     
    if (vector1.size() != vector2.size()) {
        throw std::invalid_argument("Vectors must have the same size");
    }
    
    double result = 0.0;
    for (size_t i = 0; i < vector1.size(); ++i) {
        result += vector1[i] * vector2[i];
    }
    return result;
}

int t_hash::get_index() const { return pos; };
        
// h(p) evaluation
int t_hash::evaluate(const std::vector<uint8_t>& p) {
    return static_cast<int>(std::floor(( abs(inner_product(p,v)) + noise) / window));
}

// Initiliaze the t_hash structure, passing the window
// will initiliaze a random v vector
// and noise t, for later evaluation
t_hash::t_hash(int d, int w, int index) : 
window(w), 
noise(get_double_uniform_distribution(0,window)), 
r(get_int_uniform_distribution(-5,5)) ,
pos(index)
{
    fill_random_normal_vector(v, d);
    // std::cout << "Initialized <t_hash>_" << pos << "\n";
}

t_hash::t_hash(const t_hash& copy) : 
v(copy.v),
window(copy.window), 
noise(copy.noise), 
r(copy.r),
pos(copy.pos)
{    
}

t_hash::t_hash() {
}