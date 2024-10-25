#include "headers.hpp"

g_hash::g_hash() {

};

g_hash::g_hash(const vector<t_hash>& family, int k, int family_size, int table_size)
 : K(k), family_size(family_size), table_size(table_size), M( pow(2,32)-5 ) {
    // Select h_{j} used by g[i]
    // amplified_h[i] = std::set<int>();

    assert(h_indices.empty());
    while (h_indices.size() != K) {
        int chosen = get_int_uniform_distribution(0, family_size-1);
        h_indices.insert(chosen);
    }
    
    for(const int& index : h_indices) {
        // l_hashes[index]
        t_hash temp = family[index];
        // std::cout << "Chosen h_index: " << index << " - " << temp.get_index();
        assert(index == temp.get_index());
        l_hashes.push_back(temp);
        int result = get_int_uniform_distribution(1,5); 
        r.push_back(result);
        // std::cout << " and r: " << result  << std::endl;
    }

}

t_point g_hash::evaluate(const std::vector<uint8_t>& point,const int image_id) {
    
    // Calculate the point ID
    // As the summation of r_i*h_i
    // from i = 1 to i = K
    int id = 0;

    for(long unsigned i = 0; i < K; i++) {
        int r_i = r[i];
        int h_i = l_hashes[i].evaluate(point);
        // cout << r_i << " " << h_i << endl;
        id += (r_i * h_i) % M ;
    }
    
    // Final ID is summation % M
    id = id % M;

    // The bucket it the ID % TableSize
    int g = id % table_size;
    struct t_point to_insert(id, g, point, image_id);
    
    // if point is querry point, don't insert it to bucket
    if(image_id < 0) return to_insert;
    buckets[g].push_back(to_insert);
    
    return to_insert;
}

vector<t_point> g_hash::find_in_bucket(const t_point& point,const int bucket_number ) {
    
    vector<t_point> return_points;

    vector<struct t_point> hash_bucket = buckets[bucket_number];
    
    for(t_point& p : hash_bucket) {
        if (p.get_id() == point.get_id()) {
            return_points.push_back(p);
        }
    }
    return return_points;
}