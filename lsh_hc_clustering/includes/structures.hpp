#include "headers.hpp"
#pragma once
class t_point{
    private:
        // ID(p) = Sum_{r_{i}*h_{i}} mod M
        int _id;
        // A unique image id, initialized when creating a t_point
        int _image_id;
        // On which bucket the point is saved at
        int _bucket;
        // The point itself
        std::vector<uint8_t> _point_vector;

    public:
        t_point(int id, int bucket, const std::vector<uint8_t>& point, int image_id) : _id(id),_image_id(image_id),_bucket(bucket), _point_vector(point)   {};
        t_point(int bucket, const std::vector<uint8_t>& point, int image_id) : _id(-1), _image_id(image_id), _bucket(bucket),_point_vector(point) {};
        t_point() {};

        int get_id() const { return _id; };
        int get_bucket() const { return _bucket; };
        int get_img_id() const { return _image_id; };
        std::vector<uint8_t> get_point() const { return _point_vector; };
        bool operator<(const t_point& other ) const {
            return _image_id < other._image_id;
        }
};


class t_hash {
    std::vector<double> v;
    int window;
    float noise;
    int r;
    int pos;

    double inner_product(const vector<uint8_t>&, const vector<double>&);

    public:
        int get_index() const;

        // h(p) evaluation
        /*
        *   @brief 
        *       Given a point p, evaluate its h(p) value
        *   
        *   @param point A 1D vector of binary values, representing a flattened image
        *   @return Its evaluation number, representing a bucket
        */
        int evaluate(const std::vector<uint8_t>& point);

        // Initiliaze the t_hash structure, passing the window
        // will initiliaze a random v vector
        // and noise t, for later evaluation
        t_hash(int, int, int);

        t_hash(const t_hash&);

        t_hash();

};

class g_hash {

    // Set of indeces picked from the H family to construct g
    std::set<int> h_indices;

    std::vector<int> r;
    // Set of actuall h hash functions (t_hash)
    std::vector<t_hash> l_hashes;
    // Actual hash map

    

    long unsigned int K;
    int family_size;
    int table_size;
    unsigned int M;
    
    public:
    std::map<int, std::vector<t_point>> buckets;

    public:
        g_hash();
        g_hash(const std::vector<t_hash>& ,int, int , int);

        /*
        *   @brief 
        *       Given a point p, and its ID number,
        *       evaluate its hash table position 
        *   
        *   @param point A 1D vector of binary values, representing a flattened image
        *   @param image_id A unique image identifier
        *   @return A t_point value representing this point to the fullest in the hash table
        *   saving its image id (provided as parameter), it hashing ID (sum of h_j),
        *   the point itself, and the bucket in which it was inserted in
        *   
        */
        t_point evaluate(const std::vector<uint8_t>& point,const int image_id);

        /*
        *   @brief 
        *       Given a point p, and its ID number,
        *       evaluate its hash table position 
        *   
        *   @param point A 1D vector of binary values, representing a flattened image as a query image
        *   @param bucket_number The bucket in which it was inserted
        *   @return A vector of points (t_point type) that are close to the given point
        *   , according to the ID of the bucket in which they hashed
        */
        vector<t_point> find_in_bucket(const t_point& point, const int bucket_number );
};
