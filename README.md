# Software Development for hard Algorithmic Problems
##  - 

##### Project implemented by: 

1) Georgios Zorpidis

2) [Konstantinos Fragkos](https://github.com/Konstantinos72002)

##### Overseen by: 

1) (Ioannis Emiris)[https://cgi.di.uoa.gr/~emiris/]

2) (Ioannis Chamodrakas)[http://users.uoa.gr/~ihamod/]

for the University of Athens, Department of Informatics and Telecommunications (2023-2024) 

---

This project explores solutions to Approximate Nearest Neighbors (ANN) clustering using a combination of algorithms and approaches. The main focus is on:

- **Locality Sensitive Hashing (LSH)** algorithms, including a **Hypercube variation**.
- **Graph-based clustering** using **Graph Nearest Neighbor** algorithms with **Monotonic Relative Neighborhood Graphs (MRNG)**. For more details on MRNG, refer to the research [paper](https://www.vldb.org/pvldb/vol12/p461-fu.pdf).

### Project Structure

The project consists of three primary directories:
1. `lsh_hc_clustering` — Implementations of LSH algorithms and its Hypercube variations.
2. `mrng_gnn` — Implementations of MRNG-based clustering algorithms.
3. `cnn_autoencoders` — Basic usage of Convolutional Autoencoders to train a CNN network to produce new points in a latent space (10 and 30 dimensional) and test effectiveness on this latent space.

More info about each sub-project is inside each directory, including our approach, useful comments and code structure.

### Dataset

The dataset used in this project is the well-known **MNIST** dataset, consisting of handwritten digits (0-9).

### Language

The algorithms are implemented from scratch using C++ for time effiency and high control over our code.

---