SRC_DIR := modules
OBJ_DIR := objects
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES)) $(SRC_DIR)/lsh/*.cpp
INCDIR := includes
CFLAGS := -Wall -Wextra -I$(INCDIR)
TARGET := lsh

$(TARGET): $(OBJ_FILES)
	g++ $(CFLAGS) -o $@ $^
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	g++ $(CFLAGS) -c -o $@ $<
SRC_DIR_CUBE := modules/hypercube
OBJ_DIR_CUBE := objects
SRC_FILES_CUBE := $(wildcard $(SRC_DIR_CUBE)/*.cpp)
OBJ_FILES_CUBE := $(patsubst $(SRC_DIR_CUBE)/%.cpp,$(OBJ_DIR_CUBE)/%.o,$(SRC_FILES_CUBE)) $(SRC_DIR)/input_handling.cpp $(SRC_DIR)/common.cpp $(SRC_DIR)/t_hash.cpp
TARGET_CUBE := cube
$(TARGET_CUBE): $(OBJ_FILES_CUBE)
	g++ $(CFLAGS) -o $@ $^
$(OBJ_DIR_CUBE)/%.o: $(SRC_DIR_CUBE)/%.cpp
	g++ $(CFLAGS) -c -o $@ $<

run_lsh:
	./$(TARGET) -d data/d10.bin -q data/q10.bin -o outputs/output.txt  -N 5 -R 3000 -latent Yes -id data/train-images.idx3-ubyte -iq data/t10k-images-idx3-ubyte

run_cube:
	./$(TARGET_CUBE) -d data/train-images.idx3-ubyte -q data/t10k-images-idx3-ubyte -o outputs/output.txt -N 5 -M 2000 -probes 2 -R 3500 -k 5


CLUSTER_SRC_DIR := modules
CLUSTER_OBJ_DIR := objects
CLUSTER_SRC_FILES := $(wildcard $(CLUSTER_SRC_DIR)/*.cpp)
CLUSTER_OBJ_FILES := $(patsubst $(CLUSTER_SRC_DIR)/%.cpp,$(CLUSTER_OBJ_DIR)/%.o,$(CLUSTER_SRC_FILES)) $(CLUSTER_SRC_DIR)/clustering/*.cpp $(CLUSTER_SRC_DIR)/lsh/lsh.cpp $(CLUSTER_SRC_DIR)/hypercube/lsh_hypercube.cpp
CLUSTER_INCDIR := includes
CLUSTER_CFLAGS := -Wall -Wextra -I$(INCDIR)
CLUSTER_TARGET := cluster

$(CLUSTER_TARGET): $(CLUSTER_OBJ_FILES)
	g++ $(CLUSTER_CFLAGS) -o $@ $^
$(CLUSTER_OBJ_DIR)/%.o: $(CLUSTER_SRC_DIR)/%.cpp
	g++ $(CLUSTER_CFLAGS) -c -o $@ $<

run_clustering_classic:
	./$(CLUSTER_TARGET) -i data/d10.bin -c data/cluster.conf -o outputs/output_cluster.txt -complete Yes -m Classic
run_clustering_lsh:
	./$(CLUSTER_TARGET) -i data/train-images.idx3-ubyte -c data/cluster.conf -o outputs/output_cluster.txt -complete Yes -m LSH
run_clustering_hc:
	./$(CLUSTER_TARGET) -i data/train-images.idx3-ubyte -c data/cluster.conf -o outputs/output_cluster.txt -complete Yes -m Hypercube


.PHONY: all clean
# Clean the build directory
clean:
	rm -rf $(OBJ_DIR)/*.o $(TARGET) $(CLUSTER_TARGET) $(TARGET_CUBE)



