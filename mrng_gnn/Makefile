SRC_DIR := modules/gnn
OBJ_DIR := obj
SRC_FILES := $(wildcard $(SRC_DIR)/*.cpp)
OBJ_FILES := $(patsubst $(SRC_DIR)/%.cpp,$(OBJ_DIR)/%.o,$(SRC_FILES))
INCDIR := include
CFLAGS := -Wall -Wextra -I$(INCDIR) -I$(INCDIR)/lsh_includes -Llib
LIBS := lib/new_lib.a
TARGET := graph_search

$(TARGET): $(OBJ_FILES)
	g++ $(CFLAGS) -o $@ $^ -lpthread $(LIBS)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
	g++ $(CFLAGS) -c -o $@ $<

clean:
	rm -rf  $(OBJ_DIR)/*.o graph_search

run_new: 
	./graph_search -d data/d30.bin -q data/q30.bin -o outputs/gnn_output.txt -m 1 -g data/gnn_graph.txt -N 1 -R 100 -E 50


run_gnn_construct_new:
	./graph_search -d data/d30.bin -q data/q30.bin -o outputs/gnn_output.txt -m 1 -N 1 -R 1 -E 1

run_mrng_construct_new:
	./graph_search -d data/d30.bin -q data/q30.bin -o outputs/mrng_output.txt -m 2 -N 1 -l 20

run_mrng_new:
	./graph_search -d data/d30.bin -q data/q30.bin -o outputs/mrng_output.txt -m 2 -g data/graph.txt -N 1 -l 1000