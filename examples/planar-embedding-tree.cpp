#include <sdsl/pemb.hpp>
#include <sdsl/pembSuccinct.hpp>
#include <complementary/GraphVisitor.hpp>
#include <complementary/Graph.hpp>
#include <complementary/utils.hpp>

using namespace sdsl;
using namespace std;

int main(int argc, char **argv) {
  // argv[1] is the path to a file with the planar embedding.
  // To check the input format, visit https://users.dcc.uchile.cl/~jfuentess/datasets/graphs.php
  Graph g = read_graph_from_file(argv[1]);
  pembSuccinct<> pe(g,0,new DegHeurSTGenerator());

  cout << "Size in bytes: " << size_in_bytes(pe) << " B" << endl;
  cout << "Degree of vertex 10: " << pe.degree(10) << endl;
}
