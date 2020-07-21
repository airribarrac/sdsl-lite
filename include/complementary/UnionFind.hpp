#ifndef UNIONFIND_HPP
#define UNIONFIND_HPP


// code obtained from Competitive Programming 3 Book by Steven Halim and Felix Halim

class UnionFind {                                              // OOP style
private:
  std::vector<unsigned int> p, rank;
  int numSets;
  int initialN;
public:
  UnionFind(int N) {
    initialN = N;
     numSets = N; rank.assign(N, 0);
    p.assign(N, 0); for (int i = 0; i < N; i++) p[i] = i; }
  int findSet(int i) { return (p[i] == i) ? i : (p[i] = findSet(p[i])); }
  bool isSameSet(int i, int j) { return findSet(i) == findSet(j); }
  void unionSet(int i, int j) {
    if (!isSameSet(i, j)) { numSets--;
    int x = findSet(i), y = findSet(j);
    // rank is used to keep the tree short
    if (rank[x] > rank[y]) { p[y] = x;  }
    else                   { p[x] = y;
                             if (rank[x] == rank[y]) rank[y]++; } } }
  int numDisjointSets() { return numSets; }
  //int sizeOfSet(int i) { return setSize[findSet(i)]; }
  void reset(){
    numSets = initialN;
    for (int i = 0; i < numSets; i++){
      p[i] = i;
      rank[i] = 0;
    }
  }
};

#endif
