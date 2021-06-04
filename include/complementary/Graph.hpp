#ifndef GRAPH_COMP_HPP
#define GRAPH_COMP_HPP

#include <iostream>
#include <stack>
#include <algorithm>
#include "Tree.hpp"
#include "UnionFind.hpp"
#include <sdsl/int_vector.hpp>

using namespace std;
using namespace sdsl;

#define CLOSE_PAR 0
#define OPEN_PAR 1
#define EMPTY 2
#define NONEN 3


class Graph {

private:
  Vertex *V;      // Array of vertices of the tree
  Edge *E;        // Array of edges of the tree. It is the concatenation of the
	          // adjacency lists of all nodes
  unsigned int n; // Number of vertices in the tree
  unsigned int m; // Number of edges in the tree

public:
  Graph () {
    this->n = 0;
    this->m = 0;
  }

  Graph (unsigned int n, unsigned int m) {
    this->n = n;
    this->m = m;
    this->V = new Vertex[n];
    this->E = new Edge[2*m];
  }

  unsigned int vertices() {
    return n;
  }

  unsigned int edges() {
    return m;
  }

  // E[i].src = src and E[i].tgt = tgt
  void setEdge(int i, unsigned int src, unsigned int tgt) {
    this->E[i].setSrc(src);
    this->E[i].setTgt(tgt);
    this->E[i].setCmp(-1);
    this->E[i].setOrientation(1); // Direct
  }

  // V[i].first = first and V[i].last = last
  void setVertex(int i, unsigned int first, unsigned int last) {
    this->V[i].setFirst(first);
    this->V[i].setLast(last);
  }

  // V[i].first = first
  void setVertexFirst(int i, unsigned int first) {
    this->V[i].setFirst(first);
  }

  // V[i].last = last
  void setVertexLast(int i, unsigned int last) {
    this->V[i].setLast(last);
  }

  void setEdgeSrc(int i, unsigned int s) {
    this->E[i].setSrc(s);
  }

  void setEdgeTgt(int i, unsigned int t) {
    this->E[i].setTgt(t);
  }

  void setEdgeCmp(int i, int c) {
    this->E[i].setCmp(c);
  }

  void setEdgeOrientation(int i, bool o) {
    this->E[i].setOrientation(o);
  }

  unsigned int getEdgeSrc(int i) {
    return this->E[i].getSrc();
  }

  unsigned int getEdgeTgt(int i) {
    return this->E[i].getTgt();
  }

  int getEdgeCmp(int i) {
    return this->E[i].getCmp();
  }

  bool getEdgeOrientation(int i) {
    return this->E[i].getOrientation();
  }

  // Return the next edge in counter-clockwise order
  unsigned int getNextEdgeCCW(int i) {

    uint src = this->getEdgeSrc(i);
    uint first = this->getVertexFirst(src);
    uint last = this->getVertexLast(src);

    return ((uint)(i+1) > last)? first : i+1;
  }

  // Return the next edge in clockwise order
  unsigned int getNextEdgeCW(int i) {

    uint src = this->getEdgeSrc(i);
    int first = (int)this->getVertexFirst(src);
    int last = (int)this->getVertexLast(src);

    return (i-1 < first)? last : i-1;
  }

  Vertex getVertex(int i) {
    return this->V[i];
  }

  Edge getEdge(int i) {
    return this->E[i];
  }

  unsigned int getVertexFirst(int i) {
    return this->V[i].getFirst();
  }

  unsigned int getVertexLast(int i) {
    return this->V[i].getLast();
  }

  Tree dfs_spanning_tree(unsigned int init, int *parent, unsigned int
				*count_edges, unsigned int *references) {
    unsigned int n = this->vertices();
    unsigned int m = this->edges();
    Tree t(n);


    char *visited = new char[n](); // TODO: Change to a boolean array
    unsigned int *edges = new unsigned int[2*m]();
    unsigned int num_tree_edges = 2*(n-1);

    stack <unsigned int> s;
    visited[init] = 1;
    s.push(init);
    parent[init] = -1;

    while(!s.empty()) {
      unsigned int curr = s.top(); s.pop();

      for(unsigned int i = this->V[curr].getFirst(); i <= this->V[curr].getLast(); i++)
	{
	  unsigned int tgt = this->E[i].getTgt();

	  if(visited[tgt] == 0) { // Not visited
	    visited[tgt] = 1;
	    s.push(tgt);
	    parent[tgt] = this->E[i].getCmp(); // Edge child-to-parent
	  }
	}
    }

    // Marking the edges of G that are in T
    for(unsigned int i = 0; i < init; i++) {
      edges[(unsigned int)parent[i]] = 1;
      edges[(unsigned int)this->E[parent[i]].getCmp()] = 1;
    }
    for(unsigned int i = init+1; i < n; i++) {
      edges[(unsigned int)parent[i]] = 1;
      edges[(unsigned int)this->E[parent[i]].getCmp()] = 1;
    }


    unsigned int mm = 0;
    /* Counting edges */
    for(unsigned int i = 0; i < n; i++) {
      unsigned int first = this->V[i].getFirst();
      unsigned int last = this->V[i].getLast();
      unsigned int zeros = 0;
      int last_one = -1;
      unsigned int carry_last = 0;

      for(unsigned int j = first; j <= last; j++) {
	if(edges[j] == 1) {
	  t.setEdge(mm, this->E[j].getSrc(), this->E[j].getTgt());

	  if(last_one == -1)
	    carry_last = zeros;
	  else
	    count_edges[last_one] += zeros;

	  references[mm] = j;
	  last_one = mm;
	  zeros = 0;

	  if(mm == 0)
	    t.setNodeFirst(t.getEdgeSrc(mm), mm);
	  else if(t.getEdgeSrc(mm) != t.getEdgeSrc(mm-1)) {
	    t.setNodeLast(t.getEdgeSrc(mm-1), mm-1);
	    t.setNodeFirst(t.getEdgeSrc(mm), mm);
	  }
	  mm++;
	}
	else
	  zeros++;
	edges[j] = mm;
      }

      count_edges[last_one] += zeros + carry_last;
      t.setNodeLast(t.getEdgeSrc(mm-1), mm-1);
    }

    for(unsigned int i = 0; i < init; i++)
      parent[i] = edges[parent[i]]-1;
    for(unsigned int i = init+1; i < n; i++)
      parent[i] = edges[parent[i]]-1;

    for(unsigned int i = 0; i < num_tree_edges; i++) {
      Vertex tgt = t.getNode(t.getEdgeTgt(i));

      for(unsigned int j = tgt.getFirst(); j <= tgt.getLast(); j++)
    	if(t.getEdgeTgt(j) == t.getEdgeSrc(i)) {
    	  t.setEdgeCmp(i, j);
    	  break;
    	}
    }

    return t;
  }

  bool connected_graph() {
    unsigned int n = this->vertices();
    stack <unsigned int> s;
    char *visited = new char[n]();
    int curr = -1;
    unsigned int edge;
    int first = 1;
    int num_vertices = 0;

    while(!s.empty() || first) {
      if(first) { // Root
	curr = 0;
	first = 0;
      }
      else {
	edge = s.top(); s.pop();
	curr = this->E[edge].getTgt(); // current
      }
      visited[curr] = 1;
      for(unsigned int i = this->V[curr].getFirst(); i <= this->V[curr].getLast();
	  i++) {
	if(!visited[this->E[i].getTgt()])
	  s.push(i);
      }
    }

    for(unsigned int i = 0; i < n; i++)
      if(visited[i] == 0) { // There are unvisited vertices
	num_vertices++;
      }

    cout << "unvisited vertices: " << num_vertices << ", visited vertices: " <<
      n-num_vertices << endl;
      return num_vertices==0;
  }

  int_vector<> ps_tree_encoding() {
    unsigned int n = this->vertices();
    unsigned int m = this->edges();

    int_vector<> S(4*n-5, NONEN, 8); // Output string
    char *visited = new char[2*m]();
    /*** Initial setting ***/

    // Set the external edges as visited (they are not considered in the
    // traversal)
    uint ext0 = 0; // By definition
    uint ext1 = this->getEdgeTgt(this->getVertexFirst(ext0));
    uint ext2 = this->getEdgeTgt(this->getVertexLast(ext0));
    visited[this->getVertexFirst(ext0)] = 1;
    visited[this->getVertexLast(ext0)] = 1;
    visited[this->getVertexFirst(ext1)] = 1;
    visited[this->getVertexLast(ext1)] = 1;
    visited[this->getVertexFirst(ext2)] = 1;
    visited[this->getVertexLast(ext2)] = 1;

    // For any maximal planar graph, the following entries of S are already defined
    S[0] = OPEN_PAR;       S[4*n-6] = CLOSE_PAR;
    S[1] = OPEN_PAR;       S[4*n-7] = CLOSE_PAR;
    S[2] = OPEN_PAR;       S[4*n-8] = CLOSE_PAR;
    S[3] = EMPTY;

    /*** Traversal ***/

    uint v = 0; // Starting vertex (by definition)
    uint e = this->getVertexLast(v); // Starting edge (by definition)
    uint end_e = this->getVertexFirst(v); // Stop condition
    uint ee = this->getNextEdgeCW(e);
    uint id = 4;

    while(ee != end_e) {
      if(!visited[ee] && this->getEdgeOrientation(ee)) {
    	visited[ee] = 1;
    	visited[this->getEdgeCmp(ee)] = 1;
    	S[id++] = EMPTY;
    	e = ee;
      }
      else if(!visited[ee] && !this->getEdgeOrientation(ee)) {
    	visited[ee] = 1;
    	visited[this->getEdgeCmp(ee)] = 1;
    	S[id++] = OPEN_PAR;
    	e = this->getEdgeCmp(ee);
      }
      else if(visited[ee] && this->getEdgeOrientation(ee)) {
    	S[id++] = CLOSE_PAR;
    	e = this->getEdgeCmp(ee);
      }
      else {
    	e = ee;
      }
      ee = this->getNextEdgeCW(e);

    }

    /* Convert the cw traversal into a ccw traversal */

    int mid = (4*n-5) / 2;

    // Reverse the sequence S
    for(int i=0, j=4*n-6; i < mid; i++, j--) {
      char cp = S[i];
      S[i] = S[j];
      S[j] = cp;
    }

    // Change close parentheses by open parentheses
    for(uint i=0; i < 4*n-5; i++) {
      if(S[i] == OPEN_PAR) S[i] = CLOSE_PAR;
      else if(S[i] == CLOSE_PAR) S[i] = OPEN_PAR;
    }

    return S;
  }

  void debugPrintEdges(unsigned int vertex){

    for(unsigned int i = getVertexFirst(vertex); i<=getVertexLast(vertex); i++){
      std::cerr << i <<
      " src=" << getEdgeSrc(i) <<
      " tgt=" << getEdgeTgt(i) <<
      " cmp=" << getEdgeCmp(i) <<
      " cmpSrc=" << getEdgeSrc(getEdgeCmp(i)) <<
      " cmpTgt=" << getEdgeTgt(getEdgeCmp(i)) <<
      " cmpCmp=" << getEdgeCmp(getEdgeCmp(i)) << endl;
    }
  }

  std::vector<unsigned int> externalFace(){
    std::vector<unsigned int> res;
    unsigned int start = 0;
    int iters = 0;
    unsigned int currNode = start;
    unsigned int currEdge = getVertexFirst(start);
    unsigned int startEdge = currEdge;
    unsigned int nextNode = getEdgeTgt(currEdge);
    // std::cerr << "currNode " << start << std::endl;
    // std::cerr << "currEdge " << currEdge << std::endl;
    // std::cerr << "nextNode " << nextNode << std::endl;
    do{
      // search for complementary edge, find next in counter counterclockwise
      // and continue
      res.push_back(currEdge);
      // std::cerr << "insert edge " << getEdgeSrc(currEdge) << "->" << getEdgeTgt(currEdge) << std::endl;
      //std::cerr << "currEdge " << currEdge << std::endl;
      currEdge = getEdgeCmp(currEdge);
      currEdge = currEdge == getVertexLast(nextNode) ? getVertexFirst(nextNode) : currEdge+1;
      currNode = getEdgeSrc(currEdge);
      nextNode = getEdgeTgt(currEdge);
      //std::cerr << "currNode " << currNode <<std::endl;
      //std::cerr << "nextNode " << nextNode <<std::endl;
    }while(currEdge!=startEdge);
    return res;
  }

  bool generateDfsSpanningTree(std::vector<bool> &belongsToT, unsigned int initEdge){
    unsigned int init = getEdgeSrc(initEdge);
    belongsToT.assign(2*edges(),false);
    std::vector<bool> visited(vertices(),false);
    std::stack<unsigned int> s;
    s.push(init);
    visited[init] = true;

    while(!s.empty()){
      unsigned int node = s.top();
      //std::cerr << node << std::endl;
      s.pop();
      for(unsigned int nextEdge = getVertexFirst(node); nextEdge <= getVertexLast(node); ++nextEdge ){
        unsigned int nextNode = getEdgeTgt(nextEdge);
        if(!visited[nextNode]){
          visited[nextNode] = true;
          s.push(nextNode);
          belongsToT[nextEdge] = belongsToT[getEdgeCmp(nextEdge)] = true;
        }
      }
    }
    // std::cerr << "aristas que pertenecen: " << std::endl;
    // for(unsigned int i = 0; i<belongsToT.size();i++){
    //   if(belongsToT[i]){
    //     std::cerr << getEdgeSrc(i) << "->" << getEdgeTgt(i) << std::endl;
    //   }
    // }
    return true;
  }

  bool generateLimitedDfsSpanningTree(std::vector<bool> &belongsToT, unsigned int initEdge, unsigned int maxHeight ){
    unsigned int init = getEdgeSrc(initEdge);
    belongsToT.assign(2*edges(),false);
    std::vector<bool> visited(vertices(),false);
    // std::vector<unsigned int> parentEdge(vertices(), -1);
    std::vector<unsigned int> height(vertices(), 0);
    std::stack<pair<unsigned int,unsigned int> > s;
    s.push(pair<unsigned int,unsigned int>(init,-1) );
    //visited[init] = true;

    //doing this to avoid checking if node is initial in every iteration
    // for(unsigned int nextEdge = getVertexFirst(init); nextEdge <= getVertexLast(init); ++nextEdge ){
    //   unsigned int nextNode = getEdgeTgt(nextEdge);
    //   parentEdge[nextNode] = nextEdge;
    //   s.push(nextNode);
    // }
    std::queue<pair<unsigned int,unsigned int> > frontierNodes;
    while(!s.empty()){
      unsigned int node = s.top().first;
      unsigned int parentEdge = s.top().second;
      unsigned int parentHeight = parentEdge==-1 ? 0 : height[getEdgeSrc(parentEdge)];
      //std::cerr << node << std::endl;
      s.pop();
      if(!visited[node] && parentHeight<maxHeight){
        // std::cerr << "nodo "<<node<<std::endl;
        visited[node] = true;
        height[node] = parentHeight+1;
        //s.push(nextNode);
        if(parentEdge!=-1) belongsToT[parentEdge] = belongsToT[getEdgeCmp(parentEdge)] = true;
        for(unsigned int nextEdge = getVertexFirst(node); nextEdge <= getVertexLast(node); ++nextEdge ){
          unsigned int nextNode = getEdgeTgt(nextEdge);
          // parentEdge[nextNode] = nextEdge;
          s.push(pair<unsigned int, unsigned int>(nextNode,nextEdge));
        }
      }else if(!visited[node]){
        frontierNodes.push(pair<unsigned int, unsigned int>(node, parentEdge));
      }
    }

    if(frontierNodes.empty()){
      std::cerr << "todos durante dfs" << std::endl;
    }
    while(!frontierNodes.empty()){
      unsigned int node = frontierNodes.front().first;
      //unsigned int nodeParentEdge = parentEdge[node];
      unsigned int parentEdge = frontierNodes.front().second;
      unsigned int parentNode = getEdgeSrc(parentEdge);
      //std::cerr << node << std::endl;
      frontierNodes.pop();
      if(!visited[node]){
        visited[node] = true;
        //height[node] = height[parentNode]+1;
        //s.push(nextNode);
        belongsToT[parentEdge] = belongsToT[getEdgeCmp(parentEdge)] = true;
        for(unsigned int nextEdge = getVertexFirst(node); nextEdge <= getVertexLast(node); ++nextEdge ){
          unsigned int nextNode = getEdgeTgt(nextEdge);
          // parentEdge[nextNode] = nextEdge;
          frontierNodes.push(pair<unsigned int, unsigned int>(nextNode, nextEdge));
        }
      }
    }
    // bool working = true;
    // for(int i=0; i<vertices(); i++){
    //   if(!visited[i]){
    //     working = false;
    //   }
    // }
    // if(!working){
    //   std::cerr << "faltan nodos" << std::endl;
    //   return false;
    // }
    //
    // std::cerr << "aristas que pertenecen: " << std::endl;
    // for(unsigned int i = 0; i<belongsToT.size();i++){
    //   if(belongsToT[i]){
    //     std::cerr << getEdgeSrc(i) << "->" << getEdgeTgt(i) << std::endl;
    //   }
    // }
    return true;
  }


  bool generateBfsSpanningTree(std::vector<bool> &belongsToT, unsigned int initEdge){
    unsigned int init = getEdgeSrc(initEdge);
    belongsToT.assign(2*edges(),false);
    std::vector<bool> visited(vertices(),false);
    std::queue<unsigned int> s;
    s.push(init);
    visited[init] = true;

    while(!s.empty()){
      unsigned int node = s.front();
      //std::cerr << node << std::endl;
      s.pop();
      for(unsigned int nextEdge = getVertexFirst(node); nextEdge <= getVertexLast(node); ++nextEdge ){
        unsigned int nextNode = getEdgeTgt(nextEdge);
        if(!visited[nextNode]){
          visited[nextNode] = true;
          s.push(nextNode);
          belongsToT[nextEdge] = belongsToT[getEdgeCmp(nextEdge)] = true;
        }
      }
    }
    // std::cerr << "aristas que pertenecen: " << std::endl;
    // for(unsigned int i = 0; i<belongsToT.size();i++){
    //   if(belongsToT[i]){
    //     std::cerr << getEdgeSrc(i) << "->" << getEdgeTgt(i) << std::endl;
    //   }
    // }
    return true;
  }

  // reorder edges so that selected edge is first of its source vertex
  // returns position of edge after rotation
  unsigned int rotateVertexEdges(unsigned int edge){
    // std::cerr << "rotating around edge n=" << edge << " [" <<
    //   getEdgeSrc(edge) << "," << getEdgeTgt(edge) <<endl;
    unsigned int vertex = getEdgeSrc(edge);

    unsigned int first = getVertexFirst(vertex), last = getVertexLast(vertex);

    // std::cerr << "edge first=" << first << " last=" << last << std::endl;

    unsigned int edges = last - first + 1;
    unsigned int delta = edge - first;

    // if(!delta){
    //   return first;
    // }

    unsigned int cmps[edges];
    for(unsigned int i = first; i <= last; i++){
      cmps[i-first] = getEdgeCmp(i);
    }

    for(unsigned int i=0; i< edges; i++){
      unsigned int cmp = cmps[(i+delta)%edges];
      setEdge(first + i, getEdgeTgt(cmp), getEdgeSrc(cmp));
      setEdgeCmp(first + i, cmp);
      setEdgeCmp(cmp, first + i);
    }
    return first;
  }

  //
  unsigned int markEdgesFaceEdge(unsigned int edge, std::vector<unsigned int> &edgeFace,  std::vector<unsigned int> &faceEdge, unsigned int currFaceId){
    if(edgeFace[edge]==-1){
      faceEdge.push_back(edge);
      unsigned int currFaceEdge = edge;
      do{
        edgeFace[currFaceEdge] = currFaceId;
        unsigned int cmp = getEdgeCmp(currFaceEdge);
        currFaceEdge = getNextEdgeCCW(cmp);
      }while(currFaceEdge!=edge);
      return currFaceId+1;

    }else{

      return currFaceId;
    }
  }

  //generate tree based on a BFS on dual graph
  bool generateBfsDualSpanningTree(std::vector<bool> &belongsToT, unsigned int initEdge){
      std::cerr << "construccion con dual" << std::endl;
      belongsToT.assign(2*edges(),true);
      //face to which edge i belongs
      std::vector<unsigned int> edgeFace(2*edges(), -1);
      //an edge which belongs to face i
      std::vector<unsigned int> faceEdge;
      unsigned int currFaceId = 0;
      for(unsigned int i=0; i<2*edges(); i++){
        currFaceId = markEdgesFaceEdge(i, edgeFace, faceEdge, currFaceId);
      }
      std::cerr << "graph has " << currFaceId << " faces" << std::endl;
      std::queue<unsigned int> q;
      std::vector<bool> visited(currFaceId, false);
      unsigned initFace = edgeFace[initEdge];

      q.push(initFace);
      visited[initFace] = true;

      while(!q.empty()){
        unsigned int currFace = q.front();
        q.pop();

        unsigned int initFaceEdge = faceEdge[currFace];
        unsigned int currFaceEdge = initFaceEdge;

        do {
          unsigned int cmpEdge = getEdgeCmp(currFaceEdge);
          unsigned int neighborFace = edgeFace[cmpEdge];
          if(!visited[neighborFace]){
            belongsToT[cmpEdge] = belongsToT[getEdgeCmp(cmpEdge)] = false;
            q.push(neighborFace);
            visited[neighborFace] = true;
          }
          currFaceEdge = getNextEdgeCCW(cmpEdge);
        } while(currFaceEdge!=initFaceEdge);

      }

      // std::cerr << "aristas que pertenecen: " << std::endl;
      // for(unsigned int i = 0; i<belongsToT.size();i++){
      //   if(belongsToT[i]){
      //     std::cerr << getEdgeSrc(i) << "->" << getEdgeTgt(i) << std::endl;
      //   }
      // }
      return true;
  }

  //generate tree based on a BFS on dual graph
  bool generateHeuristicBfsSpanningTree(std::vector<bool> &belongsToT, unsigned int initEdge){
    std::cerr << "construccion con heuristico" << std::endl;
    unsigned initNode = getEdgeSrc(initEdge);
      belongsToT.assign(2*edges(),true);
      std::vector<bool> usedEdge(2*edges(),false);
      //face to which edge i belongs
      std::vector<unsigned int> edgeFace(2*edges(), -1);
      //an edge which belongs to face i
      std::vector<unsigned int> faceEdge;
      unsigned int currFaceId = 0;
      for(unsigned int i=0; i<2*edges(); i++){
        currFaceId = markEdgesFaceEdge(i, edgeFace, faceEdge, currFaceId);
      }
      std::cerr << "graph tiene " << currFaceId << " faces" << std::endl;
      std::queue<unsigned int> qFaces;
      std::queue<unsigned int> qNodes;
      std::vector<unsigned int> heightNode(vertices(),0);
      std::vector<unsigned int> heightFace(currFaceId,0);
      std::vector<bool> visitedFace(currFaceId, false);
      std::vector<bool> visitedNodes(vertices(), false);
      unsigned initFace = edgeFace[initEdge];

      qFaces.push(initFace);
      qNodes.push(initNode);

      visitedFace[initFace] = true;
      visitedNodes[initNode] = true;
      /*
        basic idea, we try doing a BFS on T and T* at the same time
        to do so, we alternate between
      */
      while(!qFaces.empty()){

        //we first try expanding T as to "block" some edges to be avoided in
        // T* construction

        int currHeight = -1;
        if(!qNodes.empty()) currHeight = heightNode[qNodes.front()];
        while(!qNodes.empty() && currHeight==heightNode[qNodes.front()]){
          unsigned int currNode = qNodes.front();
          qNodes.pop();
          for(unsigned int nextEdge = getVertexFirst(currNode); nextEdge <= getVertexLast(currNode); ++nextEdge ){
            unsigned int nextNode = getEdgeTgt(nextEdge);
            if(!visitedNodes[nextNode] && !usedEdge[nextEdge]){
               // std::cerr << "bloqueo arista " << getEdgeSrc(nextEdge) << "->" << getEdgeTgt(nextEdge) << std::endl;
              usedEdge[nextEdge] = usedEdge[getEdgeCmp(nextEdge)] = true;
              visitedNodes[nextNode] = true;
              heightNode[nextNode] = currHeight+1;
              qNodes.push(nextNode);
              //belongsToT[nextEdge] = belongsToT[getEdgeCmp(nextEdge)] = true;
            }
          }
          // continue;
        }

        currHeight = qFaces.empty() ? -1 : heightFace[qFaces.front()];
        while(!qFaces.empty() && currHeight == heightFace[qFaces.front()]){

          unsigned int currFace = qFaces.front();
          qFaces.pop();
          // std::cerr << "continuo con cara " << currFace << std::endl;
          unsigned int initFaceEdge = faceEdge[currFace];
          unsigned int currFaceEdge = initFaceEdge;

          do {


            unsigned int cmpEdge = getEdgeCmp(currFaceEdge);
            if(usedEdge[currFaceEdge]){
              currFaceEdge = getNextEdgeCCW(cmpEdge);
              // std::cerr << "la arista " << getEdgeSrc(currFaceEdge) << " " << getEdgeTgt(currFaceEdge) <<  " esta bloqueada" << std::endl;
              continue;
            }

            unsigned int neighborFace = edgeFace[cmpEdge];
            if(!visitedFace[neighborFace]){
              // std::cerr << "la arista " << getEdgeSrc(currFaceEdge) << " " << getEdgeTgt(currFaceEdge) <<  " se agrega" << std::endl;
              belongsToT[cmpEdge] = belongsToT[getEdgeCmp(cmpEdge)] = false;
              usedEdge[cmpEdge] = usedEdge[getEdgeCmp(cmpEdge)] = true;
              qFaces.push(neighborFace);
              visitedFace[neighborFace] = true;
              heightFace[neighborFace] = currHeight +1;
            }
            currFaceEdge = getNextEdgeCCW(cmpEdge);
          } while(currFaceEdge!=initFaceEdge);

        }

      }

      // std::cerr << "aristas que pertenecen: " << std::endl;
      // for(unsigned int i = 0; i<belongsToT.size();i++){
      //   if(belongsToT[i]){
      //     std::cerr << getEdgeSrc(i) << "->" << getEdgeTgt(i) << std::endl;
      //   }
      // }
      return true;
  }

  bool dumpReducedGraph(std::ostream &out, int percentage){
    double totalPer = (double)percentage*0.01f;
    unsigned int edgesToKeepTotal = (unsigned int)(edges()*totalPer);

    unsigned int edgesToRemove = edges()-edgesToKeepTotal;

    if(edgesToKeepTotal<vertices()-1){
      //cant generate a connected graph with less than n-1 edges
      return false;
    }

    //unsigned int edgesToRemove = edges()-edgesToKeep;

    std::vector<bool> currMarkedForErase(2*edges(),false);
    std::vector<bool> cantTouchThis(2*edges(),false);
    //we use an UnionFind structure to create a random spanning tree
    //as to keep track of connected components
    UnionFind uf(vertices());

    //we won't try to delete the first edge
    cantTouchThis[0] = cantTouchThis[getEdgeCmp(0)] = true;
    //we create edge in UF
    uf.unionSet(getEdgeSrc(0),getEdgeTgt(0));
    //first we generate a random tree
    std::cout << "creating spanning tree "<< std::endl;
    for(unsigned int i=1; i<vertices()-1;i++){
      //we pick a random edge
      // if(i%1000==1)
      //  std::cout << "a "<<i<< std::endl;
      unsigned int randomEdge = rand()%(2*edges());
      if(cantTouchThis[randomEdge]){
        //already used, try again;
        i--;
      }else{
        unsigned int src = getEdgeSrc(randomEdge);
        unsigned int tgt = getEdgeTgt(randomEdge);
        if(uf.isSameSet(src,tgt)){
          //we can't use an edge which would create loop
          i--;
        }else{
          //edge belongs to tree, won't use for deletion
          uf.unionSet(src,tgt);
          cantTouchThis[randomEdge] = cantTouchThis[getEdgeCmp(randomEdge)] = true;
        }
      }

    }
    std::cout << "picking edgesToRemove" << std::endl;
    std::vector<unsigned int> edgeToPick(2*edges(),0);
    for(unsigned int i=0; i<2*edges(); i++){
      edgeToPick[i] = i;
    }
    unsigned int idx = 0;
    std::random_shuffle(edgeToPick.begin(), edgeToPick.end());
    for(unsigned int i=0; i<edgesToRemove; i++){
      unsigned int randomEdge = edgeToPick[idx++];
      if(cantTouchThis[randomEdge]){
        i--;
      }else{
        cantTouchThis[randomEdge] = cantTouchThis[getEdgeCmp(randomEdge)] = true;
        currMarkedForErase[randomEdge] = currMarkedForErase[getEdgeCmp(randomEdge)] = true;
      }
    }

    if(true){
      //unsigned int edgesToKeep = (unsigned int)(edges()-currentlyDeleted);
      std::cout <<edgesToKeepTotal  <<"/"<<edges() <<std::endl;
      //writing to string first is much faster;
      string output;
      output += to_string(vertices()) + "\n";
      // out << vertices() << std::endl;
      // out << edgesToKeep << std::endl;
      output += to_string(edgesToKeepTotal) + "\n";
      for(unsigned int src = 0; src< vertices(); src++){
        for(unsigned int nextEdge = getVertexFirst(src); nextEdge <= getVertexLast(src); ++nextEdge ){
          if(currMarkedForErase[nextEdge]) continue;
          unsigned int tgt = getEdgeTgt(nextEdge);
          // out << src << " " << tgt << std::endl;
          output += to_string(src) + " " + to_string(tgt) + "\n";

        }
      }
      out << output;
      return true;
    }
  }

};

#endif
