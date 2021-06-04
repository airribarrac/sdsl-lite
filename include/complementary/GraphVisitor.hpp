#ifndef GRAPH_VISITOR_COMP_HPP
#define GRAPH_VISITOR_COMP_HPP

#include "Graph.hpp"
#include <vector>
#include <stack>
#include <queue>


//interface for generating spanning trees
class GraphSpanningTreeGenerator{
public:
  void setGraph(Graph *g){ _g=g;}
  void setInitialEdge(unsigned int initialEdge){ _initialEdge=initialEdge;}
  virtual vector<bool> getSpanningTree()=0;
protected:
  Graph *_g=NULL;
  unsigned int _initialEdge=0;
};

class DfsPrimalSTGenerator : public GraphSpanningTreeGenerator{
public:


  vector<bool> getSpanningTree() override {
    unsigned int init = _g->getEdgeSrc(_initialEdge);
    vector<bool> belongsToT(2*_g->edges(),false);

    std::vector<bool> visited(_g->vertices(),false);
    std::stack<unsigned int> s;
    s.push(init);
    visited[init] = true;

    while(!s.empty()){
      unsigned int node = s.top();
      s.pop();
      for(unsigned int nextEdge = _g->getVertexFirst(node); nextEdge <= _g->getVertexLast(node); ++nextEdge ){
        unsigned int nextNode = _g->getEdgeTgt(nextEdge);
        if(!visited[nextNode]){
          visited[nextNode] = true;
          s.push(nextNode);
          belongsToT[nextEdge] = belongsToT[_g->getEdgeCmp(nextEdge)] = true;
        }
      }
    }

    return belongsToT;
  }
};

class BfsPrimalSTGenerator: public GraphSpanningTreeGenerator{
public:

  vector<bool> getSpanningTree() override {
    unsigned int init = _g->getEdgeSrc(_initialEdge);
    vector<bool> belongsToT(2*_g->edges(),false);
    std::vector<bool> visited(_g->vertices(),false);
    std::queue<unsigned int> s;
    s.push(init);
    visited[init] = true;

    while(!s.empty()){
      unsigned int node = s.front();
      //std::cerr << node << std::endl;
      s.pop();
      for(unsigned int nextEdge = _g->getVertexFirst(node); nextEdge <= _g->getVertexLast(node); ++nextEdge ){
        unsigned int nextNode = _g->getEdgeTgt(nextEdge);
        if(!visited[nextNode]){
          visited[nextNode] = true;
          s.push(nextNode);
          belongsToT[nextEdge] = belongsToT[_g->getEdgeCmp(nextEdge)] = true;
        }
      }
    }
    return belongsToT;
  }
};



class BfsDualSTGenerator: public GraphSpanningTreeGenerator{
public:

  vector<bool> getSpanningTree() override {
    unsigned int initEdge = _initialEdge;
    vector<bool> belongsToT(2*_g->edges(),true);
    //face to which edge i belongs
    std::vector<unsigned int> edgeFace(2*_g->edges(), -1);
    //an edge which belongs to face i
    std::vector<unsigned int> faceEdge;
    unsigned int currFaceId = 0;

    for(unsigned int edge=0; edge<2*_g->edges(); edge++){

      if(edgeFace[edge]==-1){
        faceEdge.push_back(edge);
        unsigned int currFaceEdge = edge;
        do{
          edgeFace[currFaceEdge] = currFaceId;
          unsigned int cmp = _g->getEdgeCmp(currFaceEdge);
          currFaceEdge = _g->getNextEdgeCCW(cmp);
        }while(currFaceEdge!=edge);
        currFaceId++;
      }

    }
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
        unsigned int cmpEdge = _g->getEdgeCmp(currFaceEdge);
        unsigned int neighborFace = edgeFace[cmpEdge];
        if(!visited[neighborFace]){
          belongsToT[cmpEdge] = belongsToT[_g->getEdgeCmp(cmpEdge)] = false;
          q.push(neighborFace);
          visited[neighborFace] = true;
        }
        currFaceEdge = _g->getNextEdgeCCW(cmpEdge);
      } while(currFaceEdge!=initFaceEdge);

    }

    return belongsToT;
  }
};

class AltSTGenerator: public GraphSpanningTreeGenerator{
public:

  vector<bool> getSpanningTree() override {
      unsigned int initEdge = _initialEdge;
      unsigned initNode = _g->getEdgeSrc(initEdge);
      vector<bool> belongsToT(2*_g->edges(),true);
      std::vector<bool> usedEdge(2*_g->edges(),false);
      //face to which edge i belongs
      std::vector<unsigned int> edgeFace(2*_g->edges(), -1);
      //an edge which belongs to face i
      std::vector<unsigned int> faceEdge;
      unsigned int currFaceId = 0;


      for(unsigned int edge=0; edge<2*_g->edges(); edge++){

        if(edgeFace[edge]==-1){
          faceEdge.push_back(edge);
          unsigned int currFaceEdge = edge;
          do{
            edgeFace[currFaceEdge] = currFaceId;
            unsigned int cmp = _g->getEdgeCmp(currFaceEdge);
            currFaceEdge = _g->getNextEdgeCCW(cmp);
          }while(currFaceEdge!=edge);
          currFaceId++;
        }

      }


      std::queue<unsigned int> qFaces;
      std::queue<unsigned int> qNodes;
      std::vector<unsigned int> heightNode(_g->vertices(),0);
      std::vector<unsigned int> heightFace(currFaceId,0);
      std::vector<bool> visitedFace(currFaceId, false);
      std::vector<bool> visitedNodes(_g->vertices(), false);
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
          for(unsigned int nextEdge = _g->getVertexFirst(currNode); nextEdge <= _g->getVertexLast(currNode); ++nextEdge ){
            unsigned int nextNode = _g->getEdgeTgt(nextEdge);
            if(!visitedNodes[nextNode] && !usedEdge[nextEdge]){
               // std::cerr << "bloqueo arista " << getEdgeSrc(nextEdge) << "->" << getEdgeTgt(nextEdge) << std::endl;
              usedEdge[nextEdge] = usedEdge[_g->getEdgeCmp(nextEdge)] = true;
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


            unsigned int cmpEdge = _g->getEdgeCmp(currFaceEdge);
            if(usedEdge[currFaceEdge]){
              currFaceEdge = _g->getNextEdgeCCW(cmpEdge);
              // std::cerr << "la arista " << getEdgeSrc(currFaceEdge) << " " << getEdgeTgt(currFaceEdge) <<  " esta bloqueada" << std::endl;
              continue;
            }

            unsigned int neighborFace = edgeFace[cmpEdge];
            if(!visitedFace[neighborFace]){
              // std::cerr << "la arista " << getEdgeSrc(currFaceEdge) << " " << getEdgeTgt(currFaceEdge) <<  " se agrega" << std::endl;
              belongsToT[cmpEdge] = belongsToT[_g->getEdgeCmp(cmpEdge)] = false;
              usedEdge[cmpEdge] = usedEdge[_g->getEdgeCmp(cmpEdge)] = true;
              qFaces.push(neighborFace);
              visitedFace[neighborFace] = true;
              heightFace[neighborFace] = currHeight +1;
            }
            currFaceEdge = _g->getNextEdgeCCW(cmpEdge);
          } while(currFaceEdge!=initFaceEdge);

        }

      }

      return belongsToT;
  }
};


class LimitedDfsSTGenerator: public GraphSpanningTreeGenerator{
public:
  LimitedDfsSTGenerator(unsigned int maxHeight) : _maxHeight(maxHeight){};

  vector<bool> getSpanningTree() override {
    unsigned int initEdge = _initialEdge;
    unsigned int maxHeight = _maxHeight;
    unsigned int init = _g->getEdgeSrc(initEdge);
    vector<bool> belongsToT(2*_g->edges(),false);
    std::vector<bool> visited(_g->vertices(),false);
    // std::vector<unsigned int> parentEdge(vertices(), -1);
    std::vector<unsigned int> height(_g->vertices(), 0);
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
      unsigned int parentHeight = parentEdge==-1 ? 0 : height[_g->getEdgeSrc(parentEdge)];
      //std::cerr << node << std::endl;
      s.pop();
      if(!visited[node] && parentHeight<maxHeight){
        // std::cerr << "nodo "<<node<<std::endl;
        visited[node] = true;
        height[node] = parentHeight+1;
        //s.push(nextNode);
        if(parentEdge!=-1) belongsToT[parentEdge] = belongsToT[_g->getEdgeCmp(parentEdge)] = true;
        for(unsigned int nextEdge = _g->getVertexFirst(node); nextEdge <= _g->getVertexLast(node); ++nextEdge ){
          unsigned int nextNode = _g->getEdgeTgt(nextEdge);
          // parentEdge[nextNode] = nextEdge;
          s.push(pair<unsigned int, unsigned int>(nextNode,nextEdge));
        }
      }else if(!visited[node]){
        frontierNodes.push(pair<unsigned int, unsigned int>(node, parentEdge));
      }
    }

    // if(frontierNodes.empty()){
    //   std::cerr << "todos durante dfs" << std::endl;
    // }
    while(!frontierNodes.empty()){
      unsigned int node = frontierNodes.front().first;
      //unsigned int nodeParentEdge = parentEdge[node];
      unsigned int parentEdge = frontierNodes.front().second;
      unsigned int parentNode = _g->getEdgeSrc(parentEdge);
      //std::cerr << node << std::endl;
      frontierNodes.pop();
      if(!visited[node]){
        visited[node] = true;
        //height[node] = height[parentNode]+1;
        //s.push(nextNode);
        belongsToT[parentEdge] = belongsToT[_g->getEdgeCmp(parentEdge)] = true;
        for(unsigned int nextEdge = _g->getVertexFirst(node); nextEdge <= _g->getVertexLast(node); ++nextEdge ){
          unsigned int nextNode = _g->getEdgeTgt(nextEdge);
          // parentEdge[nextNode] = nextEdge;
          frontierNodes.push(pair<unsigned int, unsigned int>(nextNode, nextEdge));
        }
      }
    }


      return belongsToT;
  }
private:
  unsigned int _maxHeight;
};




#endif
