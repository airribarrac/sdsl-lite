#ifndef GRAPH_VISITOR_COMP_HPP
#define GRAPH_VISITOR_COMP_HPP

#include "Graph.hpp"
#include <vector>
#include <stack>
#include <queue>

#ifdef PRINT_TREE
bool printTree = true;
#else
bool printTree = false;
#endif

//interface for generating spanning trees
class GraphSpanningTreeGenerator{
public:
  void setGraph(Graph *g){ _g=g;}
  void setInitialEdge(unsigned int initialEdge){ _initialEdge=initialEdge;}
  unsigned int getInitialEdge(){ return _initialEdge;}
  virtual vector<bool> getSpanningTree()=0;
  void printTreeEdges(vector<bool> &belongsToT){
	  if(!printTree) return;
	  std::cout << "---------BEGIN PRINT TREE----------"<<std::endl;
	  unsigned int edgesCnt = 0;
	for(unsigned int e=0; e<2*_g->edges();e++){
		if(belongsToT[e]){
			std::cout << _g->getEdgeSrc(e)<<" "<<_g->getEdgeTgt(e)<<endl;
			edgesCnt++;
		}
		if(belongsToT[e] && !belongsToT[_g->getEdgeCmp(e)]){
			std::cerr <<"WARNING: edge "<<e<<" is in tree, but not its cmp"<<std::endl;
		}

	}
	  std::cout << "----------END PRINT TREE-----------"<<std::endl;
	  std::cout << "#edges="<<edgesCnt<<std::endl;
  }

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

class DfsPrimalBadSTGenerator : public GraphSpanningTreeGenerator{
public:


  vector<bool> getSpanningTree() override {
    unsigned int init = _g->getEdgeSrc(_initialEdge);
    vector<bool> belongsToT(2*_g->edges(),false);

    std::vector<bool> visited(_g->vertices(),false);
    std::stack<unsigned int> s;

    for(unsigned int nextEdge = _g->getVertexFirst(init); nextEdge <= _g->getVertexLast(init); ++nextEdge ){
      unsigned int nextNode = _g->getEdgeTgt(nextEdge);

      if(!visited[nextNode]){
        s.push(nextEdge);
      }
    }


    visited[init] = true;

    while(!s.empty()){
      unsigned int edge = s.top();
      s.pop();
      unsigned int node = _g->getEdgeTgt(edge);
      if(visited[node]) continue;
      visited[node]=true;
      belongsToT[edge] = belongsToT[_g->getEdgeCmp(edge)] = true;
      // std::cerr << edge <<" "<<_g->getEdgeCmp(edge)<<" "<<_g->getVertexFirst(node)<<" "<<_g->getVertexLast(node)<<std::endl;


      // unsigned int nextEdge = _g->getEdgeCmp(edge);
      // do{
      //   unsigned int nextNode = _g->getEdgeTgt(nextEdge);
      //   if(!visited[nextNode]){
      //     // std::cerr << "pusheo nodo "<<nextNode<<std::endl;
      //     // std::cerr << "la arista "<<_g->getEdgeCmp(nextEdge)<<std::endl;
      //     s.push(nextEdge);
      //   }
      //
      //   if(nextEdge==_g->getVertexLast(node)){
      //     nextEdge=_g->getVertexFirst(node);
      //   }else{
      //     nextEdge++;
      //   }
      // }while(nextEdge!=_g->getEdgeCmp(edge));

      for(unsigned int nextEdge = _g->getVertexFirst(node); nextEdge <= _g->getVertexLast(node); ++nextEdge ){
        unsigned int nextNode = _g->getEdgeTgt(nextEdge);

        if(!visited[nextNode]){
          // std::cerr << "pusheo nodo "<<nextNode<<std::endl;
          // std::cerr << "la arista "<<_g->getEdgeCmp(nextEdge)<<std::endl;
          s.push(nextEdge);
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
   
    std::cerr << "ROOTED AT "<<init<<std::endl;
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

printTreeEdges(belongsToT);	
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



printTreeEdges(belongsToT);	
      return belongsToT;
  }
private:
  unsigned int _maxHeight;
};


class DegHeurSTGenerator: public GraphSpanningTreeGenerator{
public:

  vector<bool> getSpanningTree() override {
    vector<bool> belongsToT(2*_g->edges(),false);
    unsigned int init = _g->getEdgeSrc(_initialEdge);
    std::vector<bool> visited(_g->vertices(),false);
    std::priority_queue<std::pair<unsigned int, unsigned int>> pq;

    std::vector<unsigned int> externalEdges = _g->externalFace();

    for(int i=0; i<externalEdges.size()&&false; i++){
      unsigned int currEdge = externalEdges[i];
      unsigned int currNode = _g->getEdgeSrc(currEdge);
      unsigned int currDegree = _g->getVertexLast(currNode)-_g->getVertexFirst(currNode);
      unsigned int initDegree = _g->getVertexLast(init)-_g->getVertexFirst(init);

      if(currDegree>initDegree){
        init = currNode;
      }
    }

    unsigned int initEdgeCmp = _g->getVertexFirst(init);
    unsigned int initEdge = _g->getEdgeCmp(initEdgeCmp);
	
    std::cerr << "ROOTED AT"<<init<<std::endl;
    pq.push(std::make_pair(_g->getVertexLast(init)-_g->getVertexFirst(init), initEdge));
    visited[init] = true;
    initEdge = initEdgeCmp;
    _initialEdge=initEdge;
    while(!pq.empty()){
      auto pa = pq.top();
      unsigned int edge = pa.second;
      unsigned int node = _g->getEdgeTgt(edge);
      // std::cerr << node << std::endl;
      pq.pop();
      for(unsigned int nextEdge = _g->getVertexFirst(node); nextEdge <= _g->getVertexLast(node); ++nextEdge ){
        unsigned int nextNode = _g->getEdgeTgt(nextEdge);
        if(!visited[nextNode]){
          unsigned int nextNodeDeg = _g->getVertexLast(nextNode)-_g->getVertexFirst(nextNode);
          visited[nextNode] = true;
          pq.push(std::make_pair(nextNodeDeg,nextEdge));
          belongsToT[nextEdge] = belongsToT[_g->getEdgeCmp(nextEdge)] = true;
        }
      }
    }


printTreeEdges(belongsToT);	

    return belongsToT;
  }
private:
  unsigned int _maxHeight;
};

class DegHeurWithLSGenerator : public DegHeurSTGenerator{
  vector<bool> getSpanningTree() override{

    vector<bool> belongsToT=DegHeurSTGenerator::getSpanningTree();
    unsigned int init = _g->getEdgeSrc(_initialEdge);


    vector<unsigned int> parent(_g->vertices(),-1);
    vector<unsigned int> level(_g->vertices(),0);
    vector<bool> visited(_g->vertices(),false);
    vector<bool> isLeaf(_g->vertices(),true);
    isLeaf[init]=false;
    vector<unsigned int> parentEdge(_g->vertices(),-1);
    queue<unsigned int> q;
    q.push(init);
    visited[init]=true;

    while(!q.empty()){
      unsigned int u = q.front();
      q.pop();

      for(unsigned int nextEdge = _g->getVertexFirst(u); nextEdge <= _g->getVertexLast(u); ++nextEdge ){
        unsigned int v = _g->getEdgeTgt(nextEdge);
        if(belongsToT[nextEdge] && !visited[v]){
          isLeaf[u]=false;
          parent[v]=u;
          level[v]=level[u]+1;
          visited[v]=true;
          parentEdge[v]=nextEdge;
          q.push(v);
        }
      }

    }
    unsigned int improved = 0;
    for(unsigned int u=0; u<_g->vertices(); u++){
      if(u==init) continue;
      for(unsigned int nextEdge = _g->getVertexFirst(u); nextEdge <= _g->getVertexLast(u); ++nextEdge ){
        unsigned int v =  _g->getEdgeTgt(nextEdge);
        //v is a neighbor parent closer to the root
        if(level[u]>level[v]+1 && isLeaf[u] && !isLeaf[v]){
          parent[u] = v;
          level[u] = level[v]+1;
          belongsToT[parentEdge[u]]=belongsToT[_g->getEdgeCmp(parentEdge[u])]=false;
          parentEdge[u] = nextEdge;
          belongsToT[nextEdge]=belongsToT[_g->getEdgeCmp(nextEdge)]=true;
          improved++;
        }
      }
    }
    cerr << improved << " nodes were improved " << endl;

printTreeEdges(belongsToT);	
    return belongsToT;
  }
};





class BfsDegSTGenerator: public GraphSpanningTreeGenerator{
public:

  vector<bool> getSpanningTree() override {
    vector<bool> belongsToT(2*_g->edges(),false);
    unsigned int init = _g->getEdgeSrc(_initialEdge);
    std::vector<bool> visited(_g->vertices(),false);
//    std::priority_queue<std::pair<unsigned int, unsigned int>> pq;
    std::priority_queue<std::pair<std::pair<int,unsigned int>, unsigned int> > pq;
    //					   lev, degree,        edge

    std::vector<unsigned int> externalEdges = _g->externalFace();


    unsigned int initEdgeCmp = _g->getVertexFirst(init);
    unsigned int initEdge = _g->getEdgeCmp(initEdgeCmp);
	
    std::cerr << "ROOTED AT"<<init<<std::endl;
    pq.emplace(std::make_pair(0,_g->getVertexLast(init)-_g->getVertexFirst(init)), initEdge);
//    pq.push(std::make_pair(_g->getVertexLast(init)-_g->getVertexFirst(init), initEdge));
    visited[init] = true;
    initEdge = initEdgeCmp;
    _initialEdge=initEdge;
    while(!pq.empty()){
      auto pa = pq.top();
      unsigned int edge = pa.second;
      unsigned int node = _g->getEdgeTgt(edge);
      unsigned int level = (-pa.first.first);
      //std::cerr << node <<" degree="<<pa.first.second<<" level="<<level << std::endl;
      pq.pop();
      for(unsigned int nextEdge = _g->getVertexFirst(node); nextEdge <= _g->getVertexLast(node); ++nextEdge ){
        unsigned int nextNode = _g->getEdgeTgt(nextEdge);
        if(!visited[nextNode]){
          unsigned int nextNodeDeg = _g->getVertexLast(nextNode)-_g->getVertexFirst(nextNode);
          visited[nextNode] = true;
          pq.emplace(std::make_pair(-(level+1),nextNodeDeg),nextEdge);
          belongsToT[nextEdge] = belongsToT[_g->getEdgeCmp(nextEdge)] = true;
        }
      }
    }


printTreeEdges(belongsToT);	

    return belongsToT;
  }
private:
  unsigned int _maxHeight;
};



class LSMoveChildrenGenerator : public DegHeurSTGenerator{
  vector<bool> getSpanningTree() override{

    vector<bool> belongsToT=DegHeurSTGenerator::getSpanningTree();
    unsigned int init = _g->getEdgeSrc(_initialEdge);


    vector<unsigned int> parent(_g->vertices(),-1);
    vector<unsigned int> level(_g->vertices(),0);
    vector<unsigned int> children(_g->vertices(),0);
    vector<bool> visited(_g->vertices(),false);
    vector<bool> isLeaf(_g->vertices(),true);
    isLeaf[init]=false;
    vector<unsigned int> parentEdge(_g->vertices(),-1);
    queue<unsigned int> q;
    q.push(init);
    visited[init]=true;

    while(!q.empty()){
      unsigned int u = q.front();
      q.pop();

      for(unsigned int nextEdge = _g->getVertexFirst(u); nextEdge <= _g->getVertexLast(u); ++nextEdge ){
        unsigned int v = _g->getEdgeTgt(nextEdge);
        if(belongsToT[nextEdge] && !visited[v]){
          isLeaf[u]=false;
	  children[u]++;
          parent[v]=u;
          level[v]=level[u]+1;
          visited[v]=true;
          parentEdge[v]=nextEdge;
          q.push(v);
        }
      }

    }
    unsigned int leaves = 0;
    for(int rep=0; rep<500; rep++){
    priority_queue<pair<unsigned int,unsigned int> > pq;
    for(unsigned int i=0; i<_g->vertices(); i++){
	pq.emplace(-children[i],i);
	//if(children[i]==0) leaves++;
    }

    
    while(!pq.empty()){
      unsigned int u = pq.top().second;
      pq.pop();
      if(children[u]==0) continue;
      if(u==init) continue;
      vector<pair<int,int> >  replaces;
      unsigned int numChildren = 0;
      for(unsigned int nextEdge = _g->getVertexFirst(u); nextEdge <= _g->getVertexLast(u); ++nextEdge ){
        unsigned int v = _g->getEdgeTgt(nextEdge);
	if(parent[v]!=u) continue;
	if(v==init) continue;
        if(belongsToT[nextEdge]){
		numChildren++;
		// v es su hijo
		bool addedChild = false;
		unsigned int bestEdge = -1;
		unsigned int bestEdgeDegree=0;
		for(unsigned int vEdge=_g->getVertexFirst(v); vEdge<= _g->getVertexLast(v); ++vEdge){
			unsigned int w = _g->getEdgeTgt(vEdge);
			//no queremos agregarle hijos a una hoja
			if(u==w) continue;
			if(children[w]==0) continue;
			if(belongsToT[vEdge]) continue;
			bool canBeParent = true;
			unsigned int ancestor = parent[w];
			while(ancestor!=init && canBeParent){
				if(ancestor==v) canBeParent=false;
				ancestor=parent[ancestor];
			}
			if(!canBeParent) continue;
			if(printTree){
				cerr << "removed edge ("<<_g->getEdgeSrc(nextEdge)<<","<<_g->getEdgeTgt(nextEdge)<<")"<<endl;

				cerr << "added edge ("<<_g->getEdgeSrc(vEdge)<<","<<_g->getEdgeTgt(vEdge)<<")"<<endl;
				cerr << "old parent of "<<v<<" was "<<parent[v]<<" not it is "<<w<<endl;

			}
			unsigned int thisNodeDeg = _g->getVertexLast(w)-_g->getVertexFirst(w)+1;
			if(thisNodeDeg>bestEdgeDegree){

			bestEdgeDegree = thisNodeDeg;
			bestEdge = vEdge;

			addedChild=true;
			}
			//break;
		}
		if(!addedChild) break;
		unsigned int vEdge = bestEdge;
		unsigned int w = _g->getEdgeTgt(vEdge);
			belongsToT[nextEdge]=belongsToT[_g->getEdgeCmp(nextEdge)]=false;
			belongsToT[vEdge]=belongsToT[_g->getEdgeCmp(vEdge)]=true;

			replaces.emplace_back(nextEdge,vEdge);
			children[u]--;
			children[w]++;	
			parent[v]=w;
        }
      }
      if(replaces.size()!=numChildren){
	//this node can become leaf
	for(unsigned int i=0; i<replaces.size(); i++){
		unsigned int nextEdge = replaces[i].first;
		unsigned int vEdge = replaces[i].second;
		if(printTree){

			cerr << "undoing replace edge "<<nextEdge<<" with "<<vEdge<<endl;
		}
		unsigned int v = _g->getEdgeTgt(nextEdge);
		unsigned int w = _g->getEdgeTgt(vEdge);
			children[u]++;
			children[w]--;	
			parent[v]=u;
			belongsToT[nextEdge]=belongsToT[_g->getEdgeCmp(nextEdge)]=true;
			belongsToT[vEdge]=belongsToT[_g->getEdgeCmp(vEdge)]=false;
	}

      }

    }


    visited.assign(_g->vertices(),false);
    q.push(init);
    visited[init]=true;
    unsigned int visNodes =0;


    
    while(!q.empty()){
      unsigned int u = q.front();
      q.pop();

      for(unsigned int nextEdge = _g->getVertexFirst(u); nextEdge <= _g->getVertexLast(u); ++nextEdge ){
        unsigned int v = _g->getEdgeTgt(nextEdge);
        if(belongsToT[nextEdge] && !visited[v]){
          isLeaf[u]=false;
          parent[v]=u;
          level[v]=level[u]+1;
          visited[v]=true;
          parentEdge[v]=nextEdge;
          q.push(v);
        }
      }

    }
    unsigned int improved = 0;
    for(unsigned int u=0; u<_g->vertices(); u++){
      if(u==init) continue;
      for(unsigned int nextEdge = _g->getVertexFirst(u); nextEdge <= _g->getVertexLast(u); ++nextEdge ){
        unsigned int v =  _g->getEdgeTgt(nextEdge);
        //v is a neighbor parent closer to the root
        if(level[u]>level[v]+1 && isLeaf[u] && !isLeaf[v]){
	  children[parent[u]]--;
	  children[v]++;
          parent[u] = v;
          level[u] = level[v]+1;
          belongsToT[parentEdge[u]]=belongsToT[_g->getEdgeCmp(parentEdge[u])]=false;
          parentEdge[u] = nextEdge;
          belongsToT[nextEdge]=belongsToT[_g->getEdgeCmp(nextEdge)]=true;
          improved++;
        }
      }
    }

    unsigned int newLeaves =0;
    for(unsigned int i=0; i<_g->vertices(); i++){
//	pq.emplace(-children[i],i);
	if(children[i]==0) newLeaves++;
    }
    cerr <<"old leaves="<<leaves<<" new leaves="<<newLeaves<<endl;
    if(leaves == newLeaves) break;
    leaves = newLeaves;
    }

    

//    cerr <<"visNodes is "<<visNodes<<endl;
  //  cerr << improved << " nodes were improved " << endl;
    cerr <<"init edge is "<<init<<endl;
printTreeEdges(belongsToT);	
    return belongsToT;
  }
};



class DegHeurDualSTGenerator: public GraphSpanningTreeGenerator{
public:

  vector<bool> getSpanningTree() override {
    unsigned int initEdge = _initialEdge;
    vector<bool> belongsToT(2*_g->edges(),true);
    //face to which edge i belongs
    std::vector<unsigned int> edgeFace(2*_g->edges(), -1);
    //an edge which belongs to face i
    std::vector<unsigned int> faceEdge;
    std::vector<unsigned int> faceDegree;
    unsigned int currFaceId = 0;

    for(unsigned int edge=0; edge<2*_g->edges(); edge++){

      if(edgeFace[edge]==-1){
        faceEdge.push_back(edge);
	faceDegree.push_back(1);
        unsigned int currFaceEdge = edge;
        do{
          edgeFace[currFaceEdge] = currFaceId;
          unsigned int cmp = _g->getEdgeCmp(currFaceEdge);
          currFaceEdge = _g->getNextEdgeCCW(cmp);
	  faceDegree.back()++;
        }while(currFaceEdge!=edge);
        currFaceId++;
      }

    }
    std::priority_queue<std::pair<unsigned int,unsigned int >> q;
    std::vector<bool> visited(currFaceId, false);
    unsigned initFace = edgeFace[initEdge];

    q.emplace(faceDegree[initFace],initFace);
    visited[initFace] = true;

    while(!q.empty()){
      unsigned int currFace = q.top().second;
      q.pop();

      unsigned int initFaceEdge = faceEdge[currFace];
      unsigned int currFaceEdge = initFaceEdge;

      do {
        unsigned int cmpEdge = _g->getEdgeCmp(currFaceEdge);
        unsigned int neighborFace = edgeFace[cmpEdge];
        if(!visited[neighborFace]){
          belongsToT[cmpEdge] = belongsToT[_g->getEdgeCmp(cmpEdge)] = false;
          q.emplace(faceDegree[neighborFace],neighborFace);
          visited[neighborFace] = true;
        }
        currFaceEdge = _g->getNextEdgeCCW(cmpEdge);
      } while(currFaceEdge!=initFaceEdge);

    }

    return belongsToT;
  }
};


class LongTreeGenerator : public DfsPrimalBadSTGenerator{
  vector<bool> getSpanningTree() override{

    vector<bool> belongsToT=DfsPrimalBadSTGenerator::getSpanningTree();
    unsigned int init = _g->getEdgeSrc(_initialEdge);

    unsigned int numLeaves = 0;

    vector<unsigned int> parent(_g->vertices(),-1);
    vector<unsigned int> level(_g->vertices(),0);
    vector<bool> visited(_g->vertices(),false);
    vector<bool> isLeaf(_g->vertices(),true);
    vector<unsigned int> numChildren(_g->vertices(),0);
    isLeaf[init]=false;
    vector<unsigned int> parentEdge(_g->vertices(),-1);
    queue<unsigned int> q;
    q.push(init);
    visited[init]=true;

    while(!q.empty()){
      unsigned int u = q.front();
      q.pop();

      for(unsigned int nextEdge = _g->getVertexFirst(u); nextEdge <= _g->getVertexLast(u); ++nextEdge ){
        unsigned int v = _g->getEdgeTgt(nextEdge);
        if(belongsToT[nextEdge] && !visited[v]){
          isLeaf[u]=false;
	  numChildren[u]++;
          parent[v]=u;
          level[v]=level[u]+1;
          visited[v]=true;
          parentEdge[v]=nextEdge;
          q.push(v);
        }
      }

    }
    unsigned int improved = 0;

    for(unsigned int u=0; u<_g->vertices(); u++){
	numLeaves+=(numChildren[u]==0);
//	cerr << "padre de "<<u<<" ->"<<parent[u] <<endl;
    }
    cerr << "initial number of leaves "<<numLeaves <<endl;
  //  cerr << "la raiz es "<<init <<endl;
    unsigned itersBadLuck = 0;

    for(unsigned int j=0;j<200000;j++){
      if(itersBadLuck==10) break;
      if(j%1000==0)cerr << itersBadLuck << " " << numLeaves <<endl;
      int u = rand()%_g->vertices();
      while(u==init || numChildren[u]>1) u= rand()%_g->vertices();
      unsigned int parentU = parent[u];
      // avanzamos hasta llegar a un fork
      while( numChildren[parentU]<=1 && parentU!=init){
	u = parentU;
	parentU = parent[u];
      }
      

//    for(unsigned int u=0; u<_g->vertices(); u++){
      if(u==init){ itersBadLuck++;continue;}
    //  cerr << "escojo nodo " <<u<<endl;
      for(unsigned int nextEdge = _g->getVertexFirst(u); nextEdge <= _g->getVertexLast(u); ++nextEdge ){
        unsigned int v =  _g->getEdgeTgt(nextEdge);
	if(v==parent[u]) continue;
	unsigned int pp = v;
	while(pp!=init && pp!=u) pp=parent[pp];
	if(pp == u) break;
        //v is is a leaf, lets try to cover it by putting u under it
	if(numChildren[v]==0){
		numChildren[parentU]--;
		if(numChildren[parentU]==0) numLeaves++;
		numChildren[v]++;
		numLeaves --;
//		cerr << " nuevo padre va a ser " << v <<endl;
		parent[u] = v;
		belongsToT[parentEdge[u]]=belongsToT[_g->getEdgeCmp(parentEdge[u])]=false;
          	parentEdge[u] = nextEdge;
          	itersBadLuck=0;
		belongsToT[nextEdge]=belongsToT[_g->getEdgeCmp(nextEdge)]=true;
		break;


	}
      }
      itersBadLuck++;

    }
    cerr << "counted number of leaves by pred " <<numLeaves<<endl;
    numLeaves = 0;
    for(unsigned int u=0; u<_g->vertices(); u++){
	numLeaves+=(numChildren[u]==0);
    }
    cerr << "last number Of leaves "<<numLeaves <<endl;
    cerr << improved << " nodes were improved " << endl;

printTreeEdges(belongsToT);	
    return belongsToT;
  }
};

#endif
