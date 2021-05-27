#ifndef PLANAR_UTILS_COMP_HPP
#define PLANAR_UTILS_COMP_HPP

#include "../sdsl/int_vector.hpp"
#include <vector>
#include <iostream>
#include <stdlib.h>

using namespace std;
using namespace sdsl;

unsigned long long int getRandom(){
  unsigned long long a,b,c;
  a=rand();
  b=rand();
  c=rand();
  return ((a<<33) ^ (b<<1) ^ (c&1));
}

bool generateTuranMaximalGrid(vector<bool> &a, vector<bool> &b, vector<bool> &bStar,
  unsigned int n,  unsigned int removedEdges=0, unsigned int cols=1,
  unsigned int rows=1){
  if(removedEdges>2*n-3) return false;
  a.clear();
  b.clear();
  bStar.clear();
  a.reserve(7*n*rows*cols);
  b.reserve(3*n*rows*cols);
  bStar.reserve(5*n*rows*cols);
  unsigned int i;
  unsigned aIdx=0, bIdx=0, bStarIdx=0;
  b.push_back(1);
  bStar.push_back(1);
  unsigned int removeLeft = removedEdges ? getRandom()%(removedEdges+1) : 0;
  removeLeft = removeLeft>n-3 ? n-3 : removeLeft;
  unsigned int removeRight = removedEdges - removeLeft;
  removeRight = removeRight>n-2 ? n-2 : removeRight;
  removeLeft = removedEdges - removeRight;

  unsigned int toRemoveRight = removeRight;
  for(unsigned int row = 0; row<rows; ++row){
    // cerr << "starting row "<<row<<endl;
    // cerr << a.size()<<"/ "<<a.capacity()<<endl;
    if(row>0){
      a.push_back(1);
      b.push_back(1);

    }

    for(unsigned int col = 0; col<cols; ++col){
      //join with previous triangle

      unsigned int toRemoveLeft = removeLeft;
      // {([}^n-3
      for(i=0; i<n-3; ++i){

        a.push_back(1);
        b.push_back(1);
        unsigned int closedRemaining = n-3-i;
        unsigned int randomDelete = getRandom()%closedRemaining;
        if( randomDelete>= toRemoveLeft){
          // don't delete these parentheses
          a.push_back(0);
          bStar.push_back(1);
        }else{
          toRemoveLeft--;
        }

      }
      // (^2
      for(i=0; i<2; ++i){
        a.push_back(1);
        b.push_back(1);
      }

      // ]^n-3
      for(i=0; i<n-3-removeLeft; ++i){
        a.push_back(0);
        bStar.push_back(0);
      }

      if(row>0 ){
        //hacia abajo ]
        a.push_back(0);
        bStar.push_back(0);
      }

      if(col<cols-1){
        a.push_back(1);
        b.push_back(1);
      }

      //end column
    }

    for(unsigned int col = 0; col<cols; ++col){
      //join with previous
      if( row!=rows-1){
        a.push_back(0);
        bStar.push_back(1);
      }
      if(col>0){
        a.push_back(1);
        b.push_back(0);
      }
      // {[)}^n-2
      for(i=0; i<n-2; ++i){

        unsigned int closedRemaining = n-2-i;
        unsigned int randomDelete = getRandom()%closedRemaining;
        if( randomDelete>= toRemoveRight){
          // don't delete these parentheses
          a.push_back(0);
          bStar.push_back(1);
        }else{
          toRemoveRight--;
        }
        a.push_back(1);
        b.push_back(0);
      }
      // )
      a.push_back(1);
      b.push_back(0);
      // {]}^n-2
      for(i=0; i<n-2-removeRight ; ++i){
        a.push_back(0);
        bStar.push_back(0);
      }
    }

  }//end for rows
  // cerr << "ended rows" <<endl;
  for(unsigned int row = 0; row<rows-1; ++row){
    a.push_back(1);
    b.push_back(0);
  }

  b.push_back(0);
  bStar.push_back(0);
  aIdx=0;
  bIdx=1;
  bStarIdx=1;
  // cerr << "ended construction " <<endl;
  return true;
}

#endif
