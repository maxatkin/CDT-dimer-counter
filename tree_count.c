#include "tree_count.h"

void printtree(Tree *T, int level = 0){
  for(unsigned int i = 0; i < T->Root.neighbour_trees.size(); ++i){
    for(int j = 0; j < level; ++j) printf(" ");
    printf("%p\n", T->Root.neighbour_trees[i]);
    printtree(T->Root.neighbour_trees[i], level + 1);
  }
}

void Node::CopyNeighbourTrees(Tree *Neighbour){
  neighbour_trees = Neighbour->Root.neighbour_trees;
}

void Node::AddNeighbourTree(Tree *Neighbour){
  neighbour_trees.push_back(Neighbour);
}

Tree::Tree(){

}

void create_trees_size_n_fast(int n, vector<vector<Tree *> > *Trees){
  Tree *NewTree;

  for(unsigned int i = 1; i < n; ++i){
    //printf("%i,%i\n",n-i,i);
    for(int j = 0; j < Trees->at(n-i).size(); ++j){
      for(int k = 0; k < Trees->at(i).size(); ++k){
        NewTree = new Tree;
        NewTree->Root.neighbour_trees = Trees->at(n-i).at(j)->Root.neighbour_trees;
        NewTree->Root.neighbour_trees.push_back(Trees->at(i).at(k));
        Trees->at(n).push_back(NewTree);
      }
    }
  }
}

int main(){
  Tree *tmptree;
  vector<Tree *> tmp;
  vector<vector<Tree *> > Trees(40);
  tmptree = new Tree;
  tmp.push_back(tmptree);
  Trees[1] = tmp;
  for(unsigned int i = 0; i < 30; ++i){
    create_trees_size_n_fast(i, &Trees);
    printf("%i, %i\n", 2*(i-1), Trees[i].size());
  }

  printtree(Trees[13][9424]);
}