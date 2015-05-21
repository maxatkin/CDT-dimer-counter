#include "count.h"
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <algorithm>

//function for printing contents of vector
void printvec(vector<unsigned int> v){
  for(unsigned int i = 0; i < v.size() - 1; i++){
    cout << v[i] <<",";
  }
  cout << v[v.size() - 1] << endl;
}

void printdimers(vector<unsigned int> *m){
  for(unsigned int i = 0; i < m->size() - 1; i++){
//    cout << "       " <<  i->first << "       |        " << i->second << endl;
    if(i == 0) cout << m->at(i) << "+";
    else if(i == 1) cout << m->at(i) << "x" << "+";
    else cout << m->at(i) << "x^"<< i << "+";
  }
  cout << m->at(m->size() - 1) << "x^" << m->size() - 1 << ")" << endl;  
}

unsigned int fibonacci_n(unsigned int n){
  if(n == 0 || n == 1) return 1;
  //pow((1+sqrt(5))/2,n)-pow((1-sqrt(5))/2,n)
  return fibonacci_n(n-1) + fibonacci_n(n-2);
}

/*function for printing the data the creates each triangulation. The CDT skeleton data consists of 
a vector of vector, in which each vector is a spacelike layer of the CDT. The vector of a particular layer
contains numbers specifying how many down-pointing triangles are in each block of down-pointing triangles.
These blocks are necessarily separated by up-pointing triangles.*/

void printCDTskel(vector<vector<unsigned int> > cdt_skel){
  for(unsigned int i = 0; i < cdt_skel.size(); i++){
    printvec(cdt_skel[i]);
  }
  cout << endl;
}

/* partition the number n into m parts */
void partition_n_into_m(unsigned int m, unsigned int n, vector<vector<unsigned int> > *allparts, vector<unsigned int> partlist = vector<unsigned int>(0), unsigned int level = 0){
  if(level == m-1){
    partlist.push_back(n);
    allparts->push_back(partlist);
  }
  for(unsigned int i = 0; i <= n && level != m-1; i++){
    partlist.push_back(i);
    partition_n_into_m(m, n-i, allparts, partlist, level + 1);
    partlist.pop_back();
  }
}

/* partition the numbers less than or equal to n into m parts*/ 
void partition_ltn_into_m(unsigned int m, unsigned int n, vector<vector<unsigned int> > *allparts, vector<unsigned int> partlist = vector<unsigned int>(0)){
  for(unsigned int i = 1; i <= n; i++) partition_n_into_m(m, i, allparts, partlist); 
}

/* this cyclicly rotates the elements of a vector*/
void cycle_vector_once(vector<unsigned int> *temp){
  unsigned int front;
  front = temp->front();
  temp->erase(temp->begin());
  temp->push_back(front);
}

/*this removes vecotrs in next_layers which are related to one another by a cyclic permutation*/
void remove_cyclic_perms(vector<vector<unsigned int> > *next_layers){
  vector<unsigned int> temp;
  for(unsigned int i = 0; i < next_layers->size(); i++){
    temp = next_layers->at(i);
    for(unsigned int j = 0; j < temp.size() - 1; j++){
      cycle_vector_once(&temp);
      next_layers->erase(find(next_layers->begin(), next_layers->end(), temp));
    }
  }
}

/*this creates an ensemble of CDT skeleton data with each CDT of size total_size. This function requires an intial layer of the triangulation 
to already be loaded in cdt_skel*/

void create_CDT_data(vector<cdt_skeleton> *cdt_ensemble, cdt_skeleton cdt_skel, unsigned int total_size, unsigned int current_size = 0, bool second_layer = false){
  unsigned int up_pointing_tri_total = 0;
  unsigned int number_of_triangles_to_add = 0;
  vector<unsigned int> last_layer = cdt_skel.back();
  vector<vector<unsigned int> > next_layers;

  for(unsigned int i = 0; i < last_layer.size(); i++) up_pointing_tri_total += last_layer[i];
  current_size += 2*up_pointing_tri_total;
  number_of_triangles_to_add = (unsigned int)((total_size - current_size)/2);
  if(number_of_triangles_to_add == 0){
    cdt_ensemble->push_back(cdt_skel);
    return;
  }
  partition_ltn_into_m(up_pointing_tri_total, number_of_triangles_to_add, &next_layers);
  if(second_layer == true) remove_cyclic_perms(&next_layers);
  for(unsigned int i = 0; i < next_layers.size(); i++){
    cdt_skel.push_back(next_layers[i]);
    create_CDT_data(cdt_ensemble, cdt_skel, total_size, current_size, false);
    cdt_skel.pop_back();
  }
}

/*create an ensemble of CDTs skeleton data containing all CDTs of size n*/

void create_CDT_size_n(vector<cdt_skeleton> *cdt_e, unsigned int size){
  vector<unsigned int> partlist;
  cdt_skeleton cdt_skel;

  for(unsigned int i = 1; i <= size/2; i++){  
    partlist.push_back(i);
    cdt_skel.push_back(partlist);
    create_CDT_data(cdt_e, cdt_skel, size, 0, false);
    cdt_skel.pop_back();
    partlist.pop_back();
  }
}

Triangle::Triangle(){
  on_left_side = NULL;
  on_right_side = NULL;
  on_hrztl_side = NULL;
  contains_dimer = false;
}

Triangulation::Triangulation(){
  CurrentLayer = NULL;
  LastLayer = NULL;
  NextDownPointTri = NULL;
}

void Triangulation::Clear(){
  CurrentLayer = NULL;
  LastLayer = NULL;
  NextDownPointTri = NULL;
  Simplex.clear();
}

void Triangulation::StartNewLayer(){
  vector<Triangle *> *newlayer = new vector<Triangle *>;
  if(CurrentLayer != NULL){
    CurrentLayer->back()->on_right_side = CurrentLayer->front();
    CurrentLayer->front()->on_left_side = CurrentLayer->back();
  }
  Simplex.push_back(newlayer);
  LastLayer = CurrentLayer;
  CurrentLayer = Simplex.back();
}

Triangle *Triangulation::FindNextDownPointTri(){
  for(unsigned int i = 0; i < LastLayer->size(); i++){
    if((*LastLayer)[i]->orientation == DOWN && (*LastLayer)[i]->on_hrztl_side == NULL) return (*LastLayer)[i];
  }
}

void Triangulation::AddTriangle(bool orientation_){
  Triangle *CurrentTriangle = new Triangle;
  CurrentTriangle->orientation = orientation_;

  if(LastLayer != NULL){
    if(CurrentLayer->size() > 0){
      CurrentTriangle->on_left_side = CurrentLayer->back();
      if(orientation_ == UP){
        CurrentTriangle->on_hrztl_side = FindNextDownPointTri();
        CurrentTriangle->on_hrztl_side->on_hrztl_side = CurrentTriangle;
      }
      CurrentLayer->push_back(CurrentTriangle);
      CurrentTriangle->on_left_side->on_right_side = CurrentLayer->back();
    }
    else if(orientation_ == UP){
      CurrentTriangle->on_hrztl_side = FindNextDownPointTri();
      CurrentTriangle->on_hrztl_side->on_hrztl_side = CurrentTriangle;
      CurrentLayer->push_back(CurrentTriangle);  
    }
  }
  else if(orientation_ == DOWN){
    if(CurrentLayer->size() > 0){
      CurrentTriangle->on_left_side = CurrentLayer->back();
      CurrentLayer->push_back(CurrentTriangle);
      CurrentTriangle->on_left_side->on_right_side = CurrentLayer->back();
    }
    else CurrentLayer->push_back(CurrentTriangle);
  }
}

void Triangulation::Print(bool verbose = false){
  Triangle *CurrentTri;
  for(unsigned int i = 0; i < Simplex.size(); i++){
    cout << endl;
    for(unsigned int j = 0; j < Simplex[i]->size(); j++){
      CurrentTri = (*(Simplex[i]))[j];
      if(verbose){
        if(CurrentTri->orientation == DOWN) printf("(D,%p;%p,%p,%p)", CurrentTri, CurrentTri->on_left_side, CurrentTri->on_hrztl_side, CurrentTri->on_right_side);
        if(CurrentTri->orientation == UP) printf("(U,%p;%p,%p,%p)", CurrentTri, CurrentTri->on_left_side, CurrentTri->on_hrztl_side, CurrentTri->on_right_side);
      }
      else{
        if(CurrentTri->orientation == DOWN) printf("(D,%p|",CurrentTri);
        if(CurrentTri->orientation == UP) printf("(U,%p|",CurrentTri);
        if(CurrentTri->contains_dimer) printf("%p)",CurrentTri->dimer_end);
        else printf("0)");
      }
    }
  }
  cout << endl;
}

unsigned int Triangulation::Create(cdt_skeleton cdt_skel){
  for(unsigned int i = 0; i < cdt_skel.size(); i++){
    StartNewLayer();
    for(unsigned int j = 0; j < cdt_skel[i].size(); j++){
      AddTriangle(UP);
      for(unsigned int k = 0; k < cdt_skel[i][j]; k++) AddTriangle(DOWN);
    }
  }

  //Add cap to top of triangulation
  StartNewLayer();
  for(unsigned int i = 0; i < LastLayer->size(); i++){
    if((*LastLayer)[i]->orientation == DOWN) AddTriangle(UP);
  }
  StartNewLayer();
  Simplex.pop_back();
}

void pick_boxes(vector<Triangle *> *Boxes, vector<vector<Triangle *> > *BoxSelection, vector<Triangle *> *current_selection = NULL, unsigned int level = 0){
  if(current_selection == NULL) current_selection = new vector<Triangle *>(0);
  if(level == Boxes->size()){
    BoxSelection->push_back(*current_selection);
    return;
  }
  pick_boxes(Boxes, BoxSelection, current_selection, level + 1);
  current_selection->push_back((*Boxes)[level]);
  pick_boxes(Boxes, BoxSelection, current_selection, level + 1);
  current_selection->pop_back();
}

void pick_wedge_boxes(vector<Triangle *> *Boxes, vector<vector<Triangle *> > *BoxSelection, vector<Triangle *> *current_selection = NULL, unsigned int level = 0){
  vector<Triangle *>::iterator dup;
 
  pick_boxes(Boxes, BoxSelection, current_selection, level);
  //find overlapping dimers and remove these configurations
  for(unsigned int i = 0; i < BoxSelection->size(); i++){
    for(unsigned int j = 0; j < BoxSelection->at(i).size(); j++){
      if(BoxSelection->at(i).at(j)->on_left_side == BoxSelection->at(i).at(j)) continue;
      dup = find(BoxSelection->at(i).begin(), BoxSelection->at(i).end(), BoxSelection->at(i).at(j)->on_left_side);
      if(dup != BoxSelection->at(i).end()){
        BoxSelection->erase(BoxSelection->begin() + i);
        i = i - 1;
        break;
      }
    }
  }
}

void Triangulation::FindDimerBoxes(vector<Triangle *> *EmptyBoxes, unsigned int type){
  Triangle *CurrentTri;
  EmptyBoxes->clear();
  for(unsigned int i = 0; i < Simplex.size(); i++){
    for(unsigned int j = 0; j < Simplex[i]->size(); j++){
      CurrentTri = (*(Simplex[i]))[j];
      if(CurrentTri->contains_dimer) continue;

      if(type == WEDGE && CurrentTri->orientation == CurrentTri->on_left_side->orientation){
        if(CurrentTri->on_left_side->contains_dimer == true) continue;
        EmptyBoxes->push_back(CurrentTri);        
      }
      if(CurrentTri->orientation == DOWN){      
        if(type == DIAMOND){
          if(CurrentTri->on_hrztl_side->contains_dimer == true) continue;
          EmptyBoxes->push_back(CurrentTri);
        }
        if(type == UPDOWN){
          if(CurrentTri->on_left_side->contains_dimer == true || CurrentTri->on_left_side->orientation == DOWN) continue;
          EmptyBoxes->push_back(CurrentTri);
        }
        if(type == DOWNUP){
          if(CurrentTri->on_right_side->contains_dimer == true || CurrentTri->on_right_side->orientation == DOWN) continue;
          EmptyBoxes->push_back(CurrentTri);
        }
      }
    }
  }
}

void Triangulation::SetDimers(vector<Triangle *> *BoxSelection, unsigned int type){
  for(unsigned int i = 0; i < BoxSelection->size(); i++){
    (*BoxSelection)[i]->contains_dimer = true;
    if(type == DIAMOND){
      (*BoxSelection)[i]->dimer_end = (*BoxSelection)[i]->on_hrztl_side;
      (*BoxSelection)[i]->on_hrztl_side->contains_dimer = true;
      (*BoxSelection)[i]->on_hrztl_side->dimer_end = (*BoxSelection)[i];
    }
    if(type == UPDOWN){
      (*BoxSelection)[i]->dimer_end = (*BoxSelection)[i]->on_left_side;
      (*BoxSelection)[i]->on_left_side->contains_dimer = true;
      (*BoxSelection)[i]->on_left_side->dimer_end = (*BoxSelection)[i];
    }
    if(type == DOWNUP){
      (*BoxSelection)[i]->dimer_end = (*BoxSelection)[i]->on_right_side;
      (*BoxSelection)[i]->on_right_side->contains_dimer = true;
      (*BoxSelection)[i]->on_right_side->dimer_end = (*BoxSelection)[i];
    }
    if(type == WEDGE){
      (*BoxSelection)[i]->dimer_end = (*BoxSelection)[i]->on_left_side;
      (*BoxSelection)[i]->on_left_side->contains_dimer = true;
      (*BoxSelection)[i]->on_left_side->dimer_end = (*BoxSelection)[i];
    }
  }
}

void Triangulation::ClearDimers(vector<Triangle *> *EmptyBoxes, unsigned int type){
  for(unsigned int i = 0; i < EmptyBoxes->size(); i++){
    (*EmptyBoxes)[i]->contains_dimer = false;
    if(type == DIAMOND) (*EmptyBoxes)[i]->on_hrztl_side->contains_dimer = false;
    if(type == UPDOWN) (*EmptyBoxes)[i]->on_left_side->contains_dimer = false;
    if(type == DOWNUP) (*EmptyBoxes)[i]->on_right_side->contains_dimer = false;
    if(type == WEDGE) (*EmptyBoxes)[i]->on_left_side->contains_dimer = false;
  }  
}

unsigned int Triangulation::CountDimers(){
  Triangle *CurrentTri;
  unsigned int count = 0;
  unsigned int count_edge_dimer = 0;
  for(unsigned int i = 0; i < Simplex.size(); i++){
    for(unsigned int j = 0; j < Simplex[i]->size(); j++){
      CurrentTri = (*(Simplex[i]))[j];
      if(CurrentTri->contains_dimer){
        if(CurrentTri->dimer_end == CurrentTri) count_edge_dimer++;
        else count++;
      }
    }
  }
  return count_edge_dimer + count/2;
}

unsigned int Triangulation::GenerateDimerConfigs(unsigned int no_of_dimers, vector<unsigned int> *dimer_count){
  unsigned int count = 0;
  vector<Triangle *> EmptyDiamondBoxes;
  vector<vector<Triangle *> > DiamondBoxSelection;

  vector<Triangle *> EmptyUpDownBoxes;
  vector<vector<Triangle *> > UpDownBoxSelection;

  vector<Triangle *> EmptyDownUpBoxes;
  vector<vector<Triangle *> > DownUpBoxSelection;

  vector<Triangle *> EmptyWedgeBoxes;
  vector<vector<Triangle *> > WedgeBoxSelection;

  FindDimerBoxes(&EmptyDiamondBoxes, DIAMOND);
  pick_boxes(&EmptyDiamondBoxes, &DiamondBoxSelection);

  for(unsigned int i = 0; i < DiamondBoxSelection.size(); i++){
    SetDimers(&(DiamondBoxSelection[i]), DIAMOND);
    FindDimerBoxes(&EmptyUpDownBoxes, UPDOWN);
    UpDownBoxSelection.clear();
    pick_boxes(&EmptyUpDownBoxes, &UpDownBoxSelection);

    for(unsigned int j = 0; j < UpDownBoxSelection.size(); j++){
      SetDimers(&(UpDownBoxSelection[j]), UPDOWN);
      FindDimerBoxes(&EmptyDownUpBoxes, DOWNUP);
      DownUpBoxSelection.clear();
      pick_boxes(&EmptyDownUpBoxes, &DownUpBoxSelection);

      for(unsigned int k = 0; k < DownUpBoxSelection.size(); k++){
        SetDimers(&(DownUpBoxSelection)[k], DOWNUP);

        FindDimerBoxes(&EmptyWedgeBoxes, WEDGE);
        WedgeBoxSelection.clear();
        pick_wedge_boxes(&EmptyWedgeBoxes, &WedgeBoxSelection);

        for(unsigned int l = 0; l < WedgeBoxSelection.size(); l++){
          SetDimers(&(WedgeBoxSelection)[l], WEDGE);
          if(CountDimers() == no_of_dimers) count++; 
          (*dimer_count)[CountDimers()] = (*dimer_count)[CountDimers()] + 1; 
          ClearDimers(&EmptyWedgeBoxes, WEDGE);
        }
        ClearDimers(&EmptyDownUpBoxes, DOWNUP);
      }
      ClearDimers(&EmptyUpDownBoxes, UPDOWN);
    }
    ClearDimers(&EmptyDiamondBoxes, DIAMOND);
  }
  return count;
}

int main(unsigned int argc, char **argv){
  vector<unsigned int> partlist;
  cdt_skeleton cdt_skel;
  vector<cdt_skeleton> cdt_e;
  vector<unsigned int> *dimer_count;
  Triangulation TestTri;
  unsigned int tot = 0;
  unsigned int size = 16;
  unsigned int dimer_size = 4;

  if(argc > 1) size = atoi(argv[1]);

  if(argc > 2) dimer_size = atoi(argv[2]); 
  //printf("No. of Dimers  |   No. of Configurations\n", size, dimer_size);

  for(unsigned int j = 2; j <= size; j = j + 2){ 
    create_CDT_size_n(&cdt_e, j);
    dimer_count = new vector<unsigned int>((unsigned int)(j/2) + 2);
    for(unsigned int i = 0; i < cdt_e.size(); i++){
      TestTri.Clear();
      TestTri.Create(cdt_e[i]);
      TestTri.GenerateDimerConfigs(dimer_size, dimer_count);
    }
    cdt_e.clear();
    cout << "z^"<< j/2 << "(";
    printdimers(dimer_count);
    delete dimer_count;
  }
  return 0;
}