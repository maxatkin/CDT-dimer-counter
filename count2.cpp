#include "count.h"
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <algorithm>

void printvec(vector<int> v){
  for(int i = 0; i < v.size() - 1; i++){
    cout << v[i] <<",";
  }
  cout << v[v.size() - 1] << endl;
}

void printmap(map<int, int> *m){
  map<int, int>::iterator i;
  for(i = m->begin(); i != m->end(); i++){
    cout << "       " <<  i->first << "       |        " << i->second << endl;
  }
}


void printCDTskel(vector<vector<int> > cdt_skel){
  for(int i = 0; i < cdt_skel.size(); i++){
    printvec(cdt_skel[i]);
  }
  cout << endl;
}

void partition_n_into_m(int m, int n, vector<vector<int> > *allparts, vector<int> partlist = vector<int>(0), int level = 0){
  if(level == m-1){
    partlist.push_back(n);
    allparts->push_back(partlist);
  }
  for(int i = 0; i <= n && level != m-1; i++){
    partlist.push_back(i);
    partition_n_into_m(m, n-i, allparts, partlist, level + 1);
    partlist.pop_back();
  }
}

void partition_ltn_into_m(int m, int n, vector<vector<int> > *allparts, vector<int> partlist = vector<int>(0)){
  for(int i = 1; i <= n; i++) partition_n_into_m(m, i, allparts, partlist); 
}

void cycle_vector_once(vector<int> *temp){
  int front;
  front = temp->front();
  temp->erase(temp->begin());
  temp->push_back(front);
}

void remove_cyclic_perms(vector<vector<int> > *next_layers){
  vector<int> temp;
  for(int i = 0; i < next_layers->size(); i++){
    temp = next_layers->at(i);
    for(int j = 0; j < temp.size() - 1; j++){
      cycle_vector_once(&temp);
      next_layers->erase(find(next_layers->begin(), next_layers->end(), temp));
    }
  }
}

void create_CDT_data(vector<cdt_skeleton> *cdt_ensemble, cdt_skeleton cdt_skel, int total_size, int current_size = 0, bool second_layer = false){
  int up_pointing_tri_total = 0;
  int number_of_triangles_to_add = 0;
  vector<int> last_layer = cdt_skel.back();
  vector<vector<int> > next_layers;

  for(int i = 0; i < last_layer.size(); i++) up_pointing_tri_total += last_layer[i];
  current_size += 2*up_pointing_tri_total;
  number_of_triangles_to_add = (int)((total_size - current_size)/2);
  if(number_of_triangles_to_add == 0){
    cdt_ensemble->push_back(cdt_skel);
    return;
  }
  partition_ltn_into_m(up_pointing_tri_total, number_of_triangles_to_add, &next_layers);
  if(second_layer == true) remove_cyclic_perms(&next_layers);
  for(int i = 0; i < next_layers.size(); i++){
    cdt_skel.push_back(next_layers[i]);
    create_CDT_data(cdt_ensemble, cdt_skel, total_size, current_size, false);
    cdt_skel.pop_back();
  }
}

void create_CDT_size_n(vector<cdt_skeleton> *cdt_e, int size){
  vector<int> partlist;
  cdt_skeleton cdt_skel;

  for(int i = 1; i <= size/2; i++){  
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
  for(int i = 0; i < LastLayer->size(); i++){
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
  for(int i = 0; i < Simplex.size(); i++){
    cout << endl;
    for(int j = 0; j < Simplex[i]->size(); j++){
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

int Triangulation::Create(cdt_skeleton cdt_skel){
  for(int i = 0; i < cdt_skel.size(); i++){
    StartNewLayer();
    for(int j = 0; j < cdt_skel[i].size(); j++){
      AddTriangle(UP);
      for(int k = 0; k < cdt_skel[i][j]; k++) AddTriangle(DOWN);
    }
  }

  //Add cap to top of triangulation
  StartNewLayer();
  for(int i = 0; i < LastLayer->size(); i++){
    if((*LastLayer)[i]->orientation == DOWN) AddTriangle(UP);
  }
  StartNewLayer();
  Simplex.pop_back();
}

void pick_boxes(vector<Triangle *> *Boxes, vector<vector<Triangle *> > *BoxSelection, vector<Triangle *> *current_selection = NULL, int level = 0){
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

void pick_wedge_boxes(vector<Triangle *> *Boxes, vector<vector<Triangle *> > *BoxSelection, vector<Triangle *> *current_selection = NULL, int level = 0){
  vector<Triangle *>::iterator dup;
 
  pick_boxes(Boxes, BoxSelection, current_selection, level);
  //find overlapping dimers and remove these configurations
  for(int i = 0; i < BoxSelection->size(); i++){
    for(int j = 0; j < BoxSelection->at(i).size(); j++){
      dup = find(BoxSelection->at(i).begin(), BoxSelection->at(i).end(), BoxSelection->at(i).at(j)->on_left_side);
      if(dup != BoxSelection->at(i).end()){
        BoxSelection->erase(BoxSelection->begin() + i);
        i = i - 1;
      }
    }
  }
}

void Triangulation::FindDimerBoxes(vector<Triangle *> *EmptyBoxes, int type){
  Triangle *CurrentTri;
  EmptyBoxes->clear();
  for(int i = 0; i < Simplex.size(); i++){
    for(int j = 0; j < Simplex[i]->size(); j++){
      CurrentTri = (*(Simplex[i]))[j];
      if(CurrentTri->contains_dimer) continue;

      if(type == WEDGE && CurrentTri->orientation == CurrentTri->on_left_side->orientation){
        if(CurrentTri->on_left_side->contains_dimer == true || CurrentTri == CurrentTri->on_left_side) continue;
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

void Triangulation::SetDimers(vector<Triangle *> *BoxSelection, int type){
  for(int i = 0; i < BoxSelection->size(); i++){
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

void Triangulation::ClearDimers(vector<Triangle *> *EmptyBoxes, int type){
  for(int i = 0; i < EmptyBoxes->size(); i++){
    (*EmptyBoxes)[i]->contains_dimer = false;
    if(type == DIAMOND) (*EmptyBoxes)[i]->on_hrztl_side->contains_dimer = false;
    if(type == UPDOWN) (*EmptyBoxes)[i]->on_left_side->contains_dimer = false;
    if(type == DOWNUP) (*EmptyBoxes)[i]->on_right_side->contains_dimer = false;
    if(type == WEDGE) (*EmptyBoxes)[i]->on_left_side->contains_dimer = false;
  }  
}

int Triangulation::CountDimers(){
  Triangle *CurrentTri;
  int count = 0;
  for(int i = 0; i < Simplex.size(); i++){
    for(int j = 0; j < Simplex[i]->size(); j++){
      CurrentTri = (*(Simplex[i]))[j];
      if(CurrentTri->contains_dimer) count++;
    }
  }
  return count/2;
}

int Triangulation::GenerateDimerConfigs(int no_of_dimers, map <int, int> *dimer_count){
  int count = 0;
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
//  Print(true);

  for(int i = 0; i < DiamondBoxSelection.size(); i++){
    SetDimers(&(DiamondBoxSelection[i]), DIAMOND);
    FindDimerBoxes(&EmptyUpDownBoxes, UPDOWN);
    UpDownBoxSelection.clear();
    pick_boxes(&EmptyUpDownBoxes, &UpDownBoxSelection);

    for(int j = 0; j < UpDownBoxSelection.size(); j++){
      SetDimers(&(UpDownBoxSelection[j]), UPDOWN);
      FindDimerBoxes(&EmptyDownUpBoxes, DOWNUP);
      DownUpBoxSelection.clear();
      pick_boxes(&EmptyDownUpBoxes, &DownUpBoxSelection);

      for(int k = 0; k < DownUpBoxSelection.size(); k++){
        SetDimers(&(DownUpBoxSelection)[k], DOWNUP);

//        if(CountDimers() == no_of_dimers) count++; 

        FindDimerBoxes(&EmptyWedgeBoxes, WEDGE);
        WedgeBoxSelection.clear();
        pick_wedge_boxes(&EmptyWedgeBoxes, &WedgeBoxSelection);
/*
        printf("---------------\n");
        Print();
        printf("---------------\n");

        for(int ii = 0; ii < EmptyWedgeBoxes.size(); ii++){
          printf("%p\n", EmptyWedgeBoxes[ii]);
        }
        printf("gg: %i, %i\n", EmptyWedgeBoxes.size(), WedgeBoxSelection.size());
        for(int ih = 0; ih < WedgeBoxSelection.size(); ih++){
          for(int jh = 0; jh < WedgeBoxSelection[ih].size(); jh++){
            printf("%p,", WedgeBoxSelection[ih][jh]);
          }
        cout << endl;
        }
        getchar();
*/
//        cout << endl;
        for(int l = 0; l < WedgeBoxSelection.size(); l++){
          SetDimers(&(WedgeBoxSelection)[l], WEDGE);
//          Print();
//          getchar();
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

int main(int argc, char **argv){
  vector<int> partlist;
  cdt_skeleton cdt_skel;
  vector<cdt_skeleton> cdt_e;
  map <int,int> dimer_count;
  Triangulation TestTri;
  int tot = 0;
  int size = 16;
  int dimer_size = 4;
/*
  create_CDT_size_n(&cdt_e, size);
  printCDTskel(cdt_e[85]);
  TestTri.Clear();
  TestTri.Create(cdt_e[85]);
  TestTri.GenerateDimerConfigs(dimer_size);
*/
  if(argc > 1) size = atoi(argv[1]);
  create_CDT_size_n(&cdt_e, size);

  if(argc > 2) dimer_size = atoi(argv[2]); 
  printf("No. of Dimers  |   No. of Configurations\n", size, dimer_size);

  for(int i = 0; i < cdt_e.size(); i++){
    //printCDTskel(cdt_e[i]);
    TestTri.Clear();
    TestTri.Create(cdt_e[i]);
    tot += TestTri.GenerateDimerConfigs(dimer_size, &dimer_count);
  }
  printmap(&dimer_count);
  //cout << tot << endl;
  return 0;
}