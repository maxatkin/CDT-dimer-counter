#include <vector>
#include <map>
#include <string>
#include <stdio.h>
#define UP true
#define DOWN false
#define DIAMOND 0
#define UPDOWN 1
#define DOWNUP 2
#define WEDGE 3

using namespace std;
typedef vector<vector<unsigned int> > cdt_skeleton;

class Triangle{
  public:
    Triangle();
    bool orientation;
    bool contains_dimer;
    Triangle *dimer_end;
    Triangle *on_left_side;
    Triangle *on_right_side;
    Triangle *on_hrztl_side;
  private:
};

class Dimer{
  public:
    Triangle *End1;
    Triangle *End2;
    Triangle *OtherEnd(Triangle *ThisEnd);
};

class Triangulation{
  public:
    Triangulation();
    unsigned int Create(cdt_skeleton cdt_skel);
    void StartNewLayer();
    void AddTriangle(bool orientation);
    void Print(bool verbose);
    void FindDimerBoxes(vector<Triangle *> *EmptyBoxes, unsigned int type);
    void SetDimers(vector<Triangle *> *BoxSelection, unsigned int type);
    void ClearDimers(vector<Triangle *> *EmptyBoxes, unsigned int type);
    unsigned int CountDimersInWedgeBoxes();
    unsigned int GenerateDimerConfigs(unsigned int no_of_dimers, vector<unsigned int> *dimer_count);
    void Clear();
    unsigned int CountDimers();

    vector<vector<Triangle *> *> Simplex;
    vector<vector<Dimer *> > DimerConfigs;
    vector<Triangle *> EmptyDiamondBoxes;
    vector<Triangle *> UpDownBoxes;
    vector<Triangle *> DownUpBoxes;

    vector<Triangle *> *CurrentLayer;
    vector<Triangle *> *LastLayer;
    Triangle *NextDownPointTri;
    Triangle *FindNextDownPointTri();
    unsigned int size;
  private:
};

