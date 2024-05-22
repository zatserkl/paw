// Andriy_Zatserklyaniy@fnal.gov

#if !defined(__STDLIB_H)
#include <stdlib.h>
#endif

#if !defined(__MATH_H)
#include <math.h>
#endif

#if !defined(__IOSTREAM_H)
#include <iostream.h>
#endif

#include "wait.h"

struct Point;
struct SHole;

template <class T> class Queue;
class SideList;

class Tile;
Tile *MakeTile(char *name);

class Tile {
  // static members
//private:
public:
  static int count;
public:
  static int OldDrawings;

  static int get_count() {return count;};

  // semihole diameter, defined by Megatile::Dcut in Megatile constructor
  static double Dsemihole    ;

  static double Dgroove      ;
  static double DepthGroove  ;
  static double R38          ;
  static double R37          ;
  static double R36          ;
  static double R35          ;
  static double GrooveShift  ;

  // init height to start grooving
  static double P0height     ;

  static double ShiftP1      ;
  static double ShiftPlast   ;
  static double OddShiftP1   ;
  static double OddShiftP11  ;
  static double EvenShiftP1  ;
  static double EvenShiftP11 ;

  static int MinHigh;
  static char tile_file[33];

  static int MaxTile(int Nmegatile, int Nmegaside);
  static double ToleranceShift();

  // other members
public:
  char* name;
  char* nameAB;
  int makegroove;
  int megatile, megaside, number;
  SideList *contour;
  Groove *groove;

private:
  int UpperTile();
  double R38narrow();
  void MakeInner(SideList *inner);
  void ClockwiseGroove(int ns_start);
  void CounterClockwiseGroove(int ns_start);
  void MakeTile(char *tname);
  int OddLow(char* line);
  int OddHigh(char* line);
  int EvenLow(char* line);
  int EvenHigh(char* line);
public:
   Tile(char* nameInit);
   Tile(const Tile& tile);
  ~Tile();
  Tile& operator = (const Tile& tile);
  virtual void Show() const;
};

class TileSide: public Tile {
//
// Meaning of side status:
// -----------------------
// todo  - side to do
// done  - side which done but necessary for cutting sequence
// extra - side which is not necessary for cutting but close tile chain
// lift  - side which should be passed with lifted cutter
//
// Notes:
// -----
// cutter can down in every place on side 'done';
// sides which non-lift draw the tile contour
//
private:
  double Dcut;
public:
  SideList* sides;
private:
  void MakeSides();
  void Move(double dx,double dy, double dtheta, int inv) {
    sides->Move(dx,dy,dtheta,inv);
    groove->Move(dx,dy,dtheta,inv);
    contour->Move(dx,dy,dtheta,inv);
  }
  void Move(const Point& r) {
    sides->Move(r);
    groove->Move(r);
    contour->Move(r);
  }
public:
   TileSide(char *tname, double dcut);
  ~TileSide() {delete sides; sides=NULL;}
  TileSide& operator = (const TileSide& tile);
  void Show() const {cout << "Tile " << nameAB << endl; sides->Show();}
  int LetStrip(Side* alien);
  double cosZ();
  Point ClosestEnd(const Point& point);
  Point ClosestEndX(const Point& point);
  Point ClosestEndY(const Point& point);
  Point ClosestEndXY(const Point& point);
  Point ClosestEndYX(const Point& point);
  Point ContourClosestEndX(const Point& point);
  Side* ClosestSide(const Segment& tosegment);
  Side* ClosestSide(const Segment& tosegment, double& r2p1, double& r2p2);
  int Include(const Point& p);
  void AddSide(Side* NewSide) {*sides << NewSide;}
  void RemoveSide(Side *side) {
    if (sides->Count()) sides->Remove(side);
    else WAIT("RemoveSide(Side *side): Trying to remove from empty list!");
  }
  void RemoveSide(int n) {RemoveSide(sides->GetItem(n));}
  void KillSide(Side *side) {
    if (sides->Count()) sides->Kill(side);
    else WAIT("KillSide(Side *side): Trying to remove from empty list!");
  }
  void KillSide(int n) {KillSide(sides->GetItem(n));}
  void SetHead(Side *side) {sides->SetHead(side);}
  void SetHead(int n) {sides->SetHead(n);}
  void Optimize() {sides->Optimize();}
  void Rotate(const Point& O, double theta) {
    contour->Rotate(O,theta);
     groove->Rotate(O,theta);
      sides->Rotate(O,theta);
  }
  void TurnBack(double A, double B) {
    contour->TurnBack(A,B);
     groove->TurnBack(A,B);
      sides->TurnBack(A,B);
  }
};

struct TileData {
protected:
  Frame frame;
public:
  TileSide *tile;
  WorkStatus status;
  int defined_frame;
private:
  void MakeSides();
  public:
   TileData(char *tname, double Dcut,
            double x0,double y0, double thdeg, int up);
   TileData(char *tname, double Dcut, double thdeg, int up);
   TileData(char *tname, double Dcut, double thdeg);
   TileData(char *tname, double Dcut, int up);
  ~TileData() {delete tile; tile=NULL;}
  TileData& operator = (const TileData& td);
  void Move(double dx,double dy, double dtheta, int inv);
  void Move(double dx,double dy, double dtheta);
  void Move(double dx,double dy);
  void Move(const Point& r);
  void Rotate(const Point& O, double theta) {tile->Rotate(O,theta);}
  void TurnBack(double A, double B)  {tile->TurnBack(A,B);}
public:
  void Show() const {tile->Show();}
};
