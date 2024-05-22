// Andriy_Zatserklyaniy@fnal.gov

#if !defined(__STDLIB_H)
#include <stdlib.h>
#endif

#if !defined(__STDIO_H)
#include <stdio.h>
#endif

#if !defined(__STRING_H)
#include <string.h>
#endif

#if !defined(__IOSTREAM_H)
#include <iostream.h>
#endif

#if !defined(__MATH_H)
#include <math.h>
#endif

#include "wait.h"

// data structure
#include "queue.h"
#include "sidelist.h"
#include "elements.h"

#include "tile.h"

// define static variables for class Tile
int    Tile::count = 0;

int    Tile::OldDrawings = 0;

// semihole diameter. Real value - through Megatile constructor.
double Tile::Dsemihole = 3.;

// init height to start grooving
double Tile::P0height            = 0.5;

double Tile::Dgroove        = 1.35;
double Tile::DepthGroove    = 1.78;
double Tile::GrooveShift    = 3.0;
double Tile::R38            = 38.;

// for narrow down tiles
// Note: depend from current ToleranceShift!

//double Tile::R37            = 37.;
//double Tile::R37            = 36.79; // the best best fit
double Tile::R37            = 36.78; // the best fit
//double Tile::R36            = 36.;
double Tile::R36            = 35.9;  // the best fit
double Tile::R35            = 35.;

double Tile::ShiftP1        = 1.5;
double Tile::ShiftPlast     = 6.5;
double Tile::OddShiftP1     = 52.;
double Tile::OddShiftP11    = 49.;
double Tile::EvenShiftP1    = 62.;
double Tile::EvenShiftP11   = 59.;

char   Tile::tile_file[33];

int    Tile::MinHigh = 9;

//
// Implementation of the classes
//

// class Tile

Tile::Tile(char* tname) {
  char tline[32];
  strcpy(tline, tname);
  char *aname = strtok(tline, " ,;\t");
  name = new char[strlen(aname)+1];
  strcpy(name,aname);

  char *groovesign = strtok(NULL, " ,;\t");
  if (groovesign) {
    sscanf(groovesign, "%ul", &makegroove);
    switch (makegroove) {
      case 0:
      case 1: break;
      default: MESS("Unknown groove sign in line");
               STOPCR(tline);
    }
  }
  else makegroove = 1;

  nameAB = new char[strlen(name)+1];
  strcpy(nameAB,name);
  char *p = strstr(nameAB,"-");
  p++;
  if (*p == '1') *p = 'A';
  else           *p = 'B';

  contour = new SideList;
  groove  = new Groove;
  MakeTile(name);
  count++;
}

Tile::Tile(const Tile& tile)
{
  contour = new SideList;
  groove  = new Groove;
  name = new char[strlen(tile.name)+1];
  strcpy(name,tile.name);
  makegroove = tile.makegroove;
  megatile = tile.megatile;
  megaside = tile.megaside;
  number   = tile.number;
  *contour = *tile.contour;
  *groove  = *tile.groove;
  count++;
}

Tile::~Tile() {
  delete[] name;    name   = NULL;
  delete[] nameAB;  nameAB = NULL;
  delete contour; contour=NULL;
  delete groove; groove=NULL;
  count--;
}

Tile& Tile::operator = (const Tile& tile)
{
  if (this != &tile) {
    name = new char[strlen(tile.name)+1];
    strcpy(name,tile.name);
    nameAB = new char[strlen(tile.nameAB)+1];
    strcpy(nameAB,tile.nameAB);
    megatile = tile.megatile;
    megaside = tile.megaside;
    number   = tile.number;
    *contour = *tile.contour;
    *groove  = *tile.groove;
  }
  return *this;
}

void Tile::Show() const {
  int n = contour->Count();
  cout << "n = " << n << endl;
  contour->Show();
}

double Tile::ToleranceShift()
{
  // Note: ToleranceShift influes on shift in R35, R36 and R38!

  if (OldDrawings) return -0.70;
  else             return -0.75;
}

int Tile::MaxTile(int Nmegatile, int Nmegaside)
{
  // static function

  if (Nmegaside != 1 && Nmegaside != 2) STOPCR("Wrong Nmegaside");

  if (OldDrawings) goto OLD;

  switch (Nmegatile) {
    case  0: switch (Nmegaside) {
               case 1: return 12;
               case 2: return 12;
    }
    case  1: switch (Nmegaside) {
               case 1: return 16;
               case 2: return 16;
    }
    case  2: switch (Nmegaside) {
               case 1: return 16;
               case 2: return 16;
    }
    case  3: switch (Nmegaside) {
               case 1: return 16;
               case 2: return 18;
    }
    case  4: switch (Nmegaside) {
               case 1: return 18;
               case 2: return 18;
    }
    case  5: switch (Nmegaside) {
               case 1: return 18;
               case 2: return 18;
    }
    case  6: switch (Nmegaside) {
               case 1: return 18;
               case 2: return 18;
    }
    case  7: switch (Nmegaside) {
               case 1: return 18;
               case 2: return 18;
    }
    case  8: switch (Nmegaside) {
               case 1: return 18;
               case 2: return 18;
    }
    case  9: switch (Nmegaside) {
               case 1: return 16;
               case 2: return 16;
    }
    case 10: switch (Nmegaside) {
               case 1: return 16;
               case 2: return 16;
    }
    case 11: switch (Nmegaside) {
               case 1: return 16;
               case 2: return 16;
    }
    case 12: switch (Nmegaside) {
               case 1: return 16;
               case 2: return 16;
    }
    case 13: switch (Nmegaside) {
               case 1: return 16;
               case 2: return 16;
    }
    case 14: switch (Nmegaside) {
               case 1: return 14;
               case 2: return 14;
    }
    case 15: switch (Nmegaside) {
               case 1: return 14;
               case 2: return 14;
    }
    case 16: switch (Nmegaside) {
               case 1: return 12;
               case 2: return 12;
    }
    case 17: switch (Nmegaside) {
               case 1: return 12;
               case 2: return 12;
    }
    default: cerr << "MaxTile: No data for megatile " << Nmegatile << endl;
  }
  exit(1);

  OLD:
  ;
  switch (Nmegatile) {
    case  0: switch (Nmegaside) {
               case 1: return 16;
               case 2: return 16;
    }
    case  1: switch (Nmegaside) {
               case 1: return 16;
               case 2: return 16;
    }
    case  2: switch (Nmegaside) {
               case 1: return 18;
               case 2: return 18;
    }
    case  3: switch (Nmegaside) {
               case 1: return 18;
               case 2: return 18;
    }
    case  4: switch (Nmegaside) {
               case 1: return 20;
               case 2: return 20;
    }
    case  5: switch (Nmegaside) {
               case 1: return 20;
               case 2: return 20;
    }
    case  6: switch (Nmegaside) {
               case 1: return 20;
               case 2: return 20;
    }
    case  7: switch (Nmegaside) {
               case 1: return 20;
               case 2: return 20;
    }
    case  8: switch (Nmegaside) {
               case 1: return 20;
               case 2: return 20;
    }
    case  9: switch (Nmegaside) {
               case 1: return 20;
               case 2: return 20;
    }
    case 10: switch (Nmegaside) {
               case 1: return 20;
               case 2: return 20;
    }
    case 11: switch (Nmegaside) {
               case 1: return 18;
               case 2: return 18;
    }
    case 12: switch (Nmegaside) {
               case 1: return 18;
               case 2: return 18;
    }
    case 13: switch (Nmegaside) {
               case 1: return 18;
               case 2: return 18;
    }
    case 14: switch (Nmegaside) {
               case 1: return 18;
               case 2: return 18;
    }
    case 15: switch (Nmegaside) {
               case 1: return 16;
               case 2: return 16;
    }
    case 16: switch (Nmegaside) {
               case 1: return 16;
               case 2: return 16;
    }
    case 17: switch (Nmegaside) {
               case 1: return 16;
               case 2: return 16;
    }
    default: cerr << "MaxTile: No data for megatile " << Nmegatile << endl;
  }
  exit(1);
  return 0;
}

int Tile::UpperTile()
{
  if (OldDrawings) goto OLD;

  switch (megatile) {
    case  0: switch (megaside) {
               case 1: if (number >= 11) return 1; else return 0;
               case 2: if (number >= 11) return 1; else return 0;
    }
    case  1: switch (megaside) {
               case 1: if (number >= 15) return 1; else return 0;
               case 2: if (number >= 15) return 1; else return 0;
    }
    case  2: switch (megaside) {
               case 1: if (number >= 15) return 1; else return 0;
               case 2: if (number >= 15) return 1; else return 0;
    }
    case  3: switch (megaside) {
               case 1: if (number >= 15) return 1; else return 0;
               case 2: if (number >= 17) return 1; else return 0;
    }
    case  4: switch (megaside) {
               case 1: if (number >= 17) return 1; else return 0;
               case 2: if (number >= 17) return 1; else return 0;
    }
    case  5: switch (megaside) {
               case 1: if (number >= 17) return 1; else return 0;
               case 2: if (number >= 17) return 1; else return 0;
    }
    case  6: switch (megaside) {
               case 1: if (number >= 17) return 1; else return 0;
               case 2: if (number >= 17) return 1; else return 0;
    }
    case  7: switch (megaside) {
               case 1: if (number >= 17) return 1; else return 0;
               case 2: if (number >= 17) return 1; else return 0;
    }
    case  8: switch (megaside) {
               case 1: if (number >= 17) return 1; else return 0;
               case 2: if (number >= 17) return 1; else return 0;
    }
    case  9: switch (megaside) {
               case 1: if (number >= 15) return 1; else return 0;
               case 2: if (number >= 15) return 1; else return 0;
    }
    case 10: switch (megaside) {
               case 1: if (number >= 15) return 1; else return 0;
               case 2: if (number >= 15) return 1; else return 0;
    }
    case 11: switch (megaside) {
               case 1: if (number >= 15) return 1; else return 0;
               case 2: if (number >= 15) return 1; else return 0;
    }
    case 12: switch (megaside) {
               case 1: if (number >= 15) return 1; else return 0;
               case 2: if (number >= 15) return 1; else return 0;
    }
    case 13: switch (megaside) {
               case 1: if (number >= 15) return 1; else return 0;
               case 2: if (number >= 15) return 1; else return 0;
    }
    case 14: switch (megaside) {
               case 1: if (number >= 13) return 1; else return 0;
               case 2: if (number >= 13) return 1; else return 0;
    }
    case 15: switch (megaside) {
               case 1: if (number >= 13) return 1; else return 0;
               case 2: if (number >= 13) return 1; else return 0;
    }
    case 16: switch (megaside) {
               case 1: if (number >= 11) return 1; else return 0;
               case 2: if (number >= 11) return 1; else return 0;
    }
    case 17: switch (megaside) {
               case 1: if (number >= 11) return 1; else return 0;
               case 2: if (number >= 11) return 1; else return 0;
    }
    default: cerr << "UpperTile: No data for megatile " << megatile << endl;
  }
  exit(1);

  OLD:
  ;
  switch (megatile) {
    case  0: switch (megaside) {
               case 1: if (number >= 15) return 1; else return 0;
               case 2: if (number >= 15) return 1; else return 0;
    }
    case  1: switch (megaside) {
               case 1: if (number >= 15) return 1; else return 0;
               case 2: if (number >= 15) return 1; else return 0;
    }
    case  2: switch (megaside) {
               case 1: if (number >= 17) return 1; else return 0;
               case 2: if (number >= 17) return 1; else return 0;
    }
    case  3: switch (megaside) {
               case 1: if (number >= 17) return 1; else return 0;
               case 2: if (number >= 17) return 1; else return 0;
    }
    case  4: switch (megaside) {
               case 1: if (number >= 19) return 1; else return 0;
               case 2: if (number >= 19) return 1; else return 0;
    }
    case  5: switch (megaside) {
               case 1: if (number >= 19) return 1; else return 0;
               case 2: if (number >= 19) return 1; else return 0;
    }
    case  6: switch (megaside) {
               case 1: if (number >= 19) return 1; else return 0;
               case 2: if (number >= 19) return 1; else return 0;
    }
    case  7: switch (megaside) {
               case 1: if (number >= 19) return 1; else return 0;
               case 2: if (number >= 19) return 1; else return 0;
    }
    case  8: switch (megaside) {
               case 1: if (number >= 19) return 1; else return 0;
               case 2: if (number >= 19) return 1; else return 0;
    }
    case  9: switch (megaside) {
               case 1: if (number >= 19) return 1; else return 0;
               case 2: if (number >= 19) return 1; else return 0;
    }
    case 10: switch (megaside) {
               case 1: if (number >= 19) return 1; else return 0;
               case 2: if (number >= 19) return 1; else return 0;
    }
    case 11: switch (megaside) {
               case 1: if (number >= 17) return 1; else return 0;
               case 2: if (number >= 17) return 1; else return 0;
    }
    case 12: switch (megaside) {
               case 1: if (number >= 17) return 1; else return 0;
               case 2: if (number >= 17) return 1; else return 0;
    }
    case 13: switch (megaside) {
               case 1: if (number >= 17) return 1; else return 0;
               case 2: if (number >= 17) return 1; else return 0;
    }
    case 14: switch (megaside) {
               case 1: if (number >= 17) return 1; else return 0;
               case 2: if (number >= 17) return 1; else return 0;
    }
    case 15: switch (megaside) {
               case 1: if (number >= 15) return 1; else return 0;
               case 2: if (number >= 15) return 1; else return 0;
    }
    case 16: switch (megaside) {
               case 1: if (number >= 15) return 1; else return 0;
               case 2: if (number >= 15) return 1; else return 0;
    }
    case 17: switch (megaside) {
               case 1: if (number >= 15) return 1; else return 0;
               case 2: if (number >= 15) return 1; else return 0;
    }
    default: cerr << "UpperTile: No data for megatile " << megatile << endl;
  }
  exit(1);
  return 0;
}

double Tile::R38narrow()
{
  if (OldDrawings) return R38;

  if (number > 1) return R38;
  else switch (megatile) {
    case 0: return R35;
    case 1: return R36;
    case 2: return R37;
    default: return R38;
  }
}

void Tile::MakeInner(SideList *inner)
{
  const double eps=1e-7;
  int i, n=contour->Count();

  for (i=0; i<n; i++) if ((*contour)[i]->Length() < eps) {
    cout << " MakeInner: tile " << name; STOP(" has a zero side");
  }

  Point (*Inscribed)(const Point&, const Point&, const Point&, double);
  switch (megaside) {
    // Note: inverted order with respect to TileSide::MakeSides
    case 1: // contour(groove) goes clockwise
            Inscribed = &InscribedInside; break;
    case 2: // contour(groove) goes counterclockwise
            Inscribed = &InscribedOutside;  break;
    default: cerr << "MakeSides: Wrong megaside" << endl; exit(1);
  }

  // calculate tail's end point p2
  Side *side = contour->Tail();
  Point *prev = &side->p1;
  Point *p    = &side->p2;
  Point *next = &contour->Head()->p2;
  Point p2 = (*Inscribed)(*prev,*p,*next,GrooveShift);

  for (i=1; i<n; i++) {
    // p1
    Point p1 = p2;
    // p2
    prev = p;
    p    = next;
    side = (*contour)[i];
    next = &side->p2;
    // cutter center
    p2 = (*Inscribed)(*prev,*p,*next,GrooveShift);

    Side *cut = new Side(p1,p2, todo);
    // don't need: cut->CopySHoles((*contour)[i-1]);
    *inner << cut;
  }

  // make tail cut using already calculated points
  Side *cut = new Side(inner->Tail()->p2, inner->Head()->p1, todo);
  // don't need: cut->CopySHoles((*contour)[i-1]);
  *inner << cut;

  //cout << "Tile " << name << endl;
  //cout << "In tile frame: contour:" << endl; contour->Show();
  //cout << "                 inner:" << endl; inner->Show(); WAIT();
}

void Tile::ClockwiseGroove(int ns_start)
{
  SideList *inner = new SideList;
  MakeInner(inner);

  double shiftP1;
  double shiftPlast;
  // the number of turns
  int nturns;

  // cross point of start inner side with previous contour side
  Point P0 = CrossPoint(*(*inner)[ns_start], *(*contour)[ns_start-1]);
  Side start_ray = Side(P0, (*inner)[ns_start]->p2);

  Point *last = new Point;

  if (!UpperTile()) {
    shiftP1    = ShiftP1 + Dgroove/2.;
    shiftPlast = ShiftPlast;
    nturns = 3;

    // end point P[8]
    int ns_end = ns_start-1; // for clockwise groove
    // lengh of the last contour side
    double len = (*contour)[ns_end]->Length();
    Point p = (*contour)[ns_end]->AtPar(len-shiftPlast);
    *last = (*inner)[ns_end]->ClosestPoint(p);
  }
  else {
    nturns = 4;
    if (number % 2) {
      // Odd
      shiftP1    =  OddShiftP1;
      shiftPlast =  OddShiftP11 - Dgroove/2.;
    }
    else {
      // Even
      shiftP1    = EvenShiftP1;
      shiftPlast = EvenShiftP11 - Dgroove/2.;
    }
    // end point P[10]
    *last = start_ray.AtPar(shiftPlast);
  }

  // compute and store the first point P1
  *groove << new Point(start_ray.AtPar(shiftP1));

  // second point
  double L0 = sqrt(SQR(R38) - SQR(R38-DepthGroove));
  Point *second = new Point(start_ray.AtPar(shiftP1+L0));
  *groove << second;

  // initial point
  double Lh = sqrt(SQR(R38) - SQR(R38-(DepthGroove+P0height)));
  // store the initial point as a hole assigning P0height to D
  Hole *init = new Hole(start_ray.AtPar(Distance(P0,*second)-Lh),P0height);
  groove->InsertBefore(0,init);

  // turns

  Point (*Inscribed)(const Point&, const Point&, const Point&, double);
  switch (megaside) {
    // Note: inverted order with respect to TileSide::MakeSides
    case 1: // contour(groove) goes clockwise
            Inscribed = &InscribedInside; break;
    case 2: // contour(groove) goes counterclockwise
            Inscribed = &InscribedOutside;  break;
    default: cerr << "MakeSides: Wrong megaside" << endl; exit(1);
  }

  int ns=ns_start;
  // the number of turn
  int iR = 1;
  for (int i=0; i<nturns; i++)
  {
    Point prev = (*inner)[ns]->p1;
    Point p    = (*inner)[ns]->p2;
    Point next = (*inner)[ns+1]->p2;
    // The center of the turn
    double Rturn;
    if (iR <= 2) Rturn = R38narrow();
    else         Rturn = R38;

    Hole R((*Inscribed)(prev,p,next,Rturn), Rturn);
    *groove << new Point((*inner)[ns]->  ClosestPoint(R));
    *groove << new Hole(R);
    *groove << new Point((*inner)[ns+1]->ClosestPoint(R));
    ns++;
    iR++;
  }

  // store the last point
  *groove << last;

  /*
  cout << "Tile::ClockwiseGroove: tile " << name << endl;
  cout << "contour" << endl;
  contour->Show();
  cout << "inner" << endl;
  inner->Show();
  groove->Show();

  //cout << "rotate to 10" << endl;
  //groove->Rotate(Point(0.,0.),10.);
  //groove->Show();
  WAIT();
  */

  delete inner; inner=NULL;
}

void Tile::CounterClockwiseGroove(int ns_start)
{
  SideList *inner = new SideList;
  MakeInner(inner);

  double shiftP1;
  double shiftPlast;
  // the number of turns
  int nturns;

  // cross point of start inner side with next contour side
  Point P0 = CrossPoint(*(*inner)[ns_start], *(*contour)[ns_start+1]);
  Side start_ray = Side(P0, (*inner)[ns_start]->p1);

  Point *last = new Point;

  if (!UpperTile()) {
    shiftP1    = ShiftP1 + Dgroove/2.;
    shiftPlast = ShiftPlast;
    nturns = 3;

    // end point P[8]
    int ns_end = ns_start+1; // for counterclockwise groove
    Point p = (*contour)[ns_end]->AtPar(shiftPlast);
    *last = (*inner)[ns_end]->ClosestPoint(p);
  }
  else {
    nturns = 4;
    if (number % 2) {
      // Odd
      shiftP1    =  OddShiftP1;
      shiftPlast =  OddShiftP11 - Dgroove/2.;
    }
    else {
      // Even
      shiftP1    = EvenShiftP1;
      shiftPlast = EvenShiftP11 - Dgroove/2.;
    }
    // end point P[10]
    *last = start_ray.AtPar(shiftPlast);
  }

  // compute and store the start point
  *groove << new Point(start_ray.AtPar(shiftP1));

  // second point
  double L0 = sqrt(SQR(R38) - SQR(R38-DepthGroove));
  Point *second = new Point(start_ray.AtPar(shiftP1+L0));
  *groove << second;

  // initial point
  double Lh = sqrt(SQR(R38) - SQR(R38-(DepthGroove+P0height)));
  // store the initial point as a hole assigning P0height to D
  Hole *init = new Hole(start_ray.AtPar(Distance(P0,*second)-Lh),P0height);
  groove->InsertBefore(0,init);

  // turns

  Point (*Inscribed)(const Point&, const Point&, const Point&, double);
  switch (megaside) {
    // Note: order like to TileSide::MakeSides
    case 1: // contour(groove) goes counterclockwise
            Inscribed = &InscribedOutside; break;
    case 2: // contour(groove) goes clockwise
            Inscribed = &InscribedInside;  break;
    default: cerr << "MakeSides: Wrong megaside" << endl; exit(1);
  }

  int ns=ns_start;
  // the number of turn
  int iR = 1;
  for (int i=0; i<nturns; i++)
  {
    Point prev = (*inner)[ns]->p2;    // p1 --> p2
    Point p    = (*inner)[ns]->p1;    // p2 --> p1
    Point next = (*inner)[ns-1]->p1;  // ns-1 & p2 --> p1
    // The center of the turn
    double Rturn;
    if (iR <= 2) Rturn = R38narrow();
    else         Rturn = R38;

    Hole R((*Inscribed)(prev,p,next,Rturn), Rturn);
    *groove << new Point((*inner)[ns]->  ClosestPoint(R));
    *groove << new Hole(R);
    *groove << new Point((*inner)[ns-1]->ClosestPoint(R));
    ns--; // --
    iR++;
  }

  // store the last point
  *groove << last;

  delete inner; inner=NULL;
}

// class TileSide

TileSide::TileSide(char *tname, double dcut): Tile(tname) {
  sides = new SideList;
  Dcut = dcut;
  MakeSides();

  // coordinate origin now in point contour->Head()->p1
  // shift coordinate origin to point sides->Head()->p1

  Point p = sides->Head()->p1;

  //cout << "TileSide::TileSide() before Move to "; p.Show(); Show();

  Move(-p);

  //cout << "TileSide::TileSide() after Move "; Show(); WAIT();
}

TileSide& TileSide::operator = (const TileSide& tile)
{
  if (this != &tile) {
    //name = new char[strlen(tile.name)+1];
    //strcpy(name,tile.name);
    //nameAB = new char[strlen(tile.nameAB)+1];
    //strcpy(nameAB,tile.nameAB);
    //megatile = tile.megatile;
    //megaside = tile.megaside;
    //number   = tile.number;
    //*contour = *tile.contour;
    //*groove  = *tile.groove;
    Tile(tile.name);
    *sides   = *tile.sides;
  }
  return *this;
}

Point TileSide::ClosestEnd(const Point& topoint)
{
  double r2;
  Point p = sides->Head()->ClosestEnd(topoint,r2);
  for (int i=1; i<sides->Count(); i++) {
    double r2curr;
    Point   pcurr = (*sides)[i]->ClosestEnd(topoint,r2curr);
    if (r2curr < r2) {
      p  =  pcurr;
      r2 = r2curr;
    }
  }
  return p;
}

Point TileSide::ClosestEndX(const Point& topoint)
{
  Point pmin = sides->Head()->ClosestEnd(topoint);
  double dmin = fabs(pmin.X - topoint.X);
  for (int ns=1; ns<sides->Count(); ns++) {
    Point p = (*sides)[ns]->ClosestEnd(topoint);
    double d = fabs(p.X-topoint.X);
    if (d < dmin) {
      pmin = p;
      dmin = d;
    }
  }
  return pmin;
}

Point TileSide::ClosestEndY(const Point& topoint)
{
  Point pmin = sides->Head()->ClosestEnd(topoint);
  double dmin = fabs(pmin.Y - topoint.Y);
  for (int ns=1; ns<sides->Count(); ns++) {
    Point p = (*sides)[ns]->ClosestEnd(topoint);
    double d = fabs(p.Y-topoint.Y);
    if (d < dmin) {
      pmin = p;
      dmin = d;
    }
  }
  return pmin;
}

Point TileSide::ClosestEndXY(const Point& topoint)
{
  Point p = ClosestEnd(topoint);
  // this point can have large X-distance,
  // so point with the same Y and required X is
  Point X(p.Y, topoint.X);
  // and to be ensure
  p = ClosestEnd(X);
  return p;
}

Point TileSide::ClosestEndYX(const Point& topoint)
{
  Point p = ClosestEnd(topoint);
  // this point can have large Y-distance,
  // so point with the same X and required Y is
  Point Y(p.X, topoint.Y);
  // and to be ensure
  p = ClosestEnd(Y);
  return p;
}

Point TileSide::ContourClosestEndX(const Point& topoint)
{
  Point pmin = contour->Head()->ClosestEnd(topoint);
  double dmin = fabs(pmin.X - topoint.X);
  for (int ns=1; ns<contour->Count(); ns++) {
    Point p = (*contour)[ns]->ClosestEnd(topoint);
    double d = fabs(p.X-topoint.X);
    if (d < dmin) {
      pmin = p;
      dmin = d;
    }
  }
  return pmin;
}

Side* TileSide::ClosestSide(const Segment& tosegment)
{
  const double eps = 1.e-7;
  double r2p1,r2p2;
  Side *side = sides->Head();
  Point p1 = tosegment.ClosestPoint(side->p1);
  Point p2 = tosegment.ClosestPoint(side->p2);
  r2p1 = R2(p1,side->p1);
  r2p2 = R2(p2,side->p2);
  Side *pside = side;
  for (int i=1; i<sides->Count(); i++) {
    double r2p1curr,r2p2curr;
    Side *side = (*sides)[i];
    Point p1 = tosegment.ClosestPoint(side->p1);
    Point p2 = tosegment.ClosestPoint(side->p2);
    r2p1curr = R2(p1,side->p1);
    r2p2curr = R2(p2,side->p2);
    //if (r2p1curr <= r2p1 && r2p2curr <= r2p2) pside = side;
    if (r2p1-r2p1curr > -eps && r2p2-r2p2curr > -eps) {
      pside = side;
      r2p1 = r2p1curr;
      r2p2 = r2p2curr;
    }
  }
  return pside;
}

Side* TileSide::ClosestSide(const Segment& tosegment,
                            double& r2p1, double& r2p2) {
  const double eps = 1.e-7;
  Side *side = sides->Head();
  Point p1 = tosegment.ClosestPoint(side->p1);
  Point p2 = tosegment.ClosestPoint(side->p2);
  r2p1 = R2(p1,side->p1);
  r2p2 = R2(p2,side->p2);
  Side *pside = side;
  for (int i=1; i<sides->Count(); i++) {
    double r2p1curr,r2p2curr;
    Side *side = (*sides)[i];
    Point p1 = tosegment.ClosestPoint(side->p1);
    Point p2 = tosegment.ClosestPoint(side->p2);
    r2p1curr = R2(p1,side->p1);
    r2p2curr = R2(p2,side->p2);
    if (r2p1-r2p1curr > -eps && r2p2-r2p2curr > -eps) {
      pside = side;
      r2p1 = r2p1curr;
      r2p2 = r2p2curr;
    }
  }
  return pside;
}

/*
int TileSide::Include(const Point& p)
{
  const double eps = 1.e-7;
  const double pi2 = 2.*M_PI;
  double sum = 0.;
  //for (int i=0; i<sides->Count(); i++) sum += VertexAngleExclude(p,*(*sides)[i]);
  for (int i=0; i<sides->Count(); i++) {
    Side *side = (*sides)[i];
    double angle = VertexAngleExclude(p,*side);
    double adegr = angle*180./M_PI;
    sum += angle;
  }
  if (sum < pi2+eps) return 0;
  else               return 1;
}
*/
int TileSide::Include(const Point& p)
{
  const double eps = 1.e-7;
  // find the closest side to p
  Point p1   = sides->Head()->p1;
  Point p2   = sides->Head()->p2;
  Point pmin = sides->Head()->ClosestPoint(p);
  double dmin = Distance(pmin, p);
  for (int i=1; i<sides->Count(); i++) {
    Side *side = (*sides)[i];
    Point  pcur = side->ClosestPoint(p);
    double dcur = Distance(pcur, p);
    if (dcur < dmin) {
      p1 = side->p1;
      p2 = side->p2;
      pmin = pcur;
      dmin = dcur;
    }
  }

  if (dmin < eps) return 0; // p situated on the side

  // distances to ends
  if (Distance(p1,pmin) < eps || Distance(p2,pmin) < eps) return 0;

  // here pmin between the ends p1 and p2

  LinePar vs(pmin, p2);
  LinePar vp(pmin, p);

  // calc z-component of vector product (should be +1 or -1)
  double cosz = vs.cosx*vp.cosy - vs.cosy*vp.cosx;

  if (cosz*cosZ() > 0.) return 1;
  else                  return 0;
}

double TileSide::cosZ()
{
  const double eps = 1.e-7;

  Side *side1 = (*sides)[0];
  double theta1 = side1->Theta();
  Side *side2;
  for (int i=1; i<sides->Count(); i++) {
    side2 = (*sides)[i];
    if (fabs(side2->Theta() - theta1) > eps) break;
  }
  LinePar s1(*side1);
  LinePar s2(*side2);
  // calc z-component of vector product
  return s1.cosx*s2.cosy - s1.cosy*s2.cosx;
}

int TileSide::LetStrip(Side* alien)
{
  // return 1 if some tile side stripped by alien, 0 otherwise

  Side *overlap = new Side;
  overlap->status = done;

  Side *side;
  int nsides = sides->Count();
  for (int nside=0; nside<nsides; nside++)
  {
    side = (*sides)[nside];
    int overlapped = side->Overlap(alien, overlap);
    if (overlapped) goto REPLACE;
  }
  delete overlap;
  return 0;

  REPLACE:
  // let alien strip sholes from 'side'
  alien->StripSHoles(side);

  // Replace the 'side' by 'overlap' marked 'done'
  sides->ReplaceItem(side,overlap);

  LinePar This(side);
  double tmax = This.tmax;

  // parameters of ends in my frame
  double t1 = This.t(&overlap->p1);
  double t2 = This.t(&overlap->p2);

  if (t1 > t2) {
    double buf = t1;
    t1 = t2;
    t2 = buf;
  }

  double eps=1.E-6;

  if (t1 < eps) t1=0.;
  if (t2 > tmax-eps) t2=tmax;

  // if (t1 == 0 && t2 == tmax) { // side is entirely overlapped }

  if (t1 > eps) {
    // create new side left
    Side* left = new Side(side->p1, overlap->p1);
    left->StripSHoles(side);
    // add left side in the sides
    sides->InsertBefore(overlap,left);
  }

  if (fabs(tmax-t2) > eps) {
    // create new side right
    Side* right = new Side(overlap->p2, side->p2);
    right->StripSHoles(side);
    // add right side in the sides
    sides->InsertAfter(overlap,right);
  }

  // the sholes from 'side' are stripped
  delete side;

  //cout << "TileSide::LetStrip(Side* alien) before int. check"
  //     << endl;
  //sides->Show();
  return 1;
}

void TileSide::MakeSides()
{
  const double eps=1e-7;
  int i, n=contour->Count();

  for (i=0; i<n; i++) if ((*contour)[i]->Length() < eps) {
    cout << " MakeSides: tile " << name; STOP(" has a zero side");
  }

  Point (*Inscribed)(const Point&, const Point&, const Point&, double);
  switch (megaside) {
    case 1: // contour(groove) goes clockwise
            Inscribed = &InscribedOutside; break;
    case 2: // contour(groove) goes counterclockwise
            Inscribed = &InscribedInside;  break;
    default: cerr << "MakeSides: Wrong megaside" << endl; exit(1);
  }

  // calculate tail's end point p2
  Side *side = contour->Tail();
  Point *prev = &side->p1;
  Point *p    = &side->p2;
  Point *next = &contour->Head()->p2;
  //Point p2 = InscribedOutside(*prev,*p,*next,Dcut/2.);
  Point p2 = (*Inscribed)(*prev,*p,*next,Dcut/2.);

  for (i=1; i<n; i++) {
    // p1
    Point p1 = p2;
    // p2
    prev = p;
    p    = next;
    side = (*contour)[i];
    next = &side->p2;
    // cutter center
    //p2 = InscribedOutside(*prev,*p,*next,Dcut/2.);
    p2 = (*Inscribed)(*prev,*p,*next,Dcut/2.);

    Side *cut = new Side(p1,p2, todo);
    cut->CopySHoles((*contour)[i-1]);
    *sides << cut;
  }

  // make tail cut using already calculated points
  Side *cut = new Side(sides->Tail()->p2, sides->Head()->p1, todo);
  cut->CopySHoles((*contour)[i-1]);
  *sides << cut;

  //cout << "Tile " << name << endl;
  //cout << "In tile frame: contour:" << endl; contour->Show();
  //cout << "                 sides:" << endl; sides->Show(); WAIT();
}

// struct TileData

TileData::TileData(char *tname, double Dcut,
                   double dx,double dy, double dtheta, int inv)
{
  tile = new TileSide(tname, Dcut);
  status = todo;
  defined_frame = 1;
  Move(dx,dy,dtheta,inv);
}

TileData::TileData(char *tname, double Dcut, double dtheta, int inv)
{
  tile = new TileSide(tname, Dcut);
  status = todo;
  defined_frame = 0;
  Move(0.,0.,dtheta,inv);
  //cout << "After TileData constructor: "; Show(); WAIT();
}

TileData::TileData(char *tname, double Dcut, double dtheta)
{
  tile = new TileSide(tname, Dcut);
  status = todo;
  defined_frame = 0;
  Move(0.,0.,dtheta,0);
}

TileData::TileData(char *tname, double Dcut, int inv)
{
  tile = new TileSide(tname, Dcut);
  status = todo;
  defined_frame = 0;
  Move(0.,0.,0.,inv);
}

TileData& TileData::operator = (const TileData& td)
{
  if (this != &td) {
    *tile  = *td.tile;
    status = td.status;
    frame  = td.frame;
    defined_frame = td.defined_frame;
  }
  return *this;
}

void TileData::Move(double dx, double dy, double dtheta, int inv)
{
  frame.X0 += dx;
  frame.Y0 += dy;
  frame.theta += dtheta;
  tile->sides->Move(dx,dy,dtheta,inv);
  tile->contour->Move(dx,dy,dtheta,inv);
  tile->groove->Move(dx,dy,dtheta,inv);
  defined_frame = 1;
}

void TileData::Move(double dx, double dy, double dtheta)
{
  tile->sides->Move(dx,dy,dtheta,0);
  tile->contour->Move(dx,dy,dtheta,0);
  tile->groove->Move(dx,dy,dtheta,0);
  defined_frame = 1;
}

void TileData::Move(double dx, double dy)
{
  frame.X0 += dx;
  frame.Y0 += dy;
  tile->sides->Move(dx,dy,0.,0);
  tile->contour->Move(dx,dy,0.,0);
  tile->groove->Move(dx,dy,0.,0);
  defined_frame = 1;
}

void TileData::Move(const Point& r)
{
  Point dr = Point(r.X-frame.X0, r.Y-frame.Y0);
  frame.X0 = r.X;
  frame.Y0 = r.Y;
  tile->sides->Move(dr);
  tile->contour->Move(dr);
  tile->groove->Move(dr);
  defined_frame = 1;
}
