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

/*
#if !defined(__MATH_H)
#include <math.h>
#endif
*/

#if !defined(__IOSTREAM_H)
#include <iostream.h>
#endif

#include "queue.h"
#include "elements.h"

#include "sidelist.h"

#include "tile.h"
#include "megatile.h"

//class TileList: public Queue<TileData> {
//
//};

void Megatile::MakeName(char* name, int megatile, int megaside, int number)
{
  sprintf(name, "%d-%d-%d", megatile, megaside, number);
}

int Megatile::ExistName(char *name, Queue<TileData> *tds)
{
  for (int i=0; i<tds->Count(); i++)
    if (strcmp(name, (*tds)[i]->tile->name) == 0) return 1;
  return 0;
}

int Megatile::IndexName(char *name, Queue<TileData> *tds)
{
  for (int i=0; i<tds->Count(); i++)
    if (strcmp(name, (*tds)[i]->tile->name) == 0) return i;
  return -1;
}

TileData* Megatile::tdName(char *name, Queue<TileData> *tds)
{
  for (int i=0; i<tds->Count(); i++)
    if (strcmp(name, (*tds)[i]->tile->name) == 0) return (*tds)[i];
  return NULL;
}

/*
int Megatile::SheetIncludes(TileSide *tile)
  BAD ROUTINE
{
  const double eps = 1.e-6;
  const double pi2 = 2.*M_PI;
  Segment s[4];
  s[0]=Segment(SheetX1+xShift,SheetY1+yShift, SheetX1+xShift,SheetY2-yShift);
  s[1]=Segment(SheetX1+xShift,SheetY2-yShift, SheetX2-xShift,SheetY2-yShift);
  s[2]=Segment(SheetX2-xShift,SheetY2-yShift, SheetX2-xShift,SheetY1+yShift);
  s[3]=Segment(SheetX2-xShift,SheetY1+yShift, SheetX1+xShift,SheetY1+yShift);

  for (int nt=0; nt<tile->sides->Count(); nt++) {
    double sum;
    Side *side = (*tile->sides)[nt];
    int i;
    sum = 0.;
    for (i=0; i<4; i++) sum += VertexAngle(side->p1, s[i]);
    if (sum < pi2-eps) return 0;
    sum = 0.;
    for (i=0; i<4; i++) sum += VertexAngle(side->p2, s[i]);
    if (sum < pi2-eps) return 0;
  }
  return 1;
}
*/

int Megatile::SheetIncludes(TileSide *tile)
{
  const double eps = 1.e-6;
                                     //Norm   //MAHO
  double xmin = (SheetX2 > SheetX1)? SheetX1: SheetX2;
  double xmax = (SheetX2 > SheetX1)? SheetX2: SheetX1;
  double ymin = (SheetY2 > SheetY1)? SheetY1: SheetY2;
  double ymax = (SheetY2 > SheetY1)? SheetY2: SheetY1;
  xmin -= eps;
  xmax += eps;
  ymin -= eps;
  ymax += eps;
  for (int nt=0; nt<tile->contour->Count(); nt++) {
    Point p;
    Side *side = (*tile->contour)[nt];
    p = side->p1;
    if (p.X < xmin || p.X > xmax) return 0;
    if (p.Y < ymin || p.Y > ymax) return 0;
    p = side->p2;
    if (p.X < xmin || p.X > xmax) return 0;
    if (p.Y < ymin || p.Y > ymax) return 0;
  }
  return 1;
}

int Megatile::InSheetWorkArea()
{
  for (int i=0; i<tds->Count(); i++) {
    if (!SheetIncludes((*tds)[i]->tile)) return 0;
  }
  return 1;
}

int Megatile::TilesOverlap()
{
  for (int i=0; i<tds->Count(); i++) {
    TileData *td = (*tds)[i];
    for (int ns=0; ns<td->tile->sides->Count(); ns++) {
      Side *side = (*td->tile->sides)[ns];
      Point p1(side->p1);
      Point p2(side->p2);
      for(int j=0; j<tds->Count(); j++) // from 0 to avoid dipped tiles
      {
        TileData *tdcheck = (*tds)[j];
        if (tdcheck == td) continue;
        if (tdcheck->tile->Include(p1)) return 1;
        if (tdcheck->tile->Include(p2)) return 1;
      }
    }
  }
  return 0;
}

int Megatile::Fit()
{
  if (!InSheetWorkArea()) return 0;
  if (TilesOverlap())     return 0;
  return 1;
}

int Ntile(char *tname) {
  int megatile,number;
  sscanf(tname, "%d%*1s%*1s%*1s%d", &megatile,&number);
  return number;
}

void Megatile::Sort()
{
  int mtile0, mside0;

  if (NameList->Count() == 0) {
    cerr << "Empty tiles list!" << endl;
    exit(1);
  }
  else sscanf(NameList->Head(), "%d%*1s%1d%", &mtile0,&mside0);

  for (int i=1; i<NameList->Count(); i++)
  {
    int mtile,mside,num;

    char *tname = (*NameList)[i];
    sscanf(tname, "%d%*1s%1d%*1s%d", &mtile,&mside,&num);

    // check the megaside
    if (mtile != mtile0 || mside != mside0)
    {
      cerr << "Sorry, different megatiles:" << endl;
      STOPCR("Program processes tiles from the same megatile!");
    }

    // check the num <= max
    if (num > Tile::MaxTile(mtile,mside))
    {
      char nameAB[32];
      strcpy(nameAB, tname);
      // change 1 --> A, 2 --> B
      char *AB = strstr(nameAB,"-") + 1;
      switch (*AB) {
        case 'A':
        case '1': *AB = 'A'; break;
        case 'B':
        case '2': *AB = 'B'; break;
        default: cerr << "Wrong tile name "; STOPCR(nameAB);
      }
      cerr << "The number of tile " << nameAB << endl;
      cerr << "exeeds the max tile number in the megatile " << mtile << endl;
      STOPCR();
    }
  }

  int swapped = 1;
  while (swapped) {
    swapped = 0;
    for (int i=0; i<NameList->Count()-1; i++)
    {
      int n = i;
      char *nname = (*NameList)[n];
      for (int n1=n+1; n1<NameList->Count(); n1++)
      {
        if (Ntile(nname) > Ntile((*NameList)[n1]))
        {
          NameList->Swap(n,n1);
          n = n1;
          swapped = 1;
        }
      }
    }
  }
}

void Megatile::TurnBack()
{
  double A = DownRightPoint.X;
  double B =    UpLeftPoint.X;
  for (int i=0; i<tds->Count(); i++) (*tds)[i]->TurnBack(A,B);
}

void Megatile::DefineLimitPoints()
{
  // redefine the points UpLeftPoint and DownRightPoint
  Point p1,p2;

  p1 = ClosestX(tds, DownRightPoint);
  p2 = ClosestY(tds, DownRightPoint);
  // redefine the DownRightPoint
  DownRightPoint = Point(p1.X, p2.Y);

  p1 = ClosestX(tds, UpLeftPoint);
  p2 = ClosestY(tds, UpLeftPoint);
  // redefine the UpLeftPoint
  UpLeftPoint = Point(p1.X, p2.Y);
}

void Megatile::Arrange()
{
  if (NameList->Count() == 0) {
    cerr << "Empty tiles list!" << endl;
    exit(1);
  }

  /*
  cout << "Tiles to arrange:" << endl;
  for (int j=0; j<NameList->Count(); j++) cout << (*NameList)[j] << endl;
  cout << "----------------" << endl;
  */

  cout << endl << "Standard arranging. Please, wait." << endl;

  Sort();

  /*
  cout << "After tiles sorting:" << endl;
  for (int k=0; k<NameList->Count(); k++) cout << (*NameList)[k] << endl;
  cout << "----------------" << endl;
  WAIT();
  */

  int megatile,megaside,number;
  char *tname = (*NameList)[0];
  sscanf(tname, "%d%*1s%1d%*1s%d", &megatile,&megaside,&number);

  if (megaside % 2) ArrangeMegaSide1();
  else              ArrangeMegaSide2();

  if (SqueezeParts) Squeeze();

  // redefine the points UpLeftPoint and DownRightPoint
  DefineLimitPoints();

  if (BackSide)
  {
    TurnBack();
    // swap RightPart <--> LeftPart
    Queue<TileData> *BufPart = new Queue<TileData>;

    while (LeftPart->Count())  *BufPart   << LeftPart->Remove(0);
    while (RightPart->Count()) *LeftPart  << RightPart->Remove(0);
    while (BufPart->Count())   *RightPart << BufPart->Remove(0);
    delete BufPart;
  }

  if (!Fit()) {
    nopath = 1;
    cout << endl << "This tile set doesn't fit in work area!" << endl;
    cout         << "---------------------------------------" << endl;
    WAIT("Code will not be created.\n");
  }
}

double Megatile::dxPart(const Queue<TileData> *list, const Point& point)
{
  const double eps = 1.e-7;
  Segment line(DownRightPoint.X,point.Y, UpLeftPoint.X, point.Y);
  //double sign = (SheetX2 > SheetX1)? +1.: -1.;
  double sign = (point.X > SheetX1)? +1.: -1.;
  //double dxmin = SheetX2-SheetX1;
  double dxmin = point.X-SheetX1;
  // for check of horizontal side
  double dY = 5.; // mm

  for (int nt=0; nt<list->Count(); nt++) {
    TileSide *tile = (*list)[nt]->tile;
    Point p1 = tile->ClosestEndY(UpLeftPoint);
    Point p2 = tile->ClosestEndY(DownRightPoint);
    double ymin = (p1.Y < p2.Y)? p1.Y: p2.Y;
    double ymax = (p1.Y > p2.Y)? p1.Y: p2.Y;
    ymin -= eps;
    ymax += eps;
    if (point.Y < ymin || point.Y > ymax) continue;
    for (int ns=0; ns<tile->sides->Count(); ns++) {
      Segment *side = (*tile->sides)[ns];
      if (fabs(side->p1.Y - side->p2.Y) < dY) continue; // horizontal side
      Point p = CrossPoint(line,*side);
      double dx = point.X-p.X;
      if (sign*dx < sign*dxmin) dxmin=dx;
    }
  }
  return dxmin;
}

double Megatile::dxPart(const Queue<TileData> *list, const Point& point, double x0)
{
  const double eps = 1.e-7;
  Segment line(DownRightPoint.X,point.Y, UpLeftPoint.X, point.Y);
  //double sign = (SheetX2 > SheetX1)? +1.: -1.;
  double sign = (point.X > x0)? +1.: -1.;
  //double dxmin = SheetX2-SheetX1;
  double dxmin = point.X-x0;
  // for check of horizontal side
  double dY = 5.; // mm

  for (int nt=0; nt<list->Count(); nt++) {
    TileSide *tile = (*list)[nt]->tile;
    Point p1 = tile->ClosestEndY(UpLeftPoint);
    Point p2 = tile->ClosestEndY(DownRightPoint);
    double ymin = (p1.Y < p2.Y)? p1.Y: p2.Y;
    double ymax = (p1.Y > p2.Y)? p1.Y: p2.Y;
    ymin -= eps;
    ymax += eps;
    if (point.Y < ymin || point.Y > ymax) continue;
    for (int ns=0; ns<tile->sides->Count(); ns++) {
      Segment *side = (*tile->sides)[ns];
      if (fabs(side->p1.Y - side->p2.Y) < dY) continue; // horizontal side
      Point p = CrossPoint(line,*side);
      double dx = point.X-p.X;
      if (sign*dx < sign*dxmin) dxmin=dx;
    }
  }
  return dxmin;
}

double Megatile::dxRightPart(const Point& point)
{
  const double eps = 1.e-7;
  Segment line(DownRightPoint.X,point.Y, UpLeftPoint.X, point.Y);
  double sign = (SheetX2 > SheetX1)? +1.: -1.;
  //double dxmin = point.X-SheetX1;
  double dxmin = SheetX2-SheetX1;
  // for check of horizontal side
  double dY = 1.; // mm

  for (int nt=0; nt<RightPart->Count(); nt++) {
    TileSide *tile = (*RightPart)[nt]->tile;
    Point p1 = tile->ClosestEndY(UpLeftPoint);
    Point p2 = tile->ClosestEndY(DownRightPoint);
    double ymin = (p1.Y < p2.Y)? p1.Y: p2.Y;
    double ymax = (p1.Y > p2.Y)? p1.Y: p2.Y;
    ymin -= eps;
    ymax += eps;
    if (point.Y < ymin || point.Y > ymax) continue;
    for (int ns=0; ns<tile->sides->Count(); ns++) {
      Segment *side = (*tile->sides)[ns];
      if (fabs(side->p1.Y - side->p2.Y) < dY) continue; // horizontal side
      Point p = CrossPoint(line,*side);
      double dx = point.X-p.X;
      if (sign*dx < sign*dxmin) dxmin=dx;
    }
  }
  return dxmin;
}

double Megatile::dxLeftPart(const Point& point)
{
  const double eps = 1.e-7;
  Segment line(DownRightPoint.X,point.Y, UpLeftPoint.X, point.Y);
  double sign = (SheetX2 > SheetX1)? +1.: -1.;
  //double dxmin = SheetX2-point.X;
  double dxmin = SheetX2-SheetX1;
  // for check of horizontal side
  double dY = 1.; // mm

  for (int nt=0; nt<LeftPart->Count(); nt++) {
    TileSide *tile = (*LeftPart)[nt]->tile;
    Point p1 = tile->ClosestEndY(UpLeftPoint);
    Point p2 = tile->ClosestEndY(DownRightPoint);
    double ymin = (p1.Y < p2.Y)? p1.Y: p2.Y;
    double ymax = (p1.Y > p2.Y)? p1.Y: p2.Y;
    ymin -= eps;
    ymax += eps;
    if (point.Y < ymin || point.Y > ymax) continue;
    for (int ns=0; ns<tile->sides->Count(); ns++) {
      Segment *side = (*tile->sides)[ns];
      if (fabs(side->p1.Y - side->p2.Y) < dY) continue; // horizontal side
      Point p = CrossPoint(line,*side);
      double dx = p.X - point.X;
      if (sign*dx < sign*dxmin) dxmin=dx;
    }
  }
  return dxmin;
}

void Megatile::Squeeze()
{
  if (LeftPart->Count()==0 || RightPart->Count()==0) return;

  double sign = (SheetX2 > SheetX1)? +1.: -1.;

  double dxmin = UpLeftPoint.X - DownRightPoint.X;

  int nt;
  for (nt=0; nt<LeftPart->Count(); nt++)
  {
    TileSide *tile = (*LeftPart)[nt]->tile;
    for (int ns=0; ns<tile->sides->Count(); ns++) {
      double dx;
      Point p;
      Side *side = (*tile->sides)[ns];
      p = side->p1;
      dx = dxRightPart(p);
      if (sign*dx < sign*dxmin) dxmin=dx;
      p = side->p2;
      dx = dxRightPart(p);
      if (sign*dx < sign*dxmin) dxmin=dx;
    }
  }

  for (nt=0; nt<RightPart->Count(); nt++)
  {
    TileSide *tile = (*RightPart)[nt]->tile;
    for (int ns=0; ns<tile->sides->Count(); ns++) {
      double dx;
      Point p;
      Side *side = (*tile->sides)[ns];
      p = side->p1;
      dx = dxLeftPart(p);
      if (sign*dx < sign*dxmin) dxmin=dx;
      p = side->p2;
      dx = dxLeftPart(p);
      if (sign*dx < sign*dxmin) dxmin=dx;
    }
  }

  double gap = 5.e-3;         // gap 5 micrometers
  //double gap = 0.;            // gap 0 micrometers
  dxmin += gap;

  for (nt=0; nt<LeftPart->Count(); nt++) (*LeftPart)[nt]->Move(-dxmin,0.);
}

void Megatile::ArrangeMegaSide1()
{
  ArrangeMegaSide1Right();
  ArrangeMegaSide1Left();
}

void Megatile::ArrangeMegaSide1Right()
{
  // forses to UpRight corner

  if (NameList->Count() == 0) return;

  const double eps = 1.e-7;

  // extention of SheetY2 to be sure to fit entire tile set
  double LargeY2 = SheetY1 + 1000.*(SheetY2-SheetY1);

  // fixed corners of sheet
  Point DownRightCorner(SheetX1, LargeY2);
  Point    UpLeftCorner(SheetX2, SheetY1);
  // for check of overlap even high tiles
  Point DownLeftCorner(SheetX2, LargeY2);

  // variable points to arrange
  // - all tiles with number < MaxHigh and odd tiles with number > MaxHigh
  Point DownRight;
  // - even tiles with number > MaxHigh
  Point   UpLeft;
  // for overlapping check
  Point DownLeft;

  double Rx = (SheetX2 > SheetX1)? +Dcut/2.: -Dcut/2.;
  double Ry = (SheetY2 > SheetY1)? +Dcut/2.: -Dcut/2.;
  DownRight=UpLeft=Point(SheetX1+xShift-Rx, SheetY1+yShift-Ry);

  int real;
  TileData *td;

  // check the number of the first tile
  int megatile,megaside,number;
  sscanf(NameList->Head(), "%d%*1s%1d%*1s%d", &megatile,&megaside,&number);

  int odd = number%2;
  if (number > 8 && !odd) number--;

  // position in NameList
  int the_name;
  // composed name (name only)
  char name[32];
  // tile line
  char *tline;

  while (number<=Tile::MaxTile(megatile,megaside))
  {
    sprintf(name, "%d-%d-%d", megatile, megaside, number);

    real = 0;
    tline = name;
    for(the_name=0; the_name<NameList->Count(); the_name++) {
      char line[32];
      strcpy(line,(*NameList)[the_name]);
      if(strcmp(strtok(line, " ,;\t"), name) == 0) {
        real = 1;
        tline = (*NameList)[the_name];
        break;
      }
    }

    if (number < Tile::MinHigh || number%2)
    {
      // create the tile at the point DownRight
      td = new TileData(tline, Dcut, DownRight.X, DownRight.Y, 180., 0);

      DownRight = td->tile->ClosestEndYX(DownRightCorner);
        UpLeft  = td->tile->ClosestEndYX(UpLeftCorner);
      // for overlapping check
      DownLeft  = td->tile->ClosestEndYX(DownLeftCorner);
    }
    else {
      // create the tile at the point UpLeft
      td = new TileData(tline, Dcut, UpLeft.X, UpLeft.Y, 180., 0);

      double delta = td->tile->ClosestEndYX(DownLeft).X - DownLeft.X;
      if (delta > eps)
      {
        Point shift = UpLeft;
        shift.X = shift.X - delta;
        td->Move(shift);
      }
    }
    double worklen = fabs((SheetY2-yShift)-(SheetY1+yShift));
    if (fabs(DownRight.Y-(SheetY1+yShift)) > worklen) {
      // out of sheet
      delete td; td=NULL;
      goto SHIFTRIGHT;
    }
    number++;
    //td->Show();
    if (real) {
      // register in tds
      *tds << td;
      // register in RightPart to make down cut
      *RightPart << td;
      // remove the name from NameList
      NameList->Kill(the_name);
      if (NameList->Count() == 0) goto SHIFTRIGHT;
    }
    else {
      delete td; td=NULL;
    }
  }

  SHIFTRIGHT:
  ;
  Point X1Y1(SheetX1+xShift, SheetY1+yShift);
  Point p = ContourClosestX(RightPart,X1Y1);
  double dx = p.X - X1Y1.X;
  double sign = (SheetX2 > SheetX1)? +1.: -1.;
  if (+sign*dx > eps) for (int j=0; j<RightPart->Count(); j++) {
    (*RightPart)[j]->Move(-dx,0.);
  }
}

void Megatile::ArrangeMegaSide1Left()
{
  // forses to DownLeft corner and then shifts set to top

  if (NameList->Count() == 0) return;

  const double eps = 1.e-7;

  // extention of SheetY1 to be sure to fit entire tile set
  double LargeY1 = SheetY2 - 1000.*(SheetY2-SheetY1);

  // fixed corners of sheet
  Point    UpLeftCorner(SheetX2, LargeY1);
  Point DownRightCorner(SheetX1, SheetY2);
  // for overlapping check of even high tiles
  Point UpRightCorner(SheetX1, SheetY1);

  // variable points to arrange
  // - all tiles with number < MaxHigh and odd tiles with number > MaxHigh
  Point UpLeft;
  // - even tiles with number > MaxHigh
  Point DownRight;
  // for overlapping check of even high tiles
  Point UpRight;

  //UpLeft=DownRight=Point(SignSub(SheetX2,xShift), SignSub(SheetY2,yShift));
  double Rx = (SheetX2 > SheetX1)? +Dcut/2.: -Dcut/2.;
  double Ry = (SheetY2 > SheetY1)? +Dcut/2.: -Dcut/2.;
  UpLeft=DownRight=Point(SheetX2-xShift+Rx, SheetY2-yShift+Ry);

  int real;
  Point upper; // upper point of real tile

  TileData *td;

  // check the number of the first tile
  int megatile,megaside,number;
  sscanf(NameList->Head(), "%d%*1s%1d%*1s%d", &megatile,&megaside,&number);

  int odd = number%2;
  if (number > 8 && !odd) number--;

  // position in NameList
  int the_name;
  // composed name (name only)
  char name[32];
  // tile line
  char *tline;

  while (number<=Tile::MaxTile(megatile,megaside))
  {
    sprintf(name, "%d-%d-%d", megatile, megaside, number);

    real = 0;
    tline = name;
    for(the_name=0; the_name<NameList->Count(); the_name++) {
      char line[32];
      strcpy(line,(*NameList)[the_name]);
      if(strcmp(strtok(line, " ,;\t"), name) == 0) {
        real = 1;
        tline = (*NameList)[the_name];
        break;
      }
    }

    if (number < Tile::MinHigh || number%2)
    {
      // create the tile at the point UpLeft
      td = new TileData(tline, Dcut, UpLeft.X, UpLeft.Y, 0., 0);

        UpLeft  = td->tile->ClosestEndYX(UpLeftCorner);
      DownRight = td->tile->ClosestEndYX(DownRightCorner);
      // for overlapping check of even high tiles
      UpRight   = td->tile->ClosestEndYX(UpRightCorner);
    }
    else
    {
      // create the tile at the point DownRight
      td = new TileData(tline, Dcut, DownRight.X, DownRight.Y, 0., 0);

      double delta = td->tile->ClosestEndYX(UpRight).X - UpRight.X;
      if (delta < -eps)
      {
        Point shift = DownRight;
        shift.X = shift.X - delta;
        td->Move(shift);
      }
    }
    number++;
    if (real) {
      // register in tds
      *tds << td;
      // register in LeftPart to make down cut
      *LeftPart << td;
      // remove the name from NameList
      NameList->Kill(the_name);
      // recalculate the upper point
      upper = td->tile->ClosestEnd(UpLeftCorner);
      if (NameList->Count() == 0) goto SHIFTUP;
    }
    else {
      delete td; td=NULL;
    }
  }

  SHIFTUP:
  // shift real tiles to upper position
  double dy = SheetY1+yShift-Ry - upper.Y;
  for (int i=0; i<LeftPart->Count(); i++) {
    TileData *td = (*LeftPart)[i];
    //cout << "ShiftUp: Before Move: "; td->Show();
    td->Move(0.,dy);
    //cout << "ShiftUp:  After Move to dy=" << dy; td->Show(); WAIT();
  }
}

void Megatile::ArrangeMegaSide2() {
  ArrangeMegaSide2Left();
  ArrangeMegaSide2Right();
}

void Megatile::ArrangeMegaSide2Left()
{
  // forses to UpLeft corner

  if (NameList->Count() == 0) return;

  const double eps = 1.e-7;

  // extention of SheetY2 to be sure to fit entire tile set
  double LargeY2 = SheetY1 + 1000.*(SheetY2-SheetY1);

  // fixed corners of sheet
  Point DownLeftCorner(SheetX2, LargeY2);
  Point  UpRightCorner(SheetX1, SheetY1);
  // for overlapping check of even tiles
  Point DownRightCorner(SheetX1, SheetY2);

  // variable points to arrange
  // - all tiles with number < MaxHigh and odd tiles with number > MaxHigh
  Point DownLeft;
  // - even tiles with number > MaxHigh
  Point   UpRight;
  // for overlapping check of even tiles
  Point DownRight;

  double Rx = (SheetX2 > SheetX1)? +Dcut/2.: -Dcut/2.;
  double Ry = (SheetY2 > SheetY1)? +Dcut/2.: -Dcut/2.;
  DownLeft=UpRight=Point(SheetX2-xShift+Rx, SheetY1+yShift-Ry);

  int real;

  TileData *td;

  // check the number of the first tile
  char *tname = NameList->Head();
  int megatile,megaside,number;
  sscanf(tname, "%d%*1s%1d%*1s%d", &megatile,&megaside,&number);

  int odd = number%2;
  if (number > 8 && !odd) number--;

  // position in NameList
  int the_name;
  // composed name (name only)
  char name[32];
  // tile line
  char *tline;

  while (number<=Tile::MaxTile(megatile,megaside))
  {
    sprintf(name, "%d-%d-%d", megatile, megaside, number);

    real = 0;
    tline = name;
    for(the_name=0; the_name<NameList->Count(); the_name++) {
      char line[32];
      strcpy(line,(*NameList)[the_name]);
      if(strcmp(strtok(line, " ,;\t"), name) == 0) {
        real = 1;
        tline = (*NameList)[the_name];
        break;
      }
    }

    if (number < Tile::MinHigh || number%2)
    {
      // create the tile at the point DownLeft
      td = new TileData(tline, Dcut, DownLeft.X, DownLeft.Y, 180., 0);

      DownLeft  = td->tile->ClosestEndYX(DownLeftCorner);
        UpRight = td->tile->ClosestEndYX(UpRightCorner);
      // for overlapping check of even tiles
      DownRight = td->tile->ClosestEndYX(DownRightCorner);
    }
    else {
      // create the tile at the point UpRight
      td = new TileData(tline, Dcut, UpRight.X, UpRight.Y, 180., 0);

      double delta = td->tile->ClosestEndYX(DownRight).X - DownRight.X;
      if (delta < -eps)
      {
        Point shift = UpRight;
        shift.X = shift.X - delta;
        td->Move(shift);
      }
    }
    double worklen = fabs((SheetY2-yShift) - (SheetY1+yShift));
    if (fabs(DownLeft.Y-(SheetY1+yShift)) > worklen) {
      // out of sheet
      delete td; td=NULL;
      goto SHIFTLEFT;
    }
    number++;

    if (real) {
      // register in tds
      *tds << td;
      // register in LeftPart to make down cut
      *LeftPart << td;
      // remove the name from NameList
      NameList->Kill(the_name);
      if (NameList->Count() == 0) goto SHIFTLEFT;
    }
    else {
      delete td; td=NULL;
    }
  }

  SHIFTLEFT:
  ;
  Point X2Y1(SheetX2-xShift, SheetY1+yShift);
  Point p = ContourClosestX(LeftPart,X2Y1);
  double dx = p.X-X2Y1.X;
  double sign = (SheetX2 > SheetX1)? +1.: -1.;
  if (-sign*dx > eps) for (int j=0; j<LeftPart->Count(); j++) {
    (*LeftPart)[j]->Move(-dx,0.);
  }
}

void Megatile::ArrangeMegaSide2Right()
{
  // forses to DownRight corner and then shifts set to top

  if (NameList->Count() == 0) return;

  const double eps = 1.e-7;

  // extention of SheetY1 to be sure to fit entire tile set
  double LargeY1 = SheetY2 - 1000.*(SheetY2-SheetY1);

  // fixed corners of sheet
  Point    UpRightCorner(SheetX1, LargeY1);
  Point   DownLeftCorner(SheetX2, SheetY2);
  // for overlapping check of even tiles
  Point UpLeftCorner(SheetX2, LargeY1);

  // variable points to arrange
  // - all tiles with number < MaxHigh and odd tiles with number > MaxHigh
  Point UpRight;
  // - even tiles with number > MaxHigh
  Point DownLeft;
  // for overlapping check of even tiles
  Point UpLeft;

  double Rx = (SheetX2 > SheetX1)? +Dcut/2.: -Dcut/2.;
  double Ry = (SheetY2 > SheetY1)? +Dcut/2.: -Dcut/2.;
  UpRight=DownLeft=Point(SheetX1+xShift-Rx, SheetY2-yShift+Ry);

  int real;
  Point upper; // upper point of real tile

  TileData *td;

  // check the number of the first tile
  char *tname = NameList->Head();
  int megatile,megaside,number;
  sscanf(tname, "%d%*1s%1d%*1s%d", &megatile,&megaside,&number);

  int odd = number%2;
  if (number > 8 && !odd) number--;

  // position in NameList
  int the_name;
  // composed name (name only)
  char name[32];
  // tile line
  char *tline;

  while (number<=Tile::MaxTile(megatile,megaside))
  {
    sprintf(name, "%d-%d-%d", megatile, megaside, number);

    real = 0;
    tline = name;
    for(the_name=0; the_name<NameList->Count(); the_name++) {
      char line[32];
      strcpy(line,(*NameList)[the_name]);
      if(strcmp(strtok(line, " ,;\t"), name) == 0) {
        real = 1;
        tline = (*NameList)[the_name];
        break;
      }
    }

    if (number < Tile::MinHigh || number%2)
    {
      // create the tile at the point UpRight
      td = new TileData(tline, Dcut, UpRight.X, UpRight.Y, 0., 0);

        UpRight = td->tile->ClosestEndYX(UpRightCorner);
      DownLeft  = td->tile->ClosestEndYX(DownLeftCorner);
      // for overlapping check of even tiles
        UpLeft  = td->tile->ClosestEndYX(UpLeftCorner);
    }
    else
    {
      // create the tile at the point DownLeft
      td = new TileData(tline, Dcut, DownLeft.X, DownLeft.Y, 0., 0);

      double delta = td->tile->ClosestEndYX(UpLeft).X - UpLeft.X;
      if (delta > eps)
      {
        Point shift = DownLeft;
        shift.X = shift.X - delta;
        td->Move(shift);
      }
    }
    number++;

    if (real) {
      // register in tds
      *tds << td;
      // register in RightPart to make down cut
      *RightPart << td;
      // remove the name from NameList
      NameList->Kill(the_name);
      // recalculate the upper point
      upper = td->tile->ClosestEnd(UpRightCorner);
      if (NameList->Count() == 0) goto SHIFTUP;
    }
    else {
      delete td; td=NULL;
    }
  }

  SHIFTUP:
  // shift real tiles to upper position
  double dy = SheetY1+yShift-Ry - upper.Y;
  for (int i=0; i<RightPart->Count(); i++) {
    TileData *td = (*RightPart)[i];
    //cout << "ShiftUp: Before Move: "; td->Show();
    td->Move(0.,dy);
    //cout << "ShiftUp:  After Move to d=" << d; td->Show(); WAIT();
  }
}

Side* Megatile::DownSide(Queue<TileData> *tdates)
{
  const double eps = 1e-7;
  double r2p1,r2p2, r2p1curr,r2p2curr;
  //Side down(SheetX1,SheetY2,SheetX2,SheetY2);
  //double LargeY2 = SignAdd(SheetY2, 100.*fabs(SheetY2-SheetY1));
  double LargeY2 = SheetY2 + 100.*(SheetY2-SheetY1);
  Side down(SheetX1,LargeY2,SheetX2,LargeY2);

  //double dx = SheetX2-SheetX1;
  //Side down(SheetX1-100.*dx,LargeY2,SheetX2+100.*dx,LargeY2);

  //r2p1 = r2p2 = 1000.*SQR(SheetY2-SheetY1);
  r2p1 = r2p2 = 1000.*SQR(LargeY2-SheetY1);
  //Side *closest_side = tdates->Head()->tile->ClosestSide(down,r2p1,r2p2);
  Side *closest_side;
  for (int i=0; i<tdates->Count(); i++)
  {
    Side *side = (*tdates)[i]->tile->ClosestSide(down,r2p1curr,r2p2curr);
    //if (r2p1curr <= r2p1 && r2p2curr <= r2p2) {
    if (r2p1-r2p1curr > -eps && r2p2-r2p2curr > -eps) {
      closest_side = side;
      r2p1 = r2p1curr;
      r2p2 = r2p2curr;
    }
  }
  return closest_side;
}

TileData* Megatile::ClosestTileX(const Queue<TileData> *list, const Point& corner)
{
  TileData *td;
  Point  pmin = list->Head()->tile->ClosestEndX(corner);
  double dmin = fabs(pmin.X - corner.X);
  for (int nt=1; nt<list->Count(); nt++) {
    td = (*list)[nt];
    Point p = td->tile->ClosestEndX(corner);
    double d = fabs(p.X-corner.X);
    if (d < dmin) {
      pmin = p;
      dmin = d;
    }
  }
  return td;
}

Point Megatile::ClosestX(const Queue<TileData> *list, const Point& corner)
{
  Point  pmin = list->Head()->tile->ClosestEndX(corner);
  double dmin = fabs(pmin.X - corner.X);
  for (int nt=1; nt<list->Count(); nt++) {
    Point p = (*list)[nt]->tile->ClosestEndX(corner);
    double d = fabs(p.X-corner.X);
    if (d < dmin) {
      pmin = p;
      dmin = d;
    }
  }
  return pmin;
}

Point Megatile::ClosestY(const Queue<TileData> *list, const Point& corner)
{
  Point  pmin = list->Head()->tile->ClosestEndY(corner);
  double dmin = fabs(pmin.Y - corner.Y);
  for (int nt=1; nt<list->Count(); nt++) {
    Point p = (*list)[nt]->tile->ClosestEndY(corner);
    double d = fabs(p.Y-corner.Y);
    if (d < dmin) {
      pmin = p;
      dmin = d;
    }
  }
  return pmin;
}

Point Megatile::ClosestXY(const Queue<TileData> *list, const Point& corner)
{
  const double eps = 1.e-7;
  Point  pmin  = list->Head()->tile->ClosestEndXY(corner);
  double dxmin = fabs(pmin.X - corner.X);
  double dymin = fabs(pmin.Y - corner.Y);
  for (int nt=1; nt<list->Count(); nt++) {
    Point p = (*list)[nt]->tile->ClosestEndXY(corner);
    double dx = fabs(p.X-corner.X);
    if (dx < dxmin+eps) {
      // check Y distance
      double dy = fabs(p.Y-corner.Y);
      if (dy < dymin+eps) {
        pmin  = p;
        dxmin = dx;
        dymin = dy;
      }
    }
  }
  return pmin;
}

Point Megatile::ClosestYX(const Queue<TileData> *list, const Point& corner)
{
  const double eps = 1.e-7;
  Point  pmin  = list->Head()->tile->ClosestEndYX(corner);
  double dxmin = fabs(pmin.X - corner.X);
  double dymin = fabs(pmin.Y - corner.Y);
  for (int nt=1; nt<list->Count(); nt++) {
    Point p = (*list)[nt]->tile->ClosestEndYX(corner);
    double dy = fabs(p.Y-corner.Y);
    if (dy < dymin+eps) {
      // check X distance
      double dx = fabs(p.X-corner.X);
      if (dx < dxmin+eps) {
        pmin  = p;
        dxmin = dx;
        dymin = dy;
      }
    }
  }
  return pmin;
}

Point Megatile::ContourClosestX(const Queue<TileData> *list, const Point& corner)
{
  Point  pmin = list->Head()->tile->ContourClosestEndX(corner);
  double dmin = fabs(pmin.X - corner.X);
  for (int nt=1; nt<list->Count(); nt++) {
    Point p = (*list)[nt]->tile->ContourClosestEndX(corner);
    double d = fabs(p.X-corner.X);
    if (d < dmin) {
      pmin = p;
      dmin = d;
    }
  }
  return pmin;
}

void Megatile::MakeDownCut()
{
  if (down_cut == 0) return;

  Point p1,p2;

  /*
  p1 = ClosestX(tds, DownRightPoint);
  p2 = ClosestY(tds, DownRightPoint);
  // redefine the DownRightPoint
  DownRightPoint = Point(p1.X, p2.Y);

  p1 = ClosestX(tds, UpLeftPoint);
  p2 = ClosestY(tds, UpLeftPoint);
  // redefine the UpLeftPoint
  UpLeftPoint = Point(p1.X, p2.Y);
  */

  double ydown = DownRightPoint.Y;

  p1 = Point(DownRightPoint.X, ydown);
  p2 = Point(   UpLeftPoint.X, ydown);

  Point  start(SheetX1-ddown, ydown);
  Point finish(SheetX2+ddown, ydown);

  *down_cuts << new Side(start,finish);

  if (down_cut > 1)
  {
    // back to p2
    //if (Distance(start,p2) < 2.*fabs(ddown))
    //     *down_cuts << new Side(start,p2);
    //else *down_cuts << new Side(start,p2,lift);
    *down_cuts << new Side(start,p2,lift);

    p1 = p2;
    p2 = UpLeftPoint;
    *down_cuts << new Side(p1,p2);

    p1 = p2;
    p2 = Point(DownRightPoint.X, UpLeftPoint.Y);
    *down_cuts << new Side(p1,p2);

    p1 = p2;
    finish = DownRightPoint;
    *down_cuts << new Side(p1,finish);
  }

  // move cutter to start
  *down_cuts << new Side(finish,start, lift);

  // strip holes from down tiles
  for (int i=0; i<down_cuts->Count(); i++) {
    Side *side = (*down_cuts)[i];
    for (int ntd=0; ntd<tds->Count(); ntd++) {
      TileData *td = (*tds)[ntd];
      if (td->status == todo) {
        if (td->tile->LetStrip(side))
        {
          CheckStatus(td);
        }
      }
    }
  }

  if (down_todo == 0) {
    for (int i=0; i<down_cuts->Count(); i++) (*down_cuts)[i]->status=done;
  }
}

/*
void Megatile::MakeDownCut()
{
  if (down_cut == 0) return;

  Side *down1=NULL, *down2=NULL;
  if (RightPart->Count() > 0) down1 = DownSide(RightPart);
  if (LeftPart->Count()  > 0) down2 = DownSide(LeftPart);

  //double LargeY2 = SignAdd(SheetY2, 100.*fabs(SheetY2-SheetY1));
  double LargeY2 = SheetY2 + 100.*(SheetY2-SheetY1);
  //Side edge1(SignSub(SheetX1, Dcut), SignSub(SheetY1, Dcut),
  //           SignSub(SheetX1, Dcut), LargeY2);
  //Side edge2(SignAdd(SheetX2, Dcut), SignSub(SheetY1, Dcut),
  //           SignAdd(SheetX2, Dcut), LargeY2);
  double d = (ddown > 0.)? ddown+Dcut/2.: ddown-Dcut/2.;
  Side edge1(SheetX1-d, SheetY1, SheetX1-d, LargeY2);
  Side edge2(SheetX2+d, SheetY1, SheetX2+d, LargeY2);
  Point p1,p2;
  double r2p1,r2p2;

  Point start;

  if (down1 && down2)
  {
    Point p;
    // find closest down1 point to edge1
    r2p1 = R2(down1->p1,edge1.ClosestPoint(down1->p1));
    p1 = down1->p1;
    p2 = down1->p2;
    p  = p2;
    r2p2 = R2(down1->p2,edge1.ClosestPoint(down1->p2));
    if (r2p2 < r2p1) {
      p1 = down1->p2;
      p2 = down1->p1;
      p  = p1;
    }

    // start point of tile
    start = p1;

    // cut from start of sheet
    *down_cuts << new Side(edge1.ClosestPoint(p1), p1);
    // cut down1
    *down_cuts << new Side(p1,p2);
    // find closest down2 point to edge1
    r2p1 = R2(down2->p1,edge1.ClosestPoint(down2->p1));
    p1 = down2->p1;
    p2 = down2->p2;
    r2p2 = R2(down2->p2,edge1.ClosestPoint(down2->p2));
    if (r2p2 < r2p1) {
      p1 = down2->p2;
      p2 = down2->p1;
    }
    // set a bridge
    *down_cuts << new Side(p,p1);

    // cut down2
    *down_cuts << new Side(p1,p2);
    // cut to end of sheet
    *down_cuts << new Side(p2,edge2.ClosestPoint(p2));
  }
  else {
    if (down1) {
      // find closest down1 point to edge1
      r2p1 = R2(down1->p1,edge1.ClosestPoint(down1->p1));
      p1 = down1->p1;
      p2 = down1->p2;
      r2p2 = R2(down1->p2,edge1.ClosestPoint(down1->p2));
      if (r2p2 < r2p1) {
        p1 = down2->p2;
        p2 = down2->p1;
      }

      // start point of tile
      start = p1;

      // cut from start of sheet
      *down_cuts << new Side(edge1.ClosestPoint(p1), p1);
      // cut down1
      *down_cuts << new Side(p1,p2);
      // cut to end of sheet
      *down_cuts << new Side(p2, edge2.ClosestPoint(p2));
    }
    else {
      // find closest down2 point to edge1
      r2p1 = R2(down2->p1,edge1.ClosestPoint(down2->p1));
      p1 = down2->p1;
      p2 = down2->p2;
      r2p2 = R2(down2->p2,edge1.ClosestPoint(down2->p2));
      if (r2p2 < r2p1) {
        p1 = down2->p2;
        p2 = down2->p1;
      }

      // start point of tile
      start = p1;

      // cut from start of sheet
      *down_cuts << new Side(edge1.ClosestPoint(p1), p1);
      // cut down2
      *down_cuts << new Side(p1,p2);
      // cut to end of sheet
      *down_cuts << new Side(p2, edge2.ClosestPoint(p2));
    }
  }

  // clear lists RightPart and LeftPart
  while (RightPart->Count()) RightPart->Remove(0);
  while ( LeftPart->Count())  LeftPart->Remove(0);

  if (down_cut > 1)
  {
    // make a sheet cut
    double xshift = SignSub(xShift, Dcut/2.);
    double yshift = SignSub(yShift, Dcut/2.);
    // move(cut) to end of area
    Point DownLeft(SheetX2-xshift, down_cuts->Tail()->p1.Y);
    *down_cuts << new Side(down_cuts->Tail()->p2, DownLeft);
    // cut to top
    Point UpLeft(SheetX2-xshift, SheetY1+yshift);
    *down_cuts << new Side(DownLeft, UpLeft);
    // cut top
    Point UpRight(SheetX1+xshift, SheetY1+yshift);
    *down_cuts << new Side(UpLeft, UpRight);
    // cut to bottom
    // point start can be from the left tile
    *down_cuts << new Side(UpRight, Point(UpRight.X,start.Y));
  }

  // strip holes from down tiles
  for (int i=0; i<down_cuts->Count(); i++) {
    Side *side = (*down_cuts)[i];
    for (int ntd=0; ntd<tds->Count(); ntd++) {
      TileData *td = (*tds)[ntd];
      if (td->status == todo) {
        if (td->tile->LetStrip(side))
        {
          CheckStatus(td);
        }
      }
    }
  }
}
*/
