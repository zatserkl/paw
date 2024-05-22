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

#if !defined(__MATH_H)
#include <math.h>
#endif

#if !defined(__IOSTREAM_H)
#include <iostream.h>
#endif

#if !defined(__FSTREAM_H)
#include <fstream.h>
#endif

// data structure
#include "queue.h"
#include "sidelist.h"
#include "elements.h"

#include "tile.h"      /* for size of class Tile */

#include "wait.h"

#include "megatile.h"

// general-purpose functions

int ToCut(const Queue<End> *list)
{
  for (int i=0; i<list->Count(); i++) if ((*list)[i]->status==tocut) return 1;
  return 0;
}

int ToCut(const Queue<EndGroove> *list)
{
  for (int i=0; i<list->Count(); i++) if ((*list)[i]->status==tocut) return 1;
  return 0;
}

// struct End

// declaration and definition of static variable
short End::EndAnchor; // End key

End::End(const Point& p) {
  X = p.X;
  Y = p.Y;
  Dcutter = 0.;
  cutype  = cutcut;
  status  = tofly;
  tile    = NULL;
}

End::End(const Point& p, double dcutter) {
  X = p.X;
  Y = p.Y;
  Dcutter = dcutter;
  cutype  = cutcut;
  status  = tofly;
  tile    = NULL;
}

End::End(const Point& p, EndStatus thestatus) {
  X = p.X;
  Y = p.Y;
  Dcutter = 0.;
  cutype  = cutcut;
  status  = thestatus;
  tile    = NULL;
}

End::End(const Point& p, double dcutter, EndStatus stat) {
  X = p.X;
  Y = p.Y;
  Dcutter = dcutter;
  cutype  = cutcut;
  status  = stat;
  tile    = NULL;
}

End::End(const Point& p, double dcutter, EndStatus stat, TileSide *master) {
  X = p.X;
  Y = p.Y;
  Dcutter = dcutter;
  cutype  = cutcut;
  status  = stat;
  tile    = master;
}

End::End(const Point& p, double dcutter, EndCutter ecut, EndStatus stat) {
  X = p.X;
  Y = p.Y;
  Dcutter = dcutter;
  cutype  = ecut;
  status  = stat;
  tile    = NULL;
}

End::End(const Point& p, double dcutter, EndCutter ecut, EndStatus stat,
         TileSide *master)
{
  X = p.X;
  Y = p.Y;
  Dcutter = dcutter;
  cutype  = ecut;
  status  = stat;
  tile    = master;
}

void End::Show() const {
  cout.precision(3);
  cout.width(12);
  cout << X << " ";
  cout.precision(3);
  cout.width(12);
  cout << Y << " ";
  cout.precision(3);
  cout.width(7);
  cout << Dcutter << " ";
  cout << cutype;
  cout << " ";
  cout << status;
}

// struct EndGroove

// declaration and definition of static variable
short EndGroove::EndGrooveAnchor; // EndGroove key

void EndGroove::Show() const {
  End::Show();
  cout << " Ox, Oy, TableAngle" << endl;
  P->Head()->Show();
  cout << " P0, yMAHO" << endl;
  for (int i=1; i<P->Count(); i++)
  {
    (*P)[i]->Show();
    cout << " P" << i+1 << endl;
  }
}

// struct is used to select tile to cut
struct R2TileSidePoint {
  // square of distance to pmin
  double R2cut;
  // square of distance to current cutter position
  double R2cutter;
  // closest point at the nearest done cut
  Point pmin;

  TileData *td;
  Side     *side;
  Point    *p;
};

// class Megatile

// declare static variable for class Megatile //
// -------------------------------------------//
double Megatile::Dcut;
double Megatile::Dhole;
char   Megatile::task_file[33];
char   Megatile::path_file[33];
char   Megatile::groove_file[33];

double Megatile::SheetX1;
double Megatile::SheetY1;
double Megatile::SheetX2;
double Megatile::SheetY2;
// shift with respect to (SheetX1,SheetY1)
double Megatile::xShift;
double Megatile::yShift;

// is turn MAHO table to real tile groove angle
int    Megatile::TrueAngle;
int    Megatile::SqueezeParts;
int    Megatile::BackSide;

int    Megatile::down_cut  = -1; // known to be wrong initial value
int    Megatile::down_todo =  1; // means todo
// exstra space for down cut
double Megatile::ddown;

// optional show of turned contour
int Megatile::ShowTurnedContour = 0;

// safety gates
Queue<Side> *Megatile::gates;

// tiles
Queue<char> *Megatile::NameList = NULL;

// start after break
char   Megatile::continue_tile[33] = "";
int    Megatile::continue_level = 0;
// -------------------------------------------//

Megatile::Megatile(double Ox, double Oy) : O(Ox,Oy)
{
  // set static Tile::Dsemihole to Megatile::Dcut
  Tile::Dsemihole = Dcut;

  // sign of xShift
  if (SheetX2 > SheetX1) {
    xShift = +fabs(xShift);
    ddown  = +fabs(ddown);
  }
  else {
    xShift = -fabs(xShift);
    ddown  = -fabs(ddown);
  }
  // sign of yShift
  if (SheetY2 > SheetY1) yShift = +fabs(yShift);
  else                   yShift = -fabs(yShift);

  DownRightPoint = Point(SheetX1,SheetY2);
  UpLeftPoint    = Point(SheetX2,SheetY1);

  nopath = 0;

  tds = new Queue<TileData>;
  cuts      = new SideList;
  down_cuts = new SideList;

  ends        = new Queue<End>;
  down_ends   = new Queue<End>;
  groove_ends = new Queue<EndGroove>;

  RightPart = new Queue<TileData>;
  LeftPart  = new Queue<TileData>;

  continue_cut = NULL;
}

Megatile::~Megatile()
{
  // RightPart and LeftPart contain the pointers to td
  // which are use by Squeeze() and old MakeDownCut()
  if (RightPart) {
    while (RightPart->Count()) RightPart->Remove(0);
    delete RightPart; RightPart=NULL;
  }
  if (LeftPart) {
    while ( LeftPart->Count())  LeftPart->Remove(0);
    delete  LeftPart;  LeftPart=NULL;
  }

  if (tds != NULL) {
    // remove from the cuts the tiles sides
    for (int ntd=0; ntd<tds->Count(); ntd++)
    {
      TileSide *tile = (*tds)[ntd]->tile;
      for (int ns=0; ns<tile->sides->Count(); ns++)
      {
        Side *side = (*tile->sides)[ns];
        if (cuts->Exist(side)) cuts->Remove(side);
      }
    }
    delete tds; tds=NULL;
  }

  delete      cuts;      cuts=NULL;
  delete down_cuts; down_cuts=NULL;

  delete ends;        ends=NULL;
  delete down_ends;   down_ends=NULL;
  delete groove_ends; groove_ends=NULL;

  delete gates; gates=NULL;

  delete NameList; NameList=NULL;
  //MESS("~Megatile done");
}

void Megatile::WhoAlive() {
  cout << "Alive Points "   << Point::count   << endl;
  cout << "Alive Segments " << Segment::count << endl;
  cout << "Alive Sides "    << Side::count    << endl;
  cout << "Alive Tiles "    << Tile::count    << endl;
}

void Megatile::WhoAlive(char *mess) {
  if (strlen(mess)) cout << mess << endl;
  cout << "Alive Points "   << Point::count   << endl;
  cout << "Alive Segments " << Segment::count << endl;
  cout << "Alive Sides "    << Side::count    << endl;
  cout << "Alive Tiles "    << Tile::count    << endl;
}

void Megatile::Insert(char *tname, double x, double y, double thdeg, int inv)
{
  *tds << new TileData(tname, Dcut, x,y,thdeg,inv);
}

void Megatile::Insert(char *tname, double thdeg, int inv)
{
  TileData *td = new TileData(tname, Dcut, thdeg,inv);
  *tds << td;
}

void Megatile::Insert(char *tname)
{
  int megatile,megaside,number;
  sscanf(tname, "%d%*1s%1d%*1s%d", &megatile,&megaside,&number);

  if (number > Tile::MaxTile(megatile,megaside))
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
    cerr << "exeeds the max tile number in the megatile " << megatile << endl;
    STOPCR();
  }

  // store only tile's name in NameList, without creating of tile
  char *namestr = new char[strlen(tname)+1];
  strcpy(namestr,tname);

  if (NameList->Count() == 0) *NameList << namestr;
  else {
    char *name;
    int mtile,side,num;

    // check the megaside
    name = (*NameList)[0];
    sscanf(name, "%d%*1s%1d%*1s%d", &mtile,&side,&num);
    if (megatile != mtile || megaside != side)
    {
      cerr << "Sorry, different megatiles:" << endl;
      STOPCR("Program processes tiles from the same megatile!");
    }

    // store in increasing order
    for (int i=0; i<NameList->Count(); i++)
    {
      name = (*NameList)[i];
      sscanf(name, "%d%*1s%1d%*1s%d", &mtile,&side,&num);
      if (num > number) {
        NameList->InsertBefore(i,namestr);
        return;
      }
    }
    *NameList << namestr;
  }
}

void Megatile::Show() const
{
  for (int i=0; i<tds->Count(); i++) (*tds)[i]->Show();
}

void Megatile::WriteEnds()
{
  if (strlen(path_file) == 0) return;
  ofstream fout(path_file);
  fout << "* Task file " << task_file << endl;
  fout << "*                                             " << endl;
  fout << "* Ends in program frame                       " << endl;
  fout << "* Field status: 0=tofly, 1=tocut, 2=tosemihole" << endl;
  fout << "*                                             " << endl;
  fout << "*   end->X       end->Y   Dcutter status      " << endl;
  fout << "*                                             " << endl;
  int i;
  fout << "down_ends:" << endl;
  for (i=0; i<down_ends->Count(); i++) {
    End *end = (*down_ends)[i];
    fout.precision(3);
    fout.width(12);
    fout << end->X << " ";
    fout.precision(3);
    fout.width(12);
    fout << end->Y << " ";
    fout.precision(3);
    fout.width(7);
    fout << end->Dcutter << " ";
    fout << end->status << endl;
  }
  fout << "ends:" << endl;
  for (i=0; i<ends->Count(); i++) {
    End *end = (*ends)[i];
    fout.precision(3);
    fout.width(12);
    fout << end->X << " ";
    fout.precision(3);
    fout.width(12);
    fout << end->Y << " ";
    fout.precision(3);
    fout.width(7);
    fout << end->Dcutter << " ";
    fout << end->status  << " ";
    char *name = GetName(end->tile);
    if (name) fout << name;
    else      fout << "SHEET";
    fout << endl;
  }
  fout.close();
}

void Megatile::WriteGrooves()
{
  if (strlen(groove_file) == 0) return;
  ofstream fout(groove_file);
  fout << "* Task file " << task_file << endl;
  fout << "*                                            " << endl;
  fout << "* Grooves and contour in program frame (X,Y) " << endl;
  fout << "*                                            " << endl;
  for (int ntile=0; ntile<tds->Count(); ntile++)
  {
    Tile *tile = (*tds)[ntile]->tile;
    int i;
    fout << "Tile " << tile->nameAB << endl;
    fout << "contour:" << endl;
    for (i=0; i<tile->contour->Count(); i++) {
      Side *side = (*tile->contour)[i];
      fout.precision(3);
      fout << "(";
      fout.width(12);
      fout << side->p1.X << " ";
      fout.width(12);
      fout << side->p1.Y << ")";
      fout.width(7);

      fout << "--> ";

      fout << "(";
      fout.width(12);
      fout << side->p2.X << " ";
      fout.width(12);
      fout << side->p2.Y << " ";
      fout << ")";
      fout << endl;
    }
    fout << "groove:" << endl;
    int iP=0, iR=0;

    // initial point
    Hole *p = (Hole*) (*tile->groove)[0];
    fout << "P" << iP++ << " = ";
    fout << "(";
    fout.precision(3);
    fout.width(12);
    fout << p->X << " ";
    fout.width(12);
    fout << p->Y << " ";
    fout.width(7);
    fout << ")";
    fout.width(12);
    fout << " MAHO Y = ";
    fout.width(12);
    fout << p->D;
    fout << endl;

    for (i=1; i<tile->groove->Count(); i++)
    {
      Point *p = (*tile->groove)[i];
      if (p->GetKey() == Point::Key()) {
        fout << "P" << iP++ << " = ";
        fout << "(";
        fout.width(12);
        fout << p->X << " ";
        fout.width(12);
        fout << p->Y << " ";
        fout.width(7);
        fout << ")";
      }
      else {
        Hole *p = (Hole*) (*tile->groove)[i];
        fout << "R" << ++iR << " = ";
        fout << "(";
        fout.width(12);
        fout << p->X << " ";
        fout.width(12);
        fout << p->Y << " ";
        fout.width(7);
        fout << ")";
        fout.width(12);
        fout << " R = ";
        fout.width(12);
        fout << p->D;
      }
      fout << endl;
    }
  }
  fout.close();
}

void Megatile::PrintCuts() {
  cout << "The whole path (without semoholes)" << endl;
  int i;
  cout << "down cuts:" << endl;
  for (i=0; i<down_cuts->Count(); i++) {
    Side *cut = (*down_cuts)[i];
    //(*down_cuts)[i]->Show();
    cut->Show();
    cout << endl;
  }
  cout << "cuts:" << endl;
  for (i=0; i<cuts->Count(); i++) {
    (*cuts)[i]->Show();
    cout << endl;
  }
}

void Megatile::PrintEnds() {
  int i;
  cout << "down ends:" << endl;
  for (i=0; i<down_ends->Count(); i++) {
    (*down_ends)[i]->Show();
    cout << endl;
  }
  cout << "ends:" << endl;
  for (i=0; i<ends->Count(); i++) {
    (*ends)[i]->Show();
    cout << endl;
  }
}

void Megatile::PrintCutsEnds() {
  PrintCuts();
  WAIT();
  PrintEnds();
}

int Megatile::CheckStatus(TileData *td) {
  //
  // return the number of elements 'todo',
  // set status to 'done' if all elements are 'todo'
  //
  int ntodo = 0;
  int nes = td->tile->sides->Count();
  for (int ne=0; ne<nes; ne++) if ((*td->tile->sides)[ne]->status==todo) ntodo++;
  if (ntodo == 0) td->status = done;
  return ntodo;
}

int Megatile::SelectTileToCut(R2TileSidePoint& tsp) {
  // returns 1 if found tile to cut, 0 otherwise
  // additional information available thru structure tsp
  int ncut,ncuts;
  int ntile,ntds;
  int nside,nsides;
  TileData *td;
  Side *side;
  Point *p;
  Point pCutter;
  const double eps = 1.e-7;

  // last cutter position
  ncuts = cuts->Count();
  //Point pCutter = (*cuts)[ncuts-1]->p2;
  pCutter = (*cuts)[ncuts-1]->p2;

  int selected = 0;

  double R2cut, R2cutter;
  tsp.R2cut = tsp.R2cutter = 1.e100; tsp.pmin = cuts->Head()->p2;

  ntds = tds->Count();
  for (ntile=0; ntile<ntds; ntile++) {
    td = (*tds)[ntile];
    if (td->status != todo) continue;
    selected = 1;
    nsides = td->tile->sides->Count();
    for (nside=0; nside<nsides; nside++) {
      side = (*td->tile->sides)[nside];
      if (side->status != todo) continue;
      // begin cycle thru cuts
      Point pmin;
      // side->p1
      p = &side->p1;
      for (ncut=0; ncut<ncuts; ncut++) {
        pmin = (*cuts)[ncut]->ClosestPoint(*p);
        R2cut = R2(&pmin,p);
        R2cutter = R2(&pmin,&pCutter);
        if (R2cut > tsp.R2cut+eps) continue;
        if (fabs(R2cut - tsp.R2cut) < eps
            && R2cutter > tsp.R2cutter+eps) continue;
        tsp.td = td;
        tsp.side = side;
        tsp.p    = p;
        tsp.pmin = pmin;
        tsp.R2cut = R2cut;
        tsp.R2cutter = R2cutter;
      }
      // side->p2
      p = &side->p2;
      for (ncut=0; ncut<ncuts; ncut++) {
        pmin = (*cuts)[ncut]->ClosestPoint(*p);
        R2cut = R2(&pmin,p);
        R2cutter = R2(&pmin,&pCutter);
        if (R2cut > tsp.R2cut+eps) continue;
        if (fabs(R2cut - tsp.R2cut) < eps
            && R2cutter > tsp.R2cutter+eps) continue;
        tsp.td = td;
        tsp.side = side;
        tsp.p    = p;
        tsp.pmin = pmin;
        tsp.R2cut = R2cut;
        tsp.R2cutter = R2cutter;
      }
      // end of cycle thru cuts for current side of current tile
    }
  }
  return selected;
}

void Megatile::CutPath() {
  int ntds, ntile;
  int nsides, nside;
  int ncuts,  ncut;
  const double eps = 1.e-7;
  TileData *td;
  Side *side;
  Point p=0.;
  Point pCutter;

  R2TileSidePoint tsp;

  // Let's consider sheet edges as a some kind of sides

  // add to cuts the copies of gates
  for(int ngate=0;ngate<gates->Count();ngate++) {
    *cuts << new Side(*(*gates)[ngate]);
  }

  if (down_cut)
  {
    MakeDownCut();

    // add to cuts a copy of down cuts (will be removed at the end)
    for (int i=0; i<down_cuts->Count(); i++) {
      *cuts << new Side(*(*down_cuts)[i]);
    }
  }
  else {
    // define the current cutter position by gate point closest to (0,0)
    Point p0(0.,0.);
    Point  p = (*gates)[0]->ClosestPoint(p0);
    double r = R2(p,p0);
    for (int i=1; i<gates->Count(); i++) {
      Point  pi = (*gates)[i]->ClosestPoint(p0);
      double ri = R2(pi,p0);
      if (ri < r) {
        p = pi;
        r = ri;
      }
    }
    *cuts << new Side(cuts->Tail()->p2, p, lift);
  }

  // All these dummy cuts will be removed at the end

  // the number of dummy cuts
  int ndummy = cuts->Count();

  // select the nearest side

  while (SelectTileToCut(tsp))
  {
    // Now point tsp.p of side tsp.side of tile tsp.td has a minimal distance
    // to just performed cuts and (desirable) to current cutter position.

    if (continue_level > 2) {
      if (strcmp(tsp.td->tile->name, continue_tile) == 0)
      {
        continue_cut = cuts->Tail();
        if (cuts->Count() == ndummy)
        {
          // it's a joke: user try to continue from the first tile!
          cout << endl << "WARNING: "
               << "You are trying to continue from the first tile to cut!"
               << endl;
          WAIT();
          // clear continue_cut but keep the current continue_level
          continue_cut     = NULL;
        }
      }
    }

    // Set this side the Head of tile

    tsp.td->tile->SetHead(tsp.side);

    if (tsp.p == &tsp.side->p2) {
      tsp.td->tile->sides->Reverse();

      // start point now is tsp.side->p1
      //tsp.p = &tsp.side->p1;
    }

    //tsp.td->Show(); MESS("before Optimize.");

    tsp.td->tile->Optimize();

    // start point of current tile after Optimize()
    tsp.p = &tsp.td->tile->sides->Head()->p1;

    //tsp.td->Show(); WAIT(" after Optimize.");

    // set a bridge: we have fly to pmin and cut to side point
    pCutter = cuts->Tail()->p2;

    if (Distance(pCutter, tsp.pmin) > eps) {
      // fly to pmin
      *cuts << new Side(pCutter,tsp.pmin, lift);
    }
    // the first point of current tile is *tsp.p
    if (Distance(tsp.pmin, *tsp.p) > eps) {
      // cut from pmin to side point
      *cuts << new Side(tsp.pmin, *tsp.p, todo);
    }

    // put pointers to all sides of choosen tile to cuts

    nsides = tsp.td->tile->sides->Count();
    for (nside=0; nside<nsides; nside++) {
      side = (*tsp.td->tile->sides)[nside];
      if (side->status == todo)
      {
        // strip holes from overlapped sides
        for (int i=0; i<tds->Count(); i++) {
          td = (*tds)[i];
          if (td == tsp.td) continue;
          if (td->status == todo) {
            if (td->tile->LetStrip(side))
            {
              CheckStatus(td);
            }
          }
        }
      }
      // add a pointer to the side to cuts
      if (side->status != extra) *cuts << side;
    }
    tsp.td->status = done;

    cout << "Processed tile " << tsp.td->tile->nameAB << endl;
  }

  // kill the dummy cuts
  for (int i=0; i<ndummy; i++) cuts->Kill(0);
}

TileSide* Megatile::Master(Side *side)
{
  for (int ntd=0; ntd<tds->Count(); ntd++) {
    TileSide *tile = (*tds)[ntd]->tile;
    for (int ns=0; ns<tile->sides->Count(); ns++) {
      if ((*tile->sides)[ns] == side) return tile;
    }
  }
  return NULL;
}

void Megatile::AddCutToEnds(Side *cut, Queue<End> *endlist, TileSide *master)
{
  EndStatus status = (cut->status == todo)? tocut: tofly;
  if (continue_cut)
  {
    // continue from last performed cut
    status=tofly;
    if (continue_cut == cut) continue_cut = NULL;
  }
  if (status == tocut) {
    if (cut->sholes->Count())
    {
      // parameters of cut line
      LinePar s(cut);

      int nsholes=cut->sholes->Count();
      for (int i=0; i<nsholes; i++) {
        Hole *hole = (*cut->sholes)[i];

        // local variables for cutting up to hole
        EndStatus cut_status, hole_status;
        double Dcutter;

        if (hole->GetKey() == SHole::Key()) {
          cut_status  = tocut;
          hole_status = tosemihole;
          Dcutter = hole->D;
        }
        else
        if (hole->GetKey() == AHole::Key()) {
          cut_status  = tofly;
          hole_status = tohole;
          Dcutter = 0.;
        }
        else {
          MESS("Megatile::MakePath(): Unknown object:");
          hole->Show();
          STOP("Don't know what to do!");
        }

        // parameter of hole base
        double t = s.cosx*(hole->X - s.x0) + s.cosy*(hole->Y - s.y0);

        Point base;
        base.X = s.x0 + s.cosx*t;
        base.Y = s.y0 + s.cosy*t;

        // cut up to base
        *endlist << new End( base,Dcutter,cut_status,master);

        // make a hole
        *endlist << new End(*hole,Dcutter,hole_status,master);

        // return back: for back from semihole we don't use tofly
        *endlist << new End( base,Dcutter,cut_status,master);
      }
    }
    // second end
    *endlist << new End(cut->p2, Dcut, status, master);
  }
  else *endlist << new End(cut->p2, 0., status, master);
}

void Megatile::DownEnds()
{
  if (down_cut == 0) return;

  if (continue_level > 0)
  {
    // mark all down cuts as done
    for (int i=0; i<down_cuts->Count(); i++) (*down_cuts)[i]->status = done;
  }

  // fly cutter to start of the first down cut
  *down_ends << new End((*down_cuts)[0]->p1, Dcut, tofly);

  // the rest of down_cuts
  //for (int i=1; i<down_cuts->Count(); i++) {   // old MakeDownCut()
  for (int i=0; i<down_cuts->Count(); i++) {
    AddCutToEnds((*down_cuts)[i],down_ends,NULL);
  }
}

void Megatile::CuttingEnds()
{
  int ncuts = cuts->Count();
  if (ncuts == 0) return;

  // fly cutter to start of the first cut
  //*ends << new End(cuts->Head()->p1, Dcut, tofly, Master(cuts->Head()));
  *ends << new End(cuts->Head()->p1, Dcut, tofly);

  // Include in ends only side's second point

  //for (int ncut=0; ncut<ncuts; ncut++) AddCutToEnds((*cuts)[ncut],ends);
  for (int ncut=0; ncut<ncuts; ncut++) {
    Side *cut = (*cuts)[ncut];
    AddCutToEnds(cut,ends);
  }

  // Accept to bridges the names of corresponding tiles
  for (int i=0; i<ends->Count(); i++) {
    End *end = (*ends)[i];
    if (end->tile == NULL) {
      for (int j=i+1; j<ends->Count(); j++) {
        End *next = (*ends)[j];
        if (next->tile != NULL) {
          end->tile = next->tile;
          break;
        }
      }
    }
  }
}

void Megatile::MakePath()
{
  // arrange tiles
  Arrange();

  // fill down_cuts end cuts
  CutPath();

  // compute the ends
  DownEnds();
  CuttingEnds();
  GroovingEnds();

  WriteGrooves();

  //PrintCutsEnds();
}

void Megatile::Rotate(const Point& O, double theta)
{
  for (int i=0; i<tds->Count(); i++) (*tds)[i]->Rotate(O,theta);
}

// For graphics package

Point Megatile::TileUpLeftPoint(TileSide *tile)
{
  /*
  Point UpLeftCorner(SheetX2,SheetY1);
  Point UpLeft = tile->ClosestEnd(UpLeftCorner);
  // this point can be at tile's down side,
  // so point stright under the tile is
  Point top(UpLeft.X, UpLeftCorner.Y);
  // and to be ensure
  UpLeft  = tile->ClosestEnd(top);
  return UpLeft;
  */
  Point UpLeft = tile->ClosestEndYX(UpLeftPoint);
  return UpLeft;
}

char* Megatile::GetName(TileSide *tile)
{
  if (tile) return tile->nameAB;
  else      return NULL;
}

void Megatile::Area(double& xmin, double& xmax, double& ymin, double& ymax)
{
  if (ends->Count() == 0)
  {
    xmin = ymin = 0.;
    xmax = ymax = 1.;
    return;
  }

  int i;
  //int npoints = ends->Count();
  End *end;

  // initial values of xmin,ymin from ends
  i = 0;
  end = (*ends)[i];
  xmin = xmax = end->X;
  ymin = ymax = end->Y;

  // search through ends
  for (i=1; i<ends->Count(); i++) {
    end = (*ends)[i];
    if (end->X < xmin) xmin = end->X;
    if (end->X > xmax) xmax = end->X;
    if (end->Y < ymin) ymin = end->Y;
    if (end->Y > ymax) ymax = end->Y;
  }

  // search through down ends
  for (i=0; i<down_ends->Count(); i++) {
    end = (*down_ends)[i];
    if (end->X < xmin) xmin = end->X;
    if (end->X > xmax) xmax = end->X;
    if (end->Y < ymin) ymin = end->Y;
    if (end->Y > ymax) ymax = end->Y;
  }

  // search through contours
  for (int nt=0; nt<tds->Count(); nt++) {
    Queue<Side> *contour = (*tds)[nt]->tile->contour;
    for (int ns=0; ns<contour->Count(); ns++) {
      Point *end = &(*contour)[ns]->p1;
      if (end->X < xmin) xmin = end->X;
      if (end->X > xmax) xmax = end->X;
      if (end->Y < ymin) ymin = end->Y;
      if (end->Y > ymax) ymax = end->Y;
    }
  }

  // search through turned contours
  if (ShowTurnedContour) {
    for (int nt=0; nt<tds->Count(); nt++)
    {
      if ((*tds)[nt]->tile->makegroove == 0) continue;
      // to get the angle
      Groove* groove       = (*tds)[nt]->tile->groove;
      Queue<Side> *contour = (*tds)[nt]->tile->contour;
      for (int ns=0; ns<contour->Count(); ns++) {
        Side *side = (*contour)[ns];
        // rotate the side
        side->Rotate(O, +groove->angle);
        Point *end = &side->p1;
        if (end->X < xmin) xmin = end->X;
        if (end->X > xmax) xmax = end->X;
        if (end->Y < ymin) ymin = end->Y;
        if (end->Y > ymax) ymax = end->Y;
        // rotate back
        side->Rotate(O, -groove->angle);
      }
    }
  }
  /*
  // search through gates
  for (int ngate=0; ngate<gates->Count(); ngate++) {
    Side *gate = (*gates)[ngate];
    Point *p = &gate->p1;
    if (p->X < xmin) xmin = p->X;
    if (p->X > xmax) xmax = p->X;
    if (p->Y < ymin) ymin = p->Y;
    if (p->Y > ymax) ymax = p->Y;
    p = &gate->p2;
    if (p->X < xmin) xmin = p->X;
    if (p->X > xmax) xmax = p->X;
    if (p->Y < ymin) ymin = p->Y;
    if (p->Y > ymax) ymax = p->Y;
  }
  */
}

/*
int Megatile::TilesGetCount() {return tds->Count();}

int Megatile::SidesGetCount(int ntile)
{
  if (ntile < 0 || ntile >= tds->Count()) return NULL;
  return (*tds)[ntile]->tile->sides->Count();
}

Point Megatile::UpLeftPoint(int ntile)
{
  if (ntile < 0 || ntile >= tds->Count()) return NULL;
  Point UpLeftCorner(SheetX2,SheetY1);
  TileData *td = (*tds)[ntile];
  Point UpLeft = td->tile->ClosestEnd(UpLeftCorner);
  // this point can be at tile's down side,
  // so point stright under the tile is
  Point top(UpLeft.X, UpLeftCorner.Y);
  // and to be ensure
  UpLeft  = td->tile->ClosestEnd(top);
  return UpLeft;
}

char* Megatile::GetName(int ntile)
{
  if (ntile < 0 || ntile >= tds->Count()) return NULL;
  return (*tds)[ntile]->tile->nameAB;
}

Queue<Side>* Megatile::GetContour(int ntile)
{
  if (ntile < 0 || ntile >= tds->Count()) return NULL;
  return (*tds)[ntile]->tile->contour;
}

Queue<Side>* Megatile::GetTurnedContour(int ntile)
{
  if (ntile < 0 || ntile >= tds->Count()) return NULL;
  return (*tds)[ntile]->tile->turned_contour;
}

Groove* Megatile::GetGroove(int ntile)
{
  if (ntile < 0 || ntile >= tds->Count()) return NULL;
  return (*tds)[ntile]->tile->groove;
}

Queue<Side>* Megatile::GetSides(int ntile)
{
  if (ntile < 0 || ntile >= tds->Count()) return NULL;
  return (*tds)[ntile]->tile->sides;
}
*/

