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

void Megatile::GroovingEnds()
{
  // Corresponds to the current Arrange method

  Queue<TileData> *strit, *pos10, *neg10;
  strit  = new Queue<TileData>;
  pos10  = new Queue<TileData>;
  neg10  = new Queue<TileData>;

  const double limit = 1.;
  double angle, dangle;
  const double todeg = 180./M_PI;

  TileData *td;
  int nt;
  for (nt=0; nt<tds->Count(); nt++)
  {
    td = (*tds)[nt];
    //if (td->tile->makegroove == 0) continue;

    Point p1 = *(*td->tile->groove)[1];
    Point p2 = *(*td->tile->groove)[2];
    angle = todeg*Theta(p1,p2);

    //cout << "Tile " << td->tile->nameAB << ", angle = " << angle << endl;

    if (angle < 0.) angle += 180.;

    angle = 90.-angle;

    if (fabs(angle) < 5.) {
      if (fabs(angle) < 1.e-3) angle = 0.;
      if (TrueAngle) td->tile->groove->angle = angle;
      else           td->tile->groove->angle = 0.;
      *strit << td;
      dangle = fabs(angle);
    }
    else {
      if (angle > 0.) {
        if (TrueAngle) td->tile->groove->angle = angle;
        else           td->tile->groove->angle = +10.;
        *pos10 << td;
        dangle = fabs( angle-10.);
      }
      else {
        if (TrueAngle) td->tile->groove->angle = angle;
        else           td->tile->groove->angle = -10.;
        *neg10 << td;
        dangle = fabs(-angle-10.);
      }
    }
    if (dangle > limit) {
      cout.precision(2);
      cout << "Groove angle for tile " << td->tile->nameAB
           << " exeeds limit. Angle = " << angle << " degrees" << endl;
      WAIT();
      WAIT("Please, confirm");
    }
  }

  /*
  cout << "strit: " << strit->Count() << endl;
  cout << "pos10: " << pos10->Count() << endl;
  cout << "neg10: " << neg10->Count() << endl;
  */

  // cyl. groove
  AddGrooves(strit, cylcut);
  AddGrooves(pos10, cylcut);
  AddGrooves(neg10, cylcut);

  // ball groove
  AddGrooves(neg10, balcut);
  AddGrooves(pos10, balcut);
  AddGrooves(strit, balcut);

  while (strit->Count()) strit->Remove(0);
  delete strit;
  while (pos10->Count()) pos10->Remove(0);
  delete pos10;
  while (neg10->Count()) neg10->Remove(0);
  delete neg10;
}

void Megatile::AddGrooves(Queue<TileData> *list, EndCutter ecut)
{
  if (list->Count() == 0) return;

  const double todeg = 180./M_PI;

  for (int i=0; i<list->Count(); i++)
  {
    TileSide *tile = (*list)[i]->tile;

    EndStatus estatus;

    estatus = tocut;
    switch (continue_level) {
      case 0: break;
      case 3: estatus = tofly; break;
      case 1: if (ecut==cylcut)
              {
                if (strcmp(tile->name,continue_tile)==0) {
                  estatus = tocut;
                  continue_level   = 0;
                  continue_tile[0] = '\0';
                }
                else estatus = tofly;
              }
              break;
      case 2: if (ecut==cylcut) {
                estatus = tofly;
              }
              else if (strcmp(tile->name,continue_tile)==0) {
                estatus = tocut;
                continue_level   = 0;
                continue_tile[0] = '\0';
              }
              else estatus = tofly;
              break;
    }

    if (tile->makegroove == 0)
    {
      estatus = tofly;
      tile->groove->cylstatus = lift;
      tile->groove->balstatus = lift;
      continue;
    }

    if (estatus == tofly) {
      // set status field in tile->groove from todo to done
      switch (ecut) {
        case cylcut: tile->groove->cylstatus = done; break;
        case balcut: tile->groove->balstatus = done; break;
      }
    }

    Groove *groove = tile->groove;

    EndGroove *egroove = new EndGroove;
    // the center of rotation and angle
    egroove->X = O.X;
    egroove->Y = O.Y;
    egroove->Dcutter = groove->angle;
    egroove->cutype = ecut;
    egroove->status = estatus;
    // pointer to parent tile
    egroove->tile   = tile;

    //cout << "Tile " << (*list)[i]->tile->nameAB
    //     << " status = " << egroove->status << endl;

    // rotate
    groove->Rotate(O,+groove->angle);

    /*
    //DEBUG
    //cout << "contour:" << endl;
    //tile->contour->Show();
    tile->contour->Rotate(O,+groove->angle);
    cout << "sloped contour:" << endl;
    tile->contour->Show();
    tile->contour->Rotate(O,-groove->angle);
    */

    Point *p;
    Hole  *h;

    EndStatus status;

    // Field Dcutter:
    // for initial point : Tile::dh0groove (0.5 mm)
    // for first point   : vertical   +R38 or -R38
    // for common cut    : 0.
    // for arc           : horizontal +R38 or -R38

    int np;
    // initial point
    np = 0;
    h = (Hole*) (*groove)[np];
    status = tofly;
    *egroove->P << new End(Point(h->X,h->Y), Tile::P0height,ecut,status,tile);

    status = tocut;

    // vertical cut to point P2
    np = 1;
    np = 2;
    p = (*groove)[2];
    double theta = todeg * Theta(*(*groove)[1], *p);
    if (theta > 0.)
    {
      // angle +90 corresponds to counterclockwise vertical move
      *egroove->P << new End(*p, -Tile::R38,ecut,status,tile);
    }
    else {
      *egroove->P << new End(*p, +Tile::R38,ecut,status,tile);
    }

    for(np=3; np<groove->Count(); np++) {
      Point *next;
      if ((*groove)[np]->GetKey() == Point::Key()) {
        next = (*groove)[np];
        //line(Xscr(p->X),Yscr(p->Y), Xscr(next->X),Yscr(next->Y));
        status = tocut;
        *egroove->P << new End(*next, 0.,ecut,status,tile);
      }
      else {
        Hole *r = (Hole*) (*groove)[np];
        next    = (*groove)[++np];

        int ccw = RightFrame(Segment(*r, *p), Segment(*r, *next));
        //if (ccw) arc(Xscr(r->X),Yscr(r->Y),theta1,theta2,R);
        //else     arc(Xscr(r->X),Yscr(r->Y),theta2,theta1,R);
        if (ccw) *egroove->P << new End(*next, -r->D,ecut,status,tile);
        else     *egroove->P << new End(*next, +r->D,ecut,status,tile);
      }
      p = next;
    }

    // turn back
    groove->Rotate(O,-groove->angle);

    *groove_ends << egroove;
  }
}
