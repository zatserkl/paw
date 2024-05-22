// Andriy_Zatserklyaniy@fnal.gov

#if !defined(__STDLIB_H)
#include <stdlib.h>
#endif

#if !defined(__MATH_H)
#include <math.h>
#endif

#if !defined(__VALUES_H)
#include <values.h>
#endif

#if !defined(__IOSTREAM_H)
#include <iostream.h>
#endif

#include "queue.h"
#include "elements.h"
#include "wait.h"

// General-purpose routines

double Theta(const Point& p1, const Point& p2) {
  double dx = p2.X - p1.X;
  double dy = p2.Y - p1.Y;
  return (fabs(dx)+fabs(dy) > MINFLOAT)? atan2(dy,dx): 0.;
}

double VertexAngle(const Point& vertex, const Point& p1, const Point& p2)
{
  const double eps = 1.e-6;

  Segment p1p2(p1,p2);
  if (p1p2.Length() < eps) return 0.;

  if (Distance(vertex, p1p2.ClosestPoint(vertex)) < eps) return M_PI;

  LinePar v1(vertex,p1);
  LinePar v2(vertex,p2);
  // scalar product
  double scal = v1.cosx*v2.cosx + v1.cosy*v2.cosy;
  if (scal < -1.) scal=-1.;
  if (scal > +1.) scal=+1.;
  return acos(scal);
}

double VertexAngleExclude(const Point& vertex, const Point& p1, const Point& p2)
{
  // exclude point in (p1,p2)

  const double eps = 1.e-6;

  Segment p1p2(p1,p2);
  if (p1p2.Length() < eps) return 0.;

  if (Distance(vertex, p1p2.ClosestPoint(vertex)) < eps) return 0.;
  /*
  Point p = p1p2.ClosestPoint(vertex);
  double dx = Distance(vertex, p);
  if (dx < eps) return 0.;
  */

  LinePar v1(vertex,p1);
  LinePar v2(vertex,p2);
  // scalar product
  double scal = v1.cosx*v2.cosx + v1.cosy*v2.cosy;
  if (scal < -1.) scal=-1.;
  if (scal > +1.) scal=+1.;
  return acos(scal);
}

/*
double VertexAngle(const Point& vertex, const Segment& s)
{
  return VertexAngle(vertex, s.p1, s.p2);
}
*/

Point CrossPoint(const Segment& s1, const Segment& s2)
{
  // simplified version for non-parallel segments
  double eps2=1.e-12;

  double dx,dy,theta;
  // directrix cosines
  double cx1,cy1,cx2,cy2;

  dx = s1.p2.X - s1.p1.X;
  dy = s1.p2.Y - s1.p1.Y;
  theta = (fabs(dx)+fabs(dy) > MINFLOAT)? atan2(dy,dx): 0.;
  cx1 = cos(theta);
  cy1 = sin(theta);

  dx = s2.p2.X - s2.p1.X;
  dy = s2.p2.Y - s2.p1.Y;
  theta = (fabs(dx)+fabs(dy) > MINFLOAT)? atan2(dy,dx): 0.;
  cx2 = cos(theta);
  cy2 = sin(theta);

  // scalar product
  double A = cx1*cx2 + cy1*cy2;
  if (fabs(fabs(A)-1.) < eps2) {
    cerr << "CrossPoint: Parallel lines (or zero length):" << endl;
    cerr << "s1: "; s1.Show(); cerr << endl;
    cerr << "s2: "; s2.Show(); cerr << endl;
    STOP();
  }

  dx = s1.p1.X - s2.p1.X;
  dy = s1.p1.Y - s2.p1.Y;
  // Parameter of the first line that corresponds to cross point
  double u = (dx*(A*cx2-cx1) + dy*(A*cy2-cy1)) / (1-A*A);
  // simular for the second line:
  // double v = (-dx*(A*cx1-cx2) - dy*(A*cy1-cy2)) / (1-A*A);

  return Point(s1.p1.X+cx1*u, s1.p1.Y+cy1*u);
}

Point InscribedInside
      (const Point& prev, const Point& p, const Point& next, double r)
{
  // Inscribes the circle into the corner

  const double eps = 1.e-7;
  // bisectrix at point p
  //     atan2(p - prev)
  double theta1 = Theta(prev,p);
  //     atan2(p - next)
  double theta2 = Theta(next,p);
  //
  // bisectrix parametric equation
  // x = x0 + cosx*t
  // y = y0 + cosy*t
  //
  double theta = (theta2+theta1)/2.;
  double cosx = cos(theta);
  double cosy = sin(theta);

  // bisectrix parameter of circle center,
  // counting from the corner p inside the angle
  if (fabs(theta1-theta2) < eps) {
    cerr << "InscribedInside: Angle is closed to 180 degrees:" << endl;
    cerr << "prev = "; prev.Show(); cerr << endl;
    cerr << "p    = "; p.Show();    cerr << endl;
    cerr << "next = "; next.Show(); cerr << endl;
    STOP();
  }
  double t = r / sin((theta1-theta2)/2.);

  // coordinates of circle center
  return Point(p.X + cosx*t, p.Y + cosy*t);
}

Point InscribedOutside
      (const Point& prev, const Point& p, const Point& next, double r)
{
  // Inscribes the circle into the supplementary corner to given one

  const double eps = 1.e-7;
  // bisectrix at point p
  //     atan2(p - prev)
  double theta1 = Theta(prev,p);
  //     atan2(p - next)
  double theta2 = Theta(next,p);
  //
  // bisectrix parametric equation
  // x = x0 + cosx*t
  // y = y0 + cosy*t
  //
  double theta = (theta2+theta1)/2.;
  double cosx = cos(theta);
  double cosy = sin(theta);

  // bisectrix parameter of circle center,
  // counting from the corner p outside the angle
  if (fabs(theta1-theta2) < eps) {
    cerr << "InscribedOutside: Angle is closed to 180 degrees:" << endl;
    cerr << "prev = "; prev.Show(); cerr << endl;
    cerr << "p    = "; p.Show();    cerr << endl;
    cerr << "next = "; next.Show(); cerr << endl;
    STOP();
  }
  double t = r / sin((theta2-theta1)/2.);

  // coordinates of circle center
  return Point(p.X + cosx*t, p.Y + cosy*t);
}

int RightFrame(const Segment& s1, const Segment& s2)
{
  const double eps2=1.e-12;

  LinePar ray1(s1);
  LinePar ray2(s2);
  if (fabs(ray1.cosx*ray2.cosx + ray1.cosy*ray2.cosy) > 1.-eps2) {
    cerr << "RightFrame: Segmants are close to be collinear!!" << endl;
    exit(1);
  }
  if (ray1.cosx*ray2.cosy - ray1.cosy*ray2.cosx > 0.) return 1;
  else                                                return 0;
}

// class Point

// declaration and definition of static variable
short Point::PointAnchor; // Point key
int Point::count = 0;

/*
void Point::Rotate(const Point& O, double theta) {
  // coordinates with respect to O
  X -= O.X;
  Y -= O.Y;
  //double thrad = theta*acos(-1.)/180.;
  double thrad = theta*M_PI/180.;
  // turn around the O
  double X1 = X*cos(thrad) - Y*sin(thrad);
  double Y1 = X*sin(thrad) + Y*cos(thrad);
  // accept and restore the coordinate origin
  X = X1 + O.X;
  Y = Y1 + O.Y;
}
*/

// class Hole

// declaration and definition of static variable
short Hole::HoleAnchor; // Hole key

// clas LinePar

LinePar::LinePar(const Segment *segment)
{
  x0 = segment->p1.X;
  y0 = segment->p1.Y;
  double dx = segment->p2.X - x0;
  double dy = segment->p2.Y - y0;
  theta = (fabs(dx)+fabs(dy) > MINFLOAT)? atan2(dy,dx): 0.;

  cosx = cos(theta);
  cosy = sin(theta);

  pnorm  = y0*cosx - x0*cosy;

  tmax =
  (fabs(cosx)>fabs(cosy))? (segment->p2.X-x0)/cosx: (segment->p2.Y-y0)/cosy;
}

LinePar::LinePar(const Segment& segment)
{
  x0 = segment.p1.X;
  y0 = segment.p1.Y;
  double dx = segment.p2.X - x0;
  double dy = segment.p2.Y - y0;
  theta = (fabs(dx)+fabs(dy) > MINFLOAT)? atan2(dy,dx): 0.;

  cosx = cos(theta);
  cosy = sin(theta);

  pnorm  = y0*cosx - x0*cosy;

  tmax =
  (fabs(cosx)>fabs(cosy))? (segment.p2.X-x0)/cosx: (segment.p2.Y-y0)/cosy;
}

LinePar::LinePar(const Point& p1, const Point& p2)
{
  x0 = p1.X;
  y0 = p1.Y;
  double dx = p2.X - x0;
  double dy = p2.Y - y0;
  theta = (fabs(dx)+fabs(dy) > MINFLOAT)? atan2(dy,dx): 0.;

  cosx = cos(theta);
  cosy = sin(theta);

  pnorm  = y0*cosx - x0*cosy;

  tmax = (fabs(cosx) > fabs(cosy))? (p2.X - x0)/cosx: (p2.Y - y0)/cosy;
}

LinePar::LinePar(const Point *p1, const Point *p2)
{
  x0 = p1->X;
  y0 = p1->Y;
  double dx = p2->X - x0;
  double dy = p2->Y - y0;
  theta = (fabs(dx)+fabs(dy) > MINFLOAT)? atan2(dy,dx): 0.;

  cosx = cos(theta);
  cosy = sin(theta);

  pnorm  = y0*cosx - x0*cosy;

  tmax = (fabs(cosx) > fabs(cosy))? (p2->X - x0)/cosx: (p2->Y - y0)/cosy;
}

// class Segment

// declaration and definition of static variable
int Segment::count = 0;

Point Segment::AtPar(double t) const
{
  double dx = p2.X - p1.X;
  double dy = p2.Y - p1.Y;
  double theta = (fabs(dx)+fabs(dy) > MINFLOAT)? atan2(dy,dx): 0.;
  return Point(p1.X + cos(theta)*t, p1.Y + sin(theta)*t);
}

int Segment::Overlap(const Segment* segment) const
{
  // Routine for test purpose mainly

  const double eps  = 1.E-7;
  const double eps2 = 1.E-12;
  // my pars.
  LinePar This(this);
  double cosx = This.cosx;
  double cosy = This.cosy;
  double tmax = This.tmax;
  // segment pars.
  LinePar That(segment);

  //if (fabs(cosx-That.cosx) > eps || fabs(cosy-That.cosy) > eps) return 0;
  if (fabs(fabs(cosx*That.cosx + cosy*That.cosy) - 1.) > eps2) return 0;

  // segments are colinear; test distance to (0,0)

  //if (fabs(This.pnorm-That.pnorm) > eps) return 0;
  // parameter that corresponds to segment first point
  double x0 = This.x0;
  double y0 = This.y0;
  double t = cosx*(That.x0 - x0) + cosy*(That.y0 - y0);
  // distance between the lines
  double dx = x0 + cosx*t - That.x0;
  double dy = y0 + cosy*t - That.y0;
  double r = sqrt(dx*dx + dy*dy);
  if (r > eps) return 0;

  // segments lie at the same line

  double t1,t2;
  if (fabs(cosx) > fabs(cosy)) {
    t1 = (segment->p1.X - p1.X)/cosx;
    t2 = (segment->p2.X - p1.X)/cosx;
  }
  else {
    t1 = (segment->p1.Y - p1.Y)/cosy;
    t2 = (segment->p2.Y - p1.Y)/cosy;
  }
  if (t1 > t2) {
    double buf = t1;
    t1 = t2;
    t2 = buf;
  }
  if (t2 < 0 || t1 > tmax) return 0;

  // segments are overlapped!

  if (t1 < -eps)     t1=0;
  if (t2 > tmax+eps) t2=tmax;

  // same part segment -- Do not will be return!!!!!!!!!!!!!!!!
  Segment same(p1.X+cosx*t1, p1.Y+cosy*t1, p1.X+cosx*t2, p1.Y+cosy*t2);
  return 1;
}

// class SHole

// declaration and definition of static variable
short SHole::SHoleAnchor; // SHole key

// class AHole

// declaration and definition of static variable
short AHole::AHoleAnchor; // AHole key

// class Side

// declaration and definition of static variable
short Side::SideAnchor; // AHole key
int Side::count = 0;

Side::Side(const Side& side) : Segment(side.p1, side.p2) {
  status = side.status;
  sholes = new Queue<SHole>;
  //*sholes = *side.sholes; // uses Queue<T>::operator=
  ///*
  // Possible way:
  for (int i=0; i<side.sholes->Count(); i++) {
    *sholes << new SHole(*(*side.sholes)[i]);
    //Hole *shole = (*side.sholes)[i];
    //if (hole->GetKey() == Hole::Key()) *sholes << new Hole(*hole);
    //else {
    //  if   (hole->GetKey() == SHole::Key()) *sholes << new SHole(*hole);
    //  else (hole->GetKey() == AHole::Key()) *sholes << new AHole(*hole);
    //}
  }
  //*/
  count++;
}

Side& Side::operator << (SHole *shole)
{
  //*sholes << shole;
  sholes->Add(shole);
  return *this;
}

Side& Side::operator = (const Side& side) {
  if (this != &side) {
    p1 = side.p1;
    p2 = side.p1;
    status = side.status;

    //*sholes = *side.sholes;
    delete sholes;
    sholes = new Queue<SHole>;
    for (int i=0; i<side.sholes->Count(); i++)
    {
      *sholes << new SHole(*(*side.sholes)[i]);
    }
  }
  return *this;
}

void Side::InvX() {
  Segment::InvX();
  for(int i=0; i<sholes->Count(); i++) (*sholes)[i]->InvX();
}
void Side::InvY() {
  Segment::InvY();
  for(int i=0; i<sholes->Count(); i++) (*sholes)[i]->InvY();
}

void Side::Move(double dx, double dy, double dtheta, int inv)
{
  // ends
  Segment::Move(dx,dy,dtheta,inv);
  // sholes
  for (int i=0;i<sholes->Count();i++) (*sholes)[i]->Move(dx,dy,dtheta,inv);
}

void Side::Move(const Point& r)
{
  // ends
  Segment::Move(r);
  // sholes
  for (int i=0; i<sholes->Count(); i++) (*sholes)[i]->Move(r);
}

void Side::Rotate(const Point& O, double theta) {
  //p1.Rotate(O,theta);
  //p2.Rotate(O,theta);
  Segment::Rotate(O,theta);
  for (int i=0; i<sholes->Count(); i++) (*sholes)[i]->Rotate(O,theta);
}

void Side::TurnBack(double A, double B) {
  Segment::TurnBack(A,B);
  for (int i=0; i<sholes->Count(); i++) (*sholes)[i]->TurnBack(A,B);
}

void Side::StripSHoles(Side *from)
{
  if (from->sholes->Count() == 0) return;

  double eps = 1.E-5;
  // my pars.
  LinePar This(this);
  double cosx = This.cosx;
  double cosy = This.cosy;
  double tmax = This.tmax;
  double   x0 = This.x0;
  double   y0 = This.y0;

  for (int nshole=0; nshole<from->sholes->Count(); nshole++) {
    SHole *shole = (*from->sholes)[nshole];
    // parameter of base of shole in my frame
    double t = cosx*(shole->X - x0) + cosy*(shole->Y - y0);
    if (t > -eps && t < tmax + eps) {
      // have to add
      if (sholes->Count() > 0) {
        // put in right position among my sholes
        for (int i=0; i<sholes->Count(); i++) {
          SHole *curr = (*sholes)[i];
          // parameter of base of my current shole
          double tcurr = cosx*(curr->X - x0) + cosy*(curr->Y - y0);
          if (t < tcurr) {
            sholes->InsertBefore(i,shole);
            goto REMOVE;
          }
        }
      }
      // add at the end of my shole list
      *sholes << shole;

      REMOVE:
      from->RemoveSHole(nshole);
      // Note: next shole will have the same number 'nshole'!
      nshole--;
    }
  }
}

void Side::CopySHoles(Side *from)
{
  //WAIT("Side::CopySHoles");
  // briefly, using "dangerous" operator
  //*sholes = *from->sholes;

  ///*
  // clear my sholes
  delete sholes; sholes = new Queue<SHole>;

  int nsholes = from->sholes->Count();
  if (nsholes == 0) return;

  for (int nshole=0; nshole<nsholes; nshole++)
  {
    Hole *hole = (*from->sholes)[nshole];
    if (hole->GetKey() == SHole::Key()) {
      *sholes << new SHole(*( (SHole*) hole) );
      continue;
    }
    if (hole->GetKey() == AHole::Key()) {
      //*sholes << new AHole(*( (AHole*) hole) );
      //continue;
      STOP("Possibility don't realised yet!");
    }
  }
  //*/
  //cout << "Side::CopySHoles complete" << endl;
}

int Side::Overlap(Side *alien, Side *overlap)
{
  const double eps  = 1.E-7;
  const double eps2 = 1.E-12;
  // my pars.
  LinePar This(this);
  double cosx = This.cosx;
  double cosy = This.cosy;
  double tmax = This.tmax;
  double   x0 = This.x0;
  double   y0 = This.y0;
  // alien pars.
  LinePar That(alien);

  //if (fabs(cosx-That.cosx) > eps || fabs(cosy-That.cosy) > eps) return 0;
  //if (fabs(fabs(cosx*That.cosx + cosy*That.cosy) - 1.) > eps2) return 0;
  double ee = cosx*That.cosx + cosy*That.cosy;
  if (fabs(fabs(ee) - 1.) > eps2) return 0;

  // segments are colinear; test distance to (0,0)

  //if (fabs(This.pnorm-That.pnorm) > eps) return 0;
  // parameter that corresponds to segment first point
  double t = cosx*(That.x0 - x0) + cosy*(That.y0 - y0);
  // distance between the lines
  double dx = x0 + cosx*t - That.x0;
  double dy = y0 + cosy*t - That.y0;
  double r = sqrt(dx*dx + dy*dy);
  if (r > eps) return 0;

  // segments lie at the same line!

  // parameters of alien side ends in my side frame
  double t1,t2;

  if (fabs(cosx) > fabs(cosy)) {
    t1 = (alien->p1.X - x0)/cosx;
    t2 = (alien->p2.X - x0)/cosx;
  }
  else {
    t1 = (alien->p1.Y - y0)/cosy;
    t2 = (alien->p2.Y - y0)/cosy;
  }
  if (t1 > t2) {
    double buf = t1;
    t1 = t2;
    t2 = buf;
  }
  if (t2 < -eps || t1-tmax > eps) return 0;

  // segments are overlapped!

  if (t1 < eps)      t1=0;
  if (t2 > tmax-eps) t2=tmax;

  // same part segment
  overlap->p1.X = x0 + cosx*t1;
  overlap->p1.Y = y0 + cosy*t1;
  overlap->p2.X = x0 + cosx*t2;
  overlap->p2.Y = y0 + cosy*t2;

  if (t2-t1 < eps) return 0; // zero overlapping
  return 1;
}

double Side::tshole(SHole *shole)
{
  LinePar This(this);
  double cosx = This.cosx;
  double cosy = This.cosy;
  double   x0 = This.x0;
  double   y0 = This.y0;
  double t = cosx*(shole->X - x0) + cosy*(shole->Y - y0);
  return t;
}

// class Groove

// declaration and definition of static variable
short Groove::GrooveAnchor; // Groove key

Groove::Groove(const Groove& groove) : Queue<Point>()
//Groove::Groove(const Groove& groove) // possible way
{
  angle     = groove.angle;
  cylstatus = groove.cylstatus;
  balstatus = groove.balstatus;

  // add in list copies of 'queue' items
  for (Node<Point> *curr=groove.head; curr; curr=curr->Next) {
    // create a copy of Item
    //
    // It's more sutable a usage of typeid stament to create a real copy of Item!
    //
    Point* cp;
    if (curr->Item->GetKey() == Point::Key()) cp = new Point(*curr->Item);
    else {
      Hole *hole = (Hole*) curr->Item;
      cp = new Hole(*hole);
    }
    Add(cp);
  }
}

Groove& Groove::operator =(const Groove& groove)
{
  if (this != &groove) {
    angle     = groove.angle;
    cylstatus = groove.cylstatus;
    balstatus = groove.balstatus;

    Node<Point> *curr;
    // delete existing list
    while (head != NULL) {
      curr = head;
      head = curr->Next;
      // I have to check this Item in the other nodes in the list!
      delete curr;
      count--;
    }
    // add in list copies of 'queue' items
    for (curr=groove.head; curr; curr=curr->Next) {
      // create a copy of Item
      //
      // I need in typeid stament to create a real copy of Item!
      //
      Point* cp;
      if (curr->Item->GetKey() == Point::Key()) cp = new Point(*curr->Item);
      else {
        Hole *hole = (Hole*) curr->Item;
        cp = new Hole(*hole);
      }
      Add(cp);
    }
  }
  return *this;
}

void Groove::Move(double dx, double dy, double dtheta, int inv) {
  angle += dtheta;
  for (Node<Point> *curr=head; curr; curr=curr->Next)
    curr->Item->Move(dx,dy,dtheta,inv);
}
void Groove::Move(const Point& r) {
  for (Node<Point> *curr=head; curr; curr=curr->Next) curr->Item->Move(r);
}

void Groove::Rotate(const Point& O, double theta) {
  for (Node<Point> *curr=head; curr; curr=curr->Next)
    curr->Item->Rotate(O,theta);
}
void Groove::TurnBack(double A, double B) {
  for (Node<Point> *curr=head; curr; curr=curr->Next)
    curr->Item->TurnBack(A,B);
}

void Groove::Show() const {
  cout << "angle = " << angle << endl;
  int iP=0, iR=0;
  for (Node<Point> *curr=head; curr; curr=curr->Next) {
    if (curr->Item->GetKey() == Point::Key()) {
      //cout << "P" << ++iP << " = "; curr->Item->Show(); cout << endl;
      cout << "P" << iP++ << " = "; curr->Item->Show(); cout << endl;
    }
    else {
      cout << "R" << ++iR << " = "; curr->Item->Show(); cout << endl;
    }
  }
}
