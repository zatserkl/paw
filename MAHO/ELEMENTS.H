// Andriy_Zatserklyaniy@fnal.gov

#if !defined(__STDLIB_H)
#include <stdlib.h>
#endif

/*
#if !defined(__STDIO_H)
#include <stdio.h>
#endif
*/

#if !defined(__MATH_H)
#include <math.h>
#endif

#if !defined(__VALUES_H)
#include <values.h>
#endif

#if !defined(__IOSTREAM_H)
#include <iostream.h>
#endif

#if !defined(__IOMANIP_H)
#include <iomanip.h>
#endif

//// Borland macros from stdlib.h
//#define max(a,b)    (((a) > (b)) ? (a) : (b))
//#define min(a,b)    (((a) < (b)) ? (a) : (b))

////template <class T> T TipeId(T*) {return T;}

// General structure Frame

struct Frame {
  int Inv;
  double theta;
  double X0,Y0;
  Frame() {X0=Y0=theta=0.; Inv=0;}
  Frame(double x0, double y0, double thdeg, int inv)
  {
    X0=x0; Y0=y0; theta=thdeg; Inv=inv;
  }
};

// General-purpose function sign

inline int    sign(int    x) {return (x<0)? -1 : 1 ;}
inline double sign(double x) {return (x<0)? -1.: 1.;}
inline double SignAdd(double x, double dx) {return (x<0.)? x-dx: x+dx;}
inline double SignSub(double x, double dx) {return (x<0.)? x+dx: x-dx;}

struct Point;
double SQR(double a);
double R2(const Point& p1, const Point& p2);
double Distance(const Point& p1, const Point& p2);

struct Segment;
int RightFrame(const Segment& s1, const Segment& s2);

// struct Point

struct Point {
private:
  static short PointAnchor;
public:
  static void* Key() {return (void*) &PointAnchor;}

//private:
  static int count;
//public:
  double X,Y;
  Point(const Point& p) {X=p.X; Y=p.Y; count++;}
  Point() {count++;}
  Point(double xy) {X=xy; Y=xy; count++;}
  Point(double x, double y) {X=x; Y=y; count++;}
  ~Point() {count--;}
  Point& operator =(const Point& p) {X=p.X; Y=p.Y; return *this;}
  int operator==(const Point& to) {return (X == to.X && Y == to.Y)? 1: 0;}
  Point operator -() const {return Point(-X,-Y);}
  virtual void* GetKey() {return (void*) &PointAnchor;}
  void SetPoint(double x, double y) {X=x; Y=y;}
  virtual void Show() const {
    cout << "(";
    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(3);
    cout.width(12);
    cout << X << ",";
    cout.setf(ios::fixed, ios::floatfield);
    cout.precision(3);
    cout.width(12);
    cout << Y << ")";
    //printf("(%8.2g,%8.2g)", X,Y);
  }
  virtual size_t size_this() {return sizeof(*this);}
  virtual void Move(double dx, double dy, double dtheta, int inv) {
    if (inv) Y = -Y;
    //double thrad = dtheta*acos(-1.)/180.;
    double thrad = dtheta*M_PI/180.;
    // turn around the coordinate origin
    double X1 = X*cos(thrad) - Y*sin(thrad);
    double Y1 = X*sin(thrad) + Y*cos(thrad);
    // shift
    X = X1 + dx;
    Y = Y1 + dy;
  }
  virtual void Move(const Point&  r) {
    // Note: not a (const Point& ) to prevent side effects like p.Move(p)
    X += r.X;
    Y += r.Y;
  }
  void InvX() {X = -X;}
  void InvY() {Y = -Y;}
  //void Rotate(const Point& O, double theta);
  void Rotate(const Point& O, double theta) {
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
  void TurnBack(double A, double B) {
    // turn back with respect of y-axis
    X = A+B - X;
  }
};

struct Segment {
//private:
  static int count;
//public:
  Point p1,p2;

  Segment() {count++;}
  Segment(const Point& point1, const Point& point2) {
    p1 = point1; p2 = point2; count++;
  }
  Segment(double x1, double y1, double x2, double y2) {
    //p1.X = x1; p1.Y = y1; p2.X = x2; p2.Y = y2;
    p1=Point(x1,y1); p2=Point(x2,y2);
    count++;
  }
  Segment(const Segment& s) {p1=s.p1; p2=s.p2; count++;}
  ~Segment() {count--;}
  Segment& operator =(const Segment& s) {p1=s.p1; p2=s.p2; return *this;}
  Point& operator [] (int n) {if (n == 0) return p1; else return p2;}
  virtual void Move(double dx, double dy, double dtheta, int inv)
  {
    p1.Move(dx,dy,dtheta,inv);
    p2.Move(dx,dy,dtheta,inv);
  }
  virtual void Move(const Point& r)
  {
    p1.Move(r);
    p2.Move(r);
  }
  void Rotate(const Point& O, double theta) {
    p1.Rotate(O,theta);
    p2.Rotate(O,theta);
  }
  void TurnBack(double A, double B) {
    p1.TurnBack(A,B);
    p2.TurnBack(A,B);
  }
  virtual void Show() const {
    cout << "("; p1.Show();
    cout << ", ";
    p2.Show(); cout << ")";
  }
  virtual void InvX() {p1.InvX(); p2.InvX();}
  virtual void InvY() {p1.InvY(); p2.InvY();}
  virtual void Invert() {Point buf=p1; p1=p2; p2=buf;}
  Point MinX() {return (p1.X < p2.X)? p1: p2;}
  Point MaxX() {return (p1.X > p2.X)? p1: p2;}
  Point MinY() {return (p1.Y < p2.Y)? p1: p2;}
  Point MaxY() {return (p1.Y > p2.Y)? p1: p2;}
  Point MinX(double X0) {return (fabs(p1.X-X0) < fabs(p2.X-X0))? p1: p2;}
  Point MaxX(double X0) {return (fabs(p1.X-X0) > fabs(p2.X-X0))? p1: p2;}
  Point MinY(double Y0) {return (fabs(p1.Y-Y0) < fabs(p2.Y-Y0))? p1: p2;}
  Point MaxY(double Y0) {return (fabs(p1.Y-Y0) > fabs(p2.Y-Y0))? p1: p2;}
  int Overlap(const Segment* segment) const;
  double Length() const {return sqrt(SQR(p2.X-p1.X) + SQR(p2.Y-p1.Y));}
  double Theta() const
  {
    double dx = p2.X - p1.X;
    double dy = p2.Y - p1.Y;
    return (fabs(dx)+fabs(dy) > MINFLOAT)? atan2(dy,dx): 0.;
  }
  Point AtPar(double t) const;
  /*
  Point AtPar(double t) const
  {
    double dx = p2.X - p1.X;
    double dy = p2.Y - p1.Y;
    double theta = (fabs(dx)+fabs(dy) > MINFLOAT)? atan2(dy,dx): 0.;
    return Point(p1.X + cos(theta)*t, p1.Y + sin(theta)*t);
  }
  */
  Point ClosestPoint(const Point& topoint) const
  {
    double dx = p2.X - p1.X;
    double dy = p2.Y - p1.Y;
    double theta = (fabs(dx)+fabs(dy) > MINFLOAT)? atan2(dy,dx): 0.;

    double cosx = cos(theta);
    double cosy = sin(theta);

    double tmax = (fabs(cosx) > fabs(cosy))? dx/cosx: dy/cosy;

    double t = cosx*(topoint.X - p1.X) + cosy*(topoint.Y - p1.Y);
    if (t < 0.)   t=0.;
    if (t > tmax) t=tmax;
    return Point(p1.X + cosx*t, p1.Y + cosy*t);
  }
  Point ClosestEnd(const Point& topoint) const
  {
    double r2p1 = R2(p1,topoint);
    double r2p2 = R2(p2,topoint);
    //double r2p1 = Distance(p1,topoint);
    //double r2p2 = Distance(p2,topoint);
    return (r2p1 < r2p2)? p1: p2;
  }
  Point ClosestEnd(const Point& topoint, double& r2) const
  {
    r2 = R2(p1,topoint);
    double r2p2 = R2(p2,topoint);
    if (r2 < r2p2) return p1;
    else {
      r2 = r2p2;
      return p2;
    }
  }
};

// General-purpose routines

inline int min(int a, int b) {return (a<b)? a: b;}
inline int max(int a, int b) {return (a>b)? a: b;}
inline double min(double a, double b) {return (a<b)? a: b;}
inline double max(double a, double b) {return (a>b)? a: b;}

inline double SQR(double a) {return a*a;}

inline double Distance(const Point& p1, const Point& p2) {
  return sqrt(SQR(p1.X - p2.X) + SQR(p1.Y - p2.Y));
}

inline double R(const Point *p1, const Point *p2) {
  return sqrt(SQR(p1->X - p2->X) + SQR(p1->Y - p2->Y));
}

inline double R2(const Point *p1, const Point *p2) {
  return SQR(p1->X - p2->X) + SQR(p1->Y - p2->Y);
}

inline double R2(const Point& p1, const Point& p2) {
  return SQR(p1.X - p2.X) + SQR(p1.Y - p2.Y);
}

double Theta(const Point& p1, const Point& p2);
/*
inline double Theta(const Point& p1, const Point& p2) {
  double dx = p2.X - p1.X;
  double dy = p2.Y - p1.Y;
  return (fabs(dx)+fabs(dy) > MINFLOAT)? atan2(dy,dx): 0.;
}
*/

       double VertexAngle(const Point& vertex, const Point& p1, const Point& p2);
inline double VertexAngle(const Point& vertex, const Segment& s)
{
  return VertexAngle(vertex, s.p1, s.p2);
}
       double VertexAngleExclude(const Point& vertex, const Point& p1, const Point& p2);
inline double VertexAngleExclude(const Point& vertex, const Segment& s)
{
  return VertexAngleExclude(vertex, s.p1, s.p2);
}

Point CrossPoint(const Segment& s1, const Segment& s2);

Point InscribedInside
      (const Point& prev, const Point& p, const Point& next, double r);

Point InscribedOutside
      (const Point& prev, const Point& p, const Point& next, double r);

// General-purpose struct LinePar

struct LinePar
{
  double x0,y0;
  double tmax;
  double cosx,cosy;
  double theta;
  // distance from line to (0,0) like as normalized line eq.
  double pnorm;

  LinePar(const Segment *segment);
  LinePar(const Segment& segment);
  LinePar(const Point& p1, const Point& p2);
  LinePar(const Point *p1, const Point *p2);
  double t(const Point* p) {
    // parameter of my line point p
    if (fabs(cosx) > fabs(cosy)) return (p->X - x0)/cosx;
    else                         return (p->Y - y0)/cosy;
  }
};

struct tSegment {
  Point p;
  double tmax;
  double cosx,cosy;

  tSegment() {tmax=0;}
  tSegment(Point point1, Point point2) {
    p = point1;
    LinePar par(&point1,&point2);
    tmax = par.tmax;
    cosx = par.cosx;
    cosy = par.cosy;
  }
  tSegment(double x1, double y1, double x2, double y2) {
    p.X = x1; p.Y = y1;
    Point p2(x2,y2);
    LinePar par(&p,&p2);
    tmax = par.tmax;
    cosx = par.cosx;
    cosy = par.cosy;
  }
  virtual void Show() const {
    cout << "("; p.Show();
    cout << ", ";
    Point(p.X+cosx*tmax, p.Y+cosy*tmax).Show(); cout << ")";
  }
  void SetSegment(Point point1, Point point2) {
    p=point1;
    LinePar par(&point1,&point2);
    tmax = par.tmax;
    cosx = par.cosx;
    cosy = par.cosy;
  }
  virtual void Move(double dx, double dy, double dtheta, int inv) {

    // Should be something with inversion!

    p.Move(dx,dy,dtheta,inv);

    // should be invert here

    //double thrad = dtheta*acos(-1.)/180.;
    double thrad = dtheta*M_PI/180.;
    double cosx1 = cosx*cos(thrad) - cosy*sin(thrad);
    double cosy1 = cosx*sin(thrad) - cosy*cos(thrad);
    cosx = cosx1;
    cosy = cosy1;
  }
  virtual void Move(const Point& r) {
    // shift only
    p.Move(r);
  }
  virtual void Invert() {
    p.X = p.X + cosx*tmax;
    p.Y = p.Y + cosy*tmax;
    cosx = -cosx;
    cosy = -cosy;
  }
  int Overlap(const tSegment* segment) {return 1;}
};

//
// Abstract struct Element: base for independent cut elements:
// SHole, AHole and Side
//

struct Element {
  virtual void* GetKey() = 0;
  virtual void Show() const = 0;
  virtual void Move(double dx, double dy, double dtheta, int inv) = 0;
  virtual void Move(const Point& r) = 0;
  virtual void Rotate(const Point& O, double theta) = 0;
  virtual void TurnBack(double A, double B) = 0;
};

// struct Hole

struct Hole : public Element, public Point {
private:
  static short HoleAnchor;
public:
  static void* Key() {return (void*) &HoleAnchor;}

  double D;
  Hole(const Hole& hole) : Point(hole.X,hole.Y) {D=hole.D;}
  Hole() : Point() {D=0.;}
  Hole(double d) : Point() {D=d;}
  Hole(double x, double y, double d) : Point(x,y) {D=d;}
  Hole(const Point& p, double d) : Point(p) {D=d;}
  void Move(double dx, double dy, double dtheta, int inv) {
    Point::Move(dx,dy,dtheta,inv);
  }
  void Move(const Point& r) {
    Point::Move(r);
  }
  void Rotate(const Point& O, double theta) {Point::Rotate(O,theta);}
  void TurnBack(double A, double B) {Point::TurnBack(A,B);}
  virtual void* GetKey() {return (void*) &HoleAnchor;}
  void Show() const {Point::Show(); cout << " D=" << D;}
  //void Graph();
  size_t size_this() {return sizeof(*this);}
};

struct tHole : public Point {
  double t; // hole parameter with respect to Segment
  double D;
  tHole() : Point() {D=0.; t=0.;}
  tHole(double d) : Point() {D=d; t=0.;}
  tHole(const Point& p, double d, double t0) : Point(p.X,p.Y) {D=d; t=t0;}
  tHole(double x, double y, double d) : Point(x,y) {D=d; t=0.;}
  tHole(double x, double y, double d, double t0) : Point(x,y) {D=d; t=t0;}
  void Rotate(const Point& O, double theta) {Point::Rotate(O,theta);}
  void TurnBack(double A, double B) {Point::TurnBack(A,B);}
  void Show() const {cout << "(" << X << "," << Y << ") D=" << D << " t=" << t;}
  size_t size_this() {return sizeof(*this);}
};

// struct AHole: an Element, hole that stand Alone

struct AHole : public Hole {
private:
  static short AHoleAnchor;
public:
  static void* Key() {return (void*) &AHoleAnchor;}

  AHole(const AHole& ahole) : Hole(ahole) {}
  AHole() : Hole() {}
  AHole(double d) : Hole(d) {}
  AHole(double x, double y, double d) : Hole(x,y,d) {}
  virtual void* GetKey() {return (void*) &AHoleAnchor;}
  void Show() const {cout << "(" << X << "," << Y << ") D=" << D;}
  size_t size_this() {return sizeof(*this);}
};

// struct SHole: an Element, semihole that connected with a side

struct SHole : public Hole {
private:
  static short SHoleAnchor;
public:
  static void* Key() {return (void*) &SHoleAnchor;}

  SHole(const SHole& shole0) : Hole(shole0) {}
  SHole() : Hole() {}
  SHole(double d) : Hole(d) {}
  SHole(double x, double y, double d) : Hole(x,y,d) {}
  virtual void* GetKey() {return (void*) &SHoleAnchor;}
  void Show() const {cout << "(" << X << "," << Y << ") D=" << D;}
  size_t size_this() {return sizeof(*this);}
};

//
// General enum WorkStatus: will be used as a status by Side and TileData
//

enum WorkStatus {todo, done, extra, lift};

// struct Side: an Element that may contain semisholes

template <class T> class Queue; // used by class Side

struct Side : public Segment, public Element {
private:
  static short SideAnchor;
public:
  static void* Key() {return (void*) &SideAnchor;}

//private:
  static int count;
//public:
  // tile side that contains semisholes
  Queue<SHole>* sholes;
  WorkStatus status;

  Side() : Segment() {sholes=new Queue<SHole>; status=todo; count++;}
  Side(const Side& side);
  Side(const Point& p1, const Point& p2) : Segment(p1,p2) {
    sholes = new Queue<SHole>; status=todo; count++;
  }
  Side(double x1, double y1, double x2, double y2) :
      Segment(Point(x1,y1),Point(x2,y2))
  {
    sholes=new Queue<SHole>; status=todo; count++;
  }
  Side(const Point& p1, const Point& p2, WorkStatus thestatus):
      Segment(p1,p2)
  {
    sholes = new Queue<SHole>; status=thestatus; count++;
  }
  Side(double x1, double y1, double x2, double y2, WorkStatus thestatus):
      Segment(Point(x1,y1),Point(x2,y2))
  {
    sholes=new Queue<SHole>; status=thestatus; count++;
  }
  ~Side() {delete sholes; sholes=NULL; count--;}
  virtual void* GetKey() {return (void*) &SideAnchor;}
  Side& operator = (const Side& side);
  void Show() const {
    //p1.Show(); cout << " --> "; p2.Show(); cout << " " << status;
    p1.Show(); cout << " --> "; p2.Show(); cout << " " << status;
  }
  //virtual void Graph();
  void Clear() {delete sholes; sholes=NULL; status = todo;}
  //Side& Side::operator << (SHole *shole) {*sholes << shole; return *this;}
  Side& Side::operator << (SHole *shole);
  void AddSHole(SHole* NewSHole) {*sholes << NewSHole;}
  void RemoveSHole(int n) {
    // Routine can produce a side effects if n used as a cycle variable!
    if (sholes->Count()) sholes->Remove(n);
    else cout << "Side::RemoveSHole: Attempt of removing from empty list!\n";
  }
  void RemoveSHole(SHole *shole) {
    if(sholes->Count()) sholes->Remove(shole);
    else cout << "Side::RemoveSHole: Attempt of removing from empty list!\n";
  }
  void StripSHoles(Side *from);
  void  CopySHoles(Side *from);
  int Overlap(Side* alien, Side* overlap);
  virtual void Move(double dx, double dy, double dtheta, int inv);
  virtual void Move(const Point& r);
  void Rotate(const Point& O, double theta);
  void TurnBack(double A, double B);
  void InvX();
  void InvY();
  virtual void Invert() {Segment::Invert(); sholes->Invert();}
  double tshole(SHole *shole);
};

// struct Groove

struct Groove : public Element, public Queue<Point> {
private:
  static short GrooveAnchor;
public:
  static void* Key() {return (void*) &GrooveAnchor;}

  double angle;
  WorkStatus cylstatus;
  WorkStatus balstatus;

  Groove() : Queue<Point>() {angle=0.; cylstatus=balstatus=todo;};
  Groove(const Groove& groove);
  Groove& operator =(const Groove& groove);
  void* GetKey() {return (void*) &GrooveAnchor;}
  void Show() const;
  void Move(double dx, double dy, double dtheta, int inv);
  void Move(const Point& r);
  void Rotate(const Point& O, double theta);
  void TurnBack(double A, double B);
};
