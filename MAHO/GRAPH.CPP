// Andriy_Zatserklyaniy@fnal.gov

#if !defined(__STDIO_H)
#include <stdio.h>
#endif

#if !defined(__STRING_H)
#include <string.h>
#endif

#if !defined(__CTYPE_H)
#include <ctype.h>      /* for toupper */
#endif

#if !defined(__MATH_H)
#include <math.h>
#endif

#if !defined(__IOSTREAM_H)
#include <iostream.h>
#endif

#if !defined(__TIME_H)
#include <time.h>
#endif

#include <graphics.h>
#include <dos.h>        /* for delay() */
#include <conio.h>

#include "queue.h"
#include "elements.h"

#include "sidelist.h"
#include "tile.h"

#include "megatile.h"
#include "maho.h"

#include "graph.h"

#include "wait.h"

//
// static
//
int  GManager::graphics = 0;
long GManager::CuttingDelay  = 0;  /* in milliseconds */
long GManager::FlyingDelay   = 0;  /* in milliseconds */
long GManager::GroovingDelay = 0;  /* in milliseconds */

int GManager::Color           = -1;
int GManager::BkColor         = -1;
int GManager::TextColor       = -1;
int GManager::MessageColor    = -1;
int GManager::CutColor        = -1;
int GManager::FlyColor        = -1;
int GManager::ContourColor    = -1;
int GManager::SideColor       = -1;
int GManager::CylGrooveColor  = -1;
int GManager::BallGrooveColor = -1;
int GManager::GatesColor      = -1;
int GManager::SheetColor      = -1;
int GManager::WorkAreaColor   = -1;
int GManager::FinishColor     = -1;
int GManager::Rfinish         = -1;
int GManager::FinishSnow      =  0;

// globals
int YMARGIN = 0;
int XMARGIN = 0;
double minx,maxx,miny,maxy;
double Xscal,Yscal, scal;
double Xscal0,Yscal0, scal0;

void Colors() {
  int c0 = getcolor();
  setcolor(0); outtext(" 0");
  setcolor(1); outtext(" 1");
  setcolor(2); outtext(" 2");
  setcolor(3); outtext(" 3");
  setcolor(4); outtext(" 4");
  setcolor(5); outtext(" 5");
  setcolor(6); outtext(" 6");
  setcolor(7); outtext(" 7");
  setcolor(8); outtext(" 8");
  setcolor(9); outtext(" 9");
  setcolor(10); outtext(" 10");
  setcolor(11); outtext(" 11");
  setcolor(12); outtext(" 12");
  setcolor(13); outtext(" 13");
  setcolor(14); outtext(" 14");
  setcolor(15); outtext(" 15");
  setcolor(16); outtext(" 16");
  setcolor(c0);
}

// class Window
struct viewporttype     Window::VPORTCURR;
struct textsettingstype Window::TEXTSETCURR;
int                     Window::COLORCURR;

void Window::SaveCurrentWindow() {
  getviewsettings(&VPORTCURR);
  gettextsettings(&TEXTSETCURR);
  COLORCURR = getcolor();
  Xscal0 = Xscal;
  Yscal0 = Yscal;
   scal0 =  scal;
}

void Window::RestoreSavedWindow() {
  setviewport(VPORTCURR.left, VPORTCURR.top,
              VPORTCURR.right,VPORTCURR.bottom,VPORTCURR.clip);
  settextjustify(LEFT_TEXT, TOP_TEXT);
  setcolor(COLORCURR);
  Xscal = Xscal0;
  Yscal = Yscal0;
   scal =  scal0;
}

int Window::hText() {return textheight("W");}
int Window::wText() {return textwidth("W");}

Window::Window(int x1, int y1, int x2, int y2, int clip)
{
  X1=x1; Y1=y1; X2=x2; Y2=y2; Clip=clip;
  dX = X2-X1+1;
  dY = Y2-Y1+1;
}

Window::Window() {}

void Window::DefineWindow(int x1, int y1, int x2, int y2, int clip)
{
  X1=x1; Y1=y1; X2=x2; Y2=y2; Clip=clip;
  dX = X2-X1+1;
  dY = Y2-Y1+1;
}

void Window::Clear()
{
  clearviewport();
  int c0 = getcolor();
  setcolor(WHITE);
  int dx = X2-X1+1;
  int dy = Y2-Y1+1;
  line(1,    1,    dx-2, 1);
  line(dx-2, 1,    dx-2, dy-2);
  line(dx-2, dy-2, 1,    dy-2);
  line(1,    dy-2, 1,    1);
  setcolor(c0);
}

void Window::Set()
{
  SaveCurrentWindow();
  setviewport(X1,Y1, X2,Y2, Clip);
  viewporttype vp;
  getviewsettings(&vp);
  int c0 = getcolor();
  setcolor(WHITE);
  int dx = X2-X1+1;
  int dy = Y2-Y1+1;
  line(1,    1,    dx-2, 1);
  line(dx-2, 1,    dx-2, dy-2);
  line(dx-2, dy-2, 1,    dy-2);
  line(1,    dy-2, 1,    1);
  setcolor(c0);

  /*
  Xscal = (X2-X1)/(maxx-minx);
  Yscal = (Y2-Y1)/(maxy-miny);
  scal = (Xscal < Yscal)? Xscal: Yscal;
  */
  ///*

  int dX = X2-X1-XMARGIN;
  int dY = Y2-Y1-YMARGIN;

  //int smin = (X2-X1 < Y2-Y1)? X2-X1: Y2-Y1;
  int smin = (dX < dY)? dX: dY;
  double rmax = (maxx-minx > maxy-miny)? maxx-minx: maxy-miny;
  scal = smin/(rmax+.5);
  //*/
}

void Window::Show() {setviewport(X1,Y1, X2,Y2, Clip);}

MenuWindow::MenuWindow(int x1, int y1, int x2, int y2, int clip, int color):
            Window(x1,y1,x2,y2,clip) {
  items = new Queue<char>;
  TextColor = color;
  xcurr = 0;
  ycurr = 0;
}

MenuWindow::MenuWindow(): Window() {
  items = new Queue<char>;
  xcurr = 0;
  ycurr = 0;
}

MenuWindow::~MenuWindow()
{
  delete items; items=NULL;
}

void MenuWindow::DefineMenuWindow(int x1, int y1, int x2, int y2, int clip, int color)
{
  Window::DefineWindow(x1,y1,x2,y2,clip);
  TextColor = color;
  xcurr = 0;
  ycurr = 0;
}

void MenuWindow::Set() {
  Window::Set();
  Show();
}

void MenuWindow::writeln(char *item) {
  outtextxy(xcurr,ycurr,item);
  ycurr += 2*hText();
}

void MenuWindow::Show() {
  int color0 = getcolor();
  setcolor(TextColor);
  xcurr = 2*wText();
  ycurr = 2*hText();
  for (int i=0; i<items->Count(); i++) {
    char *item = (*items)[i];
    writeln(item);
  }
  setcolor(color0);
}

void MenuWindow::Add(char *Item) {
  char *item = new char[strlen(Item)+1];
  strcpy(item,Item);
  *items << item;
}

// class GManager

//GManager::GManager() : Megatile()
GManager::GManager() : MAHO()
{
  code_written = 0;
  MachineArea = 0;
  SheetArea = 0;
  Names = 0;
  Ruler = 0;
  hitCR = 0;

  work   = new Window();
  main   = new MenuWindow();
  mess   = new Window();
  entire = new Window();

  if (Color         < 0) Color        = BLACK;
  if (BkColor       < 0) BkColor      = LIGHTGRAY;
  if (TextColor     < 0) TextColor    = BLUE;
  if (MessageColor  < 0) MessageColor = RED;
  if (CutColor      < 0) CutColor     = RED;
  if (FlyColor      < 0) FlyColor     = LIGHTBLUE;
  if (ContourColor  < 0) ContourColor = WHITE;
  if (SideColor     < 0) SideColor    = RED;
  if (GatesColor    < 0) GatesColor   = LIGHTRED;
  if (SheetColor    < 0) SheetColor   = DARKGRAY;
  if (WorkAreaColor < 0) SheetColor   = BLUE;
  if (FinishColor   < 0) FinishColor  = YELLOW;
  // radius of finish circle, pels
  Rfinish = 10.;
}

GManager::~GManager()
{
  if (graphics) CloseGraph();

  //dheap("Before ~GManager");

  delete work; work=NULL;
  delete main; main=NULL;
  delete mess; mess=NULL;
  delete entire; entire=NULL;
  //cout << "~GManager done" << endl;
}

void GManager::InitGraph()
{
  // request auto detection
  int gdriver = DETECT, gmode, errorcode;

  /* register a driver that was added into graphics.lib */
  /* For information on adding the driver, see the
  /* BGIOBJ section of UTIL.DOC */
  errorcode = registerbgidriver(EGAVGA_driver);

  /* report any registration errors */
  if (errorcode < 0)
  {
    cout << "Graphics error: " << grapherrormsg(errorcode) << endl;
    cout << "Press any key to halt:";
    getch();
    exit(1); /* terminate with an error code */
  }

  // initialize graphics and local variables
  initgraph(&gdriver, &gmode, "");

  // read result of initialization
  errorcode = graphresult();
  // an error occurred
  if (errorcode != grOk)
  {
    cout << "Graphics error: " << grapherrormsg(errorcode) << endl;
    cout << "Press any key to halt:" << endl;
    getch();
    // terminate with an error code
    exit(1);
  }

  graphics = 1;

  MAXX = getmaxx();
  MAXY = getmaxy();
  MAXCOLOR = getmaxcolor() + 1;
  //setviewport(0,0,MAXX,MAXY,1);
  //clearviewport();

  YMARGIN = 3*textheight("W");
  XMARGIN = 7*textwidth("W");

  setcolor(Color);
  setbkcolor(BkColor);

  setcolor(RED);
}

void GManager::Area()
{
  Megatile::Area(minx,maxx,miny,maxy);
  if (MachineArea) {
    double zmin = MAHO_ZMIN - MAHO_Z0;
    double zmax = MAHO_ZMAX - MAHO_Z0;
    double xmin = MAHO_XMIN - MAHO_X0;
    double xmax = MAHO_XMAX - MAHO_X0;
    if (zmin < minx) minx = zmin;
    if (zmax > maxx) maxx = zmax;
    if (xmin < miny) miny = xmin;
    if (xmax > maxy) maxy = xmax;
  }

  /*
  // search through gates
  if (Gates) for (int ngate=0; ngate<gates->Count(); ngate++) {
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

  //double d = 0.15*(maxy-miny);
  double d = 0.05*(maxy-miny);

  /*
  minx = (minx < 0.)? SignAdd(minx,d): SignSub(minx,d);
  maxx = SignAdd(maxx,d);
  miny = (miny < 0.)? SignAdd(miny,d): SignSub(miny,d);
  maxy = SignAdd(maxy,d);
  */
  minx -= d;
  maxx += d;
  miny -= d;
  maxy += d;
}

void GManager::DefineWindows()
{
  Area();

  // windows
  entire->DefineWindow(0,0,MAXX,MAXY,1);

  int x1,y1,x2,y2, dx,dy,delta;
  dy = 2*Window::hText();
  delta = MAXX - MAXY - dy;
  x1 = delta;
  y1 = dy;
  x2 = MAXX;
  y2 = MAXY;
  work->DefineWindow(x1,y1,x2,y2,1);
  //work->Set();

  x1 = 0;
  y1 = dy;
  x2 = delta-1;
  y2 = MAXY;
  main->DefineMenuWindow(x1,y1,x2,y2,1,TextColor);
  //main->Set();

  x1 = 0;
  y1 = 0;
  x2 = MAXX;
  y2 = dy-1;
  mess->DefineWindow(x1,y1,x2,y2,1);
  //mess->Set();

  while(main->items->Count()) main->items->Kill(0);
  if (BackSide) {
    main->Add(" BACK SIDE");
    main->Add(" ---------");
    main->Add("");
  }
  main->Add("<CR>=start");
  main->Add("S=Slow");
  main->Add("F=Faster");
  main->Add("P=Pause");
  main->Add("G=Groove pause");
  main->Add("N=Names");
  main->Add("C=Colors");
  main->Add("A=sheet Area");
  main->Add("M=Machine Area");
  main->Add("T=Turned tiles");
  main->Add("L=ruLer");
  main->Add("H=Heap memory");
  main->Add("R=Repeat");
  main->Add("W=Write code");
  main->Add("Q=Quit");
  //main->Add("z=Zoom");
  //main->Add("Z=unZoom");
}

int GManager::Xscr(double x) {return scal*(x-minx);}
int GManager::Yscr(double y) {return YMARGIN+scal*(maxy-y);}

void GManager::Slower() {
  if (CuttingDelay == 0)  CuttingDelay   = 1;
  else                    CuttingDelay  *= 2;
  if (FlyingDelay == 0)   FlyingDelay    = 1;
  else                    FlyingDelay   *= 2;
  if (GroovingDelay == 0) GroovingDelay  = 4;
  else                    GroovingDelay += GroovingDelay/4;
}

void GManager::Faster() {
  CuttingDelay  /= 2;
  FlyingDelay   /= 2;
  GroovingDelay -= GroovingDelay/4;
}

void GManager::CloseGraph()
{
  closegraph();
  graphics = 0;
}

void GManager::Xtick(double x, double y, int h)
{
  int ix = Xscr(x);
  int iy = Yscr(y);
  line(ix,iy-h, ix,iy);
}

void GManager::Xlabl(double x, double y)
{
  int ix = Xscr(x);
  int iy = Yscr(y);
  int h = 2.*textheight("W");
  char number[16];
  if (x > 0.) sprintf(number, "%0.0f",  x);
  else        sprintf(number, "%0.0f ", x); // to center with '-'
  outtextxy(ix, iy-h, number);
}

void GManager::Ytick(double x, double y, int h)
{
  int ix = Xscr(x);
  int iy = Yscr(y);
  line(ix,iy, ix+h,iy);
}

void GManager::Ylabl(double x, double y)
{
  int ix = Xscr(x);
  int iy = Yscr(y);
  int h = 2.*textwidth("W");
  char number[16];
  sprintf(number, "%0.0f", y);
  outtextxy(ix+h, iy, number);
}

void GManager::GraphRuler()
{
  if (!Ruler) return;

  int color0 = getcolor();
  setcolor(ContourColor);

  struct textsettingstype t0;
  // get information about current text settings
  gettextsettings(&t0);

  int Xmin = 0;
  int Xmax = work->GetX2() - work->GetX1();
  int Ymin = 0;
  int Ymax = work->GetY2() - work->GetY1();
  int htick;
  double x,dx, y,dy;

  line(Xmin,Yscr(0.), Xmax,Yscr(0.));
  settextjustify(CENTER_TEXT, TOP_TEXT);

  htick=6;
  dx = 100.;
  for (x=-100.; Xscr(x)>0; x-=dx) {
    Xtick(x,0.,htick);
    Xlabl(x,0.);
  }
  if (abs(Xscr(50.)-Xscr(0.)) > 5)
  {
    dx = 50.;
    htick = 4;
    for (x=-50.; Xscr(x)>0; x-=dx) Xtick(x,0.,htick);

    if (abs(Xscr(10.)-Xscr(0.)) > 5)
    {
      dx = 10.;
      htick = 2;
      for (x=-10.; Xscr(x)>0; x-=dx) Xtick(x,0.,htick);
    }
  }

  htick=6;
  dx = 100.;
  for (x=100.; Xscr(x)<Xmax; x+=dx) {
    Xtick(x,0.,htick);
    Xlabl(x,0.);
  }
  if (abs(Xscr(50.)-Xscr(0.)) > 5)
  {
    dx = 50.;
    htick = 4;
    for (x=+50.; Xscr(x)<Xmax; x+=dx) Xtick(x,0.,htick);

    if (abs(Xscr(10.)-Xscr(0.)) > 5)
    {
      dx = 10.;
      htick = 2;
      for (x=+10.; Xscr(x)<Xmax; x+=dx) Xtick(x,0.,htick);
    }
  }

  line(Xscr(0.),Ymin, Xscr(0.),Ymax);
  settextjustify(LEFT_TEXT, CENTER_TEXT);

  htick = 6;
  dy = 100.;
  for (y=-100.; Yscr(y)<Ymax; y-=dy) {
    Ytick(0.,y,htick);
    Ylabl(0.,y);
  }
  if (abs(Yscr(50.)-Yscr(0.)) > 5)
  {
    dy = 50.;
    htick = 4;
    for (y=-50.; Yscr(y)<Ymax; y-=dy) Ytick(0.,y,htick);

    if (abs(Yscr(10.)-Yscr(0.)) > 5)
    {
      dy = 10.;
      htick = 2;
      for (y=-10.; Yscr(y)<Ymax; y-=dy) Ytick(0.,y,htick);
    }
  }

  htick = 6;
  dy = 100.;
  for (y=100.; Yscr(y)>0; y+=dy) {
    Ytick(0.,y,htick);
    Ylabl(0.,y);
  }
  if (abs(Yscr(50.)-Yscr(0.)) > 5)
  {
    dy = 50.;
    htick = 4;
    for (y=+50.; Yscr(y)>0; y+=dy) Ytick(0.,y,htick);

    if (abs(Yscr(10.)-Yscr(0.)) > 5)
    {
      dy = 10.;
      htick = 2;
      for (y=+10.; Yscr(y)>0; y+=dy) Ytick(0.,y,htick);
    }
  }

  settextstyle(t0.font, t0.direction, t0.charsize);
  settextjustify(t0.horiz, t0.vert);
  setcolor(color0);
}

//inline void GManager::GraphSheet()
void GManager::GraphSheet()
{
  if (SheetArea==0 && MachineArea==0 && ShowTurnedContour==0) return;

  // sheet area
  int color0 = getcolor();
  setcolor(SheetColor);
  line(Xscr(SheetX1),Yscr(SheetY1), Xscr(SheetX1),Yscr(SheetY2));
  line(Xscr(SheetX1),Yscr(SheetY2), Xscr(SheetX2),Yscr(SheetY2));
  line(Xscr(SheetX2),Yscr(SheetY2), Xscr(SheetX2),Yscr(SheetY1));
  line(Xscr(SheetX2),Yscr(SheetY1), Xscr(SheetX1),Yscr(SheetY1));
  setcolor(color0);

  // machine area
  setcolor(WorkAreaColor);
  int x1 = MAHO_XMIN-MAHO_X0;
  int x2 = MAHO_XMAX-MAHO_X0;
  int z1 = MAHO_ZMIN-MAHO_Z0;
  int z2 = MAHO_ZMAX-MAHO_Z0;
  line(Xscr(z1),Yscr(x1), Xscr(z1),Yscr(x2));
  line(Xscr(z1),Yscr(x2), Xscr(z2),Yscr(x2));
  line(Xscr(z2),Yscr(x2), Xscr(z2),Yscr(x1));
  line(Xscr(z2),Yscr(x1), Xscr(z1),Yscr(x1));
  setcolor(color0);

  // the center of rotation
  setcolor(SheetColor);
  circle(Xscr(MAHO_OZ-MAHO_Z0),Yscr(MAHO_OX-MAHO_X0),Rfinish);
  circle(Xscr(MAHO_OZ-MAHO_Z0),Yscr(MAHO_OX-MAHO_X0),Rfinish/2);
  circle(Xscr(MAHO_OZ-MAHO_Z0),Yscr(MAHO_OX-MAHO_X0),1);
  //putpixel(Xscr(MAHO_OZ-MAHO_Z0),Yscr(MAHO_OX-MAHO_X0),SheetColor);
  setcolor(color0);

  // point (0,0)
  setcolor(FlyColor);
  circle(Xscr(0.),Yscr(0.),Rfinish);
  circle(Xscr(0.),Yscr(0.),Rfinish/2);
  circle(Xscr(0.),Yscr(0.),1);
  //putpixel(Xscr(0.),Yscr(0.),FlyColor);
  setcolor(color0);
}

void GManager::GraphGates()
{
  int color0 = getcolor();
  setcolor(GatesColor);
  for (int i=0; i<gates->Count(); i++) {
    Side *gate = (*gates)[i];
    int x1 = Xscr(gate->p1.X);
    int y1 = Yscr(gate->p1.Y);
    int x2 = Xscr(gate->p2.X);
    int y2 = Yscr(gate->p2.Y);
    line(x1,y1,x2,y2);
  }
  setcolor(color0);
}

void GManager::GraphTiles()
{
  for (int ntd=0; ntd<tds->Count(); ntd++) GraphTile((*tds)[ntd]->tile);
}

void GManager::GraphNames()
{
  if (Names == 0) return;
  for (int ntd=0; ntd<tds->Count(); ntd++) GraphName((*tds)[ntd]->tile);
}

void GManager::GraphName(TileSide *tile)
{
  int color0 = getcolor();
  setcolor(ContourColor);
  Point upleft = TileUpLeftPoint(tile);
  int x = Xscr(upleft.X) + 5;
  int y = Yscr(upleft.Y) + 5;
  char *name = tile->nameAB;
  outtextxy(x,y,name);
  setcolor(color0);
}

void GManager::GraphTile(TileSide *tile)
{
  // graph contour
  GraphContour(tile);

  // graph turned contour
  if (ShowTurnedContour) GraphTurnedContour(tile);

  // graph sides
  GraphSides(tile);

  // graph groove
  GraphGrooveContour(tile, SheetColor);
  //GraphCylGroove(tile);
  //GraphBalGroove(tile);
}

void GManager::GraphContour(TileSide *tile)
{
  Queue<Side> *contour = tile->contour;
  int color0 = getcolor();
  setcolor(ContourColor);

  for (int i=0; i<contour->Count(); i++)
  {
    Side *side = (*contour)[i];
    line(Xscr(side->p1.X),Yscr(side->p1.Y),
         Xscr(side->p2.X),Yscr(side->p2.Y));
  }
  setcolor(color0);
}

void GManager::GraphTurnedContour(TileSide *tile)
{
  if (ShowTurnedContour)
  {
    if (tile->makegroove == 0) return;

    // to get the angle
    Groove* groove = tile->groove;

    Queue<Side> *contour = tile->contour;
    int color0 = getcolor();
    setcolor(ContourColor);

    struct linesettingstype oldline;
    getlinesettings(&oldline);
    setlinestyle(DOTTED_LINE, 1, 1);
    //setlinestyle(DASHED_LINE, 1, 1);

    for (int i=0; i<contour->Count(); i++) {
      Side *side = (*contour)[i];
      // rotate the side
      side->Rotate(O, +groove->angle);
      line(Xscr(side->p1.X),Yscr(side->p1.Y),
           Xscr(side->p2.X),Yscr(side->p2.Y));
      // rotate back
      side->Rotate(O, -groove->angle);
    }
    setlinestyle(oldline.linestyle, oldline.upattern, oldline.thickness);

    // graph turned groove
    groove->Rotate(O, +groove->angle);
    GraphGrooveContour(tile, ContourColor);
    groove->Rotate(O, -groove->angle);
    setcolor(color0);
  }
}

void GManager::GraphSides(TileSide *tile)
{
  Queue<Side>* sides = tile->sides;
  int color0 = getcolor();
  setcolor(SideColor);

  for (int i=0; i<sides->Count(); i++)
  {
    Side *side = (*sides)[i];
    if (side->status != lift) line(Xscr(side->p1.X),Yscr(side->p1.Y),
                                   Xscr(side->p2.X),Yscr(side->p2.Y));
  }
  setcolor(color0);
}

int GManager::GraphCylGroove(TileSide *tile)
{
  Groove* groove = tile->groove;
  int color0 = getcolor();

  // cyl. groove
  switch (groove->cylstatus) {
    case todo:  setcolor(CylGrooveColor);  break;
    case done:  setcolor(FlyColor);        break;
    case lift:  setcolor(SheetColor);      break;
  }
  int ireturn = GraphGroove(groove, 0);
  setcolor(color0);
  return ireturn;
}

int GManager::GraphBalGroove(TileSide *tile)
{
  Groove* groove = tile->groove;
  int color0 = getcolor();

  // ball groove
  switch (groove->balstatus) {
    case todo:  setcolor(BallGrooveColor); break;
    case done:  setcolor(FlyColor);        break;
    case lift:  setcolor(SheetColor);      break;
  }
  int ireturn = GraphGroove(groove, 10);
  setcolor(color0);
  return ireturn;
}

void GManager::GraphGrooveContour(TileSide *tile, int color)
{
  Groove* groove = tile->groove;
  int color0 = getcolor();
  setcolor(color);

  long delay0 = GroovingDelay;
  GroovingDelay = 0;
  int hitCR0 = hitCR;
  hitCR = 0;
  GraphGroove(groove, 0);
  GroovingDelay = delay0;
  hitCR = hitCR0;
  setcolor(color0);
}

void DashedArc(int x, int y, int theta1, int theta2, int R, int nparts)
{
  // should be theta2 > theta1

  if (nparts > 1) {
    if (theta1 > theta2) theta2 += 360;
    int dtheta = (theta2-theta1)/nparts;
    int th1,th2;
    th1 = theta1;
    th2 = th1+dtheta;
    for (int i=0; i<nparts; i+=2) {
      arc(x,y,th1,th2,R);
      th1 = th2 + dtheta;
      if (th1 > 360) th1 = 360-th1;
      th2 = th1 + dtheta;
      if (th2 > 360) th2 = 360-th2;
    }
  }
  else arc(x,y,theta1,theta2,R);
}

int GManager::GraphGroove(Groove* groove, int nparts)
{
  struct linesettingstype oldline;
  getlinesettings(&oldline);

  if (nparts > 1) {
    // change line style
    //setlinestyle(DOTTED_LINE, 1, 1);
    setlinestyle(DASHED_LINE, 1, 1);
  }

  int i = 1;   // skip initial point
  Point *p = (*groove)[i++];
  for(; i<groove->Count(); i++) {
    Point *next;
    if ((*groove)[i]->GetKey() == Point::Key()) {
      next = (*groove)[i];
      line(Xscr(p->X),Yscr(p->Y), Xscr(next->X),Yscr(next->Y));
    }
    else {
      Hole *r = (Hole*) (*groove)[i];
      next    = (*groove)[++i];

      const double eps   = 1.e-7;
      const double todeg = 180./M_PI;
      double dx,dy,theta1,theta2;

      dx = p->X - r->X;
      dy = p->Y - r->Y;
      theta1 = (fabs(dx)+fabs(dy) > MINFLOAT)? todeg*atan2(dy,dx): 0.;
      if (fabs(theta1) < eps) theta1 = 0.;
      if (theta1 < -eps) theta1 = 360. + theta1;

      dx = next->X - r->X;
      dy = next->Y - r->Y;
      theta2 = (fabs(dx)+fabs(dy) > MINFLOAT)? todeg*atan2(dy,dx): 0.;
      if (fabs(theta2) < eps) theta2 = 0.;
      if (theta2 < -eps) theta2 = 360. + theta2;

      // counterclockwise rotation
      Segment s1(*r, *p);
      Segment s2(*r, *next);
      int ccw = RightFrame(s1,s2);
      double R=scal*r->D;
      //if (ccw) arc(Xscr(r->X),Yscr(r->Y),theta1,theta2,R);
      //else     arc(Xscr(r->X),Yscr(r->Y),theta2,theta1,R);
      if (ccw) {
        //if (theta2 < theta1) theta2 = 360.-theta2;
        DashedArc(Xscr(r->X),Yscr(r->Y),theta1,theta2,R, nparts);
      }
      else {
        //if (theta1 < theta2) theta1 = 360.-theta1;
        DashedArc(Xscr(r->X),Yscr(r->Y),theta2,theta1,R, nparts);
      }
    }
    p = next;
  }

  // show points
  double Rpoint = 1;
  for (i=0; i<groove->Count(); i++) {
    Point *p = (*groove)[i];
    if (p->GetKey() == Point::Key()) setcolor(FinishColor);
    else                             setcolor(GatesColor);
    circle(Xscr(p->X),Yscr(p->Y),Rpoint);
  }
  setlinestyle(oldline.linestyle, oldline.upattern, oldline.thickness);
  if (kbhit()) {
    char c=getch();
    switch (c) {
      case '\0': getch(); break;
      case 'l':
      case 'L': Ruler=!Ruler; GraphRuler(); break;
      case 's':
      case 'S': Slower(); break;
      case 'f':
      case 'F': Faster(); break;
      case 'n':
      case 'N': Names=!Names; GraphNames(); break;
      case 'a':
      case 'A': SheetArea = !SheetArea; GraphSheet(); break;
      case 'w':
      case 'W': WriteCode(); break;
      case 'p':
      case 'P': if(Pause()) return 1; else break;
      case 'g':
      case 'G': hitCR = !hitCR; break;
      case 'h':
      case 'H': Heap(); break;
      case 'c':
      case 'C': Message("Possible colors:"); Colors(); break;
      case 'q':
      case 'Q': return 1;
    }
  }
  return 0;
}

void GManager::GraphCut(End *end) {
  const double eps = 1.e-6;
  int color0 = getcolor();
  setcolor(CutColor);
  LinePar cut(Point(FROMend->X,FROMend->Y),Point(end->X,end->Y));
  double r = end->Dcutter/2.;

  if (r < 1.) r=1.; // for graphics only

  double step = r/2.;
  if (fabs(cut.tmax - r) < eps) step /= 5.;
  double t=0.;
  do {
    if (kbhit()) {
      char c=getch();
      switch (c) {
        case '\0': getch(); break;
        case 'l':
        case 'L': Ruler = !Ruler; GraphRuler(); break;
        case 's':
        case 'S': Slower(); break;
        case 'f':
        case 'F': Faster(); break;
        case 'n':
        case 'N': Names=!Names; GraphNames(); break;
        case 'a':
        case 'A': SheetArea = !SheetArea; GraphSheet(); break;
        case 'w':
        case 'W': WriteCode(); delay(2000); break;
        case 'p':
        case 'P': if(Pause()) return; else break;
        case 'h':
        case 'H': Heap(); break;
        case 'g':
        case 'G': hitCR = !hitCR; break;
        case 'c':
        case 'C': Message("Possible colors:"); Colors(); break;
        case 'q':
        case 'Q': Quit(); return;
      }
    }
    delay(CuttingDelay);
    int x = Xscr(cut.x0 + cut.cosx*t);
    int y = Yscr(cut.y0 + cut.cosy*t);
    int radius = scal*r;
    if (radius < 1) radius = 1;
    circle(x,y,radius);
    t = t+step;
  } while (t <= cut.tmax);
  setcolor(color0);
}

//inline void GManager::GraphFly(End *end) {
void GManager::GraphFly(End *end) {
  int color0 = getcolor();
  setcolor(FlyColor);
  line(Xscr(FROMend->X),Yscr(FROMend->Y),Xscr(end->X),Yscr(end->Y));

  delay(FlyingDelay);

  /*
  LinePar cut(Point(FROMend->X,FROMend->Y),Point(end->X,end->Y));
  double r = end->Dcutter/2.;
  double step = r/2.;
  //if (fabs(cut.tmax - r) < eps) step /= 10.;
  double t=0.;
  do {
    if (kbhit()) return;
    //delay(CuttingDelay);
    int x = Xscr(cut.x0 + cut.cosx*t);
    int y = Yscr(cut.y0 + cut.cosy*t);
    //int radius = scal*r;
    int radius = 1;
    //if (radius < 1) radius = 1;
    circle(x,y,radius);
    t = t+step;
  } while (t <= cut.tmax);
  */
  setcolor(color0);
}

inline void GManager::GraphMakeHole(End *end)
{
  int color0 = getcolor();
  setcolor(CutColor);
  circle(Xscr(end->X),Yscr(end->Y),scal*end->Dcutter/2.);
  setcolor(color0);
}

inline void GManager::GraphLiftCutter() {}

int GManager::GraphCutting()
{
  if (ends->Count() == 0) return 1;
  int ireturn = 0;

  int idown, iends;

  // perform the down cut
  if (down_cut) {
    iends = 0;
    // initialize the FROMend from down_ends
    idown = 0;
    FROMend = (*down_ends)[idown++];

    for (; idown<down_ends->Count(); idown++) {

      if (graphics == 0) return 1; // after Pause Quit

      End *end = (*down_ends)[idown];
      switch (end->status) {
        case tocut:      GraphCut(end); break;
        case tofly:      GraphFly(end); break;
        case tosemihole: GraphCut(end); break;
        case tohole:     GraphFly(end);
                         GraphMakeHole(end);
                         GraphLiftCutter();
                         break;
      }
      FROMend = end;
    }
  }
  else {
    // initialize the FROMend from ends. Omit the first.
    iends = 1;
    FROMend = (*ends)[iends++];
  }

  // grooving
  for (int ngroove=0; ngroove<groove_ends->Count(); ngroove++) {
    EndGroove *egroove = (*groove_ends)[ngroove];
    TileSide *tile = egroove->tile;
    switch (egroove->cutype) {
      case cylcut: ireturn = GraphCylGroove(tile); break;
      case balcut: ireturn = GraphBalGroove(tile); break;
    }
    if (ireturn) return 1;
    if (hitCR) {
      if (Pause()) return 1;
    }
    else delay(GroovingDelay);
  }

  // cutting
  for (; iends<ends->Count(); iends++) {

    if (graphics == 0) return 1; // after Pause Quit

    End *end = (*ends)[iends];
    switch (end->status) {
      case tocut:      GraphCut(end); break;
      case tofly:      GraphFly(end); break;
      case tosemihole: GraphCut(end); break;
      case tohole:     GraphFly(end);
                       GraphMakeHole(end);
                       GraphLiftCutter();
                       break;
    }
    FROMend = end;
  }
  // circle at the end point
  int c0 = getcolor();
  setcolor(FinishColor);
  circle(Xscr(FROMend->X),Yscr(FROMend->Y), Rfinish);
  setcolor(c0);
  return 0;
}

void GManager::Message(char *str) {
  mess->SaveCurrentWindow();
  mess->Set();
  mess->Clear();
  int color0 = getcolor();
  setcolor(MessageColor);
  outtextxy(12,4,str);
  mess->RestoreSavedWindow();
  setcolor(color0);
}

void GManager::Message(char *str1, char *str2) {
  mess->SaveCurrentWindow();
  mess->Set();
  mess->Clear();
  int color0 = getcolor();
  setcolor(MessageColor);

  char *s = new char[strlen(str1)+strlen(str2)+1];
  strcpy(s,str1);
  strcat(s,str2);
  outtextxy(12,4,s);
  delete[] s;

  mess->RestoreSavedWindow();
  setcolor(color0);
}

void GManager::WriteCode()
{
  if (nopath) {
    Message("MAHO code can't be written");
    //wait for a second
    //delay(2000);
    return;
  }

  if (!code_written) {
    MAHO::WriteCode();
    code_written = 1;

    Message("MAHO code has been written now to file ", code_file);
    //wait for a second
    //delay(2000);
  }
  else {
    Message("MAHO code already has been written to file ", code_file);
    //wait for a second
    //delay(2000);
  }
}

int GManager::Pause()
{
  AGAIN:
  ;
  Message("Press <CR> to continue");

  char c=getch();
  switch (c) {
    case '\r': Message(""); break;
    case '\0': getch(); goto AGAIN;
    case 'l' :
    case 'L' : Ruler = !Ruler; GraphRuler(); goto AGAIN;
    case 's' :
    case 'S' : Slower(); goto AGAIN;
    case 'f' :
    case 'F' : Faster(); goto AGAIN;
    case 'n' :
    case 'N' : Names=!Names; GraphNames(); goto AGAIN;
    case 'a' :
    case 'A' : SheetArea = !SheetArea; GraphSheet(); goto AGAIN;
    case 'm' :
    case 'M' : MachineArea=!MachineArea; goto AGAIN;
    case 't' :
    case 'T' : ShowTurnedContour = !ShowTurnedContour; goto AGAIN;
    case 'c' :
    case 'C' : Message("Possible colors:"); Colors(); goto AGAIN;
    case 'w' :
    case 'W' : WriteCode(); delay(2000); goto AGAIN;
    case 'g' :
    case 'G' : hitCR = !hitCR; goto AGAIN;
    case 'h' :
    case 'H' : Heap(); goto AGAIN;
    case 'q' :
    case 'Q' : Quit(); return 1;
    default  : goto AGAIN;
  }
  Message("");
  return 0;
}

void GManager::Quit()
{
  if (graphics == 0) return;

  main->Set();
  main->Clear();

  if (!code_written)
  {
    if (nopath) {
      Message("MAHO code can't be written. Hit any key to quit");
      getch();
    }
    else {
      AGAIN:
      Message("MAHO code has not been written yet. "
              "Hit W to write or Q to quit");
      char c=getch();
      switch (c) {
        case '\0': getch(); goto AGAIN;
        case 'w' :
        case 'W' : WriteCode(); delay(2000); break;
        case 'q' :
        case 'Q' : break;
        default: goto AGAIN;
      }
    }
  }
  if (!code_written) {
    Message("Exit without saving of MAHO code. Press any key to quit.");
    getch();
  }

  if (FinishSnow) Snow();
  CloseGraph();
  return;
}

void GManager::Menu()
{
  MakePath();
  MakeCode();

  //PrintCutsEnds();

  //MAHO::WriteCode();
  //return;

  cout << endl << "To make sloped groove tile ";
  if (TrueAngle) cout << "will be turned to its real angle" << endl;
  else           cout << "will be turned to +/- 10 degrees" << endl;
  cout << endl;

  if (BackSide) cout << "Process back side of sheet" << endl << endl;

  cout << "<CR>=Visual control, Q=Quit" << endl;
  char c=getchar();
  if (c=='Q' || c=='q') goto EXIT;

  Graph();
  EXIT:;
}

void GManager::Snow()
{
  entire->Set();

  time_t t0,t,dt;
  dt = 15;
  t0 = time(NULL);
  t  = t0+dt;
  while (time(NULL) < t) {
    for (int i=0; i<16384; i++) {
      putpixel(random(MAXX),random(MAXY),random(10));
    }
    if (kbhit()) {
      getch();
      return;
    }
  }
}

void GManager::Heap()
{
  char s[64];
  sprintf(s, "Used heap memory %ld bytes", dheap());
  Message(s);
}

void GManager::Graph()
{
  // initialize the graphics system
  InitGraph();

  char c = '\0';

  if (nopath) SheetArea=1;

  INIT:
  ;
  DefineWindows();

  mess->Set();
  main->Set();
  main->Show();

  work->Set();

  REPEAT:
  ;
  work->Clear();

  if (BackSide) Message("BACK SIDE.         Press <CR> to start");
  else          Message(                   "Press <CR> to start");

  /*
  if (Ruler)             GtaphRuler();
  if (Names)             GtaphNames();
  if (MachineArea)       GraphSheet();
  if (SheetArea)         GraphSheet();
  if (ShowTurnedContour) GraphSheet();
  */

  GraphRuler();
  GraphGates();
  GraphTiles();
  GraphNames();
  GraphSheet();

  AGAIN:
  ;
  c=getch();
  switch (c) {
    case '\0': getch(); goto AGAIN;
    case '\r': Message(""); if(GraphCutting()==0) goto AGAIN; else break;
    case 'l' :
    case 'L' : Ruler = !Ruler; GraphRuler();
               if (!Ruler) goto REPEAT; else goto AGAIN;
    case 's' :
    case 'S' : Slower(); goto AGAIN;
    case 'f' :
    case 'F' : Faster(); goto AGAIN;
    case 'n' :
    case 'N' : Names=!Names; GraphNames();
               if (!Names) goto REPEAT; else goto AGAIN;
    case 'a' :
    case 'A' : SheetArea = !SheetArea; GraphSheet();
               if (!SheetArea) goto REPEAT; else goto AGAIN;
    case 'm' :
    case 'M' : MachineArea=!MachineArea; goto INIT;
    case 't' :
    case 'T' : ShowTurnedContour = !ShowTurnedContour; goto INIT;
    case 'c' :
    case 'C' : Message("Possible colors:"); Colors(); goto AGAIN;
    case 'w' :
    case 'W' : WriteCode(); goto AGAIN;
    case 'r' :
    case 'R' : goto REPEAT;
    case 'g' :
    case 'G' : hitCR = !hitCR; goto AGAIN;
    case 'h' :
    case 'H' : Heap(); goto AGAIN;
    case 'q' :
    case 'Q' : break;
    default  : goto AGAIN;
  }
  Quit();
}
