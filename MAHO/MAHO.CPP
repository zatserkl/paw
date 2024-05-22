// Andriy_Zatserklyaniy@fnal.gov

#if !defined(__STDLIB_H)
#include <stdlib.h>
#endif

#if !defined(__STDIO_H)
#include <stdio.h>
#endif

#if !defined(__MATH_H)
#include <math.h>
#endif

#if !defined(__CTYPE_H)
#include <ctype.h>
#endif

#if !defined(__STRING_H)
#include <string.h>
#endif

#if !defined(__FSTREAM_H)
#include <fstream.h>
#endif

#if !defined(__IOSTREAM_H)
#include <iostream.h>
#endif

#include "wait.h"

#include "queue.h"
#include "elements.h"

#include "megatile.h"
#include "MAHO.h"

// General-purpose function Clip

inline char* Clip(char *string) {
  char *s = new char[strlen(string)+1];
  strcpy(s,string);
  return s;
}

// class MAHO

// Declaration and defenition of static variables

int MAHO::G28on = 0;

// Machine work area in machine frame
double MAHO::MAHO_XMIN = 0.   ;
double MAHO::MAHO_XMAX = 1000.;
double MAHO::MAHO_ZMIN = 0.   ;
double MAHO::MAHO_ZMAX = 800. ;

double MAHO::MAHO_X0;
double MAHO::MAHO_Y0;
double MAHO::MAHO_Z0;

double MAHO::MAHO_OX;
double MAHO::MAHO_OZ;

double MAHO::ChangeCutter_Y;
double MAHO::ChangeCutter_Z;

int MAHO::NumberDcut;
int MAHO::NumberDsemihole;
int MAHO::NumberDhole;
int MAHO::NumberDcyl;
int MAHO::NumberDbal;
// the height of 'fly'
double MAHO::hfly;
// the height of cutter top position
double MAHO::htop;
// rotation
long MAHO::rotation    = 2000;
// direction of rotation
int MAHO::dir_rotation = 3;

Queue<char>   *MAHO::comment  = NULL;

double         MAHO::thickness;
double         MAHO::minheight;
int            MAHO::correct_depth = 1;

Queue<double> *MAHO::depth  = NULL;
Queue<long>   *MAHO::feed   = NULL;
Queue<double> *MAHO::cdepth = NULL;
Queue<long>   *MAHO::cfeed  = NULL;
Queue<double> *MAHO::bdepth = NULL;
Queue<long>   *MAHO::bfeed  = NULL;
// No. of prigram
long MAHO::ProgramNumber;
// No. of prigram
long MAHO::SubNumber;
// increment of the program line number
long MAHO::LineInc = 1;

// code file name
char MAHO::code_file[33];

// class MAHO

MAHO::MAHO() : Megatile(MAHO_OZ-MAHO_Z0, MAHO_OX-MAHO_X0)
{
  CurrentFeed   = 0;
  CurrentDepth  = hfly;
  CurrentCutter = noncut;
  CurrentAngle  = 0.;
  CurrentTile   = NULL;

  if (correct_depth)
  {
    // correct the last cut depth for minheight-thickness
    double *lastdepth = depth->Tail();
    double oldlastdepth = *lastdepth;
    *lastdepth = minheight - thickness;
    cout << endl
         <<"MAHO: replace the last cut depth " << oldlastdepth
         << " by " << *lastdepth << endl;
  }

  // initial value
  NoRotation = 1;
  NoFeed = 1;

  CutterStatus = fly;
  G0G1status   = G00;

  code = new Queue<char>;
  main = new Queue<char>;
  sub  = new Queue<char>;

  BUF[0] = '\0';
}

MAHO::~MAHO() {
  delete depth;   depth   = NULL;
  delete feed;    feed    = NULL;
  delete cdepth;  cdepth  = NULL;
  delete cfeed;   cfeed   = NULL;
  delete bdepth;  bdepth  = NULL;
  delete bfeed;   bfeed   = NULL;
  delete comment; comment = NULL;

  delete code;    code = NULL;
  delete main;    main = NULL;
  delete sub;     sub  = NULL;
  //MESS("~MAHO done");
}

//inline void MAHO::PUT()
void MAHO::PUT()
{
  char buf[256];
  int nstring;
  if (part == &main) {
    nstring = (nmain > 0)? nmain: ProgramNumber;
    nmain += LineInc;
  }
  else {
    nstring = (nsub  > 0)?  nsub: SubNumber;
    nsub  += LineInc;
  }
  sprintf(buf, "N%d,", nstring);
  strcat(buf,BUF);
  char *s = new char[strlen(buf)+1];
  strcpy(s,buf);
  **part << s;
}

// MAHO commands

//inline void MAHO::Stop()
void MAHO::Stop()
{
  if (CutterStatus == cut) LiftCutter();
  if (G0G1status == G0) sprintf(BUF,    "Y%0.3f,M30",htop);
  else                  sprintf(BUF, "G0,Y%0.3f,M30",htop);
  PUT();
  CutterStatus = fly;
  G0G1status   = G0;
}

//inline void MAHO::Pause(char *mess)
void MAHO::Pause(char *mess)
{
  LiftCutter();
  if (G0G1status == G0) sprintf(BUF,    "Y%0.3f,M0,(%s)", htop,mess);
  else                  sprintf(BUF, "G0,Y%0.3f,M0,(%s)", htop,mess);
  PUT();

  // fast down cutter to hfly (previous was G0)
  sprintf(BUF, "Y%0.3f",hfly);
  PUT();
  CutterStatus = fly;
  G0G1status   = G0;
}

void MAHO::Rotate(const double angle)
{
  if (CurrentAngle == angle) return;

  LiftCutter();

  int large_turn = (fabs(CurrentAngle-angle) > 1.);

  if (large_turn) {
    // Lift cutter at the safety height ChangeCutter_Y
    if (G0G1status == G0) sprintf(BUF,    "Y%0.3f", ChangeCutter_Y);
    else                  sprintf(BUF, "G0,Y%0.3f", ChangeCutter_Y);
    PUT();
    G0G1status   = G0;
  }

  // make a turn
  if (G0G1status == G0) sprintf(BUF,    "B%0.3f", angle);
  else                  sprintf(BUF, "G0,B%0.3f", angle);
  PUT();
  CurrentAngle = angle;

  if (large_turn) {
    // fast down cutter back at the hfly (current G0G1status == G0)
    sprintf(BUF,    "Y%0.3f", hfly); PUT();
  }
  CutterStatus = fly;
  G0G1status   = G0;
}

//inline void MAHO::LiftCutterFast()
void MAHO::LiftCutterFast()
{
  // uses G0 only

  if (G0G1status == G0) sprintf(BUF,    "Y%0.3f", hfly);
  else                  sprintf(BUF, "G0,Y%0.3f", hfly);
  PUT();

  CutterStatus = fly;
  G0G1status   = G0;
}

//--- commands which use G1 ---

//inline void MAHO::LiftCutter()
void MAHO::LiftCutter()
{
  // uses G1 only
  char RotationStr[32];
  char FeedStr[32];

  if (CutterStatus != fly) {
    if (G0G1status == G1) {
      sprintf(BUF,    "Y%0.3f", hfly);
      if (NoFeed) {
        sprintf(FeedStr, ",F%ld", CurrentFeed);
        strcat(BUF,FeedStr);
        NoFeed = 0;
      }
      if (NoRotation) {
        sprintf(RotationStr, ",S%ld,M%d",rotation,dir_rotation);
        strcat(BUF,RotationStr);
        NoRotation = 0;
      }
      PUT();
    }
    else {
      sprintf(BUF, "G1,Y%0.3f", hfly);
      if (NoFeed) {
        sprintf(FeedStr, ",F%ld", CurrentFeed);
        strcat(BUF,FeedStr);
        NoFeed = 0;
      }
      if (NoRotation) {
        sprintf(RotationStr, ",S%ld,M%d",rotation,dir_rotation);
        strcat(BUF,RotationStr);
        NoRotation = 0;
      }
      PUT();
    }

    CutterStatus = fly;
    G0G1status   = G1;
  }
}

//inline void MAHO::DownCutter()
void MAHO::DownCutter()
{
  // uses G1 only
  char RotationStr[32];
  char FeedStr[32];

  if (part == &main) {
    // Down cutter with respect to CurrentDepth
    if (G0G1status == G1) {
      sprintf(BUF,    "Y%0.3f", CurrentDepth);
    }
    else {
      sprintf(BUF, "G1,Y%0.3f", CurrentDepth);
    }

    if (NoFeed) {
      sprintf(FeedStr, ",F%ld", CurrentFeed);
      strcat(BUF,FeedStr);
      NoFeed = 0;
    }
    if (NoRotation) {
      sprintf(RotationStr, ",S%ld,M%d",rotation,dir_rotation);
      strcat(BUF,RotationStr);
      NoRotation = 0;
    }
    PUT();
  }
  else {
    // Down cutter with respect to subroutine parameter E1
    if (G0G1status == G1) {
      sprintf(BUF,    "Y=E1");
    }
    else {
      sprintf(BUF, "G1,Y=E1");
    }

    if (NoFeed) {
      sprintf(FeedStr, ",F%ld", CurrentFeed);
      strcat(BUF,FeedStr);
      NoFeed = 0;
    }
    if (NoRotation) {
      sprintf(RotationStr, ",S%ld,M%d",rotation,dir_rotation);
      strcat(BUF,RotationStr);
      NoRotation = 0;
    }
    PUT();
  }
  CutterStatus = cut;
  G0G1status = G1;
}

//inline void MAHO::MoveHeight(double h)
void MAHO::MoveHeight(double h)
{
  // uses G1 only
  char RotationStr[32];
  char FeedStr[32];

  if (G0G1status == G1) {
    sprintf(BUF,    "Y%0.3f", h);
    if (NoFeed) {
      sprintf(FeedStr, ",F%ld", CurrentFeed);
      strcat(BUF,FeedStr);
      NoFeed = 0;
    }
    if (NoRotation) {
      sprintf(RotationStr, ",S%ld,M%d",rotation,dir_rotation);
      strcat(BUF,RotationStr);
      NoRotation = 0;
    }
    PUT();
  }
  else {
    sprintf(BUF, "G1,Y%0.3f", h);
    if (NoFeed) {
      sprintf(FeedStr, ",F%ld", CurrentFeed);
      strcat(BUF,FeedStr);
      NoFeed = 0;
    }
    if (NoRotation) {
      sprintf(RotationStr, ",S%ld,M%d",rotation,dir_rotation);
      strcat(BUF,RotationStr);
      NoRotation = 0;
    }
    PUT();
  }

  CutterStatus = fly;
  G0G1status   = G1;
}

//inline void MAHO::Cut(const End *end)
void MAHO::Cut(const End *end)
{
  // Cut at CurrentDepth
  char RotationStr[32];
  char FeedStr[32];

  if (CurrentCutter != end->cutype) SetCutterAtEnd(end);

  if (CutterStatus == fly) DownCutter();

  if (G0G1status == G1) {
    sprintf(BUF,    "X%0.3f,Z%0.3f", end->Y,end->X);
    if (NoFeed) {
      sprintf(FeedStr, ",F%ld", CurrentFeed);
      strcat(BUF,FeedStr);
      NoFeed = 0;
    }
    if (NoRotation) {
      sprintf(RotationStr, ",S%ld,M%d",rotation,dir_rotation);
      strcat(BUF,RotationStr);
      NoRotation = 0;
    }
  }
  else {
    sprintf(BUF, "G1,X%0.3f,Z%0.3f,F%ld", end->Y,end->X,CurrentFeed);
    if (NoFeed) {
      sprintf(FeedStr, ",F%ld", CurrentFeed);
      strcat(BUF,FeedStr);
      NoFeed = 0;
    }
    if (NoRotation) {
      sprintf(RotationStr, ",S%ld,M%d",rotation,dir_rotation);
      strcat(BUF,RotationStr);
      NoRotation = 0;
    }
  }
  // add to BUF the end master name
  if (end->tile != CurrentTile) {
    if (end->tile == NULL) strcat(BUF,",(DOWN,CUT)");
    else {
      strcat(BUF,",(TILE,");
      strcat(BUF,GetName(end->tile));
      strcat(BUF,")");
    }
    CurrentTile = end->tile;
  }
  PUT();

  CutterStatus = cut;
  G0G1status   = G1;
}

//---

//inline void MAHO::Fly(const End *end)
void MAHO::Fly(const End *end)
{
  // Possible recurrencive call through SetCutterAtEnd

  if (CurrentCutter != end->cutype) SetCutterAtEnd(end);
  else {
    LiftCutter();
    if (G0G1status == G0) sprintf(BUF,    "X%0.3f,Z%0.3f", end->Y,end->X);
    else                  sprintf(BUF, "G0,X%0.3f,Z%0.3f", end->Y,end->X);
    // add to BUF the end master name
    if (end->tile != CurrentTile) {
      if (end->tile == NULL) strcat(BUF,",(DOWN,CUT)");
      else {
        strcat(BUF,",(TILE,");
        strcat(BUF,GetName(end->tile));
        strcat(BUF,")");
      }
      CurrentTile = end->tile;
    }
    PUT();

    CutterStatus = fly;
    G0G1status   = G0;
  }
}

//inline void MAHO::SetCutterAtEnd(const End *end)
void MAHO::SetCutterAtEnd(const End *end)
{
  if (CurrentCutter == end->cutype) return;

  LiftCutter();

  int Number;
  switch (end->cutype) {
    case cutcut: Number=NumberDcut; break;
    case cylcut: Number=NumberDcyl; break;
    case balcut: Number=NumberDbal; break;
    default: STOPCR("MAHO::SetCutterAtEnd(EndCutter cutter): Unknown cutter");
  }
  CurrentCutter = end->cutype;

  // go to change cutter point
  if (G0G1status == G0) {
    sprintf(BUF,    "Y%0.3f,Z%0.3f,M0",ChangeCutter_Y,ChangeCutter_Z);
  }
  else {
    sprintf(BUF, "G0,Y%0.3f,Z%0.3f,M0",ChangeCutter_Y,ChangeCutter_Z);
  }
  switch (end->cutype) {
    case cutcut: strcat(BUF, ",(SET,CUT,CUTTER)");         break;
    case cylcut: strcat(BUF, ",(SET,CYLINDRICAL,CUTTER)"); break;
    case balcut: strcat(BUF, ",(SET,BALL,CUTTER)");        break;
    default: STOPCR("MAHO::SetCutterAtEnd(EndCutter cutter): Unknown cutter");
  }
  PUT();
  sprintf(BUF, "G18,T%d,M67",Number); PUT();

  CutterStatus = fly;
  G0G1status   = G0;

  // fly back to the end
  Fly(end);

  // fast down to hfly
  sprintf(BUF, "Y%0.3f", hfly); PUT();

  CutterStatus = fly;
  G0G1status   = G0;
}

//inline void MAHO::MakeHole(const End *end)
#pragma warn -par
void MAHO::MakeHole(const End *end)
{
  STOP("MAHO::MakeHole(): Possibility isn't realized yet!");
}
#pragma warn +par

void MAHO::MakeCode()
{
  /*
  // Using subroutine
  MakeSubMain();
  MakeSub();
  // main
  //while (main->Count()) *code << main->Remove(0);
  // sub
  //while ( sub->Count()) *code <<  sub->Remove(0);
  */

  // Main only
  MakeMain();
  while (main->Count()) *code << main->Remove(0);
}

void MAHO::MakeSubMain()
{
  MainHeader();

  // PUT to main
  part = &main;

  // metric system
  sprintf(BUF, "G71"); PUT();
  // set the reference size system
  sprintf(BUF, "G90"); PUT();

  // Loop over npass

  for (int i=0; i<depth->Count(); i++) {
    SetDepth(*(*depth)[i]);
    // call of subroutine with parameter E1 = CurrentDepth
    sprintf(BUF, "G22,N%ld,E1=%0.3f(DEPTH%0.3f)",
                   SubNumber, CurrentDepth, CurrentDepth);
    PUT();
  }
  Stop();
}

void MAHO::SubHeader()
{
  /*
  // title of sub module
  sprintf(BUF, "%%MM");
  *sub << Clip(BUF);
  */
  // the number of subroutine
  sprintf(BUF, "N%ld,(CUTTING)", SubNumber);
  *sub << Clip(BUF);
  nsub = 10;
}

void MAHO::MakeSub()
{
  // PUT to sub
  part = &sub;
  //CutPath();
}

void MAHO::MainHeader()
{
  // title of main module
  sprintf(BUF, "%%PM");
  //sprintf(BUF, "%%");     // to satisfy to current txt --> bin convertor
  *main << Clip(BUF);

  // the number of program
  //sprintf(BUF, "N%ld,(AUTOGENERETED)", ProgramNumber);

  char *Mname = new char[strlen(GetName(ends->Head()->tile))+1];
  strcpy(Mname,GetName(ends->Head()->tile));
  char *p = strstr(Mname,"-");
  p = strstr(p+1,"-");
  *p = '\0';

  char *ftask = new char[strlen(task_file)+1];
  #pragma warn -pia
  for (int i=0; ftask[i] = toupper(task_file[i]); i++);
  #pragma warn +pia

  sprintf(BUF, "N%ld,(TASK,%s,MEGATILE,%s", ProgramNumber,ftask,Mname);
  //sprintf(BUF, "N%ld,(MEGATILE,%s", ProgramNumber,Mname);
  if (BackSide) strcat(BUF, ",BACK,SIDE)");
  else          strcat(BUF,           ")");
  *main << Clip(BUF);

  delete[] Mname;
  delete[] ftask;

  nmain = 1;

  /*
  char *ftask = new char[strlen(task_file)+1];
  #pragma warn -pia
  for (int i=0; ftask[i] = toupper(task_file[i]); i++);
  #pragma warn +pia

  sprintf(BUF, "N%d,(TASK,FILE,%s)", nmain++, ftask);
  *main << Clip(BUF);
  //delete[] ftask;
  */

  if (comment) for (int i=0; i<comment->Count(); i++) {
    sprintf(BUF, "N%d,(%s)", nmain++, (*comment)[i]);
    *main << Clip(BUF);
  }

  // WorkPiece zero
  double WorkPieceY0 = MAHO_Y0 + thickness;
  sprintf(BUF, "N%d,G93,X%0.3f,Y%0.3f,Z%0.3f",
                    nmain++, MAHO_X0,WorkPieceY0,MAHO_Z0);
  *main << Clip(BUF);

  // metric system
  sprintf(BUF, "N%d,G71", nmain++);
  *main << Clip(BUF);
  // set the reference size system
  sprintf(BUF, "N%d,G90", nmain++);
  *main << Clip(BUF);
  // select the plane
  sprintf(BUF, "N%d,G18", nmain++);
  *main << Clip(BUF);

  if (G28on) {
    sprintf(BUF, "N%d,G28", nmain++);
    *main << Clip(BUF);
  }

  // unconditional set angle 0.
  sprintf(BUF, "N%d,G0,B0", nmain++);
  *main << Clip(BUF);
  CurrentAngle = 0.;
  G0G1status   = G0;

  // Graphics options
  /*
  //End *pXmax, *pZmax;
  //End *end = (*ends)[0];
  //pXmax = pZmax = end;
  //nends = ends->Count();
  //for (int nend=1; nend<nends; nend++) {
  //  end = (*ends)[nend];
  //  if (fabs(end->X) > fabs(pXmax->X)) pXmax = end;
  //  if (fabs(end->Y) > fabs(pZmax->Y)) pZmax = end;
  //}
  //Xmax = pXmax->X;
  //Zmax = pZmax->Y;
  //Ymax = depth[npass-1]; if (Ymax > -20.) Ymax=-20.;

  double d = 10.; // extra space in mm
  double Zmax = (fabs(SheetX1) > fabs(SheetX2))? SheetX1: SheetX2;
  double Xmax = (fabs(SheetY1) > fabs(SheetY2))? SheetY1: SheetY2;
  double Ymax = *(*depth)[depth->Count()-1];

  sprintf(BUF, "N%d,G98,X%0.0f,Y%0.0f,Z%0.0f,B-30,I150,J25,K200",
          nmain++,Xmax+sign(Xmax)*d, Ymax+sign(Ymax)*d, Zmax+sign(Zmax)*d);
  sprintf(BUF, "N%d,G99,X%0.0f,Y%0.0f,Z%0.0f,B-30,I150,J25,K200",
          nmain++,Xmax,Ymax,Zmax);
  */

  /*
  // set the number of the next main line
  nmain = 10;
  */
}

void MAHO::MakeMain()
{
  MainHeader();

  // PUT to main
  part = &main;

  // initialize the CurrentFeed by first cut feed
  SetFeed(*feed->Head());

  if (continue_level) Pause("CHECK,TOOL,MEMORY");

  if (ToCut(down_ends))
  {
    // make down cut

    // to get comment "(DOWN,CUT)" set CurrentTile to some tile from the ends
    for (int i=0; i<ends->Count(); i++) {
      End *end = (*ends)[i];
      if (end->tile != NULL) {
        CurrentTile = end->tile;
        break;
      }
    }

    CutPath(down_ends);

    if (!ToCut(groove_ends)) Pause("REMOVE,EXTRA,SHEET");
  }

  // make grooves
  CutGrooves();

  // cut the tiles
  CutPath(ends);

  Stop();
}

void MAHO::CutPath(const Queue<End> *path)
{
  if (!ToCut(path)) return;

  // set the angle
  Rotate(0.);

  Queue<End> *cache = new Queue<End>;

  int i;
  // search the first cutting end
  for (i=0; i<path->Count(); i++) {
    End *end = (*path)[i];
    if (end->status != tofly) break;
    else {
      // store fly end as a head of cache
      while (cache->Count()) cache->Remove(0);
      *cache << end;
    }
  }

  if (cache->Count() == 0) STOPCR("CutPath: path should start from fly end");

  while (i<path->Count())
  {
    End *end = (*path)[i];
    if (end->status == tofly) {
      CutCache(cache);
      // Clear existing cache. The first end of new cache always will be fly
      while (cache->Count()) cache->Remove(0);
    }
    *cache << end;
    i++;
  }

  // cut the last cache
  CutCache(cache);
  // clear cache
  while (cache->Count()) cache->Remove(0);
  delete cache;
}

void MAHO::CutCache(const Queue<End> *cache)
{
  if (!ToCut(cache)) return;

  // perform the first pass
  int npass = 0;
  SetDepth(*(*depth)[0]);
  SetFeed(*(*feed)[0]);

  if (CurrentCutter != cache->Head()->cutype)
  {
    // set the right cutter and fly to the start cache end
    SetCutterAtEnd(cache->Head());
  }
  else
  {
    // simple fly to the start cache end
    Fly(cache->Head());
  }

  // down the cutter to the first depth
  DownCutter();

  // store the number of the next line
  int nline1 = NextLine();

  for (int i=1; i<cache->Count(); i++) {
    End *end = (*cache)[i];
    switch (end->status) {
      case tocut:      Cut(end); break;
      case tofly:      Fly(end); break;
      case tosemihole: Cut(end); break;
      case tohole:     Fly(end);
                       MakeHole(end);
                       LiftCutter();
                       break;
    }
  }

  // store the number of the current line (cutter cutted)
  int nline2 = CurrentLine();

  for (npass=1; npass<depth->Count(); npass++)
  {
    LiftCutter();
    SetDepth(*(*depth)[npass]);
    if (npass < feed->Count()) SetFeed(*(*feed)[npass]);
    // fly at the head of cache
    Fly(cache->Head());
    DownCutter();
    // repeat the commands in the lines nline1..nline2
    sprintf(BUF, "G14,J1,N1=%d,N2=%d", nline1,nline2); PUT();
    // statuses after the last repeated line (lift the cutter)
    CutterStatus = cut;
    G0G1status   = G1;
  }
}

void MAHO::CutGroove(const EndGroove *groove)
{
  // Do not uses Cut()

  if (groove->status == tofly) return;

  /*
  if (groove->cutype == cylcut) sprintf(BUF, "(BEGIN,OF,CYL,GROOVE)");
  else                          sprintf(BUF, "(BEGIN,OF,BALL,GROOVE)");
  PUT();
  */

  // set the angle (field Dcutter in struct EndGroove)
  double angle = groove->Dcutter;
  // direction of MAHO rotation is opposite to program one
  if (angle != 0.) {
    // prevent rotation to -0.000
    angle = -angle;
  }
  Rotate(angle);

  // depth
  Queue<double> *depths;
  double dcurr;
  // feed
  Queue<long>   *feeds;
  long   fcurr;

  switch (groove->cutype) {
    case cylcut: depths = cdepth;
                 feeds  = cfeed;
                 break;
    case balcut: depths = bdepth;
                 feeds  = bfeed;
                 break;
    default: STOPCR("Wrong groove cutter");
  }

  char RotationStr[32];
  char FeedStr[32];

  Queue<End> *P = groove->P;
  End *p;
  int i;

  // current depth
  switch (groove->cutype) {
    case cylcut: dcurr = *cdepth->Head(); break;
    case balcut: dcurr = *bdepth->Head(); break;
    default: STOPCR("Wrong groove cutter");
  }

  // current feed
  switch (groove->cutype) {
    case cylcut: fcurr = *cfeed->Head(); break;
    case balcut: fcurr = *bfeed->Head(); break;
    default: STOPCR("Wrong groove cutter");
  }

  SetFeed(fcurr);

  // fly in the groove initial point
  i = 0;
  p = (*P)[i];

  // fast move to over the point (with possible changing the cutter)
  Fly(p);

  // Field Dcutter in the 0 point represents height
  double Y0 = p->Dcutter;
  // move slowely to point
  MoveHeight(Y0);

  // vertical move

  i = 1;
  p = (*P)[i];
  if (p->Dcutter > 0.) {
    // clockwise move
    sprintf(BUF, "G2,X%0.3f,Y%0.3f,R%0.3f", p->Y, dcurr, +p->Dcutter);
  }
  else {
    // counterclockwise move
    sprintf(BUF, "G3,X%0.3f,Y%0.3f,R%0.3f", p->Y, dcurr, -p->Dcutter);
  }
  if (NoFeed) {
    sprintf(FeedStr, ",F%ld", CurrentFeed);
    strcat(BUF,FeedStr);
    NoFeed = 0;
  }
  if (NoRotation) {
    sprintf(RotationStr, ",S%ld,M%d",rotation,dir_rotation);
    strcat(BUF,RotationStr);
    NoRotation = 0;
  }
  PUT();
  CutterStatus = cut;
  G0G1status   = G00;

  // store the next line number
  int nline1 = nmain;

  for (i++; i<P->Count(); i++)
  {
    p = (*P)[i];
    if (p->Dcutter == 0.) {
      // cut to this point
      if (G0G1status == G1) sprintf(BUF,    "X%0.3f,Z%0.3f", p->Y,p->X);
      else                  sprintf(BUF, "G1,X%0.3f,Z%0.3f", p->Y,p->X);
      PUT();
      G0G1status = G1;
    }
    else {
      // circular cut
      if (p->Dcutter > 0.) {
        // clockwise
        sprintf(BUF, "G2,X%0.3f,Z%0.3f,R%0.3f", p->Y, p->X, +p->Dcutter);
      }
      else {
        // counterclockwise
        sprintf(BUF, "G3,X%0.3f,Z%0.3f,R%0.3f", p->Y, p->X, -p->Dcutter);
      }
      PUT();
      CutterStatus = cut;
      G0G1status   = G00;
    }
  }
  LiftCutter();

  // last line number
  int nline2 = nmain-LineInc;

  // repeat for each depth
  for (int npass=1; npass<depths->Count(); npass++)
  {
    // current depth
    int idepth = npass;
    switch (groove->cutype) {
      case cylcut: dcurr = *(*cdepth)[idepth]; break;
      case balcut: dcurr = *(*bdepth)[idepth]; break;
      default: STOPCR("Wrong groove cutter");
    }

    // current feed
    int ifeed = npass;
    if (ifeed >= feeds->Count()) ifeed = feeds->Count()-1;
    switch (groove->cutype) {
      case cylcut: fcurr = *(*cfeed)[ifeed]; break;
      case balcut: fcurr = *(*bfeed)[ifeed]; break;
      default: STOPCR("Wrong groove cutter");
    }

    SetFeed(fcurr);

    // fly in the groove initial point
    i = 0;
    p = (*P)[i];
    // Field Dcutter in the 0 point represents height
    double Y0 = p->Dcutter;
    // fast move to over the point
    Fly(p);
    // move slowely to point
    MoveHeight(Y0);

    // vertical move

    i = 1;
    p = (*P)[i];
    if (p->Dcutter > 0.) {
      // clockwise move
      sprintf(BUF, "G2,X%0.3f,Y%0.3f,R%0.3f", p->Y, dcurr, +p->Dcutter);
    }
    else {
      // counterclockwise move
      sprintf(BUF, "G3,X%0.3f,Y%0.3f,R%0.3f", p->Y, dcurr, -p->Dcutter);
    }
    if (NoFeed) {
      sprintf(FeedStr, ",F%ld", CurrentFeed);
      strcat(BUF,FeedStr);
      NoFeed = 0;
    }
    if (NoRotation) {
      sprintf(RotationStr, ",S%ld,M%d",rotation,dir_rotation);
      strcat(BUF,RotationStr);
      NoRotation = 0;
    }
    PUT();

    // repeat the commands in the lines nline1..nline2
    sprintf(BUF, "G14,J1,N1=%d,N2=%d", nline1,nline2); PUT();
    // statuses after the last repeated line
    CutterStatus = cut;
    G0G1status   = G00;
  }

  /*
  if (groove->cutype == cylcut) sprintf(BUF, "(END,OF,CYL,GROOVE)");
  else                          sprintf(BUF, "(END,OF,BALL,GROOVE)");
  PUT();
  */

  // clear CurrentTile
  CurrentTile = NULL;
}

void MAHO::CutGrooves()
{
  if (!ToCut(groove_ends)) return;
  for (int i=0; i<groove_ends->Count(); i++) CutGroove((*groove_ends)[i]);
}

void MAHO::PrintCode()
{
  for (int i=0; i<code->Count(); i++) cout << (*code)[i] << endl;
}

void MAHO::WriteCode()
{
  if (nopath) return;

  ofstream fout(code_file);
  for (int i=0; i<code->Count(); i++) fout << (*code)[i] << endl;
  fout.close();
}
