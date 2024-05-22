// Andriy_Zatserklyaniy@fnal.gov

#if !defined(__STDIO_H)
#include <stdio.h>
#endif

#if !defined(__MATH_H)
#include <math.h>
#endif

#if !defined(__STRING_H)
#include <string.h>
#endif

#if !defined(__IOSTREAM_H)
#include <iostream.h>
#endif

#include "queue.h"
#include "sidelist.h"
#include "elements.h"

#include "wait.h"

#include "tile.h"

// OddLow

int Tile::OddLow(char* line) {
  char h;
  int items = 0;
  int res = 0;
  int nsholes;

  double A,B,C,S,H;
  //A,B,C,S,H
  if (strstr(line+strlen(name),"-")) {
    nsholes = 1;
    items = sscanf(line+strlen(name), "%lf %lf %lf %lf %c",
                                        &A, &B, &C, &S, &h);
  }
  else {
    nsholes = 2;
    items = sscanf(line+strlen(name), "%lf %lf %lf %lf %lf",
                                        &A, &B, &C, &S, &H);
  }

  if (items != 5) {
    cerr << endl << "Error reading line for tile " << nameAB << endl;
    STOPCR();
    res = 1;
  }

  // correction for tolerance
  A += ToleranceShift();
  B += ToleranceShift();
  C += ToleranceShift();
//
// Important note: tiles from megaside=2 already inverted in this way,
// ------------------------------------------------------------------
//           namely, inverted with respect to C side.
//
  double O=0.;
  Side *side1 = new Side(O,O, O,C);
  Side *side2 = new Side(O,C, A,C);
  Side *side3 = new Side(A,C, B,O);
  Side *side4 = new Side(B,O, O,O);

  *side2 << new SHole(S,C,Dsemihole);

  switch (nsholes) {
    case 1:  break;
    case 2:  *side4 << new SHole(H,O,Dsemihole); break;
    default: cout << "Tile " << name; STOP(": wrong number of holes");
  }

  *contour << side1 << side2 << side3 << side4;

  //MakeInner();

  if (megaside==2) {
    // reinvert tile
    //cout << "Before invert" << endl; contour->Show();
    for(int i=0; i<contour->Count(); i++) (*contour)[i]->InvX();
    //cout << "After  invert" << endl; contour->Show();
    // Correct the direction going round
    //contour->Reverse();
    //contour->SetHead(1);
    //cout << "After  SideList's Reverse and SetHead" << endl; contour->Show();
  }

  // start point situates at the side 3
  int ns_start = 2;
  ClockwiseGroove(ns_start);

  /*
  cout << "Tile " << name << "   A,B,C,S,H: "
       << A << " " << B << " " << C << " " << S << " ";
  if (H) cout << H;
  else   cout << "-";
  cout << endl;
  */
  return res;
}

// OddHigh

int Tile::OddHigh(char* line) {
  char v;
  int items = 0;
  int res = 0;
  int nsholes;

  double A,B,C,R,V;
  //A,B,C,R,V
  if (strstr(line+strlen(name),"-")) {
    nsholes = 1;
    items = sscanf(line+strlen(name), "%lf %lf %lf %lf %c",
                                        &A, &B, &C, &R, &v);
  }
  else {
    nsholes = 2;
    items = sscanf(line+strlen(name), "%lf %lf %lf %lf %lf",
                                        &A, &B, &C, &R, &V);
  }

  if (items != 5) {
    cerr << endl << "Error reading line for tile " << nameAB << endl;
    STOPCR();
    res = 1;
  }

  // correction for tolerance
  A += ToleranceShift();
  B += ToleranceShift();
  C += ToleranceShift();
//
// Important note: tiles from megaside=2 already inverted in this way,
// ------------------------------------------------------------------
//           namely, inverted with respect to C side.
//
  double O=0.;
  Side *side1 = new Side(O,O, O,C);
  Side *side2 = new Side(O,C, A,C);
  Side *side3 = new Side(A,C, B,O);
  Side *side4 = new Side(B,O, O,O);

  *side4 << new SHole(R,O,Dsemihole);

  switch (nsholes) {
    case 1:  break;
    case 2:  *side2 << new SHole(V,C,Dsemihole); break;
    default: cout << "Tile " << name; STOP(": wrong number of holes");
  }

  *contour << side1 << side2 << side3 << side4;

  if (megaside==2) {
    // reinvert tile
    for(int i=0; i<contour->Count(); i++) (*contour)[i]->InvX();
    // Correct the direction going round
    //contour->Reverse();
    //contour->SetHead(1);
  }

  //  cout << " the number of holes: " <<
  /*
    cout << "Tile " << name << "   A,B,C,R,V: "
         << A << " " << B << " " << C << " " << R << " ";
    if (V) cout << V;
    else   cout << "-";
    cout << endl;
  */

  // start point situates at the side 0
  int ns_start = 0;
  CounterClockwiseGroove(ns_start);
  return res;
}

// EvenLow

int Tile::EvenLow(char* line) {
  char s,t,w;
  int items = 0;
  int res = 0;
  char *substr;
  int nsholes;

  substr = strstr(line+strlen(name),"-");

  double A,B,C,H,S,T,W;
  //A,B,C,H,S,T,W
  if (strstr(substr+1,"-")) {
    nsholes = 2;
    items = sscanf(line+strlen(name), "%lf %lf %lf %lf %lf  %c  %c",
                                        &A, &B, &C, &H, &S, &t, &w);
  }
  else {
    nsholes = 3;
    items = sscanf(line+strlen(name), "%lf %lf %lf %lf %c  %lf %lf",
                                        &A, &B, &C, &H, &s, &T, &W);
  }

  if (items != 7) {
    cerr << endl << "Error reading line for tile " << nameAB << endl;
    STOPCR();
    res = 1;
  }

  // correction for tolerance
  A += ToleranceShift();
  B += ToleranceShift();
  C += ToleranceShift();
//
// Important note: tiles from megaside=2 already inverted in this way,
// ------------------------------------------------------------------
//           namely, inverted with respect to C side.
//
  double O=0.;
  Side *side1 = new Side(O,O, O,C);
  Side *side2 = new Side(O,C, A,C);
  Side *side3 = new Side(A,C, B,O);
  Side *side4 = new Side(B,O, O,O);

  *side4 << new SHole(H,O,Dsemihole);

  switch (nsholes) {
    case 2:  *side2 << new SHole(S,C,Dsemihole); break;
    case 3:  *side2 << new SHole(T,C,Dsemihole) << new SHole(W,C,Dsemihole); break;
    default: cout << "Tile " << name; STOP(": wrong number of holes");
  }

  *contour << side1 << side2 << side3 << side4;

  if (megaside==2) {
    // reinvert tile
    for(int i=0; i<contour->Count(); i++) (*contour)[i]->InvX();
    // Correct the direction going round
    //contour->Reverse();
    //contour->SetHead(1);
  }

  /*
    cout << "Tile " << name << "   A,B,C,H,S,T,W: "
         << A << " " << B << " " << C << " " << H;
    if (S) cout << " " << S;   else cout << " " << "-";
    if (T) cout << " " << T;   else cout << " " << "-";
    if (W) cout << " " << W;   else cout << " " << "-";
    cout << endl;
  */

  // start point situates at the side 0
  int ns_start = 0;
  CounterClockwiseGroove(ns_start);
  return res;
}

// EvenHigh

int Tile::EvenHigh(char* line) {
  char u;
  int items = 0;
  int res = 0;
  int nsholes;

  double E,D,C,F,Q,U;
  // E,D,C,F,Q,U;
  if (strstr(line+strlen(name),"-")) {
    nsholes = 1;
    items = sscanf(line+strlen(name), "%lf %lf %lf %lf %lf %c",
                                        &E, &D, &C, &F, &Q, &u);
  }
  else {
    nsholes = 2;
    items = sscanf(line+strlen(name), "%lf %lf %lf %lf %lf %lf",
                                        &E, &D, &C, &F, &Q, &U);
  }

  if (items != 6) {
    cerr << endl << "Error reading line for tile " << nameAB << endl;
    STOPCR();
    res = 1;
  }

  // correction for tolerance
  F += ToleranceShift();
  D += ToleranceShift();
  E += ToleranceShift();
  C += ToleranceShift();
//
// Important note: tiles from megaside=2 already inverted in this way,
// ------------------------------------------------------------------
//           namely, inverted with respect to C side.
//
  double O=0.;
  Side *side1 = new Side(O,  O,  F-D,C);
  Side *side2 = new Side(F-D,C,  F,  C);
  Side *side3 = new Side(F,  C,  E,  O);
  Side *side4 = new Side(E,  O,  O,  O);

  *side4 << new SHole(Q,O,Dsemihole);

  switch (nsholes) {
    case 1:  break;
    case 2:  *side2 << new SHole(F-D+U,C,Dsemihole); break;
    default: cout << "Tile " << name; STOP(": wrong number of holes");
  }

  *contour << side1 << side2 << side3 << side4;

  if (megaside==2) {
    // reinvert tile
    for(int i=0; i<contour->Count(); i++) (*contour)[i]->InvX();
    // Correct the direction going round
    //contour->Reverse();
    //contour->SetHead(1);
  }

  /*
    cout << "Tile " << name << "   E,D,C,F,Q,U: "
         << E << " " << D << " " << C << " " << F << " " << Q << " ";
    if (U) cout << U;   else cout << "-";
    cout << endl;
  */

  // start point situates at the side 3
  int ns_start = 2;
  ClockwiseGroove(ns_start);
  return res;
}

void Tile::MakeTile(char *tname) {
  sscanf(tname, "%d%*1s%1d%*1s%d", &megatile,&megaside,&number);
  int len = strlen(tname);

  if (len < 5) {
    cerr << "Wrong tile name " << tname << endl;
    exit(1);
  }

  char tnameAB[16];
  strcpy(tnameAB, tname);

  if (OldDrawings == 0)
  {
    // change 1 --> A, 2 --> B
    char *AB = strstr(tnameAB,"-") + 1;
    switch (*AB) {
      case 'A':
      case '1': *AB = 'A'; break;
      case 'B':
      case '2': *AB = 'B'; break;
      default: cerr << "Wrong tile name "; STOPCR(nameAB);
    }
  }

  // Read the file to find corresponding line

  FILE *f;
	f = fopen(tile_file, "r");
  if (f == NULL) {
    cerr << "\nTile parameters file " << tile_file << " did not found.\n"
         << endl;
    STOPCR();
  }
  char line[256], token_line[256], *match;
  while (fgets(line,256,f)) {
    if (line[0] == '#' || line[0] == '*') continue;

    // look the token tname in line
    strcpy(token_line,line);
    match = strtok(token_line, " \t");

    if (strcmp(match,tnameAB) == 0) goto FOUND;
  }
  cerr << "Tile " << tnameAB << " did not found" << endl;
  fclose(f);
  STOPCR();

  FOUND:
  fclose(f);

  if (number%2 == 0) {
    // the number is even
    if (number < MinHigh) EvenLow(strstr(line,match));
    else                  EvenHigh(strstr(line,match));
  }

  else {
    // the number is odd
    if (number < MinHigh) OddLow(strstr(line,match));
    else                  OddHigh(strstr(line,match));
  }
}
