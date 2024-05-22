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

#if !defined(__CTYPE_H)
#include <ctype.h>        /* for islower */
#endif

#include "queue.h"
#include "elements.h"

#include "sidelist.h"
#include "tile.h"

#include "megatile.h"
#include "maho.h"

#include "graph.h"

#include "wait.h"

#include "textfile.h"

int TestTask(char *taskname);

int ReadTask(char *taskname, char *exename)
{
  TaskText text(taskname,"*",exename);

  char line[256];

  strcpy(MAHO::task_file, taskname);
  text.MoveRvalue("code_file", MAHO::code_file);
  text.MoveRvalue("tile_file", Tile::tile_file);
  *MAHO::path_file = '\0';
  if (text.Occurence("path_file"))   text.MoveRvalue("path_file", MAHO::path_file);
  if (text.Occurence("groove_file")) text.MoveRvalue("groove_file", MAHO::groove_file);

  text.MoveRvalue("SheetZ0", MAHO::SheetX1);
  text.MoveRvalue("SheetX0", MAHO::SheetY1);
  double xSize,ySize;
  text.MoveRvalue("SizeZ",   xSize);
  text.MoveRvalue("SizeX",   ySize);
  if (xSize < 0.) STOPCR("Wrong value of SizeZ");
  if (ySize < 0.) STOPCR("Wrong value of SizeX");
  MAHO::SheetX2 = MAHO::SheetX1 - xSize;
  MAHO::SheetY2 = MAHO::SheetY1 - ySize;

  text.MoveRvalue("G28on",           MAHO::G28on);
  text.MoveRvalue("MAHO_XMIN",       MAHO::MAHO_XMIN);
  text.MoveRvalue("MAHO_XMAX",       MAHO::MAHO_XMAX);
  text.MoveRvalue("MAHO_ZMIN",       MAHO::MAHO_ZMIN);
  text.MoveRvalue("MAHO_ZMAX",       MAHO::MAHO_ZMAX);
  text.MoveRvalue("MAHO_X0",         MAHO::MAHO_X0);
  text.MoveRvalue("MAHO_Y0",         MAHO::MAHO_Y0);
  text.MoveRvalue("MAHO_Z0",         MAHO::MAHO_Z0);
  text.MoveRvalue("MAHO_OX",         MAHO::MAHO_OX);
  text.MoveRvalue("MAHO_OZ",         MAHO::MAHO_OZ);
  text.MoveRvalue("ChangeCutter_Y",  MAHO::ChangeCutter_Y);
  text.MoveRvalue("ChangeCutter_Z",  MAHO::ChangeCutter_Z);
  text.MoveRvalue("zShift",          MAHO::xShift);
  text.MoveRvalue("xShift",          MAHO::yShift);
  text.MoveRvalue("thickness",       MAHO::thickness);
  text.MoveRvalue("minheight",       MAHO::minheight);
  text.MoveRvalue("correct_depth",   MAHO::correct_depth);
  text.MoveRvalue("TrueAngle",       MAHO::TrueAngle);
  text.MoveRvalue("SqueezeParts",    MAHO::SqueezeParts);
  text.MoveRvalue("BackSide",        MAHO::BackSide);
  text.MoveRvalue("Dcut",            MAHO::Dcut);
  text.MoveRvalue("NumberDcut",      MAHO::NumberDcut);
  /*
  text.MoveRvalue("Dhole",           MAHO::Dhole);
  text.MoveRvalue("NumberDhole",     MAHO::NumberDhole);
  */
  text.MoveRvalue("NumberDcyl",      MAHO::NumberDcyl);
  text.MoveRvalue("NumberDball",     MAHO::NumberDbal);
  text.MoveRvalue("P0height",        Tile::P0height);
  text.MoveRvalue("hfly",            MAHO::hfly);
  text.MoveRvalue("htop",            MAHO::htop);
  text.MoveRvalue("rotation",        MAHO::rotation);
  /*
  text.MoveRvalue("dir_rotation",    MAHO::dir_rotation);
  */
  text.MoveRvalue("ProgramNumber",   MAHO::ProgramNumber);
  /*
  text.MoveRvalue("SubNumber",       MAHO::SubNumber);
  */

  // down cut

  if (text.Occurence("down_cut"))
  {
    text.MoveRvalue("down_cut", line);
    char *sub = strtok(line, " ,\t");
    sscanf(sub,"%d",&MAHO::down_cut);
    if (MAHO::down_cut) {
      sub = strtok(NULL, " ,\t");
      if (sub) sscanf(sub,"%d",&MAHO::down_todo);
    }
  }
  else STOPCR("ERROR: Field down_cut did not defined!");

  text.MoveRvalue("ddown",           MAHO::ddown);

  // comments

  if (text.Occurence("comment"))
  {
    // create the MAHO::comment list
    MAHO::comment = new Queue<char>;

    #pragma warn -pia
    while (text.Occurence("comment"))
    {
      text.MoveRvalue("comment", line);
      if (strlen(line) == 0) continue;

      char *s = new char[strlen(line)+1];
      for (int i=0; s[i] = line[i]; i++) {
        if (islower(s[i])) {
          cout << "ERROR in comment line" << endl;
          cout << line << endl;
          STOPCR("Comment should be in UPPER case");
        }
        if (s[i] == '(' || s[i] == ')') {
          STOPCR("Brackets () are not allowable");
        }
        if (s[i] == ' ') s[i]=',';
      }
      if (strlen(s) > 60) {
        cerr << "Truncate comment line:" << endl;
        cerr << s;
        cerr << "to 60 symbols:" << endl;
        s[61] = '\0';
        cerr << s << endl;
        WAIT();
      }
      *MAHO::comment << s;
    }
    #pragma warn +pia
  }

  // gates

  if (text.Occurence("gates"))
  {
    // create the MAHO::gates list
    MAHO::gates = new Queue<Side>;

    #pragma warn -pia
    while (text.Occurence("gates"))
    {
      int np = 0;
      double p[4];

      text.MoveRvalue("gates", line);

      char *sub = strtok(line, " (),;\t");
      if (sub) sscanf(sub,"%lf",&p[np++]);
      while (sub=strtok(NULL, " (),;\t")) {
        sscanf(sub,"%lf",&p[np++]);
        if (np%4 == 0) {
          double x1 = p[1];
          double y1 = p[0];
          double x2 = p[3];
          double y2 = p[2];
          *MAHO::gates << new Side(x1,y1,x2,y2);
          np = 0;
        }
      }
      if (np != 0) STOPCR("Wrong gates");
    }
    #pragma warn +pia
  }
  else STOPCR("Field gates did not defined");

  // depths

  if (text.Occurence("cut_depth"))
  {
    // create the MAHO::depth list
    MAHO::depth = new Queue<double>;

    double the_depth;
    text.MoveRvalue("cut_depth", line);
    char *sub = strtok(line, " ,\t");
    if (sub) {
      sscanf(sub,"%lf",&the_depth);
      *MAHO::depth << new double(the_depth);
    }
    #pragma warn -pia
    while (sub=strtok(NULL, " ,\t")) {
      sscanf(sub,"%lf",&the_depth);
      *MAHO::depth << new double(the_depth);
    }
    #pragma warn +pia
  }
  else STOPCR("ERROR: Field cut_depth did not defined!");

  // cut feeds

  if (text.Occurence("cut_feed"))
  {
    // create the MAHO::feed list
    MAHO::feed = new Queue<long>;

    long the_feed;
    text.MoveRvalue("cut_feed", line);
    char *sub = strtok(line, " ,\t");
    if (sub) {
      sscanf(sub,"%ld",&the_feed);
      *MAHO::feed << new long(the_feed);
    }
    #pragma warn -pia
    while (sub=strtok(NULL, " ,\t")) {
      sscanf(sub,"%ld",&the_feed);
      *MAHO::feed << new long(the_feed);
    }
    #pragma warn +pia
  }
  else STOPCR("ERROR: Field cut_feed did not defined!");

  // cyl.groove depths

  if (text.Occurence("cyl_depth"))
  {
    // create the MAHO::cdepth list
    MAHO::cdepth = new Queue<double>;

    double the_depth;
    text.MoveRvalue("cyl_depth", line);
    char *sub = strtok(line, " ,\t");
    if (sub) {
      sscanf(sub,"%lf",&the_depth);
      *MAHO::cdepth << new double(the_depth);
    }
    #pragma warn -pia
    while (sub=strtok(NULL, " ,\t")) {
      sscanf(sub,"%lf",&the_depth);
      *MAHO::cdepth << new double(the_depth);
    }
    #pragma warn +pia
  }
  else STOPCR("ERROR: Field cyl_depth did not defined!");

  // ball groove depths

  if (text.Occurence("ball_depth"))
  {
    // create the MAHO::bdepth list
    MAHO::bdepth = new Queue<double>;

    double the_depth;
    text.MoveRvalue("ball_depth", line);
    char *sub = strtok(line, " ,\t");
    if (sub) {
      sscanf(sub,"%lf",&the_depth);
      *MAHO::bdepth << new double(the_depth);
    }
    #pragma warn -pia
    while (sub=strtok(NULL, " ,\t")) {
      sscanf(sub,"%lf",&the_depth);
      *MAHO::bdepth << new double(the_depth);
    }
    #pragma warn +pia
  }
  else STOPCR("ERROR: Field ball_depth did not defined!");

  // cyl.groove feeds

  if (text.Occurence("cyl_feed"))
  {
    // create the MAHO::cfeed list
    MAHO::cfeed = new Queue<long>;

    long the_feed;
    text.MoveRvalue("cyl_feed", line);
    char *sub = strtok(line, " ,\t");
    if (sub) {
      sscanf(sub,"%ld",&the_feed);
      *MAHO::cfeed << new long(the_feed);
    }
    #pragma warn -pia
    while (sub=strtok(NULL, " ,\t")) {
      sscanf(sub,"%ld",&the_feed);
      *MAHO::cfeed << new long(the_feed);
    }
    #pragma warn +pia
  }
  else STOPCR("ERROR: Field cyl_feed did not defined!");

  // ball groove feeds

  if (text.Occurence("ball_feed"))
  {
    // create the MAHO::bfeed list
    MAHO::bfeed = new Queue<long>;

    long the_feed;
    text.MoveRvalue("ball_feed", line);
    char *sub = strtok(line, " ,\t");
    if (sub) {
      sscanf(sub,"%ld",&the_feed);
      *MAHO::bfeed << new long(the_feed);
    }
    #pragma warn -pia
    while (sub=strtok(NULL, " ,\t")) {
      sscanf(sub,"%ld",&the_feed);
      *MAHO::bfeed << new long(the_feed);
    }
    #pragma warn +pia
  }
  else STOPCR("ERROR: Field ball_feed did not defined!");

  // input of tiles

  // create the tiles NameList
  if (MAHO::NameList == NULL) MAHO::NameList = new Queue<char>;

  while (text.OccurenceDigit())
  {
    text.MoveLineStartsDigit(line);
    char *item = new char[strlen(line)+1];
    strcpy(item,line);
    // change A --> 1, B --> 2
    char *AB = strstr(item,"-") + 1;
    switch (*AB) {
      case 'A':
      case '1': *AB = '1'; break;
      case 'B':
      case '2': *AB = '2'; break;
      default: cerr << "Wrong tile name "; STOPCR(item);
    }
    *MAHO::NameList << item;
  }

  // input of tile to continue broken previous cutting

  int ncontinue = text.Occurence("@");
  if (ncontinue)
  {
    if (ncontinue > 1) STOPCR("More then one line starts from @");

    text.MoveLineStarts("@", line);

    int i=1;
    switch (line[i]) {
      case '1': Megatile::continue_level = 1; break;
      case '2': Megatile::continue_level = 2; break;
      case '3': Megatile::continue_level = 3; break;
      default: STOPCR("Unknown mark after @");
    }

    // check the space after the @1
    i = 2;
    if (line[i] != ' ') STOPCR("Error in line which starts from @");

    i=3;
    while (line[i] == ' ') i++;

    char *item = new char[strlen(line+i)+1];
    strcpy(item,(line+i));

    // change A --> 1, B --> 2
    char *AB = strstr(item,"-") + 1;
    switch (*AB) {
      case 'A':
      case '1': *AB = '1'; break;
      case 'B':
      case '2': *AB = '2'; break;
      default: cerr << "Wrong tile name "; STOPCR(item);
    }

    strcpy(MAHO::continue_tile,item);
    char *blank = strstr(MAHO::continue_tile, " ");
    if (blank) *blank = '\0';

    // Note: NameList already created
    *MAHO::NameList << item;
  }

  //cout << "ReadTask: Tiles to arrange:" << endl;
  //for (int k=0; k<MAHO::NameList->Count(); k++) cout << (*MAHO::NameList)[k] << endl;
  //cout << "----------------" << endl;

  //cout << "depths:" << endl;
  //for (int idepth=0; idepth<MAHO::depth->Count(); idepth++) {
  //  cout << *(*MAHO::depth)[idepth] << endl;
  //}
  //cout << "--------" << endl;

  // colors
  text.MoveRvalue("Color",           GManager::Color);
  text.MoveRvalue("BkColor",         GManager::BkColor);
  text.MoveRvalue("TextColor",       GManager::TextColor);
  text.MoveRvalue("MessageColor",    GManager::MessageColor);
  text.MoveRvalue("CutColor",        GManager::CutColor);
  text.MoveRvalue("FlyColor",        GManager::FlyColor);
  text.MoveRvalue("ContourColor",    GManager::ContourColor);
  text.MoveRvalue("SideColor",       GManager::SideColor);
  text.MoveRvalue("CylGrooveColor",  GManager::CylGrooveColor);
  text.MoveRvalue("BallGrooveColor", GManager::BallGrooveColor);
  text.MoveRvalue("GatesColor",      GManager::GatesColor);
  text.MoveRvalue("SheetColor",      GManager::SheetColor);
  text.MoveRvalue("WorkAreaColor",   GManager::WorkAreaColor);
  text.MoveRvalue("FinishColor",     GManager::FinishColor);

  GManager::CuttingDelay = 2;  // in milliseconds
  if (text.Occurence("CuttingDelay")) {
    text.MoveRvalue("CuttingDelay",  GManager::CuttingDelay);
  }
  GManager::FlyingDelay = 200;  // in milliseconds
  if (text.Occurence("FlyingDelay")) {
    text.MoveRvalue("FlyingDelay",   GManager::FlyingDelay);
  }
  GManager::GroovingDelay = 500;  // in milliseconds
  if (text.Occurence("GroovingDelay")) {
    text.MoveRvalue("GroovingDelay", GManager::GroovingDelay);
  }
  if (text.Occurence("ShowTurnedContour")) {
    text.MoveRvalue("ShowTurnedContour", MAHO::ShowTurnedContour);
  }
  if (text.Occurence("FinishSnow")) {
    text.MoveRvalue("FinishSnow", GManager::FinishSnow);
  }

  /*
  // check the possible incorrect usage of symbols instead of @
  char *wrongline;
  #pragma warn -pia
  if (wrongline = text.LineStarts("&")) {
    cerr << "Unknown mark & in the line" << endl;
    STOPCR(wrongline);
  }
  if (wrongline = text.LineStarts("$")) {
    cerr << "Unknown mark $ in the line" << endl;
    STOPCR(wrongline);
  }
  #pragma warn +pia
  */

  if (text.Nlines()) {
    cout << endl;
    cout << "Unknown lines in task file " << taskname << ":" << endl;
    text.Show();
    return 0;
  }

  if (!TestTask(taskname)) return 0;

  return 1;
}

int TestTask(char *taskname)
{
  cout << "Test of parameters of task file " << taskname << endl << endl;

  FILE *f;

	f = fopen(Tile::tile_file, "r");
  if (f == NULL) {
    cerr << "Tile parameters file " << Tile::tile_file
         << " did not found." << endl;
    STOPCR();
  }
  else fclose(f);

	f = fopen(MAHO::code_file, "r");
  if (f != NULL) {
    fclose(f);
    cerr << "MAHO code file " << MAHO::code_file
         << " already exist"    << endl;
    STOPCR("Remove it to run the program");
    //WAIT("Remove it to run the program");     // for debugging only
  }
  if (MAHO::G28on           < 0)    STOPCR("Wrong value of G28on");
  if (MAHO::correct_depth   < 0)    STOPCR("Wrong value of correct_depth");
  if (fabs(MAHO::MAHO_XMIN) > 1.e5) STOPCR("Wrong value of MAHO_XMIN");
  if (fabs(MAHO::MAHO_XMAX) > 1.e5) STOPCR("Wrong value of MAHO_XMAX");
  if (fabs(MAHO::MAHO_ZMIN) > 1.e5) STOPCR("Wrong value of MAHO_ZMIN");
  if (fabs(MAHO::MAHO_ZMAX) > 1.e5) STOPCR("Wrong value of MAHO_ZMAX");
  if (fabs(MAHO::MAHO_X0)   > 1.e5) STOPCR("Wrong value of MAHO_X0");
  if (fabs(MAHO::MAHO_Y0)   > 1.e5) STOPCR("Wrong value of MAHO_Y0");
  if (fabs(MAHO::MAHO_Z0)   > 1.e5) STOPCR("Wrong value of MAHO_Z0");
  if (fabs(MAHO::MAHO_OX)   > 1.e5) STOPCR("Wrong value of MAHO_OX");
  if (fabs(MAHO::MAHO_OZ)   > 1.e5) STOPCR("Wrong value of MAHO_OZ");
  if (fabs(MAHO::ChangeCutter_Y) > 1.e5) STOPCR("Wrong value of ChangeCutter_Y");
  if (fabs(MAHO::ChangeCutter_Z) > 1.e5) STOPCR("Wrong value of ChangeCutter_Z");
  if (fabs(MAHO::SheetX1)   > 1.e5) STOPCR("Wrong value of SheetZ0");
  if (fabs(MAHO::SheetY1)   > 1.e5) STOPCR("Wrong value of SheetX0");
  if (     MAHO::xShift     < 0.)   STOPCR("Wrong value of zShift");
  if (     MAHO::yShift     < 0.)   STOPCR("Wrong value of xShift");
  if (fabs(MAHO::thickness) < 0.)   STOPCR("Wrong value of thickness");
  if (fabs(MAHO::minheight) < 0.)   STOPCR("Wrong value of minheight");
  if (MAHO::down_cut        < 0)    STOPCR("Wrong value of down_cut");
  if (MAHO::ddown           < 0.)   STOPCR("Wrong value of ddown");
  if (MAHO::TrueAngle       < 0)    STOPCR("Wrong value of TrueAngle");
  if (MAHO::SqueezeParts    < 0)    STOPCR("Wrong value of SqueezeParts");
  if (MAHO::BackSide        < 0)    STOPCR("Wrong value of BackSide");
  if (fabs(MAHO::Dcut)      > 1.e5) STOPCR("Wrong value of Dcut");
  if (MAHO::NumberDcut      > 100)  STOPCR("Wrong value of NumberDcut");
  /*
  if (fabs(MAHO::Dhole)     > 1.e5) STOPCR("Wrong value of Dhole");
  if (MAHO::NumberDhole     > 100)  STOPCR("Wrong value of NumberDhole");
  */
  if (MAHO::NumberDcyl      > 100)  STOPCR("Wrong value of NumberDcyl");
  if (MAHO::NumberDbal      > 100)  STOPCR("Wrong value of NumberDball");
  if (fabs(Tile::P0height)  > 1.e5) STOPCR("Wrong value of P0height");
  if (fabs(MAHO::hfly)      > 1.e5) STOPCR("Wrong value of hfly");
  if (fabs(MAHO::htop)      > 1.e5) STOPCR("Wrong value of htop");
  if (MAHO::rotation > 1000000)     STOPCR("Wrong value of rotation");
  /*
  if (MAHO::dir_rotation    > 100)  STOPCR("Wrong value of dir_rotation");
  */
  if (MAHO::ProgramNumber < 9000 || MAHO::ProgramNumber > 9999999)
     STOPCR("Wrong value of ProgramNumber");
  /*
  if (MAHO::SubNumber < 9000 || MAHO::SubNumber > 9999999)
     STOPCR("Wrong value of ProgramNumber");
  if (MAHO::SubNumber == MAHO::ProgramNumber)
     STOPCR("Error: SubNumber == ProgramNumber");
  */

  // depths
  int n;
  for (n=0; n<MAHO::depth->Count(); n++) if (fabs(*(*MAHO::depth)[n])>1.e3) {
    cerr << "Wrong value of depth " << *(*MAHO::depth)[n];
    STOPCR();
  }

  // down cut / continue broken work
  if (strlen(MAHO::continue_tile))
  {
    /*
    if (MAHO::down_cut) {
      WAIT("Continue work and down cut are requested together.");
    }
    */

    char nameAB[16];
    strcpy(nameAB, MAHO::continue_tile);
    // change 1 --> A, 2 --> B
    char *AB = strstr(nameAB,"-") + 1;
    switch (*AB) {
      case 'A':
      case '1': *AB = 'A'; break;
      case 'B':
      case '2': *AB = 'B'; break;
      default: cerr << "Wrong tile name "; STOPCR(nameAB);
    }
    cout << endl;
    cout << "Be careful!" << endl;
    cout << "-----------" << endl;
    cout << "Are you sure to continue broken ";
    switch (Megatile::continue_level) {
      case 1: cout << "cyl. groove ";
      case 2: cout << "ball groove ";
      case 3: cout << "cutting ";
    }
    cout << "from tile " << nameAB << "?" << endl << endl;
    WAIT();
    cout << "Confirm start from tile " << nameAB << endl;
    WAIT();
    WAIT("Recomendation: Remove symbol @ right away this run!");
  }
  else {
    if (MAHO::down_cut) {
      if (MAHO::down_todo == 0) {
        // down cut is already done
        cout << "The down cut is marked as already done. Is it right?" << endl;
        WAIT("Confirm that down cut is already done");
      }
      else switch (MAHO::down_cut) {
        case  1: WAIT("Are you sure to perform down cut of sheet?"); break;
        default: WAIT("Are you sure to perform cut of sheet contour?");
      }
    }
  }

  if (MAHO::BackSide)
  {
    WAIT("Are you sure to process back side of sheet?");

    if (MAHO::thickness < 7.) {
      cout << endl << "Current thickness now is only "
           << MAHO::thickness << endl;
      WAIT("Are you sure to process back side of sheet?");
    }

    if (MAHO::down_todo) {
      WAIT("Confirm the making of down cut with back side of sheet");
    }
  }

  // safety gates
  n = MAHO::gates->Count();
  if (n < 1) STOPCR("Field gates did not defined");
  /*
  cout << "TestTask: Gates:" << endl;
  int np = MAHO::gates->Count();
  for (int ns=0; ns<np; ns++) {
    Side *side = (*MAHO::gates)[ns];
    side->Show();
    cout << endl;
  }
  cout << endl;
  */

  // check tiles
  int ntiles = MAHO::NameList->Count();
  if (ntiles) {
    char *nline = new char[80];
    char *iline = new char[80];
    for (int n=0; n<ntiles; n++) {
      strcpy(nline, (*MAHO::NameList)[n]);
      char *nname = strtok(nline, " ,;\t");
      for (int i=n+1; i<ntiles; i++) {
        strcpy(iline, (*MAHO::NameList)[i]);
        char *iname = strtok(iline, " ,;\t");
        if (strcmp(nname, iname) == 0) {
          cerr << "Tile name " << nname << " occurs more then once" << endl;
          STOPCR();
        }
      }
    }
    delete nline;
    delete iline;
  }

  // colors
  if (GManager::Color < 0 || GManager::Color > 15) {
    WAIT("Wrong value of Color. Assume default.");
    GManager::Color = 1;
  }
  if (GManager::BkColor < 0 || GManager::BkColor > 15) {
    WAIT("Wrong value of BkColor. Assume default.");
    GManager::BkColor = 1;
  }
  if (GManager::TextColor < 0 || GManager::TextColor > 15) {
    WAIT("Wrong value of TextColor. Assume default.");
    GManager::TextColor = 1;
  }
  if (GManager::MessageColor < 0 || GManager::MessageColor > 15) {
    WAIT("Wrong value of MessageColor. Assume default.");
    GManager::MessageColor = 1;
  }
  if (GManager::CutColor < 0 || GManager::CutColor > 15) {
    WAIT("Wrong value of CutColor. Assume default.");
    GManager::CutColor = 1;
  }
  if (GManager::FlyColor < 0 || GManager::FlyColor > 15) {
    WAIT("Wrong value of FlyColor. Assume default.");
    GManager::FlyColor = 1;
  }
  if (GManager::ContourColor < 0 || GManager::ContourColor > 15) {
    WAIT("Wrong value of ContourColor. Assume default.");
    GManager::ContourColor = 1;
  }
  if (GManager::SideColor < 0 || GManager::SideColor > 15) {
    WAIT("Wrong value of SideColor. Assume default.");
    GManager::SideColor = 1;
  }
  if (GManager::CylGrooveColor < 0 || GManager::CylGrooveColor > 15) {
    WAIT("Wrong value of CylGrooveColor. Assume default.");
    GManager::CylGrooveColor = 1;
  }
  if (GManager::BallGrooveColor < 0 || GManager::BallGrooveColor > 15) {
    WAIT("Wrong value of BallGrooveColor. Assume default.");
    GManager::BallGrooveColor = 1;
  }
  if (GManager::GatesColor < 0 || GManager::GatesColor > 15) {
    WAIT("Wrong value of GatesColor. Assume default.");
    GManager::GatesColor = 1;
  }
  if (GManager::SheetColor < 0 || GManager::SheetColor > 15) {
    WAIT("Wrong value of SheetColor. Assume default.");
    GManager::SheetColor = 1;
  }
  if (GManager::WorkAreaColor < 0 || GManager::WorkAreaColor > 15) {
    WAIT("Wrong value of WorkAreaColor. Assume default.");
    GManager::WorkAreaColor = 1;
  }
  if (GManager::FinishColor < 0 || GManager::FinishColor > 15) {
    WAIT("Wrong value of FinishColor. Assume default.");
    GManager::FinishColor = 1;
  }

  cout << "Test is OK" << endl;
  return 1;
}
