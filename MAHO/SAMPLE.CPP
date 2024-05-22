// Andriy_Zatserklyaniy@fnal.gov

#if !defined(__IOSTREAM_H)
#include <iostream.h>
#endif

#if !defined(__FSTREAM_H)
#include <fstream.h>
#endif

#if !defined(__STDLIB_H)
#include <stdlib.h>
#endif

#if !defined(__STDIO_H)
#include <stdio.h>
#endif

#include "wait.h"

void sample(char *file, char *exename)
{
FILE *f;
f = fopen(file, "r");
if (f != NULL) {
  fclose(f);
  cerr << "\nMAHO task file " << file << " already exist" << endl;
  STOPCR("Choose another file name");
}

ofstream fout(file);

fout<<"*"<<endl;
fout<<"* Sample of task file for " << exename <<endl;
fout<<"* The order of the lines in this file is not substantial"<<endl;
fout<<"* Line that starts from symbol * is treated as a commentary"<<endl;
fout<<"*"<<endl;
fout<<"* MAHO code file"<<endl;
fout<<"code_file   = sample.cod"<<endl;
fout<<"* tiles parameters file"<<endl;
fout<<"tile_file   = sample.par"<<endl;
fout<<""<<endl;
fout<<"* path file name (optional)"<<endl;
fout<<"path_file   = sample.way"<<endl;
fout<<"* groove file name (optional)"<<endl;
fout<<"groove_file = sample.gro"<<endl;
fout<<""<<endl;
fout<<"* No. of program (9000 <= ProgramNumber <= 9999999)"<<endl;
fout<<"ProgramNumber = 10000"<<endl;
fout<<""<<endl;
fout<<"* Comment line (should be in UPPER case, no brackets (), max 60 symbols)"<<endl;
fout<<"* Note: keyword comment can be used as many times as you need"<<endl;
fout<<"comment = SOME COMMENT"<<endl;
fout<<""<<endl;
/*
fout<<"* No. of subroutine (9000 <= SubNumber <= 9999999) "
      "and != ProgramNumber)"<<endl;
fout<<"SubNumber = 10001"<<endl;
fout<<""<<endl;
*/
fout<<"*"<<endl;
fout<<"*                    SizeZ                          "<<endl;
fout<<"*          +------------------------+                 "<<endl;
fout<<"*                                    (SheetX0,SheetZ0)"<<endl;
fout<<"*          +------------------------*   +             "<<endl;
fout<<"*        | |                        |   |             "<<endl;
fout<<"* xShift | |                        |   |             "<<endl;
fout<<"*        | |      ___________       |   |             "<<endl;
fout<<"*          |     |       /   |      |   |             "<<endl;
fout<<"*          |     |      /    |      |   |             "<<endl;
fout<<"*          |     |     /     |      |   | SizeX       "<<endl;
fout<<"*          |     |    / tile |      |   |             "<<endl;
fout<<"*          |     |   /       |      |   |             "<<endl;
fout<<"*          |     |  /________|      |   |             "<<endl;
fout<<"*          |     |           |      |   |             "<<endl;
fout<<"*          |     |           |      |   |             "<<endl;
fout<<"*          |     +-----------+      |   |             "<<endl;
fout<<"*        | |                        |   |             "<<endl;
fout<<"* xShift | |                        |   |             "<<endl;
fout<<"*        | |                        |   |             "<<endl;
fout<<"*          +-----+-----------+------+   +             "<<endl;
fout<<"*          ------             ------                  "<<endl;
fout<<"*          zShift             zShift                  "<<endl;
fout<<""<<endl;
fout<<"* Function G28: 0 = off, >0 = on"<<endl;
fout<<"G28on = 0"<<endl;
fout<<""<<endl;
fout<<"* Machine work area in machine frame"<<endl;
fout<<"MAHO_XMIN = 0."<<endl;
fout<<"MAHO_XMAX = 1000."<<endl;
fout<<"MAHO_ZMIN = 0."<<endl;
fout<<"MAHO_ZMAX = 800."<<endl;
fout<<""<<endl;
fout<<"* the center of rotation in machine frame"<<endl;
fout<<"MAHO_OX = 500.099"<<endl;
fout<<"MAHO_OZ = 382.532"<<endl;
fout<<""<<endl;
fout<<"* WorkPiece zero (X,Z) in machine frame"<<endl;
fout<<"MAHO_X0 = 983.07"<<endl;
fout<<"MAHO_Z0 = 630.01"<<endl;
fout<<"* Table zero (Y) in machine frame"<<endl;
fout<<"MAHO_Y0 = 101.98"<<endl;
fout<<""<<endl;
fout<<"* Point to change the cutter in program frame"<<endl;
fout<<"ChangeCutter_Y =  200."<<endl;
fout<<"ChangeCutter_Z = -600."<<endl;
fout<<""<<endl;
fout<<"* Coordinates of sheet corner in program frame"<<endl;
fout<<"SheetZ0 = -20."<<endl;
fout<<"SheetX0 = -20."<<endl;
fout<<"* The lengths of sheet sides"<<endl;
fout<<"SizeX = 600."<<endl;
fout<<"SizeZ = 550."<<endl;
fout<<""<<endl;
fout<<"* small shift from sheet edges. Should be positive."<<endl;
fout<<"xShift  = 5."<<endl;
fout<<"zShift  = 5."<<endl;
fout<<""<<endl;
fout<<"* safety gates: (x1,z1, x2,z2) in program frame"<<endl;
fout<<"* Note: keyword gates can be used as many times as you need"<<endl;
fout<<"gates = (0., 0., 0., -400) (0. 0.,-50., 0.)"<<endl;
fout<<"gates = (-400., 0., -500., 0.)"<<endl;
fout<<""<<endl;
fout<<"* down cut of sheet: 0 = not perform, 1 = down cut, >1 = cut sheet contour"<<endl;
fout<<"* symbol 0 after the value of down_cut (like down_cut = 2 0)"<<endl;
fout<<"* means that done cut is already done"<<endl;
fout<<"down_cut = 2"<<endl;
fout<<"* extra space for down cut. Should be positive."<<endl;
fout<<"ddown = 5."<<endl;
fout<<""<<endl;
fout<<"* Squeeze of left and right sets of tiles: 0 = not perform, 1 = perform"<<endl;
/*
fout<<"* Note: Commonly squeezing of parts inconsitent with BackSide"<< endl;
*/
fout<<"SqueezeParts = 1"<<endl;
fout<<""<<endl;
fout<<"* Back side sign: 0 = normal side, >0 = back side"<<endl;
fout<<"BackSide = 0"<<endl;
fout<<""<<endl;
fout<<"* Cutter diameter for cut and semihole"<<endl;
fout<<"Dcut  = 3.175"<<endl;
fout<<"* The number of cut cutter"<<endl;
fout<<"NumberDcut = 12"<<endl;
fout<<""<<endl;
/*
fout<<"* Hole diameter"<<endl;
fout<<"Dhole = 3."<<endl;
fout<<"* The number of hole cutter"<<endl;
fout<<"NumberDhole = 0"<<endl;
fout<<""<<endl;
*/
fout<<"* The number of cyl. groove cutter"<<endl;
fout<<"NumberDcyl  = 10"<<endl;
fout<<"* The number of ball groove cutter"<<endl;
fout<<"NumberDball = 11"<<endl;
fout<<""<<endl;
fout<<"* height of cutter over groove initial point P0"<<endl;
fout<<"P0height = 0.5"<<endl;
fout<<""<<endl;
fout<<"* height of fly"<<endl;
fout<<"hfly = 6."<<endl;
fout<<""<<endl;
fout<<"* the height of cutter top position"<<endl;
fout<<"htop = 200."<<endl;
fout<<""<<endl;
fout<<"* rotation"<<endl;
fout<<"rotation = 2000"<<endl;
fout<<""<<endl;
/*
fout<<"* direction of rotation"<<endl;
fout<<"dir_rotation = 3"<<endl;
fout<<""<<endl;
*/
fout<<"* Turn table to tile groove angle:"<<endl;
fout<<"*  0: turn exactly +/- 10 degrees"<<endl;
fout<<"* >0: turn to the real angle for every tile"<<endl;
fout<<"TrueAngle = 1"<<endl;
fout<<""<<endl;
fout<<"* the total thickness (with gaskets) > 0."<<endl;
fout<<"thickness = 4.0"<<endl;
fout<<""<<endl;
fout<<"* min cutter height over table > 0."<<endl;
fout<<"minheight = 0.15"<<endl;
fout<<""<<endl;
fout<<"* correct the last cut depth for minheight-thickness: 1, no: 0"<<endl;
fout<<"correct_depth = 1"<<endl;
fout<<""<<endl;
fout<<"* the depths of the cut"<<endl;
fout<<"cut_depth = -2.0 -3.9"<<endl;
fout<<"* the feeds of cut"<<endl;
fout<<"cut_feed  = 300"<<endl;
fout<<""<<endl;
fout<<"* the depths of the cyl. groove passes"<<endl;
fout<<"cyl_depth  = -1. -1.78"<<endl;
fout<<"* the feeds of the cyl. groove passes"<<endl;
fout<<"cyl_feed   = 200 201"<<endl;
fout<<""<<endl;
fout<<"* the depths of the ball groove passes"<<endl;
fout<<"ball_depth = -1.78 -1.78"<<endl;
fout<<"* the feeds of the ball groove passes"<<endl;
fout<<"ball_feed  = 202   252"<<endl;
fout<<""<<endl;
fout<<"* Graphics options: colors"<<endl;
fout<<"Color           = 0"<<endl;
fout<<"BkColor         = 7"<<endl;
fout<<"TextColor       = 1"<<endl;
fout<<"MessageColor    = 4"<<endl;
fout<<"CutColor        = 4"<<endl;
fout<<"FlyColor        = 9"<<endl;
fout<<"ContourColor    = 15"<<endl;
fout<<"SideColor       = 4"<<endl;
fout<<"CylGrooveColor  = 2"<<endl;
fout<<"BallGrooveColor = 6"<<endl;
fout<<"GatesColor      = 12"<<endl;
fout<<"SheetColor      = 8"<<endl;
fout<<"WorkAreaColor   = 1"<<endl;
fout<<"FinishColor     = 14"<<endl;
fout<<"* delay times for graphics (in milliseconds)"<<endl;
fout<<"CuttingDelay  = 2"<<endl;
fout<<"FlyingDelay   = 200"<<endl;
fout<<"GroovingDelay = 500"<<endl;
fout<<""<<endl;
fout<<"* Show finish show (optional)"<<endl;
fout<<"FinishSnow = 0"<<endl;
fout<<""<<endl;
fout<<"* Show turned contour (optional)"<<endl;
fout<<"ShowTurnedContour = 1"<<endl;
fout<<""<<endl;
fout<<"*"<<endl;
fout<<"* Tiles to cut"<<endl;
fout<<"* Both forms (e.g. 2-1-6 and 2-A-6) are acceptable."<<endl;
fout<<"* Symbol 0 after the tile name (like 2-A-6 0) means no groove on this tile"<<endl;
fout<<"* Symbol @ (like @3 2-A-6) marks tile to continue broken work"<<endl;
fout<<"* @1 2-A-6 means continue broken cyl. groove from tile 2-A-6,"<<endl;
fout<<"* @2 2-A-6   -   continue broken ball groove from tile 2-A-6,"<<endl;
fout<<"* @3 2-A-6   -   continue broken cutting from tile 2-A-6,"<<endl;
fout<<"* Recomendation: If you use symbol @ in current run"<<endl;
fout<<"*                remove it right away the run"<<endl;
fout<<"*"<<endl;
fout<<""<<endl;
fout<<"2-A-1"<<endl;
fout<<"2-A-2"<<endl;
fout<<"2-A-3"<<endl;
fout<<"2-A-4"<<endl;
fout<<"2-A-5"<<endl;
fout<<"2-A-6"<<endl;
fout<<"2-A-7"<<endl;
fout<<"*2-A-8"<<endl;
fout<<"*2-A-9"<<endl;
fout<<"*2-A-10"<<endl;
fout<<"*2-A-11"<<endl;
fout<<"*2-A-12"<<endl;
fout<<"*2-A-13"<<endl;
fout<<"*2-A-14"<<endl;
fout<<"*2-A-15"<<endl;
fout<<"*2-A-16"<<endl;
fout<<""<<endl;
fout<<"*2-B-1"<<endl;
fout<<"*2-B-2"<<endl;
fout<<"*2-B-3"<<endl;
fout<<"*2-B-4"<<endl;
fout<<"*2-B-5"<<endl;
fout<<"*2-B-6"<<endl;
fout<<"*2-B-7"<<endl;
fout<<"*2-B-8"<<endl;
fout<<"*2-B-9"<<endl;
fout<<"*2-B-10"<<endl;
fout<<"*2-B-11"<<endl;
fout<<"*2-B-12"<<endl;
fout<<"*2-B-13"<<endl;
fout<<"*2-B-14"<<endl;
fout<<"*2-B-15"<<endl;
fout<<"*2-B-16"<<endl;
fout.close();

cout << "Sample task file " << file << " is created." << endl;
}

