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

#if !defined(__ALLOC_H)
#include <alloc.h>          /* for coreleft() */
#endif

#if !defined(__IOSTREAM_H)
#include <iostream.h>
#endif

#include "queue.h"
#include "elements.h"
#include "megatile.h"
#include "maho.h"
#include "graph.h"
#include "wait.h"

extern void sample(char *file, char *exename);
extern void tile_par(char *file, char *exename);
extern int ReadTask(char *taskname, char *exename);
void write_std(char *namestd, char *exename);
void write_par(char *namepar, char *exename);
void help(char *exename);

void main(int argc, char *argv[])
{
  char exename[256];
  char taskname[256]="";

  mark_heap();
  int ShowHeap = 0;

  long before, after;
  before = long (coreleft());

  // processing of command line argumets

  #pragma warn -pia
  for(int i=0; exename[i] = argv[0][i]; i++);

  if (argc > 1)
  {
    for(i=0; taskname[i] = argv[1][i]; i++);

    char *slash, *key;

    switch (argc) {
      case 2:
        slash = strstr(taskname,"/");
        if (slash) switch(slash[1]) {
          case '?':
          case 'H':
          case 'h':
          help(exename);
          cout << "Press <CR> to quit "; getchar();
          return;

          case 'T':
          case 't': write_std("sample.std", exename);
          cout << "Press <CR> to quit "; getchar();
          return;

          case 'P':
          case 'p': write_par("sample.par", exename);
          cout << "Press <CR> to quit "; getchar();
          return;

          case 'A':
          case 'a':
          cout << endl << "Andriy_Zatserklyaniy@fnal.gov" << endl;
          cout << "Press <CR> to quit "; getchar();
          return;

          case 'M':
          case 'm':
          cout << "Usage of /M key: " << exename << " taskfile /M" << endl;
          cout << "Press <CR> to quit "; getchar();
          return;

          default: cerr << "Unknown key " << taskname << endl;
                   cerr << "To get help run " << exename << " /?" << endl;
          cout << "Press <CR> to quit "; getchar();
          return;
        }
        break;
      case 3:
        key = argv[2];
        slash = strstr(key,"/");
        if (slash) switch(slash[1]) {
          case '?':
          case 'H':
          case 'h':
          help(exename);
          cout << "Press <CR> to quit "; getchar();
          return;

          case 'T':
          case 't': write_std(taskname, exename);
          cout << "Press <CR> to quit "; getchar();
          return;

          case 'P':
          case 'p': write_par(taskname, exename);
          cout << "Press <CR> to quit "; getchar();
          return;

          case 'M':
          case 'm': ShowHeap=1; goto START;

          default: cerr << "Unknown key " << key << endl;
                   cerr << "To get help run " << exename << " /?" << endl;
          cout << "Press <CR> to quit "; getchar();
          return;
        }
      default: cerr << "Extra parameters" << endl;
               cerr << "To get help run " << exename << " /?" << endl;
               cout << "Press <CR> to quit "; getchar();
    }
  }
  else {
    cout << "Enter task file name (extension .STD can be omitted): ";
    cin.get(taskname,'\n');
    while (taskname[0]==' ') for (int i=0; taskname[i]=taskname[i+1]; i++);
  }
  #pragma warn +pia

  START:
  ;

  if (strstr(taskname,".") == 0) strcat(taskname,".std");

  if (!ReadTask(taskname, exename)) {
    cerr << endl << "Error reading task file " << taskname << endl;
    cerr << "To get help run " << exename << " /?" << endl;
    STOPCR();
  }

  GManager manager;
  //MAHO maho;

  manager.Menu();

  if (ShowHeap) {
    dheap("Used heap memory:");
    cout << "Press <CR> to quit "; getchar();
  }

  goto TheEnd;

  #pragma warn -rch
  manager.~GManager();

  //maho.MakePath();
  //maho.MakeCode();
  //maho.~MAHO();

  Megatile::WhoAlive();

  after = long (coreleft());

  cout << "after cre/del: "
       << before << " " << after
       << ", difference = " << before-after << endl;
  cout << "Program terminated. Call of destructor" << endl;

  #pragma warn +rch
  TheEnd:;
}

void write_std(char *namestd, char *exename)
{
  FILE *f = fopen(namestd, "r");
  if (f != NULL) {
    fclose(f);
    cerr << "\nSample task file " << namestd << " already exist" << endl;
    STOPCR("Remove it to run the program");
  }
  sample(namestd, exename);
}

void write_par(char *namepar, char *exename)
{
  FILE *f = fopen(namepar, "r");
  if (f != NULL) {
    fclose(f);
    cerr << "\nSample tiles parameters file " << namepar << " already exist"
         << endl;
    STOPCR("Remove it to run the program");
  }
  tile_par(namepar,exename);
}

void help(char *exename)
{
  cout << endl << "Usage: " << exename << " [taskfile] [/parameter]" << endl;
  cout << "Possible parameter:" << endl;
  cout << "/?" << endl << "   help" << endl;
  cout << "/T" << endl << "   create sample task file" << endl;
  cout << "/P" << endl << "   create sample tiles parameters file" << endl;
  cout << "/M" << endl << "   Usage heap memory" << endl;
  cout << "/A" << endl << "   Author" << endl;
  cout << endl;

  cout << "/parameter should be separated by blank" << endl << endl;

  cout << "Example: " << exename << " maho.std" << endl;
  cout << "Extension .std can be omitted" << endl;
  cout << endl;
  cout << "Example: " << exename << " maho.std /t" << endl;
  cout << "or simple to create sample.std: " << exename << " /t" << endl;
  cout << "Example: " << exename << " sometile.par /p" << endl;
  cout << "or simple to create sample.par: " << exename << " /p" << endl;
}
