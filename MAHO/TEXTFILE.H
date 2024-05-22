// Andriy_Zatserklyaniy@fnal.gov

#if !defined(__STDLIB_H)
#include <stdlib.h>
#endif

#if !defined(__STRING_H)
#include <string.h>
#endif

#if !defined(__IOSTREAM_H)
#include <iostream.h>
#endif

template <class T> class Queue;

struct TextFile {
  char *SetName;
  char *Comments;
  Queue<char> *lines;

  TextFile(char *taskname, char *comments);
 ~TextFile() {
    delete lines; lines=NULL;
    delete SetName; SetName = NULL;
    delete Comments; Comments = NULL;
  }
  void Add(char *theline) {
    char *line = new char[strlen(theline)+1];
    strcpy(line,theline);
    *lines << line;
  }
  int Occurence(char *match);
  int OccurenceDigit();
  char* LineStarts(char *match);
  char* LineStarts(char *match, int nth);
  char* LineStartsDigit(int nth);
  char* Rvalue(char *Lvalue);
  char* Rvalue(char *Lvalue, int nth);
  int   MoveRvalue(char *Lvalue, int& var);
  int   MoveRvalue(char *Lvalue, long& var);
  int   MoveRvalue(char *Lvalue, double& var);
  int   MoveRvalue(char *Lvalue, char *s);
  int   MoveLineStarts(char *match, char *s);
  int   MoveLineStartsDigit(char *s);
  int   Nlines() {return lines->Count();}
  virtual void Show();
  virtual int  SingleOccurence(char *match);
  virtual int  NoSign(char *line);
  virtual int  NotOccures(char *Lvalue);
  virtual int  NotOccuresStart(char *match);
  virtual int  NotOccuresStartDigit();
};

struct TaskText : public TextFile
{
  // struct to supply errors handling of SingleOccurence

  char *ExeName;
  TaskText(char *taskname, char *comments, char *exename):
       TextFile(taskname,comments)
  {
    ExeName = new char[strlen(exename)+1];
    strcpy(ExeName,exename);
  }
 ~TaskText() {delete ExeName; ExeName = NULL;}
  int SingleOccurence(char *match) {
    int n = Occurence(match);
    if (n==0) {
      cerr << "Field "<< match <<" did not defined\n";
      STOPCR();
      return 0;
    }
    if (n > 1) {
      cerr << "Error: Field " << match << endl;
      cerr << "is defined " << n << " times" << endl;
      cerr << "Please, correct the task file" << endl;
      cerr << endl;
      cerr << "To get help run " << ExeName << " /?" << endl;
      STOPCR();
      return 0;
    }
    return 1;
  }
  int NoSign(char *line) {
    cerr << "ReadTask: Missed symbol = in the line" << endl;
    cerr << line << endl;
    STOPCR();
   return 0;
  }
  int NotOccures(char *Lvalue) {
    cerr << "ReadTask: Field " << Lvalue << " did not found" << endl;
    cerr << "To get help run " << ExeName << " /?" << endl;
    STOPCR();
    return 0;
  }
  int NotOccuresStart(char *match) {
    cerr << "ReadTask: Line starts from " << match
         << " did not found" << endl;
    cerr << "To get help run " << ExeName << " /?" << endl;
    STOPCR();
    return 0;
  }
  int NotOccuresStartDigit() {
    cerr << "ReadTask: Line starts from digit did not found" << endl;
    cerr << "To get help run " << ExeName << " /?" << endl;
    STOPCR();
    return 0;
  }
};
