// Andriy_Zatserklyaniy@fnal.gov

#if !defined(__STDLIB_H)
#include <stdlib.h>
#endif

#if !defined(__CTYPE_H)
#include <ctype.h>        /* for isdigit */
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

#include "queue.h"
#include "wait.h"
#include "textfile.h"

TextFile::TextFile(char *taskname, char *comments)
{
  SetName = new char[strlen(taskname)+1];
  strcpy(SetName,taskname);
  Comments = new char[strlen(comments)+1];
  strcpy(Comments,comments);

  lines = new Queue<char>;

  FILE *f;
	if ((f = fopen(SetName, "r")) == NULL) {
    cout << "\nFile " << SetName << " not found!\n";
    STOP("STOP");
  }
  // string with one character to use strstr
  char first[] = " ";

  char *line = new char[256];
  while (fgets(line, 256, f))
  {
    if (*line == '\n') continue;

    strtok(line, "\n");

    int i=0;
    while (line[i]==' ') i++;
    first[0] = line[i];
    if (line[i]=='\0' || strstr(Comments,first)) continue;
    Add(line+i);
  }
  fclose(f);
  delete line;
}

char* TextFile::LineStarts(char *match)
{
  for (int i=0; i<lines->Count(); i++)
  {
    char *line = (*lines)[i];
    if (strstr(line,match) == line) return line;
  }
  return NULL;
}

char* TextFile::LineStarts(char *match, int nth)
{
  int n=0;
  for (int i=0; i<lines->Count(); i++)
  {
    char *line = (*lines)[i];
    if (strstr(line,match) == line) n++;
    if (n == nth) return line;
  }
  return NULL;
}

char* TextFile::LineStartsDigit(int nth)
{
  int n=0;
  for (int i=0; i<lines->Count(); i++)
  {
    char *line = (*lines)[i];
    if (isdigit(*line)) n++;
    if (n == nth) return line;
  }
  return NULL;
}

char* TextFile::Rvalue(char *Lvalue)
{
  if (!SingleOccurence(Lvalue)) return NULL;
  for (int i=0; i<lines->Count(); i++)
  {
    char *line = (*lines)[i];
    if (strstr(line,Lvalue) == line)
    {
      // line starts from Lvalue
      char *ptr = strstr(line,"=");
      if (ptr == NULL) {
        cerr << "Missed symbol = in the line" << endl;
        STOPCR(line);
      }
      int j=1;
      while (ptr[j]==' ') j++;
      return ptr+j;
    }
  }
  return NULL;
}

char* TextFile::Rvalue(char *Lvalue, int nth)
{
  int n=0;
  for (int i=0; i<lines->Count(); i++)
  {
    char *line = (*lines)[i];
    if (strstr(line,Lvalue) == line) n++;
    if (n == nth)
    {
      char *ptr = strstr(line,"=");
      if (ptr == NULL) {
        cerr << "Missed symbol = in the line" << endl;
        STOPCR(line);
      }
      int j=1;
      while (ptr[j]==' ') j++;
      return ptr+j;
    }
  }
  return NULL;
}

int TextFile::Occurence(char *match)
{
  int n=0;
  for (int i=0; i<lines->Count(); i++) {
    char *line = (*lines)[i];
    if (strstr(line,match) == line) n++;
  }
  return n;
}

int TextFile::OccurenceDigit()
{
  int n=0;
  for (int i=0; i<lines->Count(); i++) {
    char *line = (*lines)[i];
    if (isdigit(*line)) n++;
  }
  return n;
}

void TextFile::Show() {
  for (int i=0; i<lines->Count(); i++) cout << (*lines)[i] << endl;
}

int TextFile::SingleOccurence(char *match)
{
  int n = Occurence(match);
  if (n == 1) return 1;
  else        return 0;
}

#pragma warn -par
int TextFile::NoSign(char *line) {return 0;}
int TextFile::NotOccures(char *Lvalue) {return 0;}
int TextFile::NotOccuresStart(char *match) {return 0;}
int TextFile::NotOccuresStartDigit() {return 0;}
#pragma warn +par

int TextFile::MoveRvalue(char *Lvalue, int& var)
{
  if (!SingleOccurence(Lvalue)) return NotOccures(Lvalue);
  for (int i=0; i<lines->Count(); i++)
  {
    char *line = (*lines)[i];
    if (strstr(line,Lvalue) == line)
    {
      // line starts from Lvalue
      char *ptr = strstr(line,"=");
      if (ptr == NULL) return NoSign(line);
      int j=1;
      while (ptr[j]==' ') j++;
      sscanf(ptr+j, "%d", &var);
      lines->Kill(i);
      return 1;
    }
  }
  return 0;
}

int TextFile::MoveRvalue(char *Lvalue, long& var)
{
  if (!SingleOccurence(Lvalue)) return 0;
  for (int i=0; i<lines->Count(); i++)
  {
    char *line = (*lines)[i];
    if (strstr(line,Lvalue) == line)
    {
      // line starts from Lvalue
      char *ptr = strstr(line,"=");
      if (ptr == NULL) return NoSign(line);
      int j=1;
      while (ptr[j]==' ') j++;
      sscanf(ptr+j, "%ld", &var);
      lines->Kill(i);
      return 1;
    }
  }
  return 0;
}

int TextFile::MoveRvalue(char *Lvalue, double& var)
{
  if (!SingleOccurence(Lvalue)) return 0;
  for (int i=0; i<lines->Count(); i++)
  {
    char *line = (*lines)[i];
    if (strstr(line,Lvalue) == line)
    {
      // line starts from Lvalue
      char *ptr = strstr(line,"=");
      if (ptr == NULL) return NoSign(line);
      int j=1;
      while (ptr[j]==' ') j++;
      sscanf(ptr+j, "%lf", &var);
      lines->Kill(i);
      return 1;
    }
  }
  return 0;
}

int TextFile::MoveRvalue(char *Lvalue, char *s)
{
  // any number of occurance
  for (int i=0; i<lines->Count(); i++)
  {
    char *line = (*lines)[i];
    if (strstr(line,Lvalue) == line)
    {
      // line starts from Lvalue
      char *ptr = strstr(line,"=");
      if (ptr == NULL) return NoSign(line);
      int j=1;
      while (ptr[j]==' ') j++;
      strcpy(s, ptr+j);
      lines->Kill(i);
      return 1;
    }
  }
  return 0;
}

int TextFile::MoveLineStarts(char *match, char *s)
{
  // any number of occurance
  for (int i=0; i<lines->Count(); i++)
  {
    char *line = (*lines)[i];
    if (strstr(line,match) == line)
    {
      // line starts from match
      strcpy(s, line);
      lines->Kill(i);
      return 1;
    }
  }
  return NotOccures(match);
}

int TextFile::MoveLineStartsDigit(char *s)
{
  // any number of occurance
  for (int i=0; i<lines->Count(); i++)
  {
    char *line = (*lines)[i];
    if (isdigit(*line))
    {
      // line starts from digit
      strcpy(s, line);
      lines->Kill(i);
      return 1;
    }
  }
  return NotOccuresStartDigit();
}

