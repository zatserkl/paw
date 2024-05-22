// Andriy_Zatserklyaniy@fnal.gov

#if !defined(__STDLIB_H)
#include <stdlib.h>
#endif

#if !defined(__STDIO_H)
#include <stdio.h>
#endif

#if !defined(__ALLOC_H)
#include <alloc.h>
#endif

/*
#define Borland

#ifdef Borland
#include <conio.h>  // For Borland C++
#endif
*/

void MESS(char *mess) {
  if (mess[0]) printf("%s\n", mess);
}

void WAIT(char *mess) {
  char prompt[] = "<CR>=Continue, Q=Quit: ";
  if (mess[0]) printf("%s\n", mess);
  printf("%s", prompt);
  char ch = getchar();
  //char ch = getche();
  if (ch=='q' || ch=='Q') exit(1);
}

void WAIT() {
  char prompt[] = "<CR>=Continue, Q=Quit: ";
  printf("%s", prompt);
  char ch = getchar();
  //char ch = getche();
  if (ch=='q' || ch=='Q') exit(1);
}

void STOP(char *mess) {
  if (mess[0]) printf("%s\n", mess);
  exit(1);
}

void STOPCR(char *mess) {
  if (mess[0]) printf("%s\n", mess);
  printf("Hit <CR> to quit\n");
  getchar();
  exit(1);
}

void STOPCR() {
  printf("Hit <CR> to quit\n");
  getchar();
  exit(1);
}

void STOP() {
  exit(1);
}

long mem() {return (long) coreleft();}

static long HEAP_MARK;

void mark_heap() {HEAP_MARK = (long) coreleft();}

long dheap() {return HEAP_MARK - (long) coreleft();}

long dheap(char *mess) {
  long after = HEAP_MARK - (long) coreleft();
  printf("%s %ld\n", mess, after);
  return after;
}
