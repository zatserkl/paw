#include <stdlib.h>
#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <iostream.h>

struct Table {
  int ns;
  unsigned char  b[257];
  unsigned short a[256];
  int            n[256];
  Table();
};

Table::Table() {
  ns = 0;
  for(int i=0; i<256; i++) {
    b[i] = '\0';
    a[i] = 0;
    n[i] = 0;
  }
  b[257] = '\0';
}

int main(int argc, char *argv[]) {
	char namebin[32], nametxt[32], namelst[32];
  char extbin[] = ".bin";
  char exttxt[] = ".txt";
  char extlst[] = ".lst";
	FILE *fbin, *ftxt, *flst;
  int i;

	if (argc > 3) {
    #pragma warn -pia
    for(i=0; namebin[i] = argv[1][i]; i++);
    for(i=0; nametxt[i] = argv[2][i]; i++);
    for(i=0; namelst[i] = argv[3][i]; i++);
    #pragma warn +pia
  }
	else {
    cout << "Usage: " << argv[0] << " binfile txtfile lstfile" << endl;
    return 1;
  }

  if (strchr(namebin, '.') == NULL) strcat(namebin,extbin);
  if (strchr(nametxt, '.') == NULL) strcat(nametxt,exttxt);
  if (strchr(namelst, '.') == NULL) strcat(namelst,extlst);

  cout << "binfile = " << namebin << endl;
  cout << "txtfile = " << nametxt << endl;
  cout << "lstfile = " << namelst << endl;

	fbin = fopen(namebin, "rb");
	if (fbin == NULL) {
    cerr << "File " << namebin << " did not found" << endl;
    return 1;
  }
	ftxt = fopen(nametxt, "rb");
	if (ftxt == NULL) {
    cerr << "File " << nametxt << " did not found" << endl;
    return 1;
  }

  Table table;


  unsigned char end = 0x84;

  unsigned char binbyte, txtbyte;

  while (binbyte = fgetc(fbin), binbyte==0x0);
  ungetc(binbyte,fbin);

  printf(  "     #   Bin               ASCII      Symbol\n");
  int n;

  int res;
  while ((res = fgetc(fbin)) != end)
  {
    if (res == EOF) {
      cerr << "Absence of terminating symbol" << endl;
      return 1;
    }
    binbyte = res;

    res = fgetc(ftxt);
    if (res == EOF) {
      cerr << "Information length of file " << nametxt
           << " is less then information length of file " << namebin
           << endl;
      return 1;
    }
    txtbyte = res;

    for(i=0; i<table.ns; i++) if (table.b[i] == binbyte) goto FOUND;

    table.b[table.ns] = binbyte;
    table.a[table.ns] = txtbyte;
    table.n[table.ns] = 1;
    n = table.ns;
    table.ns++;

    printf("%5d %6X       %6X   %6u   ", n, table.b[n],table.a[n],table.a[n]);
    if (isprint(table.a[n])) printf("   %c\n", table.a[n]);
    //else printf("\n");
    else switch (table.a[n]) {
      case 0x0A: printf("   LF\n"); break;
      case 0x0D: printf("   CR\n"); break;
      default: printf("\n");
    }

    continue;
    FOUND:
    if (table.a[i] != txtbyte) {
      printf("Different representation of binary symbol %X\n", binbyte);
      printf("with ACSII symbols %X and %X\n", table.a[i], txtbyte);
      return 0;
    }
    else table.n[i]++;
  }

  fclose(fbin);
  fclose(ftxt);

	flst = fopen(namelst, "wt");
  printf("\n\n--- Summary ---\n\n");
  printf(      "Termination symbol 0x84\n");
  fprintf(flst,"Termination symbol 0x84\n");
  printf(        "       Bin             ASCII      Symbol\n");
  fprintf(flst,  "       Bin             ASCII      Symbol\n");

  unsigned short j;
  for(j=1; j<256; j++)
  for(i=0; i<table.ns; i++) if (j == table.a[i])
  {
    printf("%6X  %6u    %6X  %6u   ", table.b[i],table.b[i],table.a[i],table.a[i]);
    fprintf(flst,"%6X %6u   %6X  %6u   ", table.b[i],table.b[i],table.a[i],table.a[i]);
    if (isprint(table.a[i])) {
      printf("%c\n", table.a[i]);
      fprintf(flst,"   %c\n", table.a[i]);
    }
    //else {
    //  printf("\n");
    //  fprintf(flst,"\n");
    //}
    else switch (table.a[i]) {
      case 0x0A: printf("   LF\n"); fprintf(flst,"   LF\n"); break;
      case 0x0D: printf("   CR\n"); fprintf(flst,"   CR\n"); break;
      default: printf("\n"); fprintf(flst,"\n");
    }
  }
  fclose(flst);

  return 0;
}
