// Andriy_Zatserklyaniy@fnal.gov

#include <stdio.h>
#include <stdlib.h>
#include <conio.h>
#include <string.h>
#include <math.h>

void main(int argc, char *argv[]) {
	char fname[256];
  char *ext = ".dat";
	double chan, ampl;
	double ped, spectrum, sum;
  double under, over;
	int nped = 8;
	FILE *f;
  int i;

  #pragma warn -pia
	if (argc-1) for(i=0; fname[i] = argv[1][i]; i++);
  #pragma warn +pia
	else {
		printf("\nEnter file name: ");
		scanf("%s", fname);
    while(getchar() != '\n'); //flushing stdin after scanf
    if (strchr(fname, '.') == NULL) strcat(fname,ext);
	}
  if (strchr(fname, '.') == NULL) strcat(fname,ext);
	printf("\nProcessing file %s\n", fname);

  #pragma warn -pia
	if (f = fopen(fname, "r")) {
    #pragma warn +pia
    printf("Enter max. channel of the pedestal (default - 8):  ");
    char string[80];
    //  For flushing stdin after scanf we can use:
    //    flushall();
    //      or
    //    fflush(stdin);
    gets(string);
    if(string[0]) nped = atoi(string);

    int maxchan = 1024;
    printf("Enter max. channel of the spectrum (default - 1024):  ");
    gets(string);
    if(string[0]) maxchan = atoi(string);
    printf("\n");

    ped = spectrum = sum = 0L;
    under = over = 0L;
    long ampl_sum=0L, ampl_chan=0L;
    char firstch = '\'';
    chan = ampl = 0;
    char line[256];
    while ((firstch=fgetc(f)) != EOF) {
      if (firstch == '\'') {
        //while (int (firstch=fgetc(f)) != 0x0A);
        printf("%s", fgets(line, 256, f));
        continue;
      }
      ungetc(firstch, f);
      fscanf(f, "%lf %lf", &chan, &ampl);
      //fflush(f);
      if (chan==1) {
        under += ampl;
        continue;
      }
      //
      sum += ampl;
      if (chan > nped) spectrum += ampl;
      else ped += ampl;
      if(chan <= maxchan) {
        ampl_sum += ampl;
        ampl_chan += ampl*chan;
      }
      else {
        over += ampl;
      }
      // in case empty next line scanf will not change previous readen values
      // therefore:
      chan = ampl = 0;
    }

    double cmass = 0.;
    if (ampl_sum) cmass = ampl_chan/ampl_sum;
    printf("\nCenter of mass in %1.1f channel\n", cmass);

    printf("\n");
    printf("sum = %1.f\n", sum);
    printf("spectrum = %1.f\n", spectrum);
    printf("pedestal = %1.f\n", ped);
    printf("\n");
    printf("under = %1.f\n", under);
    printf("over = %1.f\n", over);

    int accident = 0;
    printf("\nEnter number of accident counts (default - 0):  ");
    gets(string);
    if(string[0]) accident = atoi(string);
    printf("\n");

    // correct pedestal by accident counts
    ped = ped - accident;

    double error;
    if (ped > 0.)
      error = sqrt((ped+accident)+accident) / ped;
    else {
      printf("\nZero pedestal counts. Method doesn't work!\n");
      printf(  "==========================================\n");
      printf("\nPress any key\n");
      getch();
      exit(1);
    }

    double ped_ratio = ped/sum;
    double spectrum_ratio = spectrum/sum;
    printf("spectrum/sum ratio = %1.1f%%\n", 100.*spectrum_ratio);
    printf("pedestal/sum ratio = %1.1f%%\n", 100.*ped_ratio);
    printf("\n");

    double pe = -log(ped_ratio);
    printf
      ("Pedestal/sum ratio correspond to %1.4f +/- %1.4f photoelectrons\n",
      pe, error);
    printf("\nPress any key\n");
    char ch = getch();
    if (ch == 0) ch = getch();
  }
  else printf("\nFile %s not found!\n", fname);
}
