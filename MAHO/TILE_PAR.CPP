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

void tile_par(char *file, char *exename)
{
FILE *f;
f = fopen(file, "r");
if (f != NULL) {
  fclose(f);
  cerr << "\nTiles parameters file " << file << " already exist" << endl;
  STOPCR("Choose another file name");
}

ofstream fout(file);

fout<<"*"<<endl;
fout<<"* Sample of tiles parameters file for " << exename <<endl;
fout<<"* Line that starts from symbols * or # is treated as a commentary"
    <<endl;
fout<<"*"<<endl;
fout<<"#tile   A       B       C       S       H"<<endl;
fout<<"2-A-1   104.9   73.9    175.7   52.0    -"<<endl;
fout<<"2-A-3   146.2   122.0   137.2   72.6    60.6"<<endl;
fout<<"2-A-5   194.7   170.5   137.4   96.6    84.6"<<endl;
fout<<"2-A-7   243.1   219.0   137.2   120.6   108.6"<<endl;
fout<<"2-B-1   105.9   74.6    177.5   52.6    -"<<endl;
fout<<"2-B-3   147.8   123.3   138.7   73.3    61.2"<<endl;
fout<<"2-B-5   196.7   172.3   138.8   97.6    85.5"<<endl;
fout<<"2-B-7   245.7   221.2   138.6   121.9   109.8"<<endl;
fout<<""<<endl;
fout<<"#tile   A       B       C       R       V"<<endl;
fout<<"2-A-9   145.6   132.7   147.9   66.2    72.7"<<endl;
fout<<"2-A-11  160.0   145.6   164.5   72.7    79.9"<<endl;
fout<<"2-A-13  176.1   160.0   183.7   79.9    87.9"<<endl;
fout<<"2-A-15  196.1   176.1   228.5   87.9    -"<<endl;
fout<<"2-B-9   147.1   134.1   149.5   66.9    73.4"<<endl;
fout<<"2-B-11  161.7   147.1   166.2   73.4    80.7"<<endl;
fout<<"2-B-13  177.9   161.7   185.6   80.7    88.8"<<endl;
fout<<"2-B-15  201.1   177.9   265.3   88.8    -"<<endl;
fout<<""<<endl;
fout<<"#tile   A       B       C       H       S       T       W"<<endl;
fout<<"2-A-2   122.0   104.9   97.5    52.0    60.6    -       -"<<endl;
fout<<"2-A-4   170.5   146.2   137.5   72.6    84.6    -       -"<<endl;
fout<<"2-A-6   219.0   194.7   137.5   96.6    108.6   -       -"<<endl;
fout<<"2-A-8   267.4   243.1   137.8   120.6   -       66.2    199.7"<<endl;
fout<<"2-B-2   123.3   105.9   98.5    52.6    61.2    -       -"<<endl;
fout<<"2-B-4   172.3   147.8   138.9   73.3    85.5    -       -"<<endl;
fout<<"2-B-6   221.2   196.7   138.9   97.6    109.8   -       -"<<endl;
fout<<"2-B-8   270.2   245.7   139.2   121.9   -       66.9    201.8"<<endl;
fout<<""<<endl;
fout<<"#tile   E       D       C       F       Q       U"<<endl;
fout<<"1-A-10  131.9   144.8   144.8   157.5   65.6    72.0"<<endl;
fout<<"1-A-12  144.8   159.1   161.1   173.2   72.0    79.1"<<endl;
fout<<"1-A-14  159.1   175.1   179.9   190.8   79.1    87.0"<<endl;
fout<<"1-A-16  175.1   188.9   154.9   202.4   87.0    -"<<endl;
fout<<"1-B-10  133.3   146.3   146.4   159.2   66.3    72.8"<<endl;
fout<<"1-B-12  146.3   160.8   162.8   175.1   72.8    79.9"<<endl;
fout<<"1-B-14  160.8   177.0   181.8   192.9   79.9    88.0"<<endl;
fout<<"1-B-16  177.0   194.0   191.7   210.8   88.0    -"<<endl;
fout.close();

cout << "Tiles parameters file " << file << " is created." << endl;
}

