#include <stdio.h>
#include <string.h>
#include <iostream.h>

unsigned char bin(unsigned char ascii) {
  switch(ascii) {
    case 0x0A: return 0x0A;
    case 0x0D: return 0x8D;
    case '%':  return 0xA5;
    case '(':  return 0x28;
    case ')':  return 0xA9;
    case '*':  return 0xAA;
    case '+':  return 0x2B;
    case ',':  return 0xA0;
    case '-':  return 0x2D;
    case '.':  return 0x2E;
    case '/':  return 0xAF;
    case '0':  return 0x30;
    case '1':  return 0xB1;
    case '2':  return 0xB2;
    case '3':  return 0x33;
    case '4':  return 0xB4;
    case '5':  return 0x35;
    case '6':  return 0x36;
    case '7':  return 0xB7;
    case '8':  return 0xB8;
    case '9':  return 0x39;
    case 'A':  return 0x41;
    case 'B':  return 0x42;
    case 'C':  return 0xC3;
    case 'D':  return 0x44;
    case 'E':  return 0xC5;
    case 'F':  return 0xC6;
    case 'G':  return 0x47;
    case 'H':  return 0x48;
    case 'I':  return 0xC9;
    case 'J':  return 0xCA;
    case 'K':  return 0x4B;
    case 'L':  return 0xCC;
    case 'M':  return 0x4D;
    case 'N':  return 0x4E;
    case 'O':  return 0xCF;
    case 'P':  return 0x50;
    case 'Q':  return 0xD1;
    case 'R':  return 0xD2;
    case 'S':  return 0x53;
    case 'T':  return 0xD4;
    case 'U':  return 0x55;
    case 'V':  return 0x56;
    case 'W':  return 0xD7;
    case 'X':  return 0xD8;
    case 'Y':  return 0x59;
    case 'Z':  return 0x5A;
    default :  return 0xFF;
  }
}

unsigned char ascii(unsigned char bin) {
  switch(bin) {
    case 0x0A: return 0x0A;
    case 0x8D: return 0x0D;
    case 0xA5:  return '%';
    case 0x28:  return '(';
    case 0xA9:  return ')';
    case 0xAA:  return '*';
    case 0x2B:  return '+';
    case 0xA0:  return ',';
    case 0x2D:  return '-';
    case 0x2E:  return '.';
    case 0xAF:  return '/';
    case 0x30:  return '0';
    case 0xB1:  return '1';
    case 0xB2:  return '2';
    case 0x33:  return '3';
    case 0xB4:  return '4';
    case 0x35:  return '5';
    case 0x36:  return '6';
    case 0xB7:  return '7';
    case 0xB8:  return '8';
    case 0x39:  return '9';
    case 0x41:  return 'A';
    case 0x42:  return 'B';
    case 0xC3:  return 'C';
    case 0x44:  return 'D';
    case 0xC5:  return 'E';
    case 0xC6:  return 'F';
    case 0x47:  return 'G';
    case 0x48:  return 'H';
    case 0xC9:  return 'I';
    case 0xCA:  return 'J';
    case 0x4B:  return 'K';
    case 0xCC:  return 'L';
    case 0x4D:  return 'M';
    case 0x4E:  return 'N';
    case 0xCF:  return 'O';
    case 0x50:  return 'P';
    case 0xD1:  return 'Q';
    case 0xD2:  return 'R';
    case 0x53:  return 'S';
    case 0xD4:  return 'T';
    case 0x55:  return 'U';
    case 0x56:  return 'V';
    case 0xD7:  return 'W';
    case 0xD8:  return 'X';
    case 0x59:  return 'Y';
    case 0x5A:  return 'Z';
    default  :  return 0xFF;
  }
}

int Quit() {
  char byte = 'X';
  do {
    cout << "<CR>=Continue, Q=Quit" << endl;
    byte = getchar();
    if (byte != '\n') while (getchar() != '\n'); // flushing stdin
    if (byte == 'q' || byte == 'Q') return 1;
  }
  while (byte != '\n');
  return 0;
}

int main(int argc, char *argv[]) {
	char inname[32], outname[32];
  unsigned char byte, outbyte, null='\0', end=0x84;
	FILE *infile, *outfile;
  int txt2bin=1;
  int i, offset=-1;

	if (argc > 3) {
    cerr << "Extra parameters" << endl;
    cout << "Usage: " << argv[0] << " infile outfile" << endl;
    return 1;
  }

	if (argc == 1) {
    cout << "Usage: " << argv[0] << " infile outfile" << endl << endl;
    cout << "Enter input  file name (^C to quit): ";
    cin >> inname;
    cout << "Enter output file name: ";
    cin >> outname;
  }
  else {
    // input file name only
    #pragma warn -pia
    for(i=0; inname[i]  = argv[1][i]; i++);
    #pragma warn +pia
  }

	infile = fopen(inname, "rb");
	if (infile == NULL) {
    cerr << "File " << inname << " did not found" << endl;
    return 1;
  }

  // identify the infile
  byte = fgetc(infile);
  if (byte == '\0') txt2bin=0;
  fclose(infile);

  // output file name
  #pragma warn -pia
	if (argc == 3) for(i=0; outname[i] = argv[2][i]; i++);
  else {
    for(i=0; outname[i] = inname[i]; i++) if(inname[i]=='.') break;
    // close name to use strcat
    outname[++i] = '\0';
    if (txt2bin) strcat(outname,"bin");
    else         strcat(outname,"txt");
  }
  #pragma warn +pia

  if (txt2bin) cout << "Convert text file "   << inname
                    << " to binary file "     << outname << endl;
  else         cout << "Convert binary file " << inname
                    << " to text file "       << outname << endl;
  if (Quit()) return 0;

	outfile = fopen(outname, "rb");
	if (outfile != NULL) {
    fclose(outfile);
    cout << "File " << outname << " already exits." << endl;
    if (Quit()) return 0;
  }

  infile  = fopen(inname,  "rb");
  outfile = fopen(outname, "wb");

  if (txt2bin) {
    // write head nulls
    for (i=0; i<100; i++) fwrite(&null,1,1,outfile);

    while (fread(&byte,1,1,infile)) {
      offset++;
      cout << byte;
      outbyte = bin(byte);
      if (outbyte != 0xFF) fwrite(&outbyte,1,1,outfile);
      else {
        cerr << endl << endl
             << "Convert error: Unknown symbol with ASCII code" << endl
             << "decimal " << (int) byte << endl;
        cout.setf(ios::uppercase);
        cout << "hex     " << hex << (int) byte << endl;
        cout << "Hex offset in the input file is "  << hex << offset << endl;
        fclose(infile);
        fclose(outfile);
        remove(outname);
        cout << endl << "Output file did not written. Press <CR> to quit";
        getchar();
        return 1;
      }
    }
    // write termination symbol
    fwrite(&end,1,1,outfile);
  }
  else {
    while (fread(&byte,1,1,infile)) {
      offset++;
      if (byte == '\0') continue;
      if (byte == end)  break;

      outbyte = ascii(byte);
      cout << outbyte;
      if (outbyte != 0xFF) fwrite(&outbyte,1,1,outfile);
      else {
        cerr << endl << endl
             << "Convert error: Unknown symbol with ASCII code" << endl
             << "decimal " << (int) byte << endl;
        cout.setf(ios::uppercase);
        cout << "hex     " << hex << (int) byte << endl;
        cout << "Hex offset in the input file is "  << hex << offset << endl;
        fclose(infile);
        fclose(outfile);
        remove(outname);
        cout << endl << "Output file did not written. Press <CR> to quit";
        getchar();
        return 1;
      }
    }
    if (byte != end)
    {
      cerr << endl
           << "Convert error: Missing of termination symbol "
           << "with ASCII decimal code " << (int) end << endl;
      fclose(infile);
      fclose(outfile);
      remove(outname);
      cout << "Output file did not written. Press <CR> to quit";
      cin >> byte;
      return 1;
    }
  }

  fclose(infile);
  fclose(outfile);
  return 0;
}
