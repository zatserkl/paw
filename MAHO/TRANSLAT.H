#ifndef __TRANSLAT_H
#define __TRANSLAT_H

#ifndef __MYCONST_H
#include <heads/myconsts.h>
#endif
#ifndef __CTYPE_H
#include <ctype.h>
#endif
#ifndef __STDIO_H
#include <stdio.h>
#endif
#ifndef __STDLIB_H
#include <stdlib.h>
#endif
#if !defined(__DIR_H)
#include <dir.h>
#endif
#ifndef __STRING_H
#include <string.h>
#endif
BYTE UpToBinChar(BYTE byte)
{
BYTE y=0;
char x=toupper(char(byte));
switch(x)
{
    case 0x0A:y=0x0A; break;
    case 0x0D:y=0x0D; break;
    case 0xFF:y=0xFF; break;
    case 0x00:y=0x00; break;
    case 0x84:y=0x84; break;
    case ' ' :y=0xA0; break;
    case '!' :y= 0x21; break;
    case '%' :y= 0xA5; break;
    case '(' :y= 0x28; break;
    case ')' :y= 0xA9; break;
    case '+' :y= 0x2B; break;
    case '*' :y= 0xAA; break;
    case ',' :y= 0xAC; break;
    case '-' :y= 0x2D; break;
    case '.' :y= 0x2E; break;
    case '/' :y= 0xAF; break;
    case '0' :y= 0x30; break;
    case '1' :y= 0xB1; break;
    case '2' :y= 0xB2; break;
    case '3' :y= 0x33; break;
    case '4' :y= 0xB4; break;
    case '5' :y= 0x35; break;
    case '6' :y= 0x36; break;
    case '7' :y= 0xB7; break;
    case '8' :y= 0xB8; break;
    case '9' :y= 0x39; break;
    case ':' :y= 0x3A; break;
    case ';' :y= 0xBB; break;
    case '<' :y= 0x3C; break;
    case '=' :y= 0xBD; break;
    case '>' :y= 0xBE; break;
    case '?' :y= 0x3F; break;
    case '@' :y= 0xC0; break;
    case 'A' :y= 0x41; break;
    case 'B' :y= 0x42; break;
    case 'C' :y= 0xC3; break;
    case 'D' :y= 0x44; break;
    case 'E' :y= 0xC5; break;
    case 'F' :y= 0xC6; break;
    case 'G' :y= 0x47; break;
    case 'H' :y= 0x48; break;
    case 'I' :y= 0xC9; break;
    case 'J' :y= 0xCA; break;
    case 'K' :y= 0x4B; break;
    case 'L' :y= 0xCC; break;
    case 'M' :y= 0x4D; break;
    case 'N' :y= 0x4E; break;
    case 'O' :y= 0xCF; break;
    case 'P' :y= 0x50; break;
    case 'Q' :y= 0xD1; break;
    case 'R' :y= 0xD2; break;
    case 'S' :y= 0x53; break;
    case 'T' :y= 0xD4; break;
    case 'U' :y= 0x55; break;
    case 'V' :y= 0x56; break;
    case 'W' :y= 0xD7; break;
    case 'X' :y= 0xD8; break;
    case 'Y' :y= 0x59; break;
    case 'Z' :y= 0x5A; break;
    case '[' :y= 0xDB; break;
    case ']' :y= 0xDC; break;
    default :y=0xFF;
};
return BYTE(y);
};
FILE* UpToBin(char fname[MAXPATH])
//возвращает указатель на файл ISO
{
FILE* fin=fopen(fname,"rb");//принимает имя файла в UP
FILE* fout=fopen( "sended.iso","wb");
BYTE byte;
BYTE biso;
int i;
while(!feof(fin)){
    byte=getc(fin);
    biso=UpToBinChar(byte);
    if ((byte!=0x0D)&(biso!=0xFF))putc(biso,fout);
		 };
putc(EOFC,fout);
fclose(fin);
fclose(fout);
fout=fopen("sended.iso","rb");
return fout;
};
BYTE BinToUpChar(BYTE x)
{
BYTE y=0;
switch(x)
{   case 0x0A:y=0x0A; break;
    case 0x0D:y=0x0D; break;
    case 0xA0:y=' ' ; break;
    case 0x21:y='!' ; break;
    case 0xA5:y='%' ; break;
    case 0x28:y='(' ; break;
    case 0xA9:y=')' ; break;
    case 0x2B:y='+' ; break;
    case 0xAA:y='*' ; break;
    case 0xAC:y=',' ; break;
    case 0x2D:y='-' ; break;
    case 0x2E:y='.' ; break;
    case 0xAF:y='/' ; break;
    case 0x30:y='0' ; break;
    case 0xB1:y='1' ; break;
    case 0xB2:y='2' ; break;
    case 0x33:y='3' ; break;
    case 0xB4:y='4' ; break;
    case 0x35:y='5' ; break;
    case 0x36:y='6' ; break;
    case 0xB7:y='7' ; break;
    case 0xB8:y='8' ; break;
    case 0x39:y='9' ; break;
    case 0x3A:y=':' ; break;
    case 0xBB:y=';' ; break;
    case 0x3C:y='<' ; break;
    case 0xBD:y='=' ; break;
    case 0xBE:y='>' ; break;
    case 0x3F:y='?' ; break;
    case 0xC0:y='@' ; break;
    case 0x41:y='A' ; break;
    case 0x42:y='B' ; break;
    case 0xC3:y='C' ; break;
    case 0x44:y='D' ; break;
    case 0xC5:y='E' ; break;
    case 0xC6:y='F' ; break;
    case 0x47:y='G' ; break;
    case 0x48:y='H' ; break;
    case 0xC9:y='I' ; break;
    case 0xCA:y='J' ; break;
    case 0x4B:y='K' ; break;
    case 0xCC:y='L' ; break;
    case 0x4D:y='M' ; break;
    case 0x4E:y='N' ; break;
    case 0xCF:y='O' ; break;
    case 0x50:y='P' ; break;
    case 0xD1:y='Q' ; break;
    case 0xD2:y='R' ; break;
    case 0x53:y='S' ; break;
    case 0xD4:y='T' ; break;
    case 0x55:y='U' ; break;
    case 0x56:y='V' ; break;
    case 0xD7:y='W' ; break;
    case 0xD8:y='X' ; break;
    case 0x59:y='Y' ; break;
    case 0x5A:y='Z' ; break;
    case 0xDB:y='[' ; break;
    case 0xDC:y=']' ; break;
    default:y=0xFF;break;
    };
return y;
};
void BinToUp(char fin[MAXPATH],char fout[MAXPATH])
{
FILE* pfin=fopen(fin,"rb");//принимает имя файла в UP
FILE* pfout=fopen( fout,"wb");
BYTE biso;
BYTE btxt;
int i;
while(!feof(pfin))
     {
    biso=getc(pfin);
    btxt=BinToUpChar(biso);
    if(btxt!=0xFF){
    if(biso==0x0A)putc(0x0D,pfout);
    putc(btxt,pfout);};
     };
fclose(pfin);
fclose(pfout);
};

void BinToUpMaho(char fin[MAXPATH],char fout[MAXPATH])
{
FILE* pfin=fopen(fin,"rb");//принимает имя файла в UP
FILE* pfout=fopen( fout,"wb");
BYTE biso;
BYTE btxt;
int i;
while(getc(pfin)!=0xA5);
getc(pfin);getc(pfin);putc('%',pfout);
while(!feof(pfin))
     {
    biso=getc(pfin);
    btxt=BinToUpChar(biso);
   if(btxt!=0xFF){
	if(btxt==0x8D)putc(0x0D,pfout);
	if(btxt==0xA0)putc(0xAC,pfout);
	if(btxt==0x0A)putc(0x0D,pfout);
	if((btxt!=0x8D)&(btxt!=0x0D))putc(btxt,pfout);
		  };
     };
fclose(pfin);
fclose(pfout);
};

FILE* UpToBinMaho(char fname[MAXPATH])
//возвращает указатель на файл ISO
{
FILE* fin=fopen(fname,"rb");//принимает имя файла в UP
FILE* fout=fopen( "sended.iso","wb");
BYTE byte;
BYTE biso;
int i;
for(i=0;i<100;i++){putc(0,fout);};
i=0;
while(!feof(fin)){
    byte=getc(fin);
    biso=UpToBinChar(byte);
    if (biso!=0xFF){
       if(biso==0x0D)putc(0x8D,fout);
       if(biso==0xAC)putc(0xA0,fout);
       if((biso!=0x0D)&(biso!=0xAC))putc(biso,fout);
	  };
    if (biso==0xA5){putc(0x50,fout);putc(0x4D,fout);};
		 };
putc(0x84,fout);
for(i=0;i<100;i++){putc(0,fout);};
putc(0x2B,fout);
fclose(fin);
fclose(fout);
fout=fopen("sended.iso","rb");
return fout;
};
#endif  //TRANSLAT.H