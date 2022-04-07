
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <float.h>

/*
EXAMPLE.ine:

cube.ine
H-representation
begin
4 3 rational
-1.1 1 0
-3.1459 0 1
23 -1 0.2323
42.91238 0 -1.0
end
*/


#define USAGE printf("Usage: %s fin fout\n", argv[0]);

FILE *fin;/* input file pointer       */

FILE *fout;/* output file pointer       */

int main(int argc,char *argv[])
{
  if (argc < 3){
    USAGE;
    return -1;
  }

  long int  m,n;

  int i,j;

  long atol();

  char  buf[BUFSIZ];
  fin = fopen (argv[1], "r+");
  fout = fopen (argv[2], "w+");

  while ( fgets(buf,BUFSIZ,fin) !=NULL ){

    fputs(buf,fout);

    if (strncmp(buf,"begin",5)==0) break;

  }

  if (fscanf(fin,"%ld %ld %s",&m,&n,buf)==EOF){

    fprintf(stderr,"No begin line");

    exit(1);

  }

  fprintf(fout,"%ld %ld rational\n",m,n);

  for (i=0;i<m;i++)   {
    for(j=0;j<n;j++)	{

      char *p;
      char *frac;
      int k;

      fscanf(fin,"%s",buf);

      if ((p=strchr(buf,'.'))){

        *p=0;

        frac=&p[1];

        fprintf(fout,"%s%s/1",buf,frac);

        for (k=0; k<strlen(frac); k++)
          fprintf(fout,"0");

      }else {

        fprintf(fout,"%s",buf);
      }

      fprintf(fout," ");

    }

    fputs("\n",fout);

  }

  fgets(buf,BUFSIZ,fin);

  while (fgets(buf,BUFSIZ,fin) !=NULL ) {

    fputs(buf,fout);

  }

  return 0;
}
