#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <stdlib.h>
#include <time.h>

float random01()
{
  float scale=RAND_MAX+1.;
  float base=rand()/scale;
  float fine=rand()/scale;
  return base+fine/scale;
}

int main(int argc, char *argv[])

{
 char line[1000];
 char * temp;
 float * positions;
 char* atoms[1000];
 int nat;
 FILE * fileresult;
 FILE * outputfile;
 if (argc !=2)
   {
   printf("Usage: randomize filename");
   return -1;
   }
 srand(time(NULL));
 fileresult = fopen(argv[1],"r");
 if (!fileresult)
   {
   printf ("can't open %s\n",argv[1]);
   return -1;
   }
 strcpy(line,argv[1]);
 strcat(line,"_randomized");
 outputfile = fopen(line,"w");
 if (!outputfile)
   {
   printf ("can't open %s\n",line);
   return -1;
   }
 while (!feof(fileresult))
  {
   fgets(line,1000,fileresult);
   temp = strstr(line,"nat");
   if (temp != 0)
     {
     temp = strstr(temp,"=");
     temp = temp + 1;
     sscanf(temp , "%d", &nat); 
     printf ("found %d atoms \n", nat);  
     }
  }
 fclose(fileresult);
 int i;
 positions = (float *) malloc(3*nat*sizeof(float));
 fileresult = fopen(argv[1],"r");
 float largest = 0;
 float secondlargest =0;
 while (!feof(fileresult))
  {
   fgets(line,1000,fileresult);
   temp = strstr(line,"ATOMIC_POSITIONS");
   fputs(line,outputfile);
   if (temp != 0)
     {
     for (i = 0; i<nat;i++)
      {
      fgets(line,1000,fileresult);
      atoms[i] = (char *) malloc(4*sizeof(char));
      sscanf(line,"%s %f %f %f",(atoms[i]),(positions+i*3),(positions+1+i*3),(positions+2+i*3));
      if (*(positions+i*3) > largest)
       {
       secondlargest = largest;
       largest = *(positions+i*3);
       }
      if (*(positions+1+i*3) > largest)
       {
       secondlargest = largest;
       largest = *(positions+1+i*3);
       }
      if (*(positions+2+i*3) > largest)
       {
       secondlargest = largest;
       largest = *(positions+2+i*3);
       }      
      printf( "Read: %f %f %f ",(*(positions+i*3)),(*(positions+1+i*3)),(*(positions+2+i*3)));
      printf ("From source line :%s", line);

      }
      float difference = largest - secondlargest;
      for (i = 0; i<nat;i++)
      {
        fprintf(outputfile," %s %f %f %f \n",atoms[i],(*(positions+i*3))+random01()*difference/10.0,(*(positions+1+i*3))+random01()*difference/10.0,(*(positions+2+i*3))+random01()*difference/10.0);
      }
      printf("largest: %f   second largest %f difference: %f\n", largest, secondlargest, largest-secondlargest);
     }
  }
 fclose(fileresult);
 fclose(outputfile);
//printf("%d", argc);
return 0;
};
