/*
    PARAMETERS.C

    Calculates putative parameters of schematic 4-designs.

    Vedran Krcadinac (krcko@math.hr), 20.11.2021.
    Department of Mathematics, University of Zagreb, Croatia
*/

#include <stdio.h>
#include <stdlib.h>

#define MAXV 1000

/****************/
/* Global stuff */
/****************/

long int bin[MAXV+1][5]; /* Binomial coefficients */

long int lambda(long int v, long int k, long int l, long int i)
{ long int numer,denom;

  numer=l*bin[v-i][4-i];
  denom=bin[k-i][4-i];

  if (numer%denom) return 0;
  else return numer/denom;
}

/* Subfactorial function */
long int sub(long int x, long int i)
{ long int j,res;
  res=x;
  for (j=1; j<i; ++j) res*=x-j;
  return res;
}

/****************/
/* Main program */
/****************/

int main(int argc,char *argv[])
{ long int n,i,j,ok,count,count2,count3;
  long int v,k,l,maxl;
  long int ll[5];
  long int x,y,z;
  long int n1,n2,n3;
  long int numer,denom;
  long int A,B,C;

  /* Command line arguments */

  for(i=1; i<argc; ++i) 
  { j = 0;
    while (argv[i][j] != '\0')
    { /* Help */
      if ((argv[i][j] == 'h') || (argv[i][j] == 'H') || (argv[i][j] == '?'))
      { printf("Usage: parameters [options]\n");
	printf("Options (lowercase for yes, uppercase for no):\n");
        printf("\n");
	exit(0);
      }
      ++j;
    }
  }  

  for (i=0; i<=MAXV; ++i) bin[i][0]=1;
  for (i=0; i<5; ++i) bin[i][i]=1;
  for (i=2; i<=MAXV; ++i) for (j=1; j<5 && j<i; ++j) bin[i][j]=bin[i-1][j-1]+bin[i-1][j];
  /* for (i=0; i<=MAXV; ++i)
  { for (j=0; j<5 && j<=i; ++j) printf("%li ",bin[i][j]);
    printf("\n");
  }
  exit(0); */

  count=0;
  count2=0;
  count3=0;
  for (v=8; v<=MAXV; ++v) for (k=5; 2*k<=v; ++k) 
  { maxl=4L*bin[k][4]/((long int)(v-3));
    for (l=1; l<=maxl; ++l)
    { ok=1;
      for (i=0; i<5 && ok; ++i) ok=(ll[i]=lambda(v,k,l,i))!=0;
      if (v==23 && k==7 && l==1) ok=0;
      if (ok) 
      { ++count;
        A=ll[0]-1L;
        B=k*(ll[1]-1L);
        C=k*(k-1L)*(ll[2]-1L);	
        for (x=0; x<k-2; ++x) for (y=x+1; y<k-1; ++y) for (z=y+1; z<k; ++z)
        { numer=y*z*A + (1L-y-z)*B + C; 
	  denom=(y-x)*(z-x);
	  if ((numer%denom)==0)
	  { n1=numer/denom;
            numer=x*z*A + (1L-x-z)*B + C; 
            denom=(x-y)*(z-y);
	    if ((numer%denom)==0)
            { n2=numer/denom;
              numer=x*y*A + (1L-x-y)*B + C; 
	      denom=(x-z)*(y-z);
	      if ((numer%denom)==0)
              { n3=numer/denom;
	        ok=(n1>=0) && (n2>=0) && (n3>=0);
	        if (ok) ++count3;
	        for (i=3; ok && i<=4; ++i) ok=sub(x,i)*n1 + sub(y,i)*n2 + sub(z,i)*n3 == sub(k, i)*(ll[i] - 1);
	        if (ok)
                { ++count2;
	          printf("%li & %li & %li & %li & %li & %li & %li & %li & %li & %li \\\\ \n",count2,v,k,l,x,y,z,n1,n2,n3); 
	          fflush(stdout);
                } 
	      }	
            }
          }
        }
      } 
    } 
  }

  printf("Number of (v,k,lambda)    : %ld\n",count);
  printf("No. of (v,k,lambda,x,y,z) : %ld\n",count3);
  printf("No. of putative parameters: %ld\n",count2);

}
