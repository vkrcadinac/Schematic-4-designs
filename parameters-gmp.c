/*
    PARAMETERS-GMP.C

    Calculates putative parameters of schematic 4-designs using
    the GMP library. Assuming the library is installed, the program 
    is compiled with:

    gcc -o parameters-gmp parameters-gmp.c -lgmp

    Vedran Krcadinac (krcko@math.hr), 20.11.2021.
    Department of Mathematics, University of Zagreb, Croatia
*/

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>

#define MAXV 1000

/****************/
/* Global stuff */
/****************/

long int bin[MAXV+1][5]; /* Binomial coefficients */

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
  long int x,y,z;
  long int mx,my,mz;
  long int br,naz;

  mpz_t numer,denom,tmp,tmp2;
  mpz_t ll[5];
  mpz_t A,B,C;
  mpz_t n1,n2,n3;

  mpz_init(numer);
  mpz_init(denom);
  mpz_init(tmp);
  mpz_init(tmp2);
  for (i=0; i<5; ++i) mpz_init(ll[i]);
  mpz_init(A);
  mpz_init(B);
  mpz_init(C);
  mpz_init(n1);
  mpz_init(n2);
  mpz_init(n3);

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

      for (i=0; i<5 && ok; ++i)
      { mpz_set_si(numer,l);
        mpz_set_si(tmp,bin[v-i][4-i]);
        mpz_mul(numer,numer,tmp);
        mpz_set_si(denom,bin[k-i][4-i]);

        if (mpz_divisible_p(numer,denom)) mpz_divexact(ll[i],numer,denom);
        else ok=0;
      }
           	
      if (v==23 && k==7 && l==1) ok=0; 

      if (ok) 
      { ++count;
        mpz_sub_ui(A,ll[0],1L);
        mpz_sub_ui(B,ll[1],1L);
        mpz_set_si(tmp,k);
	mpz_mul(B,B,tmp);
	mpz_sub_ui(C,ll[2],1L);
	mpz_mul(C,C,tmp);
        mpz_set_si(tmp,k-1);
	mpz_mul(C,C,tmp);
	/* mpz_out_str(stdout,10,A);
        printf(" ");
	mpz_out_str(stdout,10,B);
        printf(" ");
	mpz_out_str(stdout,10,C);
        printf("\n"); */
        for (x=0; x<k-2; ++x) for (y=x+1; y<k-1; ++y) for (z=y+1; z<k; ++z)
        { mpz_mul_si(numer,A,y*z); 
	  mpz_mul_si(tmp,B,1L-y-z);
	  mpz_add(numer,numer,tmp);
	  mpz_add(numer,numer,C);
	  mpz_set_si(denom,(y-x)*(z-x));
          if (mpz_divisible_p(numer,denom)) 
 	  { mpz_divexact(n1,numer,denom);
            mpz_mul_si(numer,A,x*z); 
	    mpz_mul_si(tmp,B,1L-x-z);
	    mpz_add(numer,numer,tmp);
	    mpz_add(numer,numer,C);
	    mpz_set_si(denom,(x-y)*(z-y));
            if (mpz_divisible_p(numer,denom)) 
            { mpz_divexact(n2,numer,denom);
              mpz_mul_si(numer,A,x*y); 
	      mpz_mul_si(tmp,B,1L-x-y);
	      mpz_add(numer,numer,tmp);
	      mpz_add(numer,numer,C);
	      mpz_set_si(denom,(x-z)*(y-z));
              if (mpz_divisible_p(numer,denom)) 
              { mpz_divexact(n3,numer,denom);
   	        ok= (mpz_sgn(n1)>=0) && (mpz_sgn(n2)>=0) && (mpz_sgn(n3)>=0);
                if (ok) ++count3;
	        for (i=3; ok && i<=4; ++i)
                { mpz_mul_si(tmp,n1,sub(x,i));
                  mpz_addmul_ui(tmp,n2,sub(y,i));
                  mpz_addmul_ui(tmp,n3,sub(z,i));
                  mpz_sub_ui(tmp2,ll[i],1L);
		  mpz_mul_si(tmp2,tmp2,sub(k,i));
		  ok=mpz_cmp(tmp,tmp2)==0;
		}
		if (ok)
		{ ++count2;
		  printf("%li & %li & %li & %li & %li & %li & %li & ",count2,v,k,l,x,y,z);
		  mpz_out_str(stdout,10,n1);
		  printf(" & ");
		  mpz_out_str(stdout,10,n2);
		  printf(" & ");
		  mpz_out_str(stdout,10,n3);
		  printf(" \\\\ \n");
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

  mpz_clear(numer);
  mpz_clear(denom);
  mpz_clear(tmp);
  mpz_clear(tmp2);
  for (i=0; i<5; ++i) mpz_clear(ll[i]);
  mpz_clear(A);
  mpz_clear(B);
  mpz_clear(C);
  mpz_clear(n1);
  mpz_clear(n2);
  mpz_clear(n3);

}
