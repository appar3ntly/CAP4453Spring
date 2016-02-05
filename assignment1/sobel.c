#include <stdio.h>                          /* Sobel.c */
#include <math.h>

int pic[256][256];
int outpicx[256][256];
int outpicy[256][256];
int maskx[3][3] = {{-1,0,1},{-2,0,2},{-1,0,1}};
int masky[3][3] = {{1,2,1},{0,0,0},{-1,-2,-1}};
double ival[256][256],maxival;

int main(int argc,char** argv)
{
    int i,j,p,q,mr,sum1,sum2;
    int low_threshold, high_threshold;
    FILE *mag_outfile, *low_outfile, *high_outfile, *fp1, *fopen();
    char *foobar;

    argc--; argv++;
    foobar = *argv;
    fp1=fopen(foobar,"rb");

	argc--; argv++;
	foobar = *argv;
	mag_outfile=fopen(foobar,"wb");

    argc--; argv++;
	foobar = *argv;
	low_outfile=fopen(foobar,"wb");

    argc--; argv++;
	foobar = *argv;
	high_outfile=fopen(foobar,"wb");

    argc--; argv++;
	foobar = *argv;
	low_threshold = atoi(foobar);
    
    argc--; argv++;
    foobar = *argv;
    high_threshold = atoi(foobar);

    // Read in and throw away the pgm header information.
    // This is really crude and needs to be done better for more possible formats without newlines
    char dummy[20];
    fgets(dummy, sizeof(dummy), fp1);
    fgets(dummy, sizeof(dummy), fp1);
    fgets(dummy, sizeof(dummy), fp1);

    for (i=0;i<256;i++)
    { for (j=0;j<256;j++)
            {
              pic[i][j]  =  getc (fp1);
              pic[i][j]  &= 0377;
            }
    }

    // Add the PGM header formatted as:
    // P5
    // rows columns
    // 255
    // For each of the output files
    fprintf(mag_outfile, "%s\n", "P5");
    fprintf(mag_outfile, "%d %d\n", 256, 256);
    fprintf(mag_outfile, "%d\n", 255);

    fprintf(low_outfile, "%s\n", "P5");
    fprintf(low_outfile, "%d %d\n", 256, 256);
    fprintf(low_outfile, "%d\n", 255);

    fprintf(high_outfile, "%s\n", "P5");
    fprintf(high_outfile, "%d %d\n", 256, 256);
    fprintf(high_outfile, "%d\n", 255);

    mr = 1;
    for (i=mr;i<256-mr;i++)
    { for (j=mr;j<256-mr;j++)
      {
         sum1 = 0;
         sum2 = 0;
         for (p=-mr;p<=mr;p++)
         {
            for (q=-mr;q<=mr;q++)
            {
               sum1 += pic[i+p][j+q] * maskx[p+mr][q+mr];
               sum2 += pic[i+p][j+q] * masky[p+mr][q+mr];
            }
         }
         outpicx[i][j] = sum1;
         outpicy[i][j] = sum2;
      }
    }

    maxival = 0;
    for (i=mr;i<256-mr;i++)
    { for (j=mr;j<256-mr;j++)
      {
         ival[i][j]=sqrt((double)((outpicx[i][j]*outpicx[i][j]) +
                                  (outpicy[i][j]*outpicy[i][j])));
         if (ival[i][j] > maxival)
            maxival = ival[i][j];

       }
    }



    for (i=0;i<256;i++)
      { for (j=0;j<256;j++)
        {
         ival[i][j] = (ival[i][j] / maxival) * 255;            
         fprintf(mag_outfile,"%c",(char)((int)(ival[i][j])));
         
        }
      }

    // Output the file resulting from the low threshold value
    for (i=0;i<256;i++)
    {
        for (j=0;j<256;j++)
        {
            if (ival[i][j] > low_threshold)
            {
                fprintf(low_outfile, "%c", (char)255);
            }
            else
            {
                fprintf(low_outfile, "%c", (char)0);
            }
        }
    }

    // Output the file resulting from the high threshold value
    for (i=0;i<256;i++)
    {
        for (j=0;j<256;j++)
        {
            if (ival[i][j] > high_threshold)
            {
                fprintf(high_outfile, "%c", (char)255);
            }
            else
            {
                fprintf(high_outfile, "%c", (char)0);;
            }
        }
    }

}
