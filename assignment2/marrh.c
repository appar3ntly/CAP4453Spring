#include <stdio.h>                  /*  Marr-Hildreth.c  (or marrh.c) */
#include <math.h>
#include <stdlib.h>
#define  PICSIZE 256
#define  MAXMASK 100

// USE value of 4.7 to obtain proper threshold results

         int    pic[PICSIZE][PICSIZE];
         double xconv[PICSIZE][PICSIZE];
         double yconv[PICSIZE][PICSIZE];
         int    edgeflag[PICSIZE][PICSIZE];
         double cand[PICSIZE][PICSIZE];
         double xmask[MAXMASK][MAXMASK];
         double ymask[MAXMASK][MAXMASK];
         double mag[PICSIZE][PICSIZE];
         int histogram[PICSIZE];

main(argc,argv)
int argc;
char **argv;
{
        int     i,j,p,q,s,t,mr,centx,centy, more_to_do;
        double  maskval,sum,sig,maxival,slope, high_thresh, low_thresh, cutoff, area_of_tops, percent;
        FILE    *mag_outfile, *peak_outfile, *edge_outfile, *fp1, *fopen();
        char    *foobar;

        argc--; argv++;
        foobar = *argv;
        fp1=fopen(foobar,"rb");

        argc--; argv++;
        foobar = *argv;
        mag_outfile=fopen(foobar,"wb");

        argc--; argv++;
        foobar = *argv;
        peak_outfile=fopen(foobar,"wb");

        argc--; argv++;
        foobar = *argv;
        edge_outfile=fopen(foobar,"wb");


        argc--; argv++;
        foobar = *argv;
        sig = atoi(foobar);

        argc--; argv++;
        foobar = *argv;
        percent = atof(foobar);

        percent *= 0.01;

        mr = (int)(sig * 3);
        centx = (MAXMASK / 2);
        centy = (MAXMASK / 2);

        // Ignore .pgm header
        char dummy[20];
        fgets(dummy, sizeof(dummy), fp1);
        fgets(dummy, sizeof(dummy), fp1);
        fgets(dummy, sizeof(dummy), fp1);

        for (i=0;i<256;i++)
        { for (j=0;j<256;j++)
                {
                  pic[i][j]  =  getc (fp1);
                }
        }

        // Generate x and y masks for convolution from 1st Gaussian derivatives
        for (p=-mr;p<=mr;p++)
        {  for (q=-mr;q<=mr;q++)
           {
              maskval = -p*
                      (exp(-1*(((p*p)+(q*q))/(2*(sig*sig)))));
              xmask[p+centy][q+centx] = maskval;
              maskval = -q*
                      (exp(-1*(((p*p)+(q*q))/(2*(sig*sig)))));
              ymask[p+centy][q+centx] = maskval;
           }
        }

        // Convolve for x gradient
        for (i=mr;i<=255-mr;i++)
        { for (j=mr;j<=255-mr;j++)
          {
             sum = 0;
             for (p=-mr;p<=mr;p++)
             {
                for (q=-mr;q<=mr;q++)
                {
                   sum += pic[i+p][j+q] * xmask[p+centy][q+centx];
                }
             }
             xconv[i][j] = sum;
          }
        }

        // Convolve for y gradient
        for (i=mr;i<=255-mr;i++)
        { for (j=mr;j<=255-mr;j++)
          {
             sum = 0;
             for (p=-mr;p<=mr;p++)
             {
                for (q=-mr;q<=mr;q++)
                {
                   sum += pic[i+p][j+q] * ymask[p+centy][q+centx];
                }
             }
             yconv[i][j] = sum;
          }
        }

        // Combine gradient information to get gradient magnitudes
        maxival = 0;
        for (i=mr;i<256-mr;i++)
        { for (j=mr;j<256-mr;j++)
          {
             mag[i][j]=sqrt((double)((xconv[i][j]*xconv[i][j]) +
                                      (yconv[i][j]*yconv[i][j])));
             if (mag[i][j] > maxival)
                maxival = mag[i][j];

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

        fprintf(peak_outfile, "%s\n", "P5");
        fprintf(peak_outfile, "%d %d\n", 256, 256);
        fprintf(peak_outfile, "%d\n", 255);

        fprintf(edge_outfile, "%s\n", "P5");
        fprintf(edge_outfile, "%d %d\n", 256, 256);
        fprintf(edge_outfile, "%d\n", 255);



        for (i=0;i<256;i++)
          { for (j=0;j<256;j++)
            {
             mag[i][j] = (mag[i][j] / maxival) * 255;            
             fprintf(mag_outfile,"%c",(char)((int)(mag[i][j])));
             
            }
          }
    
        // I don't understand why but i and j (row, col) appears
        // to flip here from what's expected...
        for (i=mr;i<256-mr;i++)
        {
            for(j=mr;j<256-mr;j++)
            {
                if((xconv[i][j] == 0.0))
                {
                    xconv[i][j] = 0.00001;
                }
                slope = yconv[i][j]/xconv[i][j];
                if((slope <= .4142)&&(slope > -0.4142))
                {
                    if((mag[i][j]>mag[i-1][j])&&(mag[i][j]>mag[i+1][j]))
                    {
                        cand[i][j] = 255;
                    }
                }
                else if((slope <= 2.4142)&&(slope > .4142))
                {
                    if((mag[i][j]>mag[i+1][j+1])&&(mag[i][j]>mag[i-1][j-1]))
                    {
                        cand[i][j] = 255;
                    }
                }
                else if ((slope <= -.4142)&&(slope > -2.4142))
                {
                    if((mag[i][j]>mag[i-1][j+1])&&(mag[i][j]>mag[i+1][j-1]))
                    {
                        cand[i][j] = 255;
                    }
                }
                else
                {
                    if((mag[i][j]>mag[i][j-1])&&(mag[i][j]>mag[i][j+1]))
                    {
                        cand[i][j] = 255;
                    }
                }
            }
        }


        // Output Peaks image
        for (i=0;i<256;i++)
          { for (j=0;j<256;j++)
            {
             fprintf(peak_outfile,"%c",(char)((int)(cand[i][j])));
            }
          }

        // Zero out edge flag matrix
        for (i = 0; i < 256; i++)
            for (j = 0; j < 256; j++)
                edgeflag[i][j] = 0;

        // Build histogram of scaled magnitudes
        for (i = 0; i < 256; i++)
        {
            for (j = 0; j < 256; j++)
            {
                histogram[(int)mag[i][j]]++;
            }
        }

        // I don't understand this part at all
        cutoff = percent*256*256;
        area_of_tops = 0;
        for (i = 256; i > 0; i--)
        {
            area_of_tops += histogram[i];
            if (area_of_tops > cutoff)
            {
                break;
            }
        }

        // Threshold the peaks found to find edges
        printf("%d\n", i);
        high_thresh = i;
        low_thresh = 0.35*high_thresh;
 
        for (i = 0; i < 256; i++)
        {
            for (j = 0; j < 256; j++)
            {
                if (cand[i][j] == 255)
                {
                    if (mag[i][j] > high_thresh)
                    {
                        cand[i][j] = 0;
                        edgeflag[i][j] = 255;
                    }
                    else if (mag[i][j] < low_thresh)
                    {
                        cand[i][j] = 0;
                        edgeflag[i][j] = 0;
                    }
                }
            }
        }

        more_to_do = 1;
        while (more_to_do)
        {
            more_to_do = 0;
            for (i = 0; i < 256; i++)
            {
                for (j = 0; j < 256; j++)
                {
                    if (cand[i][j] == 255)
                    {
                        for (p = -1; p <= 1; p++)
                        {
                            for (q = -1; q <= 1; q++)
                            {
                                if (edgeflag[i+p][j+q] == 255)
                                {
                                    cand[i][j] = 0;
                                    edgeflag[i][j] = 255;
                                    more_to_do = 1;
                                }
                            }
                        }
                    }
                }
            }
        }

        
        // Output the edge detected image to a file
        for (i=0;i<256;i++)
          { for (j=0;j<256;j++)
            {
             fprintf(edge_outfile,"%c",(char)((int)(edgeflag[i][j])));
            }
          }


}

