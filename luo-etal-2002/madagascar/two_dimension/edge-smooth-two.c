/* One Dimensional Mean Smoothing */
#include <rsf.h>
#include "array_math.h"

int main (int argc, char* argv[])
{
    int n1, r1, n2, r2;
	int flag;
    sf_file in, out;

    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");

    /* Get sampling dimensions */
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1=");
    if (!sf_histint(in,"n2",&n2)) sf_error("No n2=");
    /* Get smoothing radii */
    if (!sf_getint("r1",&r1)) sf_error("Need r1=");
	if (!sf_getint("r2",&r2)) sf_error("Need r2=");
    /* smoothing radius */
	
	if(!sf_getint("flag",&flag)) flag=1;
	/* reproducibility flag.  0 changes the stencil to 4x4 so we may reproduce the paper's results */
	

    /* allocate input array */
    float* input = sf_floatalloc(n1*n2);
	/* allocate output array */
	float* output = sf_floatalloc(n1*n2);
	/* read array from file */
	sf_floatread(input,n1*n2,in);
	/* allocate array for storing mean */
	float* mean = sf_floatalloc(n1*n2);
	/* allocate array for storing variance */
	float* vari = sf_floatalloc(n1*n2);
    /* loop through data, calculate mean */
	for (int i2 = 0; i2 < n2 ; i2++){
	    for (int i1 = 0; i1 < n1 ; i1++){
			 /* calculate mean value */
	         mean[i2*n1+i1] = window_mean(input,i1,i2,r1,r2,n1,n2,flag);		
			 /* and the corresponding variance */
			 vari[i2*n1+i1] = window_variance(input,i1,i2,r1,r2,n1,n2,mean[i2*n1+i1],flag);
	    }
	}
	/* find the mean value corresponding to minimum variance value in each window */
	for (int i2 = 0 ; i2 < n2 ; i2++){
		for (int i1 = 0; i1 < n1 ; i1++){
			output [ i2*n1 + i1 ] = mean[ find_min_index(vari,i1,i2,r1,r2,n1,n2,flag) ];
		}
	}
    /* write to file */
	sf_floatwrite( output, n1*n2, out );
    /* exit program */
    exit (0);
}
