/* One Dimensional Edge Preserving Smoothing */
#include <rsf.h>

float calc_mean (float* array, int i1, int r1, int n1){
	/* takes an array with an index and a smoothing radius, outputs mean value */
	int counter = 0;
	float tracker = 0;
	/* loop through smoothing window */
	for (int j = 0 ; j < 2*r1+1 ; j++){
		int i = j - r1 + i1;
		/* make sure in bounds */
		if (i < 0){continue;}
		if (i >= n1){continue;}
		/* increment counter */
		counter += 1;
		/* add array value to tracker */
		tracker += array[i];
	}
	/* return mean */
	return tracker / (float)counter;
}

int main (int argc, char* argv[])
{
    int n1, r1;
    sf_file in, out;
    sf_init (argc,argv);
    in = sf_input("in");
    out = sf_output("out");
    /* Get sampling dimension */
    if (!sf_histint(in,"n1",&n1)) sf_error("No n1=");
    /* Get smoothing radius */
    if (!sf_getint("r1",&r1))   sf_error("Need r1=");
    /* smoothing radius */
    /* allocate input array */
    float* input  = sf_floatalloc(n1);
	/* allocate output array */
	float* output = sf_floatalloc(n1);
	/* read array from file */
	sf_floatread(input,n1,in);
    /* loop through data, calculate mean */
	for (int i1 = 0; i1 < n1 ; i1++){
		output[i1] = calc_mean(input,i1,r1,n1);
		}
    /* write to file */
	sf_floatwrite(output,n1,out);
    /* exit program */
    exit (0);
}
