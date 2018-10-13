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
float calc_variance (float* array, int i1, int r1, int n1, float mean){
	/* takes an array with an index and a smoothing radius, outputs mean value */
	int counter = 0;
	float tracker = 0;
	/* loop through smoothing window */
	for (int j = 0 ; j < 2*r1+1; j++){
		int i = j-r1 + i1;
		/* make sure in bounds */
		if (i < 0){continue;}
		if (i >= n1){continue;}
		/* increment counter */
		counter += 1;
		/* add array value to tracker */
		tracker += (array[i]-mean)*(array[i]-mean);
	}
	/* return variance */
	return tracker / (float)counter;
}
/* returns the minimum of two integers, x and y */
int min(int x, int y){
	int min;
	if (x < y){ min = x; }
	else      { min = y; }
	return min;	
}
/* returns the maximum of two integers, x and y */
int max(int x, int y){
	int max;
	if (x > y){ max = x; }
	else      { max = y; }
	return max;	
}
/* finds the minimal variance index in an input array */
int find_mindex(float* array, int i1, int r1, int n1){
    /* determine the first and last indexes in the window */
	int first = max(i1 - r1, 0 );
	int last  = min(i1 + r1, n1);
	/* initialize minimal index */
	int least = first;
	float minval = array[first];
    /* loop through window to find lowest value */
	for (int ind = first+1 ; ind < last ; ind++){
		if (array[ind] < minval){
			minval = array[ind];
			least = ind;
		}
	}
	return least;
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
	/* allocate array for storing mean values */
	float* mean   = sf_floatalloc(n1);
	/* allocate array for storing variance values */
	float* var    = sf_floatalloc(n1);
	/* read array from file */
	sf_floatread(input,n1,in);
    /* loop through data, calculate mean */
	for (int i1 = 0; i1 < n1 ; i1++){
		/* calculate mean */
	    mean[i1] = calc_mean(input,i1,r1,n1);
		/* calculate variance */
		var[i1]  = calc_variance(input,i1,r1,n1,mean[i1]);
	}
	/* loop through points and select the index with the minimal variance */
	for (int i1 = 0; i1 < n1 ; i1++){
		output[i1] = mean[find_mindex(var,i1,r1,n1)];
		}
    /* write to file */
	sf_floatwrite(output,n1,out);
    /* exit program */
    exit (0);
}
