#include <rsf.h>
/* performs array operations for mean and edge preserving smoothing on 2d arrays */
float window_mean (float* array, int i1, int i2, int r1, int r2, int n1, int n2, int flag){
	/* takes an array with an index and a smoothing radius, outputs mean value */
	int counter = 0;
	float tracker = 0;
	/* loop through smoothing window */
	for (int j2 = i2-r2 ; j2 < i2 + r2 + flag ; j2++){
		if (j2 < 0 ) {continue;}
		if (j2 >= n2) {continue;}
		for (int j1 = i1-r1 ; j1 < i1 + r1 +flag ; j1++){
			/* make sure in bounds */
			if (j1 < 0){continue;}
			if (j1 >= n1){continue;}
			/* increment counter */
			counter += 1;
			/* add array value to tracker */
			tracker += array[ j2*n1 + j1 ];
		}
	}
	/* return mean */
	return tracker / (float)counter;
}
/* performs array operations for mean and edge preserving smoothing on 2d arrays */
float window_variance (float* array, int i1, int i2, int r1, int r2, int n1, int n2, float mean, int flag){
	/* takes an array with an index and a smoothing radius, outputs mean value */
	int counter = 0;
	float tracker = 0;
	/* loop through smoothing window */
	for (int j2 = i2 - r2 ; j2 < i2 + r2 + flag ; j2++){
		if (j2 < 0 ) {continue;}
		if (j2 >= n2) {continue;}
		for (int j1 = i1-r1 ; j1 < i1 + r1 + flag; j1++){
			/* make sure in bounds */
			if (j1 < 0){continue;}
			if (j1 >= n1){continue;}
			/* increment counter */
			counter += 1;
			/* add array value to tracker */
			tracker += sqrt(( array[ j2*n1 + j1 ] - mean ) * ( array[ j2*n1 + j1 ] - mean ));
		}
	}
	/* return mean */
	return (tracker / (float) counter);
}
int find_min_index(float* array, int i1, int i2, int r1, int r2, int n1, int n2, int flag){
	/* initialize minval trackers */
	float mini = -1;
	int mindex = -1;
	for (int j2 = i2-r2 ; j2 < i2 + r2 + flag; j2++){
		if (j2 < 0 ) {continue;}
		if (j2 >= n2) {continue;}
		for (int j1 = i1-r1 ; j1 < i1 + r1 + flag; j1++){
			/* make sure in bounds */
			if (j1 < 0){continue;}
			if (j1 >= n1){continue;}
			/* see if need to seed comparison */
			if (mindex < 0){
				mindex = j2*n1 + j1;
				mini   = array[ mindex];
			}
			/* compare to see if current value is minimum */
			if ( array[j2*n1 + j1] < mini ){
				mindex = j2*n1 + j1;
				mini   = array[mindex];
			}
		}
	}
	return mindex ;
}
