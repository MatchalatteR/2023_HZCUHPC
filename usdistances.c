#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<sys/time.h>

double** distances; 

// Calculate sum of distance while combining different pivots. Complexity : O( n^2 )
double SumDistance(const int k, const int n, const int dim, double* coord, int* pivots){
    int i;
    // Rebuild coordinates. New coordinate of one point is its distance to each pivot.
    #pragma omp parallel for
    for(i=0; i<n; i++){
        int ki;
        for(ki=0; ki<k; ki++){
            distances[i][ki] = 0;
            int pivoti = pivots[ki];
            int j;
            for(j=0; j<dim; j++){
                distances[i][ki] += pow(coord[pivoti*dim + j] - coord[i*dim + j] ,2);
            }
            distances[i][ki] = sqrt(distances[i][ki]);
        }
    }

    // Calculate the sum of Chebyshev distance with rebuilt coordinates between every points
    double chebyshevSum = 0;
    for(i=0; i<n; i++){
        int j;
        for(j=i+1; j<n; j++){
            double chebyshev = 0;
            int ki;
            for(ki=0; ki<k; ki++){
                double dis = fabs(distances[i][ki] - distances[j][ki]);
                chebyshev = dis>chebyshev ? dis : chebyshev;
            }
            chebyshevSum += chebyshev;
        }
    }
    chebyshevSum*=2;
    return chebyshevSum;
}

// Recursive function Combination() : combine pivots and calculate the sum of distance while combining different pivots.
// ki  : current depth of the recursion
// k   : number of pivots
// n   : number of points
// dim : dimension of metric space
// M   : number of combinations to store
// coord  : coordinates of points
// pivots : indexes of pivots
// maxDistanceSum  : the largest M distance sum
// maxDisSumPivots : the top M pivots combinations
// minDistanceSum  : the smallest M distance sum
// minDisSumPivots : the bottom M pivots combinations
void Combination(int ki, int k, int n, int dim, double* coord, int* pivots, double* maxDistanceSum, int* maxDisSumPivots, double* minDistanceSum, int* minDisSumPivots){
    // Initialize stack and pivots
    int stack[k+1];

    int i;
    for (i=0; i<k; i++) {
        pivots[i] = i;
    }
    stack[0] = 0;
    int top = 0;

    while (top >= 0) {
        int ki = stack[top];
        top--;

        if (ki == k-1) {
            int j;
            for (j=pivots[ki-1]+1; j<n; j++) {
                pivots[ki] = j;
                double distanceSum = SumDistance(k, n, dim, coord, pivots);
                maxDistanceSum[0] = fmax(maxDistanceSum[0], distanceSum);
                minDistanceSum[0] = fmin(minDistanceSum[0], distanceSum);
                int kj;
                for (kj=0; kj<k; kj++) {
                    maxDisSumPivots[kj] = pivots[kj];
                    minDisSumPivots[kj] = pivots[kj];
                }
                for (kj=0; kj<k; kj++) {
                    maxDisSumPivots[kj] = pivots[kj];
                    minDisSumPivots[kj] = pivots[kj];
                }
            }
        } else {
            int j;
            for (j=pivots[ki-1]+1; j<n; j++) {
                pivots[ki] = j;
                top++;
                stack[top] = ki+1;
            }
        }
    }
}

int main(int argc, char* argv[]){
    // filename : input file namespace
    char* filename = (char*)"uniformvector-2dim-5h.txt";
    if( argc==2 ) {
        filename = argv[1];
    }  else if(argc != 1) {
        printf("Usage: ./pivot <filename>\n");
        return -1;
    }
    // M : number of combinations to store
    const int M = 1000;
    // dim : dimension of metric space
    int dim;
    // n : number of points
    int n;
    // k : number of pivots
    int k;

    // Read parameter
    FILE* file = fopen(filename, "r");
    if( file == NULL ) {
        printf("%s file not found.\n", filename);
        return -1;
    }
    fscanf(file, "%d", &dim);
    fscanf(file, "%d", &n);
    fscanf(file, "%d", &k);
    printf("dim = %d, n = %d, k = %d\n", dim, n, k);

    // Start timing
    struct timeval start;
    gettimeofday(&start, NULL);

    // Read Data
    double* coord = (double*)malloc(sizeof(double) * dim * n);
    int i;
    for(i=0; i<n; i++){
        int j;
        for(j=0; j<dim; j++){
            fscanf(file, "%lf", &coord[i*dim + j]);
        }
    }
    fclose(file);

    // maxDistanceSum : the largest M distance sum
    double* maxDistanceSum = (double*)malloc(sizeof(double) * (M+1));
    for(i=0; i<M; i++){
        maxDistanceSum[i] = 0;
    }
    // maxDisSumPivots : the top M pivots combinations
    int* maxDisSumPivots = (int*)malloc(sizeof(int) * k * (M+1));
    for(i=0; i<M; i++){
        int ki;
        for(ki=0; ki<k; ki++){
            maxDisSumPivots[i*k + ki] = 0;
        }
    }
    // minDistanceSum : the smallest M distance sum
    double* minDistanceSum = (double*)malloc(sizeof(double) * (M+1));
    for(i=0; i<M; i++){
        minDistanceSum[i] = __DBL_MAX__;
    }
    // minDisSumPivots : the bottom M pivots combinations
    int* minDisSumPivots = (int*)malloc(sizeof(int) * k * (M+1));
    for(i=0; i<M; i++){
        int ki;
        for(ki=0; ki<k; ki++){
            minDisSumPivots[i*k + ki] = 0;
        }
    }

    // temp : indexes of pivots with dummy array head
    int* temp = (int*)malloc(sizeof(int) * (k+1));
    temp[0] = -1;

distances = (double **) malloc(n * sizeof(double *));
    for(int i = 0; i < n; i++)
        distances[i] = (double *) malloc(k * sizeof(double));



    // Main loop. Combine different pivots with recursive function and evaluate them. Complexity : O( n^(k+2) )
    Combination(0, k, n, dim, coord, &temp[1], maxDistanceSum, maxDisSumPivots, minDistanceSum, minDisSumPivots);

    // End timing
    struct timeval end;
    gettimeofday (&end, NULL);
    printf("Using time : %f ms\n", (end.tv_sec-start.tv_sec)*1000.0+(end.tv_usec-start.tv_usec)/1000.0);

    // Store the result
    FILE* out = fopen("result.txt", "w");
    for(i=0; i<M; i++){
        int ki;
        for(ki=0; ki<k-1; ki++){
            fprintf(out, "%d ", maxDisSumPivots[i*k + ki]);
        }
        fprintf(out, "%d\n", maxDisSumPivots[i*k + k-1]);
    }
    for(i=0; i<M; i++){
        int ki;
        for(ki=0; ki<k-1; ki++){
            fprintf(out, "%d ", minDisSumPivots[i*k + ki]);
        }
        fprintf(out, "%d\n", minDisSumPivots[i*k + k-1]);
    }
    fclose(file);

    // Log
    int ki;
    printf("max : ");
    for(ki=0; ki<k; ki++){
        printf("%d ", maxDisSumPivots[ki]);
    }
    printf("%lf\n", maxDistanceSum[0]);
    printf("min : ");
    for(ki=0; ki<k; ki++){
        printf("%d ", minDisSumPivots[ki]);
    }
    printf("%lf\n", minDistanceSum[0]);

        for(int i = 0; i < n; i++)
        free(distances[i]);
    free(distances);
    return 0;
}