#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>  //All vectors have been declared as gsl vectors for consistency as some are required to
                             //are required to solve the tridiagonal matrix.

#define N 100                //Radial resolution of simulation
#define M 10000.0f           //Temporal resolution of simulation
#define k 20000000.0f        //Fuel rod constant
#define a 25.0f              //Radius of fuel rod
#define RC 100.0f            //Outer radius of simulation (where T is assumed to be fixed at 300K)

void set_initial_conditions(gsl_vector *r, gsl_vector *time, gsl_vector *T[])
{
    int i;

    for(i = 0; i < N; i++)
    {
        gsl_vector_set(r, i, i*RC / (N-1));   //Sets the discrete radial intervals between 0 and RC

        gsl_vector_set(T[0], i, 300);         //Sets the initial environmental temperature to 300K
        gsl_vector_set(T[1], i, 300);
        gsl_vector_set(T[2], i, 300);
    }

    gsl_vector_set(time, 0, 0);               //Sets the start time equal to 0
    gsl_vector_set(time, 1, 0);
    gsl_vector_set(time, 2, 0);
}

void solve_matrix(double duration, gsl_vector *r, gsl_vector *time, gsl_vector *T[])
{
    int i, j;
    double s, timestep;

    gsl_vector *source, *diag, *above_diag, *below_diag, *b;        //Gsl_vector format is required to use the
                                                                    //non symmetric tridiagonal matrix solver

    source = gsl_vector_alloc(N);               //Define the lengths of each gsl_vector
    diag = gsl_vector_alloc(N);
    b = gsl_vector_alloc(N);
    above_diag = gsl_vector_alloc(N-1);
    below_diag = gsl_vector_alloc(N-1);


    timestep = duration/M;
    s = k * timestep / (gsl_vector_get(r, 1)*gsl_vector_get(r, 1));    //s is a required variable for the tridiagonal matrix as seen
                                                                       //in the solution to this problem on page 45 of the booklet

    for(i = 0; i < N-1; i++)                                           //The solution is of the form Ax = b, where A is a matrix,
    {                                                                  //b is a vector and x is a vector of unknowns.
        gsl_vector_set(diag, i, 1 + 2*s);                          //This is the diagonal vector of matrix A
        gsl_vector_set(above_diag, i, -s - s/(2*(i+1)));           //This is the super-diagonal vector of matrix A
        gsl_vector_set(below_diag, i, -s + s/(2*(i+2)));           //This is the sub-diagonal vector of matrix A

    }

    gsl_vector_set(diag, N-1, 1 + 2*s);        //Sets the final diagonal element as the loop ends one short due to super/sub vectors being N-1 in length.
    gsl_vector_set(diag, 0, 1 + s + s/2);      //Heat cannot flow into position r[0], so an exception for this limit is made and approximated as r[1]

    for(j = 0; j < M; j++) //Loops until duration is reached
    {
        gsl_vector_set(time, 1, j*timestep);  //Jumps time forward by the timestep.

        for(i = 0; i < N; i++)
        {
            gsl_vector_set(T[1], i, gsl_vector_get(T[2], i));  //Redefines the current temperature T[1] as
                                                               //the previously stored future temperature T[2]
            if(gsl_vector_get(r, i) <= a)
            {
                gsl_vector_set(source, i, exp(-gsl_vector_get(time, 1)/100)/(a*a)); //Recalculates the heat source at the new time
            }
            else
            {
                gsl_vector_set(source, i, 0);  //When position r lies beyond the edge of the rod (a), heat source is 0.
            }

            gsl_vector_set(b, i, gsl_vector_get(T[1], i) + k * timestep * gsl_vector_get(source, i)); //Recalculates b based on the updated current
                                                                                                      //temperature T[1] and the updated source.

            if(i == N-1)        //An exception to b is added for the final coordinate due to the temperature at r = RC
            {                   //being fixed at the environmental constant 300K
                gsl_vector_set(b, i, (gsl_vector_get(T[1], i) + k * timestep * gsl_vector_get(source, i)) + 300*(s + s/(2*N)));
            }
        }
                    //The matrix is now set up in the format at the bottom of page 45, the new
                    //temperature T[2] is the unknown in Ax = b.

                    //the GSL scientific library has an inbuilt function that is able to solve
                    //such a tridiagonal matrix efficiently.

        gsl_linalg_solve_tridiag(diag, above_diag, below_diag, b, T[2]); //This calculates and sets the new temperature T[2]
    }
}


void plotdata()
{
    char command[4096];
    char *p = command;
    FILE *fp;
    int i, rc;


    fp = fopen("output.gp", "w");               //This will be the gnuplot script file
    if (fp == NULL) {
        printf("Could not open the script file for writing\n");  //Check to make sure file has been opened
        exit(2);
    }
    else
    {
        //Write the following commands to script file, this commands will be run in the gnuplot terminal
        //in order to plot the graphical output

        fprintf(fp, "set   autoscale                        # scale axes automatically\n");
        fprintf(fp, "set title \"How The Temperature Distribution of Nuclear Fuel Rod Varies Over Time\"\n");
        fprintf(fp, "set xlabel \"Radius(cm)\"\n");
        fprintf(fp, "set ylabel \"Temperature(K)\"\n");
        fprintf(fp, "set xr [0:100]\n");
        fprintf(fp, "plot \"output.dat\" using 1:2 title 'Time = 1 year' with lines, ");        //plot all the data on the same graph
        fprintf(fp, " \"output.dat\" using 1:3 title 'Time = 10 years' with lines, ");
        fprintf(fp, " \"output.dat\" using 1:4 title 'Time = 50 years' with lines, ");
        fprintf(fp, " \"output.dat\" using 1:5 title 'Time = 100 years' with lines ,");

        fprintf(fp, "\npause -1\n");
        fclose(fp);
    }

    /* Call gnuplot to use the above script file */
    p += sprintf(p, "cmd /K \"C:\\Program Files (x86)\\gnuplot\\bin\\gnuplot.exe\"  output.gp"); //IMPORTANT!!! THIS CALLS GNUPLOT WHEN
                                                                                                 //INSTALLED ON WINDOWS, IF GNUPLOT CANNOT BE
    printf("command: [%s]\n", command);                                                          //CALLED, CHANGE THIS PATH APPROPRIATELY
    rc = system(command);
    printf("command returned %d\n", rc);
};

void save_output(gsl_vector *T_final[], gsl_vector *r)
{
    int i, j;
    FILE* output = fopen("output.dat", "w+");

    if (output == NULL)
    {
        printf("Could not open output data file!\n");  //Check to make sure output file has been opened or created.
        exit(2);                                       //return error message in the event of a problem
    }
    else
    {
        for(i = 0; i < N; i++)
        {
            fprintf(output, "%lf ", gsl_vector_get(r, i)); //Writes the r coordinates to the output file

            for(j = 0; j < 4; j++)
            {
                fprintf(output, "%lf ",                     //Writes a column for each T_Final vector.
                        gsl_vector_get(T_final[j], i));
            }

            if(i != N-1)
            {
                fprintf(output, "\n");
            }
        }
        fclose(output);
    }
    printf("Output file has been generated.\n");
}

void store_data(gsl_vector *T[], gsl_vector *T_final[], gsl_vector *r, gsl_vector *time, int j)
{
    int i;
    for(i = 0; i < N ; i++)
    {
        gsl_vector_set(T_final[j], i, gsl_vector_get(T[2], i));      //Due to the iterative temperatures T[1] and T[2] needing to be
    }                                                                //re-used for each scenario, the results must be stored as a more
    set_initial_conditions(r, time, T);                              //permanent vector T_Final before the next simulation is run.
}

int main()
{
    int i;
    gsl_vector *r, *T[3], *time, *T_final[4];

    r = gsl_vector_alloc(N);                        //Define the lengths of each vector.
    time = gsl_vector_alloc(3);
    T[0] = gsl_vector_alloc(N);
    T[1] = gsl_vector_alloc(N);
    T[2] = gsl_vector_alloc(N);
    T_final[0] = gsl_vector_alloc(N);
    T_final[1] = gsl_vector_alloc(N);
    T_final[2] = gsl_vector_alloc(N);
    T_final[3] = gsl_vector_alloc(N);

    set_initial_conditions(r, time, T);

    solve_matrix(0.01, r, time, T);              //Finds the Temperature as a function of r after a 1 year simulation
                                                //using the method described on page 45 to find the temperature at the
    store_data(T, T_final, r, time, 0);         //next time interval, looped until duration is reached
                                            //Once calculated, this final heat distribution
                                            //vector is stored permanently as a T_final vector.
    solve_matrix(10.0, r, time, T);

    store_data(T, T_final, r, time, 1);              //Process is repeated to generate results for duration = 1/10/50/100 years

    solve_matrix(50.0, r, time, T);

    store_data(T, T_final, r, time, 2);

    solve_matrix(100.0, r, time, T);

    store_data(T, T_final, r, time, 3);

    save_output(T_final, r);        //Once all calculations have been made and each set of results stored,
                                    //prints the position vector r and each T_final result as columns to an output file output.dat

    while(i != 0 && i != 1)     //Allows the user the choice to plot the output data file using gnuplot, due to
                                //some incompatibility across different setups, if gnuplot cannot be called the program can still be run
                                //without producing the graphical output.
    {
        printf("\nDo you wish to plot this data?\n\n"
               "Enter 1 to plot.\n"
               "Enter 0 to continue.\n");
        scanf("%d", &i);
    }
    if(i == 1)
    {
        plotdata();     //The exciting bit!!!
    }

    return 0;
}
