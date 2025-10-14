/*  
ThinFilmMassSource_HB UDF

    This UDF computes the heat generated/lost by the thin film to add as a 
    source term to the ANSYS Fluent Simulation

    Developed By: Saaras Pakanati, Ayaaz Yasin - 2025

    At The University of Cincinnati Lab for Interfacial Dynamics
    


    
    Inputs******************Description*********************************************format

    DEFINE_SOURCE Function
    --------------------------------------------------------------------------------------
    lengthT                 Number of data points the thin film int. temp. provides double
*/

#include <udf.h>
#include <math.h>
#include <stdio.h>

extern double Sm_Global;
double Q_Global;

// DEFINE_SOURCE FUNCTION
//  THE MAIN FUNCTION.
//      This function contains the wall temperature curveFit code and in general is the function that communicates with Ansys.
DEFINE_SOURCE(ThinFilmMassSource_HB,c,t,dS,eqn)
{

    int iter  = N_ITER;				// Get the current iteration number
	int niter = iter + 1;           // Get the new iteration number
    
	// Cell position
	double pos[ND_ND], xPosVap, yPosVap;
	C_CENTROID(pos,c,t);			// Cell centroid location
	xPosVap = pos[0];				// x-position of cell
	yPosVap = pos[1];				// y-position of cell

    // Initialize the heat source term
    double Q;

    if (iter > 100) {
        // Start running the UDF at the 105th iteration.
        if (iter % 5 == 0.0) {

            // Start of reading wall temperatures from Ansys Domain
            long long xAddTemp = floor(yPosVap*1e12);
            long long xAdd = (long long)xAddTemp;


            // Define the maximum number of rows and columns
            int maxRows = 100000;
            int maxCols = 2;
            int rows = 0;
            int cols = 0;
            int i = 0;


            // Declare an array to store the data
            long long xLoc = -1;					    // Initialization the xLoc.

            // CALLING THE THIN FILM SOLUTION.
            double vol = C_VOLUME(c,t);			        // cell volume
            double Sm = Sm_Global;                      // Mass source from the thin film calculations

            // Declare an array to store the data
            double Ti_data[maxRows][maxCols];
            rows = 0;
            cols = 0;
            double Ti = -1;								// Initializing the Wall Temperature.


            // Read data from the file
            char temp_Ti[] = "Ti.txt";					// Temp filename
            char file_Ti[20];

            sprintf(file_Ti, "%d%s", iter, temp_Ti);	// Concatenating iteration number 
            FILE *fp_Ti = fopen(file_Ti, "r");			// Open File
            
            if (fp_Ti == NULL) {						// File Opening Error Check
                printf("Unable to open the file, %s.\n", file_Ti);
            }
            
            while (fscanf(fp_Ti, "%lf\t%lf\t%lf", &Ti_data[rows][0], &Ti_data[rows][1], &Ti_data[rows][2]) == 3) {
                rows++;
                
                if (rows >= maxRows) {
                    printf("Array size exceeded. Increase maxRows: %s\n", file_Ti);
                    break;
                }
            }
            fclose(fp_Ti);								// Close the file


            int lengthT = 20000;                        // # of cells reading.

            // Initializing interfacial temperature calculation related variables
            double yT[lengthT];
            double yM[lengthT];
            double yMadd[lengthT];
            double Tcount = 0.0;
            double Mcount = 0.0;
            

            i = 0;
            while (i < lengthT){
                yMadd[i] = Ti_data[i][2];  // Extract the Interface Temperatures
                Mcount = Mcount + yMadd[i];
                i = i + 1;                 // increment
            }

            i = 0;
            while (i < lengthT){
                yT[i] = Ti_data[i][1];   // Extract the Interface Temperatures
                yM[i] = Ti_data[i][2];   // Extract the Interface Temperatures
                Tcount = Tcount + ((yT[i] * yM[i]) / Mcount);
                i = i + 1;               // increment
            }
            double Ti_avg = Tcount;

            

            // Calculating Heat Loss/Gain.
            double cp = C_CP(c,t);                      // [J/kg-K] Specific Heat.
            double TV = C_T(c,t);			            // [K] Vapor Temperature.
            double TVsat = Ti_avg; 				        // [K] Vapor Saturation Temperature.
            
            if (Sm < 0.0) {
                Q = (-Sm * cp * (298.15 - TV));         // if evaporation rate is negative, remove mass at current cell temperature
            } else {
                Q = (-Sm * cp * (298.15 - TVsat));      // for zero or positive evaporation rate
            }

            // Save the heat energy for future iterations.
            Q_Global = Q;

            // Save Computed Curve Fit Coefficient Values and Mass Source Term to external file
            FILE *file_data = fopen("_data_ThinFilmMassSource_HB.txt", "a");
            if (file_data == NULL) {
                Message("Error: Could not open output file (ThinFilmMassSource_HB).\n");
            }
            fprintf(file_data, "%d\t%0.16f\t%0.16f\t%0.16f\t%0.16f\n", iter, Sm, Q, Ti_avg, Tcount);
            fclose(file_data);
        } else {
            Q = Q_Global;
        }
    } else if (iter == 1) {
        Q_Global = 0;
        Q = Q_Global;
    } else {
        Q = 0;
    }
    
    // Save Values for next iteration Interfacial TEMPERATURE (Ti)
    FILE *fp_nQ = fopen("_data_Q.txt", "a");
    if (fp_nQ == NULL) {
        Message("Error: Could not open output file, Q.\n");
    }
    fprintf(fp_nQ, "%d\t%0.16f\t%0.16f\n", iter, Q_Global, Q);
    fclose(fp_nQ);

    // Returning the source term to the Ansys Domain.
	dS[eqn] = Q;
	return Q;
}