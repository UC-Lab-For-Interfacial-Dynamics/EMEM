/*  
ThinFilm_MassSource.c

    This User Defined Function (UDF) computes the mass source term for the computed thin film solution.

    Developed By: Saaras Pakanati, Ayaaz Yasin, and Kishan Bellur - October 2025
    At The University of Cincinnati Lab for Interfacial Dynamics

----------------------------------------------------------------------------------------------------

    General Instructions:

        (1) Thermo-Physical properties are curve fit for H2 from 15K-35K.
            - Change curvefits, if not appropriate to use case.

        (2) Refer to the properties mentioned below and make the appropriate changes.
        
        (3) Kinetic Theory of Gases Model
            - Choose the appropriate model by uncommenting the code in the preCondition function.
              [Default - Wayner's model]

            - Available Models:
                1. Wayner's KTG model
                      Wayner, P.; Kao, Y.; LaCroix, L. The Interline Heat-Transfer Coefficient of an 
                      Evaporating Wetting Film. International Journal of Heat and Mass Transfer 1976, 
                      19, 487–492. https://doi.org/10.1016/0017-9310(76)90161-7. 
                2. Bellur's KTG model
                      Bellur, K.; Médici, E.F.; Hermanson, J.C.; Choi, C.K.; Allen, J.S. Modeling 
                      Liquid–Vapor Phase Change Experiments: Cryogenic Hydrogen and Methane. Colloids 
                      and Surfaces A: Physicochemical and Engineering Aspects 2023, 131932. 
                      https://doi.org/10.1016/j.colsurfa.2023.131932.

        (4) Adsorbed film stopping conditions.
            - The framework to terminate the thin film solution to a non-evaporating flat film exists
              for the conditions laid out below. Refer to the thinFilmSolve function for more information.
                    1. Adsorbed film height reaches a certain value.
                    2. Adsorbed film height's first derivative becomes zero.
                    3. Adsorbed film height's second derivative becomes zero.
                    4. Mass flux becomes zero.

                    [Default - Zero Mass flux]

----------------------------------------------------------------------------------------------------
    Description of Properties:
    
    Property****************Description*****************************************************units

    getProperties Function
    -----------------------------------------------------------------------------------------------
    KL                      Liquid Thermal Conductivity Curve Fit                           W/m*K
    KL                      Vapor Thermal Conductivity Curve Fit                            W/m*K
    PV                      Vapor Pressure                                                  MPa
    sigma                   Surface Tension                                                 N/m
    RL                      Liquid Density                                                  kg/m^3
    RV                      Vapor Density                                                   kg/m^3
    Hfg                     Latent Heat of Vaporization                                     J/kg
    mu                      Viscosity                                                       Pa*s

    preCondition Function
    -----------------------------------------------------------------------------------------------
    HaC                     Hamacker Constant                                               -
    zeta                    Under-relaxation factor                                         -
    R                       Radius                                                          m
    M                       Molar Mass                                                      g/mol
    mB                      Molecular Mass                                                  g

    evolution Function
    -----------------------------------------------------------------------------------------------
    R                       Radius                                                          m
    HaC                     Hamacker Constant                                               -
    gamma                   Temperature derivative of Surface Tension                       -

    thinFilmSolve Function
    -----------------------------------------------------------------------------------------------
    x                       Initial location                                                m
    L                       Interface length                                                m   
    dx                      Spatial discretization                                          m
    h0                      Initial height                                                  m
    h1                      First Height Derivative                                         m/m
    h2                      Second Height Derivative                                        1/m
    M_lastEV                Inlet Mass Flow                                                 kg/s
    R                       Radius                                                          m
    H_ad                    Adsorbed Film Height [Optional stopping condition]              m
   
    DEFINE_SOURCE Function
    -----------------------------------------------------------------------------------------------
    length                  Number of cells the wall temperature UDF provides               #
    Tv_main                 Saturation Temperature                                          K
    gammaTv                 Gamma Formulation for Knudsen Layer Reduction                   -

----------------------------------------------------------------------------------------------------

*/

#include <udf.h>
#include <math.h>
#include <stdio.h>

// Initialize global variable for mass source term to communicate current mass source to the ThinFilm_HeatBal code.
double Sm_Global;

// getProperties FUNCTION
//      This function provides the thermoPhysical properties for a given temperature.
//      - Curvefit generated from NIST Standard Reference Database.
//        Linstrom, P. NIST Chemistry WebBook, NIST Standard Reference Database 69, 1997. https://doi.org/10.18434/T4D303.
double getProperties(double T, double* KL, double* KV, double* PV, double* sigma, double* RL, double* RV, double* Hfg, double* mu)
{
    *KL     = -0.000124979154984 * pow(T, 2)        + 0.005335130035383 * T         + 0.047038794468053;
    *KV     = 0.000001479775739 * pow(T, 4)         - 0.000126505585113 * pow(T, 3) + 0.004045030999534 * pow(T, 2) - 0.055949424433763 * T     + 0.293374322089212;
    *PV     = 121405;
    *sigma  = 0.000000099194193 * pow(T, 3)         - 0.000005761321932 * pow(T, 2) - 0.000063021125245 * T         + 0.004748749018374;
    *RL     = -0.0044690793821 * pow(T, 3)          + 0.2482366185492 * pow(T, 2)   - 5.6135094695879 * T           + 119.8947608515605;
    *RV     = 0.004725701257276 * pow(T, 3)         + -0.262910551267955 * pow(T, 2)+ 5.109903042079692 * T         - 33.539922830008848;
    *Hfg    = 445571.00;
    *mu     = 1e-3 * (0.000000112504984 * pow(T, 4) + -0.000014093589989 * pow(T, 3)+ 0.000657338803965 * pow(T, 2) + -0.014156307035266 * T    + 0.128818861971111);
}

// getAlpha FUNCTION
//      Transition state theory based Accomodation coefficient.
//      - Nagayama, G.; Tsuruta, T. A General Expression for the Condensation Coefficient Based on Transition State Theory and Molecular 
//        Dynamics Simulation. Journal of Chemical Physics 2003, 118, 1392–1399. https://doi.org/10.1063/1.1528192.
double getAlpha(double RV, double RL)
{
    double ratio = RV/RL;               // "l" variable initialization
    double a = 1;                       // "alpha"  variable initialization

    // Error handling to prevent SIGSEV errors.
    // To debug: Ensure that the curvefit density values are correct.
    if (ratio < 0.0) {
        // Computing the Accommodation Coefficients
        double l    = -1 * pow(-1 * RV/RL, 1.0/3.0);
        double linv = -1 * pow(-1 * RL/RV, 1.0/3.0);
        a = (1.0 - l) * exp((-0.5) * (1.0 /(linv - 1.0)));
    }
    
    if (ratio >= 0.0) {
        // Computing the Accommodation Coefficients
        double l    = pow(RV/RL, 1.0/3.0);
        double linv = pow(RL/RV, 1.0/3.0);
        a = (1 - l) * exp((-0.5) * (1.0 /(linv - 1)));
    }

    // a = 0.58;                        // <---------- Uncomment to override TST based alpha calculations.

    return a;                           // Return the accommodation coefficient.
}


// preCondition Function
//      This function provides the MassFlux and Interfacial temperature at a given x location on the thin film.
double preCondition(double A, double B, double x, int iter, double h, double h1, double h2, double T_lastEv, double Tv_main, double gammaTv, double* Ti, double* mFlux) 
{ 
    int maxIter = 1e4;                              // Maximum Inner Itterations to find the Interfacial Temperature.
    int minIter = 30;                               // Minimum Inner Itterations to find the Interfacial Temperature.

    double Tw = B + x*A;                            // [K] Wall Temperature.
    double Tv = Tv_main * (1 - gammaTv);               // [K] Vapor Temperature.

    double T = T_lastEv;                            // [K] Interfacial Temperature.
    
    
    double a        = 10000;                        // [-] Initializing Accomodation Coefficint Value.
    double mflux    = 10000;                        // [kg/m^2-s] Initializing Mass Flux.

    for (int i = 0; i < maxIter; ++i)
    {
        // Constants
        double HaC  = 5.1100e-21;                   // Hamaker Constant
        double zeta = 1e-1;                         // Underrelaxation factor
        double Pd   = HaC / pow(h, 3);              // [kPa] Disjoining Pressure
        double R    = 0.005;                        // [m] Radius of the CryoStat.
        double pi   = 3.141592653589;		        // pi
        double M    = 2.01588e-3; 		            // [kg/mol], molar mass of H2
        double Rc   = 8.314462618153240;            // [J/K-mol], gas constant
        double kB   = 1.380649e-23;                 // [kg-m^2/s^2-K], Boltzmann's Constant
        double mB   = 3.34744749e-27;               // [kg], Molecular Mass of H2

        // Getting Fluid Properties
        double KL, KV, PV, sigma, RL, RV, Hfg, mu;
        getProperties(T, &KL, &KV, &PV, &sigma, &RL, &RV, &Hfg, &mu);
        double kappa    = (1.0 / ((R-h) * pow((1 + pow(h1, 2)), 0.5))) + (h2 / (pow((1 + pow(h1, 2)), 1.5)));       // [] Curvature.
        double Pc       = sigma * kappa;                                                                            // [kPa] Capillary Pressure.
        double Vl       = 1.0 / RV;                                                                                 // [m^3/kg] Specific Volume of Liquid.

        a = getAlpha(RV, RL);                                                               //  Accomodation Coefficint Value
        

        // Calculating Mass Flux: Wayner
        double waynerCoeff = (2 * a/(2 - a)) * pow((M / (2 * pi * Rc * T)), 0.5);           // Wayner Coefficient
        double term1  = waynerCoeff * (PV * M * Hfg / (Rc * Tv * T * 1000))*(T - Tv);       // Wayner Thermal Term
        double term2  = waynerCoeff * ((Vl * PV)/(Rc * T * 1000)) * (Pd + Pc);              // Wayner Pressure Term
        mflux  = term1  -  term2;                                                           // mass flux
        
        // Calculating Mass Flux: Beta-Form Bellur2023
        //double cRB = pow(2.0*kB*Tv/mB, 0.5);                                                // Most Probable Velocity
        //double alphaCoeffB = 2.0*a/(2.0-a);                                                 // Alpha term
        //double S_coeffB = (PV)/(cRB * 1000.0 * pow(pi, 0.5));                               // [s/m]
        //double beta = 1.0;                                                                  // Beta
        //double tempRatio = pow(Tv/T, 0.5);                                                  // Vapor and Interfacial Temperature Ratio

        //double term1B = 1.0;
        //double term2AB = 1.0 - (Tv/T);
        //double term2BB = RV * Hfg / (PV);
        //double term3AB = Tv/T;
        //double term3BB = RV/RL;
        //double term3CB = (Pc+Pd)/(PV);
        //double presRatio  = term1B + (term2AB * term2BB) + (term3AB * term3BB * term3CB);   // Pressure Ratio
        //mflux = alphaCoeffB * S_coeffB * ((beta * presRatio * tempRatio) - 1.0);            // mass flux

        // Interface Temperature
        double Ti_star = -((Hfg / KL) * mflux * (R - h) * log(R / (R - h))) + Tw;           // NEW Interface Temperature
        double Ti_new  = zeta * Ti_star + (1 - zeta) * T;                                   // Relaxation
        
        if (fabs(T - Ti_star) < 1e-4 && i > minIter) {
            break;                                                                          // Stopping Criterion
        }
        T = Ti_new;                                                                         // Variable shift for next iteration

    }

    // Returning Final Converged Interfacial Temperature and Mass Flux.
    *Ti = T;
    *mFlux = mflux;
    
}











// evolution FUNCTION
//      Computes the h''' value [Solving the Evolution Equation].
double evolution(double dx, double h, double h1, double h2, double k, double A, double B, double x, int iter, double Tv_main, double gammaTv, double *Hppp, double *T_lastEv, double *M_lastEv, double *mFlux) 
{
    
    // Initializing and Computing the Interfacial Temperature and Mass Flux.
    double Ti;
    preCondition(A, B, x, iter, h, h1, h2, *T_lastEv, Tv_main, gammaTv, &Ti, mFlux);

    // Getting Fluid Properties
    double KL, KV, PV, sigma, RL, RV, Hfg, mu;
    getProperties(Ti, &KL, &KV, &PV, &sigma, &RL, &RV, &Hfg, &mu);
    double alpha = getAlpha(RV, RL);
    double R = 0.0050;                                              // [m] Radius of the cryoStat.
    double HaC = 5.1100e-21;                                        // [-] Hamacker Constant

    double gamma = -1.645615366357070e-04;                          // [N/m-K] Surface Tension Temperature Gradient
    double pi  = 3.141592653589;                                    // [-] Pi - Constant

    // Derivatives 
    double dT_dx    = (Ti - *T_lastEv) / dx;                        // [K/m] Spatial Temperature Gradient
    double dPd_dx   =  -3 * HaC * (h1) / pow(h, 4);                 // [kPa/m] Disjoining Pressure Gradient
    double ds_dx    = gamma * dT_dx;                                // [N/m] Surface Tension Gradient

    // Computing Liquid Pressure Gradient.
    double GammaM = *M_lastEv - *mFlux * (pi * (pow(R, 2) - pow((R - h), 2)));
    double num1 = (8.0 * mu * GammaM) / (pi * RL * (pow(R,2) - pow((R - h), 2)));
    double num2 = (R - h) * ds_dx * (4.0 - 8.0 * log(R / (R - h)));
    double den1 = (pow(R,2) - pow((R - h), 2)) - 2.0 * pow(R, 2) + pow((R - h), 2) * (4.0 * log(R / (R - h)) - 2.0);
    double dPl_dx = (num1 - num2) / den1;

    // Saving Interfacial Temperature and Mass Flux.
    *T_lastEv = Ti;
    *M_lastEv = GammaM;

    // THE EVOLUTION EQUATION.
    *Hppp = (3.0 * h1 * pow(h2, 2) / alpha) + (h1*h2/(R-h)) - (h1*alpha/pow((R-h), 2))  -  ((alpha/(R-h)) + h2)*(gamma/sigma)*dT_dx - ((pow(alpha, 1.5))/sigma)*(dPl_dx + dPd_dx);

}


// thinFilmSolve FUNCTION
//      This function contains the RangeKutta Code to solve the thinFilm
double thinFilmSolve(double A, double B, double vol, int iter, double Tv_main, double gammaTv) 
{
        
    int initer = 0;

    // **************** The Spatial Specifications of the Thin Film ****************
    double x    = 0.0;                                           // [m] Initial Location
    double L    = 200e-6;                                        // [m] Length of the Interface.
    double dx   = 1e-8;                                          // [m] Discretization Length.
    double n    = ((L - x) / dx) + 1.0;                          // [#] Number of steps.

    double Tw   = B + x*A;                                       // [K] Wall Temperature
    double Tv   = Tv_main * (1 - gammaTv);                       // [K] Vapor Temperature

    double Ti_last = (Tv + Tw) / 2.0;                            // [K] Initial Interface Temperature

    double h0   = 10e-6;                                         // [m] Initial Height
    double h1   = -1.265280644375428e-01;                        // [-] Initial Height First Derivative
    double h2   = 7.811371253212997e+02;                         // [1/m] Initial Height Second Derivative
    double h3   = 1;                                             // [1/m^2] Initial Height Third Derivative [Not Required, however helps with debugging]

    double H[]  = {h0, h1, h2};                                  // H vector definition. H = [h, h', h'']
    double Hppp;                                                 // Initializing the third derivative.

    double T_lastEv = Ti_last;                                   // Initializing the interfacial Temperature.
    double Ti, mFlux, M_lastEv;                                  // Additional Initializations.

    // Obtaining the interfacial temperature and mass flux.
    preCondition(A, B, x, iter, h0, h1, h2, T_lastEv, Tv_main, gammaTv, &Ti, &mFlux);
    T_lastEv = Ti;
    Ti_last = Ti;

    // Initializing Mass flow.
    M_lastEv = 1e-8;                                             // Inlet Mass Flow
    double M_last = M_lastEv;

    // Mass Source Term
    double Sm   = 0;
    double R    = 0.0050;                                        // [m] Radius of the cryoStat.
    double pi   = 3.141592653589;                                // [-] Pi - Constant

    // Adsorbed Region.
    double H_ad = 100.0e-9;                                      // [m] Adsorbed Film Height
    double H_check[] = {1.0, 1.0, 1.0};                          // [m] Initializing Condition for adsorbed height check       

    // THE Runge-Kutta Loop
    for (int  j = 0; j < n; j++) {
        
        // Check for zero mass flux stopping condition to define the adsorbed region.
        // Can change to "if (H_check[0] > 1e-10) {" to terminate the solution at a specific height (in example 1e-10 m).
        if (mFlux > 1e-5) {
            // Initializer Step
            evolution(dx, H[0], H[1], H[2], 0.0, A, B, x, iter, Tv_main, gammaTv, &Hppp, &T_lastEv, &M_lastEv, &mFlux);
            T_lastEv = Ti_last;
            M_lastEv = M_last;
            double k1[] = {H[1], H[2], Hppp};

            // Predictor Step
            double H1[] = {H[0] + 0.5*dx*k1[0], H[1] + 0.5*dx*k1[1], H[2] + 0.5*dx*k1[2]};
            evolution(dx, H[0], H[1], H[2], 0.0, A, B, x+(dx/2), iter, Tv_main, gammaTv, &Hppp, &T_lastEv, &M_lastEv, &mFlux);
            T_lastEv = Ti_last;
            M_lastEv = M_last;
            double k2[] = {H1[1], H1[2], Hppp};

            // Corrector Step
            double H2[] = {H[0] + 0.5*dx*k2[0], H[1] + 0.5*dx*k2[1], H[2] + 0.5*dx*k2[2]};
            evolution(dx, H[0], H[1], H[2], 0.0, A, B, x+(dx/2), iter, Tv_main, gammaTv, &Hppp, &T_lastEv, &M_lastEv, &mFlux);
            T_lastEv = Ti_last;
            M_lastEv = M_last;
            double k3[] = {H2[1], H2[2], Hppp};

            // Extrapolation Step
            double H3[] = {H[0] + dx*k3[0], H[1] + dx*k3[1], H[2] + dx*k3[2]};
            evolution(dx, H[0], H[1], H[2], 0.0, A, B, x+dx, iter, Tv_main, gammaTv, &Hppp, &T_lastEv, &M_lastEv, &mFlux);
            double k4[] = {H3[1], H3[2], Hppp};

            // Obtaining the H values.
            double kc[] = {k1[0] + 2.0*k2[0] + 2.0*k3[0] + k4[0],     k1[1] + 2.0*k2[1] + 2.0*k3[1] + k4[1],      k1[2] + 2.0*k2[2] + 2.0*k3[2] + k4[2]};

            // Integrating step
            for (int idx = 0; idx < 3; idx++) {
                H[idx] = H[idx] + kc[idx] * dx / 6.0;
                H_check[idx] = H[idx];
            }
        } else {
            // If the film is defined to be in the adsorbed region, make the height the same as the previous iteration, mass flux zero and, mass flow zero.
            for (int idx = 0; idx < 3; idx++) {
                H[idx] = H_check[idx];
                M_lastEv = 0.0;
                mFlux = 0.0;
            }
        }

        // Save Values for next iteration Interfacial TEMPERATURE (Ti)
        char temp_nTi[] = "Ti.txt";
        char file_nTi[20];
        sprintf(file_nTi, "%d%s", iter, temp_nTi);	                                            // Concatenating iteration number
        FILE *fp_nTi = fopen(file_nTi, "a");
        if (fp_nTi == NULL) {
            Message("Error: Could not open output file, %s.\n", file_nTi);
        }
        // Message("Check \n", file_nTi);
        fprintf(fp_nTi, "%d\t%0.16f\t%0.16f\n", iter, Ti_last, mFlux);
        fclose(fp_nTi);

        // Saving Interfacial Temperature and Mass Flux after a completed RK step.
        Ti_last = T_lastEv;
        M_last = M_lastEv;

        // Extracting thermoPhysical properties for reference.
        double KL, KV, PV, sigma, RL, RV, Hfg, mu;
        getProperties(Ti_last, &KL, &KV, &PV, &sigma, &RL, &RV, &Hfg, &mu);
        double alpha = getAlpha(RV, RL);

        // Update Mass Source Term.
        Sm = Sm + mFlux * (pi * (2 * (R - H[0]) + H[1] * dx) * dx) / vol;               // Volumetric Mass Source Term

        // Update spatial location for the next iteration
        x = x + dx;

        // Save Read Values to external file
        FILE *file_dataG = fopen("_data_Global.txt", "a");
        if (file_dataG == NULL) {
            Message("Error: Could not open output file (Global).\n");
        }
        fprintf(file_dataG, "%d\t%d\t%0.16f\t%0.16f\t%0.16f\t%0.16f\t%0.16f\t%0.16f\t%0.16f\t%0.16f\t%0.16f\t%0.16f\n", iter, j, x, H[0], H[1], H[2], alpha, vol, A, B, Ti_last, mFlux);
        fclose(file_dataG);
    }

    return Sm;
}

















// DEFINE_SOURCE FUNCTION
//  THE MAIN FUNCTION.
//      This function contains the wall temperature curve fit code and in general is the function that communicates with Ansys.
DEFINE_SOURCE(ThinFilmMassSource,c,t,dS,eqn)
{

    int iter  = N_ITER;				// Get the current iteration number
	int niter = iter + 1;           // Calculate next iteration number
    
	// Cell position
	double pos[ND_ND], xPosVap, yPosVap;
	C_CENTROID(pos,c,t);			// Cell centroid location
	xPosVap = pos[0];				// x-position of cell
	yPosVap = pos[1];				// y-position of cell

    // Start of reading wall temperatures from Ansys Domain
	long long xAddTemp = floor(yPosVap*1e12);
	long long xAdd = (long long)xAddTemp;


	// Define the maximum number of rows and columns
    int maxRows = 100000;
    int maxCols = 2;
	int rows    = 0;
    int cols    = 0;
	int i       = 0;


    // Declare an array to store the data
    double Tw_data[maxRows][maxCols];
	double Tw = -1;					            // Initializing the Wall Temperature.
    long long xLoc = -1;			            // Initialization the xLoc.


    // Read data from the file
	char temp_Tw[] = "Tw.txt";		            // Wall Temperature filename
	char file_Tw[20];

	sprintf(file_Tw, "%d%s", iter, temp_Tw);	// Concatenating iteration number 
    FILE *fp_Tw = fopen(file_Tw, "r");			// Open File
    
    if (fp_Tw == NULL) {						// File Opening Error Check
        printf("Unable to open the file, %s.\n", file_Tw);
    }
    
    while (fscanf(fp_Tw, "%lf\t%lf", &Tw_data[rows][0], &Tw_data[rows][1]) == 2) {
        rows++;
        
        if (rows >= maxRows) {
            printf("Array size exceeded. Increase maxRows: %s\n", file_Tw);
            break;
        }
    }
    fclose(fp_Tw);								// Close the file
	// End of reading wall temperatures from Ansys Domain

    // Wall temperature CurveFit section
    double Sigma_X = 0;                         // n
    double Sigma_X1 = 0;                        // sum X
    double Sigma_X2 = 0;                        // sum X^2
    double Sigma_XY = 0;                        // sum X*Y
    double Sigma_Y = 0;                         // sum Y

    int length = 19;                            // # of cells reading.

    // Variable Initialization
    double x[length];
    double y[length];

    for (int i = 0; i < length; ++i)
    {
        x[i] = (Tw_data[i][0] / 1e12) - (Tw_data[0][0] / 1e12);            // Extract the Locations.
        y[i] = Tw_data[i][1];                                               // Extract the Temperatures.
    }

    for (int i = 0; i < length; ++i)
    {
        Sigma_X += 1;                           // Sigma x^0
        Sigma_X1 += x[i];                       // Sigma x^1
        Sigma_X2 += x[i] * x[i];                // Sigma x^2
        Sigma_XY += y[i] * x[i];                // Sigma xy
        Sigma_Y += y[i];                        // Sigma y
    }

    double M1[2][2] = { {Sigma_X2, Sigma_X1}, {Sigma_X1, Sigma_X} };
    double M2[2] = { Sigma_XY, Sigma_Y };
    

    // CALLING THE THIN FILM SOLUTION.
    double vol      = C_VOLUME(c,t);			// get cell volume from CFD
    double Sm       = 0.0;                      // Initializing the mass source term variable
    double Tv_main  = 21.0;          			// [K] Vapor Temperature.

    // Gamma Formulation
    //      A. Yasin and K. Bellur, “Computational modeling of evaporation without tuning coefficients,” Applied Thermal Engineering, vol. 276, 
    //      p. 126807, Oct. 2025, doi: 10.1016/j.applthermaleng.2025.126807.
    double gammaTv = 0.00032;
	
    // Numerical Treatment for convergence
    if (iter > 100) {
        // Start running the UDF at the 105th iteration.
        if (iter % 5 == 0.0) {
            // Recalculate the thin film every 5 iterations
            Sm = thinFilmSolve(A, B, vol, iter, Tv_main, gammaTv);
            Sm_Global = Sm;
        } else {
            // If on a non-5th iteration, maintain the same value as the previous iteration.
            Sm = Sm_Global;
        }
    } else if (iter == 1){
        // If before the 105th iteration, maintain a zero output. i.e., no microscale evaporation effect.
        Sm_Global = 0;
        Sm = Sm_Global;
    } else {
        Sm = 0.0;
    }

    // Save Values for next iteration Interfacial TEMPERATURE (Ti)
    FILE *fp_nSm = fopen("_data_Sm.txt", "a");
    if (fp_nSm == NULL) {
        Message("Error: Could not open output file, Sm.\n");
    }
    fprintf(fp_nSm, "%d\t%0.16f\t%0.16f\n", iter, Sm_Global, Sm);
    fclose(fp_nSm);


    // Returning the source term to the Ansys Domain.
	dS[eqn] = Sm;
	return Sm;
}