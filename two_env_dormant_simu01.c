#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>
#include <math.h>
#include <time.h>

// numerical calculation of extinction probability in system of age-structured cells in two repeating environments with division, death, and switching to phenotype that does not divide nor die 
// compile this file in UNIX command line with 'gcc -o two_env_dormancy_simu01.o two_env_dormancy_simu01.c -O3 -I/usr/include -L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas -lm' 
// run this program in command line with flags to set parameters, for example a reasonable choice at extinction resonance is
// './two_env_dormancy_simu01.o -L 2000 -n 20 -S 2.3 -T 2.3 -W 1.0 -w 0.55555 -X 1.0 -x 0.55555 -Y 1.0 -y 0.00001 -Z 1000.0 -z 500.0 -O 1.0 -o 0.0 -P 1.0 -p 10.0 -Q 1.0 -q 10.0 -R 1.0 -r 0.0 > test_results.dat'
// All the flags are explained below in the main function

// NOTE that the variable and array names ending on env1 and env2 are due to historical reasons swapped in their meaning, e.g. a_birth_env1 corresponds to the 2nd environment. The output labels the environments correctly though

// overview of the parameters of the stochastic process
// the 4 involved processes are division, death, switch from ND state to D state, and switch from D state to ND state
// there are 2 alternating environments and each process can change their parameters between the 2 environments
// each process in each environment is determined by two parameters of a Gamma distribution using the alpha (a) and beta (b) convention, such that the mean event time is alpha/beta
// tauND refers to the switch from ND state to the D state, and tauD refers to the switch from the D state to the ND state
// only in the ND state cells can divide or die. in the D state they can neither divide nor die 
double a_birth_env1 = 1.0;
double b_birth_env1 = 3.0;
double a_birth_env2 = 1.0;
double b_birth_env2 = 1.0;
double a_death_env1 = 1.0;
double b_death_env1 = 1.0;
double a_death_env2 = 1.0;
double b_death_env2 = 0.5;
double a_tauND_env1 = 1.0;
double b_tauND_env1 = 1.0;
double a_tauND_env2 = 1.0;
double b_tauND_env2 = 1.0;
double a_tauD_env1 = 1.0;
double b_tauD_env1 = 1.0;
double a_tauD_env2 = 1.0;
double b_tauD_env2 = 1.0;

// this parameters is needed to keep track of old cells that are older than an environmnet. All cells older than environment length * intervalextension are approximated by cells of that maximum age.
int intervalextension = 40;

// returns the birth rate for the two environments
double birth_rate(double age,int env){
    if (env == 1){
	    double rate = gsl_ran_gamma_pdf(age,a_birth_env1,1.0/b_birth_env1)/gsl_cdf_gamma_Q(age,a_birth_env1,1.0/b_birth_env1);
	    if (isnan(rate) || rate > b_birth_env1){
	    	return b_birth_env1;
	    }
	    else{
        	return rate;
	    }
    }
    else{
	    double rate = gsl_ran_gamma_pdf(age,a_birth_env2,1.0/b_birth_env2)/gsl_cdf_gamma_Q(age,a_birth_env2,1.0/b_birth_env2);
	    if (isnan(rate) || rate > b_birth_env2){
	    	return b_birth_env2;
	    }
	    else{
        	return rate;
    	    }
    }
}

// returns the death rate for the two environments
double death_rate(double age,int env){
    if (env == 1){
	    double rate = gsl_ran_gamma_pdf(age,a_death_env1,1.0/b_death_env1)/gsl_cdf_gamma_Q(age,a_death_env1,1.0/b_death_env1);
	    if (isnan(rate) || rate > b_death_env1){
	    	return b_death_env1;
	    }
	    else{
        	return rate;
	    }
    }
    else{
	    double rate = gsl_ran_gamma_pdf(age,a_death_env2,1.0/b_death_env2)/gsl_cdf_gamma_Q(age,a_death_env2,1.0/b_death_env2);
	    if (isnan(rate) || rate > b_death_env2){
	    	return b_death_env2;
	    }
	    else{
        	return rate;
	    }
    }
}

// returns the rate to switch from state ND to state D for the two environments
double tauND_rate(double age,int env){
    if (env == 1){
	    double rate = gsl_ran_gamma_pdf(age,a_tauND_env1,1.0/b_tauND_env1)/gsl_cdf_gamma_Q(age,a_tauND_env1,1.0/b_tauND_env1);
	    if (isnan(rate) || rate > b_tauND_env1){
	    	return b_tauND_env1;
	    }
	    else{
        	return rate;
	    }
    }
    else{
	    double rate = gsl_ran_gamma_pdf(age,a_tauND_env2,1.0/b_tauND_env2)/gsl_cdf_gamma_Q(age,a_tauND_env2,1.0/b_tauND_env2);
	    if (isnan(rate) || rate > b_tauND_env2){
	    	return b_tauND_env2;
	    }
	    else{
        	return rate;
	    }
    }
}

// returns the rate to switch from state D to state ND for the two environments
double tauD_rate(double age,int env){
    if (env == 1){
	    double rate = gsl_ran_gamma_pdf(age,a_tauD_env1,1.0/b_tauD_env1)/gsl_cdf_gamma_Q(age,a_tauD_env1,1.0/b_tauD_env1);
	    if (isnan(rate) || rate > b_tauD_env1){
	    	return b_tauD_env1;
	    }
	    else{
        	return rate;
	    }
    }
    else{
	    double rate = gsl_ran_gamma_pdf(age,a_tauD_env2,1.0/b_tauD_env2)/gsl_cdf_gamma_Q(age,a_tauD_env2,1.0/b_tauD_env2);
	    if (isnan(rate) || rate > b_tauD_env2){
	    	return b_tauD_env2;
	    }
	    else{
        	return rate;
	    }
    }
}
// calculates the integral (Newton sum) from age 0 to end over rates depending on species and environment
double hazards(double end,int step_no,int env,int species){
    if (end == 0){
        return 0.0;
    }
    else{
        double hazards_return = 0.0;
        int i = 0;
        double step_size = end/((double) step_no);
	if (species == 0){ // non-dormant species
            for (i=0;i<step_no;i++){
                hazards_return += step_size*(birth_rate(i*step_size,env)+death_rate(i*step_size,env)+tauND_rate(i*step_size,env));
            }
	}
	if (species == 1){ // dormant species
            for (i=0;i<step_no;i++){
                hazards_return += step_size*tauD_rate(i*step_size,env);
            }
	}
        return hazards_return;
    }
}

// calculates a G0 type function, see theory. This is also used to calculate a H0 function 
void G0(int step_no,double g0[2][step_no],double new_g0[2][step_no],double hT[2][intervalextension*step_no],double time_interval,int env,double hazards_speedup[2][2][(intervalextension+1)*step_no],double birth_speedup[2][(intervalextension+1)*step_no],double death_speedup[2][(intervalextension+1)*step_no],double tauND_speedup[2][(intervalextension+1)*step_no],double tauD_speedup[2][(intervalextension+1)*step_no]){
    int i = 0;
    int j = 0;
    double integral_sum_ND = 0.0;
    double integral_sum_D = 0.0;
    for (i=0;i<step_no;i++){
        integral_sum_ND = 0.0;
        integral_sum_D = 0.0;
        for (j=0;j<i;j++){
            integral_sum_ND += exp(-hazards_speedup[0][env-1][j])*(birth_speedup[env-1][j]*(g0[0][i-1-j]*g0[0][i-1-j])+death_speedup[env-1][j]+tauND_speedup[env-1][j]*g0[1][i-1-j]);
            integral_sum_D += exp(-hazards_speedup[1][env-1][j])*(tauD_speedup[env-1][j]*g0[0][i-1-j]);
        }
        new_g0[0][i] = exp(-hazards_speedup[0][env-1][i])*hT[0][i]+(time_interval/((double)step_no))*integral_sum_ND;
        new_g0[1][i] = exp(-hazards_speedup[1][env-1][i])*hT[1][i]+(time_interval/((double)step_no))*integral_sum_D;
        if (new_g0[0][i] > 1.0){
            new_g0[0][i] = 1;
        }
        if (new_g0[1][i] > 1.0){
            new_g0[1][i] = 1;
        }
    }
    return;
}

// calculates a GT type function, see theory. This is also used to calculate a HT function
void GT(int step_no,double g0[2][step_no],double gT[2][intervalextension*step_no],double hT[2][intervalextension*step_no],double time_interval,int env,double hazards_speedup[2][2][(intervalextension+1)*step_no],double birth_speedup[2][(intervalextension+1)*step_no],double death_speedup[2][(intervalextension+1)*step_no],double tauND_speedup[2][(intervalextension+1)*step_no],double tauD_speedup[2][(intervalextension+1)*step_no]){
    
    int j = 0;
    int k = 0;
	gT[0][0] = g0[0][step_no-1];
    gT[1][0] = g0[1][step_no-1];
    double integral_sum_ND = 0;
    double integral_sum_D = 0;
    for (j=1;j<intervalextension*step_no;j++){
        integral_sum_ND = 0;
        integral_sum_D = 0;
        for (k=0;k<step_no;k++){
            integral_sum_ND += exp(-hazards_speedup[0][env-1][j+k]+hazards_speedup[0][env-1][j])*(birth_speedup[env-1][j+k]*(g0[0][step_no-1-k]*g0[0][step_no-1-k])+death_speedup[env-1][j+k]+tauND_speedup[env-1][j+k]*g0[1][step_no-1-k]);
            integral_sum_D += exp(-hazards_speedup[1][env-1][j+k]+hazards_speedup[1][env-1][j])*(tauD_speedup[env-1][j+k]*g0[0][step_no-1-k]);
		}
        if (j-1 < step_no*(intervalextension-1)){
            gT[0][j] = exp(-hazards_speedup[0][env-1][step_no-1+j]+hazards_speedup[0][env-1][j])*hT[0][step_no-1+j]+(time_interval/((double)step_no))*integral_sum_ND;
            gT[1][j] = exp(-hazards_speedup[1][env-1][step_no-1+j]+hazards_speedup[1][env-1][j])*hT[1][step_no-1+j]+(time_interval/((double)step_no))*integral_sum_D;
        }
        else{
            gT[0][j] = exp(-hazards_speedup[0][env-1][step_no-1+j]+hazards_speedup[0][env-1][j])*hT[0][intervalextension*step_no-1]+(time_interval/((double)step_no))*integral_sum_ND;
            gT[1][j] = exp(-hazards_speedup[1][env-1][step_no-1+j]+hazards_speedup[1][env-1][j])*hT[1][intervalextension*step_no-1]+(time_interval/((double)step_no))*integral_sum_D;
        }
        if (gT[0][j] > 1){
            gT[0][j] = 1;
        }
        if (gT[1][j] > 1){
            gT[1][j] = 1;
        }
    }
    return;
}

int main(int argc, char *argv[]){

    time_t seconds1 = time(NULL); // used to record how long the calculation took, start time
    int length = 200; // total number of time steps to reach finite time T
    double time1 = 1.0; // length of time in environment 1
    double time2 = 1.0; // length of time in environment 2
    double accuracy = 0.00001; // accurancy for convergence of G(T,0) and H(T,0) equations
    double init_hT = 0.0; // initial guess for hT
    double init_gT = 0.0; // initial guess for gT
    int max_iterations = 200; // upper bound for the number of iterations for the convergence of G(T,0) and H(T,0)
    
    // flags for changing parameters from command line
    int ch;
    while ((ch = getopt(argc, argv, "L:S:T:I:J:A:M:R:W:w:X:x:Y:y:Z:z:O:o:P:p:Q:q:R:r:n:")) != -1) {
        switch (ch) {
            case 'L': // number of points per environment
            length=atoi(optarg);
            break;
            case 'S': // length of 2nd environment
            time1=strtod(optarg, NULL);
            break;
            case 'T': // length of 1st environment
            time2=strtod(optarg, NULL);
            break;
            case 'I': // initial guess for extinction probability if started in 2nd environment
            init_hT=strtod(optarg, NULL);
            break;
            case 'J': // initial guess for extinction probability if started in 1st environment
            init_gT=strtod(optarg, NULL);
            break;
            case 'A': // accuracy threshold determines when G0 and G1 are considered to have converged and iterations ares stopped
            accuracy=strtod(optarg, NULL);
            break;
            case 'M': // maximum number of iterations
            max_iterations=atoi(optarg);
            break;
            case 'W': // alpha parameter for division in 2nd environment
            a_birth_env1=strtod(optarg, NULL);
            break;
            case 'w': // beta parameter for division in 2nd environment
            b_birth_env1=strtod(optarg, NULL);
            break;
            case 'X': // alpha parameter for division in  1st environment
            a_birth_env2=strtod(optarg, NULL);
            break;
            case 'x': // beta parameter for division in 1st environment
            b_birth_env2=strtod(optarg, NULL);
            break;
            case 'Y': // alpha parameter for death in 2nd environment
            a_death_env1=strtod(optarg, NULL);
            break;
            case 'y': // beta parameter for death in 2nd environment
            b_death_env1=strtod(optarg, NULL);
            break; 
            case 'Z': // alpha parameter for death in 1st environment
            a_death_env2=strtod(optarg, NULL);
            break;
            case 'z': // beta parameter for death in 1st environment
            b_death_env2=strtod(optarg, NULL);
            break;
            case 'O': // alpha parameter for switching to dormancy in 2nd environment
            a_tauND_env1=strtod(optarg, NULL);
            break;
            case 'o': // beta parameter for switching to dormancy in 2nd environment
            b_tauND_env1=strtod(optarg, NULL);
            break;
            case 'P': // alpha parameter for switching to dormancy in 1st environment
            a_tauND_env2=strtod(optarg, NULL);
            break;
            case 'p': // beta parameter for switching to dormancy in 1st environment
            b_tauND_env2=strtod(optarg, NULL);
            break;
            case 'Q': // alpha parameter for waking up in 2nd environment
            a_tauD_env1=strtod(optarg, NULL);
            break;
            case 'q': // beta parameter for waking up in 2nd environment
            b_tauD_env1=strtod(optarg, NULL);
            break;
            case 'R': // alpha parameter for waking up in 1st environment
            a_tauD_env2=strtod(optarg, NULL);
            break;
            case 'r': // alpha parameter for wkaing up in 1st environment
            b_tauD_env2=strtod(optarg, NULL);
            break;
	    case 'n': // the maximum considered age is n*time1 or n*time2
            intervalextension=atoi(optarg);
            break;
            default:
            break;
        }
    }
    
    // results are printed to standard out
    printf("# simulation of structured birth-death process in two environments with dormancy\n");
    time_t t;
    t = time(NULL);
    struct tm tm = *localtime(&t);
    printf("# Date of simulation start (DD/MM/YYYY): %d/%d/%d\n", tm.tm_mday, tm.tm_mon+1, tm.tm_year+1900);
    printf("# length of G0 and H0 arrays = %d \n",length);
    printf("# time length of first interval = %f\n",time1);
    printf("# time length of second interval = %f\n",time2);
    printf("# accuracy = %f\n",accuracy);
    printf("# max_iterations = %d\n",max_iterations);
    printf("# interval extensions for G(T,zeta) and H(T,zeta), %i*T\n",intervalextension);
    printf("# initialisation of G0 = %f\n",init_gT);
    printf("# initialisation of H0 = %f\n",init_hT);
    printf("# Gamma function parameters for environment 1: birth = (%f,%f), death = (%f,%f)\n",a_birth_env2,b_birth_env2,a_death_env2,b_death_env2);
    printf("# Gamma function parameters for environment 2: birth = (%f,%f), death = (%f,%f)\n",a_birth_env1,b_birth_env1,a_death_env1,b_death_env1);
    printf("# Gamma function parameters for environment 1: tauND = (%f,%f), tauD = (%f,%f)\n",a_tauND_env2,b_tauND_env2,a_tauD_env2,b_tauD_env2);
    printf("# Gamma function parameters for environment 2: tauND = (%f,%f), tauD = (%f,%f)\n",a_tauND_env1,b_tauND_env1,a_tauD_env1,b_tauD_env1);
    fflush(stdout);
    // initialsation of function g0, gT, h0, hT and creation of lookup tables (ending on '_speedup') for rates and their integrals
    // the first argument in the array below is for non-dormant (entry 0) and dormant (entry 1) cells
    double hazards_speedup[2][2][(intervalextension+1)*length];
    double birth_speedup[2][(intervalextension+1)*length];
    double death_speedup[2][(intervalextension+1)*length];
    double tauND_speedup[2][(intervalextension+1)*length];
    double tauD_speedup[2][(intervalextension+1)*length];
    double g0[2][length];
    double gT[2][intervalextension*length];
    double h0[2][length];
    double hT[2][intervalextension*length];
    double inter_g0[2][length];
    double out_g0[2][length];
    double inter_h0[2][length];
    double out_h0[2][length];
    int k = 0;
    int l = 0;
    for (k=0;k<length;k++){
		for (l=0;l<2;l++){
        	g0[l][k] = init_gT;
        	gT[l][k] = init_gT;
        	h0[l][k] = init_hT;
        	hT[l][k] = init_hT;
        	inter_g0[l][k] = init_gT;
        	inter_h0[l][k] = init_hT;
        	out_g0[l][k] = init_gT;
        	out_h0[l][k] = init_hT;
			hazards_speedup[l][0][k] = hazards(k*time1/((double)length),length,1,l);
			hazards_speedup[l][1][k] = hazards(k*time2/((double)length),length,2,l);
			birth_speedup[l][k] = birth_rate(k*time1/((double)length),l+1);
			death_speedup[l][k] = death_rate(k*time1/((double)length),l+1);
			tauND_speedup[l][k] = tauND_rate(k*time1/((double)length),l+1);
			tauD_speedup[l][k] = tauD_rate(k*time1/((double)length),l+1);
		}
	}
    for (k=length;k<intervalextension*length;k++){
		for (l=0;l<2;l++){
        	gT[l][k] = init_gT;
        	hT[l][k] = init_hT;
			hazards_speedup[l][0][k] = hazards(k*time1/((double)length),length,1,l);
			hazards_speedup[l][1][k] = hazards(k*time2/((double)length),length,2,l);
			birth_speedup[l][k] = birth_rate(k*time1/((double)length),l+1);
			death_speedup[l][k] = death_rate(k*time1/((double)length),l+1);
			tauND_speedup[l][k] = tauND_rate(k*time1/((double)length),l+1);
			tauD_speedup[l][k] = tauD_rate(k*time1/((double)length),l+1);
    	}
	}
    for (k=intervalextension*length;k<(intervalextension+1)*length;k++){
		for (l=0;l<2;l++){
			hazards_speedup[l][0][k] = hazards(k*time1/((double)length),length,1,l);
			hazards_speedup[l][1][k] = hazards(k*time2/((double)length),length,2,l);
			birth_speedup[l][k] = birth_rate(k*time1/((double)length),l+1);
			death_speedup[l][k] = death_rate(k*time1/((double)length),l+1);
			tauND_speedup[l][k] = tauND_rate(k*time1/((double)length),l+1);
			tauD_speedup[l][k] = tauD_rate(k*time1/((double)length),l+1);
    	}
	}
    // the following values are used to determine convergence
    double prev_value_0ND = 0.0;
    double new_value_0ND = 0.0;
    double prev_value_0D = 0.0;
    double new_value_0D = 0.0;
    // counting iterations
    int outside_iteration = 1;
    int inside_iteration = 0;
    // iterations start here
    while (outside_iteration<max_iterations){
        prev_value_0ND = g0[0][length-1];
        prev_value_0D = g0[1][length-1];
	G0(length,g0,inter_g0,hT,time1,1,hazards_speedup,birth_speedup,death_speedup,tauND_speedup,tauD_speedup);
        G0(length,inter_g0,g0,hT,time1,1,hazards_speedup,birth_speedup,death_speedup,tauND_speedup,tauD_speedup);
        new_value_0ND = g0[0][length-1];
        new_value_0D = g0[1][length-1];
        inside_iteration = 0;
	// iterate until G0 converged
        while ((prev_value_0ND-new_value_0ND>accuracy || new_value_0ND-prev_value_0ND>accuracy || prev_value_0D-new_value_0D>accuracy || new_value_0D-prev_value_0D>accuracy) && inside_iteration<max_iterations){
            prev_value_0ND = g0[0][length-1];
            prev_value_0D = g0[1][length-1];
            G0(length,g0,inter_g0,hT,time1,1,hazards_speedup,birth_speedup,death_speedup,tauND_speedup,tauD_speedup);
            G0(length,inter_g0,g0,hT,time1,1,hazards_speedup,birth_speedup,death_speedup,tauND_speedup,tauD_speedup);
            new_value_0ND = g0[0][length-1];
            new_value_0D = g0[1][length-1];
            inside_iteration += 1;
        }
        // printf("# number of iterations for G0 %d\n",inside_iteration);
        // calculate GT based on converged G0 
	GT(length,g0,gT,hT,time1,1,hazards_speedup,birth_speedup,death_speedup,tauND_speedup,tauD_speedup);
        
	// repeat for H0 and HT
	prev_value_0ND = h0[0][length-1];
        prev_value_0D = h0[1][length-1];
        G0(length,h0,inter_h0,gT,time2,2,hazards_speedup,birth_speedup,death_speedup,tauND_speedup,tauD_speedup);
        G0(length,inter_h0,h0,gT,time2,2,hazards_speedup,birth_speedup,death_speedup,tauND_speedup,tauD_speedup);
        new_value_0ND = h0[0][length-1];

        new_value_0D = h0[1][length-1];
        inside_iteration = 0;
	// iterate until H0 converged
        while ((prev_value_0ND-new_value_0ND>accuracy || new_value_0ND-prev_value_0ND>accuracy || prev_value_0D-new_value_0D>accuracy || new_value_0D-prev_value_0D>accuracy) && inside_iteration<max_iterations){
            prev_value_0ND = h0[0][length-1];
            prev_value_0D = h0[1][length-1];
            G0(length,h0,inter_h0,gT,time2,2,hazards_speedup,birth_speedup,death_speedup,tauND_speedup,tauD_speedup);
            G0(length,inter_h0,h0,gT,time2,2,hazards_speedup,birth_speedup,death_speedup,tauND_speedup,tauD_speedup);
            new_value_0ND = h0[0][length-1];
            new_value_0D = h0[1][length-1];
            inside_iteration += 1;
        }
       // printf("# number of iterations for H0 %d\n",inside_iteration);
       // calculate HT based on converged H0
        GT(length,h0,hT,gT,time2,2,hazards_speedup,birth_speedup,death_speedup,tauND_speedup,tauD_speedup);
        
	// every 10 iterations check for convergence
	if (outside_iteration%10==0){
            for (k=0;k<length;k++){
                if (g0[0][k]-out_g0[0][k]>accuracy || out_g0[0][k]-g0[0][k]>accuracy || h0[0][k]-out_h0[0][k]>accuracy || out_h0[0][k]-h0[0][k]>accuracy || g0[1][k]-out_g0[1][k]>accuracy || out_g0[1][k]-g0[1][k]>accuracy || h0[1][k]-out_h0[1][k]>accuracy || out_h0[1][k]-h0[1][k]>accuracy){
                    break;
                }
            }
            if (k==length){
                    printf("# stopped outside iterations at %d\n",outside_iteration);
                    outside_iteration = max_iterations;        
            }
            for (k=0;k<length;k++){
                out_g0[0][k]=g0[0][k];
                out_h0[0][k]=h0[0][k];
                out_g0[1][k]=g0[1][k];
                out_h0[1][k]=h0[1][k];
            }
        }
        outside_iteration += 1;
    }
    //printf("# number of total iterations for G0-GT-H0-HT %d\n",outside_iteration);
    time_t seconds2 = time(NULL); // stop and print time
    printf("# program took %ld min, %ld sec\n",(seconds2-seconds1)/60,(seconds2-seconds1)-((seconds2-seconds1)/60)*60);
    
    // print converged functions
    printf("# columns: time argument for h, hT[0], time argument for g, gT[0],hT[1][k],gT[1][k],h0[0],g0[0],h0[1][k],g0[1][k]\n");
	for (k=0;k<5;k++){ // for initial age 0, only the second line is relevant (k=1), for higher ages, change the upper bound of counter k up to intervalextension*length
        printf("%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n",k*time2/((double) length),hT[0][k],k*time1/((double) length),gT[0][k],hT[1][k],gT[1][k],h0[0][k],g0[0][k],h0[1][k],g0[1][k]);
    }
    //for (k=length;k<intervalextension*length;k++){
    //    printf("%f\t%f\t%f\t%f\t%f\t%f\n",k*time2/((double) length),hT[0][k],k*time1/((double) length),gT[0][k],hT[1][k],gT[1][k]);
    //}
    //for (k=0;k<length;k++){
    //    printf("%f\t%f\t%f\t%f\t%f\n",k*time1/((double) length),g0[k],gT[k],h0[k],hT[k]);
    //}
    fflush(stdout);

}
