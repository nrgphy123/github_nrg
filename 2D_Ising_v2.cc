/* 
Homework 6 Problem 1: Markov chain Monte Carlo method for 2D Ising model
*/

#include <iostream>
#include <stdlib.h> //To use command abs (Absolute value)
#include <stdio.h>
#include <vector>
//#include <time.h> // To assist in generating a random number between 0 and 1.
#include <math.h>       // To use the function exp
#include <random>

using namespace std;

void Initial(vector<vector<double> > &s, int N);

double Magnetization(vector<vector<double> > &s, int N);

double Mag_avg(vector<double> &Mag);

double Suscep_avg(vector<double> &Mag, int N);

void Metropolis(vector<vector<double> > &s, double x, int N);


double fourth_order_cumulant(vector<double> &Mag, int N);

double random_number();


int main(int argc, char **argv){

	int N = 50; //Number of lattice sites in either x or y-direction
	
	int Nsteps = 100000; // Number of Monte Carlo steps
	
	int equil_step = 1000; // Equilibration step size
	
	
	cout << "#Number of lattice sites = " << N*N << endl;
	
	
	//double Magnetization_total = 0;
	
	//double Magnetization_avg;
	
	//vector<double> Mag; //create a vector to store the Magnetization for steps
	
	
	double Mag_total = 0;
	
	//double Mag_avg;
	
	vector<double> Mag; //create a vector to store the Magnetization for steps
	
	
	double step_Mag;
	
	
	// create a array to represent the Square Lattice
	vector<vector<double> > s(N, vector<double> (N));
	
	
	// Initial spin configuration of the lattice
	Initial(s,N);
	
	//Note x = T/Tc.
	
	cout << "\n#x     Magnetization per site (ensemble average)     Susceptibility per site (ensemble average)     Fourth order cumulant"; 
	
	//loop through temperatures and calculate Average Magnetization for each temeprature
	for(double x=0.000; x <= 2.00; x = x + 0.001){
		
		for(int i=0; i < Nsteps; i++){
		
			//set the sampling steps
			if(fmod((double)(i),(double)(equil_step))==0){
				
				Metropolis(s,x,N); // Calculate the energy and spin flip of each spin
				
				step_Mag = Magnetization(s,N); // Calculate total spin
						
				Mag.push_back(double(step_Mag)); //store total spin in the vector Mag
			}
					
		}
		cout << "\n" << x << "\t\t\t\t" << Mag_avg(Mag)/N/N << "\t\t\t\t" << Suscep_avg(Mag, N) << "\t\t\t\t" << fourth_order_cumulant(Mag, N);
		

		Mag.clear(); //Clear the values of the array Mag to start fresh for the value of T/Tc = x

	}
	cout << endl;	
	
	return 0;
}

// Initial spin configuration
void Initial(vector<vector<double> > &s, int N){
	for(int i=0; i<N; i++){
		for(int j=0; j<N; j++){
			s[i][j] = 1;
			}
		}
	}



// Magnetization
double Magnetization(vector<vector<double> > &s, int N){
	
	double Mag = 0;
	
	for( int i=0; i < N; i++){
		for(int j=0; j<N; j++){
			Mag += s[i][j];
		}
	}
	return abs(Mag);
}

// Average Magnetization
double Mag_avg(vector<double> &Mag){
	
	double Mag_total = 0;
	
	for(int i=0; i<Mag.size(); i++){
		Mag_total += Mag[i];
	}

	//average_Mag = Mag_total/Mag.size();
	
	return  Mag_total/Mag.size();
}

// Average Susceptibility

double Suscep_avg(vector<double> &Mag, int N){
	double Mag_square_total = 0;
	
	double suscep; 
	
	for(int i=0; i<Mag.size(); i++){
		Mag_square_total += Mag[i]*Mag[i];
	}
	suscep = (Mag_square_total/Mag.size() - Mag_avg(Mag)*Mag_avg(Mag))/N/N;
	return suscep;
}


double fourth_order_cumulant(vector<double> &Mag, int N){
	double Mag_square_total = 0;
	double Mag_4thpower_total = 0;
	
	double cumulant;
	
	for(int i=0; i<Mag.size(); i++){
		Mag_square_total += Mag[i]*Mag[i];
		Mag_4thpower_total += Mag[i]*Mag[i]*Mag[i]*Mag[i];
	}
	cumulant = 1.0 - (Mag_4thpower_total/Mag.size())/3.0/(Mag_square_total*Mag_square_total/Mag.size()/Mag.size());
	
	return cumulant;
	
}

// Change in energy and Metropolis update
void Metropolis(vector<vector<double> > &s, double x, int N){
	double delH;
	
/*
	J/kB.T = J/kB.(Tc * T)/Tc = J/(kB * Tc *x) = (1/2/x) ln(1 + sqrt(2)) 
*/
	
	double current_spin;
	
	for( int i=0; i < N; i++){
		for(int j=0; j<N; j++){
			
			// Spin at the lattice site under consideration
			current_spin = s[i][j];
		
			int right, left, top, bottom; //Nearest Neighbours 
			
			//Periodic boundary condtions
			if(i==0){left=N-1;}
			else {left=i-1;}

			if(i==N-1){right=0;}
			else {right=i+1;}

			if(j==0){top=N-1;}
			else {top=j-1;}

			if(j==N-1){bottom=0;}
			else {bottom=j+1;}			
		
			//calculate the total spin of the neighbouring lattice points
			double total_neighbour_spin = (s[left][j] + s[right][j] + s[i][top] + s[i][bottom]);
		
			//Change in energy
			// replace J/kB by Tc*1/2(ln(1+sqrt(2))),then use T/Tc = x so that x goes from 0 to 4
			double delH = 2.0*current_spin*(total_neighbour_spin*0.5*log(1.0+sqrt(2.0))/x);
			
			//Metropolis update
			
			if(delH <=0){
				s[i][j] = -1*current_spin;
				}
			else{
				double W = min(1.0, exp(-delH)); 
				double r=random_number();
				if(r<W){
					s[i][j] = -1*current_spin;
					}
				}	
			}
		}
	}

//Random Number Generator
double random_number(){	
	double r = (double)(rand()) / (double) (RAND_MAX);
	return r;
}




