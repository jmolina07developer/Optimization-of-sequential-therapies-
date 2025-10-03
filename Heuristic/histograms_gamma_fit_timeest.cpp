#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <numeric>
#include <cmath>

// Parameters
double mu1 = 1e-3;
//double mu1 = 1e-4;
double mu2 = 1e-2;
//double mu2 = 1e-3;
double kA = 0.5;
double kB = 0.5;
double CS = 0.05;
double CR = 1;
double gammaV = 0.01;

// crear la funcion de linspace
template <typename T>
std::vector<T> linspace(double start, double end, double num)
{
    std::vector<T> linspaced;

    if (0 != num)
    {
        if (1 == num) 
        {
            linspaced.push_back(static_cast<T>(start));
        }
        else
        {
            double delta = (end - start) / (num - 1);

            for (auto i = 0; i < (num - 1); ++i)
            {
                linspaced.push_back(static_cast<T>(start + delta * i));
            }
            // ensure that start and end are exactly the same as the input
            linspaced.push_back(static_cast<T>(end));
        }
    }
    return linspaced;
}


// Function to calculate N (total population)
double calculate_N(const std::vector<double> &x) {
    return x[0] + x[1] + x[2] + x[3];
}

// REPLICATION RATES SEASON A
std::vector<double> propensities_A(const std::vector<double> &x) {
  double N = calculate_N(x);
  return {kA*x[0], x[1], CS*kA*x[2], CR*x[3]};
}
std::vector<double> propensities_B(const std::vector<double> &x) {
  double N = calculate_N(x);
  return {kB*x[0], CS*kB*x[1], x[2], CR*x[3]};
}
std::vector<std::vector<double> > mutations={{1-2*mu1,mu1,mu1,0},
					     {mu2,1-mu1-mu2,0,mu1},
					     {mu2,0,1-mu1-mu2,mu1},
					     {0,mu2,mu2,1-2*mu2}};

//SAMPLING MUTACIONES:
std::vector<int> multinomialSampling(const std::vector<double>& probabilities, int num_samples, std::mt19937 &rng) {
    // Check if probabilities sum to 1
    double sum = std::accumulate(probabilities.begin(), probabilities.end(), 0.0);
    if (std::abs(sum - 1.0) > 1e-6) {
        throw std::invalid_argument("Probabilities must sum to 1.");
    }
    // Create the discrete distribution
    std::discrete_distribution<> dist(probabilities.begin(), probabilities.end());

    // Perform multinomial sampling
    std::vector<int> counts(probabilities.size(), 0);
    for (int i = 0; i < num_samples; ++i) {
        int sampled_index = dist(rng);
        counts[sampled_index]++;
    }
    return counts;
}

// Tau-leaping step
void tau_leap_step(std::vector<double> &x, double tau, std::vector<double> (*propensity_func)(const std::vector<double> &), const std::vector<std::vector<double> > &mutations, std::mt19937 &rng) {
    std::vector<double> propensities = propensity_func(x);
    std::poisson_distribution<int> poisson;
    //BIRTHS
    for (size_t i = 0; i < propensities.size(); ++i){
      poisson.param(std::poisson_distribution<int>::param_type(propensities[i] * tau));
      int births=poisson(rng);//births
      auto mut=mutations[i];
      auto smuts=multinomialSampling(mut,births,rng);
      for (size_t j = 0; j < x.size(); ++j) {
	x[j] += smuts[j];
      }
    }
    //DEATHS
    auto N=calculate_N(x);
    for (size_t i = 0; i < propensities.size(); ++i){
      poisson.param(std::poisson_distribution<int>::param_type(gammaV*N*x[i]*tau));
      int deaths=poisson(rng);//births
      x[i] -= deaths;
      if (x[i]<0) x[i]=0;
    }      
}

// Determine the current season
std::vector<double> (*current_propensity_function(double t, double season_length, double t_end))(const std::vector<double> &) {
    double period_A = season_length;
    double period_B = season_length;
    double total_period = period_A + period_B;
    double time_in_cycle = fmod(t, total_period);

    if (t < t_end) {
        if (time_in_cycle < period_A) {
        return propensities_A; // Season A
        } else {
        return propensities_B; // Season B
        }
    }
    else {

        if (fmod(t_end + 0.1, total_period) < period_A) {
        return propensities_A; // Season A
        } else {
        return propensities_B; // Season B
        }

    }

    
}

int main() {
    // Simulation parameters
    double t_start = 0.0;
    double t_end = 1050.0;
    double t_end_treatment = 1000;
    double tau = 0.1; // Time step for tau-leaping
    int num_trajectories = 10000; // Number of trajectories per season_length
    double K = 1 / gammaV;

    //auto season_lengths = linspace<double>(1, 100, 100); // Tiempos de tratamiento
    auto season_lengths = linspace<double>(1000, 1000, 1); // Tiempos de tratamiento

    // Random number generator
    std::random_device rd;
    std::mt19937 rng(rd());
    std::normal_distribution<double> gaussian(0.0, 1.0);

    // Output extinction results
    std::ofstream out("time_ext_x2.csv");
    out << "Trajectory,FinalTime" << std::endl;

    
    // Output trajectories
    //std::ofstream traj_out("trajectories_n10000.csv");
    //traj_out << "SeasonLength,Trajectory,Time,x0,x1,x2,x3" << std::endl;

    // Output complete trajectories
    //std::ofstream traj_plot("trajectories_plot.csv");
    //traj_plot << "SeasonLength,Trajectory,Time,x0,x1,x2,x3" << std::endl;
    

    
        for (double season_length : season_lengths) { // Bucle para recorrer los valores de season_length
            int extinction_count = 0;
            double avgN = 0.0;

            for (int traj = 0; traj < num_trajectories; ++traj) {
                // Initial state
                std::vector<double> x = {0.0, 0.0, 100.0, 0.0};
                double t = t_start;
                int extinct = 0;
                int Ntot = 0;

                // Run trajectory
                while (t < t_end) {
                    // Determine the current propensity function based on the season
                    auto propensity_func = current_propensity_function(t, season_length, t_end_treatment); // Incluye CS
            
                    tau_leap_step(x, tau, propensity_func, mutations, rng);
                    t += tau;

                    auto N = calculate_N(x);
                    Ntot += N;

                
                    // Check for extinction
                    if (calculate_N(x) <= 0.05 * K) {
                        extinct = 1;
                        std::cout << "Traj: " << traj << std::endl;
                        out << traj << "," << t << std::endl;
                        break;
                    }
                
                }
                
                // ADD AVERAGES
                extinction_count += extinct;
                avgN += Ntot;
            }
            avgN /= t_end;
            avgN *= tau;
            avgN /= num_trajectories;

            
            // Calculate extinction rate
            
            double extinction_rate = static_cast<double>(extinction_count) / num_trajectories;
            
        }

    out.close();
    //traj_out.close();
    //traj_plot.close();
    //return 0;
}
