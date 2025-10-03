#include <iostream>
#include <vector>
#include <fstream>
#include <random>
#include <numeric>
#include <cmath>

// Parámetros
double mu1 = 1e-3;
double mu2 = 1e-2;
double kA = 0.5;
double kB = 0.5;
double CS = 0.05;
double CR = 1;
double gammaV = 0.01;

template <typename T>
std::vector<T> linspace(double start, double end, double num) {
    std::vector<T> linspaced;
    if (num != 0) {
        if (num == 1) {
            linspaced.push_back(static_cast<T>(start));
        } else {
            double delta = (end - start) / (num - 1);
            for (auto i = 0; i < (num - 1); ++i) {
                linspaced.push_back(static_cast<T>(start + delta * i));
            }
            linspaced.push_back(static_cast<T>(end));
        }
    }
    return linspaced;
}

double calculate_N(const std::vector<double> &x) {
    return x[0] + x[1] + x[2] + x[3];
}

std::vector<double> propensities_A(const std::vector<double> &x) {
    return {kA * x[0], x[1], CS * kA * x[2], CR * x[3]};
}

std::vector<double> propensities_B(const std::vector<double> &x) {
    return {kB * x[0], CS * kB * x[1], x[2], CR * x[3]};
}

std::vector<std::vector<double>> mutations = {{1 - 2 * mu1, mu1, mu1, 0},
                                               {mu2, 1 - mu1 - mu2, 0, mu1},
                                               {mu2, 0, 1 - mu1 - mu2, mu1},
                                               {0, mu2, mu2, 1 - 2 * mu2}};

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

void tau_leap_step(std::vector<double> &x, double tau, std::vector<double> (*propensity_func)(const std::vector<double> &), std::mt19937 &rng) {
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

void gillespie_step(std::vector<double> &x, std::vector<double> (*propensity_func)(const std::vector<double> &), std::mt19937 &rng, double &t) {
    std::vector<double> propensities = propensity_func(x);
    double a0 = std::accumulate(propensities.begin(), propensities.end(), 0.0);
    double N = calculate_N(x);
    double death_rate = gammaV * N;
    a0 += death_rate * N;
    
    if (a0 == 0) return;
    
    std::exponential_distribution<double> exp_dist(a0);
    t += exp_dist(rng);
    
    propensities.push_back(death_rate * N);
    std::discrete_distribution<int> reaction_dist(propensities.begin(), propensities.end());
    int reaction = reaction_dist(rng);
    
    if (reaction < x.size()) {
        x[reaction] += 1;
        std::vector<double> mutation_probs = mutations[reaction];
        std::discrete_distribution<int> mutation_dist(mutation_probs.begin(), mutation_probs.end());
        int new_state = mutation_dist(rng);
        x[reaction] -= 1;
        x[new_state] += 1;
    } else {
        std::discrete_distribution<int> death_dist(x.begin(), x.end());
        int victim = death_dist(rng);
        if (x[victim] > 0) x[victim] -= 1;
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
    double t_start = 0.0, tau = 0.1, K = 1 / gammaV;

    //double t_end = 150.0;
    //double t_end_treatment = 100;
    double t_end = 150.0;
    double t_end_treatment = 100;
    std::random_device rd;
    std::mt19937 rng(rd());

    std::ofstream final_time_out("Gillespie_extrate_t100_boundarytheshold5_gillespie15.csv");
    final_time_out << "SeasonLength,Trajectory,FinalTime" << std::endl;
    
    //std::ofstream traj_out("trajectories_e5_tf_1000.csv");
    //traj_out << "CS,SeasonLength,Trajectory,Time,x0,x1,x2,x3" << std::endl;

    auto season_lengths = linspace<double>(5, 100, 5);
    //auto season_lengths = linspace<double>(50, 500, 5);
    for (double season_length : season_lengths) {
        for (int traj = 0; traj < 10000; ++traj) {
            std::vector<double> x = {50, 0.0, 0.0, 0.0};
            double t = t_start;

            while (t < t_end) {

                auto propensity_func = current_propensity_function(t, season_length, t_end_treatment); // Incluye CS

                if (calculate_N(x) > 15) {
                    tau_leap_step(x, tau, propensity_func, rng);
                    t += tau;
                } else {
                    gillespie_step(x, propensity_func, rng, t);
                }
                if (calculate_N(x) <= 0.05*K) break;
                
            }

            final_time_out << season_length << "," << traj << "," << t << std::endl;
        }
        std::cout << "Season Length: " << season_length << std::endl;
    }
    final_time_out.close();
    //traj_out.close();
    return 0;
}
