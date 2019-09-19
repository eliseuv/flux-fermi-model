#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>

// Gini Coefficient
double giniIndex(const std::vector<double>& vec) {

    const size_t N = vec.size();

    // Simple sum
    double sum = 0;
    for (double val : vec) sum += val;

    // Sort vector
    std::vector<double> vec_sorted(vec);
    std::sort(vec_sorted.begin(), vec_sorted.end());

    // Operation on sorted vector
    double sum2 = 0, y = 0, y_old = 0;
    for (double val : vec_sorted) {
        y = y_old + val;
        sum2 += y_old + y;
        y_old = y;
    }

    return (1.0 - (sum2/sum)/N);
}

int main()
{
    // System dimensions
    const int l_x = 128;
    const int l_y = 1;

    const int it_max = 1;
    const int idelta_max = 1;
    const int idt = 1000;
    const double epsilon = 1e-4;
    const double erro = 1e-1;
    const double dif = 0;
    const double deq = 1;
    const double alpha = 100;

    // Number of particles
    int n_a = 1;
    int n_b = 1;
    int n_total = n_a + n_b;

    // Density of particles
    double rho_a = static_cast<double>(n_a)/(l_x*l_y);
    double rho_b = static_cast<double>(n_b)/(l_x*l_y);
    double rho_total = static_cast<double>(n_total)/(l_x*l_y);

    //
    const int n_run = 1000/n_total + 200;
    int rand_seed = 793;
    std::mt19937 engine(rand_seed);
    std::uniform_real_distribution<double> uniform_dist(0.0, 1.0);

    // 
    std::vector<double> phi_l_x(l_y), phi_l_y(l_x);

    // 
    std::vector<double> vec_a(l_x*l_y), vec_b(l_x*l_y);

    //
    std::vector<int> iv_x(l_x+2), iv_y(l_y+2);

    // iv_x = (l_x, 1, 2, ... l_x, 1)
    for (int i = 1; i < l_x+1; i++) iv_x[i] = i;
    iv_x[0] = iv_x[l_x];
    iv_x[l_x+1] = iv_x[1];

    // iv_y = (l_y, 1, 2, ... l_y, 1)
    for (int i = 1; i < l_y+1; i++) iv_y[i] = i;
    iv_y[0] = iv_y[l_y];
    iv_y[l_y+1] = iv_y[1];

    // Movement of the particles
    std::vector<int> i_x(5, 0), i_y(5, 0);

    // Movement probabilities
    std::vector<double> p_x(5, 0);

    // i_x = (+1, 0, 0, -1, 0)
    i_x[0] = +1; i_x[3] = -1;

    // i_y = (0, +1, -1, 0, 0)
    i_y[1] = +1; i_y[2] = -1;

    // External loop (changes number of particles)
    while (rho_total < 4.0) {

        // Particle positions
        std::vector<int> in_x(n_total), in_y(n_total);

        // Measured parameters
        double  avg_mobility = 0, var_mobility = 0,                             // Mobility
                avg_phi_x = 0, avg_phi_y = 0, var_phi_x = 0, var_phi_y = 0,     // Order parameter
                avg_j_a = 0, avg_j_b = 0, var_j_a = 0, var_j_b = 0,             // Current density
                avg_gini_a = 0, avg_gini_b = 0, var_gini_a = 0, var_gini_b = 0, // Gini index
                p_lane = 0, p_clog = 0, p_else = 0;                             // Probabilities assign to each scenario

        // Sample loop
        for (int i_run = 0; i_run < n_run; i_run++) {

            // New seed for PRNG
            rand_seed += 35;
            engine.seed(rand_seed);

            // Lattice matrices
            std::vector<std::vector<int>> ic_a(l_x+1, std::vector<int>(l_y+1, 0)),
                                          ic_b(l_x+1, std::vector<int>(l_y+1, 0));
        
            // Random initial distribution of particles
            for (int j = 1; j <= n_total; j++) {
                // Generate coordinates
                int ixl = static_cast<int>(l_x*uniform_dist(engine)) + 1;
                int iyl = static_cast<int>(l_y*uniform_dist(engine)) + 1;
                // Fill domain with A and the with B particles
                if (j <= n_a) ic_a[ixl][iyl]++;
                else ic_b[ixl][iyl]++;
                // Store coordinates in position vectors
                in_x[j] = ixl;
                in_y[j] = iyl;
            }

            // Aux variables
            int mc_iter = 0;                            // Monte Carlo iteration
            int i_count = 1, i_count2 = 0;              // Stopping criterium interval index and time step
            int ism=0, isx=0, isy=0, isga=0, islo=0;    // Aux vars
            double  a=0, b=0, c=0, d=0;                 // Aux vars
            double  d_mom_mobility = 0,                 // Time average of mobility over interval i_dt
                    d_mom_phi_x = 0, d_mom_phi_y = 0,   // Time average of order parameter over interval i_dt
                    d_mom_j_a = 0, d_mom_j_b = 0,       // Time average of current density over interval i_dt
                    d_mom_gini_a = 0, d_mom_gini_b = 0; // Time average of Gini index over interval i_dt

            // Time loop
            while (ism == 0) { // Until mobility stationary state

                mc_iter++;
                i_count2++;
                rand_seed += 24;
                engine.seed(rand_seed);

                int in_a = 0, in_b = 0;
                double  phi_x = 0, phi_y = 0,                       
                        gini_a = 0, gini_b = 0;

                // Sequential Dynamics: n_total particles are selected
                for (int i_particle = 0; i_particle < n_total; i_particle++) {

                    int i = static_cast<int>(n_total*uniform_dist(engine));

                    if (i <= n_a) { // A particles first
                        
                        // sigma^A(j+1,k) da exponencial de Fermi
                        double aux1 = dif*ic_b[iv_x[in_x[i]+i_x[0]]][iv_y[in_y[i]+i_y[0]]] + deq*ic_a[iv_x[in_x[i]+i_x[0]]][iv_y[in_y[i]+i_y[0]]];
                        
                        // sigma^A(j,k+1) da exponencial de Fermi
                        double aux2 = dif*ic_b[iv_x[in_x[i]+i_x[1]]][iv_y[in_y[i]+i_y[1]]] + deq*ic_a[iv_x[in_x[i]+i_x[1]]][iv_y[in_y[i]+i_y[1]]];
                        
                        // sigma^A(j,k-1) da exponencial de Fermi
                        double aux3 = dif*ic_b[iv_x[in_x[i]+i_x[2]]][iv_y[in_y[i]+i_y[2]]] + deq*ic_a[iv_x[in_x[i]+i_x[2]]][iv_y[in_y[i]+i_y[2]]];
                    
                        // Calculate movement probabilities
                        p_x[0] = 1.0/(1.0 + std::exp(alpha*(aux1-idelta_max)));
                        p_x[1] = 1.0/(1.0 + std::exp(alpha*(aux2-idelta_max)))*(1.0 - 1.0/(1.0 + std::exp(alpha*(aux1-idelta_max))));
                        p_x[2] = 1.0/(1.0 + std::exp(alpha*(aux3-idelta_max)))*(1.0 - 1.0/(1.0 + std::exp(alpha*(aux1-idelta_max))));
                        p_x[3] = 0;
                        p_x[4] = 0;

                        // Normalization of probabilities
                        double anorm = 0;
                        for (double prob : p_x) anorm += prob;
                        if (anorm < 1) p_x[4] = 1 - anorm;
                        else if (anorm > 1) for (double& prob : p_x) prob /= anorm;

                        // Random particle displacement
                        double z = uniform_dist(engine);
                        int k = 0;
                        double acum = 0;
                        while (acum < z) {
                            acum += p_x[k];
                            k++;
                        }

                        if (k == 1) // x displacement
                            in_a++;
                        
                        // Update lattice
                        ic_a[in_x[i]][in_y[i]]--;
                        in_x[i] = iv_x[in_x[i]+i_x[k]];
                        in_y[i] = iv_y[in_y[i]+i_y[k]];                     
                        ic_a[in_x[i]][in_y[i]]++;

                    }       // A particles first
                    else {  // and thn B particles

                        // sigma^B(j-1,k) da exponencial de Fermi
                        double aux1 = deq*ic_b[iv_x[in_x[i]-i_x[0]]][iv_y[in_y[i]+i_y[0]]] + dif*ic_a[iv_x[in_x[i]-i_x[0]]][iv_y[in_y[i]+i_y[0]]];
                        
                        // sigma^B(j,k+1) da exponencial de Fermi
                        double aux2 = deq*ic_b[iv_x[in_x[i]-i_x[1]]][iv_y[in_y[i]+i_y[1]]] + dif*ic_a[iv_x[in_x[i]-i_x[1]]][iv_y[in_y[i]+i_y[1]]];
                        
                        // sigma^B(j,k-1) da exponencial de Fermi
                        double aux3 = deq*ic_b[iv_x[in_x[i]-i_x[2]]][iv_y[in_y[i]+i_y[2]]] + dif*ic_a[iv_x[in_x[i]-i_x[2]]][iv_y[in_y[i]+i_y[2]]];
 
                        // Calculate movement probabilities
                        p_x[0] = 1.0/(1.0 + std::exp(alpha*(aux1-idelta_max)));
                        p_x[1] = 1.0/(1.0 + std::exp(alpha*(aux2-idelta_max)))*(1.0 - 1.0/(1.0 + std::exp(alpha*(aux1-idelta_max))));
                        p_x[2] = 1.0/(1.0 + std::exp(alpha*(aux3-idelta_max)))*(1.0 - 1.0/(1.0 + std::exp(alpha*(aux1-idelta_max))));
                        p_x[3] = 0;
                        p_x[4] = 0;

                        // Normalization of probabilities
                        double anorm = 0;
                        for (double prob : p_x) anorm += prob;
                        if (anorm < 1) p_x[4] = 1 - anorm;
                        else if (anorm > 1) for (double& prob : p_x) prob /= anorm;

                        // Random particle displacement
                        double z = uniform_dist(engine);
                        int k = 0;
                        double acum = 0;
                        while (acum < z) {
                            acum += p_x[k];
                            k++;
                        }

                        if (k == 1) // x displacement
                            in_b++;
                        
                        // Update lattice
                        ic_b[in_x[i]][in_y[i]]--;
                        in_x[i] = iv_x[in_x[i]-i_x[k]];
                        in_y[i] = iv_y[in_y[i]+i_y[k]];                     
                        ic_b[in_x[i]][in_y[i]]++;

                    } // and then B particles
                } // Sequential dynamics

                // Update vectors and matrices
                for (int i = 0; i < l_x; i++) phi_l_y[i] = 0;
                for (int i = 0; i < l_y; i++) phi_l_x[i] = 0;
                
                for (int l = 0; l < l_y; l++) {
                    for (int k = 0; k < l_x; k++) {
                        phi_l_x[l] += ic_a[k][l] - ic_b[k][l];
                        phi_l_y[k] += ic_a[k][l] - ic_b[k][l];
                        vec_a[k+l_x*(l-1)] = ic_a[k][l];
                        vec_b[k+l_x*(l-1)] = ic_b[k][l];
                    }
                    phi_x += std::abs(phi_l_x[l]);
                }

                for (int l = 0; l < l_y; l++)
                    phi_y += std::abs(phi_l_y[l]);

                gini_a = giniIndex(vec_a);
                gini_b = giniIndex(vec_b);

                // Intantaneous order parameter, mobility and current density
                phi_x /= n_total;
                phi_y /= n_total;
                double d_mobility = static_cast<double>(in_a+in_b)/n_total;
                double d_j_a = static_cast<double>(in_a)/(l_x*l_y);
                double d_j_b = static_cast<double>(in_b)/(l_x*l_y);

                // Integrate over time
                d_mom_mobility += d_mobility;
                d_mom_phi_x += phi_x;
                d_mom_phi_y += phi_y;
                d_mom_j_a += d_j_a;
                d_mom_j_b += d_j_b;
                d_mom_gini_a += gini_a;
                d_mom_gini_b += gini_b;

                // Linear fit
                a += i_count2*d_mobility;
                b += i_count2*phi_x;
                c += i_count2*phi_y;
                d += i_count2*gini_a;

                double di = static_cast<double>(idt*(idt+1))/2;
                double di2 = static_cast<double>((2*idt+1)*di)/3;

                // Check if the system has reached stationary state
                if (mc_iter == i_count*idt) {

                }

            } // Time loop

        } // Sample loop

    } // External loop

    return 0;
}