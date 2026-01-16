/**
 * Numerical solution of the 1D transient heat conduction equation
 * in a solid wall with internal heat generation.
 *
 * Domain:
 *   z ∈ [0, L]
 *
 * Governing equation:
 *   ρ c_p ∂T/∂t = k ∂²T/∂z² + Q
 *
 * Boundary conditions:
 *   - z = 0 : zero heat flux (adiabatic)
 *             ∂T/∂z = 0
 *
 *   - z = L : prescribed temperature
 *             T = 300 K
 *
 * Initial condition:
 *   - Uniform temperature field
 *             T(z, 0) = 300 K
 *
 * Source term:
 *   - Uniform, constant volumetric heat generation
 *             Q = const [W/m³]
 */

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <array>
#include <omp.h>

#include <tdma.h>

int main() {

    // ===================================================================
    //                    NUMERICAL / PHYSICAL SETUP
    // ===================================================================

	constexpr int N = 100;                          // Number of wall cells
	constexpr double L = 1.0;                       // Wall length [m]
	constexpr double dz = L / N;                    // Cell size [m]
	constexpr double dt = 1;                     // Time step [s]
	constexpr int time_iter = 1000;                 // Number of time iterations

	constexpr double k = 20.0;                      // Steel thermal conductivity [W/mK]
    constexpr double rho = 7850.0;                  // Steel density [kg/m3]
	constexpr double cp = 500.0;                    // Steel specific heat [J/kgK]
    constexpr double T_amb = 300.0;                 // Ambient temperature [K]

	std::vector<double> T_w_bulk(N, 300.0);         // Initial wall temperature [K]
	std::vector<double> T_w_bulk_old;               // Old temperature [K]
	std::vector<double> Q(N, 0.0);                  // Heat pipe volumetric source term [W/m3]

	// Tridiagonal matrix coefficients
    std::vector<double> aTW(N, 0.0);
    std::vector<double> bTW(N, 0.0);
    std::vector<double> cTW(N, 0.0);
    std::vector<double> dTW(N, 0.0);

	// Output file
    std::ofstream file("wall_numerical.dat");

    double start = omp_get_wtime();

    // Time loop
    for (int n = 0; n < time_iter; ++n) {

        const double time = (n + 1) * dt;

        // Store previous temperature
        T_w_bulk_old = T_w_bulk;

        // Assembly loop
        for (int i = 1; i < N - 1; ++i) {

            aTW[i] = -k / (rho * cp * dz * dz);
            bTW[i] = 1.0 / dt + (k + k) / (rho * cp * dz * dz);
            cTW[i] = -k / (rho * cp * dz * dz);

            dTW[i] =
                T_w_bulk_old[i] / dt
                + Q[i] / (rho * cp);
        }

        // Wall temperature BCs
        aTW[0] = 0.0;
        bTW[0] = 1.0;
        cTW[0] = 0.0;
        dTW[0] = 350.0;

        aTW[N - 1] = 0.0;
        bTW[N - 1] = 1.0;
        cTW[N - 1] = 0.0;
        dTW[N - 1] = T_amb;

        // Solve
        T_w_bulk = solveTridiagonal(aTW, bTW, cTW, dTW);

        // Output
        for (int i = 0; i < N; ++i)
            file << T_w_bulk[i] << ", ";

        file << "\n";
        file.flush();
    }

    double end = omp_get_wtime();
    std::cout << "Execution time: " << end - start;

    file.close();

    return 0;
}
