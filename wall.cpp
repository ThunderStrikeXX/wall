#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <array>

// =======================================================================
//                        [TDMA ALGORITHM]
// =======================================================================

std::vector<double> solveTridiagonal(const std::vector<double>& a,
    const std::vector<double>& b,
    const std::vector<double>& c,
    const std::vector<double>& d) {

    const int n = static_cast<int>(b.size());

    std::vector<double> c_star(n, 0.0);
    std::vector<double> d_star(n, 0.0);
    std::vector<double> x(n, 0.0);

    c_star[0] = c[0] / b[0];
    d_star[0] = d[0] / b[0];

    for (int i = 1; i < n; ++i) {
        const double m = b[i] - a[i] * c_star[i - 1];
        c_star[i] = c[i] / m;
        d_star[i] = (d[i] - a[i] * d_star[i - 1]) / m;
    }

    x[n - 1] = d_star[n - 1];

    for (int i = n - 2; i >= 0; --i)
        x[i] = d_star[i] - c_star[i] * x[i + 1];

    return x;
}

int main() {

    // ===================================================================
    //                    NUMERICAL / PHYSICAL SETUP
    // ===================================================================

	constexpr int N = 100;                          // Number of wall cells
	constexpr double L = 1.0;                       // Wall length [m]
	constexpr double dz = L / N;                    // Cell size [m]
	constexpr double dt = 1e-1;                     // Time step [s]
	constexpr int time_iter = 1000;                 // Number of time iterations

	constexpr double k = 20.0;                      // Steel thermal conductivity [W/mK]
    constexpr double rho = 7850.0;                  // Steel density [kg/m3]
	constexpr double cp = 500.0;                    // Steel specific heat [J/kgK]
    constexpr double T_amb = 300.0;                 // Ambient temperature [K]

	std::vector<double> T_w_bulk(N, 300.0);         // Initial wall temperature [K]
	std::vector<double> T_w_bulk_old;               // Old temperature [K]
	std::vector<double> Q(N, 1e6);                // Heat pipe volumetric source term [W/m3]

	// Tridiagonal matrix coefficients
    std::vector<double> aTW(N, 0.0);
    std::vector<double> bTW(N, 0.0);
    std::vector<double> cTW(N, 0.0);
    std::vector<double> dTW(N, 0.0);

	// Output file
    std::ofstream file("T_wall.dat");

    // Time loop
    for (int n = 0; n < time_iter; ++n) {

        const double time = (n + 1) * dt;

        // Store previous temperature
        T_w_bulk_old = T_w_bulk;

        // ===================================================================
        //                      ASSEMBLY LOOP
        // ===================================================================

        for (int i = 1; i < N - 1; ++i) {

            aTW[i] = -k / (rho * cp * dz * dz);
            bTW[i] = 1.0 / dt + (k + k) / (rho * cp * dz * dz);
            cTW[i] = -k / (rho * cp * dz * dz);

            dTW[i] =
                T_w_bulk_old[i] / dt
                + Q[i] / (rho * cp);
        }

        // Wall temperature BCs (zero flux at first and last node)
        aTW[0] = 0.0;
        bTW[0] = 1.0;
        cTW[0] = -1.0;
        dTW[0] = 0.0;

        aTW[N - 1] = 0.0;
        bTW[N - 1] = 1.0;
        cTW[N - 1] = 0.0;
        dTW[N - 1] = T_amb;

        // ===================================================================
        //                          SOLVE
        // ===================================================================

        T_w_bulk = solveTridiagonal(aTW, bTW, cTW, dTW);

        // ===================================================================
        //                          OUTPUT
        // ===================================================================


        for (int i = 0; i < N; ++i)
            file << T_w_bulk[i] << " ";

        file << "\n";
        file.flush();
    }

    file.close();

    return 0;
}
