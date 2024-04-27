#include <iostream>
#include <fstream>
using namespace std;

// Lorenz system parameters
const double sigma = 10.0; // Prandtl number
const double rho = 28.0; // Rayleigh number
const double beta = 8.0 / 3.0; // Aspect ratio constant

// Function to calculate the derivatives of the Lorenz system
void lorenzSystem(double x, double y, double z, double& dxdt, double& dydt, double& dzdt) {
    dxdt = sigma * (y - x);
    dydt = x * (rho - z) - y;
    dzdt = x * y - beta * z;
}

// Fourth-order Runge-Kutta method to integrate the Lorenz system
void rungeKutta(double& x, double& y, double& z, double dt) {
    // Variables to hold the intermediate steps of the Runge-Kutta method
    double k1x, k1y, k1z;
    double k2x, k2y, k2z;
    double k3x, k3y, k3z;
    double k4x, k4y, k4z;

    // Calculate the first set of intermediate values (k1)
    lorenzSystem(x, y, z, k1x, k1y, k1z);
    k1x *= dt;
    k1y *= dt;
    k1z *= dt;

    // Calculate the second set of intermediate values (k2)
    lorenzSystem(x + k1x / 2, y + k1y / 2, z + k1z / 2, k2x, k2y, k2z);
    k2x *= dt;
    k2y *= dt;
    k2z *= dt;

    // Calculate the third set of intermediate values (k3)
    lorenzSystem(x + k2x / 2, y + k2y / 2, z + k2z / 2, k3x, k3y, k3z);
    k3x *= dt;
    k3y *= dt;
    k3z *= dt;

    // Calculate the fourth set of intermediate values (k4)
    lorenzSystem(x + k3x, y + k3y, z + k3z, k4x, k4y, k4z);
    k4x *= dt;
    k4y *= dt;
    k4z *= dt;

    // Update x, y, and z using the Runge-Kutta method
    x += (k1x + 2 * k2x + 2 * k3x + k4x) / 6;
    y += (k1y + 2 * k2y + 2 * k3y + k4y) / 6;
    z += (k1z + 2 * k2z + 2 * k3z + k4z) / 6;
}

// Function to simulate the Lorenz system and save data to a file
void simulateLorenz(double x0, double y0, double z0, double dt, double duration) {
    // Initial conditions
    double x = x0;
    double y = y0;
    double z = z0;

    // Open a file to save the data
    ofstream dataFile("lorenz_data.txt");
    dataFile << "Time X Y Z\n"; // Write headers

    // Simulate the Lorenz system over the specified duration
    for (double t = 0.0; t <= duration; t += dt) {
        // Save the current time, x, y, and z values to the file
        dataFile << t << " " << x << " " << y << " " << z << "\n";

        // Integrate the Lorenz system using the Runge-Kutta method
        rungeKutta(x, y, z, dt);
    }

    // Close the file
    dataFile.close();
}

// Main function
int main() {
    // Define initial conditions
    double x0 = 1.0;
    double y0 = 1.0;
    double z0 = 1.0;

    // Define time step and duration for the simulation
    double dt = 0.01; // Time step (s)
    double duration = 50.0; // Duration of the simulation (s)

    // Simulate the Lorenz system and save data to a file
    simulateLorenz(x0, y0, z0, dt, duration);

    return 0;
}
