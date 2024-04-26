#include <iostream>
#include <vector>
#include <fstream>

class SIR_model {

private:

    // Initiating variables
    double S;
    double I;
    double R;
    double gamma = 0.1;
    double beta = 0.2;
    double N = 1000;

public:
    // Defining inputs
    SIR_model(double S0, double I0, double R0) : S(S0), I(I0), R(R0) {}

    // Defining derivatives
    double dSdt() {
        return - beta * I * S / N;
    }
    double dIdt() {
        return beta * I * S / N - gamma * I;
    }
    double dRdt() {
        return gamma * I;
    }

    // Implementing the forward Euler method
    std::vector<std::vector<double>> forward_Euler(double dt, double max_t) {
        std::vector<std::vector<double>> result_vec;
        
        result_vec.push_back({0, S, I, R});

        for(double t = 0; t<=max_t; t+=dt) {
            S += dt * dSdt();
            I += dt * dIdt();
            R += dt * dRdt();

            result_vec.push_back({t, S, I, R});
            
        }
        return result_vec;
    }
};

int main() {

    // Initial conditions
    double S0 = 699;
    double R0 = 300;
    double I0 = 1;

    // Performing simulation
    SIR_model my_SIR_model{S0, I0, R0};
    std::vector<std::vector<double>> result = my_SIR_model.forward_Euler(0.1, 200);

    // Creating and filling data file
    std::ofstream outFile("sir_output.txt");

    if (outFile.is_open()) {
        // Writing the header
        outFile << "Time S I R\n";

        // Writing data
        for (const auto& row : result) {
            for (const auto& elem : row) {
                outFile << elem << " ";
            }
            outFile << "\n";
        }

        outFile.close();
    } else {
        std::cerr << "Unable to open file";
        return 1;
    }
    
    return 0;
}
