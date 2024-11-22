#include <iostream>
#include "bhjet_class.h"

int main() {
    try {
        // Create an instance of BhJetClass
        BhJetClass model;

        // Load parameters from a file
        std::cout << "Loading parameters from file...\n";
        model.read_params("test_params.dat");

        // Display loaded parameters
        std::vector<double> params = model.get_parameters();
        std::cout << "Loaded parameters:\n";
        for (size_t i = 0; i < params.size(); ++i) {
            std::cout << "Param[" << i << "]: " << params[i] << '\n';
        }

        // Modify a single parameter
        std::cout << "Modifying parameter 0 (Mbh)...\n";
        model.set_parameter(0, 30.0); // Example: Set Mbh to 20.0

        // Display updated parameters
        params = model.get_parameters();
        std::cout << "Updated parameters:\n";
        for (size_t i = 0; i < params.size(); ++i) {
            std::cout << "Param[" << i << "]: " << params[i] << '\n';
        }

        // Run the model
        std::cout << "Running the BHJet model...\n";
        model.run();

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << '\n';
    }

    return 0;
}