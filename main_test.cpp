#include "bhjet_class.hh"
#include <iostream>

int main() {
    try {
        BhJetClass jet;

        // Read parameters from a file
        jet.load_params("params.dat");

        // Set basic run parameters
        jet.set_basic_run_params();

        // Run the model
        jet.run();

        // Access the output if needed
        // const JetOutput& output = jet.get_output();

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
    }

    return 0;
}