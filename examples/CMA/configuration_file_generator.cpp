#include <iostream>
#include <fstream>
#include <string>

struct Parameters
{
    std::string EXPERIMENT_NAME{};
    std::string RESULTS_DIR{};
    double DELTA{};
    double DELTA_DELTA_T{};
    std::size_t CELL_NX{};
    std::size_t CELL_NY{};
    double CELL_DX{};
    double CELL_DY{};
    bool GHOSTED{};
    std::size_t NUM_MT{};
    double MT_MIN_SEGMENT_INIT{};
    double MT_MAX_SEGMENT_INIT{};
    std::size_t NUM_STEPS{};
    double LENGTH_DIV_FACTOR{};
    double DIV_LENGTH{};
    double DIV_LENGTH_RETRACT{};
    double V_PLUS{};
    double V_MINUS{};
    double SIGMOID_K{};
    double TOTAL_TIME{};
    double MAXIMAL_REACTION_RADIUS{};
    double DELTA_T_MIN{};
    double RHO_TEST_RATE{}; //a tunable test parameter for MT dynamics
};

void serialize(Parameters& settings, std::string filename)
{
    std::ofstream outfile(filename);

    // Check if the file was opened successfully
    if (!outfile.is_open()) {
        std::cerr << "Error opening file for writing." << std::endl;
    }

    outfile << "{\n";
    outfile << "\t\"EXPERIMENT\": \"" <<  settings.EXPERIMENT_NAME << "\",\n";
    outfile << "\t\"RESULTS_DIR\": \"" <<  settings.RESULTS_DIR << "\",\n";
    outfile << "\t\"TOTAL_TIME\": " << settings.TOTAL_TIME << ",\n";
    outfile << "\t\"CELL_NX\": " << settings.CELL_NX << ",\n";
    outfile << "\t\"CELL_NY\": " << settings.CELL_NY << ",\n";
    outfile << "\t\"CELL_DX\": " << settings.CELL_DX << ",\n";
    outfile << "\t\"CELL_DY\": " << settings.CELL_DY << ",\n";
    settings.GHOSTED ? outfile << "\t\"GHOSTED\": true,\n" : outfile << "\t\"GHOSTED\": false,\n";
    outfile << "\t\"NUM_MT\": " << settings.NUM_MT << ",\n";
    outfile << "\t\"MT_MIN_SEGMENT_INIT\": " << settings.MT_MIN_SEGMENT_INIT << ",\n";
    outfile << "\t\"MT_MAX_SEGMENT_INIT\": " << settings.MT_MAX_SEGMENT_INIT << ",\n";
    outfile << "\t\"LENGTH_DIV_FACTOR\": " << settings.LENGTH_DIV_FACTOR << ",\n";
    outfile << "\t\"DIV_LENGTH\": " << settings.DIV_LENGTH << ",\n";
    outfile << "\t\"DIV_LENGTH_RETRACT\": " << settings.DIV_LENGTH_RETRACT << ",\n";
    outfile << "\t\"V_PLUS\": " << settings.V_PLUS << ",\n";
    outfile << "\t\"V_MINUS\": " << settings.V_MINUS << ",\n";
    outfile << "\t\"MAXIMAL_REACTION_RADIUS\": " << settings.MAXIMAL_REACTION_RADIUS << ",\n";
    outfile << "\t\"RHO_TEST_RATE\": " << settings.RHO_TEST_RATE << ",\n";
    outfile << "\t\"DELTA\": " << settings.DELTA << ",\n";
    outfile << "\t\"SIGMOID_K\": " << settings.SIGMOID_K << "\n";
    outfile << "}\n";

    outfile.close();
}

//function to set the default parameters that we will modifiy for experiments
void set_default(Parameters& settings)
{
    settings.EXPERIMENT_NAME = "my_test";
    settings.RESULTS_DIR = settings.EXPERIMENT_NAME + "_results";
    //1x1 micrometer domain
    settings.CELL_NX = 1;//1;
    settings.CELL_NY = 1;//1;
    settings.CELL_DX = 2.5;//;0.5;//1.0;
    settings.CELL_DY = 1.0;//0.5;//1.0;

    //non ghosted complex
    settings.GHOSTED = false;

    //number of microtubules in the simulation
    settings.NUM_MT = 0;//18;

    //starting size of the MTs
    settings.MT_MIN_SEGMENT_INIT = 0.005;
    settings.MT_MAX_SEGMENT_INIT = 0.01;

    settings.LENGTH_DIV_FACTOR = 1.2;

    //actual MTs are 23 to 27 nm in diameter and up to 50 micrometers (um) long
    settings.DIV_LENGTH = 3*0.025; //3*0.025 //MT segments are approximated as 25 nm long
    settings.DIV_LENGTH_RETRACT = 0.0025;

    //growing and shrinking velocities in micrometers per second
    settings.V_PLUS = 4*0.0615; // 3.69um/min to um/s, shaw et al. (2003)
    settings.V_MINUS = 4*0.00883; // 0.53um/min to um/s, shaw et al. (2003)

    settings.SIGMOID_K = 10.0;

    //0.05 micrometers = 50 nanometers
    settings.MAXIMAL_REACTION_RADIUS = 2*0.05; // 2*0.05

    //simulation time in seconds
    settings.TOTAL_TIME = 20.0;//25.0;//20.0;
    settings.DELTA = 0.5/8.0; //unit of seconds

    //The internal step of the solver should be at least smaller than delta
    settings.DELTA_DELTA_T = settings.DELTA / 20.0;
    settings.DELTA_T_MIN = settings.DELTA_DELTA_T;
    settings.NUM_STEPS = settings.TOTAL_TIME / settings.DELTA;

    settings.RHO_TEST_RATE = 10.0;
}

int main() {

    std::cout << "Generating all the configuration files for the experiment" << std::endl;

    Parameters settings;
    set_default(settings);
    serialize(settings, "test.json");

    return 0;
}
