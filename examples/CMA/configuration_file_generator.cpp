#include <iostream>
#include <fstream>
#include <string>
#include <filesystem>
#include <vector>

struct Parameters
{
    // Core parameter settings
    std::string EXPERIMENT_NAME;
    std::string RESULTS_DIR;
    double TOTAL_TIME{};
    double DELTA{};
    double MIN_DELTA_STEPS{};
    double DELTA_DELTA_T{};
    double DELTA_T_MIN{};
    std::size_t NUM_STEPS{};
    double MAXIMAL_REACTION_RADIUS{};
    std::size_t CHECKPOINT_FREQUENCY{};

    // Expanded cell complex settings
    std::size_t CELL_NX{};
    std::size_t CELL_NY{};
    double CELL_DX{};
    double CELL_DY{};
    bool GHOSTED{};

    // Initialization settings
    std::size_t NUM_MT{};
    double MT_MIN_SEGMENT_INIT{};
    double MT_MAX_SEGMENT_INIT{};

    // Growth rule settings
    bool ENABLE_GROWTH{};
    double WITH_GROWTH_RATE_FACTOR{};
    double V_PLUS{};
    double LENGTH_DIV_FACTOR{};
    double DIV_LENGTH{};
    bool ENABLE_WOBBLE{};
    double WOBBLE_ANGLE{};

    // Retraction rule settings
    bool ENABLE_RETRACTION{};
    double WITH_RETRACTION_RATE_FACTOR{};
    double DIV_LENGTH_RETRACT{};
    double V_MINUS{};

    // Boundary rule settings
    bool ENABLE_STANDARD_BOUNDARY_CATASTROPHE{};
    double STANDARD_BOUNDARY_CATASTROPHE_RATE{};
    bool ENABLE_CLASP_BOUNDARY_CATASTROPHE{};
    double CLASP_BOUNDARY_CATASTROPHE_RATE{};

    // Catastrophe rule settings
    bool ENABLE_INTERMEDIATE_CIC{};
    double INTERMEDIATE_CIC_RATE{};
    bool ENABLE_POSITIVE_CIC{};
    double POSITIVE_CIC_RATE{};
    bool ENABLE_NEGATIVE_CIC{};
    double NEGATIVE_CIC_RATE{};
    double CATASTROPHE_ANGLE{};

    // Zippering rule settings
    bool ENABLE_ZIPPERING{};
    double ZIPPERING_HIT_RATE{};
    double ZIPPERING_GUARD_RATE{};
    double ZIPPERING_RETRACTION_RATE{};
    double CRITICAL_ANGLE{};
    double SEPARATION_DISTANCE{};

    // Crossover rule settings
    bool ENABLE_CROSSOVER{};
    double CROSSOVER_RATE{};
    double CROSSOVER_ANGLE{};
    bool ENABLE_UNCROSSOVER{};
    double UNCROSSOVER_RATE{};

    // Clasp rule settings
    bool CLASP_ENABLE_ENTRY{};
    double CLASP_ENTRY_RATE{};
    double CLASP_ENTRY_ANGLE{};
    bool CLASP_ENABLE_EXIT{};
    double CLASP_EXIT_RATE{};
    double CLASP_EXIT_ANGLE{};
    bool CLASP_ENABLE_CAT{};
    double CLASP_CAT_RATE{};
    bool CLASP_ENABLE_DETACHMENT{};
    double CLASP_DETACHMENT_RATE{};

    // Destruction rule settings
    bool ENABLE_MT_DESTRUCTION{};
    double MT_DESTRUCTION_RATE{};

    // Creation rule settings
    bool ENABLE_CREATION{};
    double CREATION_RATE{};
    double CREATION_FACTOR{};

    // Recovery rule settings
    bool ENABLE_RECOVERY{};
    double RECOVERY_RATE{};
    double RECOVERY_FACTOR{};

    // Other
    double RHO_TEST_RATE{}; //a tunable test parameter for MT dynamics
    double SIGMOID_K{};
    double COLLISION_DISTANCE{};

    //function to set the default parameters that we will modifiy for experiments
    void set_default()
    {
        // -----------------------
        // Core Parameter Settings
        // -----------------------

        EXPERIMENT_NAME = "my_test";
        auto name = EXPERIMENT_NAME;
        RESULTS_DIR = "my_test_results";
        TOTAL_TIME = 3600;
        //Delta should be big, but not to big. In this case, the maximum amount of time it would
        //take one MT to grow a single unit of MT
        //e.g. something like: 0.25*settings.MAXIMAL_REACTION_RADIUS / std::max(settings.V_PLUS, settings.V_MINUS);
        DELTA = 0.4;//0.4;//0.0625;
        MIN_DELTA_STEPS = 5;
        //The internal step of the solver should be at least smaller than delta
        DELTA_DELTA_T = DELTA / MIN_DELTA_STEPS;
        DELTA_T_MIN = DELTA_DELTA_T;
        NUM_STEPS = TOTAL_TIME / DELTA;
        MAXIMAL_REACTION_RADIUS = 0.1;
        CHECKPOINT_FREQUENCY = 30;

        //std::cout << "Core parameter settings parsed...\n";

        // ------------------------------
        // Expanded cell complex settings
        // ------------------------------
        CELL_NX = 3;
        CELL_NY = 3;
        CELL_DX = 1.66;//5.5;//1.666;
        CELL_DY = 1.66;//2.0;//1.666;
        GHOSTED = false;
        //std::cout << "Expanded cell complex settings parsed...\n";

        // -----------------------
        // Initialization settings
        // -----------------------
        NUM_MT = 0;
        MT_MIN_SEGMENT_INIT = 0.005;
        MT_MAX_SEGMENT_INIT = 0.01;
        //std::cout << "Initialization settings parsed...\n";

        // ----------------
        // Grammar settings
        // ----------------
        // Growth rules
        ENABLE_GROWTH = true;
        WITH_GROWTH_RATE_FACTOR = 100.0;
        V_PLUS = 0.0615;
        //actual MTs are 23 to 27 nm in diameter and up to 50 micrometers (um) long
        DIV_LENGTH = 0.075;
        LENGTH_DIV_FACTOR = 1.2;
        ENABLE_WOBBLE = false;
        WOBBLE_ANGLE = 8.0;
        //std::cout << "Growth rule settings parsed...\n";

        // Retraction rules
        ENABLE_RETRACTION = true;
        WITH_RETRACTION_RATE_FACTOR = 10.0;
        V_MINUS = 0.00883;
        DIV_LENGTH_RETRACT = 0.0025;
        //std::cout << "Retraction rule settings parsed...\n";

        // Boundary rules
        ENABLE_STANDARD_BOUNDARY_CATASTROPHE =  true;
        STANDARD_BOUNDARY_CATASTROPHE_RATE = 40000.0;
        ENABLE_CLASP_BOUNDARY_CATASTROPHE = true;
        CLASP_BOUNDARY_CATASTROPHE_RATE = 40000.0;
        //std::cout << "Boundary rule settings parsed...\n";

        // Catastrophe rule settings
        ENABLE_INTERMEDIATE_CIC = true;
        INTERMEDIATE_CIC_RATE = 4000.0;
        ENABLE_POSITIVE_CIC = true;
        POSITIVE_CIC_RATE = 1000.0;
        ENABLE_NEGATIVE_CIC = true;
        NEGATIVE_CIC_RATE = 4000.0;
        CATASTROPHE_ANGLE = 50.0;
        //std::cout << "Catastrophe rule settings parsed...\n";

        // Zippering rule settings
        ENABLE_ZIPPERING = true;
        ZIPPERING_HIT_RATE = 4000.0;
        ZIPPERING_GUARD_RATE = 4000.0;
        ZIPPERING_RETRACTION_RATE = 10.0;
        CRITICAL_ANGLE = 50.0;
        SEPARATION_DISTANCE = 0.025;
        //std::cout << "Zippering rule settings parsed...\n";

        // Crossover rule settings
        ENABLE_CROSSOVER = true;
        CROSSOVER_RATE = 40.0;
        CROSSOVER_ANGLE = 50.0;
        ENABLE_UNCROSSOVER = true;
        UNCROSSOVER_RATE = 0.01; // a rate of one => occurs once per unit of time
        //in this case, 1 => onces per second on average
        //std::cout << "Crossover rule settings parsed...\n";

        // Clasp rule settings
        CLASP_ENABLE_ENTRY = false;
        CLASP_ENTRY_RATE = 0.0005;
        CLASP_ENTRY_ANGLE = 15.0;
        CLASP_ENABLE_EXIT = false;
        CLASP_EXIT_RATE = 40000.0;
        CLASP_EXIT_ANGLE = 15.0;
        CLASP_ENABLE_CAT = true;
        CLASP_CAT_RATE = 1000.0;
        CLASP_ENABLE_DETACHMENT = false;
        CLASP_DETACHMENT_RATE = 0.01;
        //std::cout << "Clasp rule settings parsed...\n";

        // Destruction rule settings
        ENABLE_MT_DESTRUCTION = true;
        MT_DESTRUCTION_RATE = 10.0;
        //std::cout << "Destruction rule settings parsed...\n";

        // Creation rule settings
        ENABLE_CREATION = true;
        CREATION_RATE = 0.0026;
        CREATION_FACTOR = 1.0;
        //std::cout << "Creation rule settings parsed...\n";

        // Recovery rule settings
        ENABLE_RECOVERY = true;
        RECOVERY_RATE = 0.016;
        RECOVERY_FACTOR = 1.0;
        //std::cout << "Recovery rule settings parsed...\n";

        //std::cout << "Grammar settings parsed...\n";

        RHO_TEST_RATE = 10.0;
        SIGMOID_K = 10.0;
        COLLISION_DISTANCE = 0.025;
        //std::cout << "Other settings parsed...\n";
    }
};

void serialize(Parameters& settings, std::string filename)
{
    std::ofstream outfile(filename);

    // Check if the file was opened successfully
    if (!outfile.is_open()) {
        std::cerr << "Error opening file for writing." << std::endl;
    }

    //Begin the json file
    outfile << "{\n";

    // Core paramaters
    outfile << "\t\"CORE\": {\n";
    outfile << "\t\t\"EXPERIMENT\": \"" <<  settings.EXPERIMENT_NAME << "\",\n";
    outfile << "\t\t\"RESULTS_DIR\": \"" <<  settings.RESULTS_DIR << "\",\n";
    outfile << "\t\t\"TOTAL_TIME\": " << settings.TOTAL_TIME << ",\n";
    outfile << "\t\t\"DELTA\": " << settings.DELTA << ",\n";
    outfile << "\t\t\"MIN_DELTA_STEPS\": " << settings.MIN_DELTA_STEPS << ",\n";
    outfile << "\t\t\"MAXIMAL_REACTION_RADIUS\": " << settings.MAXIMAL_REACTION_RADIUS << ",\n";
    outfile << "\t\t\"CHECKPOINT_FREQUENCY\": " << settings.CHECKPOINT_FREQUENCY << "\n";
    outfile << "\t},\n";

    // Expanded cell com settings
    outfile << "\t\"EXPANDED_CELL_COMPLEX\": {\n";
    outfile << "\t\t\"CELL_NX\": " << settings.CELL_NX << ",\n";
    outfile << "\t\t\"CELL_NY\": " << settings.CELL_NY << ",\n";
    outfile << "\t\t\"CELL_DX\": " << settings.CELL_DX << ",\n";
    outfile << "\t\t\"CELL_DY\": " << settings.CELL_DY << ",\n";
    settings.GHOSTED ? outfile << "\t\t\"GHOSTED\": true\n" : outfile << "\t\t\"GHOSTED\": false\n";
    outfile << "\t},\n";

    //initialization settings
    outfile << "\t\"INITIALIZATION\": {\n";
    outfile << "\t\t\"NUM_MT\": " << settings.NUM_MT << ",\n";
    outfile << "\t\t\"MT_MIN_SEGMENT_INIT\": " << settings.MT_MIN_SEGMENT_INIT << ",\n";
    outfile << "\t\t\"MT_MAX_SEGMENT_INIT\": " << settings.MT_MAX_SEGMENT_INIT << "\n";
    outfile << "\t},\n";

    //Grammar settings
    outfile << "\t\"GRAMMAR\": {\n";

    //Growth rules
    outfile << "\t\t\"GROWTH_RULES\": {\n";
    settings.ENABLE_GROWTH ? outfile << "\t\t\t\"ENABLE_GROWTH\": true,\n" : outfile << "\t\t\t\"ENABLE_GROWTH\": false,\n";
    outfile << "\t\t\t\"WITH_GROWTH_RATE_FACTOR\": " << settings.WITH_GROWTH_RATE_FACTOR << ",\n";
    outfile << "\t\t\t\"GROWTH_VELOCITY\": " << settings.V_PLUS << ",\n";
    outfile << "\t\t\t\"DIV_LENGTH\": " << settings.DIV_LENGTH << ",\n";
    outfile << "\t\t\t\"LENGTH_DIV_FACTOR\": " << settings.LENGTH_DIV_FACTOR << ",\n";
    settings.ENABLE_WOBBLE ? outfile << "\t\t\t\"ENABLE_WOBBLE\": true,\n" : outfile << "\t\t\t\"ENABLE_WOBBLE\": false,\n";
    outfile << "\t\t\t\"WOBBLE_ANGLE\": " << settings.WOBBLE_ANGLE << "\n";
    outfile << "\t\t},\n";

    //Retraction rule
    outfile << "\t\t\"RETRACTION_RULES\": {\n";
    settings.ENABLE_RETRACTION ? outfile << "\t\t\t\"ENABLE_RETRACTION\": true,\n" : outfile << "\t\t\t\"ENABLE_RETRACTION\": false,\n";
    outfile << "\t\t\t\"WITH_RETRACTION_RATE_FACTOR\": " << settings.WITH_RETRACTION_RATE_FACTOR << ",\n";
    outfile << "\t\t\t\"RETRACTION_VELOCITY\": " << settings.V_MINUS << ",\n";
    outfile << "\t\t\t\"DIV_LENGTH_RETRACT\": " << settings.DIV_LENGTH_RETRACT << "\n";
    outfile << "\t\t},\n";

    //Boundary rules
    outfile << "\t\t\"BOUNDARY_RULES\": {\n";
    settings.ENABLE_STANDARD_BOUNDARY_CATASTROPHE ? outfile << "\t\t\t\"ENABLE_STANDARD\": true,\n" : outfile << "\t\t\t\"ENABLE_STANDARD\": false,\n";
    outfile << "\t\t\t\"STANDARD_RATE\": " << settings.STANDARD_BOUNDARY_CATASTROPHE_RATE << ",\n";
    settings.ENABLE_CLASP_BOUNDARY_CATASTROPHE ? outfile << "\t\t\t\"ENABLE_CLASP\": true,\n" : outfile << "\t\t\t\"ENABLE_CLASP\": false,\n";
    outfile << "\t\t\t\"CLASP_RATE\": " << settings.CLASP_BOUNDARY_CATASTROPHE_RATE << "\n";
    outfile << "\t\t},\n";

    //Catastrophe rules
    outfile << "\t\t\"CATASTROPHE_RULES\": {\n";
    settings.ENABLE_INTERMEDIATE_CIC ? outfile << "\t\t\t\"ENABLE_INTERMEDIATE_CIC\": true,\n" : outfile << "\t\t\t\"ENABLE_INTERMEDIATE_CIC\": false,\n";
    outfile << "\t\t\t\"INTERMEDIATE_CIC_RATE\": " << settings.INTERMEDIATE_CIC_RATE << ",\n";
    settings.ENABLE_POSITIVE_CIC ? outfile << "\t\t\t\"ENABLE_POSITIVE_CIC\": true,\n" : outfile << "\t\t\t\"ENABLE_POSITIVE_CIC\": false,\n";
    outfile << "\t\t\t\"POSITIVE_CIC_RATE\": " << settings.POSITIVE_CIC_RATE << ",\n";
    settings.ENABLE_NEGATIVE_CIC ? outfile << "\t\t\t\"ENABLE_NEGATIVE_CIC\": true,\n" : outfile << "\t\t\t\"ENABLE_NEGATIVE_CIC\": false,\n";
    outfile << "\t\t\t\"NEGATIVE_CIC_RATE\": " << settings.NEGATIVE_CIC_RATE << ",\n";
    outfile << "\t\t\t\"CATASTROPHE_ANGLE\": " << settings.CATASTROPHE_ANGLE << "\n";
    outfile << "\t\t},\n";

    //Zippering rules
    outfile << "\t\t\"ZIPPERING_RULES\": {\n";
    settings.ENABLE_ZIPPERING ? outfile << "\t\t\t\"ENABLE_ZIPPERING\": true,\n" : outfile << "\t\t\t\"ENABLE_ZIPPERING\": false,\n";
    outfile << "\t\t\t\"ZIPPERING_HIT_RATE\": " << settings.ZIPPERING_HIT_RATE << ",\n";
    outfile << "\t\t\t\"ZIPPERING_GUARD_RATE\": " << settings.ZIPPERING_GUARD_RATE << ",\n";
    outfile << "\t\t\t\"ZIPPERING_RETRACTION_RATE\": " << settings.ZIPPERING_RETRACTION_RATE << ",\n";
    outfile << "\t\t\t\"CRITICAL_ANGLE\": " << settings.CRITICAL_ANGLE << ",\n";
    outfile << "\t\t\t\"SEPARATION_DISTANCE\": " << settings.SEPARATION_DISTANCE << "\n";
    outfile << "\t\t},\n";

    //Crossover rules
    outfile << "\t\t\"CROSSOVER_RULES\": {\n";
    settings.ENABLE_CROSSOVER ? outfile << "\t\t\t\"ENABLE_CROSSOVER\": true,\n" : outfile << "\t\t\t\"ENABLE_CROSSOVER\": false,\n";
    outfile << "\t\t\t\"CROSSOVER_RATE\": " << settings.CROSSOVER_RATE << ",\n";
    outfile << "\t\t\t\"CROSSOVER_ANGLE\": " << settings.CROSSOVER_ANGLE << ",\n";
    settings.ENABLE_UNCROSSOVER ? outfile << "\t\t\t\"ENABLE_UNCROSSOVER\": true,\n" : outfile << "\t\t\t\"ENABLE_UNCROSSOVER\": false,\n";
    outfile << "\t\t\t\"UNCROSSOVER_RATE\": " << settings.UNCROSSOVER_RATE << "\n";
    outfile << "\t\t},\n";

    //Clasp rules
    outfile << "\t\t\"CLASP_RULES\": {\n";
    settings.CLASP_ENABLE_ENTRY ? outfile << "\t\t\t\"CLASP_ENABLE_ENTRY\": true,\n" : outfile << "\t\t\t\"CLASP_ENABLE_ENTRY\": false,\n";
    outfile << "\t\t\t\"CLASP_ENTRY_RATE\": " << settings.CLASP_ENTRY_RATE << ",\n";
    outfile << "\t\t\t\"CLASP_ENTRY_ANGLE\": " << settings.CLASP_ENTRY_ANGLE << ",\n";
    settings.CLASP_ENABLE_EXIT ? outfile << "\t\t\t\"CLASP_ENABLE_EXIT\": true,\n" : outfile << "\t\t\t\"CLASP_ENABLE_EXIT\": false,\n";
    outfile << "\t\t\t\"CLASP_EXIT_RATE\": " << settings.CLASP_EXIT_RATE << ",\n";
    outfile << "\t\t\t\"CLASP_EXIT_ANGLE\": " << settings.CLASP_EXIT_ANGLE << ",\n";
    settings.CLASP_ENABLE_CAT ? outfile << "\t\t\t\"CLASP_ENABLE_CAT\": true,\n" : outfile << "\t\t\t\"CLASP_ENABLE_CAT\": false,\n";
    outfile << "\t\t\t\"CLASP_CAT_RATE\": " << settings.CLASP_CAT_RATE << ",\n";
    settings.CLASP_ENABLE_DETACHMENT ? outfile << "\t\t\t\"CLASP_ENABLE_DETACHMENT\": true,\n" : outfile << "\t\t\t\"CLASP_ENABLE_DETACHMENT\": false,\n";
    outfile << "\t\t\t\"CLASP_DETACHMENT_RATE\": " << settings.CLASP_DETACHMENT_RATE << "\n";
    outfile << "\t\t},\n";

    //Destruction rules
    outfile << "\t\t\"DESTRUCTION_RULES\": {\n";
    settings.ENABLE_MT_DESTRUCTION ? outfile << "\t\t\t\"ENABLE_MT_DESTRUCTION\": true,\n" : outfile << "\t\t\t\"ENABLE_MT_DESTRUCTION\": false,\n";
    outfile << "\t\t\t\"MT_DESTRUCTION_RATE\": " << settings.MT_DESTRUCTION_RATE << "\n";
    outfile << "\t\t},\n";

    //Creation rules
    outfile << "\t\t\"CREATION_RULES\": {\n";
    settings.ENABLE_CREATION ? outfile << "\t\t\t\"ENABLE_CREATION\": true,\n" : outfile << "\t\t\t\"ENABLE_CREATION\": false,\n";
    outfile << "\t\t\t\"CREATION_RATE\": " << settings.CREATION_RATE << ",\n";
    outfile << "\t\t\t\"CREATION_FACTOR\": " << settings.CREATION_FACTOR << "\n";
    outfile << "\t\t},\n";

    //Recovery rules
    outfile << "\t\t\"RECOVERY_RULES\": {\n";
    settings.ENABLE_RECOVERY ? outfile << "\t\t\t\"ENABLE_RECOVERY\": true,\n" : outfile << "\t\t\t\"ENABLE_RECOVERY\": false,\n";
    outfile << "\t\t\t\"RECOVERY_RATE\": " << settings.RECOVERY_RATE << ",\n";
    outfile << "\t\t\t\"RECOVERY_FACTOR\": " << settings.RECOVERY_FACTOR << "\n";
    outfile << "\t\t}\n";

    //end grammar rules
    outfile << "\t},\n";


    //other settings
    outfile << "\t\"RHO_TEST_RATE\": " << settings.RHO_TEST_RATE << ",\n";
    outfile << "\t\"SIGMOID_K\": " << settings.SIGMOID_K << ",\n";
    outfile << "\t\"COLLISION_DISTANCE\": " << settings.COLLISION_DISTANCE << "\n";

    // end the json file
    outfile << "}\n";

    outfile.close();
}

void create_slurm_script(std::string path, std::string filename, std::size_t n)
{
    std::string slurm_name = filename + "_slurm.sh";
    std::ofstream outfile(path + "/"+slurm_name);
    outfile << "#!/bin/bash\n";

    // Check if the file was opened successfully
    if (!outfile.is_open()) {
        std::cerr << "Error opening file for writing." << std::endl;
    }
    outfile << "#SBATCH --job-name " << filename << "_ensemble   ## name that will show up in the queue\n";
    outfile << "#SBATCH --output slurm-%j.out   ## filename of the output; the %j is equal to jobID; default is slurm-[jobID].out\n";
    outfile << "#SBATCH --ntasks=" << n << " ## number of tasks (analyses) to run\n";
    outfile << "#SBATCH --cpus-per-task=1  ## the number of threads allocated to each task\n";
    outfile << "#SBATCH --mem-per-cpu=5000M   # memory per CPU core\n";
    outfile << "#SBATCH --partition=ai4science.p  ## the partitions to run in (comma seperated)\n";
    outfile << "#SBATCH --nodelist=blazar-1 ## Specifies the node list\n";
    outfile << "#SBATCH --time=1-00:00:00  ## time for analysis (day-hour:min:sec)\n";
    outfile << "#SBATCH --hint=nomultithread\n";
    outfile << "## Load modules\n";

    outfile << "## Your commands here\n";
    outfile << "mkdir output\n";

    outfile << "## Insert code, and run your programs here (use 'srun').\n";
    for(auto i = 1; i <= n; i++)
        outfile << "srun --ntasks=1 --cpus-per-task=1 ../../mt_dgg_simulator " << filename << "_" << i <<".json > output/out" << i <<".txt &\n";

    outfile << "wait\n";

    outfile.close();

}

//will set the file names, we just need set the other settings we want beforehand
void create_experiment(Parameters& settings, std::string root_dir, std::string experiment_name, std::size_t num_runs)
{
    std::string experiment_dir = root_dir + "/" + experiment_name;
    std::filesystem::create_directory(root_dir + "/" + experiment_name);

    for(int i = 1; i <= num_runs; i++) {
        std::string exp_name = experiment_name+"_"+std::to_string(i);
        settings.EXPERIMENT_NAME = exp_name;
        std::string results_dir = exp_name+"_results";
        settings.RESULTS_DIR = results_dir;
        std::string config_name = exp_name + ".json";
        std::string path = root_dir + "/" + experiment_name + "/" + config_name;
        serialize(settings, path);
    }
    create_slurm_script(experiment_dir, experiment_name, num_runs);
}

void create_main_bash(std::string root_dir, std::vector<std::string>& filenames)
{
    std::ofstream outfile(root_dir+"/"+"run_all.sh");
    // Check if the file was opened successfully
    if (!outfile.is_open()) {
        std::cerr << "Error opening file for writing." << std::endl;
    }
    outfile << "#!/bin/bash\n";

    for(auto& filename : filenames) {
        std::string slurm_name = filename + "_slurm.sh";
        outfile << "cd " << filename << "\n";
        outfile << "chmod +x " << filename << "\n";
        outfile << "sbatch " << slurm_name << "\n";
        outfile << "cd ..\n";
    }

    outfile.close();
}

int main() {

    std::cout << "Generating all the configuration files for the experiment" << std::endl;

    std::vector<std::string> all_filenames;
    Parameters settings;
    settings.set_default();

    std::cout << "Creating the main experiment directory and removing if it exists...\n";
    std::string root_dir = "exp";
    std::filesystem::remove_all(root_dir);
    std::filesystem::create_directory(root_dir);

    std::size_t n = 16; //number of runs for the ensemble
    create_experiment(settings, root_dir, "square_no_clasp",n);
    all_filenames.push_back("square_no_clasp");

    settings.CREATION_FACTOR = 0.2;
    create_experiment(settings, root_dir, "square_no_clasp_low_creation",n);
    all_filenames.push_back("square_no_clasp_low_creation");

    settings.CREATION_FACTOR = 1.0;
    settings.ENABLE_CROSSOVER = true;
    settings.CROSSOVER_RATE = 4000.0;
    create_experiment(settings, root_dir, "square_no_clasp_high_cross",n);
    all_filenames.push_back("square_no_clasp_high_cross");

    settings.set_default();
    settings.CLASP_ENABLE_ENTRY = true;
    settings.CLASP_ENTRY_RATE = 0.0005;
    settings.CLASP_ENABLE_EXIT = true;
    for(int i = 1; i <= 4; i ++)
    {
        settings.CLASP_EXIT_ANGLE = (double)(i)*15.0;
        std::string name = "square_with_clasp_angle_"+std::to_string((int)settings.CLASP_EXIT_ANGLE);
        create_experiment(settings, root_dir, name,n);
        all_filenames.push_back(name);
    }
    settings.CLASP_EXIT_ANGLE = 15.0;
    settings.CLASP_ENTRY_RATE = 0.005;
    create_experiment(settings, root_dir, "square_with_clasp_angle_"+std::to_string((int)settings.CLASP_EXIT_ANGLE)+"_influx",n);
    all_filenames.push_back("square_with_clasp_angle_"+std::to_string((int)settings.CLASP_EXIT_ANGLE)+"_influx");

    //resets to default
    settings.set_default();
    settings.CELL_NX = 1;
    settings.CELL_NY = 1;
    settings.CELL_DX = 8.3;
    settings.CELL_DY = 3.0;
    create_experiment(settings, root_dir, "rectangle_no_clasp", n);
    all_filenames.push_back("rectangle_no_clasp");
    settings.ENABLE_CROSSOVER = true;
    settings.CROSSOVER_RATE = 4000.0;
    create_experiment(settings, root_dir, "rectangle_no_clasp_high_cross",n);
    all_filenames.push_back("rectangle_no_clasp_high_cross");

    settings.CROSSOVER_RATE = 40.0;
    settings.CLASP_ENABLE_ENTRY = true;
    settings.CLASP_ENTRY_RATE = 0.0005;
    settings.CLASP_ENABLE_EXIT = true;
    for(int i = 1; i <= 4; i ++)
    {
        settings.CLASP_EXIT_ANGLE = (double)(i)*15.0;
        std::string name = "rectangle_with_clasp_angle_"+std::to_string((int)settings.CLASP_EXIT_ANGLE);
        create_experiment(settings, root_dir, name,n);
        all_filenames.push_back(name);
    }
    settings.CLASP_EXIT_ANGLE = 15.0;
    settings.CLASP_ENTRY_RATE = 0.005;
    create_experiment(settings, root_dir, "rectangle_with_clasp_angle_"+std::to_string((int)settings.CLASP_EXIT_ANGLE)+"_influx",n);
    all_filenames.push_back("rectangle_with_clasp_angle_"+std::to_string((int)settings.CLASP_EXIT_ANGLE)+"_influx");

    create_main_bash(root_dir, all_filenames);
    return 0;
}
