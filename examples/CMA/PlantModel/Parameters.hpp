#ifndef DGGML_PARAMETERS_HPP
#define DGGML_PARAMETERS_HPP
#include <string>
#include "simdjson.h"

//TODO: user and system params should be split, so that the user subclasses from a core params class
struct Parameters
{
    std::string EXPERIMENT_NAME;
    std::size_t CHECKPOINT_FREQUENCY;
    double DELTA;
    double DELTA_DELTA_T;
    double MIN_DELTA_STEPS;
    std::size_t CELL_NX;
    std::size_t CELL_NY;
    double CELL_DX;
    double CELL_DY;
    bool GHOSTED;
    std::size_t NUM_MT;
    double MT_MIN_SEGMENT_INIT;
    double MT_MAX_SEGMENT_INIT;
    std::size_t NUM_STEPS;
    double LENGTH_DIV_FACTOR;
    double DIV_LENGTH;
    double DIV_LENGTH_RETRACT;
    double V_PLUS;
    double V_MINUS;
    double SIGMOID_K;
    double TOTAL_TIME;
    double MAXIMAL_REACTION_RADIUS;
    double DELTA_T_MIN;
    double RHO_TEST_RATE; //a tunable test parameter for MT dynamics
    std::string RESULTS_DIR;
    bool CLASP_ENTRY;
    bool CLASP_EXIT;
//    bool ENABLE_SOLVING;
//    bool ENABLE_RETRACTION;
//    bool ENABLE_ZIPPERING;
//    bool ENABLE_CROSSOVER;
//    bool ENABLE_CATASTROPE;


    void set_parameters(simdjson::ondemand::document& interface)
    {
        // -----------------------
        // Core Parameter Settings
        // -----------------------
        std::string_view temp = interface["CORE"]["EXPERIMENT"];
        EXPERIMENT_NAME = static_cast<std::string>(temp);
        auto name = EXPERIMENT_NAME;
        temp = interface["CORE"]["RESULTS_DIR"];
        RESULTS_DIR = static_cast<std::string>(temp);
        TOTAL_TIME = double(interface["CORE"]["TOTAL_TIME"]);
        //Delta should be big, but not to big. In this case, the maximum amount of time it would
        //take one MT to grow a single unit of MT
        //e.g. something like: 0.25*settings.MAXIMAL_REACTION_RADIUS / std::max(settings.V_PLUS, settings.V_MINUS);
        DELTA = double(interface["CORE"]["DELTA"]);
        MIN_DELTA_STEPS = double(interface["CORE"]["MIN_DELTA_STEPS"]);
        //The internal step of the solver should be at least smaller than delta
        DELTA_DELTA_T = DELTA / MIN_DELTA_STEPS;
        DELTA_T_MIN = DELTA_DELTA_T;
        NUM_STEPS = TOTAL_TIME / DELTA;
        MAXIMAL_REACTION_RADIUS = double(interface["CORE"]["MAXIMAL_REACTION_RADIUS"]);
        CHECKPOINT_FREQUENCY = double(interface["CORE"]["CHECKPOINT_FREQUENCY"]);

        std::cout << "Core parameter settings parsed...\n";

        // ------------------------------
        // Expanded cell complex settings
        // ------------------------------
        CELL_NX = int64_t(interface["EXPANDED_CELL_COMPLEX"]["CELL_NX"]);
        CELL_NY = int64_t(interface["EXPANDED_CELL_COMPLEX"]["CELL_NY"]);
        CELL_DX = double(interface["EXPANDED_CELL_COMPLEX"]["CELL_DX"]);
        CELL_DY = double(interface["EXPANDED_CELL_COMPLEX"]["CELL_DY"]);
        GHOSTED = bool(interface["EXPANDED_CELL_COMPLEX"]["GHOSTED"]);
        std::cout << "Expanded cell complex settings parsed...\n";

        // -----------------------
        // Initialization settings
        // -----------------------
        NUM_MT = int64_t(interface["INITIALIZATION"]["NUM_MT"]);
        MT_MIN_SEGMENT_INIT = double(interface["INITIALIZATION"]["MT_MIN_SEGMENT_INIT"]);
        MT_MAX_SEGMENT_INIT = double(interface["INITIALIZATION"]["MT_MAX_SEGMENT_INIT"]);
        std::cout << "Initialization settings parsed...\n";

        // ----------------
        // Grammar settings
        // ----------------
        // Growth rules
        V_PLUS = double(interface["GRAMMAR"]["GROWTH_RULES"]["GROWTH_RATE"]);
        LENGTH_DIV_FACTOR = double(interface["GRAMMAR"]["GROWTH_RULES"]["LENGTH_DIV_FACTOR"]);
        //actual MTs are 23 to 27 nm in diameter and up to 50 micrometers (um) long
        DIV_LENGTH = double(interface["GRAMMAR"]["GROWTH_RULES"]["DIV_LENGTH"]);
        std::cout << "Growth rule settings parsed...\n";

        // Retraction rules
        DIV_LENGTH_RETRACT = double(interface["GRAMMAR"]["RETRACTION_RULES"]["DIV_LENGTH_RETRACT"]);
        V_MINUS = double(interface["GRAMMAR"]["RETRACTION_RULES"]["RETRACTION_RATE"]);
        std::cout << "Retraction rule settings parsed...\n";

        // Clasp rules
        CLASP_ENTRY = bool(interface["GRAMMAR"]["CLASP_RULES"]["ENABLE_ENTRY"]);
        CLASP_EXIT = bool(interface["GRAMMAR"]["CLASP_RULES"]["ENABLE_EXIT"]);
        std::cout << "Clasp rule settings parsed...\n";

        std::cout << "Grammar settings parsed...\n";
        SIGMOID_K = double(interface["SIGMOID_K"]);
        RHO_TEST_RATE = double(interface["RHO_TEST_RATE"]);
    }
};

#endif //DGGML_PARAMETERS_HPP
