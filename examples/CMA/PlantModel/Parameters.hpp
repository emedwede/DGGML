#ifndef DGGML_PARAMETERS_HPP
#define DGGML_PARAMETERS_HPP
#include <string>
#include "simdjson.h"

//TODO: user and system params should be split, so that the user subclasses from a core params class
struct Parameters
{
    // Core parameter settings
    std::string EXPERIMENT_NAME;
    std::string RESULTS_DIR;
    double TOTAL_TIME;
    double DELTA;
    double MIN_DELTA_STEPS;
    double DELTA_DELTA_T;
    double DELTA_T_MIN;
    std::size_t NUM_STEPS;
    double MAXIMAL_REACTION_RADIUS;
    std::size_t CHECKPOINT_FREQUENCY;

    // Expanded cell complex settings
    std::size_t CELL_NX;
    std::size_t CELL_NY;
    double CELL_DX;
    double CELL_DY;
    bool GHOSTED;

    // Initialization settings
    std::size_t NUM_MT;
    double MT_MIN_SEGMENT_INIT;
    double MT_MAX_SEGMENT_INIT;

    // Growth rule settings
    bool ENABLE_GROWTH;
    double WITH_GROWTH_RATE_FACTOR;
    double V_PLUS;
    double LENGTH_DIV_FACTOR;
    double DIV_LENGTH;

    // Retraction rule settings
    bool ENABLE_RETRACTION;
    double WITH_RETRACTION_RATE_FACTOR;
    double DIV_LENGTH_RETRACT;
    double V_MINUS;

    // Boundary rule settings
    bool ENABLE_STANDARD_BOUNDARY_CATASTROPHE;
    double STANDARD_BOUNDARY_CATASTROPHE_RATE;
    bool ENABLE_CLASP_BOUNDARY_CATASTROPHE;
    double CLASP_BOUNDARY_CATASTROPHE_RATE;

    // Catastrophe rule settings
    bool ENABLE_INTERMEDIATE_CIC;
    double INTERMEDIATE_CIC_RATE;
    bool ENABLE_POSITIVE_CIC;
    double POSITIVE_CIC_RATE;
    bool ENABLE_NEGATIVE_CIC;
    double NEGATIVE_CIC_RATE;
    double CATASTROPHE_ANGLE;

    // Crossover rule settings
    bool ENABLE_CROSSOVER;
    double CROSSOVER_RATE;
    double CROSSOVER_ANGLE;

    // Zippering rule settings
    bool ENABLE_ZIPPERING;
    double ZIPPERING_HIT_RATE;
    double ZIPPERING_GUARD_RATE;
    double ZIPPERING_RETRACTION_RATE;
    double CRITICAL_ANGLE;
    double SEPARATION_DISTANCE;

    // Clasp rule settings
    bool CLASP_ENABLE_ENTRY;
    double CLASP_ENTRY_RATE;
    double CLASP_ENTRY_ANGLE;
    bool CLASP_ENABLE_EXIT;
    double CLASP_EXIT_RATE;
    double CLASP_EXIT_ANGLE;
    bool CLASP_ENABLE_CAT;
    double CLASP_CAT_RATE;

    // Destruction rule settings
    bool ENABLE_MT_DESTRUCTION;
    double MT_DESTRUCTION_RATE;

    // Creation rule settings
    bool ENABLE_CREATION;
    double CREATION_RATE;
    double CREATION_FACTOR;

    // Recovery rule settings
    bool ENABLE_RECOVERY;
    double RECOVERY_RATE;
    double RECOVERY_FACTOR;

    // Other
    double SIGMOID_K;
    double RHO_TEST_RATE; //a tunable test parameter for MT dynamics
    double COLLISION_DISTANCE;


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
        ENABLE_GROWTH = bool(interface["GRAMMAR"]["GROWTH_RULES"]["ENABLE_GROWTH"]);
        WITH_GROWTH_RATE_FACTOR = double(interface["GRAMMAR"]["GROWTH_RULES"]["WITH_GROWTH_RATE_FACTOR"]);
        V_PLUS = double(interface["GRAMMAR"]["GROWTH_RULES"]["GROWTH_VELOCITY"]);
        LENGTH_DIV_FACTOR = double(interface["GRAMMAR"]["GROWTH_RULES"]["LENGTH_DIV_FACTOR"]);
        //actual MTs are 23 to 27 nm in diameter and up to 50 micrometers (um) long
        DIV_LENGTH = double(interface["GRAMMAR"]["GROWTH_RULES"]["DIV_LENGTH"]);
        std::cout << "Growth rule settings parsed...\n";

        // Retraction rules
        ENABLE_RETRACTION = bool(interface["GRAMMAR"]["RETRACTION_RULES"]["ENABLE_RETRACTION"]);
        WITH_RETRACTION_RATE_FACTOR = double(interface["GRAMMAR"]["RETRACTION_RULES"]["WITH_RETRACTION_RATE_FACTOR"]);
        DIV_LENGTH_RETRACT = double(interface["GRAMMAR"]["RETRACTION_RULES"]["DIV_LENGTH_RETRACT"]);
        V_MINUS = double(interface["GRAMMAR"]["RETRACTION_RULES"]["RETRACTION_VELOCITY"]);
        std::cout << "Retraction rule settings parsed...\n";

        // Boundary rules
        ENABLE_STANDARD_BOUNDARY_CATASTROPHE =  bool(interface["GRAMMAR"]["BOUNDARY_RULES"]["ENABLE_STANDARD"]);
        STANDARD_BOUNDARY_CATASTROPHE_RATE = double(interface["GRAMMAR"]["BOUNDARY_RULES"]["STANDARD_RATE"]);
        ENABLE_CLASP_BOUNDARY_CATASTROPHE = bool(interface["GRAMMAR"]["BOUNDARY_RULES"]["ENABLE_CLASP"]);
        CLASP_BOUNDARY_CATASTROPHE_RATE = double(interface["GRAMMAR"]["BOUNDARY_RULES"]["CLASP_RATE"]);
        std::cout << "Boundary rule settings parsed...\n";

        // Catastrophe rule settings
        ENABLE_INTERMEDIATE_CIC = bool(interface["GRAMMAR"]["CATASTROPHE_RULES"]["ENABLE_INTERMEDIATE_CIC"]);
        INTERMEDIATE_CIC_RATE = double(interface["GRAMMAR"]["CATASTROPHE_RULES"]["INTERMEDIATE_CIC_RATE"]);
        ENABLE_POSITIVE_CIC = bool(interface["GRAMMAR"]["CATASTROPHE_RULES"]["ENABLE_POSITIVE_CIC"]);
        POSITIVE_CIC_RATE = double(interface["GRAMMAR"]["CATASTROPHE_RULES"]["POSITIVE_CIC_RATE"]);
        ENABLE_NEGATIVE_CIC = bool(interface["GRAMMAR"]["CATASTROPHE_RULES"]["ENABLE_NEGATIVE_CIC"]);
        NEGATIVE_CIC_RATE = double(interface["GRAMMAR"]["CATASTROPHE_RULES"]["NEGATIVE_CIC_RATE"]);
        CATASTROPHE_ANGLE = double(interface["GRAMMAR"]["CATASTROPHE_RULES"]["CATASTROPHE_ANGLE"]);
        std::cout << "Catastrophe rule settings parsed...\n";

        // Crossover rule settings
        ENABLE_CROSSOVER = bool(interface["GRAMMAR"]["CROSSOVER_RULES"]["ENABLE_CROSSOVER"]);
        CROSSOVER_RATE = double(interface["GRAMMAR"]["CROSSOVER_RULES"]["CROSSOVER_RATE"]);
        CROSSOVER_ANGLE = double(interface["GRAMMAR"]["CROSSOVER_RULES"]["CROSSOVER_ANGLE"]);
        std::cout << "Crossover rule settings parsed...\n";

        // Zippering rule settings
        ENABLE_ZIPPERING = bool(interface["GRAMMAR"]["ZIPPERING_RULES"]["ENABLE_ZIPPERING"]);
        ZIPPERING_HIT_RATE = double(interface["GRAMMAR"]["ZIPPERING_RULES"]["ZIPPERING_HIT_RATE"]);
        ZIPPERING_GUARD_RATE = double(interface["GRAMMAR"]["ZIPPERING_RULES"]["ZIPPERING_GUARD_RATE"]);
        ZIPPERING_RETRACTION_RATE = double(interface["GRAMMAR"]["ZIPPERING_RULES"]["ZIPPERING_RETRACTION_RATE"]);
        CRITICAL_ANGLE = double(interface["GRAMMAR"]["ZIPPERING_RULES"]["CRITICAL_ANGLE"]);
        SEPARATION_DISTANCE = double(interface["GRAMMAR"]["ZIPPERING_RULES"]["SEPARATION_DISTANCE"]);
        std::cout << "Zippering rule settings parsed...\n";

        // Clasp rule settings
        CLASP_ENABLE_ENTRY = bool(interface["GRAMMAR"]["CLASP_RULES"]["CLASP_ENABLE_ENTRY"]);
        CLASP_ENTRY_RATE = double(interface["GRAMMAR"]["CLASP_RULES"]["CLASP_ENTRY_RATE"]);
        CLASP_ENTRY_ANGLE = double(interface["GRAMMAR"]["CLASP_RULES"]["CLASP_ENTRY_ANGLE"]);
        CLASP_ENABLE_EXIT = bool(interface["GRAMMAR"]["CLASP_RULES"]["CLASP_ENABLE_EXIT"]);
        CLASP_EXIT_RATE = double(interface["GRAMMAR"]["CLASP_RULES"]["CLASP_EXIT_RATE"]);
        CLASP_EXIT_ANGLE = double(interface["GRAMMAR"]["CLASP_RULES"]["CLASP_EXIT_ANGLE"]);
        CLASP_ENABLE_CAT = bool(interface["GRAMMAR"]["CLASP_RULES"]["CLASP_ENABLE_CAT"]);
        CLASP_CAT_RATE = double(interface["GRAMMAR"]["CLASP_RULES"]["CLASP_CAT_RATE"]);
        std::cout << "Clasp rule settings parsed...\n";

        // Destruction rule settings
        ENABLE_MT_DESTRUCTION = bool(interface["GRAMMAR"]["DESTRUCTION_RULES"]["ENABLE_MT_DESTRUCTION"]);
        MT_DESTRUCTION_RATE = double(interface["GRAMMAR"]["DESTRUCTION_RULES"]["MT_DESTRUCTION_RATE"]);
        std::cout << "Destruction rule settings parsed...\n";

        // Creation rule settings
        ENABLE_CREATION = bool(interface["GRAMMAR"]["CREATION_RULES"]["ENABLE_CREATION"]);
        CREATION_RATE = double(interface["GRAMMAR"]["CREATION_RULES"]["CREATION_RATE"]);
        CREATION_FACTOR = double(interface["GRAMMAR"]["CREATION_RULES"]["CREATION_FACTOR"]);
        std::cout << "Creation rule settings parsed...\n";

        // Recovery rule settings
        ENABLE_RECOVERY = bool(interface["GRAMMAR"]["RECOVERY_RULES"]["ENABLE_RECOVERY"]);
        RECOVERY_RATE = double(interface["GRAMMAR"]["RECOVERY_RULES"]["RECOVERY_RATE"]);
        RECOVERY_FACTOR = double(interface["GRAMMAR"]["RECOVERY_RULES"]["RECOVERY_FACTOR"]);
        std::cout << "Recovery rule settings parsed...\n";

        std::cout << "Grammar settings parsed...\n";

        SIGMOID_K = double(interface["SIGMOID_K"]);
        RHO_TEST_RATE = double(interface["RHO_TEST_RATE"]);
        COLLISION_DISTANCE = double(interface["COLLISION_DISTANCE"]);
        std::cout << "Other settings parsed...\n";
    }
};

#endif //DGGML_PARAMETERS_HPP
