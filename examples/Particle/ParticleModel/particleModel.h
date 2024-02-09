#ifndef DGGML_PARTICLE_MODEL_H
#define DGGML_PARTICLE_MODEL_H

#include "DGGML.h"
#include "ParticleUtils.hpp"
#include "MathUtils.hpp"
#include "simdjson.h"
#include "ExpandedComplex2D.hpp"
#include "YAGL_Algorithms.hpp"

//molecular dynamics (MD)
namespace MD {

    template <typename DataType>
    void print_numpy_array_stats(DataType& data, std::string var_name)
    {
        std::cout << var_name << " = np.asarray([";
        for(auto i = 0; i < data.size(); i++)
        {
            if(i != 0 && i % 20 == 0) std::cout << "\n";
            if(i != data.size() - 1)
                std::cout << data[i] << ", ";
            else
                std::cout << data[i];
        }
        std::cout << "]);\n";
    }

    //TODO: user and system params should be split
    // This also may be more appropriate for the base class
    struct Parameters
    {
        std::string EXPERIMENT_NAME;
        double DELTA;
        double DELTA_DELTA_T;
        std::size_t CELL_NX;
        std::size_t CELL_NY;
        double CELL_DX;
        double CELL_DY;
        bool GHOSTED;
        std::size_t NUM_PARTICLES;
        std::size_t NUM_STEPS;
        double V_PLUS;
        double TOTAL_TIME;
        double MAXIMAL_REACTION_RADIUS;
        double DELTA_T_MIN;
    };

    //checks if we're outside a box boundary
    bool boundary_check_2D(Parameters& settings, double x, double y, double padding = 0.0)
    {
        auto min_x = 0.0+padding, min_y = 0.0+padding;
        auto max_x = settings.CELL_NX*settings.CELL_DX-padding;
        auto max_y = settings.CELL_NY*settings.CELL_DY-padding;

        return (x >= min_x && x <= max_x && y >= min_y && y <= max_y) ? false : true;
    }

    using graph_grammar_t = DGGML::Grammar<MD::graph_type>;
    class particleModel : public DGGML::Model<graph_grammar_t> {
    public:
        void initialize() override
        {
            std::cout << "Initializing the particle simulation\n";

            set_parameters();
            name = settings.EXPERIMENT_NAME;
            std::cout << "Creating the grammar\n";
            using GT = MD::graph_type;

            GT boundary_bounce_lhs;
            boundary_bounce_lhs.addNode({1, {MD::Molecule{}}});

            GT boundary_bounce_rhs;
            boundary_bounce_rhs.addNode({1, {MD::Molecule{}}});


            DGGML::WithRule<GT> boundary_bounce("boundary_bounce", boundary_bounce_lhs, boundary_bounce_rhs,
         [&](auto& lhs, auto& m)
            {
                //determine if we're outside the interior of a padded boundary
                if(boundary_check_2D(settings, lhs[m[1]].position[0], lhs[m[1]].position[1], settings.MAXIMAL_REACTION_RADIUS))
                    return 1000.0;
                return 0.0;
            },
            [&](auto& lhs, auto& rhs, auto& m1, auto& m2)
            {
                //ensures particles that move out of bounds are put back on the boundary
                auto padding = settings.MAXIMAL_REACTION_RADIUS;
                auto min_x = 0.0+padding, min_y = 0.0+padding;
                auto max_x = settings.CELL_NX*settings.CELL_DX-padding;
                auto max_y = settings.CELL_NY*settings.CELL_DY-padding;
                rhs[m2[1]].position[0] = (rhs[m2[1]].position[0] < min_x) ? min_x : rhs[m2[1]].position[0];
                rhs[m2[1]].position[0] = (rhs[m2[1]].position[0] > max_x) ? max_x : rhs[m2[1]].position[0];
                rhs[m2[1]].position[1] = (rhs[m2[1]].position[1] < min_y) ? min_y: rhs[m2[1]].position[1];
                rhs[m2[1]].position[1] = (rhs[m2[1]].position[1] > max_y) ? max_y : rhs[m2[1]].position[1];
                std::get<MD::Molecule>(rhs[m2[1]].data).unit_vec[0] = -std::get<MD::Molecule>(lhs[m1[1]].data).unit_vec[0];
                std::get<MD::Molecule>(rhs[m2[1]].data).unit_vec[1] = -std::get<MD::Molecule>(lhs[m1[1]].data).unit_vec[1];
                std::get<MD::Molecule>(rhs[m2[1]].data).unit_vec[2] = -std::get<MD::Molecule>(lhs[m1[1]].data).unit_vec[2];
            });
            gamma.addRule(boundary_bounce);



            GT destroy_molecule_lhs;
            destroy_molecule_lhs.addNode({1, {MD::Molecule{}}});

            GT destroy_molecule_rhs;

            DGGML::WithRule<GT> destroy_molecule("destroy_molecule", destroy_molecule_lhs, destroy_molecule_rhs,
                    [](auto& lhs, auto& m) { return 0.2; },
                    [](auto& lhs, auto& rhs, auto& m1, auto& m2) {});

            gamma.addRule(destroy_molecule);

            GT particle_motion_lhs;
            particle_motion_lhs.addNode({1, {MD::Molecule{}}});
            //TODO: I think I need to add velocity back in, and make the growing solve a function of the two ODEs
            DGGML::SolvingRule<GT> particle_motion("solving_particle_motion", particle_motion_lhs, particle_motion_lhs, 3,
                    [](auto& lhs, auto& m1, auto& varset)
                    {
                        //TODO: fix and account for proper verlet integration
                        //bind the variables involved
                        varset.insert(&lhs[m1[1]].position[0]);
                        varset.insert(&lhs[m1[1]].position[1]);
                        varset.insert(&lhs[m1[1]].position[2]);
                    },
                    [&](auto& lhs, auto& m1, auto y, auto ydot, auto& varmap) {
                        auto v_plus = settings.V_PLUS;

                        auto& data1 = std::get<MD::Molecule>(lhs[m1[1]].data);
                        for(auto i = 0; i < 3; i++) {
                            if (auto search = varmap.find(&data1.unit_vec[i]); search != varmap.end())
                                NV_Ith_S(ydot, varmap[&lhs[m1[1]].position[i]]) += v_plus * data1.unit_vec[i]*NV_Ith_S(y, search->second);
                            else
                                NV_Ith_S(ydot, varmap[&lhs[m1[1]].position[i]]) += v_plus * data1.unit_vec[i];
                        }
                    });

            gamma.addRule(particle_motion);

            std::cout << "Generating the expanded cell complex\n";
            geoplex2D.init(settings.CELL_NX,
                           settings.CELL_NY,
                           settings.CELL_DX,
                           settings.CELL_DY,
                           settings.GHOSTED,
                           settings.MAXIMAL_REACTION_RADIUS); //ghosted
            //std::cout << geoplex2D;

            std::cout << "Initializing the system graph\n";
            MD::particle_uniform_scatter(system_graph, geoplex2D, settings, gen);
        }

        struct MyMetrics {
            std::vector<std::size_t> total_nodes;
            std::vector<std::size_t> type_counts;
            std::vector<double> time_count;

            std::size_t molecule_count;

            void reset_count()
            {
                molecule_count = 0;
            }
        } metrics;

        void collect() override
        {
            metrics.reset_count();
            //metrics.con_com.push_back(YAGL::connected_components(system_graph));

            for(auto iter = system_graph.node_list_begin();
                iter != system_graph.node_list_end(); iter++) {
                auto itype = iter->second.getData();
                if(std::holds_alternative<MD::Molecule>(itype.data))
                    metrics.molecule_count++;
            }
            metrics.type_counts.push_back(metrics.molecule_count);

            metrics.total_nodes.push_back(metrics.molecule_count);
        }

        void print_metrics() override
        {
            print_numpy_array_stats(metrics.type_counts, "molecule");
            print_numpy_array_stats(metrics.total_nodes, "total_nodes");
            print_numpy_array_stats(metrics.time_count, "time_count");
        }

        void set_parameters() {
            settings.EXPERIMENT_NAME = "box";

            //1x1 micrometer domain
            settings.CELL_NX = 2;//1;
            settings.CELL_NY = 2;//1;
            settings.CELL_DX = 1.5;//;0.5;//1.0;
            settings.CELL_DY = 1.5;//0.5;//1.0;

            //non ghosted complex
            settings.GHOSTED = false;

            //number of microtubules in the simulation
            settings.NUM_PARTICLES = 1000;

            settings.V_PLUS = 0.0615;

            settings.MAXIMAL_REACTION_RADIUS = 2*0.05;

            //simulation time in seconds
            settings.TOTAL_TIME = 25.0;//25.0;//20.0;
            settings.DELTA = 0.5/8.0; //unit of seconds

            //The internal step of the solver should be at least smaller than delta
            settings.DELTA_DELTA_T = settings.DELTA / 20.0;
            settings.DELTA_T_MIN = settings.DELTA_DELTA_T;
            settings.NUM_STEPS = settings.TOTAL_TIME / settings.DELTA;
        }

        void set_parameters(simdjson::ondemand::document& interface)
        {
            std::string_view temp = interface["META"]["EXPERIMENT"];
            settings.EXPERIMENT_NAME = static_cast<std::string>(temp);

            name = settings.EXPERIMENT_NAME;

            std::cout << settings.EXPERIMENT_NAME << "+++\n";
            settings.CELL_NX = int64_t(interface["SETTINGS"]["CELL_NX"]);
            settings.CELL_NY = int64_t(interface["SETTINGS"]["CELL_NY"]);

            settings.CELL_DX = double(interface["SETTINGS"]["CELL_DX"]);
            settings.CELL_DY = double(interface["SETTINGS"]["CELL_DY"]);
            settings.GHOSTED = bool(interface["SETTINGS"]["GHOSTED"]);

            settings.NUM_PARTICLES = int64_t(interface["SETTINGS"]["NUM_PARTICLES"]);

            settings.V_PLUS = double(interface["SETTINGS"]["V_PLUS"]);

            settings.MAXIMAL_REACTION_RADIUS = double(interface["SETTINGS"]["REACTION_RADIUS"]);

            //Simulate until the specified unit time
            settings.TOTAL_TIME = double(interface["SETTINGS"]["TOTAL_TIME"]);
            settings.DELTA = double(interface["SETTINGS"]["DELTA"]);
            //The internal step of the solver should be at least smaller than delta
            settings.DELTA_DELTA_T = settings.DELTA / 20.0;// / settings.NUM_INTERNAL_STEPS;
            settings.DELTA_T_MIN = settings.DELTA_DELTA_T;
            settings.NUM_STEPS = settings.TOTAL_TIME / settings.DELTA;
        }
        Parameters settings;
    };
}


#endif //DGGML_CMAMODEL_H
