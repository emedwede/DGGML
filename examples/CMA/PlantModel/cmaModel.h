#ifndef DGGML_CMAMODEL_H
#define DGGML_CMAMODEL_H

#include "DGGML.h"
#include "PlantGrammar.hpp"
#include "simdjson.h"
#include "ExpandedComplex2D.hpp"
#include "YAGL_Algorithms.hpp"

namespace CMA {

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
        int NUM_INTERNAL_STEPS;
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
    };

    using graph_grammar_t = DGGML::Grammar<Plant::graph_type>;
    class cmaModel : public DGGML::Model<graph_grammar_t> {
    public:
        void initialize() override
        {
            std::cout << "Initializing the plant model simulation\n";

            settings = {
                    "treadmilling",
                        10.0,
                        2.0,
                        10,
                        2,
                        2,
                        8.0,
                        8.0,
                        true,
                        64,
                        0.5,
                        0.8,
                        10,
                        1.2,
                        1.0,
                        0.125,
                        1.0, 0.25,
                        10.0,
                        5.0,
                        1.0,
                        1.0,
                        10.0
            };

            name = settings.EXPERIMENT_NAME;
            std::cout << "Creating the grammar\n";
            using GT = Plant::graph_type;

            //stochastic growing rule
            GT g1;
            g1.addNode({1, {Plant::Intermediate{}}});
            g1.addNode({2, {Plant::Positive{}}});
            g1.addEdge(1, 2);

            GT g2;
            g2.addNode({1, {Plant::Intermediate{}}});
            g2.addNode({3, {Plant::Intermediate{}}});
            g2.addNode({2, {Plant::Positive{}}});
            g2.addEdge(1, 3);
            g2.addEdge(3, 2);

            DGGML::WithRule<GT> r1;
            r1.name("growing").lhs(g1).rhs(g2)
                .with([](auto& lhs, auto& m) { return 0.5; })
                .where([](auto& lhs, auto& rhs, auto& m) { std::cout << "updating growing rule\n"; });

            gamma.addRule(r1);

            //stochastic retraction rule
            GT g3;
            g3.addNode({1, {Plant::Negative{}}});
            g3.addNode({2, {Plant::Intermediate{}}});
            g3.addNode({3, {Plant::Intermediate{}}});
            g3.addEdge(1, 2);
            g3.addEdge(2, 3);

            GT g4;
            g4.addNode({1, {Plant::Negative{}}});
            g4.addNode({3, {Plant::Intermediate{}}});
            g4.addEdge(1, 3);

            DGGML::WithRule<GT> r2;
            r2.name("retraction").lhs(g3).rhs(g4)
                    .with([](auto& lhs, auto& m) { std::cout << "retraction propensity\n"; return 0.5; })
                    .where([](auto& lhs, auto& rhs, auto& m) { std::cout << "updating retraction rule\n"; });

            gamma.addRule(r2);

            //collision catastrophe rule
            GT g5;
            g5.addNode({1, {Plant::Intermediate{}}});
            g5.addNode({2, {Plant::Positive{}}});
            g5.addNode({3, {Plant::Intermediate{}}});
            g5.addNode({4, {Plant::Intermediate{}}});
            g5.addNode({5, {Plant::Intermediate{}}});
            g5.addEdge(1, 2);
            g5.addEdge(3, 4);
            g5.addEdge(4, 5);

            GT g6;
            g6.addNode({1, {Plant::Intermediate{}}});
            g6.addNode({2, {Plant::Negative{}}});
            g6.addNode({3, {Plant::Intermediate{}}});
            g6.addNode({4, {Plant::Intermediate{}}});
            g6.addNode({5, {Plant::Intermediate{}}});
            g6.addEdge(1, 2);
            g6.addEdge(3, 4);
            g6.addEdge(4, 5);


            DGGML::WithRule<GT> r3;
            r3.name("catastrophe").lhs(g5).rhs(g6)
                    .with([](auto& lhs, auto& m) { std::cout << "catastrophe propensity\n"; return 7.5; })
                    .where([](auto& lhs, auto& rhs, auto& m) { std::cout << "updating catastrophe rule\n"; });

            gamma.addRule(r3);

            //Testable Interaction Rule
            //collision catastrophe rule
            GT g7;
            g7.addNode({1, {Plant::Negative{}}});
            g7.addNode({2, {Plant::Intermediate{}}});
            g7.addNode({3, {Plant::Positive{}}});
            g7.addNode({4, {Plant::Negative{}}});
            g7.addNode({5, {Plant::Intermediate{}}});
            g7.addNode({6, {Plant::Positive{}}});
            g7.addEdge(1, 2);
            g7.addEdge(2, 3);
            g7.addEdge(4, 5);
            g7.addEdge(5, 6);

            GT g8;

            DGGML::WithRule<GT> r4;
            r4.name("interaction").lhs(g7).rhs(g8)
                    .with([](auto& lhs, auto& m) { std::cout << "interaction propensity\n"; return 7.5; })
                    .where([](auto& lhs, auto& rhs, auto& m) { std::cout << "updating interaction rule\n"; });

            gamma.addRule(r4);

            std::cout << "Generating the expanded cell complex\n";
            geoplex2D.init(settings.CELL_NX,
                           settings.CELL_NY,
                           settings.CELL_DX,
                           settings.CELL_DY,
                           settings.GHOSTED,
                           settings.MAXIMAL_REACTION_RADIUS); //ghosted
            //std::cout << geoplex2D;

            std::cout << "Initializing the system graph\n";
            Plant::microtubule_uniform_scatter(system_graph, geoplex2D, settings);
        }

        struct MyMetrics {
            std::vector<std::size_t> con_com;
            std::vector<std::size_t> total_nodes;
            std::vector<std::size_t> type_counts[5];
            std::vector<double> time_count;

            std::size_t junction_count;
            std::size_t positive_count;
            std::size_t negative_count;
            std::size_t zipper_count;
            std::size_t intermediate_count;

            void reset_count()
            {
                junction_count = 0;
                positive_count = 0;
                negative_count = 0;
                zipper_count = 0;
                intermediate_count = 0;
            }
        } metrics;

        void collect() override
        {
            metrics.reset_count();
            metrics.con_com.push_back(YAGL::connected_components(system_graph));

            for(auto iter = system_graph.node_list_begin();
                iter != system_graph.node_list_end(); iter++) {
                auto itype = iter->second.getData();
                if(std::holds_alternative<Plant::Negative>(itype.data))
                    metrics.negative_count++;
                if(std::holds_alternative<Plant::Positive>(itype.data))
                    metrics.positive_count++;
                if(std::holds_alternative<Plant::Intermediate>(itype.data))
                    metrics.intermediate_count++;
                if(std::holds_alternative<Plant::Junction>(itype.data))
                    metrics.junction_count++;
                if(std::holds_alternative<Plant::Zipper>(itype.data))
                    metrics.zipper_count++;
            }
            metrics.type_counts[0].push_back(metrics.negative_count);
            metrics.type_counts[1].push_back(metrics.positive_count);
            metrics.type_counts[2].push_back(metrics.intermediate_count);
            metrics.type_counts[3].push_back(metrics.junction_count);
            metrics.type_counts[4].push_back(metrics.zipper_count);

            metrics.total_nodes.push_back(metrics.negative_count + metrics.positive_count +
            metrics.intermediate_count + metrics.junction_count + metrics.zipper_count);
        }

        void print_metrics() override
        {
            print_numpy_array_stats(metrics.con_com, "con_com");
            print_numpy_array_stats(metrics.type_counts[0], "negative");
            print_numpy_array_stats(metrics.type_counts[1], "positive");
            print_numpy_array_stats(metrics.type_counts[2], "intermediate");
            print_numpy_array_stats(metrics.type_counts[3], "junction");
            print_numpy_array_stats(metrics.type_counts[4], "zipper");
            print_numpy_array_stats(metrics.total_nodes, "total_nodes");
            print_numpy_array_stats(metrics.time_count, "time_count");
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

            settings.NUM_MT = int64_t(interface["SETTINGS"]["NUM_MT"]);
            settings.MT_MIN_SEGMENT_INIT = double(interface["SETTINGS"]["MT_MIN_SEGMENT_INIT"]);
            settings.MT_MAX_SEGMENT_INIT = double(interface["SETTINGS"]["MT_MAX_SEGMENT_INIT"]);

            settings.LENGTH_DIV_FACTOR = double(interface["SETTINGS"]["LENGTH_DIV_FACTOR"]);
            settings.DIV_LENGTH = double(interface["SETTINGS"]["DIV_LENGTH"]);
            settings.DIV_LENGTH_RETRACT = double(interface["SETTINGS"]["DIV_LENGTH_RETRACT"]);
            settings.V_PLUS = double(interface["SETTINGS"]["V_PLUS"]);
            settings.V_MINUS = double(interface["SETTINGS"]["V_MINUS"]);

            settings.SIGMOID_K = double(interface["SETTINGS"]["SIGMOID_K"]);

            settings.MAXIMAL_REACTION_RADIUS = settings.DIV_LENGTH*2.0;

            //Simulate until the specified unit time
            settings.TOTAL_TIME = double(interface["SETTINGS"]["TOTAL_TIME"]);
            settings.NUM_INTERNAL_STEPS = 5;
            //Delta should be big, but not to big. In this case, the maximum amount of time it would
            //take one MT to grow a single unit of MT
            settings.DELTA =
                    0.25*settings.MAXIMAL_REACTION_RADIUS / std::max(settings.V_PLUS, settings.V_MINUS);
            //The internal step of the solver should be at least this small
            settings.DELTA_DELTA_T = settings.DELTA / settings.NUM_INTERNAL_STEPS;
            settings.DELTA_T_MIN = settings.DELTA_DELTA_T;
            settings.NUM_STEPS = settings.TOTAL_TIME / settings.DELTA;

            settings.RHO_TEST_RATE = double(interface["EXPERIMENTAL"]["RHO_TEST_RATE"]);
        }
        Parameters settings;
    };
}


#endif //DGGML_CMAMODEL_H
