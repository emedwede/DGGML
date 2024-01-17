#ifndef DGGML_CMAMODEL_H
#define DGGML_CMAMODEL_H

#include "DGGML.h"
#include "PlantGrammar.hpp"
#include "PlantUtils.hpp"
#include "MathUtils.hpp"
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

    bool boundary_check_2D(Parameters& settings, double x, double y)
    {
        auto min_x = 0.0, min_y = 0.0;
        auto max_x = settings.CELL_NX*settings.CELL_DX;
        auto max_y = settings.CELL_NY*settings.CELL_DY;

        return (x > min_x && x < max_x && y > min_y && y < max_y) ? false : true;
    }

    using graph_grammar_t = DGGML::Grammar<Plant::graph_type>;
    class cmaModel : public DGGML::Model<graph_grammar_t> {
    public:
        void initialize() override
        {
            std::cout << "Initializing the plant model simulation\n";

            set_parameters();
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


            //TODO: I should make it so that any solving/propensity functions that need access to parameters
            // are actually passed as functors with states!
            DGGML::WithRule<GT> r1("growing", g1, g2,
                                   [&](auto& lhs, auto& m)
                {
                    //return 0.0*2.0;
                    auto& node_i_data = lhs.findNode(m[1])->second.getData();
                    auto& node_j_data = lhs.findNode(m[2])->second.getData();
                    auto len = DGGML::calculate_distance(node_i_data.position, node_j_data.position);
                    double propensity = 10.0*DGGML::heaviside(len, settings.DIV_LENGTH);
                    //double propensity = DGGML::sigmoid((len/settings.DIV_LENGTH) - 1.0, settings.SIGMOID_K);
                    return propensity;
                }, [](auto& lhs, auto& rhs, auto& m1, auto& m2) {
                        //in this function we are responsible for setting all new parameters
                        //first set the position
                        //modified this rule s.t. new node is placed just behind the growing end not halfway
//                        rhs[m2[3]].position[0] = (lhs[m1[2]].position[0] + lhs[m1[1]].position[0])/2.0;
//                        rhs[m2[3]].position[1] = (lhs[m1[2]].position[1] + lhs[m1[1]].position[1])/2.0;
//                        rhs[m2[3]].position[2] = (lhs[m1[2]].position[2] + lhs[m1[1]].position[2])/2.0;
                        rhs[m2[3]].position[0] = lhs[m1[2]].position[0] - (lhs[m1[2]].position[0] - lhs[m1[1]].position[0])/100.0;
                        rhs[m2[3]].position[1] = lhs[m1[2]].position[1] - (lhs[m1[2]].position[1] - lhs[m1[1]].position[1])/100.0;
                        rhs[m2[3]].position[2] = lhs[m1[2]].position[2] - (lhs[m1[2]].position[2] - lhs[m1[1]].position[2])/100.0;
                        //next set the unit vector
                        std::get<Plant::Intermediate>(rhs[m2[3]].data).unit_vec[0] = std::get<Plant::Intermediate>(lhs[m1[1]].data).unit_vec[0];
                        std::get<Plant::Intermediate>(rhs[m2[3]].data).unit_vec[1] = std::get<Plant::Intermediate>(lhs[m1[1]].data).unit_vec[1];
                        std::get<Plant::Intermediate>(rhs[m2[3]].data).unit_vec[2] = std::get<Plant::Intermediate>(lhs[m1[1]].data).unit_vec[2];
                    });

            gamma.addRule(r1);

            //TODO: reduce the copy paste coding
            //collision catastrophe rule case 1
            GT catastrophe_lhs_graph1;
            catastrophe_lhs_graph1.addNode({1, {Plant::Intermediate{}}});
            catastrophe_lhs_graph1.addNode({2, {Plant::Positive{}}});
            catastrophe_lhs_graph1.addNode({3, {Plant::Intermediate{}}});
            catastrophe_lhs_graph1.addNode({4, {Plant::Intermediate{}}});
            catastrophe_lhs_graph1.addEdge(1, 2);
            catastrophe_lhs_graph1.addEdge(3, 4);

            GT catastrophe_rhs_graph1;
            catastrophe_rhs_graph1.addNode({1, {Plant::Intermediate{}}});
            catastrophe_rhs_graph1.addNode({2, {Plant::Negative{}}});
            catastrophe_rhs_graph1.addNode({3, {Plant::Intermediate{}}});
            catastrophe_rhs_graph1.addNode({4, {Plant::Intermediate{}}});
            catastrophe_rhs_graph1.addEdge(1, 2);
            catastrophe_rhs_graph1.addEdge(3, 4);

            DGGML::WithRule<GT> catastrophe_case1("catastrophe_case1", catastrophe_lhs_graph1, catastrophe_rhs_graph1,
                    [&](auto& lhs, auto& m) {
                        //find all the node data
                        auto& dat1 = lhs.findNode(m[1])->second.getData();
                        auto& dat2 = lhs.findNode(m[2])->second.getData();
                        auto& dat3 = lhs.findNode(m[3])->second.getData();
                        auto& dat4 = lhs.findNode(m[4])->second.getData();

                        //get references to position vector
                        auto& pos1 = dat1.position;
                        auto& pos2 = dat2.position;
                        auto& pos3 = dat3.position;
                        auto& pos4 = dat4.position;

                        //get references to unit vectors
                        auto& u1 = std::get<Plant::Intermediate>(dat1.data).unit_vec;
                        auto& u2 = std::get<Plant::Positive>(dat2.data).unit_vec;
                        auto& u3 = std::get<Plant::Intermediate>(dat3.data).unit_vec;
                        auto& u4 = std::get<Plant::Intermediate>(dat4.data).unit_vec;

                        double theta = DGGML::unit_dot_product(u2, u3);
                        double propensity = 0.0;
                        double sol[2];
                        if(theta != 1.0 && theta != -1.0)
                        {
                            DGGML::paramaterized_intersection(pos2, pos4, pos3, u2, sol);

                            if(sol[0] > 0.0 && sol[1] >= 0.0 && sol[1] <= 1.0)
                            {
                                //TODO: distance should be closest point to the growing end not pos3
                                propensity = 1000.0*exp(-pow(DGGML::calculate_distance(pos2, pos3), 2.0) / pow(0.5*settings.DIV_LENGTH, 2.0));
                            }
                        }
                       //if(propensity != 0) {std::cout << "cat prop: " << propensity << "\n"; std::cin.get();}
                        return propensity;
                    },
                    [](auto& lhs, auto& rhs, auto& m1, auto& m2)
                    {
                        std::get<Plant::Negative>(rhs[m2[2]].data).unit_vec[0] = -std::get<Plant::Positive>(lhs[m1[2]].data).unit_vec[0];
                        std::get<Plant::Negative>(rhs[m2[2]].data).unit_vec[1] = -std::get<Plant::Positive>(lhs[m1[2]].data).unit_vec[1];
                        std::get<Plant::Negative>(rhs[m2[2]].data).unit_vec[2] = -std::get<Plant::Positive>(lhs[m1[2]].data).unit_vec[2];
                    });

            gamma.addRule(catastrophe_case1);

            //collision catastrophe rule case 2
            GT catastrophe_lhs_graph2;
            catastrophe_lhs_graph2.addNode({1, {Plant::Intermediate{}}});
            catastrophe_lhs_graph2.addNode({2, {Plant::Positive{}}});
            catastrophe_lhs_graph2.addNode({3, {Plant::Intermediate{}}});
            catastrophe_lhs_graph2.addNode({4, {Plant::Positive{}}});
            catastrophe_lhs_graph2.addEdge(1, 2);
            catastrophe_lhs_graph2.addEdge(3, 4);

            GT catastrophe_rhs_graph2;
            catastrophe_rhs_graph2.addNode({1, {Plant::Intermediate{}}});
            catastrophe_rhs_graph2.addNode({2, {Plant::Negative{}}});
            catastrophe_rhs_graph2.addNode({3, {Plant::Intermediate{}}});
            catastrophe_rhs_graph2.addNode({4, {Plant::Positive{}}});
            catastrophe_rhs_graph2.addEdge(1, 2);
            catastrophe_rhs_graph2.addEdge(3, 4);

            DGGML::WithRule<GT> catastrophe_case2("catastrophe_case2", catastrophe_lhs_graph2, catastrophe_rhs_graph2,
                                                  [&](auto& lhs, auto& m) {
                                                      //find all the node data
                                                      auto& dat1 = lhs.findNode(m[1])->second.getData();
                                                      auto& dat2 = lhs.findNode(m[2])->second.getData();
                                                      auto& dat3 = lhs.findNode(m[3])->second.getData();
                                                      auto& dat4 = lhs.findNode(m[4])->second.getData();

                                                      //get references to position vector
                                                      auto& pos1 = dat1.position;
                                                      auto& pos2 = dat2.position;
                                                      auto& pos3 = dat3.position;
                                                      auto& pos4 = dat4.position;

                                                      //get references to unit vectors
                                                      auto& u1 = std::get<Plant::Intermediate>(dat1.data).unit_vec;
                                                      auto& u2 = std::get<Plant::Positive>(dat2.data).unit_vec;
                                                      auto& u3 = std::get<Plant::Intermediate>(dat3.data).unit_vec;
                                                      auto& u4 = std::get<Plant::Positive>(dat4.data).unit_vec;

                                                      double theta = DGGML::unit_dot_product(u2, u3);
                                                      double propensity = 0.0;
                                                      double sol[2];
                                                      if(theta != 1.0 && theta != -1.0)
                                                      {
                                                          DGGML::paramaterized_intersection(pos2, pos4, pos3, u2, sol);

                                                          if(sol[0] > 0.0 && sol[1] >= 0.0 && sol[1] <= 1.0)
                                                          {
                                                              //TODO: distance should be closest point to the growing end not pos3
                                                              propensity = 1000.0*exp(-pow(DGGML::calculate_distance(pos2, pos3), 2.0) / pow(0.5*settings.DIV_LENGTH, 2.0));
                                                          }
                                                      }
                                                      //if(propensity != 0) {std::cout << "cat prop: " << propensity << "\n"; std::cin.get();}
                                                      return propensity;
                                                  },
                                                  [](auto& lhs, auto& rhs, auto& m1, auto& m2)
                                                  {
                                                      std::get<Plant::Negative>(rhs[m2[2]].data).unit_vec[0] = -std::get<Plant::Positive>(lhs[m1[2]].data).unit_vec[0];
                                                      std::get<Plant::Negative>(rhs[m2[2]].data).unit_vec[1] = -std::get<Plant::Positive>(lhs[m1[2]].data).unit_vec[1];
                                                      std::get<Plant::Negative>(rhs[m2[2]].data).unit_vec[2] = -std::get<Plant::Positive>(lhs[m1[2]].data).unit_vec[2];
                                                  });

            gamma.addRule(catastrophe_case2);

            //collision catastrophe rule case 3
            GT catastrophe_lhs_graph3;
            catastrophe_lhs_graph3.addNode({1, {Plant::Intermediate{}}});
            catastrophe_lhs_graph3.addNode({2, {Plant::Positive{}}});
            catastrophe_lhs_graph3.addNode({3, {Plant::Intermediate{}}});
            catastrophe_lhs_graph3.addNode({4, {Plant::Negative{}}});
            catastrophe_lhs_graph3.addEdge(1, 2);
            catastrophe_lhs_graph3.addEdge(3, 4);

            //TODO: FIX ME! The rhs of this rule is messing with the catastrophe solving rule
            GT catastrophe_rhs_graph3;
            catastrophe_rhs_graph3.addNode({1, {Plant::Intermediate{}}});
            catastrophe_rhs_graph3.addNode({2, {Plant::Negative{}}});
            catastrophe_rhs_graph3.addNode({3, {Plant::Intermediate{}}});
            catastrophe_rhs_graph3.addNode({4, {Plant::Negative{}}});
            catastrophe_rhs_graph3.addEdge(1, 2);
            catastrophe_rhs_graph3.addEdge(3, 4);

            DGGML::WithRule<GT> catastrophe_case3("catastrophe_case3", catastrophe_lhs_graph3, catastrophe_rhs_graph3,
                                                  [&](auto& lhs, auto& m) {

                                                      //find all the node data
                                                      auto& dat1 = lhs.findNode(m[1])->second.getData();
                                                      auto& dat2 = lhs.findNode(m[2])->second.getData();
                                                      auto& dat3 = lhs.findNode(m[3])->second.getData();
                                                      auto& dat4 = lhs.findNode(m[4])->second.getData();

                                                      //get references to position vector
                                                      auto& pos1 = dat1.position;
                                                      auto& pos2 = dat2.position;
                                                      auto& pos3 = dat3.position;
                                                      auto& pos4 = dat4.position;

                                                      //get references to unit vectors
                                                      auto& u1 = std::get<Plant::Intermediate>(dat1.data).unit_vec;
                                                      auto& u2 = std::get<Plant::Positive>(dat2.data).unit_vec;
                                                      auto& u3 = std::get<Plant::Intermediate>(dat3.data).unit_vec;
                                                      auto& u4 = std::get<Plant::Negative>(dat4.data).unit_vec;

                                                      double theta = DGGML::unit_dot_product(u2, u3);
                                                      double propensity = 0.0;
                                                      double sol[2];
                                                      if(theta != 1.0 && theta != -1.0)
                                                      {
                                                          DGGML::paramaterized_intersection(pos2, pos4, pos3, u2, sol);

                                                          if(sol[0] > 0.0 && sol[1] >= 0.0 && sol[1] <= 1.0)
                                                          {
                                                              //TODO: distance should be closest point to the growing end not pos3
                                                              propensity = 1000.0*exp(-pow(DGGML::calculate_distance(pos2, pos3), 2.0) / pow(0.5*settings.DIV_LENGTH, 2.0));
                                                          }
                                                      }

                                                      //if(propensity != 0) {std::cout << "cat prop: " << propensity << "\n"; std::cin.get();}
                                                      return propensity;
                                                  },
                                                  [](auto& lhs, auto& rhs, auto& m1, auto& m2)
                                                  {

                                                      std::get<Plant::Negative>(rhs[m2[2]].data).unit_vec[0] = -std::get<Plant::Positive>(lhs[m1[2]].data).unit_vec[0];
                                                      std::get<Plant::Negative>(rhs[m2[2]].data).unit_vec[1] = -std::get<Plant::Positive>(lhs[m1[2]].data).unit_vec[1];
                                                      std::get<Plant::Negative>(rhs[m2[2]].data).unit_vec[2] = -std::get<Plant::Positive>(lhs[m1[2]].data).unit_vec[2];

                                                  });

            gamma.addRule(catastrophe_case3);

            //destruction rule case 1: depolymerization on both ends
            GT destruction_lhs_graph1;
            destruction_lhs_graph1.addNode({1, {Plant::Negative{}}});
            destruction_lhs_graph1.addNode({2, {Plant::Intermediate{}}});
            destruction_lhs_graph1.addNode({3, {Plant::Negative{}}});
            destruction_lhs_graph1.addEdge(1, 2);
            destruction_lhs_graph1.addEdge(2, 3);

            GT destruction_rhs_graph1;

            DGGML::WithRule<GT> destruction_case1("destruction_case1", destruction_lhs_graph1, destruction_rhs_graph1,
                    [](auto& lhs, auto& m) { return 1; },
                    [](auto& lhs, auto& rhs, auto& m1, auto& m2) {});

            gamma.addRule(destruction_case1);

            //destruction rule case 2: depolymerization on one end polymerization on the other
            GT destruction_lhs_graph2;
            destruction_lhs_graph2.addNode({1, {Plant::Negative{}}});
            destruction_lhs_graph2.addNode({2, {Plant::Intermediate{}}});
            destruction_lhs_graph2.addNode({3, {Plant::Positive{}}});
            destruction_lhs_graph2.addEdge(1, 2);
            destruction_lhs_graph2.addEdge(2, 3);

            GT destruction_rhs_graph2;

            DGGML::WithRule<GT> destruction_case2("destruction_case2", destruction_lhs_graph2, destruction_rhs_graph2,
                                                  [](auto& lhs, auto& m) { return 1; },
                                                  [](auto& lhs, auto& rhs, auto& m1, auto& m2) {});

            //gamma.addRule(destruction_case2);


            //TODO: I think I need to add velocity back in, and make the growing solve a function of the two ODEs
            DGGML::SolvingRule<GT> r5("solving_grow", g1, g1, 3,
                    [](auto& lhs, auto& m1, auto& varset)
                    {
                        //TODO: fix and account for verlet integration
                        //bind the variables involved
                        varset.insert(&lhs[m1[2]].position[0]);
                        varset.insert(&lhs[m1[2]].position[1]);
                        varset.insert(&lhs[m1[2]].position[2]);
                    },
                    [&](auto& lhs, auto& m1, auto y, auto ydot, auto& varmap) {
                        //std::cout << "solving the grow rule\n";
                        //unless we know a variable wasn't used before, it's current value from it's solving
                        //ode must be checked and used i.e. if(varmap[&lhs[m[1]].position[0]]->second = false) do
                        //otherwise we don't have to check, but if the user is wrong, undefined behavior may ensue
                        //growth rule params
                        auto v_plus = settings.V_PLUS;
                        auto d_l = settings.DIV_LENGTH;
                        double l = 0.0;
                        for(auto i = 0; i < 3; i++)
                        {
                            //TODO: set a constraint that stops solving if a boundary is reached
                            double diff = NV_Ith_S(y, varmap[&lhs[m1[2]].position[i]]);
                            if(auto search = varmap.find(&lhs[m1[1]].position[i]); search != varmap.end())
                                diff -= NV_Ith_S(y, search->second);
                            else
                                diff -= lhs[m1[1]].position[i];
                            l += diff*diff;
                        }
                        l = sqrt(l);
                        double length_limiter = 1.0;//(1.0 - (l/d_l));
                        auto& data1 = std::get<Plant::Intermediate>(lhs[m1[1]].data);
                        for(auto i = 0; i < 3; i++) {
                            if (auto search = varmap.find(&data1.unit_vec[i]); search != varmap.end())
                                NV_Ith_S(ydot, varmap[&lhs[m1[2]].position[i]]) += v_plus * data1.unit_vec[i]*NV_Ith_S(y, search->second) * length_limiter;
                            else
                                NV_Ith_S(ydot, varmap[&lhs[m1[2]].position[i]]) += v_plus * data1.unit_vec[i] * length_limiter;
                        }

                        //TODO: see if we can make this internal to the algorithm and not user controlled
                        // so that deactivated ODEs can be removed from the system
                        //boundary check
                        //TODO: fix and account for verlet integration
                        auto x_plus_dx = NV_Ith_S(y, varmap[&lhs[m1[2]].position[0]]) + NV_Ith_S(ydot, varmap[&lhs[m1[2]].position[0]]);
                        auto y_plus_dy = NV_Ith_S(y, varmap[&lhs[m1[2]].position[1]]) + NV_Ith_S(ydot, varmap[&lhs[m1[2]].position[1]]);
                        bool out_of_bounds = boundary_check_2D(settings, x_plus_dx, y_plus_dy);
                        if(out_of_bounds) {
                            for (auto i = 0; i < 3; i++) {
                                NV_Ith_S(ydot, varmap[&lhs[m1[2]].position[i]]) = 0.0;
                            }
                        }
                    });

            gamma.addRule(r5);

            GT r6g1;
            r6g1.addNode({1, {Plant::Negative{}}});
            r6g1.addNode({2, {Plant::Intermediate{}}});
            r6g1.addEdge(1, 2);
            //TODO: I think I need to add velocity back in, and make the growing solve a function of the two ODEs
            DGGML::SolvingRule<GT> r6("solving_retraction", r6g1, r6g1, 3,
                                      [](auto& lhs, auto& m1, auto& varset)
                                      {
                                          //std::cout << "ic of the grow rule\n";
                                          //bind the variables involved
                                          varset.insert(&lhs[m1[1]].position[0]);
                                          varset.insert(&lhs[m1[1]].position[1]);
                                          varset.insert(&lhs[m1[1]].position[2]);
                                      },
                                      [&](auto& lhs, auto& m1, auto y, auto ydot, auto& varmap) {
                                          //std::cout << "here start\n";
                                          auto d_l_r = settings.DIV_LENGTH_RETRACT;
                                          auto v_minus = settings.V_MINUS;
                                          auto d_l = settings.DIV_LENGTH;
                                          double l = 0.0;
                                          for(auto i = 0; i < 3; i++)
                                          {
                                              double diff = NV_Ith_S(y, varmap[&lhs[m1[1]].position[i]]);
                                              if(auto search = varmap.find(&lhs[m1[2]].position[i]); search != varmap.end())
                                                  diff -= NV_Ith_S(y, search->second);
                                              else
                                                  diff -= lhs[m1[2]].position[i];
                                              l += diff*diff;
                                          }
                                          l = sqrt(l);
                                          double length_limiter = l/d_l;
                                          //if(length_limiter <= d_l_r) length_limiter = 0.0;
                                          if(l <= d_l_r/2.0) length_limiter = 0.0;
                                          else length_limiter = 1.0;
                                          auto& data1 = std::get<Plant::Negative>(lhs[m1[1]].data);
                                          for(auto i = 0; i < 3; i++) {
                                              if (auto search = varmap.find(&data1.unit_vec[i]); search != varmap.end())
                                                  NV_Ith_S(ydot, varmap[&lhs[m1[1]].position[i]]) += v_minus * NV_Ith_S(y, search->second) * length_limiter;
                                              else
                                                  NV_Ith_S(ydot, varmap[&lhs[m1[1]].position[i]]) += v_minus * data1.unit_vec[i] * length_limiter;
                                          }
                                          //std::cout << "here end\n";
                                      });

            gamma.addRule(r6);

            //stochastic retraction rule
            GT r7g1;
            r7g1.addNode({1, {Plant::Negative{}}});
            r7g1.addNode({2, {Plant::Intermediate{}}});
            r7g1.addNode({3, {Plant::Intermediate{}}});
            r7g1.addEdge(1, 2);
            r7g1.addEdge(2, 3);

            GT r7g2;
            r7g2.addNode({1, {Plant::Negative{}}});
            r7g2.addNode({3, {Plant::Intermediate{}}});
            r7g2.addEdge(1, 3);

            DGGML::WithRule<GT> r7("retraction", r7g1, r7g2,
                                   [&](auto& lhs, auto& m)
                                   {
                                       auto& node_i_data = lhs.findNode(m[1])->second.getData();
                                       auto& node_j_data = lhs.findNode(m[2])->second.getData();
                                       auto len = DGGML::calculate_distance(node_i_data.position, node_j_data.position);
                                       double propensity = 10.0*DGGML::heaviside(settings.DIV_LENGTH_RETRACT, len);
                                       return propensity;
                                   }, [](auto& lhs, auto& rhs, auto& m1, auto& m2) {
                                        //TODO: reset the unit vector if we have wobble rules
                    });

            gamma.addRule(r7);

//            DGGML::WithRule<GT> r7alt("retraction_alt", r7g1, r7g2,
//                                   [&](auto& lhs, auto& m)
//                                   {
//                                       return 10.0;
//                                   }, [](auto& lhs, auto& rhs, auto& m1, auto& m2) {
//                        //TODO: reset the unit vector if we have wobble rules
//                        rhs[m2[1]].position[0] = (lhs[m1[2]].position[0]);
//                        rhs[m2[1]].position[1] = (lhs[m1[2]].position[1]);
//                        rhs[m2[1]].position[2] = (lhs[m1[2]].position[2]);
//                    });

            //gamma.addRule(r7alt);

            std::cout << "Generating the expanded cell complex\n";
            geoplex2D.init(settings.CELL_NX,
                           settings.CELL_NY,
                           settings.CELL_DX,
                           settings.CELL_DY,
                           settings.GHOSTED,
                           settings.MAXIMAL_REACTION_RADIUS); //ghosted
            //std::cout << geoplex2D;

            std::cout << "Initializing the system graph\n";
            Plant::microtubule_uniform_scatter(system_graph, geoplex2D, settings, gen);
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

        void set_parameters() {
            settings.EXPERIMENT_NAME = "treadmilling";

            //1x1 micrometer domain
            settings.CELL_NX = 1;//1;
            settings.CELL_NY = 1;//1;
            settings.CELL_DX = 0.5;//1.0;
            settings.CELL_DY = 0.5;//1.0;

            //non ghosted complex
            settings.GHOSTED = false;

            //number of microtubules in the simulation
            settings.NUM_MT = 20;

            //starting size of the MTs
            settings.MT_MIN_SEGMENT_INIT = 0.005;
            settings.MT_MAX_SEGMENT_INIT = 0.01;

            settings.LENGTH_DIV_FACTOR = 1.2;
            settings.DIV_LENGTH = 0.026;
            settings.DIV_LENGTH_RETRACT = 0.0026;

            //growing and shrinking velocities in micrometers per second
            settings.V_PLUS = 0.0583;
            settings.V_MINUS = 0.00883;

            settings.SIGMOID_K = 10.0;

            //0.05 micrometers = 50 nanometers
            settings.MAXIMAL_REACTION_RADIUS = 0.05;

            //simulation time in seconds
            settings.TOTAL_TIME = 20.0;
            settings.DELTA = 0.5/8.0; //unit of seconds

            //The internal step of the solver should be at least smaller than delta
            settings.DELTA_DELTA_T = settings.DELTA / 20.0;
            settings.DELTA_T_MIN = settings.DELTA_DELTA_T;
            settings.NUM_STEPS = settings.TOTAL_TIME / settings.DELTA;

            settings.RHO_TEST_RATE = 10.0;
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
            //Delta should be big, but not to big. In this case, the maximum amount of time it would
            //take one MT to grow a single unit of MT
            //e.g. something like: 0.25*settings.MAXIMAL_REACTION_RADIUS / std::max(settings.V_PLUS, settings.V_MINUS);
            settings.DELTA = double(interface["SETTINGS"]["DELTA"]);
            //The internal step of the solver should be at least smaller than delta
            settings.DELTA_DELTA_T = settings.DELTA / 20.0;// / settings.NUM_INTERNAL_STEPS;
            settings.DELTA_T_MIN = settings.DELTA_DELTA_T;
            settings.NUM_STEPS = settings.TOTAL_TIME / settings.DELTA;

            settings.RHO_TEST_RATE = double(interface["EXPERIMENTAL"]["RHO_TEST_RATE"]);
        }
        Parameters settings;
    };
}


#endif //DGGML_CMAMODEL_H
