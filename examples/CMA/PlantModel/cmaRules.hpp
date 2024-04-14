#ifndef DGGML_CMARULES_HPP
#define DGGML_CMARULES_HPP

#include "PlantTypes.hpp"
#include "Parameters.hpp"
//checks if we're outside a box boundary
namespace CMA {
    using GT = Plant::graph_type;

    bool boundary_check_2D(Parameters &settings, double x, double y, double padding = 0.0) {
        auto min_x = 0.0 + padding, min_y = 0.0 + padding;
        auto max_x = settings.CELL_NX * settings.CELL_DX - padding;
        auto max_y = settings.CELL_NY * settings.CELL_DY - padding;

        return (x > min_x && x < max_x && y > min_y && y < max_y) ? false : true;
    }

//checks if we're outside a circular boundary
    bool boundary_check_circle(Parameters &settings, double x, double y) {
        double padding = 2.0 * settings.MAXIMAL_REACTION_RADIUS;
        double diameter = settings.CELL_NX * settings.CELL_DX - 2 * padding;
        double radius = diameter / 2.0;
        double center_x = settings.CELL_NX * settings.CELL_DX / 2.0;
        double center_y = settings.CELL_NY * settings.CELL_DY / 2.0;
        double d = sqrt((x - center_x) * (x - center_x) + (y - center_y) * (y - center_y));
        return (d < radius) ? false : true;
    }

    void create_with_mt_growing_rule(DGGML::Grammar<Plant::graph_type> &gamma, Plant::graph_type &system_graph,
                              Parameters &settings) {
        //stochastic growing rule
        GT mt_growth_lhs;
        mt_growth_lhs.addNode({1, {Plant::Intermediate{}}});
        mt_growth_lhs.addNode({2, {Plant::Positive{}}});
        mt_growth_lhs.addEdge(1, 2);

        GT mt_growth_rhs;
        mt_growth_rhs.addNode({1, {Plant::Intermediate{}}});
        mt_growth_rhs.addNode({3, {Plant::Intermediate{}}});
        mt_growth_rhs.addNode({2, {Plant::Positive{}}});
        mt_growth_rhs.addEdge(1, 3);
        mt_growth_rhs.addEdge(3, 2);


        //TODO: I should make it so that any solving/propensity functions that need access to parameters
        // are actually passed as functors with states!
        DGGML::WithRule<GT> stochastic_mt_growth("stochastic_mt_growth", mt_growth_lhs, mt_growth_rhs,
                                                 [&](auto &lhs, auto &m) {
                                                     //return 0.0*2.0;
                                                     auto &node_i_data = lhs.findNode(m[1])->second.getData();
                                                     auto &node_j_data = lhs.findNode(m[2])->second.getData();
                                                     auto len = DGGML::calculate_distance(node_i_data.position,
                                                                                          node_j_data.position);
                                                     double propensity =
                                                             settings.WITH_GROWTH_RATE_FACTOR * DGGML::heaviside(len, settings.DIV_LENGTH);
                                                     //double propensity = DGGML::sigmoid((len/settings.DIV_LENGTH) - 1.0, settings.SIGMOID_K);
                                                     return propensity;
                                                 }, [&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
                    rhs[m2[3]].position[0] =
                            lhs[m1[2]].position[0] - (lhs[m1[2]].position[0] - lhs[m1[1]].position[0]) / 100.0;
                    rhs[m2[3]].position[1] =
                            lhs[m1[2]].position[1] - (lhs[m1[2]].position[1] - lhs[m1[1]].position[1]) / 100.0;
                    rhs[m2[3]].position[2] =
                            lhs[m1[2]].position[2] - (lhs[m1[2]].position[2] - lhs[m1[1]].position[2]) / 100.0;
                    //next set the unit vector
                    auto& rhs_node3 = rhs[m2[3]];
                    auto& rhs_node3_data = std::get<Plant::Intermediate>(rhs_node3.data);
                    auto& u3 = rhs_node3_data.unit_vec;
                    u3[0] = std::get<Plant::Intermediate>(lhs[m1[1]].data).unit_vec[0];
                    u3[1] = std::get<Plant::Intermediate>(lhs[m1[1]].data).unit_vec[1];
                    u3[2] = std::get<Plant::Intermediate>(lhs[m1[1]].data).unit_vec[2];

                    //randomly wobble the growth direction of u2 on update
                    auto wobble_angle = settings.WOBBLE_ANGLE;
                    if(settings.ENABLE_WOBBLE) {
                        auto &rhs_node2 = rhs[m2[2]];
                        auto &rhs_node2_data = std::get<Plant::Positive>(rhs_node2.data);
                        auto &u2 = rhs_node2_data.unit_vec;
                        //rotate
                        std::random_device random_device;
                        std::mt19937 random_engine(random_device());
                        std::uniform_real_distribution<double> distribution_angle(-wobble_angle * (3.14 / 180.0),
                                                                                  wobble_angle * (3.14 / 180.0));
                        double angle = distribution_angle(random_engine);
                        double x = u2[0];
                        double y = u2[1];
                        u2[0] = x * cos(angle) - y * sin(angle);
                        u2[1] = x * sin(angle) + y * cos(angle);
                    }
                });

        gamma.addRule(stochastic_mt_growth);
    }

    void create_ode_mt_growing_rule(DGGML::Grammar<Plant::graph_type> &gamma, Plant::graph_type &system_graph,
                                       Parameters &settings) {
        GT mt_growth_lhs;
        mt_growth_lhs.addNode({1, {Plant::Intermediate{}}});
        mt_growth_lhs.addNode({2, {Plant::Positive{}}});
        mt_growth_lhs.addEdge(1, 2);

        GT mt_growth_rhs;
        mt_growth_rhs.addNode({1, {Plant::Intermediate{}}});
        mt_growth_rhs.addNode({3, {Plant::Intermediate{}}});
        mt_growth_rhs.addNode({2, {Plant::Positive{}}});
        mt_growth_rhs.addEdge(1, 3);
        mt_growth_rhs.addEdge(3, 2);
        //TODO: I think I need to add velocity back in, and make the growing solve a function of the two ODEs
        DGGML::SolvingRule<GT> ode_mt_growth("ode_mt_growth", mt_growth_lhs, mt_growth_lhs, 3,
                                             [](auto &lhs, auto &m1, auto &varset) {
                                                 //TODO: fix and account for verlet integration
                                                 //bind the variables involved
                                                 varset.insert(&lhs[m1[2]].position[0]);
                                                 varset.insert(&lhs[m1[2]].position[1]);
                                                 varset.insert(&lhs[m1[2]].position[2]);
                                             },
                                             [&](auto &lhs, auto &m1, auto y, auto ydot, auto &varmap) {
                                                 //std::cout << "solving the grow rule\n";
                                                 //unless we know a variable wasn't used before, it's current value from it's solving
                                                 //ode must be checked and used i.e. if(varmap[&lhs[m[1]].position[0]]->second = false) do
                                                 //otherwise we don't have to check, but if the user is wrong, undefined behavior may ensue
                                                 //growth rule params
                                                 auto v_plus = settings.V_PLUS;
                                                 auto d_l = settings.DIV_LENGTH;
                                                 double l = 0.0;
                                                 for (auto i = 0; i < 3; i++) {
                                                     //TODO: set a constraint that stops solving if a boundary is reached
                                                     double diff = NV_Ith_S(y, varmap[&lhs[m1[2]].position[i]]);
                                                     if (auto search = varmap.find(&lhs[m1[1]].position[i]); search !=
                                                                                                             varmap.end())
                                                         diff -= NV_Ith_S(y, search->second);
                                                     else
                                                         diff -= lhs[m1[1]].position[i];
                                                     l += diff * diff;
                                                 }
                                                 l = sqrt(l);
                                                 double length_limiter = 1.0;//(1.0 - (l/d_l));
                                                 auto &data2 = std::get<Plant::Intermediate>(lhs[m1[1]].data);
                                                 auto &data1 = std::get<Plant::Positive>(lhs[m1[2]].data);
                                                 for (auto i = 0; i < 3; i++) {
                                                     if (auto search = varmap.find(&data1.unit_vec[i]); search !=
                                                                                                        varmap.end())
                                                         NV_Ith_S(ydot, varmap[&lhs[m1[2]].position[i]]) +=
                                                                 v_plus * data1.unit_vec[i] *
                                                                 NV_Ith_S(y, search->second) * length_limiter;
                                                     else
                                                         NV_Ith_S(ydot, varmap[&lhs[m1[2]].position[i]]) +=
                                                                 v_plus * data1.unit_vec[i] * length_limiter;
                                                 }

                                                 //TODO: see if we can make this internal to the algorithm and not user controlled
                                                 // so that deactivated ODEs can be removed from the system
                                                 //boundary check
                                                 //TODO: fix and account for verlet integration
                                                 auto x_plus_dx = NV_Ith_S(y, varmap[&lhs[m1[2]].position[0]]) +
                                                                  NV_Ith_S(ydot, varmap[&lhs[m1[2]].position[0]]);
                                                 auto y_plus_dy = NV_Ith_S(y, varmap[&lhs[m1[2]].position[1]]) +
                                                                  NV_Ith_S(ydot, varmap[&lhs[m1[2]].position[1]]);
                                                 bool out_of_bounds = boundary_check_2D(settings, x_plus_dx, y_plus_dy);
                                                 if (false) {//out_of_bounds) {
                                                     for (auto i = 0; i < 3; i++) {
                                                         NV_Ith_S(ydot, varmap[&lhs[m1[2]].position[i]]) = 0.0;
                                                     }
                                                 }
                                             });
            gamma.addRule(ode_mt_growth);
    }

    void create_with_mt_retraction_rule(DGGML::Grammar<Plant::graph_type> &gamma, Plant::graph_type &system_graph,
                                 Parameters &settings)
    {
        //stochastic retraction rule
        GT mt_retraction_lhs1;
        mt_retraction_lhs1.addNode({1, {Plant::Negative{}}});
        mt_retraction_lhs1.addNode({2, {Plant::Intermediate{}}});
        mt_retraction_lhs1.addNode({3, {Plant::Intermediate{}}});
        mt_retraction_lhs1.addEdge(1, 2);
        mt_retraction_lhs1.addEdge(2, 3);

        GT mt_retraction_rhs1;
        mt_retraction_rhs1.addNode({1, {Plant::Negative{}}});
        mt_retraction_rhs1.addNode({3, {Plant::Intermediate{}}});
        mt_retraction_rhs1.addEdge(1, 3);

        DGGML::WithRule<GT> mt_stochastic_retraction("mt_stochastic_retraction", mt_retraction_lhs1, mt_retraction_rhs1,
                                                     [&](auto &lhs, auto &m) {
                                                         auto &node_i_data = lhs.findNode(m[1])->second.getData();
                                                         auto &node_j_data = lhs.findNode(m[2])->second.getData();
                                                         auto len = DGGML::calculate_distance(node_i_data.position,
                                                                                              node_j_data.position);
                                                         double propensity = settings.WITH_RETRACTION_RATE_FACTOR * DGGML::heaviside(
                                                                 settings.DIV_LENGTH_RETRACT, len);
                                                         return propensity;
                                                     }, [](auto &lhs, auto &rhs, auto &m1, auto &m2) {
                    //reset the unit vector
                    DGGML::set_unit_vector(rhs[m2[3]].position, rhs[m2[1]].position,
                                           std::get<Plant::Negative>(rhs[m2[1]].data).unit_vec);
                });


        gamma.addRule(mt_stochastic_retraction);
    }

    void create_ode_mt_retraction_rule(DGGML::Grammar<Plant::graph_type> &gamma, Plant::graph_type &system_graph,
                                 Parameters &settings)
    {
        GT mt_retraction_lhs2;
        mt_retraction_lhs2.addNode({1, {Plant::Negative{}}});
        mt_retraction_lhs2.addNode({2, {Plant::Intermediate{}}});
        mt_retraction_lhs2.addEdge(1, 2);
        //TODO: I think I need to add velocity back in, and make the growing solve a function of the two ODEs
        DGGML::SolvingRule<GT> mt_ode_retraction("mt_ode_retraction", mt_retraction_lhs2, mt_retraction_lhs2, 3,
                                                 [](auto &lhs, auto &m1, auto &varset) {
                                                     //std::cout << "ic of the grow rule\n";
                                                     //bind the variables involved
                                                     varset.insert(&lhs[m1[1]].position[0]);
                                                     varset.insert(&lhs[m1[1]].position[1]);
                                                     varset.insert(&lhs[m1[1]].position[2]);
                                                 },
                                                 [&](auto &lhs, auto &m1, auto y, auto ydot, auto &varmap) {
                                                     //std::cout << "here start\n";
                                                     auto d_l_r = settings.DIV_LENGTH_RETRACT;
                                                     auto v_minus = settings.V_MINUS;
                                                     auto d_l = settings.DIV_LENGTH;
                                                     double l = 0.0;
                                                     for (auto i = 0; i < 3; i++) {
                                                         double diff = NV_Ith_S(y, varmap[&lhs[m1[1]].position[i]]);
                                                         if (auto search = varmap.find(&lhs[m1[2]].position[i]);
                                                                 search != varmap.end())
                                                             diff -= NV_Ith_S(y, search->second);
                                                         else
                                                             diff -= lhs[m1[2]].position[i];
                                                         l += diff * diff;
                                                     }
                                                     l = sqrt(l);
                                                     double length_limiter = l / d_l;
                                                     //if(length_limiter <= d_l_r) length_limiter = 0.0;
                                                     if (l <= d_l_r / 2.0) length_limiter = 0.0;
                                                     else length_limiter = 1.0;
                                                     auto &data1 = std::get<Plant::Negative>(lhs[m1[1]].data);
                                                     for (auto i = 0; i < 3; i++) {
                                                         if (auto search = varmap.find(&data1.unit_vec[i]); search !=
                                                                                                            varmap.end())
                                                             NV_Ith_S(ydot, varmap[&lhs[m1[1]].position[i]]) +=
                                                                     v_minus * NV_Ith_S(y, search->second) *
                                                                     length_limiter;
                                                         else
                                                             NV_Ith_S(ydot, varmap[&lhs[m1[1]].position[i]]) +=
                                                                     v_minus * data1.unit_vec[i] * length_limiter;
                                                     }
                                                     //std::cout << "here end\n";
                                                 });
        if(settings.ENABLE_RETRACTION)
            gamma.addRule(mt_ode_retraction);

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
    }

    //case one is where an MT his the normal boundary
    void create_with_standard_boundary_catastrophe(DGGML::Grammar<Plant::graph_type> &gamma, Plant::graph_type &system_graph,
                               Parameters &settings) {
        //alternatively to boundary checking we can have rules if the boundary is a graph
        GT boundary_catastrophe_lhs1;
        boundary_catastrophe_lhs1.addNode({1, {Plant::Intermediate{}}});
        boundary_catastrophe_lhs1.addNode({2, {Plant::Positive{}}});
        boundary_catastrophe_lhs1.addEdge(1, 2);
        boundary_catastrophe_lhs1.addNode({3, {Plant::Boundary{}}});
        boundary_catastrophe_lhs1.addNode({4, {Plant::Boundary{}}});
        boundary_catastrophe_lhs1.addEdge(3, 4);
        GT boundary_catastrophe_rhs1;
        boundary_catastrophe_rhs1.addNode({1, {Plant::Intermediate{}}});
        boundary_catastrophe_rhs1.addNode({2, {Plant::Negative{}}});
        boundary_catastrophe_rhs1.addEdge(1, 2);
        boundary_catastrophe_rhs1.addNode({3, {Plant::Boundary{}}});
        boundary_catastrophe_rhs1.addNode({4, {Plant::Boundary{}}});
        boundary_catastrophe_rhs1.addEdge(3, 4);

        DGGML::WithRule<GT> boundary_catastrophe1("boundary_catastrophe1", boundary_catastrophe_lhs1,
                                                  boundary_catastrophe_rhs1,
                                                  [&](auto &lhs, auto &m) {
                                                      auto &dat1 = lhs.findNode(m[1])->second.getData();
                                                      auto &dat2 = lhs.findNode(m[2])->second.getData();
                                                      auto &dat3 = lhs.findNode(m[3])->second.getData();
                                                      auto &dat4 = lhs.findNode(m[4])->second.getData();

                                                      //get references to position vector
                                                      auto &pos1 = dat1.position;
                                                      auto &pos2 = dat2.position;
                                                      auto &pos3 = dat3.position;
                                                      auto &pos4 = dat4.position;

                                                      //get references to unit vectors
                                                      auto &u1 = std::get<Plant::Intermediate>(dat1.data).unit_vec;
                                                      auto &u2 = std::get<Plant::Positive>(dat2.data).unit_vec;
                                                      //auto& u3 = std::get<Plant::Boundary>(dat3.data).unit_vec;
                                                      //auto& u4 = std::get<Plant::Boundary>(dat4.data).unit_vec;

                                                      //distance from positive node to line segment
                                                      auto d = DGGML::distanceToLineSegment(pos3[0], pos3[1], pos4[0],
                                                                                            pos4[1], pos2[0], pos2[1]);
                                                      if (d <= settings.COLLISION_DISTANCE)
                                                          return settings.STANDARD_BOUNDARY_CATASTROPHE_RATE;
                                                      else
                                                          return 0.0;
                                                  },
                                                  [](auto &lhs, auto &rhs, auto &m1, auto &m2) {
                                                      std::get<Plant::Negative>(
                                                              rhs[m2[2]].data).unit_vec[0] = -std::get<Plant::Positive>(
                                                              lhs[m1[2]].data).unit_vec[0];
                                                      std::get<Plant::Negative>(
                                                              rhs[m2[2]].data).unit_vec[1] = -std::get<Plant::Positive>(
                                                              lhs[m1[2]].data).unit_vec[1];
                                                      std::get<Plant::Negative>(
                                                              rhs[m2[2]].data).unit_vec[2] = -std::get<Plant::Positive>(
                                                              lhs[m1[2]].data).unit_vec[2];
                                                  });
        gamma.addRule(boundary_catastrophe1);
    }

    //case one is where an MT his the clasp boundary
    void create_with_clasp_boundary_catastrophe(DGGML::Grammar<Plant::graph_type> &gamma, Plant::graph_type &system_graph,
                                                   Parameters &settings) {
        //temp rule for clasp
        GT boundary_catastrophe_lhs2;
        boundary_catastrophe_lhs2.addNode({1, {Plant::Intermediate{}}});
        boundary_catastrophe_lhs2.addNode({2, {Plant::Positive{}}});
        boundary_catastrophe_lhs2.addEdge(1, 2);
        boundary_catastrophe_lhs2.addNode({3, {Plant::Boundary{}}});
        boundary_catastrophe_lhs2.addNode({4, {Plant::Intermediate{}}});
        boundary_catastrophe_lhs2.addEdge(3, 4);
        GT boundary_catastrophe_rhs2;
        boundary_catastrophe_rhs2.addNode({1, {Plant::Intermediate{}}});
        boundary_catastrophe_rhs2.addNode({2, {Plant::Negative{}}});
        boundary_catastrophe_rhs2.addEdge(1, 2);
        boundary_catastrophe_rhs2.addNode({3, {Plant::Boundary{}}});
        boundary_catastrophe_rhs2.addNode({4, {Plant::Intermediate{}}});
        boundary_catastrophe_rhs2.addEdge(3, 4);

        DGGML::WithRule<GT> boundary_catastrophe2("boundary_catastrophe2", boundary_catastrophe_lhs2,
                                                  boundary_catastrophe_rhs2,
                                                  [&](auto &lhs, auto &m) {
                                                      auto &dat1 = lhs.findNode(m[1])->second.getData();
                                                      auto &dat2 = lhs.findNode(m[2])->second.getData();
                                                      auto &dat3 = lhs.findNode(m[3])->second.getData();
                                                      auto &dat4 = lhs.findNode(m[4])->second.getData();

                                                      //get references to position vector
                                                      auto &pos1 = dat1.position;
                                                      auto &pos2 = dat2.position;
                                                      auto &pos3 = dat3.position;
                                                      auto &pos4 = dat4.position;

                                                      //get references to unit vectors
                                                      auto &u1 = std::get<Plant::Intermediate>(dat1.data).unit_vec;
                                                      auto &u2 = std::get<Plant::Positive>(dat2.data).unit_vec;
                                                      //auto& u3 = std::get<Plant::Boundary>(dat3.data).unit_vec;
                                                      //auto& u4 = std::get<Plant::Boundary>(dat4.data).unit_vec;

                                                      //distance from positive node to line segment
                                                      auto d = DGGML::distanceToLineSegment(pos3[0], pos3[1], pos4[0],
                                                                                            pos4[1], pos2[0], pos2[1]);
                                                      if (d <= settings.COLLISION_DISTANCE)
                                                          return settings.CLASP_BOUNDARY_CATASTROPHE_RATE;
                                                      else
                                                          return 0.0;
                                                  },
                                                  [](auto &lhs, auto &rhs, auto &m1, auto &m2) {
                                                      std::get<Plant::Negative>(
                                                              rhs[m2[2]].data).unit_vec[0] = -std::get<Plant::Positive>(
                                                              lhs[m1[2]].data).unit_vec[0];
                                                      std::get<Plant::Negative>(
                                                              rhs[m2[2]].data).unit_vec[1] = -std::get<Plant::Positive>(
                                                              lhs[m1[2]].data).unit_vec[1];
                                                      std::get<Plant::Negative>(
                                                              rhs[m2[2]].data).unit_vec[2] = -std::get<Plant::Positive>(
                                                              lhs[m1[2]].data).unit_vec[2];
                                                  });
        gamma.addRule(boundary_catastrophe2);
    }
    // intermediate segment hit for collision induced catastrophe
    void create_with_intermediate_cic(DGGML::Grammar<Plant::graph_type> &gamma, Plant::graph_type &system_graph,
                                Parameters &settings) {
        GT catastrophe2_lhs_graph1;
        catastrophe2_lhs_graph1.addNode({1, {Plant::Intermediate{}}});
        catastrophe2_lhs_graph1.addNode({2, {Plant::Positive{}}});
        catastrophe2_lhs_graph1.addNode({3, {Plant::Intermediate{}}});
        catastrophe2_lhs_graph1.addNode({4, {Plant::Intermediate{}}});
        catastrophe2_lhs_graph1.addEdge(1, 2);
        catastrophe2_lhs_graph1.addEdge(3, 4);

        GT catastrophe2_rhs_graph1;
        catastrophe2_rhs_graph1.addNode({1, {Plant::Intermediate{}}});
        catastrophe2_rhs_graph1.addNode({2, {Plant::Negative{}}});
        catastrophe2_rhs_graph1.addNode({3, {Plant::Intermediate{}}});
        catastrophe2_rhs_graph1.addNode({4, {Plant::Intermediate{}}});
        catastrophe2_rhs_graph1.addEdge(1, 2);
        catastrophe2_rhs_graph1.addEdge(3, 4);

        DGGML::WithRule<GT> catastrophe2_case1("catastrophe2_case1", catastrophe2_lhs_graph1, catastrophe2_rhs_graph1,
                                               [&](auto &lhs, auto &m) {
                                                   //find all the node data
                                                   auto &dat1 = lhs.findNode(m[1])->second.getData();
                                                   auto &dat2 = lhs.findNode(m[2])->second.getData();
                                                   auto &dat3 = lhs.findNode(m[3])->second.getData();
                                                   auto &dat4 = lhs.findNode(m[4])->second.getData();

                                                   //get references to position vector
                                                   auto &pos1 = dat1.position;
                                                   auto &pos2 = dat2.position;
                                                   auto &pos3 = dat3.position;
                                                   auto &pos4 = dat4.position;

                                                   //get references to unit vectors
                                                   auto &u1 = std::get<Plant::Intermediate>(dat1.data).unit_vec;
                                                   auto &u2 = std::get<Plant::Positive>(dat2.data).unit_vec;
                                                   auto &u3 = std::get<Plant::Intermediate>(dat3.data).unit_vec;
                                                   auto &u4 = std::get<Plant::Intermediate>(dat4.data).unit_vec;

                                                   //distance from positive node to line segment
                                                   auto d = DGGML::distanceToLineSegment(pos3[0], pos3[1], pos4[0],
                                                                                         pos4[1], pos2[0], pos2[1]);
                                                   if (d <= settings.COLLISION_DISTANCE) {
                                                       auto theta = DGGML::compute_theta(u2, dat2.position, u3);
                                                       if (theta < settings.CATASTROPHE_ANGLE) {
                                                           double sol[2];
                                                           DGGML::paramaterized_intersection(pos2, pos4, pos3, u2, sol);

                                                           //TODO: determine how strict we need to be with the collision point
                                                           if (sol[0] > 0.0) {// && sol[1] >= 0.0 && sol[1] <= 1.0) {
                                                               return settings.INTERMEDIATE_CIC_RATE;
                                                           }
                                                       }
                                                   }
                                                   return 0.0;
                                               },
                                               [](auto &lhs, auto &rhs, auto &m1, auto &m2) {
                                                   std::get<Plant::Negative>(
                                                           rhs[m2[2]].data).unit_vec[0] = -std::get<Plant::Positive>(
                                                           lhs[m1[2]].data).unit_vec[0];
                                                   std::get<Plant::Negative>(
                                                           rhs[m2[2]].data).unit_vec[1] = -std::get<Plant::Positive>(
                                                           lhs[m1[2]].data).unit_vec[1];
                                                   std::get<Plant::Negative>(
                                                           rhs[m2[2]].data).unit_vec[2] = -std::get<Plant::Positive>(
                                                           lhs[m1[2]].data).unit_vec[2];
                                               });

        gamma.addRule(catastrophe2_case1);
    }

    // intermediate segment hit for collision induced catastrophe
    void create_with_positive_cic(DGGML::Grammar<Plant::graph_type> &gamma, Plant::graph_type &system_graph,
                                      Parameters &settings) {
        GT catastrophe2_lhs_graph2;
        catastrophe2_lhs_graph2.addNode({1, {Plant::Intermediate{}}});
        catastrophe2_lhs_graph2.addNode({2, {Plant::Positive{}}});
        catastrophe2_lhs_graph2.addNode({3, {Plant::Intermediate{}}});
        catastrophe2_lhs_graph2.addNode({4, {Plant::Positive{}}});
        catastrophe2_lhs_graph2.addEdge(1, 2);
        catastrophe2_lhs_graph2.addEdge(3, 4);

        GT catastrophe2_rhs_graph2;
        catastrophe2_rhs_graph2.addNode({1, {Plant::Intermediate{}}});
        catastrophe2_rhs_graph2.addNode({2, {Plant::Negative{}}});
        catastrophe2_rhs_graph2.addNode({3, {Plant::Intermediate{}}});
        catastrophe2_rhs_graph2.addNode({4, {Plant::Positive{}}});
        catastrophe2_rhs_graph2.addEdge(1, 2);
        catastrophe2_rhs_graph2.addEdge(3, 4);

        DGGML::WithRule<GT> catastrophe2_case2("catastrophe2_case2", catastrophe2_lhs_graph2, catastrophe2_rhs_graph2,
                                               [&](auto &lhs, auto &m) {
                                                   //find all the node data
                                                   auto &dat1 = lhs.findNode(m[1])->second.getData();
                                                   auto &dat2 = lhs.findNode(m[2])->second.getData();
                                                   auto &dat3 = lhs.findNode(m[3])->second.getData();
                                                   auto &dat4 = lhs.findNode(m[4])->second.getData();

                                                   //get references to position vector
                                                   auto &pos1 = dat1.position;
                                                   auto &pos2 = dat2.position;
                                                   auto &pos3 = dat3.position;
                                                   auto &pos4 = dat4.position;

                                                   //get references to unit vectors
                                                   auto &u1 = std::get<Plant::Intermediate>(dat1.data).unit_vec;
                                                   auto &u2 = std::get<Plant::Positive>(dat2.data).unit_vec;
                                                   auto &u3 = std::get<Plant::Intermediate>(dat3.data).unit_vec;
                                                   auto &u4 = std::get<Plant::Positive>(dat4.data).unit_vec;

                                                   //distance from positive node to line segment
                                                   auto d = DGGML::distanceToLineSegment(pos3[0], pos3[1], pos4[0],
                                                                                         pos4[1], pos2[0], pos2[1]);
                                                   if (d <= settings.COLLISION_DISTANCE / 2.0)
                                                       return settings.POSITIVE_CIC_RATE;
                                                   else return 0.0;
                                               },
                                               [](auto &lhs, auto &rhs, auto &m1, auto &m2) {
                                                   std::get<Plant::Negative>(
                                                           rhs[m2[2]].data).unit_vec[0] = -std::get<Plant::Positive>(
                                                           lhs[m1[2]].data).unit_vec[0];
                                                   std::get<Plant::Negative>(
                                                           rhs[m2[2]].data).unit_vec[1] = -std::get<Plant::Positive>(
                                                           lhs[m1[2]].data).unit_vec[1];
                                                   std::get<Plant::Negative>(
                                                           rhs[m2[2]].data).unit_vec[2] = -std::get<Plant::Positive>(
                                                           lhs[m1[2]].data).unit_vec[2];
                                               });

        gamma.addRule(catastrophe2_case2);
    }

    // intermediate segment hit for collision induced catastrophe
    void create_with_negative_cic(DGGML::Grammar<Plant::graph_type> &gamma, Plant::graph_type &system_graph,
                                      Parameters &settings) {
        GT catastrophe2_lhs_graph3;
        catastrophe2_lhs_graph3.addNode({1, {Plant::Intermediate{}}});
        catastrophe2_lhs_graph3.addNode({2, {Plant::Positive{}}});
        catastrophe2_lhs_graph3.addNode({3, {Plant::Intermediate{}}});
        catastrophe2_lhs_graph3.addNode({4, {Plant::Negative{}}});
        catastrophe2_lhs_graph3.addEdge(1, 2);
        catastrophe2_lhs_graph3.addEdge(3, 4);

        GT catastrophe2_rhs_graph3;
        catastrophe2_rhs_graph3.addNode({1, {Plant::Intermediate{}}});
        catastrophe2_rhs_graph3.addNode({2, {Plant::Negative{}}});
        catastrophe2_rhs_graph3.addNode({3, {Plant::Intermediate{}}});
        catastrophe2_rhs_graph3.addNode({4, {Plant::Negative{}}});
        catastrophe2_rhs_graph3.addEdge(1, 2);
        catastrophe2_rhs_graph3.addEdge(3, 4);

        DGGML::WithRule<GT> catastrophe2_case3("catastrophe2_case3", catastrophe2_lhs_graph3, catastrophe2_rhs_graph3,
                                               [&](auto &lhs, auto &m) {
                                                   //find all the node data
                                                   auto &dat1 = lhs.findNode(m[1])->second.getData();
                                                   auto &dat2 = lhs.findNode(m[2])->second.getData();
                                                   auto &dat3 = lhs.findNode(m[3])->second.getData();
                                                   auto &dat4 = lhs.findNode(m[4])->second.getData();

                                                   //get references to position vector
                                                   auto &pos1 = dat1.position;
                                                   auto &pos2 = dat2.position;
                                                   auto &pos3 = dat3.position;
                                                   auto &pos4 = dat4.position;

                                                   //get references to unit vectors
                                                   auto &u1 = std::get<Plant::Intermediate>(dat1.data).unit_vec;
                                                   auto &u2 = std::get<Plant::Positive>(dat2.data).unit_vec;
                                                   auto &u3 = std::get<Plant::Intermediate>(dat3.data).unit_vec;
                                                   auto &u4 = std::get<Plant::Negative>(dat4.data).unit_vec;

                                                   //distance from positive node to line segment
                                                   auto d = DGGML::distanceToLineSegment(pos3[0], pos3[1], pos4[0],
                                                                                         pos4[1], pos2[0], pos2[1]);
                                                   if (d <= settings.COLLISION_DISTANCE / 2.0)
                                                       return settings.NEGATIVE_CIC_RATE;
                                                   else return 0.0;
                                               },
                                               [](auto &lhs, auto &rhs, auto &m1, auto &m2) {
                                                   std::get<Plant::Negative>(
                                                           rhs[m2[2]].data).unit_vec[0] = -std::get<Plant::Positive>(
                                                           lhs[m1[2]].data).unit_vec[0];
                                                   std::get<Plant::Negative>(
                                                           rhs[m2[2]].data).unit_vec[1] = -std::get<Plant::Positive>(
                                                           lhs[m1[2]].data).unit_vec[1];
                                                   std::get<Plant::Negative>(
                                                           rhs[m2[2]].data).unit_vec[2] = -std::get<Plant::Positive>(
                                                           lhs[m1[2]].data).unit_vec[2];
                                               });

        gamma.addRule(catastrophe2_case3);
    }

    void create_with_crossover_rule(DGGML::Grammar<Plant::graph_type> &gamma, Plant::graph_type &system_graph,
                                    Parameters &settings) {
        //            GT crossover_lhs_graph;
//            crossover_lhs_graph.addNode({1, {Plant::Intermediate{}}});
//            crossover_lhs_graph.addNode({2, {Plant::Positive{}}});
//            crossover_lhs_graph.addNode({3, {Plant::Intermediate{}}});
//            crossover_lhs_graph.addNode({4, {Plant::Intermediate{}}});
//            crossover_lhs_graph.addEdge(1, 2);
//            crossover_lhs_graph.addEdge(3, 4);
//
//            GT crossover_rhs_graph;
//            crossover_rhs_graph.addNode({1, {Plant::Intermediate{}}});
//            crossover_rhs_graph.addNode({2, {Plant::Junction{}}});
//            crossover_rhs_graph.addNode({5, {Plant::Positive{}}});
//            crossover_rhs_graph.addNode({3, {Plant::Intermediate{}}});
//            crossover_rhs_graph.addNode({4, {Plant::Intermediate{}}});
//            crossover_rhs_graph.addEdge(1, 2);
//            crossover_rhs_graph.addEdge(3, 2);
//            crossover_rhs_graph.addEdge(4, 2);
//            crossover_rhs_graph.addEdge(5, 2);
//
//            DGGML::WithRule<GT> crossover("crossover", crossover_lhs_graph, crossover_rhs_graph,
//                                                   [&](auto& lhs, auto& m) {
//                                                       //find all the node data
//                                                       auto& dat1 = lhs.findNode(m[1])->second.getData();
//                                                       auto& dat2 = lhs.findNode(m[2])->second.getData();
//                                                       auto& dat3 = lhs.findNode(m[3])->second.getData();
//                                                       auto& dat4 = lhs.findNode(m[4])->second.getData();
//
//                                                       //get references to position vector
//                                                       auto& pos1 = dat1.position;
//                                                       auto& pos2 = dat2.position;
//                                                       auto& pos3 = dat3.position;
//                                                       auto& pos4 = dat4.position;
//
//                                                       //get references to unit vectors
//                                                       auto& u1 = std::get<Plant::Intermediate>(dat1.data).unit_vec;
//                                                       auto& u2 = std::get<Plant::Positive>(dat2.data).unit_vec;
//                                                       auto& u3 = std::get<Plant::Intermediate>(dat3.data).unit_vec;
//                                                       auto& u4 = std::get<Plant::Intermediate>(dat4.data).unit_vec;
//
//                                                       //distance from positive node to line segment
//                                                       auto d = DGGML::distanceToLineSegment(pos3[0], pos3[1], pos4[0], pos4[1], pos2[0], pos2[1]);
//                                                       if(d <= 0.025) {
//                                                           auto theta = DGGML::compute_theta(u2, dat2.position, u3);
//                                                           if(theta < 45.0)
//                                                           {
//                                                               double sol[2];
//                                                               DGGML::paramaterized_intersection(pos2, pos4, pos3, u2, sol);
//
//                                                               //TODO: determine how strict we need to be with the collision point
//                                                               if(sol[0] > 0.0){// && sol[1] >= 0.0 && sol[1] <= 1.0) {
//                                                                   return 1000.0;
//                                                               }
//                                                           }
//                                                       }
//                                                       return 0.0;
//                                                   },
//                                                   [](auto& lhs, auto& rhs, auto& m1, auto& m2)
//                                                   {
//                                                       for(int i = 0; i < 3; i++)
//                                                        rhs[m2[2]].position[i] = (lhs[m2[4]].position[i] + lhs[m2[3]].position[i])/2.0;
//                                                       DGGML::set_unit_vector(rhs[m2[2]].position, rhs[m2[1]].position, std::get<Plant::Junction>(rhs[m2[2]].data).unit_vec);
//                                                       DGGML::set_unit_vector(rhs[m2[2]].position, rhs[m2[1]].position, std::get<Plant::Intermediate>(rhs[m2[1]].data).unit_vec);
//                                                       for(int i = 0; i < 3; i++)
//                                                        rhs[m2[5]].position[i] = rhs[m1[2]].position[i]+0.01*std::get<Plant::Junction>(rhs[m1[2]].data).unit_vec[i];
//                                                       DGGML::set_unit_vector(rhs[m2[5]].position, rhs[m2[2]].position, std::get<Plant::Positive>(rhs[m2[5]].data).unit_vec);
//
//                                                   });
//
//            gamma.addRule(crossover);
        GT crossover_lhs_graph;
        crossover_lhs_graph.addNode({1, {Plant::Intermediate{}}});
        crossover_lhs_graph.addNode({2, {Plant::Positive{}}});
        crossover_lhs_graph.addNode({3, {Plant::Intermediate{}}});
        crossover_lhs_graph.addNode({4, {Plant::Intermediate{}}});
        crossover_lhs_graph.addEdge(1, 2);
        crossover_lhs_graph.addEdge(3, 4);

        GT crossover_rhs_graph;
        crossover_rhs_graph.addNode({1, {Plant::Intermediate{}}});
        crossover_rhs_graph.addNode({2, {Plant::Junction{}}});
        crossover_rhs_graph.addNode({5, {Plant::Intermediate{}}});
        crossover_rhs_graph.addNode({6, {Plant::Positive{}}});
        crossover_rhs_graph.addNode({3, {Plant::Intermediate{}}});
        crossover_rhs_graph.addNode({4, {Plant::Intermediate{}}});
        crossover_rhs_graph.addEdge(1, 2);
        crossover_rhs_graph.addEdge(3, 2);
        crossover_rhs_graph.addEdge(4, 2);
        crossover_rhs_graph.addEdge(5, 2);
        crossover_rhs_graph.addEdge(5, 6);

        DGGML::WithRule<GT> crossover("crossover", crossover_lhs_graph, crossover_rhs_graph,
                                      [&](auto &lhs, auto &m) {
                                          //find all the node data
                                          auto &dat1 = lhs.findNode(m[1])->second.getData();
                                          auto &dat2 = lhs.findNode(m[2])->second.getData();
                                          auto &dat3 = lhs.findNode(m[3])->second.getData();
                                          auto &dat4 = lhs.findNode(m[4])->second.getData();

                                          //get references to position vector
                                          auto &pos1 = dat1.position;
                                          auto &pos2 = dat2.position;
                                          auto &pos3 = dat3.position;
                                          auto &pos4 = dat4.position;

                                          //get references to unit vectors
                                          auto &u1 = std::get<Plant::Intermediate>(dat1.data).unit_vec;
                                          auto &u2 = std::get<Plant::Positive>(dat2.data).unit_vec;
                                          auto &u3 = std::get<Plant::Intermediate>(dat3.data).unit_vec;
                                          auto &u4 = std::get<Plant::Intermediate>(dat4.data).unit_vec;

                                          //distance from positive node to line segment
                                          auto d = DGGML::distanceToLineSegment(pos3[0], pos3[1], pos4[0], pos4[1],
                                                                                pos2[0], pos2[1]);
                                          if (d <= settings.COLLISION_DISTANCE) {
                                              auto theta = DGGML::compute_theta(u2, dat2.position, u3);
                                              if (theta < settings.CROSSOVER_ANGLE) {
                                                  double sol[2];
                                                  DGGML::paramaterized_intersection(pos2, pos4, pos3, u2, sol);

                                                  //TODO: determine how strict we need to be with the collision point
                                                  if (sol[0] > 0.0 && sol[1] >= 0.0 && sol[1] <= 1.0) {
                                                      return settings.CROSSOVER_RATE;
                                                  }
                                              }
                                          }
                                          return 0.0;
                                      },
                                      [](auto &lhs, auto &rhs, auto &m1, auto &m2) {
                                          auto &dat2 = lhs.findNode(m1[2])->second.getData();
                                          auto &dat3 = lhs.findNode(m1[3])->second.getData();
                                          auto &dat4 = lhs.findNode(m1[4])->second.getData();
                                          auto &pos2 = dat2.position;
                                          auto &pos3 = dat3.position;
                                          auto &pos4 = dat4.position;
                                          auto &u2 = std::get<Plant::Positive>(dat2.data).unit_vec;
                                          auto &u3 = std::get<Plant::Intermediate>(dat3.data).unit_vec;
                                          double sol[2];
                                          DGGML::paramaterized_intersection(pos2, pos4, pos3, u2, sol);
                                          for (int i = 0; i < 3; i++)
                                              rhs[m2[2]].position[i] = lhs[m2[4]].position[i] +
                                                                       (lhs[m2[3]].position[i] -
                                                                        lhs[m2[4]].position[i]) * sol[1];
                                          DGGML::set_unit_vector(rhs[m2[2]].position, rhs[m2[1]].position,
                                                                 std::get<Plant::Junction>(rhs[m2[2]].data).unit_vec);
                                          DGGML::set_unit_vector(rhs[m2[2]].position, rhs[m2[1]].position,
                                                                 std::get<Plant::Intermediate>(
                                                                         rhs[m2[1]].data).unit_vec);
                                          for (int i = 0; i < 3; i++)
                                              rhs[m2[5]].position[i] = rhs[m1[2]].position[i] + 0.01 *
                                                                                                std::get<Plant::Junction>(
                                                                                                        rhs[m1[2]].data).unit_vec[i];
                                          DGGML::set_unit_vector(rhs[m2[5]].position, rhs[m2[2]].position,
                                                                 std::get<Plant::Intermediate>(
                                                                         rhs[m2[5]].data).unit_vec);
                                          for (int i = 0; i < 3; i++)
                                              rhs[m2[6]].position[i] = rhs[m1[2]].position[i] + 0.011 *
                                                                                                std::get<Plant::Junction>(
                                                                                                        rhs[m1[2]].data).unit_vec[i];
                                          DGGML::set_unit_vector(rhs[m2[6]].position, rhs[m2[2]].position,
                                                                 std::get<Plant::Positive>(rhs[m2[6]].data).unit_vec);

                                      });

        gamma.addRule(crossover);
    }

    void create_with_zippering_rules(DGGML::Grammar<Plant::graph_type> &gamma, Plant::graph_type &system_graph,
                                Parameters &settings) {
        GT zippering_lhs_graph1;
        zippering_lhs_graph1.addNode({1, {Plant::Intermediate{}}});
        zippering_lhs_graph1.addNode({2, {Plant::Positive{}}});
        zippering_lhs_graph1.addNode({3, {Plant::Intermediate{}}});
        zippering_lhs_graph1.addNode({4, {Plant::Intermediate{}}});
        zippering_lhs_graph1.addEdge(1, 2);
        zippering_lhs_graph1.addEdge(3, 4);

        GT zippering_rhs_graph1;
        zippering_rhs_graph1.addNode({1, {Plant::Intermediate{}}});
        zippering_rhs_graph1.addNode({2, {Plant::Zipper{}}});
        zippering_rhs_graph1.addNode({6, {Plant::Intermediate{}}});
        zippering_rhs_graph1.addNode({5, {Plant::Positive{}}});
        zippering_rhs_graph1.addNode({3, {Plant::Intermediate{}}});
        zippering_rhs_graph1.addNode({4, {Plant::Intermediate{}}});
        zippering_rhs_graph1.addEdge(1, 2);
        zippering_rhs_graph1.addEdge(2, 6);
        zippering_rhs_graph1.addEdge(6, 5);
        zippering_rhs_graph1.addEdge(3, 4);

        DGGML::WithRule<GT> zippering_case1("zippering_case1", zippering_lhs_graph1, zippering_rhs_graph1,
                                            [&](auto &lhs, auto &m) {
                                                //find all the node data
                                                auto &dat1 = lhs.findNode(m[1])->second.getData();
                                                auto &dat2 = lhs.findNode(m[2])->second.getData();
                                                auto &dat3 = lhs.findNode(m[3])->second.getData();
                                                auto &dat4 = lhs.findNode(m[4])->second.getData();

                                                //get references to position vector
                                                auto &pos1 = dat1.position;
                                                auto &pos2 = dat2.position;
                                                auto &pos3 = dat3.position;
                                                auto &pos4 = dat4.position;

                                                //get references to unit vectors
                                                auto &u1 = std::get<Plant::Intermediate>(dat1.data).unit_vec;
                                                auto &u2 = std::get<Plant::Positive>(dat2.data).unit_vec;
                                                auto &u3 = std::get<Plant::Intermediate>(dat3.data).unit_vec;
                                                auto &u4 = std::get<Plant::Intermediate>(dat4.data).unit_vec;

                                                //distance from positive node to line segment
                                                auto d = DGGML::distanceToLineSegment(pos3[0], pos3[1], pos4[0],
                                                                                      pos4[1], pos2[0], pos2[1]);
                                                if (d <= settings.SEPARATION_DISTANCE) {
                                                    //double theta = std::acos(DGGML::unit_dot_product(u2, u3))*(180.0/3.14159265);
                                                    //theta = std::min(180.0 - theta, theta);
                                                    //if(theta > 2.0 && theta < 80.0)
                                                    auto theta = DGGML::compute_theta(u2, dat2.position, u3);
                                                    if (theta >= settings.CRITICAL_ANGLE && theta < 88.0) {
                                                        double sol[2];
                                                        DGGML::paramaterized_intersection(pos2, pos4, pos3, u2, sol);

                                                        //TODO: determine how strict we need to be with the collision point
                                                        if (sol[0] > 0.0) {// && sol[1] >= 0.0 && sol[1] <= 1.0) {
                                                            return settings.ZIPPERING_HIT_RATE;
                                                        }
                                                    }
                                                }
                                                return 0.0;
                                            },
                                            [](auto &lhs, auto &rhs, auto &m1, auto &m2) {
                                                auto &dat2 = lhs.findNode(m1[2])->second.getData();
                                                auto &dat3 = lhs.findNode(m1[3])->second.getData();
                                                auto &u2 = std::get<Plant::Positive>(dat2.data).unit_vec;
                                                auto &u3 = std::get<Plant::Intermediate>(dat3.data).unit_vec;
                                                double par_u[3];
                                                //DGGML::perfect_deflection(u2, dat2.position, u3, par_u);
                                                auto theta = DGGML::compute_theta(u2, dat2.position, u3);
                                                //std::cout << theta << "\n";
                                                //std::cin.get();
                                                DGGML::parallel_deflection(u2, dat2.position, u3, par_u);
                                                std::get<Plant::Zipper>(
                                                        rhs[m2[2]].data).unit_vec[0] = std::get<Plant::Intermediate>(
                                                        lhs[m1[3]].data).unit_vec[0];
                                                std::get<Plant::Zipper>(
                                                        rhs[m2[2]].data).unit_vec[1] = std::get<Plant::Intermediate>(
                                                        lhs[m1[3]].data).unit_vec[1];
                                                std::get<Plant::Zipper>(
                                                        rhs[m2[2]].data).unit_vec[2] = std::get<Plant::Intermediate>(
                                                        lhs[m1[3]].data).unit_vec[2];

                                                std::get<Plant::Intermediate>(rhs[m2[6]].data).unit_vec[0] = par_u[0];
                                                std::get<Plant::Intermediate>(rhs[m2[6]].data).unit_vec[1] = par_u[1];
                                                std::get<Plant::Intermediate>(rhs[m2[6]].data).unit_vec[2] = par_u[2];

                                                rhs[m2[6]].position[0] = lhs[m1[2]].position[0] + 0.005 *
                                                                                                  std::get<Plant::Intermediate>(
                                                                                                          rhs[m2[6]].data).unit_vec[0];
                                                rhs[m2[6]].position[1] = lhs[m1[2]].position[1] + 0.005 *
                                                                                                  std::get<Plant::Intermediate>(
                                                                                                          rhs[m2[6]].data).unit_vec[1];
                                                rhs[m2[6]].position[2] = lhs[m1[2]].position[2] + 0.005 *
                                                                                                  std::get<Plant::Intermediate>(
                                                                                                          rhs[m2[6]].data).unit_vec[2];

                                                std::get<Plant::Positive>(rhs[m2[5]].data).unit_vec[0] = par_u[0];
                                                std::get<Plant::Positive>(rhs[m2[5]].data).unit_vec[1] = par_u[1];
                                                std::get<Plant::Positive>(rhs[m2[5]].data).unit_vec[2] = par_u[2];

                                                rhs[m2[5]].position[0] = lhs[m1[2]].position[0] + 0.01 *
                                                                                                  std::get<Plant::Positive>(
                                                                                                          rhs[m2[5]].data).unit_vec[0];
                                                rhs[m2[5]].position[1] = lhs[m1[2]].position[1] + 0.01 *
                                                                                                  std::get<Plant::Positive>(
                                                                                                          rhs[m2[5]].data).unit_vec[1];
                                                rhs[m2[5]].position[2] = lhs[m1[2]].position[2] + 0.01 *
                                                                                                  std::get<Plant::Positive>(
                                                                                                          rhs[m2[5]].data).unit_vec[2];

                                                //DGGML::set_unit_vector(rhs[m2[5]].position, rhs[m2[2]].position, std::get<Plant::Positive>(rhs[m2[5]].data).unit_vec);

                                            });

        gamma.addRule(zippering_case1);

        GT zipper_catastrophe_lhs_graph1;
        zipper_catastrophe_lhs_graph1.addNode({1, {Plant::Intermediate{}}});
        zipper_catastrophe_lhs_graph1.addNode({2, {Plant::Positive{}}});
        zipper_catastrophe_lhs_graph1.addNode({3, {Plant::Zipper{}}});
        zipper_catastrophe_lhs_graph1.addNode({4, {Plant::Intermediate{}}});
        zipper_catastrophe_lhs_graph1.addNode({5, {Plant::Intermediate{}}});
        zipper_catastrophe_lhs_graph1.addEdge(3, 4);
        zipper_catastrophe_lhs_graph1.addEdge(3, 5);
        //zipper_catastrophe_lhs_graph1.addNode({4, {Plant::Intermediate{}}});
        zipper_catastrophe_lhs_graph1.addEdge(1, 2);
        //zipper_catastrophe_lhs_graph1.addEdge(3, 4);

        GT zipper_catastrophe_rhs_graph1;
        zipper_catastrophe_rhs_graph1.addNode({1, {Plant::Intermediate{}}});
        zipper_catastrophe_rhs_graph1.addNode({2, {Plant::Negative{}}});
        zipper_catastrophe_rhs_graph1.addNode({3, {Plant::Zipper{}}});
        zipper_catastrophe_rhs_graph1.addNode({4, {Plant::Intermediate{}}});
        zipper_catastrophe_rhs_graph1.addNode({5, {Plant::Intermediate{}}});
        zipper_catastrophe_rhs_graph1.addEdge(3, 4);
        zipper_catastrophe_rhs_graph1.addEdge(3, 5);
        //zipper_catastrophe_rhs_graph1.addNode({4, {Plant::Intermediate{}}});
        zipper_catastrophe_rhs_graph1.addEdge(1, 2);
        //zipper_catastrophe_rhs_graph1.addEdge(3, 4);

        DGGML::WithRule<GT> zipper_catastrophe("zipper_catastrophe", zipper_catastrophe_lhs_graph1,
                                               zipper_catastrophe_rhs_graph1,
                                               [&](auto &lhs, auto &m) {
                                                   //find all the node data
                                                   auto &dat1 = lhs.findNode(m[1])->second.getData();
                                                   auto &dat2 = lhs.findNode(m[2])->second.getData();
                                                   auto &dat3 = lhs.findNode(m[3])->second.getData();
                                                   //auto& dat4 = lhs.findNode(m[4])->second.getData();

                                                   //get references to position vector
                                                   auto &pos1 = dat1.position;
                                                   auto &pos2 = dat2.position;
                                                   auto &pos3 = dat3.position;
                                                   //auto& pos4 = dat4.position;

                                                   //get references to unit vectors
                                                   auto &u1 = std::get<Plant::Intermediate>(dat1.data).unit_vec;
                                                   auto &u2 = std::get<Plant::Positive>(dat2.data).unit_vec;
                                                   auto &u3 = std::get<Plant::Zipper>(dat3.data).unit_vec;
                                                   //auto& u4 = std::get<Plant::Intermediate>(dat4.data).unit_vec;

                                                   //distance from positive node to line segment
                                                   //auto d = DGGML::distanceToLineSegment(pos3[0], pos3[1], pos4[0], pos4[1], pos2[0], pos2[1]);
                                                   auto d = DGGML::calculate_distance(pos2, pos3);
                                                   if (d <= settings.SEPARATION_DISTANCE) {
                                                       //auto theta = DGGML::compute_theta(u2, dat2.position, u3);
                                                       //if(theta < 45.0)
                                                       //{
                                                       //double sol[2];
                                                       //DGGML::paramaterized_intersection(pos2, pos4, pos3, u2, sol);

                                                       //TODO: determine how strict we need to be with the collision point
                                                       //if(sol[0] > 0.0){// && sol[1] >= 0.0 && sol[1] <= 1.0) {
                                                       return settings.ZIPPERING_GUARD_RATE;
                                                       //}
                                                       //}
                                                   }
                                                   return 0.0;
                                               },
                                               [](auto &lhs, auto &rhs, auto &m1, auto &m2) {
                                                   std::get<Plant::Negative>(
                                                           rhs[m2[2]].data).unit_vec[0] = -std::get<Plant::Positive>(
                                                           lhs[m1[2]].data).unit_vec[0];
                                                   std::get<Plant::Negative>(
                                                           rhs[m2[2]].data).unit_vec[1] = -std::get<Plant::Positive>(
                                                           lhs[m1[2]].data).unit_vec[1];
                                                   std::get<Plant::Negative>(
                                                           rhs[m2[2]].data).unit_vec[2] = -std::get<Plant::Positive>(
                                                           lhs[m1[2]].data).unit_vec[2];
                                               });

        gamma.addRule(zipper_catastrophe);

        //TODO: the retracting end could switch back to grow state, leaving an open funnel point
        // this rule needs more work to guard against that
        //stochastic retraction rule
        GT zipper_retraction_lhs_graph1;
        zipper_retraction_lhs_graph1.addNode({1, {Plant::Negative{}}});
        zipper_retraction_lhs_graph1.addNode({2, {Plant::Intermediate{}}});
        zipper_retraction_lhs_graph1.addNode({3, {Plant::Zipper{}}});
        zipper_retraction_lhs_graph1.addEdge(1, 2);
        zipper_retraction_lhs_graph1.addEdge(2, 3);

        GT zipper_retraction_rhs_graph1;
        zipper_retraction_rhs_graph1.addNode({1, {Plant::Negative{}}});
        zipper_retraction_rhs_graph1.addNode({3, {Plant::Intermediate{}}});
        zipper_retraction_rhs_graph1.addEdge(1, 3);

        DGGML::WithRule<GT> zipper_retraction("zipper_retraction", zipper_retraction_lhs_graph1,
                                              zipper_retraction_rhs_graph1,
                                              [&](auto &lhs, auto &m) {
                                                  auto &node_i_data = lhs.findNode(m[1])->second.getData();
                                                  auto &node_j_data = lhs.findNode(m[2])->second.getData();
                                                  auto len = DGGML::calculate_distance(node_i_data.position,
                                                                                       node_j_data.position);
                                                  double propensity =
                                                          settings.ZIPPERING_RETRACTION_RATE * DGGML::heaviside(settings.DIV_LENGTH_RETRACT, len);
                                                  return propensity;
                                              }, [](auto &lhs, auto &rhs, auto &m1, auto &m2) {
                    //reset the unit vector
                    DGGML::set_unit_vector(rhs[m2[3]].position, rhs[m2[1]].position,
                                           std::get<Plant::Negative>(rhs[m2[1]].data).unit_vec);
                });

        gamma.addRule(zipper_retraction);
    }

    void create_with_clasp_entry_rule(DGGML::Grammar<Plant::graph_type> &gamma, Plant::graph_type &system_graph,
                                      Parameters &settings)
    {
        GT clasp_creation_lhs_graph1;
        clasp_creation_lhs_graph1.addNode({1, {Plant::Boundary{}}});
        clasp_creation_lhs_graph1.addNode({2, {Plant::Boundary{}}});
        clasp_creation_lhs_graph1.addEdge(1, 2);

        GT clasp_creation_rhs_graph1;
        clasp_creation_rhs_graph1.addNode({1, {Plant::Boundary{}}});
        clasp_creation_rhs_graph1.addNode({2, {Plant::Boundary{}}});
        clasp_creation_rhs_graph1.addNode({3, {Plant::Intermediate{}}});
        clasp_creation_rhs_graph1.addNode({4, {Plant::Positive{}}});
        clasp_creation_rhs_graph1.addEdge(1, 3);
        clasp_creation_rhs_graph1.addEdge(3, 2);
        clasp_creation_rhs_graph1.addEdge(3, 4);

        DGGML::WithRule<GT> clasp_creation_case1("clasp_creation_case1", clasp_creation_lhs_graph1,
                                                 clasp_creation_rhs_graph1,
                                                 [&](auto &lhs, auto &m) { return settings.CLASP_ENTRY_RATE; },
                                                 [&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
                                                     auto &lhs_node1 = lhs.findNode(m1[1])->second.getData();
                                                     auto &lhs_node2 = rhs.findNode(m1[2])->second.getData();
                                                     auto &rhs_node3 = rhs.findNode(m2[3])->second.getData();
                                                     auto &rhs_node4 = rhs.findNode(m2[4])->second.getData();

                                                     auto &lhs_pos1 = lhs_node1.position;
                                                     auto &lhs_pos2 = lhs_node2.position;
                                                     auto &rhs_pos3 = rhs_node3.position;
                                                     auto &rhs_pos4 = rhs_node4.position;

                                                     //compute the normal
                                                     auto n_x = lhs_pos2[0] - lhs_pos1[0];
                                                     auto n_y = lhs_pos2[1] - lhs_pos1[1];
                                                     auto n_mag = std::sqrt(n_x * n_x + n_y * n_y);
                                                     n_x /= n_mag;
                                                     n_y /= n_mag;
                                                     //make orthogonal
                                                     std::swap(n_x, n_y);
                                                     n_x = -n_x;
                                                     std::random_device random_device;
                                                     std::mt19937 random_engine(random_device());
                                                     //std::uniform_real_distribution<double>
                                                     //        distribution_local(settings.MT_MIN_SEGMENT_INIT, settings.MT_MAX_SEGMENT_INIT);
                                                     auto radian_angle = settings.CLASP_ENTRY_ANGLE * 3.14 / 180.0;
                                                     std::uniform_real_distribution<double>
                                                             distribution_angle(-radian_angle, radian_angle);

                                                     auto theta = distribution_angle(random_engine);
                                                     auto seg_len = 0.027;//distribution_local(random_engine);
                                                     double x_c, y_c, z_c;
                                                     x_c = (lhs_pos1[0] + lhs_pos2[0]) / 2.0;
                                                     y_c = (lhs_pos1[1] + lhs_pos2[1]) / 2.0;
                                                     z_c = 0;
                                                     auto x_s = n_x * seg_len;//0.0;
                                                     auto y_s = n_y * seg_len;
                                                     auto x_r_t = x_s * cos(theta) + y_s * sin(theta);
                                                     auto y_r_t = -x_s * sin(theta) + y_s * cos(theta);
                                                     auto x_r = x_c + x_r_t;
                                                     auto y_r = y_c + y_r_t;
                                                     auto z_r = 0.0;

                                                     x_r = x_c + x_r_t;//n_x*x_r_t;//seg_len;
                                                     y_r = y_c + y_r_t;//n_y*y_r_t;//seg_len;

                                                     if (boundary_check_2D(settings, x_r, y_r,
                                                                           settings.MAXIMAL_REACTION_RADIUS / 2.0)) {
                                                         x_r = x_c - x_r_t;//n_x*x_r_t;//seg_len;
                                                         y_r = y_c - y_r_t;//n_y*y_r_t;//seg_len;
                                                     }

                                                     //auto x_rot = x_r * cos(theta) + y_r * sin(theta);
                                                     //auto y_rot = -x_r * sin(theta) + y_r * cos(theta);

                                                     //compute dist and unit vector
                                                     rhs_pos3[0] = x_c;
                                                     rhs_pos3[1] = y_c;
                                                     rhs_pos3[2] = z_c;
                                                     rhs_pos4[0] = x_r;
                                                     rhs_pos4[1] = y_r;
                                                     rhs_pos4[2] = z_r;

                                                     //get references to unit vectors for right
                                                     auto &u3 = std::get<Plant::Intermediate>(rhs_node3.data).unit_vec;
                                                     auto &u4 = std::get<Plant::Positive>(rhs_node4.data).unit_vec;
//
                                                     DGGML::set_unit_vector(rhs_pos4, rhs_pos3, u3);
                                                     DGGML::set_unit_vector(rhs_pos4, rhs_pos3, u4);

                                                 });

        gamma.addRule(clasp_creation_case1);
    }

    void create_with_clasp_exit_rule(DGGML::Grammar<Plant::graph_type> &gamma, Plant::graph_type &system_graph,
                                     Parameters &settings)
    {
        GT clasp_boundary_cross_lhs_graph1;
        clasp_boundary_cross_lhs_graph1.addNode({1, {Plant::Intermediate{}}});
        clasp_boundary_cross_lhs_graph1.addNode({2, {Plant::Positive{}}});
        clasp_boundary_cross_lhs_graph1.addNode({3, {Plant::Boundary{}}});
        clasp_boundary_cross_lhs_graph1.addNode({4, {Plant::Boundary{}}});
        clasp_boundary_cross_lhs_graph1.addEdge(1, 2);
        clasp_boundary_cross_lhs_graph1.addEdge(3, 4);

        GT clasp_boundary_cross_rhs_graph1;
        clasp_boundary_cross_rhs_graph1.addNode({1, {Plant::Intermediate{}}});
        clasp_boundary_cross_rhs_graph1.addNode({2, {Plant::Intermediate{}}});
        clasp_boundary_cross_rhs_graph1.addNode({3, {Plant::Boundary{}}});
        clasp_boundary_cross_rhs_graph1.addNode({4, {Plant::Boundary{}}});
        clasp_boundary_cross_rhs_graph1.addEdge(1, 2);
        clasp_boundary_cross_rhs_graph1.addEdge(2, 3);
        clasp_boundary_cross_rhs_graph1.addEdge(2, 4);

        DGGML::WithRule<GT> clasp_boundary_cross_case1("clasp_boundary_cross_case1", clasp_boundary_cross_lhs_graph1,
                                                       clasp_boundary_cross_rhs_graph1,
                                                       [&](auto &lhs, auto &m) {
                                                           //find all the node data
                                                           auto &dat1 = lhs.findNode(m[1])->second.getData();
                                                           auto &dat2 = lhs.findNode(m[2])->second.getData();
                                                           auto &dat3 = lhs.findNode(m[3])->second.getData();
                                                           auto &dat4 = lhs.findNode(m[4])->second.getData();

                                                           //get references to position vector
                                                           auto &pos1 = dat1.position;
                                                           auto &pos2 = dat2.position;
                                                           auto &pos3 = dat3.position;
                                                           auto &pos4 = dat4.position;

                                                           //get references to unit vectors
                                                           auto &u1 = std::get<Plant::Intermediate>(dat1.data).unit_vec;
                                                           auto &u2 = std::get<Plant::Positive>(dat2.data).unit_vec;
                                                           double u3[3];
                                                           DGGML::set_unit_vector(pos3, pos4, u3);

                                                           //distance from positive node to line segment
                                                           auto d = DGGML::distanceToLineSegment(pos3[0], pos3[1],
                                                                                                 pos4[0], pos4[1],
                                                                                                 pos2[0], pos2[1]);
                                                           if (d <= settings.COLLISION_DISTANCE) {
                                                               auto theta = DGGML::compute_theta(u2, dat2.position, u3);
                                                               //std::cout << "Distance, Theta : { " << d << ", " << theta << " }\n";
                                                               if (theta < settings.CLASP_EXIT_ANGLE) {
                                                                   double sol[2];
                                                                   DGGML::paramaterized_intersection(pos2, pos4, pos3,
                                                                                                     u2, sol);
                                                                   //std::cout << "Solution: { " << sol[0] << " " << sol[1] << " }\n";
                                                                   //std::cin.get();
                                                                   //TODO: determine how strict we need to be with the collision point
                                                                   if (sol[0] > 0.0 && sol[1] >= 0.0 &&
                                                                       sol[1] <= 1.0) {
                                                                       //std::cout << "passed\n";
                                                                       //std::cin.get();
                                                                       return settings.CLASP_EXIT_RATE;
                                                                   }
                                                               }
                                                           }
                                                           return 0.0;
                                                       },
                                                       [&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
                                                           auto &lhs_node2 = lhs.findNode(m1[2])->second.getData();
                                                           auto &lhs_node3 = lhs.findNode(m1[3])->second.getData();
                                                           auto &lhs_node4 = rhs.findNode(m1[4])->second.getData();
                                                           auto &rhs_node1 = rhs.findNode(m2[1])->second.getData();
                                                           auto &rhs_node2 = rhs.findNode(m2[2])->second.getData();

                                                           auto &lhs_pos2 = lhs_node2.position;
                                                           auto &lhs_pos3 = lhs_node3.position;
                                                           auto &lhs_pos4 = lhs_node4.position;
                                                           auto &rhs_pos1 = rhs_node1.position;
                                                           auto &rhs_pos2 = rhs_node2.position;

                                                           auto &u3 = std::get<Plant::Intermediate>(
                                                                   rhs_node1.data).unit_vec;
                                                           auto &u4 = std::get<Plant::Intermediate>(
                                                                   rhs_node2.data).unit_vec;

                                                           // place in center
                                                           //for(int i = 0; i < 3; i++)
                                                           //    rhs_pos2[i] = (lhs_pos3[i] + lhs_pos4[i]) / 2.0;

                                                           // place at parameterized intersection point
                                                           auto &u2 = std::get<Plant::Positive>(
                                                                   lhs_node2.data).unit_vec;
                                                           double sol[2];
                                                           DGGML::paramaterized_intersection(lhs_pos2, lhs_pos4,
                                                                                             lhs_pos3, u2, sol);
                                                           for (int i = 0; i < 3; i++)
                                                               rhs_pos2[i] = lhs_pos4[i] +
                                                                             (lhs_pos3[i] - lhs_pos4[i]) * sol[1];

                                                           DGGML::set_unit_vector(rhs_pos2, rhs_pos1, u4);
                                                           DGGML::set_unit_vector(rhs_pos2, rhs_pos1, u3);

//                                                            double sol[2];
//                                                            DGGML::paramaterized_intersection(pos2, pos4, pos3, u2, sol);
//                                                            for(int i = 0; i < 3; i++)
//                                                                rhs[m2[2]].position[i] = lhs[m2[4]].position[i]+(lhs[m2[3]].position[i] - lhs[m2[4]].position[i])*sol[1];
//                                                            DGGML::set_unit_vector(rhs[m2[2]].position, rhs[m2[1]].position, std::get<Plant::Junction>(rhs[m2[2]].data).unit_vec);
//                                                            DGGML::set_unit_vector(rhs[m2[2]].position, rhs[m2[1]].position, std::get<Plant::Intermediate>(rhs[m2[1]].data).unit_vec);
//                                                            for(int i = 0; i < 3; i++)
//                                                                rhs[m2[5]].position[i] = rhs[m1[2]].position[i]+0.01*std::get<Plant::Junction>(rhs[m1[2]].data).unit_vec[i];
//                                                            DGGML::set_unit_vector(rhs[m2[5]].position, rhs[m2[2]].position, std::get<Plant::Intermediate>(rhs[m2[5]].data).unit_vec);
//                                                            for(int i = 0; i < 3; i++)
//                                                                rhs[m2[6]].position[i] = rhs[m1[2]].position[i]+0.011*std::get<Plant::Junction>(rhs[m1[2]].data).unit_vec[i];
//                                                            DGGML::set_unit_vector(rhs[m2[6]].position, rhs[m2[2]].position, std::get<Plant::Positive>(rhs[m2[6]].data).unit_vec);


                                                       });
        gamma.addRule(clasp_boundary_cross_case1);
    }


    void create_with_clasp_catastrophe(DGGML::Grammar<Plant::graph_type> &gamma, Plant::graph_type &system_graph,
                            Parameters &settings) {

        GT clasp_catastrophe_lhs_graph1;
        clasp_catastrophe_lhs_graph1.addNode({1, {Plant::Boundary{}}});
        clasp_catastrophe_lhs_graph1.addNode({2, {Plant::Boundary{}}});
        clasp_catastrophe_lhs_graph1.addNode({3, {Plant::Intermediate{}}});
        clasp_catastrophe_lhs_graph1.addNode({4, {Plant::Negative{}}});
        clasp_catastrophe_lhs_graph1.addEdge(1, 3);
        clasp_catastrophe_lhs_graph1.addEdge(3, 2);
        clasp_catastrophe_lhs_graph1.addEdge(3, 4);

        GT clasp_catastrophe_rhs_graph1;
        clasp_catastrophe_rhs_graph1.addNode({1, {Plant::Boundary{}}});
        clasp_catastrophe_rhs_graph1.addNode({2, {Plant::Boundary{}}});
        clasp_catastrophe_rhs_graph1.addEdge(1, 2);

        DGGML::WithRule<GT> clasp_catastrophe_case1("clasp_catastrophe_case1", clasp_catastrophe_lhs_graph1,
                                                    clasp_catastrophe_rhs_graph1,
                                                    [&](auto &lhs, auto &m) {
                                                        auto &dat3 = lhs.findNode(m[3])->second.getData();
                                                        auto &dat4 = lhs.findNode(m[4])->second.getData();
                                                        auto &pos3 = dat3.position;
                                                        auto &pos4 = dat4.position;

                                                        auto d = DGGML::calculate_distance(pos3, pos4);
                                                        auto d2 = sqrt((pos3[0] - pos4[0]) * (pos3[0] - pos4[0]) +
                                                                       (pos3[1] - pos4[1]) * (pos3[1] - pos4[1]));
                                                        //std::cin.get(); std::cout << d << " " << d2 << "\n";
                                                        auto &t3 = std::get<Plant::Intermediate>(dat3.data);
                                                        auto &t4 = std::get<Plant::Negative>(dat4.data);
                                                        //std::cin.get();
                                                        if (d <= 0.01) {
                                                            //std::cout << "here\n";
                                                            //std::cin.get();
                                                            return settings.CLASP_CAT_RATE;
                                                        } else
                                                            return 0.0;
                                                    },
                                                    [&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
                                                    });
        gamma.addRule(clasp_catastrophe_case1);
    }

    void create_with_clasp_detachment(DGGML::Grammar<Plant::graph_type> &gamma, Plant::graph_type &system_graph,
                                       Parameters &settings) {

        GT clasp_detachment_lhs_graph1;
        clasp_detachment_lhs_graph1.addNode({1, {Plant::Boundary{}}});
        clasp_detachment_lhs_graph1.addNode({2, {Plant::Boundary{}}});
        clasp_detachment_lhs_graph1.addNode({3, {Plant::Intermediate{}}});
        clasp_detachment_lhs_graph1.addNode({4, {Plant::Intermediate{}}});
        clasp_detachment_lhs_graph1.addEdge(1, 3);
        clasp_detachment_lhs_graph1.addEdge(3, 2);
        clasp_detachment_lhs_graph1.addEdge(3, 4);

        GT clasp_detachment_rhs_graph1;
        clasp_detachment_rhs_graph1.addNode({1, {Plant::Boundary{}}});
        clasp_detachment_rhs_graph1.addNode({2, {Plant::Boundary{}}});
        clasp_detachment_rhs_graph1.addNode({3, {Plant::Negative{}}});
        clasp_detachment_rhs_graph1.addNode({4, {Plant::Intermediate{}}});
        clasp_detachment_rhs_graph1.addEdge(1, 2);
        clasp_detachment_rhs_graph1.addEdge(3, 4);

        DGGML::WithRule<GT> clasp_detachment_case1("clasp_detachment_case1", clasp_detachment_lhs_graph1,
                                                    clasp_detachment_rhs_graph1,
                                                    [&](auto &lhs, auto &m) {
                                                        return settings.CLASP_DETACHMENT_RATE;
                                                    },
                                                    [&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
                                                        auto &dat3_lhs = lhs.findNode(m1[3])->second.getData();
                                                        auto &dat4_lhs = lhs.findNode(m1[4])->second.getData();
                                                        auto &dat3_rhs = rhs.findNode(m2[3])->second.getData();

                                                        //get references to position vector
                                                        auto &pos3_lhs = dat3_lhs.position;
                                                        auto &pos4_lhs = dat4_lhs.position;
                                                        auto &pos3_rhs = dat3_rhs.position;

                                                        pos3_rhs[0] = pos3_lhs[0] - (pos3_lhs[0] - pos4_lhs[0]) / 20.0;
                                                        pos3_rhs[1] = pos3_lhs[1] - (pos3_lhs[1] - pos4_lhs[1]) / 20.0;
                                                        pos3_rhs[2] = pos3_lhs[2] - (pos3_lhs[2] - pos4_lhs[2]) / 20.0;

                                                        //get references to unit vectors
                                                        auto &u3 = std::get<Plant::Negative>(dat3_rhs.data).unit_vec;
                                                        DGGML::set_unit_vector(pos4_lhs, pos3_rhs, u3);
                                                    });
        gamma.addRule(clasp_detachment_case1);
    }

    void create_with_destruction_rule(DGGML::Grammar<Plant::graph_type> &gamma, Plant::graph_type &system_graph,
                                  Parameters &settings) {
        //destruction rule case 1: depolymerization on both ends
        GT destruction_lhs_graph1;
        destruction_lhs_graph1.addNode({1, {Plant::Negative{}}});
        destruction_lhs_graph1.addNode({2, {Plant::Intermediate{}}});
        destruction_lhs_graph1.addNode({3, {Plant::Negative{}}});
        destruction_lhs_graph1.addEdge(1, 2);
        destruction_lhs_graph1.addEdge(2, 3);

        GT destruction_rhs_graph1;

        DGGML::WithRule<GT> destruction_case1("destruction_case1", destruction_lhs_graph1, destruction_rhs_graph1,
                                              [&](auto &lhs, auto &m) { return settings.MT_DESTRUCTION_RATE; },
                                              [](auto &lhs, auto &rhs, auto &m1, auto &m2) {});

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
                                              [](auto &lhs, auto &m) { return 1; },
                                              [](auto &lhs, auto &rhs, auto &m1, auto &m2) {});

        //gamma.addRule(destruction_case2);
    }

    void create_with_creation_rule(DGGML::Grammar<Plant::graph_type> &gamma, Plant::graph_type &system_graph,
                               Parameters &settings) {
        GT creation_lhs_graph1;
        creation_lhs_graph1.addNode({1, {Plant::Nucleator{}}});

        GT creation_rhs_graph1;
        creation_rhs_graph1.addNode({1, {Plant::Nucleator{}}});
        creation_rhs_graph1.addNode({2, {Plant::Negative{}}});
        creation_rhs_graph1.addNode({3, {Plant::Intermediate{}}});
        creation_rhs_graph1.addNode({4, {Plant::Positive{}}});
        creation_rhs_graph1.addEdge(2, 3);
        creation_rhs_graph1.addEdge(3, 4);

        DGGML::WithRule<GT> creation_case1("creation_case1", creation_lhs_graph1, creation_rhs_graph1,
                                           [&](auto &lhs, auto &m) { return settings.CREATION_FACTOR*settings.CREATION_RATE; },
                                           [&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
                                               //std::cin.get();
                                               //find all the node data for the left
                                               auto &lhs_node1 = lhs.findNode(m1[1])->second.getData();
                                               //find all the node data for the right
                                               auto &rhs_node2 = rhs.findNode(m2[2])->second.getData();
                                               auto &rhs_node3 = rhs.findNode(m2[3])->second.getData();
                                               auto &rhs_node4 = rhs.findNode(m2[4])->second.getData();

                                               //get references to position vectors for left
                                               auto &lhs_pos1 = lhs_node1.position;
                                               //get references to position vectors for right
                                               auto &rhs_pos2 = rhs_node2.position;
                                               auto &rhs_pos3 = rhs_node3.position;
                                               auto &rhs_pos4 = rhs_node4.position;

                                               std::random_device random_device;
                                               std::mt19937 random_engine(random_device());
                                               std::uniform_real_distribution<double>
                                                       distribution_local(settings.MT_MIN_SEGMENT_INIT,
                                                                          settings.MT_MAX_SEGMENT_INIT);
                                               std::uniform_real_distribution<double>
                                                       distribution_angle(0.0, 2.0 * 3.14);

                                               auto theta = distribution_angle(random_engine);
                                               auto seg_len = distribution_local(random_engine);
                                               double x_c, y_c, z_c;
                                               x_c = lhs_pos1[0];
                                               y_c = lhs_pos1[1];
                                               z_c = 0;
                                               auto x_s = 0.0;
                                               auto y_s = seg_len;
                                               auto x_r_t = x_s * cos(theta) + y_s * sin(theta);
                                               auto y_r_t = -x_s * sin(theta) + y_s * cos(theta);
                                               auto x_r = x_c + x_r_t;
                                               auto y_r = y_c + y_r_t;
                                               auto z_r = 0.0;
                                               auto x_l = x_c - (x_r - x_c);
                                               auto y_l = y_c - (y_r - y_c);
                                               auto z_l = 0.0;

                                               //compute dist and unit vector
                                               rhs_pos2[0] = x_l;
                                               rhs_pos2[1] = y_l;
                                               rhs_pos2[2] = z_l;
                                               rhs_pos3[0] = x_c;
                                               rhs_pos3[1] = y_c;
                                               rhs_pos3[2] = z_c;
                                               rhs_pos4[0] = x_r;
                                               rhs_pos4[1] = y_r;
                                               rhs_pos4[2] = z_r;

                                               //get references to unit vectors for right
                                               auto &u2 = std::get<Plant::Negative>(rhs_node2.data).unit_vec;
                                               auto &u3 = std::get<Plant::Intermediate>(rhs_node3.data).unit_vec;
                                               auto &u4 = std::get<Plant::Positive>(rhs_node4.data).unit_vec;
//                                                      for(int i = 0; i < 3; i++) {
//                                                          rhs_pos2[i] = lhs_pos1[i]+settings.MT_MIN_SEGMENT_INIT;
//                                                          rhs_pos3[i] = lhs_pos1[i];
//                                                          rhs_pos4[i] = lhs_pos1[i]-+settings.MT_MIN_SEGMENT_INIT;
//                                                          if(i == 2)
//                                                          {
//                                                              rhs_pos2[i] = 0.0;
//                                                              rhs_pos3[i] = 0.0;
//                                                              rhs_pos4[i] = 0.0;
//                                                          }
//                                                      }
                                               DGGML::set_unit_vector(rhs_pos4, rhs_pos3, u2);
                                               DGGML::set_unit_vector(rhs_pos4, rhs_pos3, u3);
                                               DGGML::set_unit_vector(rhs_pos4, rhs_pos3, u4);

                                           });

        gamma.addRule(creation_case1);
    }


    void create_with_recovery_rule(DGGML::Grammar<Plant::graph_type> &gamma, Plant::graph_type &system_graph,
                               Parameters &settings) {
        //TODO: bug, flipping the ordering of the lhs to {1: negative, 2: intermediate}
        // and rhs to {1: positive, 2: intermediate} causes a bad variant access
        GT rescue_lhs;
        rescue_lhs.addNode({1, {Plant::Intermediate{}}});
        rescue_lhs.addNode({2, {Plant::Negative{}}});
        rescue_lhs.addEdge(1, 2);

        GT rescue_rhs;
        rescue_rhs.addNode({1, {Plant::Intermediate{}}});
        rescue_rhs.addNode({2, {Plant::Positive{}}});
        rescue_rhs.addEdge(1, 2);

        DGGML::WithRule<GT> rescue("rescue", rescue_lhs, rescue_rhs,
                                   [&](auto &lhs, auto &m) {
                                       return settings.RECOVERY_FACTOR*settings.RECOVERY_RATE;
                                   }, [](auto &lhs, auto &rhs, auto &m1, auto &m2) {
                    std::get<Plant::Positive>(rhs[m2[2]].data).unit_vec[0] = -std::get<Plant::Negative>(
                            lhs[m1[2]].data).unit_vec[0];
                    std::get<Plant::Positive>(rhs[m2[2]].data).unit_vec[1] = -std::get<Plant::Negative>(
                            lhs[m1[2]].data).unit_vec[1];
                    std::get<Plant::Positive>(rhs[m2[2]].data).unit_vec[2] = -std::get<Plant::Negative>(
                            lhs[m1[2]].data).unit_vec[2];

                });

        gamma.addRule(rescue);
    }

    void experimental_rules(DGGML::Grammar<Plant::graph_type> &gamma, Plant::graph_type &system_graph,
                            Parameters &settings) {
        GT periodic_replicate_lhs;
        periodic_replicate_lhs.addNode({1, {Plant::Intermediate{}}});
        periodic_replicate_lhs.addNode({2, {Plant::Capture{}}});
        periodic_replicate_lhs.addEdge(1, 2);
        GT periodic_replicate_rhs;
        periodic_replicate_rhs.addNode({1, {Plant::Intermediate{}}});
        periodic_replicate_rhs.addNode({2, {Plant::Holding{}}});
        periodic_replicate_rhs.addNode({3, {Plant::Intermediate{}}});
        periodic_replicate_rhs.addNode({4, {Plant::Positive{}}});
        periodic_replicate_rhs.addEdge(1, 2);
        periodic_replicate_rhs.addEdge(3, 4);

        //TODO: replicated MTs show up on the other side and seem not to be added back into the match data struct
        DGGML::WithRule<GT> periodic_replicate("periodic_replicate", periodic_replicate_lhs, periodic_replicate_rhs,
                                               [&](auto &lhs, auto &m) {
                                                   double len_x = settings.CELL_DX * settings.CELL_NX;
                                                   double len_y = settings.CELL_DY * settings.CELL_NY;
                                                   double r = settings.MAXIMAL_REACTION_RADIUS;
                                                   double period_x = len_x - r;
                                                   double period_y = len_y - r;
                                                   if (lhs[m[2]].position[0] >= len_x - r) //went out the right side
                                                   {
                                                       return 1000.0;
                                                   } else if (lhs[m[2]].position[0] <= r) {
                                                       return 1000.0;
                                                   }
                                                   return 0.0;
                                               },
                                               [&](auto &lhs, auto &rhs, auto &m1, auto &m2) {
                                                   double len_x = settings.CELL_DX * settings.CELL_NX;
                                                   double len_y = settings.CELL_DY * settings.CELL_NY;
                                                   double r = settings.MAXIMAL_REACTION_RADIUS;
                                                   double period_x = len_x - r;
                                                   double period_y = len_y - r;
                                                   if (lhs[m1[2]].position[0] >= len_x - r) //went out the right side
                                                   {
                                                       //set the intermediate
                                                       rhs[m2[3]].position[0] =
                                                               lhs[m1[1]].position[0] - len_x + 2.0 * r;
                                                       rhs[m2[3]].position[1] = lhs[m1[1]].position[1];
                                                       rhs[m2[3]].position[2] = 0.0;
                                                       //set the positive
                                                       rhs[m2[4]].position[0] =
                                                               lhs[m1[2]].position[0] - len_x + 2.1 * r;
                                                       rhs[m2[4]].position[1] = lhs[m1[2]].position[1];
                                                       rhs[m2[4]].position[2] = 0.0;
                                                       std::get<Plant::Positive>(
                                                               rhs[m2[4]].data).unit_vec[0] = std::get<Plant::Capture>(
                                                               lhs[m1[2]].data).unit_vec[0];
                                                       std::get<Plant::Positive>(
                                                               rhs[m2[4]].data).unit_vec[1] = std::get<Plant::Capture>(
                                                               lhs[m1[2]].data).unit_vec[1];
                                                       std::get<Plant::Positive>(
                                                               rhs[m2[4]].data).unit_vec[2] = std::get<Plant::Capture>(
                                                               lhs[m1[2]].data).unit_vec[2];
                                                   } else if (lhs[m1[2]].position[0] <= r) //went out the left side
                                                   {
                                                       //set the intermediate
                                                       rhs[m2[3]].position[0] =
                                                               lhs[m1[1]].position[0] + len_x - 2.0 * r;
                                                       rhs[m2[3]].position[1] = lhs[m1[1]].position[1];
                                                       rhs[m2[3]].position[2] = 0.0;
                                                       //set the positive
                                                       rhs[m2[4]].position[0] =
                                                               lhs[m1[2]].position[0] + len_x - 2.1 * r;
                                                       rhs[m2[4]].position[1] = lhs[m1[2]].position[1];
                                                       rhs[m2[4]].position[2] = 0.0;
                                                       std::get<Plant::Positive>(
                                                               rhs[m2[4]].data).unit_vec[0] = -std::get<Plant::Capture>(
                                                               lhs[m1[2]].data).unit_vec[0];
                                                       std::get<Plant::Positive>(
                                                               rhs[m2[4]].data).unit_vec[1] = -std::get<Plant::Capture>(
                                                               lhs[m1[2]].data).unit_vec[1];
                                                       std::get<Plant::Positive>(
                                                               rhs[m2[4]].data).unit_vec[2] = -std::get<Plant::Capture>(
                                                               lhs[m1[2]].data).unit_vec[2];
                                                   }
                                               });
        //gamma.addRule(periodic_replicate);
    }

    // old boundary caputre rule, that puts the growing end into a pause state of sorts
    void create_boundary_rules(DGGML::Grammar<Plant::graph_type> &gamma, Plant::graph_type &system_graph,
                               Parameters &settings) {

        GT boundary_capture_lhs;
        boundary_capture_lhs.addNode({1, {Plant::Intermediate{}}});
        boundary_capture_lhs.addNode({2, {Plant::Positive{}}});
        boundary_capture_lhs.addEdge(1, 2);
        GT boundary_capture_rhs;
        boundary_capture_rhs.addNode({1, {Plant::Intermediate{}}});
        boundary_capture_rhs.addNode({2, {Plant::Capture{}}});
        boundary_capture_rhs.addEdge(1, 2);

        DGGML::WithRule<GT> boundary_capture("boundary_capture", boundary_capture_lhs, boundary_capture_rhs,
                                             [&](auto &lhs, auto &m) {
                                                 //determine if we're outside the interior of a padded boundary
                                                 if (boundary_check_2D(settings, lhs[m[2]].position[0],
                                                                       lhs[m[2]].position[1],
                                                                       settings.MAXIMAL_REACTION_RADIUS))
                                                     return 1000.0;
                                                 //if(boundary_check_circle(settings, lhs[m[2]].position[0], lhs[m[2]].position[1])) return 1000.0;
                                                 return 0.0;
                                             },
                                             [](auto &lhs, auto &rhs, auto &m1, auto &m2) {
                                                 std::get<Plant::Capture>(
                                                         rhs[m2[2]].data).unit_vec[0] = std::get<Plant::Positive>(
                                                         lhs[m1[2]].data).unit_vec[0];
                                                 std::get<Plant::Capture>(
                                                         rhs[m2[2]].data).unit_vec[1] = std::get<Plant::Positive>(
                                                         lhs[m1[2]].data).unit_vec[1];
                                                 std::get<Plant::Capture>(
                                                         rhs[m2[2]].data).unit_vec[2] = std::get<Plant::Positive>(
                                                         lhs[m1[2]].data).unit_vec[2];

                                                 //std::get<Plant::Negative>(rhs[m2[2]].data).unit_vec[0] = -std::get<Plant::Positive>(lhs[m1[2]].data).unit_vec[0];
                                                 //std::get<Plant::Negative>(rhs[m2[2]].data).unit_vec[1] = -std::get<Plant::Positive>(lhs[m1[2]].data).unit_vec[1];
                                                 //std::get<Plant::Negative>(rhs[m2[2]].data).unit_vec[2] = -std::get<Plant::Positive>(lhs[m1[2]].data).unit_vec[2];
                                             });
        //gamma.addRule(boundary_capture);
    }

    // old collision rules
    void create_collision_rules(DGGML::Grammar<Plant::graph_type> &gamma, Plant::graph_type &system_graph,
                                Parameters &settings) {
//            //TODO: reduce the copy paste coding
//            //collision catastrophe rule case 1
//            GT catastrophe_lhs_graph1;
//            catastrophe_lhs_graph1.addNode({1, {Plant::Intermediate{}}});
//            catastrophe_lhs_graph1.addNode({2, {Plant::Positive{}}});
//            catastrophe_lhs_graph1.addNode({3, {Plant::Intermediate{}}});
//            catastrophe_lhs_graph1.addNode({4, {Plant::Intermediate{}}});
//            catastrophe_lhs_graph1.addEdge(1, 2);
//            catastrophe_lhs_graph1.addEdge(3, 4);
//
//            GT catastrophe_rhs_graph1;
//            catastrophe_rhs_graph1.addNode({1, {Plant::Intermediate{}}});
//            catastrophe_rhs_graph1.addNode({2, {Plant::Negative{}}});
//            catastrophe_rhs_graph1.addNode({3, {Plant::Intermediate{}}});
//            catastrophe_rhs_graph1.addNode({4, {Plant::Intermediate{}}});
//            catastrophe_rhs_graph1.addEdge(1, 2);
//            catastrophe_rhs_graph1.addEdge(3, 4);
//
//            DGGML::WithRule<GT> catastrophe_case1("catastrophe_case1", catastrophe_lhs_graph1, catastrophe_rhs_graph1,
//                                                  [&](auto& lhs, auto& m) {
//                                                      //find all the node data
//                                                      auto& dat1 = lhs.findNode(m[1])->second.getData();
//                                                      auto& dat2 = lhs.findNode(m[2])->second.getData();
//                                                      auto& dat3 = lhs.findNode(m[3])->second.getData();
//                                                      auto& dat4 = lhs.findNode(m[4])->second.getData();
//
//                                                      //get references to position vector
//                                                      auto& pos1 = dat1.position;
//                                                      auto& pos2 = dat2.position;
//                                                      auto& pos3 = dat3.position;
//                                                      auto& pos4 = dat4.position;
//
//                                                      //get references to unit vectors
//                                                      auto& u1 = std::get<Plant::Intermediate>(dat1.data).unit_vec;
//                                                      auto& u2 = std::get<Plant::Positive>(dat2.data).unit_vec;
//                                                      auto& u3 = std::get<Plant::Intermediate>(dat3.data).unit_vec;
//                                                      auto& u4 = std::get<Plant::Intermediate>(dat4.data).unit_vec;
//
//                                                      double theta = DGGML::unit_dot_product(u2, u3);
//                                                      double propensity = 0.0;
//                                                      double sol[2];
//                                                      if(theta != 1.0 && theta != -1.0)
//                                                      {
//                                                          DGGML::paramaterized_intersection(pos2, pos4, pos3, u2, sol);
//
//                                                          if(sol[0] > 0.0 && sol[1] >= 0.0 && sol[1] <= 1.0)
//                                                          {
//                                                              //TODO: distance should be closest point to the growing end not pos3
//                                                              propensity = 1000.0*exp(-pow(DGGML::calculate_distance(pos2, pos3), 2.0) / pow(0.5*settings.DIV_LENGTH, 2.0));
//                                                          }
//                                                      }
//                                                      //if(propensity != 0) {std::cout << "cat prop: " << propensity << "\n"; std::cin.get();}
//                                                      return propensity;
//                                                  },
//                                                  [](auto& lhs, auto& rhs, auto& m1, auto& m2)
//                                                  {
//                                                      std::get<Plant::Negative>(rhs[m2[2]].data).unit_vec[0] = -std::get<Plant::Positive>(lhs[m1[2]].data).unit_vec[0];
//                                                      std::get<Plant::Negative>(rhs[m2[2]].data).unit_vec[1] = -std::get<Plant::Positive>(lhs[m1[2]].data).unit_vec[1];
//                                                      std::get<Plant::Negative>(rhs[m2[2]].data).unit_vec[2] = -std::get<Plant::Positive>(lhs[m1[2]].data).unit_vec[2];
//                                                  });
//
//            //gamma.addRule(catastrophe_case1);
//
//            //collision catastrophe rule case 2
//            GT catastrophe_lhs_graph2;
//            catastrophe_lhs_graph2.addNode({1, {Plant::Intermediate{}}});
//            catastrophe_lhs_graph2.addNode({2, {Plant::Positive{}}});
//            catastrophe_lhs_graph2.addNode({3, {Plant::Intermediate{}}});
//            catastrophe_lhs_graph2.addNode({4, {Plant::Positive{}}});
//            catastrophe_lhs_graph2.addEdge(1, 2);
//            catastrophe_lhs_graph2.addEdge(3, 4);
//
//            GT catastrophe_rhs_graph2;
//            catastrophe_rhs_graph2.addNode({1, {Plant::Intermediate{}}});
//            catastrophe_rhs_graph2.addNode({2, {Plant::Negative{}}});
//            catastrophe_rhs_graph2.addNode({3, {Plant::Intermediate{}}});
//            catastrophe_rhs_graph2.addNode({4, {Plant::Positive{}}});
//            catastrophe_rhs_graph2.addEdge(1, 2);
//            catastrophe_rhs_graph2.addEdge(3, 4);
//
//            DGGML::WithRule<GT> catastrophe_case2("catastrophe_case2", catastrophe_lhs_graph2, catastrophe_rhs_graph2,
//                                                  [&](auto& lhs, auto& m) {
//                                                      //find all the node data
//                                                      auto& dat1 = lhs.findNode(m[1])->second.getData();
//                                                      auto& dat2 = lhs.findNode(m[2])->second.getData();
//                                                      auto& dat3 = lhs.findNode(m[3])->second.getData();
//                                                      auto& dat4 = lhs.findNode(m[4])->second.getData();
//
//                                                      //get references to position vector
//                                                      auto& pos1 = dat1.position;
//                                                      auto& pos2 = dat2.position;
//                                                      auto& pos3 = dat3.position;
//                                                      auto& pos4 = dat4.position;
//
//                                                      //get references to unit vectors
//                                                      auto& u1 = std::get<Plant::Intermediate>(dat1.data).unit_vec;
//                                                      auto& u2 = std::get<Plant::Positive>(dat2.data).unit_vec;
//                                                      auto& u3 = std::get<Plant::Intermediate>(dat3.data).unit_vec;
//                                                      auto& u4 = std::get<Plant::Positive>(dat4.data).unit_vec;
//
//                                                      double theta = DGGML::unit_dot_product(u2, u3);
//                                                      double propensity = 0.0;
//                                                      double sol[2];
//                                                      if(theta != 1.0 && theta != -1.0)
//                                                      {
//                                                          DGGML::paramaterized_intersection(pos2, pos4, pos3, u2, sol);
//
//                                                          if(sol[0] > 0.0 && sol[1] >= 0.0 && sol[1] <= 1.0)
//                                                          {
//                                                              //TODO: distance should be closest point to the growing end not pos3
//                                                              propensity = 1000.0*exp(-pow(DGGML::calculate_distance(pos2, pos3), 2.0) / pow(0.5*settings.DIV_LENGTH, 2.0));
//                                                          }
//                                                      }
//                                                      //if(propensity != 0) {std::cout << "cat prop: " << propensity << "\n"; std::cin.get();}
//                                                      return propensity;
//                                                  },
//                                                  [](auto& lhs, auto& rhs, auto& m1, auto& m2)
//                                                  {
//                                                      std::get<Plant::Negative>(rhs[m2[2]].data).unit_vec[0] = -std::get<Plant::Positive>(lhs[m1[2]].data).unit_vec[0];
//                                                      std::get<Plant::Negative>(rhs[m2[2]].data).unit_vec[1] = -std::get<Plant::Positive>(lhs[m1[2]].data).unit_vec[1];
//                                                      std::get<Plant::Negative>(rhs[m2[2]].data).unit_vec[2] = -std::get<Plant::Positive>(lhs[m1[2]].data).unit_vec[2];
//                                                  });
//
//            //gamma.addRule(catastrophe_case2);
//
//            //collision catastrophe rule case 3
//            GT catastrophe_lhs_graph3;
//            catastrophe_lhs_graph3.addNode({1, {Plant::Intermediate{}}});
//            catastrophe_lhs_graph3.addNode({2, {Plant::Positive{}}});
//            catastrophe_lhs_graph3.addNode({3, {Plant::Intermediate{}}});
//            catastrophe_lhs_graph3.addNode({4, {Plant::Negative{}}});
//            catastrophe_lhs_graph3.addEdge(1, 2);
//            catastrophe_lhs_graph3.addEdge(3, 4);
//
//            //TODO: FIX ME! The rhs of this rule is messing with the catastrophe solving rule
//            GT catastrophe_rhs_graph3;
//            catastrophe_rhs_graph3.addNode({1, {Plant::Intermediate{}}});
//            catastrophe_rhs_graph3.addNode({2, {Plant::Negative{}}});
//            catastrophe_rhs_graph3.addNode({3, {Plant::Intermediate{}}});
//            catastrophe_rhs_graph3.addNode({4, {Plant::Negative{}}});
//            catastrophe_rhs_graph3.addEdge(1, 2);
//            catastrophe_rhs_graph3.addEdge(3, 4);
//
//            DGGML::WithRule<GT> catastrophe_case3("catastrophe_case3", catastrophe_lhs_graph3, catastrophe_rhs_graph3,
//                                                  [&](auto& lhs, auto& m) {
//
//                                                      //find all the node data
//                                                      auto& dat1 = lhs.findNode(m[1])->second.getData();
//                                                      auto& dat2 = lhs.findNode(m[2])->second.getData();
//                                                      auto& dat3 = lhs.findNode(m[3])->second.getData();
//                                                      auto& dat4 = lhs.findNode(m[4])->second.getData();
//
//                                                      //get references to position vector
//                                                      auto& pos1 = dat1.position;
//                                                      auto& pos2 = dat2.position;
//                                                      auto& pos3 = dat3.position;
//                                                      auto& pos4 = dat4.position;
//
//                                                      //get references to unit vectors
//                                                      auto& u1 = std::get<Plant::Intermediate>(dat1.data).unit_vec;
//                                                      auto& u2 = std::get<Plant::Positive>(dat2.data).unit_vec;
//                                                      auto& u3 = std::get<Plant::Intermediate>(dat3.data).unit_vec;
//                                                      auto& u4 = std::get<Plant::Negative>(dat4.data).unit_vec;
//
//                                                      double theta = DGGML::unit_dot_product(u2, u3);
//                                                      double propensity = 0.0;
//                                                      double sol[2];
//                                                      if(theta != 1.0 && theta != -1.0)
//                                                      {
//                                                          DGGML::paramaterized_intersection(pos2, pos4, pos3, u2, sol);
//
//                                                          if(sol[0] > 0.0 && sol[1] >= 0.0 && sol[1] <= 1.0)
//                                                          {
//                                                              //TODO: distance should be closest point to the growing end not pos3
//                                                              propensity = 1000.0*exp(-pow(DGGML::calculate_distance(pos2, pos3), 2.0) / pow(0.5*settings.DIV_LENGTH, 2.0));
//                                                          }
//                                                      }
//
//                                                      //if(propensity != 0) {std::cout << "cat prop: " << propensity << "\n"; std::cin.get();}
//                                                      return propensity;
//                                                  },
//                                                  [](auto& lhs, auto& rhs, auto& m1, auto& m2)
//                                                  {
//
//                                                      std::get<Plant::Negative>(rhs[m2[2]].data).unit_vec[0] = -std::get<Plant::Positive>(lhs[m1[2]].data).unit_vec[0];
//                                                      std::get<Plant::Negative>(rhs[m2[2]].data).unit_vec[1] = -std::get<Plant::Positive>(lhs[m1[2]].data).unit_vec[1];
//                                                      std::get<Plant::Negative>(rhs[m2[2]].data).unit_vec[2] = -std::get<Plant::Positive>(lhs[m1[2]].data).unit_vec[2];
//
//                                                  });
//
//            //gamma.addRule(catastrophe_case3);
    }
}

#endif //DGGML_CMARULES_HPP
