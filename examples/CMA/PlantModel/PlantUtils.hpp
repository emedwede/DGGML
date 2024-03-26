#ifndef PLANT_UTILS_HPP
#define PLANT_UTILS_HPP

#include <random>
#include <algorithm>
#include <vector>

#include "PlantTypes.hpp"
#include "Utlities/MathUtils.hpp"

#include "CartesianGrid2D.hpp"

namespace DGGML
{
namespace Plant 
{
    template <typename GraphType, typename CplexType, typename ParamType>
    void microtubule_unit_scatter(GraphType& graph, CplexType& cplex, ParamType& settings)
    {
        double epsilon_min = settings.MT_MIN_SEGMENT_INIT;
        double epsilon_max = settings.MT_MAX_SEGMENT_INIT;
        std::random_device random_device;
        std::mt19937 random_engine(random_device());
        auto& grid = cplex.getCoarseGrid();
        
        std::uniform_real_distribution<double> 
            distribution_global_x(cplex.min_x+epsilon_max, 
                    cplex.max_x-epsilon_max);

        std::uniform_real_distribution<double>
            distribution_global_y(cplex.min_y+epsilon_max, 
                    cplex.max_y-epsilon_max);

        std::uniform_real_distribution<double> 
            distribution_local(epsilon_min, epsilon_max);
        
        std::uniform_real_distribution<double>
            distribution_angle(0.0, 2.0*3.14);
        using node_type = typename GraphType::node_type;

        std::size_t segments = 3;
        for(auto i = 0; i < settings.NUM_MT; i++) 
        {
            auto x_c = distribution_global_x(random_engine);
            auto y_c = distribution_global_y(random_engine);
            auto z_c = 0.0; 

            auto theta = distribution_angle(random_engine);
            auto seg_len = distribution_local(random_engine);

            auto x_s = 0.0;
            auto y_s = seg_len;
            auto x_r_t = x_s*cos(theta) + y_s*sin(theta);
            auto y_r_t = -x_s*sin(theta) +y_s*cos(theta);
            auto x_r = x_c + x_r_t;
            auto y_r = y_c + y_r_t;
            auto z_r = 0.0;

            auto x_l = x_c - (x_r - x_c);
            auto y_l = y_c - (y_r - y_c);
            auto z_l = 0.0;
            
            //compute dist and unit vector
            double p1[3] = {x_r - x_c, y_r - y_c, 0.0};
            double p2[3] = {0.0, 0.0, 0.0};
            double p3[3] = {x_c - x_l, y_c - y_l, 0.0};
            double u1[3], u2[3];

            set_unit_vector(p1, p2, u1);
            set_unit_vector(p3, p2, u2);
            node_type node_l(i*segments, 
                    {{x_l, y_l, z_l}, 
                    {0.0, 0.0, 0.0}, 
                    negative, 
                    {-1, -1, -1}, 
                    {u2[0], u2[1], u2[2]}});

            node_type node_c(i*segments+1, 
                    {{x_c, y_c, z_c}, 
                    {0.0, 0.0, 0.0}, 
                    intermediate, 
                    {-1, -1, -1}, 
                    {u1[0], u1[1], u1[2]}});
            
            node_type node_r(i*segments+2, 
                    {{x_r, y_r, z_r}, 
                    {0.0, 0.0, 0.0}, 
                    positive, 
                    {-1, -1, -1}, 
                    {u1[0], u1[1], u1[2]}});

            graph.addNode(node_l);
            graph.addNode(node_c);
            graph.addNode(node_r);

            graph.addEdge(node_l, node_c);
            graph.addEdge(node_r, node_c);
        }

    }
    
    template <typename GraphType, typename CplexType, typename ParamType>
    void microtubule_uniform_scatter(GraphType& graph, CplexType& cplex, ParamType& settings)
    {
        double epsilon_min = settings.MT_MIN_SEGMENT_INIT;
        double epsilon_max = settings.MT_MAX_SEGMENT_INIT;
        std::random_device random_device;
        std::mt19937 random_engine(random_device());
        auto& grid = cplex.getCoarseGrid();
        
        //first create a grid that needs to fit MTs of two segments
        //without intial overlap
        double max_nx = 
            std::floor((cplex.max_x - cplex.min_x) / (2.25*settings.MT_MAX_SEGMENT_INIT));
        double max_ny = 
            std::floor((cplex.max_y - cplex.min_y) / (2.25*settings.MT_MAX_SEGMENT_INIT));
        
        CartesianGrid2D uniform_grid;
        uniform_grid.init(cplex.min_x, cplex.min_y, 
                cplex.max_x, cplex.max_y,
                max_nx, max_ny);
        std::cout << uniform_grid << "\n";
        std::size_t max_mts = uniform_grid.totalNumCells();
        
        if(settings.NUM_MT > max_mts)
        {
            std::cout << "The number of " << settings.NUM_MT 
                <<  " MTs is not supported for the reqested grid size.\n"
                << "The maximum number allowed and will be used is " << max_mts << "\n"; 
            settings.NUM_MT = max_mts;
        }
        
        //step 1 create an ordered number of selectable cells 
        std::vector<std::size_t> selectable_cells;
        for(auto i = 0; i < max_mts; i++) selectable_cells.push_back(i);
        
        //step 2 randomly select the cells to initialize into
        std::vector<std::size_t> selected_cells;
        
        for(auto i = 0; i < settings.NUM_MT; i++)
        {
            std::uniform_int_distribution<std::size_t>
                resizing_distribution(0, selectable_cells.size()-1);
            auto selected = resizing_distribution(random_engine);
            selected_cells.push_back(selectable_cells[selected]);
            selectable_cells.erase(selectable_cells.begin()+selected);
        }
        
        std::uniform_real_distribution<double> 
            distribution_global_x(cplex.min_x+3*epsilon_max, 
                    cplex.max_x-3*epsilon_max);

        std::uniform_real_distribution<double>
            distribution_global_y(cplex.min_y+3*epsilon_max, 
                    cplex.max_y-3*epsilon_max);

        std::uniform_real_distribution<double> 
            distribution_local(epsilon_min, epsilon_max);
        
        std::uniform_real_distribution<double>
            distribution_angle(0.0, 2.0*3.14);
        using node_type = typename GraphType::node_type;

        std::size_t segments = 3;
        for(auto i = 0; i < settings.NUM_MT; i++) 
        {
            double x_c, y_c;
            uniform_grid.cardinalCellToPoint(x_c, y_c, selected_cells[i]);
            auto z_c = 0.0; 

            auto theta = distribution_angle(random_engine);
            auto seg_len = distribution_local(random_engine);

            auto x_s = 0.0;
            auto y_s = seg_len;
            auto x_r_t = x_s*cos(theta) + y_s*sin(theta);
            auto y_r_t = -x_s*sin(theta) +y_s*cos(theta);
            auto x_r = x_c + x_r_t;
            auto y_r = y_c + y_r_t;
            auto z_r = 0.0;

    
            auto x_l = x_c - (x_r - x_c);
            auto y_l = y_c - (y_r - y_c);
            auto z_l = 0.0;
            
            //compute dist and unit vector
            double p1[3] = {x_r - x_c, y_r - y_c, 0.0};
            double p2[3] = {0.0, 0.0, 0.0};
            double p3[3] = {x_c - x_l, y_c - y_l, 0.0};
            double u1[3], u2[3];

            set_unit_vector(p1, p2, u1);
            set_unit_vector(p3, p2, u2);
            node_type node_l(i*segments, 
                    {{x_l, y_l, z_l}, 
                    {0.0, 0.0, 0.0}, 
                    negative, 
                    {-1, -1, -1}, 
                    {u2[0], u2[1], u2[2]}});

            node_type node_c(i*segments+1, 
                    {{x_c, y_c, z_c}, 
                    {0.0, 0.0, 0.0}, 
                    intermediate, 
                    {-1, -1, -1}, 
                    {u1[0], u1[1], u1[2]}});
            
            node_type node_r(i*segments+2, 
                    {{x_r, y_r, z_r}, 
                    {0.0, 0.0, 0.0}, 
                    positive, 
                    {-1, -1, -1}, 
                    {u1[0], u1[1], u1[2]}});

            graph.addNode(node_l);
            graph.addNode(node_c);
            graph.addNode(node_r);

            graph.addEdge(node_l, node_c);
            graph.addEdge(node_r, node_c);
        }

    }
} //end namespace Plant
} //end namespace DGGML

namespace Plant
{
    template <typename GraphType, typename CplexType, typename ParamType, typename GenType>
    void microtubule_uniform_scatter(GraphType& graph, CplexType& cplex, ParamType& settings, GenType& gen) {
        using node_type = typename GraphType::node_type;
        using key_type = typename node_type::key_type;
        double epsilon_min = settings.DIV_LENGTH;//settings.MT_MIN_SEGMENT_INIT;
        double epsilon_max = settings.DIV_LENGTH;//settings.MT_MAX_SEGMENT_INIT;
        std::random_device random_device;
        std::mt19937 random_engine(random_device());
        auto &grid = cplex.getCoarseGrid();

        DGGML::CartesianGrid2D& reaction_grid = cplex.reaction_grid;

        //create the boundary
        //bottom, left to right
        key_type prev_key; //only the first step doesn't have a prev
        key_type first_key; //needed to complete the loop around the boundary
        for(auto i = 0; i < reaction_grid._nx; i++)
        {
            auto cardinal = reaction_grid.cardinalCellIndex(i, 0);
            double px, py;
            reaction_grid.cardinalCellToPoint(px, py, cardinal);
            key_type curr_key = gen.get_key();
            node_type node_n = {curr_key, {Plant::Boundary{}, px, py, 0.0}};
            graph.addNode(node_n);
            //connect to previous node or its the first
            if(i >= 1) graph.addEdge(prev_key, curr_key);
            else first_key = curr_key;
            prev_key = curr_key;
        }
        //note: the previous gets carried over!
        //right side interior, bottom to top
        for(auto j = 1; j < reaction_grid._ny-1; j++)
        {
            auto cardinal = reaction_grid.cardinalCellIndex(reaction_grid._nx-1, j);
            double px, py;
            reaction_grid.cardinalCellToPoint(px, py, cardinal);
            key_type curr_key = gen.get_key();
            node_type node_n = {curr_key, {Plant::Boundary{}, px, py, 0.0}};
            graph.addNode(node_n);
            //connect to previous node
            graph.addEdge(prev_key, curr_key);
            prev_key = curr_key;
        }
        //note: the previous gets carried over!
        //top, right to left
        for(auto i = reaction_grid._nx-1; i >= 0; i--)
        {
            auto cardinal = reaction_grid.cardinalCellIndex(i, reaction_grid._ny-1);
            double px, py;
            reaction_grid.cardinalCellToPoint(px, py, cardinal);
            key_type curr_key = gen.get_key();
            node_type node_n = {curr_key, {Plant::Boundary{}, px, py, 0.0}};
            graph.addNode(node_n);
            //connect to previous node
            graph.addEdge(prev_key, curr_key);
            prev_key = curr_key;
        }
        //note: the previous gets carried over!
        //left side interior, bottom to top
        for(auto j = reaction_grid._nx-2; j > 0; j--)
        {
            auto cardinal = reaction_grid.cardinalCellIndex(0, j);
            double px, py;
            reaction_grid.cardinalCellToPoint(px, py, cardinal);
            key_type curr_key = gen.get_key();
            node_type node_n = {curr_key, {Plant::Boundary{}, px, py, 0.0}};
            graph.addNode(node_n);
            //connect to previous node
            graph.addEdge(prev_key, curr_key);
            prev_key = curr_key;
        }
        //complete the loop with the first
        graph.addEdge(prev_key, first_key);

        //create the nucleators
        for(auto i = 1; i < reaction_grid._nx-1; i++)
        {
            for(auto j = 1; j < reaction_grid._ny-1; j++)
            {
                ///if(i % 3 == 0 && j % 3 == 0) {
                    auto cardinal = reaction_grid.cardinalCellIndex(i, j);
                    double px, py;
                    reaction_grid.cardinalCellToPoint(px, py, cardinal);
                    node_type node_n = {gen.get_key(), {Plant::Nucleator{}, px, py, 0.0}};
                    graph.addNode(node_n);
                //}
            }
        }
        //first create a grid that needs to fit MTs of two segments without initial overlap
        double max_nx =
                std::floor((cplex.max_x - cplex.min_x) / (4.0 * settings.DIV_LENGTH));//settings.MT_MAX_SEGMENT_INIT));
        double max_ny =
                std::floor((cplex.max_y - cplex.min_y) / (4.0 * settings.DIV_LENGTH));//settings.MT_MAX_SEGMENT_INIT));

        DGGML::CartesianGrid2D uniform_grid;
        uniform_grid.init(cplex.min_x, cplex.min_y,
                          cplex.max_x, cplex.max_y,
                          max_nx, max_ny);
        std::size_t max_mts = uniform_grid.totalNumCells();
        //std::cout << "Max mts " << max_mts << "\n"; std::cin.get();
        if (settings.NUM_MT > max_mts) {
            std::cout << "The number of " << settings.NUM_MT
                      << " MTs is not supported for the requested grid size.\n"
                      << "The maximum number allowed and will be used is " << max_mts << "\n";
            settings.NUM_MT = max_mts;
        }

        //step 1 create an ordered number of selectable cells
        std::vector<std::size_t> selectable_cells;
        for (auto i = 0; i < max_mts; i++) selectable_cells.push_back(i);

        //step 2 randomly select the cells to initialize into
        std::vector<std::size_t> selected_cells;

        for (auto i = 0; i < settings.NUM_MT; i++) {
            std::uniform_int_distribution<std::size_t>
                    resizing_distribution(0, selectable_cells.size() - 1);
            auto selected = resizing_distribution(random_engine);
            selected_cells.push_back(selectable_cells[selected]);
            selectable_cells.erase(selectable_cells.begin() + selected);
        }

        std::uniform_real_distribution<double>
                distribution_global_x(cplex.min_x + 3 * epsilon_max,
                                      cplex.max_x - 3 * epsilon_max);

        std::uniform_real_distribution<double>
                distribution_global_y(cplex.min_y + 3 * epsilon_max,
                                      cplex.max_y - 3 * epsilon_max);

        std::uniform_real_distribution<double>
                distribution_local(epsilon_min, epsilon_max);

        std::uniform_real_distribution<double>
                distribution_angle(0.0, 2.0 * 3.14);

        std::size_t segments = 3;
        for (auto i = 0; i < settings.NUM_MT; i++) {
            double x_c, y_c;
            uniform_grid.cardinalCellToPoint(x_c, y_c, selected_cells[i]);
            auto z_c = 0.0;

            auto theta = distribution_angle(random_engine);
            auto seg_len = distribution_local(random_engine);

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
            double p1[3] = {x_r - x_c, y_r - y_c, 0.0};
            double p2[3] = {0.0, 0.0, 0.0};
            double p3[3] = {x_c - x_l, y_c - y_l, 0.0};
            double u1[3], u2[3];

            DGGML::set_unit_vector(p1, p2, u1);
            DGGML::set_unit_vector(p3, p2, u2);
            Plant::graph_type tg;
            node_type node_l = {gen.get_key(),//i*segments,
                        {Plant::Negative{0.0,0.0,0.0,
                                             u2[0], u2[1], u2[2]},
                         x_l, y_l, z_l}};

            node_type node_c ={gen.get_key(),//i*segments+1,
                        {Plant::Intermediate{0.0,0.0,0.0,
                                             u1[0], u1[1], u1[2]},
                         x_c, y_c, z_c}};

            node_type node_r ={gen.get_key(),//i*segments+2,
                        {Plant::Positive{0.0,0.0,0.0,
                                             u2[0], u2[1], u2[2]},
                         x_r, y_r, z_r}};

            graph.addNode(node_l);
            graph.addNode(node_c);
            graph.addNode(node_r);

            graph.addEdge(node_l, node_c);
            graph.addEdge(node_r, node_c);
        }
    }
}

#endif
