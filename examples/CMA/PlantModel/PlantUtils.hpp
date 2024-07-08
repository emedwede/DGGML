#ifndef PLANT_UTILS_HPP
#define PLANT_UTILS_HPP

#include <random>
#include <algorithm>
#include <vector>

#include "PlantTypes.hpp"
#include "Utlities/MathUtils.hpp"

#include "CartesianGrid2D.hpp"

namespace CMA {

    template<typename StreamType, typename DataType>
    void print_numpy_array_stats(StreamType &out, DataType &data, std::string var_name) {
        out << var_name << " = np.asarray([";
        for (auto i = 0; i < data.size(); i++) {
            if (i != 0 && i % 20 == 0) out << "\n";
            if (i != data.size() - 1)
                out << data[i] << ", ";
            else
                out << data[i];
        }
        out << "]);\n";
    }

    template<typename StreamType, typename DataType>
    void print_json_array_stats(StreamType &out, DataType &data, std::string var_name) {
        out << "\t\"" << var_name << "\": [";
        for (auto i = 0; i < data.size(); i++) {
            if (i != data.size() - 1)
                out << data[i] << ", ";
            else
                out << data[i];
        }
        out << "]";
    }
}
namespace Plant {
    template<typename GraphType, typename CplexType, typename ParamType, typename GenType>
    void microtubule_uniform_scatter(GraphType &graph, CplexType &cplex, ParamType &settings, GenType &gen) {
        using node_type = typename GraphType::node_type;
        using key_type = typename node_type::key_type;
        double epsilon_min = settings.MT_MIN_SEGMENT_INIT; //settings.DIV_LENGTH;
        double epsilon_max = settings.MT_MAX_SEGMENT_INIT; //settings.DIV_LENGTH;
        std::random_device random_device;
        std::mt19937 random_engine(random_device());
        auto &grid = cplex.getCoarseGrid();

        DGGML::CartesianGrid2D &reaction_grid = cplex.reaction_grid;

        //TODO: make it so the end points are aligned with the cell corners, but have custom spacing
        // between the points
        //create the boundary
        //bottom, left to right
//        key_type prev_key; //only the first step doesn't have a prev
//        key_type first_key; //needed to complete the loop around the boundary
//        for (auto i = 0; i < reaction_grid._nx; i++) {
//            auto cardinal = reaction_grid.cardinalCellIndex(i, 0);
//            double px, py;
//            reaction_grid.cardinalCellToPoint(px, py, cardinal);
//            key_type curr_key = gen.get_key();
//            node_type node_n = {curr_key, {Plant::Boundary{}, px, py, 0.0}};
//            graph.addNode(node_n);
//            //connect to previous node or its the first
//            if (i >= 1) graph.addEdge(prev_key, curr_key);
//            else first_key = curr_key;
//            prev_key = curr_key;
//        }
//        //note: the previous gets carried over!
//        //right side interior, bottom to top
//        for (auto j = 1; j < reaction_grid._ny - 1; j++) {
//            auto cardinal = reaction_grid.cardinalCellIndex(reaction_grid._nx - 1, j);
//            double px, py;
//            reaction_grid.cardinalCellToPoint(px, py, cardinal);
//            key_type curr_key = gen.get_key();
//            node_type node_n = {curr_key, {Plant::Boundary{}, px, py, 0.0}};
//            graph.addNode(node_n);
//            //connect to previous node
//            graph.addEdge(prev_key, curr_key);
//            prev_key = curr_key;
//        }
//        //note: the previous gets carried over!
//        //top, right to left
//        for (auto i = reaction_grid._nx - 1; i >= 0; i--) {
//            auto cardinal = reaction_grid.cardinalCellIndex(i, reaction_grid._ny - 1);
//            double px, py;
//            reaction_grid.cardinalCellToPoint(px, py, cardinal);
//            key_type curr_key = gen.get_key();
//            node_type node_n = {curr_key, {Plant::Boundary{}, px, py, 0.0}};
//            graph.addNode(node_n);
//            //connect to previous node
//            graph.addEdge(prev_key, curr_key);
//            prev_key = curr_key;
//        }
//        //note: the previous gets carried over!
//        //left side interior, bottom to top
//        for (auto j = reaction_grid._ny - 2; j > 0; j--) {
//            auto cardinal = reaction_grid.cardinalCellIndex(0, j);
//            double px, py;
//            reaction_grid.cardinalCellToPoint(px, py, cardinal);
//            key_type curr_key = gen.get_key();
//            node_type node_n = {curr_key, {Plant::Boundary{}, px, py, 0.0}};
//            graph.addNode(node_n);
//            //connect to previous node
//            graph.addEdge(prev_key, curr_key);
//            prev_key = curr_key;
//        }
//        //complete the loop with the first
//        graph.addEdge(prev_key, first_key);

        if(settings.ENABLE_CREATION) {
            //create the nucleators
            for (auto i = 1; i < reaction_grid._nx - 1; i++) {
                for (auto j = 1; j < reaction_grid._ny - 1; j++) {
                    ///if(i % 3 == 0 && j % 3 == 0) {
                    auto cardinal = reaction_grid.cardinalCellIndex(i, j);
                    double px, py;
                    reaction_grid.cardinalCellToPoint(px, py, cardinal);
                    node_type node_n = {gen.get_key(), {Plant::Nucleator{}, px, py, 0.0}};
                    graph.addNode(node_n);
                    //}
                }
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
                                {Plant::Negative{0.0, 0.0, 0.0,
                                                 u2[0], u2[1], u2[2]},
                                 x_l, y_l, z_l}};

            node_type node_c = {gen.get_key(),//i*segments+1,
                                {Plant::Intermediate{0.0, 0.0, 0.0,
                                                     u1[0], u1[1], u1[2]},
                                 x_c, y_c, z_c}};

            node_type node_r = {gen.get_key(),//i*segments+2,
                                {Plant::Positive{0.0, 0.0, 0.0,
                                                 u2[0], u2[1], u2[2]},
                                 x_r, y_r, z_r}};

            graph.addNode(node_l);
            graph.addNode(node_c);
            graph.addNode(node_r);

            graph.addEdge(node_l, node_c);
            graph.addEdge(node_r, node_c);
        }
    }


    double compute_correlation(const std::vector<double> &x, const std::vector<double> &y) {
        // Compute means
        double mean_x = std::accumulate(x.begin(), x.end(), 0.0) / x.size();
        double mean_y = std::accumulate(y.begin(), y.end(), 0.0) / y.size();

        // Compute variances
        double var_x = std::accumulate(x.begin(), x.end(), 0.0,
                                       [mean_x](double acc, double xi) {
                                           return acc + (xi - mean_x) * (xi - mean_x);
                                       }) / (x.size() - 1);
        double var_y = std::accumulate(y.begin(), y.end(), 0.0,
                                       [mean_y](double acc, double yi) {
                                           return acc + (yi - mean_y) * (yi - mean_y);
                                       }) / (y.size() - 1);

        // Check for zero variance
        if (var_x == 0 || var_y == 0) {
            return 0; // or any other value you choose to return in this case
        }

        // Compute covariance
        double cov = 0;
        for (size_t i = 0; i < x.size(); ++i) {
            cov += (x[i] - mean_x) * (y[i] - mean_y);
        }
        cov /= x.size() - 1;

        // Compute correlation coefficient
        double corr = cov / std::sqrt(var_x * var_y);

        return corr;
    }


    struct Point {
        double p_x, p_y;
        double u_x, u_y;
    };

    double dot_product(Point &p1, Point &p2) {
        return p1.u_x * p2.u_x + p1.u_y * p2.u_y;
    }

    double angle_ref(Point &p1, Point &p2) {
        double dot_prod = dot_product(p1, p2);
        return std::acos(std::min(std::abs(dot_prod), 1.0));
    }

    template<typename GraphType>
    double compute_correlation(GraphType &system_graph) {
        std::vector<double> x;
        std::vector<double> y;
        for (auto &node: system_graph.getNodeSetRef()) {
            auto &u = node.second.getData().unit_vec;
            double u0 = u[0];
            double u1 = u[1];

            x.push_back(u0);
            y.push_back(u1);
        }

        return compute_correlation(x, y);
    }

    std::vector<double> two_point_correlation_alpha(std::vector<Point> &points, double bin_size, double max_distance) {
        std::size_t num_bins = std::ceil(max_distance / bin_size);
        std::vector<double> correlation(num_bins, 0.0);

        std::vector<int> num_pairs(num_bins, 0);

        auto num_points = points.size();
        for (int i = 0; i < num_points; i++) {
            Point p1 = points[i];
            for (int j = i + 1; j < num_points; j++) {
                Point p2 = points[j];
                double distance = hypot(p1.p_x - p2.p_x, p1.p_y - p2.p_y);
                if (distance <= max_distance) {
                    //double angle = angle_ref(p1, p2);
                    //if (angle >= alpha1 && angle <= alpha2) {
                    int bin_index = floor(distance / bin_size);
                    auto abs_dot = dot_product(p1, p2) * dot_product(p1, p2);//std::abs(dot_product(p1, p2));
                    auto dot = dot_product(p1, p2);
                    if (dot > 1.0001 || dot < -1.0001) {std::cout << dot << "\n"; std::cin.get();}
                    if (abs_dot > 1.0) abs_dot = 1.0;
//                    if(bin_index < 0 || bin_index >= num_bins) {
//                        std::cout << "bin index error " << bin_index << "\n"; std::cin.get();
//                    }
                    correlation[bin_index] += abs_dot;//std::acos(abs_dot); //std::cos(angle);
                    num_pairs[bin_index]++;
                    //}
                }
            }
        }

        for (int i = 0; i < num_bins; i++) {
            if (num_pairs[i] > 0) {
                correlation[i] /= num_pairs[i];
            }
        }
        return correlation;
    }

    template<typename GraphType, typename ParamType>
    std::vector<double> compute_two_point_correlation_alpha(GraphType &system_graph, ParamType &settings) {
        std::vector<Point> points;
        for(auto& node : system_graph.getNodeSetRef())
        {
            auto& node_data = node.second.getData();
            if(std::holds_alternative<Plant::Intermediate>(node_data.data)) {
                auto& p = node_data.position;
                auto& u = std::get<Plant::Intermediate>(node_data.data).unit_vec;

                points.push_back({p[0], p[1], u[0], u[1]});
            }
            else if(std::holds_alternative<Plant::Positive>(node_data.data)) {
                auto& p = node_data.position;
                auto& u = std::get<Plant::Positive>(node_data.data).unit_vec;

                points.push_back({p[0], p[1], u[0], u[1]});
            }
            else if(std::holds_alternative<Plant::Negative>(node_data.data)) {
                auto& p = node_data.position;
                auto& u = std::get<Plant::Negative>(node_data.data).unit_vec;

                points.push_back({p[0], p[1], u[0], u[1]});
            }
        }

        double bin_size = settings.MAXIMAL_REACTION_RADIUS;
        double max_distance = std::max(settings.CELL_NX*settings.CELL_DX, settings.CELL_NY*settings.CELL_DY);
        return two_point_correlation_alpha(points, bin_size, max_distance);
    }

    template<typename GraphType, typename ParamType>
    std::vector<std::size_t> compute_orientation_histogram(GraphType &system_graph, ParamType &settings, std::vector<double>& angle_vec) {
        std::vector<Point> points;
        for(auto& node : system_graph.getNodeSetRef())
        {
            auto& node_data = node.second.getData();
            if(std::holds_alternative<Plant::Intermediate>(node_data.data)) {
                auto& p = node_data.position;
                auto& u = std::get<Plant::Intermediate>(node_data.data).unit_vec;

                points.push_back({p[0], p[1], u[0], u[1]});
            }
            else if(std::holds_alternative<Plant::Positive>(node_data.data)) {
                auto& p = node_data.position;
                auto& u = std::get<Plant::Positive>(node_data.data).unit_vec;

                points.push_back({p[0], p[1], u[0], u[1]});
            }
            else if(std::holds_alternative<Plant::Negative>(node_data.data)) {
                auto& p = node_data.position;
                auto& u = std::get<Plant::Negative>(node_data.data).unit_vec;

                points.push_back({p[0], p[1], u[0], u[1]});
            }
        }
        //build bins in 10 degree increments
        std::size_t num_bins = 180/10;
        std::cout << "Point size " << points.size() << "\n";
        std::vector<std::size_t> orientation_histogram(num_bins);
        for(auto& p : points)
        {
            //TODO: fix this so that we are in the right half of the plane!
            //flip to quadrant 1, compute the angle with the terminal unit vector
            double angle = std::acos(std::abs(p.u_x))*(180.0/3.14);
            //if we are in originally in quadrant 2, we must create a mirror image in quadrant 4, so subtract 90
            if(p.u_x < 0 && p.u_y > 0 || p.u_x > 0 && p.u_y < 0)
                angle = -angle;
            //if we are in quadrant 3, we must create a mirror image in quadrant 1
            //Note: anything in quadrant 1 and 4 are as is
            //auto angle = p.u_x == 0 ? 0.0 : std::atan(p.u_y/std::abs(p.u_x))*(180.0/3.14);
            //angle += 90.0;
            //std::cout << angle << " ";
            //scale the angle
            angle += 90.0;
            angle = angle > 180.0 ? 180.0 : angle; //ensures precision errors don't break the code
            angle = angle < 0.0 ? 0.0 : angle; //ensures precision errors don't break the code
            int bid = angle == 0.0 ? 0 : int(std::ceil(angle/10.0))-1;
            if(bid >= 0 && bid < num_bins) {
                orientation_histogram[bid]++;
                angle_vec.push_back(angle);
            }
        }
        //std::cout << "\n";
        auto total = std::accumulate(orientation_histogram.begin(), orientation_histogram.end(), 0);
        std::cout << "\nhist total: " << total << "\n";
        return std::move(orientation_histogram);
    }
}
#endif
