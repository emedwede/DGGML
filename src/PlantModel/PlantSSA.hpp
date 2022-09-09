#ifndef CAJETE_PLANT_SSA_HPP
#define CAJETE_PLANT_SSA_HPP 

#include "PlantTypes.hpp"
#include "PlantGrammar.hpp"

#include "ExpandedComplex2D.hpp"
#include "YAGL_Graph.hpp"

#include "SundialsUtils.hpp"

#include "RuleSystem.hpp"

#include <chrono>
#include <random>
#include <algorithm>
#include <vector>
#include <iostream>
#include <math.h>
#include <map>

#include <arkode/arkode_erkstep.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h> 
#include <sundials/sundials_math.h>

#if defined(SUNDIALS_EXTENDED_PRECISION)
#define GSYM "Lg"
#define ESYM "Le"
#define FSYM "Lf"
#else
#define GSYM "g"
#define ESYM "e"
#define FSYM "f"
#endif
namespace Cajete 
{
namespace Plant 
{

//TODO: Move to it's own random class or header
auto RandomRealsBetween = [](double low, double high)
{
    auto randomFunc = 
        [distribution_ = std::uniform_real_distribution<double>(low, high),
         random_engine_ = std::mt19937{std::random_device{}() }]() mutable
       {
           return distribution_(random_engine_);
       };
    return randomFunc;
};

auto RandomIntsBetween = [](int low, int high)
{
    auto randomFunc = 
        [distribution_ = std::uniform_int_distribution<int>(low, high),
         random_engine_ = std::mt19937{std::random_device{}() }]() mutable
       {
           return distribution_(random_engine_);
       };
    return randomFunc;
};

template <typename GraphType, typename MatchSetType>
void filter_matches(MatchSetType* unfiltered, MatchSetType* filtered, std::size_t N, GraphType& system_graph, std::size_t k, std::size_t dim)
{
    for(auto i = 0; i < N; i++) 
    {
        for(auto match : unfiltered[i])
        {
            bool bad_match = false;
            //check match integrity
            for(auto key : match)
            {
                auto dtag = system_graph.findNode(key)->second.getData().tagND[dim];
                if(dtag != k)
                {
                     bad_match = true;
                     break;
                }
            } 

            if(!bad_match)
            {
                filtered[i].push_back(match);
            }
        }
        unfiltered[i].clear();
    }
}

int rhs_func(realtype t, N_Vector y, N_Vector ydot, void* user_data) 
{
    
    return 0;
}

int root_func(realtype t, N_Vector y, realtype* gout, void* user_data)
{
    return 0;
}

template <typename UserDataType>
struct Solver 
{
    int flag; //generic reusable flag 
    
    //Create the suncontext 
    SUNContext ctx;
    void* arkode_mem;

    realtype t_start, t_final, t, tout, dt_out;
     
    realtype reltol;
    realtype abstol;
    
    static constexpr int num_roots = 1;
    int roots_found[num_roots];
    int root_flag;

    sunindextype num_eq;
    
    N_Vector y;
   
    //TODO: should this be a pointer?
    UserDataType user_data;
    
    Solver(UserDataType _user_data, sunindextype _num_eq = 0) 
        : user_data(_user_data)
    {
        num_eq = _num_eq;
        t_start = RCONST(0.0);
        dt_out = RCONST(user_data.settings.DELTA_T_MIN);
        t = t_start;
        tout = t_start + dt_out;

        flag = SUNContext_Create(NULL, &ctx);
        if(!Cajete::SundialsUtils::check_flag(&flag, "SUNContext_Create", 1))
            std::cout << "Passed the error check, suncontext created\n";
        
        num_eq = 0;
        for(const auto& r : user_data.rule_map[user_data.k])
        {
            num_eq += 2*DIM3D*user_data.rule_system[r].size();
        }
        
        //num_eq = 2;
        std::cout << "*****NumEQ="<<num_eq<<"*****\n";
        y = NULL; 
        y = N_VNew_Serial(num_eq, ctx);
        
        auto j = 0;
        for(const auto& r : user_data.rule_map[user_data.k])
        {
            for(const auto& m : user_data.rule_system[r])
            {
                auto& node_data = user_data.system_graph.findNode(m)->second.getData();
                for(int i = 0; i < DIM3D; i++)
                    NV_Ith_S(y, j++) = node_data.velocity[i];
                for(int i = 0; i < DIM3D; i++)
                    NV_Ith_S(y, j++) = node_data.position[i];
            } 
        }
        //NV_Ith_S(y, 0) = 50.0;
        //NV_Ith_S(y, 1) = 0.0;

        //std::cout << "J: " << j << "\n";
        arkode_mem = ERKStepCreate(my_func, t_start, y, ctx); 
        //if (!Cajete::SundialsUtils::check_flag((void *)arkode_mem, "ERKStepCreate", 1))
        //    std::cout << "Passed the error check, stepper initialized\n";

        reltol = 1.0e-6;
        abstol = 1.0e-10;

        flag = ERKStepSStolerances(arkode_mem, reltol, abstol);
        
        //set optional inputs
        flag = ERKStepSetUserData(arkode_mem, &user_data); 
    
        // specify the root finding function having num_roots  
        //flag = ERKStepRootInit(arkode_mem, num_roots, root_func); 
        
        //ERKStepSetMinStep(arkode_mem, user_data.settings.DELTA_T_MIN/10.0);
        ERKStepSetMaxStep(arkode_mem, dt_out);//dt_out/10.0);
    }
    
    void step()
    {
        if(num_eq == 0) return;
        flag = ERKStepEvolve(arkode_mem, tout, y, &t, ARK_NORMAL);
        if(flag == ARK_ROOT_RETURN)
        {
            root_flag = ERKStepGetRootInfo(arkode_mem, roots_found);
            std::cout << "A root has been found\n";
        }
        
        //successful solve
        if(flag >= 0)
        {
            tout += dt_out;
        }
    }
    
    static int my_func(realtype t, N_Vector y, N_Vector ydot, void* user_data) 
    {
        UserDataType* local_ptr = (UserDataType *)user_data;
        auto k = local_ptr->k;
        auto i = 0;
        auto v_plus = local_ptr->settings.V_PLUS;
        auto d_l = local_ptr->settings.DIV_LENGTH;

        //NV_Ith_S(ydot, 0) = NV_Ith_S(y, 1);
        //NV_Ith_S(ydot, 1) = -9.8;

        for(const auto& r : local_ptr->rule_map[k])
        {
            auto& instance = local_ptr->rule_system[r];
            if(instance.type == Rule::G)
            {
                auto id = instance[0];
                auto& node_i_data = local_ptr->system_graph.findNode(id)->second.getData();
                double l = 0.0;
                for(auto j = 0; j < 3; j++)
                {
                    double diff = NV_Ith_S(y, i+j+3) - node_i_data.position[j];
                    l += diff*diff;
                }
                l = sqrt(l);
                double length_limiter = (1.0 - (l/d_l));
                std::cout << d_l << " " << length_limiter << "\n";
                //std::cout << Rule::G << "\n";
                NV_Ith_S(ydot, i) = v_plus*node_i_data.unit_vec[0]*length_limiter;
                NV_Ith_S(ydot, i+1) = v_plus*node_i_data.unit_vec[1]*length_limiter;
                NV_Ith_S(ydot, i+2) = v_plus*node_i_data.unit_vec[2]*length_limiter;
                NV_Ith_S(ydot, i+3) = v_plus*node_i_data.unit_vec[0]*length_limiter;
                NV_Ith_S(ydot, i+4) = v_plus*node_i_data.unit_vec[1]*length_limiter;
                NV_Ith_S(ydot, i+5) = v_plus*node_i_data.unit_vec[2]*length_limiter;
                //NV_Ith_S(ydot, i+3) = NV_Ith_S(y, i);
                //NV_Ith_S(ydot, i+4) = NV_Ith_S(y, i+1);
                //NV_Ith_S(ydot, i+5) = NV_Ith_S(y, i+2);
                for(int j = i+6; j < i+12; j++)
                    NV_Ith_S(ydot, j) = 0.0;
            }
            if(local_ptr->rule_system[r].type == Rule::R)
            {
                for(int j = i; j < i+12; j++)
                    NV_Ith_S(ydot, j) = 0.0;
                //std::cout << Rule::R << "\n";
            }
            i += 12;
        }
        //std::cout << "In my function for cell " << k 
        //    << ", local rule size: " << local_ptr->rule_map[k].size() << "\n";  
        local_ptr = nullptr;
        return 0;
    }
    
    void copy_back()
    {
        auto j = 0;
        for(const auto& r : user_data.rule_map[user_data.k])
        {
            for(const auto& m : user_data.rule_system[r])
            {
                auto& node_data = user_data.system_graph.findNode(m)->second.getData();
                for(int i = 0; i < DIM3D; i++)
                    node_data.velocity[i] = NV_Ith_S(y, j++);
                for(int i = 0; i < DIM3D; i++)
                    node_data.position[i] = NV_Ith_S(y, j++);
            } 
        }

    }
    void print_stats()
    {
        long int nst, nst_a, nfe, netf;
        /* Print some final statistics */
        flag = ERKStepGetNumSteps(arkode_mem, &nst);
        flag = ERKStepGetNumStepAttempts(arkode_mem, &nst_a);
        flag = ERKStepGetNumRhsEvals(arkode_mem, &nfe);
        flag = ERKStepGetNumErrTestFails(arkode_mem, &netf);

        printf("\nFinal Solver Statistics:\n");
        printf("   Internal solver steps = %li (attempted = %li)\n", nst, nst_a);
        printf("   Total RHS evals = %li\n", nfe);
        printf("   Total number of error test failures = %li\n\n", netf);
    }
    ~Solver()
    {
        N_VDestroy(y);
        ERKStepFree(&arkode_mem); //free the solver memory
        SUNContext_Free(&ctx); //always call prior to MPI_Finalize
        std::cout << "Destroyed the solver\n";
    
        /* Step 13: Finalilze MPI, if used */ 
    }
};

using RuleSystemType = Cajete::RuleSystem<Plant::mt_key_type>;
template <typename B, typename T, typename U, typename GeoplexType, typename GraphType, typename ParamType>
void plant_model_ssa_new(RuleSystemType& rule_system, B& k, T& rule_map, U& cell_list, GeoplexType& geoplex2D, GraphType& system_graph, ParamType& settings)
{
    std::cout << "Cell " << k << " has " << rule_map[k].size() << " rules\n";

    double delta_t, exp_sample, tau, geocell_propensity;
    std::size_t events, steps; 
    delta_t = 0.0; events = 0; steps = 0; geocell_propensity = 0.0;
    
    //need an aggregation
    typedef struct 
    {
        RuleSystemType& rule_system;
        GraphType& system_graph;
        B& k;
        T& rule_map;
        U& cell_list;
        ParamType& settings;
    } PackType;
    
    Solver ode_system(PackType{rule_system, system_graph, k, rule_map, cell_list, settings}, 2);

    while(delta_t < settings.DELTA) 
    {
        //reset tau
        tau = 0.0;

        double uniform_sample = RandomRealsBetween(0.0, 1.0)();
        
        //sample the exponential variable
        exp_sample = -log(1-uniform_sample);
        
        //the sums of the particular rules propensities
            
        //the propensities calculated for a particular rule
        
        while(delta_t < settings.DELTA && tau < exp_sample)
        {
                    
            // STEP(1) : sum all of the rule propensities
            //zero the geocell_propensity
            geocell_propensity = 0.0;
            
           
            //the step adapts based on propensity or systems fastest dynamic
            //settings.DELTA_T_MIN = 
            //std::min(1.0/(10.0*geocell_propensity), settings.DELTA_DELTA_T);
            
            // STEP(2) : solve the system of ODES 
            ode_system.step();
            
            // STEP(3) : use forward euler to solve the TAU ODE
            tau += geocell_propensity*settings.DELTA_T_MIN; 
            
            // STEP(4) : advance the loop timer
            delta_t += settings.DELTA_T_MIN; //TODO: make delta_t adaptive
            
            //std::cout << "tout: " << ode_system.t << " , DT: " << delta_t << "\n";
            
            steps++;
        }
    

        // if we get over our threshold an event can occur
        if(tau > exp_sample) 
        {
            //determine which rule to file and fire it
        }
    }
    std::cout << "Total steps taken: " << steps << "\n";
    ode_system.copy_back();
    ode_system.print_stats();
}


//template <typename BucketType>
//void plant_model_ssa(
//        BucketType& bucket,
//        ExpandedComplex2D<>& geoplex2D, 
//        YAGL::Graph<mt_key_type, MT_NodeData>& system_graph) 
template <typename BucketType, typename GeoplexType, typename GraphType, typename ParamType>
void plant_model_ssa(BucketType& bucket, GeoplexType& geoplex2D, GraphType& system_graph, ParamType& settings)
{
    double delta_t, exp_sample, tau, geocell_propensity;
    std::size_t events, steps; 
    delta_t = 0.0; events = 0; steps = 0; geocell_propensity = 0.0;

    /* Begin the generic sundials skeleton */ 
    
    /* Step 1: initalize the parallel env */
    
    // For now, serial, so skip
    
    /* Step 2: Create the sun context */
    
    int flag; //generic reusable flag 

    //Create the suncontext 
    SUNContext ctx;
    flag = SUNContext_Create(NULL, &ctx);
    if(!Cajete::SundialsUtils::check_flag(&flag, "SUNContext_Create", 1))
        std::cout << "Passed the error check, suncontext created\n";
    
    /* Step 3: set the problem dimensions */
    realtype t_start, t_final, t, tout;
    sunindextype num_eq;
    
    /* Step 4: set the vector of initial values */ 
    N_Vector y = NULL; 
    y = N_VNew_Serial(num_eq, ctx);

    /* Step 5: create the explicit stepper object */ 
    void* arkode_mem = NULL;
    arkode_mem = ERKStepCreate(rhs_func, t_start, y, ctx); 

    /* Step 6: specify the integration tolerances */ 
    realtype reltol = 1.0e-6;
    realtype abstol = 1.0e-10;
    flag = ERKStepSStolerances(arkode_mem, reltol, abstol);

    /* Step 7: set any optional inputs */
    //flag = ERKStepSetUserData(arkode_mem, user_data); 
    
    /* Step 8: specify an optional root finding problem to solve */ 
    int num_roots = 1;
    int roots_found[num_roots];
    int root_flag;
    // specify the root finding function having num_roots  
    flag = ERKStepRootInit(arkode_mem, num_roots, root_func); 

    /* Step 9: advance the solution in time */ 

    /* Step 10: get optional outputs */ 
    
    /* Step 11: deallocate memory for solution vector */
    N_VDestroy(y);

    /* Step 12: free the solver memory and the context if not reused */ 
    ERKStepFree(&arkode_mem); //free the solver memory
    SUNContext_Free(&ctx); //always call prior to MPI_Finalize
    std::cout << "Freeing the suncontext\n";
    
    /* Step 13: Finalilze MPI, if used */ 

    /* End the generic sundials skeleton */

    while(delta_t < settings.DELTA) 
    {
        //reset tau
        tau = 0.0;
        auto k = bucket.first;
        std::size_t num_patterns = 5;
        std::vector<std::vector<mt_key_type>> rule_matches[num_patterns];
        std::vector<std::vector<mt_key_type>> unfiltered_matches[num_patterns];
        unfiltered_matches[0] = microtubule_growing_end_matcher(system_graph, bucket.second); 
        unfiltered_matches[1] = microtubule_retraction_end_matcher(system_graph, bucket.second); 
        unfiltered_matches[2] = microtubule_retraction_end_two_intermediate_matcher(system_graph, bucket.second);
        unfiltered_matches[3] = three_intermediate_matcher(system_graph, bucket.second);
        std::size_t dim = geoplex2D.getGraph().findNode(k)->second.getData().type;
        auto total_matches = 0; for(auto i = 0; i < num_patterns; i++) total_matches += unfiltered_matches[i].size(); 
        //std::cout << "Found " << total_matches << " unfiltered matches\n";
        filter_matches(unfiltered_matches, rule_matches, num_patterns, system_graph, k, dim);
        total_matches = 0; for(auto i = 0; i < num_patterns; i++) total_matches += rule_matches[i].size(); 
        //std::cout << "Found " << total_matches << " filtered matches\n";

        //std::cout << "Total growing ends: " << rule_matches[0].size() << "\n";
        //std::cout << "Total wildcard_matches: " << rule_matches[3].size() << "\n";

        collision_match_refinement(system_graph, 3.0*settings.DIV_LENGTH, rule_matches[0], rule_matches[3], rule_matches[4]);

        //std::cout << "Total potential collisions: " << rule_matches[4].size() << "\n";

        double uniform_sample = RandomRealsBetween(0.0, 1.0)();
        
        //sample the exponential variable
        exp_sample = -log(1-uniform_sample);
        
        //the sums of the particular rules propensities
        double rule_propensities[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
            
        //the propensities calculated for a particular rule
        std::vector<double> prop[8];
        prop[0].reserve(rule_matches[0].size());
        prop[1].reserve(rule_matches[2].size());
        prop[2].reserve(rule_matches[4].size());
        prop[3].reserve(rule_matches[0].size());
        prop[4].reserve(rule_matches[1].size());
        prop[5].reserve(rule_matches[4].size());
        prop[6].reserve(rule_matches[4].size());
        prop[7].reserve(rule_matches[3].size());

        std::vector<std::vector<typename GraphType::node_type>> prev_state[2];
        for(auto i = 0; i < 2; i++)
        {
            for(auto& match : rule_matches[i])
            {
                prev_state[i].push_back({system_graph.findNode(match[0])->second, 
                                         system_graph.findNode(match[1])->second});
            }
        }
        
        while(delta_t < settings.DELTA && tau < exp_sample)
        {
            // STEP(0) : store old state and find all the matches 
            //do a deep copy? TODO: only need to copy old state parameters
            auto& system_graph_old = system_graph;
                    
            // STEP(1) : sum all of the rule propensities
            //zero the geocell_propensity
            geocell_propensity = 0.0;
            
            //reset
            for(auto i = 0; i < 8; i++)
            {
                rule_propensities[i] = 0.0; 
                prop[i].clear();
            }

            for(auto& match : rule_matches[0])
            {
                auto rho = 
                    microtubule_growing_end_polymerize_propensity(system_graph, match, settings);
                rule_propensities[0] += rho;
                prop[0].push_back(rho);

                rho = microtubule_positive_to_negative_propensity(system_graph, match, settings);
                rule_propensities[3] += rho;
                prop[3].push_back(rho); 
                
            }
            
            for(auto& match : rule_matches[1])
            {
                auto rho =
                    microtubule_negative_to_positive_propensity(system_graph, match, settings);
                rule_propensities[4] += rho;
                prop[4].push_back(rho);
            }

            for(auto& match : rule_matches[2]) 
            {
                auto rho = 
                    microtubule_retraction_end_depolymerize_propensity(system_graph, match, settings);
                rule_propensities[1] += rho;
                prop[1].push_back(rho);
             }
            
            for(auto& match : rule_matches[3])
            {
                auto rho =
                    microtubule_katanin_severing_propensity(system_graph, match, settings);
                rule_propensities[7] += rho;
                prop[7].push_back(rho);
            }
            for(auto& match : rule_matches[4])
            {
                auto rho = settings.RHO_TEST_RATE*microtubule_collision_crossover_propensity(system_graph, match, settings);
                rule_propensities[2] += rho;
                rule_propensities[5] += rho;
                rule_propensities[6] += rho;
                prop[2].push_back(rho);
                prop[5].push_back(rho);
                prop[6].push_back(rho);

            } 
            
            geocell_propensity += rule_propensities[0] + rule_propensities[1] + rule_propensities[2] + rule_propensities[3] + rule_propensities[4] + rule_propensities[5] + rule_propensities[6] + rule_propensities[7];
           
            //the step adapts based on propensity or systems fastest dynamic
            settings.DELTA_T_MIN = std::min(1.0/(10.0*geocell_propensity), settings.DELTA_DELTA_T);
            
            // STEP(2) : solve the system of ODES 
            microtubule_ode_solver(rule_matches, system_graph, prev_state, system_graph_old, k, settings); 

            // STEP(3) : use forward euler to solve the TAU ODE
            tau += geocell_propensity*settings.DELTA_T_MIN; //TODO: we need to be careful not to oversolve
            
                        // STEP(4) : advance the loop timer
            delta_t += settings.DELTA_T_MIN; //TODO: make delta_t adaptive

            steps++;
        }
    

        // if we get over our threshold an event can occur
        if(tau > exp_sample) 
        {
            //determine which rule to file and fire it
            microtubule_rule_firing(rule_matches, system_graph, bucket, prop, settings);
        }
    }
    std::cout << "Total steps taken: " << steps << "\n";
}

//be careful with solving, could lead to a segmentation fault 
//in the sorting phase if a parameter is solved out to out of bounds
template <typename MatchSetType, typename GraphType, typename StateType, typename KeyType, typename ParamType>
void microtubule_ode_solver(MatchSetType* all_matches, 
        GraphType& system_graph,
        StateType& prev_state,
        GraphType& system_graph_old,
        KeyType& k,
        ParamType& settings)
{
    for(auto i = 0; i < 2; i++) 
    {
        for(auto j = 0; j < all_matches[i].size(); j++)
        {
           if(i == 0)
                microtubule_growing_end_polymerize_solve(system_graph, prev_state[i][j], system_graph_old, all_matches[i][j], settings);
           if(i == 1)
               microtubule_retraction_end_depolymerize_solve(system_graph, prev_state[i][j], system_graph_old, all_matches[i][j], settings);
    }
    }
}

template <typename MatchType, typename GraphType, typename BucketType, typename PropensityType, typename ParamType>
void microtubule_rule_firing(MatchType* all_matches, GraphType& system_graph, BucketType& bucket, PropensityType& prop, ParamType& settings)
{
    auto k = bucket.first;
    //std::vector<std::vector<mt_key_type>> good_matches;
    //std::size_t total_size = all_matches[0].size() + all_matches[2].size() + all_matches[4].size();
    //if(total_size == 0) return;
    
    double sums[8] = 
        {std::accumulate(prop[0].begin(), prop[0].end(), 0.0),
         std::accumulate(prop[1].begin(), prop[1].end(), 0.0),
         std::accumulate(prop[2].begin(), prop[2].end(), 0.0),
         std::accumulate(prop[3].begin(), prop[3].end(), 0.0),
         std::accumulate(prop[4].begin(), prop[4].end(), 0.0),
         std::accumulate(prop[5].begin(), prop[5].end(), 0.0),
         std::accumulate(prop[6].begin(), prop[6].end(), 0.0),
         std::accumulate(prop[7].begin(), prop[7].end(), 0.0)};
    
    double globSum = sums[0] + sums[1] + sums[2] + sums[3] + sums[4] + sums[5] + sums[6] + sums[7];

    if(globSum == 0.0) return;

    int ruleFired = 0; 
    auto sample = RandomRealsBetween(0.0, 1.0)();
    double progress = sample*globSum - sums[ruleFired];
    
    while(progress > 0.0)
    {
        ruleFired++;
        progress -= sums[ruleFired];
    }
    //std::cout << "Rule " << ruleFired << " has been selected to fire\n";

    if(ruleFired == 0)  //fire rule 0
    {
        //std::cout << "Sum 0: " << sums[ruleFired] << "\n";
        double local_sample = RandomRealsBetween(0.0, 1.0)();
        int eventFired = 0;
        double local_progress = local_sample*sums[ruleFired] - prop[ruleFired][eventFired];
        while(local_progress > 0.0) 
        {
            eventFired++;
            local_progress -= prop[ruleFired][eventFired];
        }
        microtubule_growing_end_polymerize_rewrite(system_graph, all_matches[0][eventFired], bucket); 
    }       
    if(ruleFired == 1)
    {
        //std::cout << "Sum 1: " << sums[ruleFired] << "\n";
        double local_sample = RandomRealsBetween(0.0, 1.0)();
        int eventFired = 0;
        double local_progress = local_sample*sums[ruleFired] - prop[ruleFired][eventFired];
        while(local_progress > 0.0) 
        {
            eventFired++;
            local_progress -= prop[ruleFired][eventFired];
        }

        microtubule_retraction_end_depolymerize_rewrite(system_graph, all_matches[2][eventFired], bucket);
    }
    if(ruleFired == 2)
    {
        //std::cout << "Sum 2: " << sums[ruleFired] << "\n";
        double local_sample = RandomRealsBetween(0.0, 1.0)();
        int eventFired = 0;
        double local_progress = local_sample*sums[ruleFired] - prop[ruleFired][eventFired];
        while(local_progress > 0.0) 
        {
            eventFired++;
            local_progress -= prop[ruleFired][eventFired];
        }
        microtubule_zippering_rewrite(system_graph, all_matches[4][eventFired], bucket, settings);
    }
    if(ruleFired == 3)
    {
        //std::cout << "Sum 3: " << sums[ruleFired] << "\n";
        double local_sample = RandomRealsBetween(0.0, 1.0)();
        int eventFired = 0;
        double local_progress = local_sample*sums[ruleFired] - prop[ruleFired][eventFired];
        while(local_progress > 0.0) 
        {
            eventFired++;
            local_progress -= prop[ruleFired][eventFired];
        }
        microtubule_positive_to_negative_rewrite(system_graph, all_matches[0][eventFired], bucket);
    }
    if(ruleFired == 4)
    {
        //std::cout << "Sum 4: " << sums[ruleFired] << "\n";
        double local_sample = RandomRealsBetween(0.0, 1.0)();
        int eventFired = 0;
        double local_progress = local_sample*sums[ruleFired] - prop[ruleFired][eventFired];
        while(local_progress > 0.0) 
        {
            eventFired++;
            local_progress -= prop[ruleFired][eventFired];
        }
        microtubule_negative_to_positive_rewrite(system_graph, all_matches[1][eventFired], bucket);
    }
    if(ruleFired == 5)
    {
        //std::cout << "Sum 5: " << sums[ruleFired] << "\n";
        double local_sample = RandomRealsBetween(0.0, 1.0)();
        int eventFired = 0;
        double local_progress = local_sample*sums[ruleFired] - prop[ruleFired][eventFired];
        while(local_progress > 0.0) 
        {
            eventFired++;
            local_progress -= prop[ruleFired][eventFired];
        }
        microtubule_crossover_rewrite(system_graph, all_matches[4][eventFired], bucket, settings);
    }
    if(ruleFired == 6)
    {
        //std::cout << "Sum 6: " << sums[ruleFired] << "\n";
        double local_sample = RandomRealsBetween(0.0, 1.0)();
        int eventFired = 0;
        double local_progress = local_sample*sums[ruleFired] - prop[ruleFired][eventFired];
        while(local_progress > 0.0) 
        {
            eventFired++;
            local_progress -= prop[ruleFired][eventFired];
        }
        microtubule_catastrophe_rewrite(system_graph, all_matches[4][eventFired], bucket, settings);
    }
    if(ruleFired == 7)
    {
        //std::cout << "Sum 7: " << sums[ruleFired] << "\n";
        double local_sample = RandomRealsBetween(0.0, 1.0)();
        int eventFired = 0;
        double local_progress = local_sample*sums[ruleFired] - prop[ruleFired][eventFired];
        while(local_progress > 0.0) 
        {
            eventFired++;
            local_progress -= prop[ruleFired][eventFired];
        }
        microtubule_katanin_severing_rewrite(system_graph, all_matches[3][eventFired], bucket, settings);
    }

}

} //end namespace plant 
} //end namespace cajete

#endif
