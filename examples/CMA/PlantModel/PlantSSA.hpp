#ifndef DGGML_PLANT_SSA_HPP
#define DGGML_PLANT_SSA_HPP

#include "PlantTypes.hpp"
#include "PlantGrammar.hpp"

#include "ExpandedComplex2D.hpp"
#include "YAGL_Graph.hpp"
#include "YAGL_Algorithms.hpp"

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
namespace DGGML
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
        t_final = RCONST(user_data.settings.DELTA);

        flag = SUNContext_Create(NULL, &ctx);
        if(!DGGML::SundialsUtils::check_flag(&flag, "SUNContext_Create", 1))
            std::cout << "Passed the error check, suncontext created\n";
        
        // the number of equations is defined by the parameters of the nodes
        // involved in all the ode rule types
        num_eq = 0;
        for(const auto& r : user_data.rule_map[user_data.k])
        {
            const auto& inst = user_data.rule_system[r];
            if(inst.type == Rule::G || inst.type == Rule::R)
                num_eq += DIM3D;
        }
        num_eq++; //one extra equation for tau 

        std::cout << "*****NumEQ="<<num_eq<<"*****\n";
        y = NULL; 
        y = N_VNew_Serial(num_eq, ctx);
        
        auto j = 0;
        for(const auto& r : user_data.rule_map[user_data.k])
        {
            for(const auto& m : user_data.rule_system[r])
            {
                if(user_data.rule_system[r].type != Rule::G && user_data.rule_system[r].type != Rule::R)
                    continue;
                auto& node_data = user_data.system_graph.findNode(m)->second.getData();
                if(node_data.type == positive || node_data.type == negative)
                {
                    for(int i = 0; i < DIM3D; i++)
                        NV_Ith_S(y, j++) = node_data.position[i];
                } 
            }
        }
        NV_Ith_S(y, j) = user_data.tau;

        //NV_Ith_S(y, 0) = 50.0;
        //NV_Ith_S(y, 1) = 0.0;

        //std::cout << "J: " << j << "\n";
        arkode_mem = ERKStepCreate(rhs, t_start, y, ctx); 
        //if (!DGGML::SundialsUtils::check_flag((void *)arkode_mem, "ERKStepCreate", 1))
        //    std::cout << "Passed the error check, stepper initialized\n";

        reltol = 1.0e-4;//1.0e-6;
        abstol = 1.0e-8;//1.0e-10;

        flag = ERKStepSStolerances(arkode_mem, reltol, abstol);
        
        //set optional inputs
        flag = ERKStepSetUserData(arkode_mem, &user_data); 
    
        // specify the root finding function having num_roots  
        flag = ERKStepRootInit(arkode_mem, num_roots, root_func); 
        
        //ERKStepSetMinStep(arkode_mem, user_data.settings.DELTA_T_MIN/10.0);
        ERKStepSetMaxStep(arkode_mem, dt_out);//dt_out/10.0);
    }
    
    void step()
    {
        if(num_eq == 0) return;
        flag = ERKStepEvolve(arkode_mem, tout, y, &t, ARK_NORMAL);
        std::cout << "t: " << t << ", tout: " << tout << "\n";
        if(flag == ARK_ROOT_RETURN)
        {
            root_flag = ERKStepGetRootInfo(arkode_mem, roots_found);
            std::cout << "A root has been found\n";
        } 
        //successful solve
        else if(flag >= 0)
        {
            tout += dt_out;
            tout = (tout > t_final) ? t_final : tout;
        }
    }
    
    static int root_func(realtype t, N_Vector y, realtype* gout, void* user_data)
    {
        UserDataType* local_ptr = (UserDataType *)user_data;
        auto last = N_VGetLength(y) - 1;
        
        //check for when the tau equation crosses the sample 
        gout[0] = NV_Ith_S(y, last) - local_ptr->waiting_time;
        
        return 0;
    }

    static int rhs(realtype t, N_Vector y, N_Vector ydot, void* user_data) 
    {
        UserDataType* local_ptr = (UserDataType *)user_data;
        auto k = local_ptr->k;
        auto i = 0;
        //growth rule params
        auto v_plus = local_ptr->settings.V_PLUS;
        auto d_l = local_ptr->settings.DIV_LENGTH;
        
        //retraction rule params 
        auto l_d_f = local_ptr->settings.LENGTH_DIV_FACTOR;
        auto d_l_r = local_ptr->settings.DIV_LENGTH_RETRACT;
        auto v_minus = local_ptr->settings.V_MINUS;

        //NV_Ith_S(ydot, 0) = NV_Ith_S(y, 1);
        //NV_Ith_S(ydot, 1) = -9.8;
        
        //TODO: move the rules to their separate functions
        for(const auto& r : local_ptr->rule_map[k])
        {
            auto& instance = local_ptr->rule_system[r];
            if(instance.type == Rule::G)
            {
                auto id = instance[1];
                auto& node_i_data = local_ptr->system_graph.findNode(id)->second.getData();
                double l = 0.0;
                for(auto j = 0; j < 3; j++)
                {
                    double diff = NV_Ith_S(y, i+j) - node_i_data.position[j];
                    l += diff*diff;
                }
                l = sqrt(l);
                double length_limiter = (1.0 - (l/d_l));
                NV_Ith_S(ydot, i) = v_plus*node_i_data.unit_vec[0]*length_limiter;
                NV_Ith_S(ydot, i+1) = v_plus*node_i_data.unit_vec[1]*length_limiter;
                NV_Ith_S(ydot, i+2) = v_plus*node_i_data.unit_vec[2]*length_limiter;
            }
            if(local_ptr->rule_system[r].type == Rule::R)
            {
                auto id = instance[1];
                auto& node_i_data = local_ptr->system_graph.findNode(id)->second.getData();
                double l = 0.0;
                for(auto j = 0; j < 3; j++)
                {
                    double diff = NV_Ith_S(y, i+j) - node_i_data.position[j];
                    l += diff*diff;
                }
                l = sqrt(l);
                
                double length_limiter = l/d_l;
                if(length_limiter <= d_l_r) length_limiter = 0.0;

                for(int j = i; j < i+3; j++)
                    NV_Ith_S(ydot, j) = v_minus*node_i_data.unit_vec[j-i]*length_limiter;
                //std::cout << Rule::R << "\n";
            }
            i += 3;
        }
        
        NV_Ith_S(ydot, i) = local_ptr->geocell_propensity;

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
                if(user_data.rule_system[r].type != Rule::G && user_data.rule_system[r].type != Rule::R)
                    continue;
                auto& node_data = user_data.system_graph.findNode(m)->second.getData();
                if(node_data.type != positive && node_data.type != negative)
                    continue;
                //TODO: update the current velocity here?
                for(int i = 0; i < DIM3D; i++)
                    node_data.position[i] = NV_Ith_S(y, j++);
            } 
        }
        user_data.tau = NV_Ith_S(y, j);
    }

    //TODO: make this function reset the whole y vector
    //since it may change in size
    void reinit()
    {
        auto pre_eq = num_eq;
        auto num_eq = 0;

        for(const auto& r : user_data.rule_map[user_data.k])
        {
            const auto& inst = user_data.rule_system[r];
            if(inst.type == Rule::G || inst.type == Rule::R)
                num_eq += DIM3D;
        }
        num_eq++; //one extra equation for tau 

        std::cout << "*****NumEQ="<<num_eq<<"*****\n";
        N_VDestroy(y); 
        y = N_VNew_Serial(num_eq, ctx);
        
        auto j = 0;
        for(const auto& r : user_data.rule_map[user_data.k])
        {
            for(const auto& m : user_data.rule_system[r])
            {
                if(user_data.rule_system[r].type != Rule::G && user_data.rule_system[r].type != Rule::R)
                    continue;
                auto& node_data = user_data.system_graph.findNode(m)->second.getData();
                if(node_data.type == positive || node_data.type == negative)
                {
                    for(int i = 0; i < DIM3D; i++)
                        NV_Ith_S(y, j++) = node_data.position[i];
                } 
            }
        }
        
        //std::cout << num_eq << " " << pre_eq << "\n"; std::cin.get();
        user_data.tau = 0.0;
        NV_Ith_S(y, num_eq-1) = 0.0;
        ERKStepReInit(arkode_mem, rhs, t, y);
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
        std::cout << "Starting solver destruction\n";
        N_VDestroy(y);
        std::cout << "Destroyed y vector\n";
        ERKStepFree(&arkode_mem); //free the solver memory
        std::cout << "Destroyed arkode_mem\n";
        SUNContext_Free(&ctx); //always call prior to MPI_Finalize
        std::cout << "Destroyed the solver\n";
    
        /* Step 13: Finalilze MPI, if used */ 
    }
};

using RuleSystemType = DGGML::RuleSystem<Plant::mt_key_type>;
template <typename B, typename T, typename U, typename GeoplexType, typename GraphType, typename ParamType>
std::pair<double, double> plant_model_ssa_new(RuleSystemType& rule_system, B& k, T& rule_map, 
        U& anchor_list, GeoplexType& geoplex2D, GraphType& system_graph, 
        ParamType& settings, std::pair<double, double> geocell_progress)
{
    std::cout << "Cell " << k << " has " << rule_map[k].size() << " rules\n";
    std::cout << "{ ";
    for(auto& item : rule_map[k])
        std::cout << item << " ";
    std::cout << "}\n";
    double delta_t, geocell_propensity;
    std::size_t events, steps; 
    delta_t = 0.0; events = 0; steps = 0, geocell_propensity = 0.0;
    double tau = geocell_progress.first;
    double exp_sample = geocell_progress.second;
    if(exp_sample <= 0.0)
    {
        double uniform_sample = RandomRealsBetween(0.0, 1.0)();
        exp_sample = -log(1-uniform_sample);
        std::cout << "Warped wating time: " << exp_sample << "\n"; 
    }

    //need an aggregation
    typedef struct 
    {
        RuleSystemType& rule_system;
        GraphType& system_graph;
        B& k;
        T& rule_map;
        ParamType& settings;
        double& geocell_propensity;
        double& tau;
        double waiting_time;
    } PackType;
    
    //TODO: need a way to map the set of nodes and associated parameters 
    //from the rules to be solved into a flattened nvector, this is important 
    //for different rules that contribute to the same parameter 
    Solver ode_system(PackType{rule_system, system_graph, k, rule_map, settings, geocell_propensity, tau, exp_sample}, 2);
    
    std::vector<double> propensity_space;

    while(delta_t < settings.DELTA) 
    {
        //the propensities calculated for a particular rule
        std::vector<std::pair<std::size_t, double>> propensity_space;
            
        //while(delta_t < settings.DELTA && tau < exp_sample)
        if(tau < exp_sample)
        {
                    
            // STEP(1) : sum all of the rule propensities
            //geocell_propensity = RandomRealsBetween(0.0, 0.01)(); 
            auto rho = 0.0;
            for(const auto& r : rule_map[k])
            {
                auto& instance = rule_system[r];
                if(instance.type == Rule::G)
                {
                    propensity_space.push_back({ r,
                        100*microtubule_growing_end_polymerize_propensity(system_graph, 
                            instance.match, settings)}
                    );
                    rho += propensity_space.back().second;
                }
            }
            geocell_propensity = rho;
            //the step adapts based on propensity or systems fastest dynamic
            //settings.DELTA_T_MIN = 
            //std::min(1.0/(10.0*geocell_propensity), settings.DELTA_DELTA_T);
            
            // STEP(2) : solve the system of ODES 
            ode_system.step();
            ode_system.copy_back(); 
            // STEP(3) : use forward euler to solve the TAU ODE
            //tau += geocell_propensity*settings.DELTA_T_MIN; 
            
            // STEP(4) : advance the loop timer
            delta_t = ode_system.t;
            //delta_t += settings.DELTA_T_MIN; //TODO: make delta_t adaptive
            
            //std::cout << "tout: " << ode_system.t << " , DT: " << delta_t << "\n";
            
            steps++;
        }
    

        // if we get over our threshold an event can occur
        if(tau >= exp_sample) 
        {
            //determine which rule to fire and fire it
            std::cout << "Size of sample space: " << propensity_space.size() << "\n"; 
            
            auto local_sample = RandomRealsBetween(0.0, 1.0)();
            
            auto eventFired = 0;
            
            auto local_progress = local_sample*geocell_propensity - propensity_space[0].second;
            while(local_progress > 0.0) 
            {
                eventFired++;
                local_progress -= propensity_space[eventFired].second;
            }
            
            auto fired_id = propensity_space[eventFired].first;
            std::cout << "Selected rule id " << fired_id << " of type " << rule_system[fired_id].type << "\n"; 
            //need to update the cell list and the rule_map as well 
            
            ///*
            auto [invalidations, inducers] = 
                test_rewrite_growth(system_graph, rule_system[fired_id].match);
            std::cout << "Before Invalidations:\n";
            std::cout << "{ Rule Map Size: " << rule_map[k].size() << " }, { Anchor list size: " << anchor_list.size() 
                << " }, { Rule System Size: " << rule_system.size() << " }\n";

            //TODO: create a list of future invalidations?
            std::set<std::size_t> invalid_rules; // a set since some rules may be repeated
            std::set<std::size_t> future_invalidations; 
            //for each invalidated node 
            for(const auto& i : invalidations)
            {
                //find all the rules that node partcipates in
                const auto& rules = rule_system.inverse_index.find(i);
                //if the node participates in any rules
                if(rules != rule_system.inverse_index.end())
                {
                    //for each of the rules a node participates in
                    for(const auto& item : rules->second)
                    {
                        //search for the rule in the rule set for this cell
                        auto loc = std::find(rule_map[k].begin(), rule_map[k].end(), item);
                        //if we find the rule in the rule set for this cell
                        if(loc != rule_map[k].end())
                        {
                            //add the rule to the invalidation set
                            invalid_rules.insert(item);
                            std::cout << "Adding rule " << item << " to invalidation set\n";
                        }
                        //this should mean a lower dim cell is responsible 
                        //for the invalidation?
                        else
                        {
                            std::cout << "Rule " << item << " not found\n";
                            future_invalidations.insert(item);    
                        }
                    }
                }
            }
            std::cout << "Rules to invalidate now:   " << invalid_rules.size() << "\n";
            std::cout << "Rules to invalidate later: " << future_invalidations.size() << "\n";

            //only invalidate and erase rules a cell owns
            for(auto& item : invalid_rules)
            {
                rule_map[k].erase(std::find(rule_map[k].begin(), rule_map[k].end(), item));
                anchor_list.erase(rule_system[item].anchor);
                rule_system.invalidate_rule(item);
            }
            std::cout << "After invalidations:\n";
            std::cout << "{ Rule Map Size: " << rule_map[k].size() << " }, { Anchor list size: " << anchor_list.size() 
                << " }, { Rule System Size: " << rule_system.size() << " }\n";


            std::cout << "Find: ";
            for(const auto& i : inducers)
            {
                std::cout << i << " "; 
                anchor_list.insert({i, k});
            }
            std::cout << "\n";

            auto temp = inducers;
            for(const auto& i : temp)
            {
                auto res = YAGL::recursive_dfs(system_graph, i, 2);
                auto n_t = system_graph.findNode(i)->second.getData().type;
                if(n_t == negative) std::cout << "N ";
                if(n_t == positive) std::cout << "P ";
                if(n_t == intermediate) std::cout << "I ";

                std::cout << "Size at " << i << ": " << res.size() << "\n";
                for(auto& j : res)
                    inducers.insert(j);
            }

            for(auto& item : inducers)
            {
                std::cout << anchor_list[item] << " ";
                auto n_t = system_graph.findNode(item)->second.getData().type;
                if(n_t == negative)
                    std::cout << "N\n";
                else if(n_t == positive)
                    std::cout << "P\n";
                else if(n_t == intermediate)
                    std::cout << "I\n";
                else std::cout << "O\n";
            }

            std::cout << "Number of Inducers: " << inducers.size() << "\n";
            auto subgraph = YAGL::induced_subgraph(system_graph, inducers);
            std::cout << "Created induced subgraph\n";
            std::vector<DGGML::Plant::mt_key_type> sub_bucket;

            for(const auto& [key, value] : subgraph.getNodeSetRef())
                sub_bucket.push_back(key);
            std::cout << "SubBucket Size: " << sub_bucket.size() << "\n";
            auto incremental_matches = 
                DGGML::Plant::microtubule_growing_end_matcher(system_graph, sub_bucket);

            for(auto& item : incremental_matches)
            {
                //only add back in matches the cell actually owns
                //since the ones it didn't own weren't invalidated 
                //--idea so far: widen local search and then a sieve--
                bool valid = true;
                for(auto& v : item)
                    if(anchor_list[v] != k) { valid = false; std::cout << "Assigned cell mismatch rule G\n"; break; }
                if(!valid) continue;
                rule_system.push_back({std::move(item), DGGML::Rule::G});
                rule_map[k].push_back(rule_system.key_gen.current_key-1);
            }

            incremental_matches = 
                DGGML::Plant::microtubule_retraction_end_matcher(system_graph, sub_bucket);
            for(auto& item : incremental_matches)
            {
                //only add back in matches the cell actually owns
                bool valid = true;
                for(auto& v : item)
                    if(anchor_list[v] != k) { valid = false; std::cout << "Assigned cell mismatch rule R\n"; break; }
                if(!valid) continue;
                rule_system.push_back({std::move(item), DGGML::Rule::R});
                rule_map[k].push_back(rule_system.key_gen.current_key-1); 
            }
            std::cout << "After Reinsertion:\n";
            std::cout << "{ Rule Map Size: " << rule_map[k].size() << " }, { Anchor list size: " << anchor_list.size() 
                << " }, { Rule System Size: " << rule_system.size() << " }\n";



            //*/
            //zero out tau since a rule has fired 
            tau = 0.0;
            //reset the exp waiting_time sample 
            double uniform_sample = RandomRealsBetween(0.0, 1.0)();
            //sample the exponential variable
            exp_sample = -log(1-uniform_sample);
            ode_system.user_data.waiting_time = exp_sample; 
            std::cout << "Warped wating time: " << exp_sample << "\n"; 
            
            ode_system.print_stats();
            ode_system.reinit();
        }
    }
    std::cout << "Total steps taken: " << steps << "\n";
    return {tau, exp_sample};
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
    //if(!DGGML::SundialsUtils::check_flag(&flag, "SUNContext_Create", 1))
        //std::cout << "Passed the error check, suncontext created\n";
    
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
    //std::cout << "Freeing the suncontext\n";
    
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
                auto rho1 = settings.RHO_TEST_RATE*microtubule_zippering_propensity(system_graph, match, settings);
                auto rho2 = settings.RHO_TEST_RATE*microtubule_collision_crossover_propensity(system_graph, match, settings);
                auto rho3 = settings.RHO_TEST_RATE*microtubule_catastrophe_propensity(system_graph, match, settings);
                rule_propensities[2] += rho1;
                rule_propensities[5] += rho2;
                rule_propensities[6] += rho3;
                prop[2].push_back(rho1);
                prop[5].push_back(rho2);
                prop[6].push_back(rho3);

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
} //end namespace DGGML

#endif
