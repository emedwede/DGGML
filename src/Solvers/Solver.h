#ifndef DGGML_SOLVER_H
#define DGGML_SOLVER_H

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

#include "SundialsUtils.hpp"
#include "HelperFunctions.hpp"

// TODO: Need to generalize the sovler code and we need to add in the root finding fix from
// CajeteCMA repo
namespace DGGML {

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
                : user_data(_user_data), num_eq(_num_eq)
        {
            t_start = RCONST(0.0);
            dt_out = RCONST(user_data.model->settings.DELTA_T_MIN);
            t = t_start;
            tout = t_start + dt_out;
            t_final = RCONST(user_data.model->settings.DELTA);

            flag = SUNContext_Create(NULL, &ctx);
            if(!DGGML::SundialsUtils::check_flag(&flag, "SUNContext_Create", 1))
                std::cout << "Passed the error check, suncontext created\n";

            std::cout << "*****NumEQ="<<num_eq<<"*****\n";
            y = NULL;
            y = N_VNew_Serial(num_eq, ctx);

            std::size_t offset = 0;
            //set the initial conditions
            for(auto& item : user_data.solving_rules)
            {
                auto& inst = user_data.rule_matches[item];
                auto& name = inst.name;
                auto& ode = user_data.grammar_analysis.solving_rules.find(name)->second;
                //TODO: find a more efficient way
                //need to construct a graph from the match keys and
                //the mapping from the orignal LHS to the match
                auto induced_graph = induce_from_set(inst, user_data.component_matches, user_data.model->system_graph);
                std::map<std::size_t, std::size_t> lhs_vertex_map;
                construct_grammar_match_map(inst, user_data.grammar_analysis, lhs_vertex_map, user_data.component_matches);
                std::cout << "for rule " << item << " we have a map: ";
                for(auto& [k, v] : lhs_vertex_map)
                    std::cout << "{ " << k << " -> " << v << " } ";
                std::cout << "\n";
                ode.ic(induced_graph, lhs_vertex_map, y, offset);
                offset += ode.num_eq;
            }
            NV_Ith_S(y, offset) = user_data.tau;

            for(auto i = 0; i < num_eq; i++)
                std::cout << NV_Ith_S(y, i) << "\n";

//            arkode_mem = ERKStepCreate(rhs, t_start, y, ctx);
//            //if (!DGGML::SundialsUtils::check_flag((void *)arkode_mem, "ERKStepCreate", 1))
//            //    std::cout << "Passed the error check, stepper initialized\n";
//
//            reltol = 1.0e-4;//1.0e-6;
//            abstol = 1.0e-8;//1.0e-10;
//
//            flag = ERKStepSStolerances(arkode_mem, reltol, abstol);
//
//            //set optional inputs
//            flag = ERKStepSetUserData(arkode_mem, &user_data);
//
//            // specify the root finding function having num_roots
//            flag = ERKStepRootInit(arkode_mem, num_roots, root_func);
//
//            //ERKStepSetMinStep(arkode_mem, user_data.settings.DELTA_T_MIN/10.0);
//            ERKStepSetMaxStep(arkode_mem, dt_out);//dt_out/10.0);
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
                    //if(node_data.type != positive && node_data.type != negative)
                    //    continue;
                    //TODO: update the current velocity here?
                    for(int i = 0; i < 3; i++)
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
                    num_eq += 3;
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
//                    if(node_data.type == positive || node_data.type == negative)
//                    {
//                        for(int i = 0; i < 3; i++)
//                            NV_Ith_S(y, j++) = node_data.position[i];
//                    }
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
//            std::cout << "Starting solver destruction\n";
//            N_VDestroy(y);
//            std::cout << "Destroyed y vector\n";
//            ERKStepFree(&arkode_mem); //free the solver memory
//            std::cout << "Destroyed arkode_mem\n";
//            SUNContext_Free(&ctx); //always call prior to MPI_Finalize
//            std::cout << "Destroyed the solver\n";

            /* Step 13: Finalilze MPI, if used */
        }
    };
} // DGGML

#endif //DGGML_SOLVER_H
