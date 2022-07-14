#include <iostream>

#include "catch.hpp"
#include <arkode/arkode_erkstep.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

static int check_flag(void *flagvalue, const char *funcname, int opt)
{
    int *errflag;

    /* Check if SUNDIALS function returned NULL pointer - no memory allocated */
    if (opt == 0 && flagvalue == NULL) {
      fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed - returned NULL pointer\n\n",
              funcname);
      return 1; }

    /* Check if flag < 0 */
    else if (opt == 1) {
      errflag = (int *) flagvalue;
      if (*errflag < 0) {
        fprintf(stderr, "\nSUNDIALS_ERROR: %s() failed with flag = %d\n\n",
                funcname, *errflag);
        return 1; }}

    /* Check if function returned NULL pointer - no memory allocated */
    else if (opt == 2 && flagvalue == NULL) {
      fprintf(stderr, "\nMEMORY_ERROR: %s() failed - returned NULL pointer\n\n",
              funcname);
      return 1; }

    return 0;
}


TEST_CASE("Running SUNDIALS ARKODE Test", "[ARKODE test]")
{    
    std::cout << "Running SUNDIALS ARKODE test\n";
    
    int flag; //general resuable error flag 

    N_Vector y = NULL; //empty vector for storing the solution
    void* arkode_mem = NULL; //empty arkode memory structure 
    realtype t, tout;
    int NEQ = 10; //number of equations

    //create the sundials suncontext
    SUNContext ctx;
    flag = SUNContext_Create(NULL, &ctx);
    REQUIRE(check_flag(&flag, "SUNContext_Create", 1) == false);
    std::cout << "SUNContext created\n";

    //initialize the data structures 
    y = N_VNew_Serial(NEQ, ctx); //number of equations and owning suncontext 
    REQUIRE(check_flag((void *)y, "N_VNew_Serial", 0) == false);
    std::cout << "NVector created\n";

    //free the vector
    N_VDestroy(y);
    std::cout << "NVector destroyed\n";

    //free the sundcontext 
    flag = SUNContext_Free(&ctx);
    REQUIRE((check_flag(&flag, "SUNContext_Free", 1)) == false);
    std::cout << "SUNContext destroyed\n";
}
