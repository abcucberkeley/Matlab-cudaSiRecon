#include "cudaSireconImpl.h"
#include "SIM_reconstructor.hpp"


void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{
    try {
        int argc = nrhs-2;
        char** argv = (char**)mxCalloc(argc, sizeof(char*));
        for (int i = 0; i < argc; i++){
            /*
             * if( !mxIsChar( prhs[i] ) ){
             * mexErrMsgTxt("Input must be of type char.");
             * return;
             * }
             */
            argv[i] = mxArrayToString(prhs[i]);
        }
        
        std::cout << "Parameters converted from MATLAB" << std::endl;
        SIM_Reconstructor myreconstructor(argc, argv, prhs[nrhs-2], prhs[nrhs-1]);
        //mexErrMsgIdAndTxt("zarr:threadError","Testing Outside");
        
        /*
         * for (int it = 0; it < myreconstructor.getNTimes(); ++it) {
         * for (int iw = 0; iw < 1; ++iw) {
         * myreconstructor.setFile(it, iw);
         * myreconstructor.loadAndRescaleImage(it, iw);
         * myreconstructor.setCurTimeIdx(it);
         * if (myreconstructor.processOneVolume())
         * myreconstructor.writeResult(it, iw);
         * }
         * }
         */
        myreconstructor.setFile(0, 0);
        myreconstructor.loadAndRescaleImage(0, 0);
        myreconstructor.setCurTimeIdx(0);
        if (myreconstructor.processOneVolume()){
            plhs[0] = myreconstructor.passResult(0, 0);
        }
        
        
        myreconstructor.closeFiles();
    }
    catch (std::exception &e) {
        std::cerr << "\n!!Error occurred: " << e.what() << std::endl;
        mexErrMsgIdAndTxt("cudsirecon:mainError","An error occured");
        //return 0;
    }
//return 0;
}
