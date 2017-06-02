%%****************************************************
%%% llr_in : LLR input 
%%% HH     : sparse matrix
%%% THH    : transpose of matrix in sparse form ( sparse(full(HH)') )
%%% WS_in  : worK State input 
%%% MT     : max number of iterations
%%% llr_out: LLR output
%%% WS_out : work state output
%%*****************************************************

mex -v -largeArrayDims LDPC_MPmaxLog_mex.c  %generate mex file
[llr_out,WS_out] = LDPC_MPmaxLog_mex(llr_in,WS_in,HH,THH,MT);