
 /* LDPC maxlog decoder: A.Benyouss ver1.2017  */

#include <stdio.h>
#include <math.h>
#include "mex.h"

int i,j,l,n,l,k,sign,iter,temp;
double minimum;

typedef struct  {
    double   inValue;  // lambda's
    double   outValue;  // Mu's
    double   summ;  // sum of values
    
} constraint;

static void LDPC_MP(
        double  WorkState_in[], /*the Tanner Graph initial state */
        double  WorkState_out[], /*the Tanner Graph final state */
        double	LLR_in[],  /* intput  LLR - coded bits */
        double	LLR_out[],  /* output  LLR - doded bits */
        mwIndex *JJC,  /* pointer to columns o the sparse matrix */
        mwIndex *TJJC,  /* pointer to columns o the sparse matrix transpse */
        mwIndex	*IIR,  /* pointer to the rows of sparse matrix */
        int     max_iter,  // max number of iteration
        mwIndex nnz, //  number of none zero elements
        mwSize  columns, //  number columns of the sparse matrix
        mwSize  Tcolumns, //  number columns of the sparse matrix transpose
        int      N
        )
{
    int Tindx[nnz];
    double *p_WorkState_out;
    constraint *cstr; // create a pointer for constraint
    cstr   =  mxMalloc(nnz*sizeof(constraint));   // allocate memory for number none zero elements of matrix H
    p_WorkState_out   =  mxMalloc(nnz*sizeof(double));  
    
    p_WorkState_out = WorkState_in;
    for( i=0; i<nnz; i++) (*(cstr + i)).outValue = *(p_WorkState_out+i);
    
    int indx = 0;
    for(l = 0; l<Tcolumns; l++) {
        for( i=0; i<nnz; i++){
            if ( (*(IIR+i)) == l){  // (*(IIR+i)) = (*(cstr + i)).Rdx
                Tindx[indx] = i;
                indx++;
            }
        }
    }
    
    
    iter = 0;
    
    do{
           for (j=0; j<N ; j++){
            temp = *(JJC+j+1) - (*(JJC+j));
            float result = (*(LLR_in+j));
            for (l=0; l<temp; l++) result += (*(cstr + l + (*(JJC+j)))).outValue;
            for (l=0; l<temp; l++) (*(cstr + l + (*(JJC+j)))).summ = result;
        }
        
        for (n=0; n<nnz; n++) (*(cstr+n)).inValue =  (*(cstr+n)).summ - (*(cstr + n)).outValue;
        
        for(j = 0; j<Tcolumns; j++){
            
            temp = *(TJJC+j+1) - (*(TJJC+j));
            for (l=0; l<temp; l++){
                mwIndex I = l + (*(TJJC+j));
                minimum = 10^6;
                sign = 1;
                for (k=0; k<temp; k++){
                    mwIndex J = k + (*(TJJC+j));
                    if ( k != l){
                        if( (*(cstr + Tindx[J] )).inValue < 0) sign = sign*(-1);
                        if( fabs((*(cstr + Tindx[J])).inValue) < minimum) minimum = fabs((*(cstr+Tindx[J] )).inValue);
                    }
                }
                (*(cstr + Tindx[I] )).outValue = sign*minimum;
                (*(WorkState_out + Tindx[I])) = (*(cstr + Tindx[I])).outValue;
            }
            
        }
        
        iter++;
    }  while(iter < max_iter);
     
         for (j=0; j<N ; j++){
            temp = *(JJC+j+1) - (*(JJC+j));
            float result = (*(LLR_in+j));
            for (l=0; l<temp; l++) result += (*(cstr + l + (*(JJC+j)))).outValue;
            (*(LLR_out+j)) = result;
        }
   
    return;
 
    mxFree(cstr);
    mxFree(p_WorkState_out);
    
}



/* Input Arguments */
#define LLR_in_IN	prhs[0]  /* LLRs  */
#define WorkState_in_IN	prhs[1]  /* Work state (initial Mu's)   */
#define Hmtrix_IN 	prhs[2]   /* Parity check matrix sparse form  */
#define THmtrix_IN 	prhs[3]   /* Transpose of Parity check matrix sparse form  */
#define max_iter_IN prhs[4]   /* max number of iterations  */


/* Output Arguments */
#define LLR_out_OUT	plhs[0]
#define WorkState_out_OUT plhs[1]
/*****************************************************/
void mexFunction( int nlhs,	 mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    int N, max_iter ; // columns
    double	*LLR_out,*LLR_in, *WorkState_out,*WorkState_in;
    mwSize  columns, Tcolumns ;
    mwIndex nnz;
    mwIndex *JJC, *IIR, *TJJC;
    
    N = mxGetN(LLR_in_IN);
    columns = mxGetN(prhs[2]);
    IIR = mxGetIr(prhs[2] );
    JJC = mxGetJc(prhs[2] );
    nnz = *(mxGetJc(prhs[2]) + columns);
    
    TJJC = mxGetJc(prhs[3] );
    Tcolumns = mxGetN(prhs[3]);
    
    
    max_iter = (int)(*mxGetPr(max_iter_IN));
    LLR_out_OUT = mxCreateDoubleMatrix(N, 1, mxREAL);
    LLR_out = mxGetPr( LLR_out_OUT );
    WorkState_out_OUT = mxCreateDoubleMatrix(nnz,1, mxREAL);
    WorkState_out = mxGetPr( WorkState_out_OUT  );
    LLR_in = mxGetPr( LLR_in_IN );
    WorkState_in = mxGetPr( WorkState_in_IN );
    
    LDPC_MP(WorkState_in,WorkState_out,LLR_in,LLR_out,JJC,TJJC,IIR,max_iter,nnz,columns,Tcolumns,N);
    
    
    return;
}
