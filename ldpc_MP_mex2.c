1

#include <math.h>
#include "mex.h"


int   i,j,k,l,iter,sign, I1,I2, J;

typedef struct  {
    int     RowIdx,ColIdx;   // (i,j) = (row,col) where 1 exist
    float   inValue;  // lambda's
    float   outValue;  // Mu's
    
} constraint;


static void LDPC_MP(
        double  WorkState_out[], /*Store the Tanner Graph state */ 
        double  WorkState_in[], /* initial Work state to initialize the tanner gragh*/ 
        double	LLR_out[],  /* output  LLR - information bits */
        double	LLR_in[],  /* input  LLR - information bits */
        double  *Hm,  // partiy check matrix
        int		N,      /* sequence length */
        int		NK,     /* N-K : parity length */
        int     max_iter,  //  number of iteration
        int     NOnes //  number of ones in Hm
        )
{

    double  *p_Hm, *p_LLR_in, *p_LLR_out, *p_WorkState_out;
    float *llr_init; 
    int numberOfonesinH =  NOnes;
    llr_init   = mxMalloc(N*sizeof(float));
    for (j=0; j<N ; j++) (*(llr_init+j))= (float) *(LLR_in+j) ;
    p_LLR_out  =  mxMalloc(N*sizeof(double));  // allocate memory for LLR_out
    p_WorkState_out  =  mxMalloc(numberOfonesinH*sizeof(double));  // allocate memory for LLR_out
    p_LLR_out  =  LLR_in;  // initialization of p_LLR_out
    constraint *cstr; // create a pointer for constrainte
    cstr       =  mxMalloc(numberOfonesinH*sizeof(constraint));   // allocate memory for number none zero elements of matrix H
    p_Hm = Hm;
    iter = 0;
    
    p_WorkState_out = WorkState_in;
    
    do{
        /****************************  Start Operation in Check node level *******************************/
        for (i = 0; i<NK; ++i) {   // run over check nodes
            J=0;
            /*------------- Compute in values to check nodes from each code node (Lamda's) -------------*/
            for(j = 0; j<N; j++) {   // run over code nodes
                if( *(p_Hm+i+NK*j)){
                    (*(cstr+J)).ColIdx =  j;
                    (*(cstr+J)).RowIdx =  i;
                    (*(cstr+J)).inValue =  *(p_LLR_out+j) - *(p_WorkState_out+J);
                    J++;
                }
            }  
            /*------------- Compute out values to each code node from check node i (Mu's) -------------*/
            float minimum;
            
            for (I1=0; I1<J; I1++){
                sign = 1;
                minimum = 10^6;
                for(I2=0; I2<J; I2++){
                    if( I2!=I1) {
                        if( (*(cstr+I2)).inValue < 0) sign = sign*(-1);
                        if( fabsf((*(cstr+I2)).inValue) < minimum) minimum = fabsf((*(cstr+I2)).inValue);
                    }
                }
                (*(cstr+I1)).outValue = sign*minimum;
            }    
            cstr +=J;  // put the pointer in the next check node (next row)
            WorkState_out +=J;
        }
        /*************************** Start Operation in Code node level **********************************/
        cstr -= numberOfonesinH;  // initilize pointer to cstr
        WorkState_out -= numberOfonesinH;
        
        for (j=0; j<N; j++) {   // run over code nodes
            float result = (*(llr_init+j));
//             printf(" LLR_in= %f \n " , *(llr_init+j));
            for (i=0; i<numberOfonesinH; i++){
                if ( (*(cstr+i)).ColIdx == j) result += (*(cstr+i)).outValue;
            }
            (*(p_LLR_out+j)) = result; 
        } 
        /*************************** End of operations in Code node level *******************************/
        for (i=0; i<numberOfonesinH; i++)  (*(p_WorkState_out+i)) =  (*(cstr+i)).outValue;
        iter++; 
        
    }while(iter<max_iter);
    for (j=0; j<N; j++)   *(LLR_out+j) = *(p_LLR_out+j);
    for (i=0; i<numberOfonesinH; i++)  (*(WorkState_out+i)) =  (*(cstr+i)).outValue;
//   mxFree(cstr);   
//   mxFree(llr_init);
//   mxFree(p_LLR_out);
    // close procedure
    return;
}





/* Input Arguments */
#define LLR_in_IN	prhs[0]  /* LLRs  */
#define WorkState_in_IN	prhs[1]  /* Work state (initial Mu's)   */
#define Hmtrix_IN 	prhs[2]   /* Parity check matrix binary form  */
#define N_IN prhs[3]          /* code Length */
#define NK_IN prhs[4]          /* check Length */
#define max_iter_IN prhs[5]   /* max number of iterations  */
#define NOnes_IN prhs[6]   /* Number of ones in Hm  */

/* Output Arguments */
#define LLR_out_OUT	plhs[0]
#define WorkState_out_OUT plhs[1]
/*****************************************************/
void mexFunction( int nlhs,	 mxArray *plhs[], int nrhs, const mxArray *prhs[] )
{
    int N, NK,max_iter, NOnes ; // NH=columns of H, NH=rows of H,
    double	*LLR_out,*LLR_in, *Hm, *WorkState_out,*WorkState_in;
    
    NOnes = (int)(*mxGetPr(NOnes_IN));
    N = (int)(*mxGetPr(N_IN));
    NK = (int)(*mxGetPr(NK_IN));
    max_iter = (int)(*mxGetPr(max_iter_IN));
    LLR_out_OUT = mxCreateDoubleMatrix(N, 1, mxREAL);
    WorkState_out_OUT = mxCreateDoubleMatrix(NOnes,1, mxREAL);
    LLR_in = mxGetPr( LLR_in_IN ); 
    WorkState_in = mxGetPr( WorkState_in_IN ); 
    LLR_out = mxGetPr( LLR_out_OUT );
    WorkState_out = mxGetPr( WorkState_out_OUT  );
    Hm = mxGetPr(Hmtrix_IN);
    LDPC_MP(WorkState_out, WorkState_in, LLR_out, LLR_in, Hm,N,NK,max_iter,NOnes);
    
    return;
}
