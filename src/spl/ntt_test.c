/* sporbiter generated code 
   created 01/12/21*/ 
void ntt_8(long* X, long* Y){
  /*no locals */ 
  Y[0]=X[0]+X[1]+X[2]+X[3]+X[4]+X[5]+X[6]+X[7];
  Y[1]=X[0]+9*X[1]+13*X[2]+15*X[3]+16*X[4]+8*X[5]+4*X[6]+2*X[7];
  Y[2]=X[0]+13*X[1]+16*X[2]+4*X[3]+X[4]+13*X[5]+16*X[6]+4*X[7];
  Y[3]=X[0]+15*X[1]+4*X[2]+9*X[3]+16*X[4]+2*X[5]+13*X[6]+8*X[7];
  Y[4]=X[0]+16*X[1]+X[2]+16*X[3]+X[4]+16*X[5]+X[6]+16*X[7];
  Y[5]=X[0]+8*X[1]+13*X[2]+2*X[3]+16*X[4]+9*X[5]+4*X[6]+15*X[7];
  Y[6]=X[0]+4*X[1]+16*X[2]+13*X[3]+X[4]+4*X[5]+16*X[6]+13*X[7];
  Y[7]=X[0]+2*X[1]+4*X[2]+8*X[3]+16*X[4]+15*X[5]+13*X[6]+9*X[7];
}
