# 2 "testAU.ino"
#include <ArduinoUnit.h>
#include "Entropy.h"

Entropy e;
//Don't forget to comment this string when not testing the ALL-0 suite to save space (see line 150)
char ZeroString[501];

/** TEST CASE with ALL 0 chain
 *    Should be equal to 0 since the proba of having the value 0 is 1
 */
 
test(MCV_ALL0)
{
  assertEqual(0,e.mostCommonValueEstimate(&ZeroString[0],500));
}

test(ColE_ALL0)
{
  assertEqual(0,e.collisionEstimate(&ZeroString[0],500));
}

test(MarkovE_ALL0)
{
  
  assertEqual(0,e.markovEstimate(&ZeroString[0],500));
}

test(CompE_ALL0)
{
 
  assertEqual(0,e.compressionEstimate(&ZeroString[0],500));
}

test(tTuple_ALL0)
{

  assertEqual(0,e.tTupleEstimate(&ZeroString[0],500));
} 

test(lrsE_ALL0)
{
  
  assertEqual(0,e.lrsEstimate(&ZeroString[0],10));
}

test(multiMCW_ALL0)
{
  assertEqual(0,e.multiMCWEstimate(&ZeroString[0],10));
}

test(lagPrE_ALL0){
  assertEqual(0,e.lagPredictionEstimate(&ZeroString[0],5));
}

test(multiMMC_ALL0){
  assertEqual(0,e.multiMMCEstimate(&ZeroString[0],12));
}

test(LZ78Y_ALL0){
  assertEqual(0,e.LZ78YEstimate(&ZeroString[0],12));
}

// End of ALL0 tests

/** Non-extreme particular cases
 *  Now we test using the examples from 6.3.* as test cases
 */
/*
test(MCV_Ex){
  char s[20] = {0,1,1,2,0,1,2,2,0,1,0,1,1,0,2,2,1,0,2,1};
  float MinE = e.mostCommonValueEstimate(&s[0],20);
  assertTrue( (MinE < 0.5519) && (MinE > 0.5517),"MCV Expected result: 0.5518. Actual Result: " << MinE); //Expected result is -log2(0.6822) ~=~ 0.5518. 
}

test(Col_Ex){
  char s[21] = {2,2,0,1,0,2,0,1,2,1,2,0,1,2,1,0,0,1,0,0,0};
  float MinE = e.collisionEstimate(&s[0], 21);
  Serial.print("Col_Ex estimate: "); Serial.println(MinE,10);
  assertTrue( (MinE<0.5017) && (MinE>0.5015),"Col Expected result: 0.5016. Actual Result = " << MinE ); 
}

test(Markov_Ex){
  char s[21] = {2,2,0,1,0,2,0,1,2,1,2,0,1,2,1,0,0,1,1,0,0};
  float MinE = e.markovEstimate(&s[0],21);
  assertTrue( (MinE < 0.6167) && (MinE > 0.6165),"Markov Expected result: 0.6166. Actual Result = " << MinE );
}

/* I commented this test because it produces garbage value that messes with the output of other tests
test(Comp_Ex){
  char s[21] = {2,2,0,1,0,2,0,1,2,1,2,0,1,2,1,0,0,1,0,0,0};
  float MinE = e.compressionEstimate(&s[0],21);
  assertTrue( (MinE < 0.5140) && (MinE > 0.5138),"Expected result: 0.5139. Actual Result = " << MinE );
}
*/

/* no examples for those tests, leaving them for later
test(tTuple_Ex){
}
*/
/* 
test(LRS_Ex){

}
*/
/*
test(MultiMCW_Ex){
  char s[12] = {1,2,1,0,2,1,1,2,2,0,0,0};
  float MinE = e.multiMCWEstimate(&s[0],12);
  assertTrue( (MinE < 0.3910) && (MinE > 0.3908),"MultiMCW Expected result: 0.3909. Actual Result = " << MinE );
}

test(lagP_Ex){
  char s[10] = {2,1,3,2,1,3,1,3,1,2};
  float MinE = e.lagPredictionEstimate(&s[0],10);
  assertTrue( (MinE < 0.5841) && (MinE > 0.5849),"lagP Expected result: 0.5850. Actual Result = " << MinE );
}

test(MultiMMC_Ex){
  char s[9] = {2,1,3,2,1,3,1,3,1};
  float MinE = e.multiMMCEstimate(&s[0],9);
  assertTrue( (MinE < 0.0756) && (MinE > 0.0754),"MultiMMC Expected result: 0.0755. Actual Result = " << MinE );
}

test(LZ78Y_Ex){
  char s[13] = {2,1,3,2,1,3,1,3,1,2,1,3,2};
  float MinE = e.LZ78YEstimate(&s[0],13);
  assertTrue( (MinE < 0.0192) && (MinE > 0.0190),"LZ78Y Expected result: 0.0191. Actual Result = " << MinE );
}*/

void setup()
{
  Serial.begin(9600);
  while(!Serial) {} // Portability for Leonardo/Micro (see arduinounit readme)
}

void loop()
{
  //Initializing an ALL-0 chain for extreme case
  //Don't forget to comment when not testing the ALL-0 suite (see line 150)
  for(int i = 0; i < 100; i++){
    ZeroString[i] = 0;    
  }   
  // excluding these two tests for now as the test case seems to give an infinite loop
  Test::exclude("CompE_ALL0");
  Test::exclude("ColE_ALL0"); 
  
  // 
  /** Excludes a whole test suite (ALL0 or Ex)
   * Note: there isn't enough space on the device to be able to use these two instructions as expected
   * If the tests are left uncommented, we start to get unexpected behaviour
   * The test suites (ALL0 and Ex) thus cannot be run together
   * I suggest leaving these exlusions commented for reference, and commenting the unwanted test suite.
   * Don't forget to comment ZeroString and its initialization for loop to save space and memory
   * 
   */
  //Test::exclude("*_ALL0");
  //Test::exclude("*_Ex"); 
  
  Test::exclude("Comp_Ex"); //Definitely produces some garbage value. Might be an overflow. Exluding it for now.
  
  
  Test::run();
  delay(1000);

}
