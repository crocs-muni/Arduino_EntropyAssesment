#ifndef Entropy_h
#define Entropy_h

#include "Arduino.h"
#include <algorithm>
//Shared constants:
#define NB_SYM 2
#define SIZE 1000000


//Constants for MultiMCW
//Sizes for MultiMCW window
#define W_ONE 63
#define W_TWO 255
#define W_THREE 1023
#define W_FOUR 4096

//Constants for Compression Estimate:
//Number of observations 
#define NB_OBS 1000

//Constants for LZ78Y estimate
#define MAX_DICT_SIZE  65536
#define B 16

//Max alphabet size should be 256
class Entropy {

public:
  //convenience functions for input formatting
  //maybe not necesary, depends on expected input
  //void convertToArray(String input, char* output, size_t* len);

  //run all tests and output the lowest value (safe estimate for min-entropy)
  float allTests(char* input, size_t len);

  //individual tests
  //6.3.1
  float mostCommonValueEstimate(char* input, size_t len);
  //6.3.2
  float collisionEstimate(char* input, size_t len);
  //6.3.3
  float markovEstimate(char* input, size_t len);
  //6.3.4
  float compressionEstimate(char* input, size_t len);
  //6.3.5
  float tTupleEstimate(char* input, size_t len);
  //6.3.6
  float lrsEstimate(char* input, size_t len);
  //6.3.7
  float multiMCWEstimate(char* input, size_t len);
  //6.3.8
  float lagPredictionEstimate(char* input, size_t len);
  //6.3.9
  float multiMMCEstimate(char* input, size_t len);
  //6.3.10
  float LZ78YEstimate(char* input, size_t len);
  
private:
  //number of bits used to represent sample output values
  //int bits_per_symbol = 8;
  void sortArray(char* input, size_t len);
  float log2( float n );  
  // function for solving for P
  size_t solveForP(float mu_bar,float* p);
  // function for solving for P
  float calcEpS(float p_c);
  
  //Functions for 6.3.4
  int solve_for_p_634(float mu_bar, size_t n, size_t v, size_t d, float *p);
  float EppM(float p, size_t n, size_t v, size_t d);
  float func_G(float p, size_t v, size_t d);

  //Functions for 6.3.7-6.3.10
  char Entropy:: mostComman(char* data, size_t length);
  size_t Entropy::findMaxRun(size_t *correct, size_t N);
  float Entropy::calcRun(size_t r, size_t N);
  float Entropy::calc_qn(float p, size_t r, size_t n);
  float Entropy::find_root(float p, size_t r);
};


#endif
