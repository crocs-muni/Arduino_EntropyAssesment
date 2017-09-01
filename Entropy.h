#ifndef Entropy_h
#define Entropy_h

#include <cstddef>


//Max alphabet size should be 256

class Entropy {

public:
  //convenience functions for input formatting
  //maybe not necesary, depends on expected input
  //void convertToArray(String input, char* output, size_t* len);

  //run all tests and output the lowest value (safe estimate for min-entropy)
  float allTests(char* input, size_t len);

  //individual tests
  float mostCommonValueEstimate(char* input, size_t len);
  float collisionEstimate(char* input, size_t len);
  float markovEstimate(char* input, size_t len);
  float compressionEstimate(char* input, size_t len);
  float tTupleEstimate(char* input, size_t len);
  float lrsEstimate(char* input, size_t len);

private:
  //number of bits used to represent sample output values
  int bits_per_symbol = 8;
  void sortArray(char* input, size_t len);

};


#endif
