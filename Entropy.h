#ifndef Entropy_h
#define Entropy_h

//Max alphabet size should be 256

class Entropy {

public:
  //convenience functions for input formatting
  //maybe not necesary, depends on expected input
  void convertToArray(char* input, char* output);

  //run all tests and output the lowest value (safe estimate for min-entropy)
  float allTests(char* input);

  //individual tests
  float mostCommonValueEstimate(char* input);
  float collisionEstimate(char* input);
  float markovEstimate(char* input);
  float compressionEstimate(char* input);
  float tTupleEstimate(char* input);
  float lrsEstimate(char* input);

private:
  //number of bits used to represent sample output values
  int bits_per_symbol = 8;

}


#endif
