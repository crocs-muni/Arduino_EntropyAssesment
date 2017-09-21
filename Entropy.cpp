#include "Entropy.h"
#include <string.h>
#include <math.h>

#define MIN(a,b) (((a)<(b))?(a):(b))


/*
void Entropy::convertToArray(String input, char* output){

}
*/
float Entropy::allTests(char* input, size_t len){

  return 0;
}

//individual tests
float Entropy::mostCommonValueEstimate(char* input, size_t len){

  char data[len];
  memcpy(data, input, len*sizeof(char));

  //sort for easy counting
  sortArray(data, len);

  char current = data[0];
  int countCurrent = 1;
  int totalCount = 1;
  int cmax = 1;
  for(int i = 1; i < len; i++){
    if(data[i] == current){
      countCurrent = countCurrent + 1;
      if(countCurrent > cmax){
        cmax = countCurrent;
      }
    } else {
      current = data[i];
      totalCount = totalCount + 1;
      countCurrent = 1;
    }
  }
  float pmax = (cmax * 1.0) / (len * 1.0);
  float ubound = pmax + 2.576 * sqrt(pmax * (1.0-pmax)/(len * 1.0));
  float pu = MIN(1,ubound);
  float entropy = - log2( pu );
  return entropy;
}

float Entropy::collisionEstimate(char* input, size_t len){

  return 0;
}

float Entropy::markovEstimate(char* input, size_t len){

  return 0;
}

float Entropy::compressionEstimate(char* input, size_t len){

  return 0;
}

float Entropy::tTupleEstimate(char* input, size_t len){

  return 0;
}

float Entropy::lrsEstimate(char* input, size_t len){

  return 0;
}

float Entropy::log2( float n ){  
  // log(n)/log(2) is log2.
  return log( n ) / log( 2 );
}


void Entropy::sortArray(char* input, size_t len){
  if (len < 2) return;

  char pivot = input[len / 2];

  int i, j;
  for (i = 0, j = len - 1; ; i++, j--) {
    while (input[i] < pivot) i++;
    while (input[j] > pivot) j--;

    if (i >= j) break;

    int temp = input[i];
    input[i] = input[j];
    input[j] = temp;
  }

  sortArray(input, i);
  sortArray(input + i, len - i);
}
