#include "Entropy.h"

/*
void Entropy::convertToArray(String input, char* output){

}
*/
float Entropy::allTests(char* input, size_t len){

  return 0;
}

//individual tests
float Entropy::mostCommonValueEstimate(char* input, size_t len){

  return 0;
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
