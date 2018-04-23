#include "Arduino.h"
#include "utils.h"
#include "Entropy.h"
#include <algorithm>
using namespace std;
void print_usage(){
	cout << "Usage is: ./non_iid_main.out <file_name> <bits_per_word> [-v]" << endl;
	cout << "\t <file_name>: Must be relative path to a binary file with at least 1 million entries (words)." << endl;
	cout << "\t <bits_per_word>: Must be between 1-8, inclusive." << endl;
	cout << "\t -v: Optional verbosity flag for more output." << endl;
	cout << endl;
}
int main(int argc, char** argv){
  Entropy e;

  char* dataset = (char*)calloc(SIZE,sizeof(char));
  int word_size = 8;
  bool verbose = false;
  char* file_path;

  // Parse args
  if(argc != 3 && argc != 4){
    cout << "Incorrect usage." << endl;
    print_usage();
    exit(-1);
  }else{

    // Gather filename
    file_path = argv[1];

    // Gather bits per word
    word_size = atoi(argv[2]);
    if(word_size < 1 || word_size > 8){
      cout << "Invalid bits per word." << endl;
      print_usage();
      exit(-1);
    }

    const char verbose_flag = 'v';
    if(argc == 4){
      verbose = (argv[3][1] == verbose_flag);
    }
  }

  if(verbose){
    cout << "Opening file: " << file_path << endl;
  }

  if(!read_file(file_path, dataset, word_size)){
    cout << "Error reading file. Need 1 million entries for the tests to work." << endl;
    print_usage();
    exit(-1);
  }

  // The maximum min-entropy is -log2(1/SIZE)
  double min_entropy = log2(SIZE);
  double H_min;

  if(verbose){
    cout << endl << "Running non-IID tests..." << endl;
  }

  cout << endl << "Entropic statistic estimates:" << endl;

  // Section 6.3.1 - Estimate entropy with Most Common Value
  H_min = e.mostCommonValueEstimate(dataset,SIZE);
  if(verbose){
    cout << "Most Common Value Estimate = " << H_min << endl;
  }
  min_entropy = min(min_entropy, H_min);

  // Section 6.3.2 - Estimate entropy with Collision Test
  //this segfaults because of a variable length array that is almost certainly not supported in c++
  H_min = e.collisionEstimate(dataset,SIZE);
  if(verbose){
    cout << "Collision Test Estimate = " << H_min << endl;
  }
  min_entropy = min(min_entropy, H_min);

  // Section 6.3.3 - Estimate entropy with Markov Test
  if(word_size > 6){
    char mapped_down[SIZE];
    for(int i = 0; i < SIZE; i++){
      mapped_down[i] = dataset[i] & 63;
    }

    H_min = e.markovEstimate(mapped_down, SIZE);
  }else{
    H_min = e.markovEstimate(dataset,SIZE);
  }
  if(verbose){
    cout << "Markov Test Estimate = " << H_min << endl;
  }
  min_entropy = min(min_entropy, H_min);

  // Section 6.3.4 - Estimate entropy with Compression Test
  H_min = e.compressionEstimate(dataset,SIZE);
  if(verbose){
    cout << "Compression Test Estimate = " << H_min << endl;
  }
  min_entropy = min(min_entropy, H_min);

  // Section 6.3.5 - Estimate entropy with t-Tuple Test
  H_min = e.tTupleEstimate(dataset,SIZE);
  if(verbose){
    cout << "t-Tuple Test Estimate = " << H_min << endl;
  }
  min_entropy = min(min_entropy, H_min);

  // Section 6.3.6 - Estimate entropy with Longest Repeated Substring Test (LRS)
  // Slow, needs speedup, wrap LRS functions into a single call
  H_min = e.lrsEstimate(dataset,SIZE);
  if(verbose){
    cout << "Longest Reapeated Substring Test Estimate = " << H_min << endl;
  }
  min_entropy = min(min_entropy, H_min);

  cout << endl << "Predictor estimates:" << endl;

  // Section 6.3.7 - Estimate entropy with Multi Most Common in Window Test
  H_min = e.multiMCWEstimate(dataset, SIZE);
  if(verbose){
    cout << "Multi Most Common in Window (Multi MCW) Test = " << H_min << endl;
  }
  min_entropy = min(min_entropy, H_min);

  // Section 6.3.8 - Estimate entropy with Lag Prediction Test
  H_min = e.lagPredictionEstimate(dataset, SIZE);
  if(verbose){
    cout << "Lag Prediction Test = " << H_min << endl;
  }
  min_entropy = min(min_entropy, H_min);

  // Section 6.3.9 - Estimate entropy with Multi Markov Model with Counting Test (MultiMMC)
  H_min = e.multiMMCEstimate(dataset, SIZE);
  if(verbose){
    cout << "Multi Markov Model with Counting (MultiMMC) Prediction Test = " << H_min << endl;
  }
  min_entropy = min(min_entropy, H_min);

  // Section 6.3.10 - Estimate entropy with LZ78Y Test
  // Not super fast, just a touch longer than the python but with -O3 is super fast
  H_min = e.LZ78YEstimate(dataset, SIZE);
  free(dataset);
  if(verbose){
    cout << "LZ78Y Prediction Test = " << H_min << endl;
  }
  min_entropy = min(min_entropy, H_min);

  cout << "Min Entropy: " << min_entropy << endl;
}
