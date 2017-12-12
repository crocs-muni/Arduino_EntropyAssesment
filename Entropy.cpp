#include "Entropy.h"
#include <string.h>
#include <math.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


/*
  void Entropy::convertToArray(String input, char* output){

  }
*/

float Entropy::allTests(char* input, size_t len) {

  return 0;
}

//individual tests
float Entropy::mostCommonValueEstimate(char* input, size_t len) {

  char data[len];
  memcpy(data, input, len * sizeof(char));

  //sort for easy counting
  sortArray(data, len);

  char current = data[0];
  int countCurrent = 1;
  int totalCount = 1;
  int cmax = 1;
  for (int i = 1; i < len; i++) {
    if (data[i] == current) {
      countCurrent = countCurrent + 1;
      if (countCurrent > cmax) {
        cmax = countCurrent;
      }
    } else {
      current = data[i];
      totalCount = totalCount + 1;
      countCurrent = 1;
    }
  }
  float pmax = (cmax * 1.0) / (len * 1.0);
  float ubound = pmax + 2.576 * sqrt(pmax * (1.0 - pmax) / (len * 1.0));
  float pu = MIN(1, ubound);
  float entropy = - log2( pu );
  return entropy;
}

float Entropy::collisionEstimate(char* input, size_t len) {

  char data[len];
  memcpy(data, input, len * sizeof(char));

  //void collistionEstimate(char *input, size_t len,  float *p, float *min_entropy)
  //{

  /*1. Set v = 1, index = 1
    2. Beginning with s_index, step through the dataset until any observed
     value is repeated; i.e., find smallest j s.t.si = sj for some
     i with 1 <= i < j
    3. Set t_v = j - index + 1, v = v + 1, and index = j + 1
    4. Repeat steps 2 - 3 until the end of the dataset is reached
    5. set v = v - 1
    6. If v < 1000, the noise source outputs will be mapped down based on
     the ranking provided, and the data will be retested.*/


  size_t v = 0;  // initializing the variable v mu_bar = mu -(2.576 * sigma / math.sqrt(v))
  size_t index = 0;
  size_t tv[len]; //
  size_t i, j, valid; // variable for loops
  size_t count = 0;//  variable for storing intermediate count value
  float mu = 0, mu_bar = 0, sum = 0; // The mean value of the tv
  float sigma = 0;
  float p = 0; // local copy of output p in this function
  float temp = 0;
  float entropy = 0;

  for (index = 0; index < len; index++)
  {

    for (j = index + 1; j < len; j++)
    {
      if (input[index] == input[j])  // searching the same value in the remaining set
      {
        tv[v] = j - index + 1;
        v = v + 1;
        index = j ;
        break;
      }
    }
  }

  // The step-5 is bypassed for testing.
  //v = v -1; // Step : 5
  /*if (v < 1000)   // Step:6
    {
    p = 0.0;
    min_entropy = 0.0;
    exit(1);
    }*/
  //Calculating mean mu
  for (i = 0; i < v; i++)   // Step :7
  {
    sum = sum + tv[i];
  }
  mu = sum / v;
  sum = 0;
  //Calculating variance Sigma
  for (i = 0; i < v; i++)
  {
    sum = sum + ((tv[i] - mu) * (tv[i] - mu));
  }
  sum = sum / v;
  sigma = sqrt(sum);
  sum = 0;
  //8. Compute the lower - bound of the confidence interval for the mean
  // based on a normal distribution with confidence level alpha = 0.95

  mu_bar = mu - (2.576 * sigma / sqrt(v));
  
  // Serial.println(mu_bar);
  valid = solveForP(mu_bar, &p);
  //Test vecor
  //valid = 1;
  //p = 0.7063;
  if (valid == 0)
  { // No solution to equation.Assume max min - entropy.

    p = 1.0 / (float) len;
    entropy = log2((float)len);
  }

  else
  {
    entropy = -log2(p);
  }

  return entropy;
}

float Entropy::markovEstimate(char* input, size_t len) {

  char data[len];
  memcpy(data, input, len * sizeof(char));
  Serial.println("step0");
  int O_i, O_j; // Intrim variable
  int d = 128;
  int k = 3; // Number of symbols in input string. Need to be set before calling
  int index;
  int Oi[k];
  float L;
  float Pi[k];     // P - probalility
  float Ppi[k];    // P prime
  float epsilon_i[k];
  float Oij[k][k];   // 2-D array for Oij
  float Tij[k][k];   // 2-D array for Tij
  float hc[k];

  float temp = 0.0;
  float alpha = 0.0;
  float epsilon = 0.0;
  float Pmax = 0.0;
  float entropy = 0.0;
  /* 1. The confidence level to be alpha = min(alpha^k ^ 2, alpha^d)
        where k ^ 2 is the number of terms in the transition matrix and d = 128 is
      the assumed length of the Markov chain.*/
  L = (float)len;
  temp = (float)k * k;
  alpha = MIN(pow(0.99, temp), pow(0.99, (float)d));

  // Probe-1
  //printf("\n Value of alpha = %f \n", alpha);
  Serial.println("step1");
  /*2. Estimate the initial state probability distribution, P, with:
       Pi = min(1, o_i / L + epsilon)cv */
  for (int i = 0; i < k; i++)
  {
    Oi[i] = 0;
    Ppi[i] = 0.0;
    hc[i] = 0.0;
  }

  // Obtaining count of occurance of symbol invariable Oi
  for (int i = 0; i < k; i++)
  {
    for (int j = 0; j < len; j++)
    {
      //printf("\n value of i = %d, value of ip[%d]=%d\n",i,j,(int)(data[j]));
      if (i == (int)(data[j] - '0'))
      {
        Oi[i] = Oi[i] + 1;
        //printf("Oi[%d]=%d\n",i,Oi[i]);
      }
    }
    //printf("%d\n", Oi[i]);
  }
  // printf("\n");
  //Probe-2
  Serial.println("step2");
  float epsilon_term = log2(1 / (1 - alpha));

  epsilon = sqrt(epsilon_term / (2 * (float)len));

  //printf("\n epsilon = %f\n", epsilon);

  for (int i = 0; i < k; i++)
  {
    Pi[i] = MIN(1.0, (Oi[i] / L + epsilon));

    //printf("%d,\t %f\n", Oi[i], Pi[i]);
  }

  /*3.Let o_s_L = o_s_L - 1
      Need to subtract 1 from count of last symbol for
    bounding matrix construction.Nothing follows the last occurance,
    therefore it should not be included in transition proportion. */
  Oi[k - 1] = Oi[k - 1] - 1;
  //printf("\n%d,\t %d", Oi[k-1], Oi[2]);

  Serial.println("step3");
  /*4. Estimate the probabilities in the transition matrix T, overestimating where...
     Ti, j = 1 if o_i = 0
     min(1, oi, j + eps_i) otherwise */
  //Initializing the Matrix
  for (int i = 0; i < k; i++)
  {
    for (int j = 0; j < k; j++)
    {
      Oij[i][j] = 0.0;
      //Tij[i][j] = 1.0;
    }
  }

  O_i = (int)(data[0] - '0');
  for (int j = 1; j < len; j++)
  {
    //O_i = (int)(data[i]-'0');
    //for (int j = 0; j < len; j++)
    //{
    O_j = (int)(data[j] - '0');
    Oij[O_i][O_j] ++;
    O_i = O_j;
    //Tij[i][j] = 1.0;
    //}

  }

  for (int i = 0; i < 3; i++)
  {
    epsilon_i[i] = sqrt(epsilon_term / (2 * (float) Oi[i]));
    //printf("\n%f\t %f\t%d\n", epsilon_i[i], epsilon_term,Oi[i]);
  }

  for (int i = 0; i < k; i++)
  {
    for (int j = 0; j < k; j++)
    {
      if (Oi[i] == 0)
        Tij[i][j] = 1;
      else
        //temp = Oij[i][j] / Oi[i] + epsilon_i[i];
        Tij[i][j] = MIN(1.0, (Oij[i][j] / Oi[i] + epsilon_i[i]));
      // printf("\t %f", Tij[i][j]);

    }
    // printf("\n");
  }
  Serial.println("step4");
  for (int j = 0; j < d - 1; j++)
  {
    for (int c = 0; c < k; c++)
    {
      float h_temp = 0.0;
      for (int i = 0; i < k; i++)
      {
        Ppi[i] = Pi[i] * Tij[i][c];
        h_temp = MAX(Ppi[i], h_temp);
      }
      hc[c] = h_temp;
      //Pi[c] = hc[c];

    }

    for (int i = 0; i < k; i++)
    {
      Pi[i] = hc[i];
    }
  }

  for (int i = 0; i < k - 1; i++)
  {
    Pmax = MAX(Pi[i], Pi[i + 1]);

  }

  Serial.println("step5");
  //printf("\n Pmax =%f",Pmax);

  // The value of Pmax is reduced to 0.0000 after 128 iteration. therfore minimum value is assigned
  //Pmax=MAX(0.001,Pmax);
  Pmax = 0.00001;
  Serial.println("Pmax");
  //Serial.println(Pmax);
  entropy = -log2(Pmax) / d;

  return entropy;
}

// 6.3.4
float Entropy::compressionEstimate(char* input, size_t len) {

  //char* data;
  //data = (char*)malloc(sizeof(char) * len); // Visual Studio
  char data[len];
  memcpy(data, input, len * sizeof(char));

  size_t k = 3;
  size_t d = 10;
  size_t v = 0;
  size_t D[15]; // variable for storing the D
  float b = 0.0;
  size_t max_sk = 0;
  float mean_x = 0.0;
  float prime_mean_x = 0.0;
  float c = 0.0;
  float std_dev = 0.0;

  float p = 0.0;
  float min_entropy = 0.0;

  int valid;


  // Calculating v
  v = len - d;

  // 2. Initializing the dictionary
  int dict[3];
  for (size_t i = 0; i < k; i++)
  {
    dict[i] = 0;
  }

  for (size_t i = 0; i < d; i++)
  {
    for (size_t j = 0; j < k; j++)
    {
      if (j == (int)(data[i] - '0'))
      {
        dict[j] = i + 1;
      }
    }

  }

  /*for (size_t j = 0; j < k; j++)
    {

    printf("value of dictionary obained are dict[%d]=%d \n",j,  dict[j]);

    }*/
  for (size_t i = d; i < len; i++)
  {
    for (size_t j = 0; j < k; j++)
    {
      if ((dict[j] != 0) && (j == (int)(data[i] - '0')))
      {
        D[i - d] = i - dict[j] + 1;
        dict[j] = i + 1;
      }

      if ((dict[j] == 0) && (j == (int)(data[i] - '0')))
      {
        D[i - d] = i + 1;
        dict[j] = i + 1;
      }

    }
    //printf("Value of D[%d] = %d\n", (i - d), D[i - d]);

  }
  // Step 4 -
  for (size_t i = 0; i < len - 1; i++)
  {
    max_sk = max(((int)(data[i] - '0')), max_sk);
  }

  //printf("\n max sk = %d", max_sk);
  b = log2((float)max_sk) + 1.0;  // The Max(x1,x2,x3..xk)
  //printf("\n b = %f", b);
  for (size_t i = 0; i < v; i++)
  {
    mean_x += log2(D[i]);
  }
  mean_x = mean_x / (float)v;
  //printf("\n Mean = %f", mean_x);

  float term = (-3.0) / b;
  c = 0.7 - (0.8 / b) + (4.0 + (32.0 / b)) * pow((float)v, term) / (float)v;
  //printf("\n c = %f", c);

  term = 0.0;
  float term1 = 0.0;
  for (size_t i = 0; i < v; i++)
  {
    term += (log2(D[i]) * log2(D[i]));
  }
  term1 = (term / v) - (mean_x * mean_x);

  std_dev = c * sqrt(term1);

  //Step -5 Calculating the mean_x_prime

  prime_mean_x = mean_x - (2.57 * std_dev) / sqrt((float)v);
  //printf("\n X prime = %f", prime_mean_x);
  valid = solve_for_p_634(prime_mean_x, k, v, d, &p);

  if (valid == 0)
  { // No solution to equation.Assume max min - entropy.

    p = 1.0 / (float)len;
    min_entropy = log2((float)len);
    //printf("\n No valid solution p= %f", p);

  }

  else
  {
    min_entropy = -log2(p);
    //printf("\n valid solution p= %f", p);

  }

  return (min_entropy);

}

//6.3.5
float Entropy::tTupleEstimate(char* input, size_t len) {

  //char* data; // For VS2015
  //data = (char*)malloc(sizeof(char) * len); // For VS2015
  char data[len];
  memcpy(data, input, len * sizeof(char));

  float min_entropy = 0.0;
  size_t L = 21;
  float Q[21];      // should be optimized for L/2
  float P[21];    // should be optimized for L/2
  float P_max[21];  // should be optimized for L/2
  float Pmax = 0.0;
  float temp = 0.0;
  int tuple_found = 0;
  // counter variables
  size_t tuple_size = 0;
  size_t step = 0;
  size_t i = 0;
  size_t j = 0;

  // Initializing the Q[i]
  for (i = 0; i < L; i++)
  {
    Q[i] = 0.0;
    P[i] = 0.0;
    P_max[i] = 0.0;
  }
  // step -1 is bypassed due to small data sample

  // Step -2 for finding Q[i]
  // first loop for forming tuples form 1 to L/2
  for (tuple_size = 1; tuple_size <= (size_t)(L / 2); tuple_size++)
  {
    // Loop for selecting tuple of given size
    for (step = 0; step < L; step = step + tuple_size)
    {
      temp = 1.0; // The minimum occurance of any tuple is 1
      // loop for searching similar tuple
      //for (i = step + tuple_size; i < L - tuple_size; i++)
      for (i = step + tuple_size; i <= L - tuple_size; i++)
      {
        tuple_found = 1;
        // loop for matching the element of tuple
        for (j = 0; j < tuple_size; j++)
        {
          //if (data[i-tuple_size+j]!=data[i+j])
          if (data[step + j] != data[i + j])
          {
            tuple_found = 0;
            break;
          }

        } // end of tuple element matching
        if (tuple_found)
        {
          temp++;
          //printf("\n Tuple match found"); //Probe
        }

      } // end of tuple count search
      Q[tuple_size] = max(Q[tuple_size], temp);
      //printf("Value of temp = %f",temp);  // Probe
      //printf("\n Q[%d]= %f",tuple_size,Q[tuple_size]); // Probe
    }
  }

  // Step -3
  for (i = 1; i <= L / 2; i++)
  {
    P[i] = Q[i] / (float)(L - i + 1);
    P_max[i] = pow(P[i], (1.0 / (float)i));
    Pmax = max(P_max[i], Pmax);
  }

  //printf("\n Pmax =%f", Pmax);
  min_entropy = -log2(Pmax);
  return (min_entropy);
  //return 0;
}
//6.3.6 lrsEstimate
float Entropy::lrsEstimate(char* input, size_t len) {
  //char* data;
  //data = (char*)malloc(sizeof(char) * len);
  char data[len];
  memcpy(data, input, len * sizeof(char));

  float min_entropy = 0.0;
  size_t L = len;
  float Q[21];      // should be optimized for L/2
  //float P[21];    // should be optimized for L/2
  //float P_max[21];  // should be optimized for L/2
  //float Pmax = 0.0;
  float temp = 0.0;

  float P_w[21];
  float Pmax_w[21];
  float Pmaxw = 0.0;

  int tuple_found = 0;
  // counter variables
  size_t tuple_size = 0;
  size_t step = 0;
  size_t i = 0;
  size_t j = 0;
  size_t u = 0;
  size_t v = 0;

  // Initializing the Q[i]
  for (i = 0; i < L; i++)
  {
    Q[i] = 0.0;
    P_w[i] = 0.0;
    Pmax_w[i] = 0.0;

  }
  // step -1 is by passed due to small data sample

  // Step -2 for finding Q[i]
  // first loop for forming tuples form 1 to L/2
  for (tuple_size = 1; tuple_size <= (size_t)(L / 2); tuple_size++)
  {
    // Loop for selecting tuple of given size
    for (step = 0; step < L; step = step + tuple_size)
    {
      temp = 1.0; // The minimum occurance of any tuple is 1
      // loop for searching similar tuple
      //for (i = step + tuple_size; i < L - tuple_size; i++)
      for (i = step + tuple_size; i <= L - tuple_size; i++)
      {
        tuple_found = 1;
        // loop for matching the element of tuple
        for (j = 0; j < tuple_size; j++)
        {
          //if (data[i-tuple_size+j]!=data[i+j])
          if (data[step + j] != data[i + j])
          {
            tuple_found = 0;
            break;
          }

        } // end of tuple element matching
        if (tuple_found)
        {
          temp++;
          //printf("\n Tuple match found"); //Probe
        }

      } // end of tuple count search
      Q[tuple_size] = max(Q[tuple_size], temp);
      //printf("Value of temp = %f",temp);  // Probe
      //printf("\n Q[%d]= %f",tuple_size,Q[tuple_size]); // Probe
    }
  }

  //Step-1/2 Finding u and v

  for (i = 1; i <= L / 2; i++)
  {
    //printf("\n Q[%d]=%f", i, Q[i]);
    if ((Q[i] < 10.0) && (Q[i] > 2.0))
    {
      u = i;

    }

    if ((Q[i] >= 2.0) && (Q[i + 1] == 1.0))
    {
      v = i;
    }

  }

  //printf("\n Value of u=%d and v=%d", u, v);
  // Step -3
  for (i = u; i <= v; i++)
  {
    P_w[i] = (Q[i] * (Q[i] - 1.0) / 2.0) / ((float)(L - i + 1) * ((float)(L - i)) / 2.0);
    Pmax_w[i] = pow(P_w[i], (1.0 / (float)i));
    Pmaxw = max(Pmaxw, Pmax_w[i]);

  }

  // Step -4

  //printf("\n Pmax =%f", Pmaxw);
  min_entropy = -log2(Pmaxw);
  return (min_entropy);
  //return 0;
}
//3.6.7 multiMCW
float Entropy::multiMCWEstimate(char* input, size_t len) {
  //char* data; // VS 2015
  // data = (char*)malloc(sizeof(char) * len); // VS 2015
  char data[len];
  memcpy(data, input, len * sizeof(char));
  
  float min_entropy=0.0;
  size_t L = len;
  size_t wlen;   // Variable for window lenght used in finding most common value
  char prediction; // Prediction variable
  size_t C = 0; // Variable for storing sum of correct[]
  //size_t s[3] = { 0,1,2 }; // array of input symbols
  size_t i, j; //Loop counters
  float Pavg = 0.0;  // Prime of Pglobal 
  float Pglob = 0.0; // is Pglobal
  float Prun = 0.0; // is Plocal
  float Pmax = 0.0; // The Claculated Plocal
  size_t r = 0; // Variable for storing the maximun run of correct
  // Step -1
  // Window size w1=63, w2=255,w3=1023,w4=4095
  // For our experiment w1=3,w2=5,w3=7,w4=9
  // size_t w[4] = {63,255,1023,4095};

  size_t w[4] = {3,5,7,9 };
  size_t N;
  N = L - w[0];
  size_t correct[12];  // VLA is not permitted in VS2015
  char mdata[12]; // Buffer to store temporary string for finding the mostCommon
  //size_t *correct  = (size_t*)alloca(sizeof(size_t) * N);
  // initializing the correct arra
  for (i = 0; i < N; i++)
  {
    correct[i] = 0;
  }

  //Step -2
  size_t scoreboard[4] = {0,0,0,0};
  char frequent[4] = { 0 };
  size_t winner = 0;
  
  size_t counters[4] = { 0,0,0,0 }; // Additonal variable

  //Step -3

  for (i = w[0] ; i < L; i++)
  {
    //printf("\n i = %d\n",i+1);
    for (j = 0; j < 4; j++)
    {
      if (i > w[j]-1)
      {
        size_t count = 0;
        for (size_t k = i - w[j]; k < i ; k++)
        {
          mdata[count] = data[k];
          //printf("%c", mdata[count]);
          count++;
        }
        
        wlen = w[j];
        frequent[j] = mostComman(mdata, wlen);
      }
      else
      {
        frequent[j] = 0;
      }
      //printf(" frequent [%d]= %c", j, frequent[j]);

    }
    //printf("\n");


    //Step- 3b

    prediction = frequent[winner];

    //Step- 3c
    if (prediction == data[i] )
    {
      correct[i - w[0]] = 1;
    }

    //Step- 3d


    //printf("\n ");
    for (j = 0; j < 4; j++)
    {
      if (frequent[j] == data[i])
      {
        scoreboard[j] += 1;
        if (scoreboard[j] >= scoreboard[winner])
        {
          winner = j;
        }
      }
      //printf("scoreboard[%d] = %d", j, scoreboard[j]);
    }//end of loop for j

    //printf("\n Winner = %d, Prediction = %c, Correct[%d] = %d", winner+1, prediction, i,correct[i-w[0]]);

        
  } // end of i loop

  // Step - 4

  for (i = 0; i < N; i++)
  {
    C = C + correct[i];
  }
    
  //Step - 5
  Pglob = (float)C / (float)N;
  Pavg = Pglob+2.576*sqrt((Pglob*(1-Pglob))/(float)(N-1));

  //printf("\nPavg = %f", Pavg);


  //step - 6
  //Prun = calcRun(correct);
  r = findMaxRun(correct, N);
  r =r+1;
  //printf("Value of r = %d",r);
  Prun = calcRun(r, N);

  //step - 7
  Pmax = max(Pavg, Prun);
  //printf("\n Pmax =%f",Pmax);
  //printf("\n Pmax =%f", Prun);
  //Pmax = 0.76;
  min_entropy = -log2(Pmax);
  return (min_entropy);
  //return 0;


}
//3.6.8 lagPredictionEstimate
float Entropy::lagPredictionEstimate(char* input, size_t len)
{
  //char* data; // for VS 2015
  //data = (char*)malloc(sizeof(char) * len); // for VS 2015
  char data[len];
  memcpy(data, input, len * sizeof(char));
  
  float min_entropy=0.0;
  size_t L = len;
  //size_t wlen;   // Variable for window lenght used in finding most common value
  char prediction; // Prediction variable
  size_t C = 0; // Variable for storing sum of correct[]
  //size_t s[3] = { 0,1,2 }; // array of input symbols
  size_t i, j, d; //Loop counters
  float Pavg = 0.0;  // Prime of Pglobal 
  float Pglob = 0.0; // is Pglobal
  float Prun = 0.0; // is Plocal
  float Pmax = 0.0; // The Claculated Plocal
  size_t r = 0; // Variable for storing the maximun run of correct
  // Step -1
  size_t N;
  size_t D;
  N = L - 1;
  D = 3;
  size_t correct[12];  // VLA is not permitted in VS2015
  //char mdata[12]; // Buffer to store temporary string for finding the mostCommon
  //size_t *correct  = (size_t*)alloca(sizeof(size_t) * N);
  // initializing the correct arra
  for (i = 0; i < N; i++)
  {
    correct[i] = 0;
  }

  //Step -2
  size_t scoreboard[3] = {0,0,0};
  char lag[3] = {0};
  size_t winner = 0;
  
  //size_t counters[4] = { 0,0,0,0 }; // Additonal variable

  //Step -2

  for (i = 1 ; i < L; i++)
  {
    //printf("\n i = %d\n",i+1);
    //Step -2a
    for (d = 0; d < D; d++)
    {
      if (d <i)
      {
        lag[d] = data[i-d-1];
      }
      else
      {
        lag[d] = 0;
      }
      //printf(" lag [%d]= %c", d, lag[d]);

    }
    //printf("\n");


    //Step- 2b

    prediction = lag[winner];

    //Step- 2c
    if (prediction == data[i] )
    {
      correct[i - 1] = 1;
    }

    //Step- 2d
    
    //printf("\n ");
    for (d = 0; d < D; d++)
    {
      if (lag[d] == data[i])
      {
        scoreboard[d] += 1;
        if (scoreboard[d] >= scoreboard[winner])
        {
          winner = d;
        }
      }
      //printf("scoreboard[%d] = %d", d, scoreboard[d]);
    }//end of loop for j

    //printf("\n Winner = %d, Prediction = %c, Correct[%d] = %d", winner+1, prediction, i,correct[i-1]);

        
  } // end of i loop

  // Step - 4

  for (i = 0; i < N; i++)
  {
    C = C + correct[i];
  }
    
  //Step - 5
  Pglob = (float)C / (float)N;
  Pavg = Pglob+2.576*sqrt((Pglob*(1.0-Pglob))/(float)(N-1));

  //printf("\nPavg = %f", Pavg);


  //step - 6
  //Prun = calcRun(correct);
  r = findMaxRun(correct, N);
  r =r+1;
  //printf("Value of r = %d",r);
  Prun = calcRun(r, N);

  //step - 7
  Pmax = max(Pavg, Prun);
  //printf("\n Pmax =%f",Pmax);
  //printf("\n Prun =%f", Prun);
  //Pmax = 0.76;
  min_entropy = -log2(Pmax);
  return (min_entropy);
  //return 0;
}

//3.6.9 multiMMCEstimate
float Entropy::multiMMCEstimate(char* input, size_t len)
{
  //char* data; // for VS2015
  //data = (char*)malloc(sizeof(char) * len); // for VS2015
  char data[len];
  memcpy(data, input, len * sizeof(char));
  
  float min_entropy=0.0;
  size_t L = len;
  //size_t N = L - 2;
  //size_t D = 3;
  char prediction; // Prediction variable
  size_t C = 0; // Variable for storing sum of correct[]
  //size_t s[3] = { 0,1,2 }; // array of input symbols
  size_t i, j,k,l,m,n, d; //Loop counters
  size_t matched = 0;
  float Pavg = 0.0;  // Prime of Pglobal 
  float Pglob = 0.0; // is Pglobal
  float Prun = 0.0; // is Plocal
  float Pmax = 0.0; // The Claculated Plocal
  size_t r = 0; // Variable for storing the maximun run of correct
  // Step -1
  size_t N;
  size_t D;
  N = L - 2;
  D = 3;
  size_t correct[12];  // VLA is not permitted in VS2015
  char subpredict[3] = {0};
  char y = 0;
  char ymax = 0;
  char prdiction;
  
  // initializing the correct array
  for (i = 0; i < N; i++)
  {
    correct[i] = 0;
  }
  
    
  //Step -2
  char MMC[10] = { 0 };
  char Md[3][4] = { 0 };

  // step -3
  size_t scoreboard[3] = { 0,0,0 };
  size_t winner = 0;
  // Step -4
  for (i = D - 1; i < L; i++)
  {
    //printf("\n i = %d\n", i + 1);
    //Step -4a
    for (d = 0; d < D; d++)
    {
      if (d < i-1)
      {
        for (k = 0; k < i ; k++)
        {
          MMC[k] = data[k];
          //printf("%c", MMC[k]);
        }
        //printf("\n");
      }


    } // end of d loop
    //printf("\n");
    //Step - 4b

    for (d = 0; d < D; d++)
    {
      l = 0;
      if (d < i - 1)
      {
        for (k = i - d - 1; k < i; k++)
        {
          Md[d][l] = data[k];
          //printf("Md[%d][%d] = %c", d, l, Md[d][l]);
          //printf("\n k ,i = %d,%d", k, i);
          l++;
        }
        //printf("\n");
      }
      
      
    } // end of d loop

    // Loop for searching ymax
    for (d = 0; d < D; d++)
    {
      //size_t l = 0;
      y = 0;
      ymax = 0;

      for (k = 0 ; k < i-d; k++)
      {
        matched = 0;
        y = 0;
        for (n = 0; n <= d; n++)
        {
          if (Md[d][n] == MMC[k + n])
          {
            matched ++;
            //break;
          }
        }
        if (matched == d+1)
        {
          y = MMC[k + n] ;
          if (y > ymax)
            ymax = y;

        }
        
      }

      subpredict[d] = ymax;
      //printf("\n");
      //printf("subpredict [%d] = %c", d, subpredict[d]);

    } // end of d loop
 
    //Step- 4c

    prediction = subpredict[winner];

    //Step- 4d
    if (prediction == data[i] )
    {
      correct[i - 2] = 1;
    }

    //Step- 4e
    
    //printf("\n ");
    for (d = 0; d < D; d++)
    {
      if (subpredict[d] == data[i])
      {
        scoreboard[d] += 1;
        if (scoreboard[d] >= scoreboard[winner])
        {
          winner = d;
        }
      }
      //printf("scoreboard[%d] = %d", d, scoreboard[d]);
    }//end of loop for d

    //printf("\n Winner = %d, Prediction = %c, Correct[%d] = %d", winner+1, prediction, i,correct[i-2]);

        
  } // end of i loop (master loop)

  // Step - 5

  for (i = 0; i < N; i++)
  {
    C = C + correct[i];
  }
    
  //Step - 6
  Pglob = (float)C / (float)N;
  Pavg = Pglob+2.576*sqrt((Pglob*(1.0-Pglob))/(float)(N-1));

  //printf("\nPavg = %f", Pavg);


  //step - 7
  //Prun = calcRun(correct);
  r = findMaxRun(correct, N);
  r =r+1;
  //printf("Value of r = %d",r);
  Prun = calcRun(r, N);

  //step - 8
  Pmax = max(Pavg, Prun);
  //printf("\n Pmax =%f",Pmax);
  //printf("\n Prun =%f", Prun);
  min_entropy = -log2(Pmax);
  return (min_entropy);
  //return 0;
}


//3.6.10 LZ78YEstimate
float Entropy::LZ78YEstimate(char* input, size_t len)
{
  //char* data;
  //data = (char*)malloc(sizeof(char) * len);
  char data[len];
  memcpy(data, input, len * sizeof(char));
  
  float min_entropy=0.0;
  size_t L = len;
   
  char prediction = 0; // Prediction variable
  size_t C = 0; // Variable for storing sum of correct[]
  int i, j,k,l,m,n; //Loop counters
  size_t matched = 0;
  float Pavg = 0.0;  // Prime of Pglobal 
  float Pglob = 0.0; // is Pglobal
  float Prun = 0.0; // is Plocal
  float Pmax = 0.0; // The Claculated Plocal
  size_t r = 0; // Variable for storing the maximun run of correct
  
  // Step -1
  size_t N;
  size_t B;
  B = 4;
  N = L - B - 1;
  size_t correct[len];  // VLA is not permitted in VS2015
  
  char y = 0;
  char maxCount = 0;
  
  
  // initializing the correct array
  for (i = 0; i < len; i++)
  {
    correct[i] = 0;
  }
  
  size_t maxDictionarySize = 65536;
  
  //Step -2
  size_t dictSize[4] = { 0 }; // Vairable to keep track of dictionary sizes
  char D[4][100] = { 0 };        // Initialzing all the 4 dictionaries with null
  char prev[5] = { 0 };

   
  for (i = B + 1; i < L; i++)
  {
    //printf("\n i = %d\n", i + 1); // for VS2015
    //Step -4a
    for (j = B - 1; j >= 0; j--)
    {
      for (k = 0; k < j + 2; k++)
      {
        D[j][dictSize[j] + k + 1] = data[i - j - 2 + k];
       // printf("%c", D[j][dictSize[j] + k + 1]); // for VS2015
      }
      dictSize[j] = dictSize[j] + j + 2;
      //printf("\n"); // for VS2015
      //printf("size of dictionary %d\n", dictSize[j]); // for VS2015

    } // end of j loop
    //printf("\n"); // for VS2015
        
    // Loop for searching maxCount
    for (j = B-1; j>=0; j--)
    {
      //size_t l = 0;
      y = 0;
      maxCount = 0;
      //printf("\nPrev ="); // for VS2015
      // 3b(i)creating prev array
      for (k = 0; k < j + 1; k++)
      {
        prev[k] = data[i - (j + 1) + k];
        //printf("%c",prev[k]); // for VS2015
      }
      //Step - 3b(ii)
      for (k = 0; k < dictSize[j] - (j+2); k=k+j+2)
      {
       // printf("\n k =%d", k); // for VS2015
        matched = 0;
        y = 0;
        for (n = 0; n < j+1; n++)
        {
          if (D[j][k+n+1] == prev[n])
          {
            matched++;
            
          }
        }
        if (matched == (j + 1))
        {
          y = D[j][k + n+1];
          //printf("\n matched\n"); // for VS2015
        }
         //Step - 3b(iii)
        if (y > maxCount)
        {
          prediction = y;
          maxCount = y;
        }
          
        //printf("\nMaxD = %c",maxCount); // for VS2015
      }   

    } // end of j loop

      
    //Step - 4
    if (prediction == data[i] )
    {
      correct[i - B - 1] = 1;
    }
    
   //printf("\n Prediction = %c, Correct[%d] = %d,Si = %c", prediction, i-B-1,correct[i-B-1],data[i]); // for VS2015
        
  } // end of i loop (master loop)

  // Step - 5

  for (i = 0; i < N; i++)
  {
    C = C + correct[i];
    //printf("%d", correct[i]); //for VS 2015
  }
    
  //Step - 6
  Pglob = (float)C / (float)N;
  Pavg = Pglob+2.576*sqrt((Pglob*(1.0-Pglob))/(float)(N-1));

  //printf("\nPavg = %f", Pavg); //for VS 2015


  //step - 7
  r = findMaxRun(correct, N);
  r =r+1;
  //printf("Value of r = %d",r); //for VS 2015
  Prun = calcRun(r, N);

  //step - 8
  Pmax = max(Pavg, Prun);
  //printf("\n Pmax =%f",Pmax); //for VS 2015
  //printf("\n Prun =%f", Prun); //for VS 2015
  min_entropy = -log2(Pmax);
  return (min_entropy);
  //return 0;
}



float Entropy::log2( float n ) {
  // log(n)/log(2) is log2.
  return log( n ) / log( 2 );
}


void Entropy::sortArray(char* input, size_t len) {
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

size_t Entropy::solveForP(float mu_bar, float* p)
{
  float minp = 0;
  float p1 = 0;
  float adj = 0;
  float Ep = 0;
  float Ep_maxvalid = 0;
  size_t k = 10;
  //float pSol = 0;
  minp = 1.0 / (float)k;
  p1 = (1 - minp) / 2.0 + minp;
  adj = (1 - minp);
  Ep = calcEpS(p1);
  Ep_maxvalid = calcEpS(minp); //The function is modified for global use of k
  if (mu_bar > Ep_maxvalid)
  {
    *p = (float)0.0;
    Serial.println(*p);
    return 0;
  }
  else
  {
    while (fabs(mu_bar - Ep) > 0.0001)
    {
      adj /= 2.0;
      if (mu_bar < Ep)
      {
        p1 += adj;
      }

      //caclEps will crash if p_c is exactly 1.
      if (p1 == 1.0)
      {
        p1 -= 0.0001;
      }
      else
      {
        p1 -= adj;
      }

      //#occasionally dips below lowest possible pmax.This is to fix that
      if (p1 < minp)
      {
        p1 = minp;

      }

      Ep = calcEpS(p1); // Func is modified for global read of k

    }
    *p = p1;
    Serial.println(*p);
    return 1;
  }

}
float Entropy::calcEpS(float p1)
{
  //printf("\n p_c calcEpS %f", p_c);
  float q = 0;
  float p2 = 0;
  float i_k = 0;
  float  ip = 0;
  float iq = 0;
  float iq2 = 0;
  float Ep = 0;
  float denom = 0;
  float q_inv = 0;
  float F_q = 0;
  int i = 0;
  size_t k = 10;
  q = (1.0 - p1) / (float)(k - 1);
  // implementing F calulation inline
  q_inv = 1.0 / q;
  denom = 1.0 + (float)k / q_inv;
  for (i = 1; i <= k; i++)
  {
    denom = q_inv + (float)(-k) / denom;
    denom = 1.0 + (float)(k - i) / denom;
  }
  denom = q_inv + (float)(-k) / denom;
  F_q = (1.0) / denom;

  i_k = 1.0 / (float)k;
  ip = 1.0 / p1;
  iq = 1.0 / q;
  iq2 = iq * iq;
  Ep = p1 * iq2 * (1.0 + i_k * (ip - iq)) * F_q - (p1 * iq * i_k * (ip - iq));
  return (Ep);
}
// functions of 6.3.4
int Entropy::solve_for_p_634(float mu_bar, size_t n, size_t v, size_t d, float *p)
{
  //*p = 0.1;
  //*min_entropy = 0.5;

  float tollerance = 0.0001;
  float min_p = 1 / (float)n;
  float p_1 = 0.0;
  float Ep = 0.0;
  p_1 = (1.0 - min_p) / (2.0 + min_p);
  float adj = 1.0 - min_p;
  //float Ep_maxvalid =0.0;
  float Ep_maxvalid = EppM(min_p, n, v, d);
  //printf("\n Value of mu_bar is %f , Ep_max =%f", mu_bar, Ep_maxvalid);
  if (mu_bar > Ep_maxvalid)
  {
    //printf("\n Value of mu_bar is %f , Ep_max =%f",mu_bar, Ep_maxvalid);
    *p = 0.0;
    return 0;
  }
  Ep = EppM(p_1, n, v, d);

  //printf("\n absolte of mu_bar -Ep = %f", fabs(mu_bar-Ep));
  while (fabs(mu_bar - Ep) > tollerance)
  {
    //printf("\nEntered in to while loop");
    adj = adj / 2.0;
    if (mu_bar > Ep)
    {
      p_1 = p_1 - adj;
    }
    else
    {
      p_1 = p_1 + adj;
    }
    Ep = EppM(p_1, n, v, d);
    //printf("\nEp = %f, p = %f", Ep, p_1);
  }

  *p = p_1;
  return (1);
}

float Entropy::EppM(float p, size_t n, size_t v, size_t d)
{
  //printf("\nvalue of p, n and v received is %f, %d, %d", p, n, v);
  float q = 0.0;
  float Ep = 0.0;
  q = (1.0 - p) / (n - 1.0);

  Ep = func_G(p, v, d) + (float)(n - 1) * func_G(q, v, d);
  //printf("\n Ep in EppM =%f", Ep);
  return (Ep);
}

float Entropy::func_G(float p, size_t v, size_t d)
{
  float Gsum = 0.0;
  float N = (float) v + (float)d;
  float inSum = 0.0;
  float tempSum = 0.0;
  float st[100];
  float temp = 0.0;

  // The inner sum for s = 1->d...compute once...use in every iteration of outer loop
  //printf("\n value of p/q, v  and N is %f, %d, %f", p, v, N);
  for (int i = 1; i <= d; i++)
  {
    tempSum = tempSum + (log2((float)i) * pow((1.0 - p), ((float)i - 1.0)));
  }
  inSum = p * p * tempSum;
  //printf("\n inSum in funcG =%f", inSum);
  //Serial.println(tempSum);
  //inSum = p * p * sum([log(s, 2.0) * pow(1.0 - p, s - 1) for s in range(1, (d + 1))])
  Gsum = inSum * (N - (float)d);
  //printf("\n Gsum in funcG =%f", Gsum);

  //Serial.println(Gsum);
  for (int i = d + 1; i <= (int)N ; i++)
  {
    st[i - (d + 1)] = log2((float)i) * pow((1.0 - p), ((float)i - 1.0));

    //printf("\n ST[%d] = %f",(i-(d+1)), st[i - (d + 1)]);
  }
  //Serial.println(st[10]);
  //# The additional 's < t' sum term for s = (d + 1)->N - 1

  tempSum = 0.0;
  for (int i = 0; i <= (int)N - (d + 1); i++)
  {
    tempSum = tempSum + ((N - (float)i - (float)(d + 1.0)) * (st[i])) ;
    //printf("\nTempsum = %f, N = %f, i= %d, ST = %f", tempSum, N, i, st[i]);
  }
  //Serial.println(tempSum);
  Gsum += p * p * tempSum;
  //printf("\n Gsum after ST addition in EppM =%f", Gsum);
  //# The 's = t' term for s = (d + 1)->N

  tempSum = 0.0;

  for (int i = 0; i <= (int)N - (d + 1); i++)
  {
    tempSum = tempSum + st[i];
  }

  Gsum += p * tempSum;
  //Serial.println(tempSum);
  //printf("\n Final Gsum =%f", Gsum);
  temp = Gsum / (float)v;
  //Serial.println(temp);
  return (temp);
  //return(0.0);

}

// Functions for 6.3.7 to 6.3.10
//Function for finding the most common value in the array
char Entropy:: mostComman(char* data, size_t length)
{
  size_t maxCount = 0;  // Initializing the maxCount Variable
  size_t count = 0;
  char commonValue;
  for (size_t i = 0; i < length; i++)
  {
    count = 0;   // re-initalizing the count for next common value
    for (size_t j = i; j < length; j++)
    {
      if (data[i] == data[j])  // searching the same value in the remaining set
      {
        count++;
      }
    }// End of for loop for j

     //printf("\n Value = %c, Count = %d",input[i],count); // Probe for displaying the result

     // Comparing the count value with previous maxCount value, if heigher replace it.

    if (count >= maxCount)
    {
      commonValue = data[i];
      maxCount = count;
      //printf("\n most Common value %c\n", commonValue);
    }


  } // End of for loop for i
  
  return(commonValue);
}
// function to find maximum run in correct
size_t Entropy::findMaxRun(size_t *correct, size_t N)
{
  size_t run = 0;
  size_t maxrun = 0;
  for (size_t i = 0; i < N; i++)
  {
    if (correct[i] == 0)
    {
      if (run > maxrun)
      {
        maxrun = run;
      }
      run = 0;
    }
    if (correct[i] == 1)
    {
      run += 1;
    }
    
  }
  if (run > maxrun)
  {
    maxrun = run;
  }
  return(maxrun);
}

// fucntion to calculate calcRun

float Entropy::calcRun(size_t r, size_t N)
{

  float alpha = 0.99;
  float p = 0.8;
  float adj = 0.5;
  float qn;
  qn = calc_qn(p, r, N);
  for (size_t i = 0; i < 30; i++)
  {
    adj = adj / 2.0;
    if (qn > alpha)
    {
      p = p + adj;
    }
    else
    {
      p = p - adj;
    }
    qn = calc_qn(p, r, N);
    if (fabs(qn - alpha) < 0.01)
    {
      break;
    }
  }
  //printf("\n Value of P= %f", p);
  return(p);
}
//Function to caluculate calc_qn
float Entropy::calc_qn(float p, size_t r, size_t n)
{
  float x;
  float q;
  float qn = 0.0;
  q = 1.0 - p;
  x = find_root(p, r);
  qn = (1.0 - p*x) / (((float)r + 1.0 - (float)r*x)*q);
  qn = qn / (pow(x,((float)n+1.0)));
  //printf("\n value of qn = %f",qn);
  return(qn);
}

// Funciton to calculate find_root

float Entropy::find_root(float p, size_t r)
{
  float q = 0.0;
  float s = 1.0;
  q = 1 - p;
  for (size_t i = 0; i < 10; i++)
  {
    s = 1 + q*pow(p, (float)r)*pow(s, ((float)r + 1.0));
  }
  //printf("\n Value of s =%f", s);
  return(s);
}

