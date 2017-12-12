#include "Entropy.h"

Entropy e;

char string[501];
unsigned long time1;
unsigned long time2;
size_t k=255;
void setup() {
  Serial.begin(9600);
  
}

void loop() {
  // put your main code here, to run repeatedly:
  Serial.println("Program ready");
  fillString();
  Serial.println("String ready");
  time1 = millis();
  //float minE=0;
  //6.3.1
  //float minE = e.mostCommonValueEstimate(&string[0], 500);
  //6.3.2
  //float minE = e.collisionEstimate(&string[0], 500);
  //6.3.3
  //float minE = e.markovEstimate(&string[0],100);
  //6.3.4
  //float minE = e.compressionEstimate(&string[0],100);
  //6.3.5
  //float minE = e.tTupleEstimate(&string[0],100);
  //6.3.6
  //float minE = e.lrsEstimate(&string[0],10);
  //6.3.7
  //float minE = e.multiMCWEstimate(&string[0],10);
  //6.3.8
  //float minE = e.lagPredictionEstimate(&string[0],5);
  //6.3.9
  //float minE = e.multiMMCEstimate(&string[0],12);
  //6.3.10
  float minE = e.LZ78YEstimate(&string[0],12);
  
  time2 = millis();
  Serial.println(minE);
  Serial.println(time2-time1);

  for(;;){
    delay(1000);
  }
}

void fillString(){
   for(int i = 0; i < 100; i++){
      string[i] = random(3);    
      if(i % 10 == 0){
         Serial.println(i);
      }
   }   
}
