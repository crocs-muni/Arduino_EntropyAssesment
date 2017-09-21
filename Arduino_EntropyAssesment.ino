#include "Entropy.h"

Entropy e;

char string[501];
unsigned long time1;
unsigned long time2;

void setup() {
  Serial.begin(9600);
  
}

void loop() {
  // put your main code here, to run repeatedly:
  Serial.println("Program ready");
  fillString();
  Serial.println("String ready");
  time1 = millis();
  float minE = e.mostCommonValueEstimate(&string[0], 500);
  time2 = millis();
  Serial.println(minE);
  Serial.println(time2-time1);

  for(;;){
    delay(1000);
  }
}

void fillString(){
   for(int i = 0; i < 500; i++){
      string[i] = random(255);    
      if(i % 100 == 0){
         Serial.println(i);
      }
   }   
}
