#include <Arduino.h>
#include <iostream>
#include <cstring>

int Serial::println(char* buf){
  cout << buf << endl;
  return strlen(buf);
}
