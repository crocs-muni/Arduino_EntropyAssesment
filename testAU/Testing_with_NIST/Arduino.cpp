#include "Arduino.h"

using namespace std;


int Serial::println(myType buf){
  cout << buf << endl;
  return strlen(buf);
}
