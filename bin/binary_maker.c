#include <stdio.h>
#include <stdint.h>
#define LEN 58

int main() {
  // According to 6.3.5:
  uint8_t x[LEN] = {2,2,0,1,0,2,0,1,2,1,2,0,1,2,1,0,0,1,0,0,0,1,2,0,2,1,0,1,2,0,1,0,2,0,1,2,0,1,0,2,1,2,0,1,1,0,2,1,2,1,0,0,1,1,1,0};
  FILE *fh = fopen ("sample_for_tuple.bin", "wb");
  if (fh != NULL) {
    for(int i=0;i<LEN;i++){
      fwrite (&x[i], sizeof(uint8_t), 1, fh);
    }
    fclose (fh);
  }
  return 0;
}
