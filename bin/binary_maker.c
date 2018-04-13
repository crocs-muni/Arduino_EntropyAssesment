#include <stdio.h>

int main() {
  // According to 6.3.5:
  char x[59] = {2,2,0,1,0,2,0,1,2,1,2,0,1,2,1,0,0,1,0,0,0,1,2,0,2,1,0,1,2,0,1,0,2,0,1,2,0,1,0,2,1,2,0,1,1,0,2,1,2,1,0,0,1,1,1,1,0};
  FILE *fh = fopen ("sample_for_tuple.bin", "wb");
  if (fh != NULL) {
    for(int i=0;i<59;i++){
      fwrite (&x[i], sizeof (char), 1, fh);
    }
    fclose (fh);
  }
  return 0;
}
