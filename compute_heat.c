// TODO: Run with valgrind!
/*
 * Combine humidity and teperature to come up with a heat risk category
 * Receives as input filenames for humidity and temperature data.
 * Assumes these files are sorted and matching exactly on their coordinates, i.e.,
 * line i has the same lat/long in both files.
 * Outputs in CSV format to standard output.
 *
 * Compile with: gcc -O3 -o compute_heat compute_heat.c
 */

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NINDICES 4
#define NPERCENTILES 20

// Lookup table: Each row represents increasing amount of humidity, from 5% to
// 100% in 5% increments. Each column represents the minimum temperature to
// reach the next heat index.
const float g_thresholds[NPERCENTILES][NINDICES] = {
 { 84, 98, 114, 999 },  // 5% humidity
 { 83, 96, 111, 999 },  // 10%
 { 82, 95, 108, 125 },  // 15%
 { 82, 94, 106, 121 },  // 20%
 { 82, 93, 104, 117 },  // 25%
 { 81, 92, 102, 114 },  // 30%
 { 80, 91, 100, 111 },  // 35%
 { 80, 90, 99, 109 },   // 40%
 { 80, 89, 97, 106 },   // 45%
 { 80, 88, 95, 104 },   // 50%
 { 80, 87, 94, 102 },   // 55%
 { 80, 86, 93, 101 },   // 60%
 { 80, 85, 91, 99 },    // 65%
 { 80, 84, 90, 97 },    // 70%
 { 80, 84, 89, 96 },    // 75%
 { 80, 83, 88, 95 },    // 80%
 { 80, 82, 87, 93 },    // 85%
 { 80, 82, 86, 92 },    // 90%
 { 80, 82, 86, 91 },    // 95%
 { 80, 81, 85, 90 }     // 100%

};

int compute_heat_index(float humidity, float temp)
{
//  assert(humidity >= 0 && humidity <= 1);
  assert(temp >= -60 && temp <= 60);

  // The formula for rounding humidity values takes all values
  // below 0.025 to 5%, and every other humidity value is rounded
  // to the nearest 5%. The first index (row = 0) is 5%.
  int row = (int)(0.5 + humidity * NPERCENTILES) - 1;
  if (row < 0) row = 0;  // For those humidities below 2.5%

  // For temperature, we just convert from C to F and round:
  temp = (int)(temp * 9. / 5. + 32 + 0.5);

  assert(row < NPERCENTILES);

  int i = 0;
  for (; i < NINDICES; ++i) {
    if (temp < g_thresholds[row][i]) {
      break;
    }
  }
  return i;
}

void compute_row(char* humline, char* templine)
{
  char *savehum, *savetemp, *tstr, *hstr;

  // Get coordinates from both rows, assert equal
  float lontemp = atof(strtok_r(templine, ",", &savetemp));
  float lattemp = atof(strtok_r(NULL, ",", &savetemp));
  float lathum = atof(strtok_r(humline, ",", &savehum));
  float lonhum = atof(strtok_r(NULL, ",", &savehum));
  assert(lontemp == lonhum);
  assert(lattemp == lathum);
  printf("%f,%f", lontemp, lattemp);

// Loop over rest of values up to templine max
  while ((tstr = strtok_r(NULL, ",", &savetemp))) {
    hstr = strtok_r(NULL, ",", &savehum);
    assert(hstr);

    // compute heat index for each value and output:
    const int heat_index = compute_heat_index(atof(hstr), atof(tstr));
    printf(",%d", heat_index);
  }

  puts("\n");
}

int main(int argc, char** argv)
{
  if (argc < 3) {
    fprintf(stderr, "Requred parameters: humidity CSV filename and temp CSV filename\n");
    exit(-1);
  }

  FILE *fp_hum = fopen(argv[1], "r");
  FILE *fp_temp = fopen(argv[2], "r");
  if (NULL == fp_hum || NULL == fp_temp) {
    fprintf(stderr, "Can't open input file %s or %s\n", argv[1], argv[2]);
    exit(-2);
  }

  char *hum_line = NULL, *temp_line = NULL;
  size_t len_hum, len_temp, dummy;

  // Output first line of temp file as CSV header. discard first line of hum file
  dummy = getline(&hum_line, &len_hum, fp_hum);
  dummy = getline(&temp_line, &len_temp, fp_temp);
  puts(temp_line);

  while (getline(&hum_line, &len_hum, fp_hum) != -1) {
    ssize_t read = getline(&temp_line, &len_temp, fp_temp);
    compute_row(hum_line, temp_line);
    assert(read > 0);
  }

  assert(getline(&temp_line, &len_temp, fp_temp) == -1);
  if (hum_line) free(hum_line);
  if (temp_line) free(temp_line);

  return 0;
}

