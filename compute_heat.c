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
// reach the next heat index. So g_thresholds[3,2] == 35 means that for a
// relative humidity of 15% and a temperature between 35 deg C and < 42.2C,
// The heat index is at level 2 of [0:4].
const float g_thresholds[NPERCENTILES][NINDICES] = {
 {  28.89, 36.67, 45.56, 999.9 },  // 5% humidity
 {  28.33, 35.56, 43.89, 999.9 },  // 10%
 {  27.78, 35.00, 42.22, 51.67 },  // 15%
 {  27.78, 34.44, 41.11, 49.44 },  // 20%
 {  27.78, 33.89, 40.00, 47.22 },  // 25%
 {  27.22, 33.33, 38.89, 45.56 },  // 30%
 {  26.67, 32.78, 37.78, 43.89 },  // 35%
 {  26.67, 32.22, 37.22, 42.78 },  // 40%
 {  26.67, 31.67, 36.11, 41.11 },  // 45%
 {  26.67, 31.11, 35.00, 40.00 },  // 50%
 {  26.67, 30.56, 34.44, 54.44 },  // 55%
 {  26.67, 30.00, 33.89, 38.33 },  // 60%
 {  26.67, 29.44, 32.78, 37.22 },  // 65%
 {  26.67, 28.89, 32.22, 36.11 },  // 70%
 {  26.67, 28.89, 31.67, 35.56 },  // 75%
 {  26.67, 28.33, 31.11, 35.00 },  // 80%
 {  26.67, 27.78, 30.56, 33.89 },  // 85%
 {  26.67, 27.78, 30.00, 33.33 },  // 90%
 {  26.67, 27.78, 30.00, 32.78 },  // 95%
 {  26.67, 27.22, 29.44, 32.22 }   // 100%
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

