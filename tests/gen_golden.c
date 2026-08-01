/**
 * @file      gen_golden.c
 * @brief     Write tests/golden.dat from the CURRENT build.
 *
 * This locks 2D/3D behaviour before the dimension refactor. Only ever run it
 * from a build you have reason to trust, and review the diff: a regenerated
 * golden file is a claim that the change in behaviour was intended.
 *
 *   ./gen_golden > tests/golden.dat
 */

#include "sweep.h"

static void emit(const double *rec, void *ctx)
{
  int i;
  (void) ctx;
  for (i = 0; i < VOFI_NREC; i++)
    printf("%s%.17g", i ? " " : "", rec[i]);
  printf("\n");
}

int main(void)
{
  int c;

  printf("# vofi golden data: %d cases, %d values per cell\n",
         VOFI_NCASE_23 + VOFI_NCASE_4, VOFI_NREC);
  for (c = 0; c < VOFI_NCASE_23; c++) {
    printf("# case %s\n", vofi_cases_23[c].name);
    vofi_run_case(&vofi_cases_23[c], emit, NULL);
  }
  for (c = 0; c < VOFI_NCASE_4; c++) {
    printf("# case %s\n", vofi_cases_4[c].name);
    vofi_run_case(&vofi_cases_4[c], emit, NULL);
  }
  return 0;
}
