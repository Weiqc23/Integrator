#include "untitled.h"

ExtU_untitled_T untitled_U;
ExtY_untitled_T untitled_Y;
static RT_MODEL_untitled_T untitled_M_;
RT_MODEL_untitled_T *const untitled_M = &untitled_M_;
void untitled_step(void)
{
  if (untitled_U.In1 > 0.0) {
    untitled_Y.Out1 = 1.5 * untitled_U.In1;
  } else {
    untitled_Y.Out1 = 2.0 * untitled_U.In1;
  }
}

void untitled_initialize(void)
{
}

void untitled_terminate(void)
{
}
