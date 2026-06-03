#ifndef untitled_h_
#define untitled_h_
#ifndef untitled_COMMON_INCLUDES_
#define untitled_COMMON_INCLUDES_
#include "rtwtypes.h"
#endif

#include "untitled_types.h"

#ifndef rtmGetErrorStatus
#define rtmGetErrorStatus(rtm)         ((rtm)->errorStatus)
#endif

#ifndef rtmSetErrorStatus
#define rtmSetErrorStatus(rtm, val)    ((rtm)->errorStatus = (val))
#endif

typedef struct {
  real_T In1;
} ExtU_untitled_T;

typedef struct {
  real_T Out1;
} ExtY_untitled_T;

struct tag_RTM_untitled_T {
  const char_T * volatile errorStatus;
};

extern ExtU_untitled_T untitled_U;
extern ExtY_untitled_T untitled_Y;
extern void untitled_initialize(void);
extern void untitled_step(void);
extern void untitled_terminate(void);
extern RT_MODEL_untitled_T *const untitled_M;

#endif

