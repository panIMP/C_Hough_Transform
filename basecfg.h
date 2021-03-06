/*   Do *not* directly modify this file.  It was    */
/*   generated by the Configuration Tool; any  */
/*   changes risk being overwritten.                */

/* INPUT base.cdb */

/*  Include Header Files  */
#include <std.h>
#include <prd.h>
#include <hst.h>
#include <swi.h>
#include <tsk.h>
#include <log.h>
#include <sem.h>
#include <mbx.h>
#include <sts.h>

#ifdef __cplusplus
extern "C" {
#endif

extern far PRD_Obj prdStack;
extern far HST_Obj RTA_fromHost;
extern far HST_Obj RTA_toHost;
extern far SWI_Obj PRD_swi;
extern far SWI_Obj KNL_swi;
extern far TSK_Obj TSK_idle;
extern far TSK_Obj TSK_vpss;
extern far TSK_Obj TSK_ndk;
extern far TSK_Obj TSK_msge;
extern far TSK_Obj TSK_uart;
extern far TSK_Obj TSK_can;
extern far TSK_Obj TSK_sendTime;
extern far TSK_Obj TSK_sendTemper;
extern far LOG_Obj LOG_system;
extern far LOG_Obj trace;
extern far LOG_Obj DVTEvent_Log;
extern far SEM_Obj msgeSem;
extern far SEM_Obj vidSem;
extern far SEM_Obj oneshotSem;
extern far SEM_Obj uartSem;
extern far SEM_Obj canSem;
extern far SEM_Obj vpfeSem;
extern far SEM_Obj vpbeSem;
extern far SEM_Obj timeSem;
extern far SEM_Obj temperSem;
extern far SEM_Obj uartSourceSem;
extern far MBX_Obj msgMbx;
extern far MBX_Obj netMbx;
extern far MBX_Obj canMbx;
extern far MBX_Obj algMbx;
extern far STS_Obj IDL_busyObj;

extern IRAM;
extern DDR2;
#ifdef __cplusplus
}
#endif /* extern "C" */
