#ifndef __KPM_ECODES__
#define __KPM_ECODES__

#include <pthread.h>

namespace kpmeans { namespace errors {
    class ecodes {
        public:
            static void test_code(int ecode) {
                if (ecode == EINVAL) {
                    printf("\n\nERROR code: EINVAL\n\n");
                } else if (ecode == EBUSY) {
                    printf("\n\nERROR code: EBUSY\n\n");
                } else if (ecode == EAGAIN) {
                    printf("\n\nERROR code: EAGAIN\n\n");
                } else if (ecode == EDEADLK) {
                    printf("\n\nERROR code: EDEADLK\n\n");
                } else if (ecode == EPERM) {
                    printf("\n\nERROR code: EPERM\n\n");
                } else if (ecode == EINTR) {
                    printf("\n\nERROR code: EINTR\n\n");
                } else{
                    printf("\n\nUNKNOWN ERROR code!\n\n");
                }
            }
    };
}} // End namespace kpmeans::errors
#endif
