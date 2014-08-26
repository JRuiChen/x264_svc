#include "common/common.h"


/* Add by chenjie */
void x264_threadpool_up_sampling();

struct x264_thread_up_sampling_arg_t
{
    x264_threadpool_t *pThreadpool;
    x264_t *pH;
};
struct x264_thread_up_sampling_arg_t up_sampling_arg;

int i_slice_encode_threads_finished = 0;

x264_pthread_mutex_t mutex_slice_encode_threads_finished;
x264_pthread_cond_t cv_start_up_sampling;

x264_pthread_mutex_t mutex_slice_encode_thread_wakeup;
x264_pthread_cond_t cv_slice_encode_thread_wakeup;
