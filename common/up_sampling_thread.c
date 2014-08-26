/* Add by chenjie */
#include "up_sampling_thread.h"

void x264_threadpool_up_sampling()
{
	while( !up_sampling_arg.pThreadpool->exit )
	{
		x264_pthread_mutex_lock( &mutex_slice_encode_threads_finished );
		x264_pthread_cond_wait( &cv_start_up_sampling, &mutex_slice_encode_threads_finished );
		printf("---------------------------------\n");
		printf("---------------------------------\n");
		printf("---------------------------------\n");
		printf("---------------------------------\n");
		x264_pthread_mutex_unlock( &mutex_slice_encode_threads_finished );
		x264_pthread_cond_broadcast( &cv_slice_encode_thread_wakeup );
	}
}