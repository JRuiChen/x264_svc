/*****************************************************************************
 * frame.c: frame handling
 *****************************************************************************
 * Copyright (C) 2003-2014 x264 project
 *
 * Authors: Laurent Aimar <fenrir@via.ecp.fr>
 *          Loren Merritt <lorenm@u.washington.edu>
 *          Jason Garrett-Glaser <darkshikari@gmail.com>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02111, USA.
 *
 * This program is also available under a commercial proprietary license.
 * For more information, contact us at licensing@x264.com.
 *****************************************************************************/

#include "common.h"
//#include <malloc.h>
//#include<string.h>

static int align_stride( int x, int align, int disalign )
{
    x = ALIGN( x, align );
    if( !(x&(disalign-1)) )
        x += align;
    return x;
}

static int align_plane_size( int x, int disalign )
{
    if( !(x&(disalign-1)) )
        x += 128;
    return x;
}

static int x264_frame_internal_csp( int external_csp )
{
    switch( external_csp & X264_CSP_MASK )
    {
        case X264_CSP_NV12:
        case X264_CSP_I420:
        case X264_CSP_YV12:
            return X264_CSP_NV12;
        case X264_CSP_NV16:
        case X264_CSP_I422:
        case X264_CSP_YV16:
  //      case X264_CSP_V210:
            return X264_CSP_NV16;
        case X264_CSP_I444:
        case X264_CSP_YV24:
        case X264_CSP_BGR:
        case X264_CSP_BGRA:
        case X264_CSP_RGB:
            return X264_CSP_I444;
        default:
            return X264_CSP_NONE;
    }
}

static x264_frame_t *x264_frame_new( x264_t *h, int b_fdec )
{
    x264_frame_t *frame;
    int i_csp = x264_frame_internal_csp( h->param.i_csp );
    int i_mb_count = h->mb.i_mb_count;
    int i_stride, i_width, i_lines, luma_plane_count;
    int i_padv = PADV << PARAM_INTERLACED;
    int align = 16;
#if ARCH_X86 || ARCH_X86_64
    if( h->param.cpu&X264_CPU_CACHELINE_64 )
        align = 64;
    else if( h->param.cpu&X264_CPU_CACHELINE_32 || h->param.cpu&X264_CPU_AVX2 )
        align = 32;
#endif
#if ARCH_PPC
    int disalign = 1<<9;
#else
    int disalign = 1<<10;
#endif

    CHECKED_MALLOCZERO( frame, sizeof(x264_frame_t) );
    PREALLOC_INIT

    /* allocate frame data (+64 for extra data for me) */
    i_width  = h->mb.i_mb_width*16;
    i_lines  = h->mb.i_mb_height*16;
    i_stride = align_stride( i_width + 2*PADH, align, disalign );
    dst_s = align_stride( i_width*2 + 2*PADH, align, disalign );  //add by chenjie

    if( i_csp == X264_CSP_NV12 || i_csp == X264_CSP_NV16 )
    {
        luma_plane_count = 1;
        frame->i_plane = 2;
        for( int i = 0; i < 2; i++ )
        {
            frame->i_width[i] = i_width >> i;
            frame->i_lines[i] = i_lines >> (i && i_csp == X264_CSP_NV12);
            frame->i_stride[i] = i_stride;
        }
    }
    else if( i_csp == X264_CSP_I444 )
    {
        luma_plane_count = 3;
        frame->i_plane = 3;
        for( int i = 0; i < 3; i++ )
        {
            frame->i_width[i] = i_width;
            frame->i_lines[i] = i_lines;
            frame->i_stride[i] = i_stride;
        }
    }
    else
        goto fail;

    frame->i_csp = i_csp;
    frame->i_width_lowres = frame->i_width[0]/2;
    frame->i_lines_lowres = frame->i_lines[0]/2;
    frame->i_stride_lowres = align_stride( frame->i_width_lowres + 2*PADH, align, disalign<<1 );

    for( int i = 0; i < h->param.i_bframe + 2; i++ )
        for( int j = 0; j < h->param.i_bframe + 2; j++ )
            PREALLOC( frame->i_row_satds[i][j], i_lines/16 * sizeof(int) );

    frame->i_poc = -1;
    frame->i_type = X264_TYPE_AUTO;
    frame->i_qpplus1 = X264_QP_AUTO;
    frame->i_pts = -1;
    frame->i_frame = -1;
    frame->i_frame_num = -1;
    frame->i_lines_completed = -1;
    frame->b_fdec = b_fdec;
    frame->i_pic_struct = PIC_STRUCT_AUTO;
    frame->i_field_cnt = -1;
    frame->i_duration =
    frame->i_cpb_duration =
    frame->i_dpb_output_delay =
    frame->i_cpb_delay = 0;
    frame->i_coded_fields_lookahead =
    frame->i_cpb_delay_lookahead = -1;

    frame->orig = frame;

    if( i_csp == X264_CSP_NV12 || i_csp == X264_CSP_NV16 )
    {
        int chroma_padv = i_padv >> (i_csp == X264_CSP_NV12);
        int chroma_plane_size = (frame->i_stride[1] * (frame->i_lines[1] + 2*chroma_padv));
        PREALLOC( frame->buffer[1], chroma_plane_size * sizeof(pixel) );
        if( PARAM_INTERLACED )
            PREALLOC( frame->buffer_fld[1], chroma_plane_size * sizeof(pixel) );
    }

    /* all 4 luma planes allocated together, since the cacheline split code
     * requires them to be in-phase wrt cacheline alignment. */

    for( int p = 0; p < luma_plane_count; p++ )
    {
        int luma_plane_size = align_plane_size( frame->i_stride[p] * (frame->i_lines[p] + 2*i_padv), disalign );
        if( h->param.analyse.i_subpel_refine && b_fdec )
        {
            /* FIXME: Don't allocate both buffers in non-adaptive MBAFF. */
            PREALLOC( frame->buffer[p], 4*luma_plane_size * sizeof(pixel) );
            if( PARAM_INTERLACED )
                PREALLOC( frame->buffer_fld[p], 4*luma_plane_size * sizeof(pixel) );
        }
        else
        {
            PREALLOC( frame->buffer[p], luma_plane_size * sizeof(pixel) );
            if( PARAM_INTERLACED )
                PREALLOC( frame->buffer_fld[p], luma_plane_size * sizeof(pixel) );
        }
    }

    frame->b_duplicate = 0;

    if( b_fdec ) /* fdec frame */
    {
        PREALLOC( frame->mb_type, i_mb_count * sizeof(int8_t) );
        PREALLOC( frame->mb_partition, i_mb_count * sizeof(uint8_t) );
        PREALLOC( frame->mv[0], 2*16 * i_mb_count * sizeof(int16_t) );
        PREALLOC( frame->mv16x16, 2*(i_mb_count+1) * sizeof(int16_t) );
        PREALLOC( frame->ref[0], 4 * i_mb_count * sizeof(int8_t) );
        if( h->param.i_bframe )
        {
            PREALLOC( frame->mv[1], 2*16 * i_mb_count * sizeof(int16_t) );
            PREALLOC( frame->ref[1], 4 * i_mb_count * sizeof(int8_t) );
        }
        else
        {
            frame->mv[1]  = NULL;
            frame->ref[1] = NULL;
        }
        PREALLOC( frame->i_row_bits, i_lines/16 * sizeof(int) );
        PREALLOC( frame->f_row_qp, i_lines/16 * sizeof(float) );
        PREALLOC( frame->f_row_qscale, i_lines/16 * sizeof(float) );
        if( h->param.analyse.i_me_method >= X264_ME_ESA )
            PREALLOC( frame->buffer[3], frame->i_stride[0] * (frame->i_lines[0] + 2*i_padv) * sizeof(uint16_t) << h->frames.b_have_sub8x8_esa );
        if( PARAM_INTERLACED )
            PREALLOC( frame->field, i_mb_count * sizeof(uint8_t) );
        if( h->param.analyse.b_mb_info )
            PREALLOC( frame->effective_qp, i_mb_count * sizeof(uint8_t) );
    }
    else /* fenc frame */
    {
        if( h->frames.b_have_lowres )
        {
            int luma_plane_size = align_plane_size( frame->i_stride_lowres * (frame->i_lines[0]/2 + 2*PADV), disalign );

            PREALLOC( frame->buffer_lowres[0], 4 * luma_plane_size * sizeof(pixel) );

            for( int j = 0; j <= !!h->param.i_bframe; j++ )
                for( int i = 0; i <= h->param.i_bframe; i++ )
                {
                    PREALLOC( frame->lowres_mvs[j][i], 2*h->mb.i_mb_count*sizeof(int16_t) );
                    PREALLOC( frame->lowres_mv_costs[j][i], h->mb.i_mb_count*sizeof(int) );
                }
            PREALLOC( frame->i_propagate_cost, (i_mb_count+7) * sizeof(uint16_t) );
            for( int j = 0; j <= h->param.i_bframe+1; j++ )
                for( int i = 0; i <= h->param.i_bframe+1; i++ )
                    PREALLOC( frame->lowres_costs[j][i], (i_mb_count+3) * sizeof(uint16_t) );

        }
        if( h->param.rc.i_aq_mode )
        {
            PREALLOC( frame->f_qp_offset, h->mb.i_mb_count * sizeof(float) );
            PREALLOC( frame->f_qp_offset_aq, h->mb.i_mb_count * sizeof(float) );
            if( h->frames.b_have_lowres )
                PREALLOC( frame->i_inv_qscale_factor, (h->mb.i_mb_count+3) * sizeof(uint16_t) );
        }
    }

    PREALLOC_END( frame->base );

    if( i_csp == X264_CSP_NV12 || i_csp == X264_CSP_NV16 )
    {
        int chroma_padv = i_padv >> (i_csp == X264_CSP_NV12);
        frame->plane[1] = frame->buffer[1] + frame->i_stride[1] * chroma_padv + PADH;
        if( PARAM_INTERLACED )
            frame->plane_fld[1] = frame->buffer_fld[1] + frame->i_stride[1] * chroma_padv + PADH;
    }

    for( int p = 0; p < luma_plane_count; p++ )
    {
        int luma_plane_size = align_plane_size( frame->i_stride[p] * (frame->i_lines[p] + 2*i_padv), disalign );
        if( h->param.analyse.i_subpel_refine && b_fdec )
        {
            for( int i = 0; i < 4; i++ )
            {
                frame->filtered[p][i] = frame->buffer[p] + i*luma_plane_size + frame->i_stride[p] * i_padv + PADH;
                frame->filtered_fld[p][i] = frame->buffer_fld[p] + i*luma_plane_size + frame->i_stride[p] * i_padv + PADH;
            }
            frame->plane[p] = frame->filtered[p][0];
            frame->plane_fld[p] = frame->filtered_fld[p][0];
        }
        else
        {
            frame->filtered[p][0] = frame->plane[p] = frame->buffer[p] + frame->i_stride[p] * i_padv + PADH;
            frame->filtered_fld[p][0] = frame->plane_fld[p] = frame->buffer_fld[p] + frame->i_stride[p] * i_padv + PADH;
        }
    }

    if( b_fdec )
    {
        M32( frame->mv16x16[0] ) = 0;
        frame->mv16x16++;

        if( h->param.analyse.i_me_method >= X264_ME_ESA )
            frame->integral = (uint16_t*)frame->buffer[3] + frame->i_stride[0] * i_padv + PADH;
    }
    else
    {
        if( h->frames.b_have_lowres )
        {
            int luma_plane_size = align_plane_size( frame->i_stride_lowres * (frame->i_lines[0]/2 + 2*PADV), disalign );
            for( int i = 0; i < 4; i++ )
                frame->lowres[i] = frame->buffer_lowres[0] + (frame->i_stride_lowres * PADV + PADH) + i * luma_plane_size;

            for( int j = 0; j <= !!h->param.i_bframe; j++ )
                for( int i = 0; i <= h->param.i_bframe; i++ )
                    memset( frame->lowres_mvs[j][i], 0, 2*h->mb.i_mb_count*sizeof(int16_t) );

            frame->i_intra_cost = frame->lowres_costs[0][0];
            memset( frame->i_intra_cost, -1, (i_mb_count+3) * sizeof(uint16_t) );

            if( h->param.rc.i_aq_mode )
                /* shouldn't really be initialized, just silences a valgrind false-positive in x264_mbtree_propagate_cost_sse2 */
                memset( frame->i_inv_qscale_factor, 0, (h->mb.i_mb_count+3) * sizeof(uint16_t) );
        }
    }

    if( x264_pthread_mutex_init( &frame->mutex, NULL ) )
        goto fail;
    if( x264_pthread_cond_init( &frame->cv, NULL ) )
        goto fail;

#if HAVE_OPENCL
    frame->opencl.ocl = h->opencl.ocl;
#endif

    return frame;

fail:
    x264_free( frame );
    return NULL;
}

void x264_frame_delete( x264_frame_t *frame )
{
    /* Duplicate frames are blank copies of real frames (including pointers),
     * so freeing those pointers would cause a double free later. */
    if( !frame->b_duplicate )
    {
        x264_free( frame->base );

        if( frame->param && frame->param->param_free )
            frame->param->param_free( frame->param );
        if( frame->mb_info_free )
            frame->mb_info_free( frame->mb_info );
        if( frame->extra_sei.sei_free )
        {
            for( int i = 0; i < frame->extra_sei.num_payloads; i++ )
                frame->extra_sei.sei_free( frame->extra_sei.payloads[i].payload );
            frame->extra_sei.sei_free( frame->extra_sei.payloads );
        }
        x264_pthread_mutex_destroy( &frame->mutex );
        x264_pthread_cond_destroy( &frame->cv );
#if HAVE_OPENCL
        x264_opencl_frame_delete( frame );
#endif
    }
    x264_free( frame );
}

static int get_plane_ptr( x264_t *h, x264_picture_t *src, uint8_t **pix, int *stride, int plane, int xshift, int yshift )
{
    int width = h->param.i_width >> xshift;
    int height = h->param.i_height >> yshift;
    *pix = src->img.plane[plane];
    *stride = src->img.i_stride[plane];
    if( src->img.i_csp & X264_CSP_VFLIP )
    {
        *pix += (height-1) * *stride;
        *stride = -*stride;
    }
    if( width > abs(*stride) )
    {
        x264_log( h, X264_LOG_ERROR, "Input picture width (%d) is greater than stride (%d)\n", width, *stride );
        return -1;
    }
    return 0;
}

#define get_plane_ptr(...) do{ if( get_plane_ptr(__VA_ARGS__) < 0 ) return -1; }while(0)

int x264_frame_copy_picture( x264_t *h, x264_frame_t *dst, x264_picture_t *src )
{
    int i_csp = src->img.i_csp & X264_CSP_MASK;
    if( dst->i_csp != x264_frame_internal_csp( i_csp ) )
    {
        x264_log( h, X264_LOG_ERROR, "Invalid input colorspace\n" );
        return -1;
    }

#if HIGH_BIT_DEPTH
    if( !(src->img.i_csp & X264_CSP_HIGH_DEPTH) )
    {
        x264_log( h, X264_LOG_ERROR, "This build of x264 requires high depth input. Rebuild to support 8-bit input.\n" );
        return -1;
    }
#else
    if( src->img.i_csp & X264_CSP_HIGH_DEPTH )
    {
        x264_log( h, X264_LOG_ERROR, "This build of x264 requires 8-bit input. Rebuild to support high depth input.\n" );
        return -1;
    }
#endif

 /*   if( BIT_DEPTH != 10 && i_csp == X264_CSP_V210 )
    {
        x264_log( h, X264_LOG_ERROR, "v210 input is only compatible with bit-depth of 10 bits\n" );
        return -1;
    }*/

    dst->i_type     = src->i_type;
    dst->i_qpplus1  = src->i_qpplus1;
    dst->i_pts      = dst->i_reordered_pts = src->i_pts;
    dst->param      = src->param;
    dst->i_pic_struct = src->i_pic_struct;
    dst->extra_sei  = src->extra_sei;
    dst->opaque     = src->opaque;
    dst->mb_info    = h->param.analyse.b_mb_info ? src->prop.mb_info : NULL;
    dst->mb_info_free = h->param.analyse.b_mb_info ? src->prop.mb_info_free : NULL;

    uint8_t *pix[3];
    int stride[3];
 /*   if( i_csp == X264_CSP_V210 )
    {
         stride[0] = src->img.i_stride[0];
         pix[0] = src->img.plane[0];

         h->mc.plane_copy_deinterleave_v210( dst->plane[0], dst->i_stride[0],
                                             dst->plane[1], dst->i_stride[1],
                                             (uint32_t *)pix[0], stride[0]/sizeof(uint32_t), h->param.i_width, h->param.i_height );
    }*/
    
//    else if( i_csp >= X264_CSP_BGR )
    if( i_csp >= X264_CSP_BGR )
    {
         stride[0] = src->img.i_stride[0];
         pix[0] = src->img.plane[0];
         if( src->img.i_csp & X264_CSP_VFLIP )
         {
             pix[0] += (h->param.i_height-1) * stride[0];
             stride[0] = -stride[0];
         }
         int b = i_csp==X264_CSP_RGB;
         h->mc.plane_copy_deinterleave_rgb( dst->plane[1+b], dst->i_stride[1+b],
                                            dst->plane[0], dst->i_stride[0],
                                            dst->plane[2-b], dst->i_stride[2-b],
                                            (pixel*)pix[0], stride[0]/sizeof(pixel), i_csp==X264_CSP_BGRA ? 4 : 3, h->param.i_width, h->param.i_height );
    }
    else
    {
        int v_shift = CHROMA_V_SHIFT;
        get_plane_ptr( h, src, &pix[0], &stride[0], 0, 0, 0 );//获得Y分量指针
        h->mc.plane_copy( dst->plane[0], dst->i_stride[0], (pixel*)pix[0],
                          stride[0]/sizeof(pixel), h->param.i_width, h->param.i_height );//复制数据
        if( i_csp == X264_CSP_NV12 || i_csp == X264_CSP_NV16 )
        {
            get_plane_ptr( h, src, &pix[1], &stride[1], 1, 0, v_shift );
            h->mc.plane_copy( dst->plane[1], dst->i_stride[1], (pixel*)pix[1],
                              stride[1]/sizeof(pixel), h->param.i_width, h->param.i_height>>v_shift );
        }
        else if( i_csp == X264_CSP_I420 || i_csp == X264_CSP_I422 || i_csp == X264_CSP_YV12 || i_csp == X264_CSP_YV16 )
        {
            int uv_swap = i_csp == X264_CSP_YV12 || i_csp == X264_CSP_YV16;//UV存储顺序是否交换
            get_plane_ptr( h, src, &pix[1], &stride[1], uv_swap ? 2 : 1, 1, v_shift );
            get_plane_ptr( h, src, &pix[2], &stride[2], uv_swap ? 1 : 2, 1, v_shift );
            h->mc.plane_copy_interleave( dst->plane[1], dst->i_stride[1],
                                         (pixel*)pix[1], stride[1]/sizeof(pixel),
                                         (pixel*)pix[2], stride[2]/sizeof(pixel),
                                         h->param.i_width>>1, h->param.i_height>>v_shift );
        }
        else //if( i_csp == X264_CSP_I444 || i_csp == X264_CSP_YV24 )
        {
            get_plane_ptr( h, src, &pix[1], &stride[1], i_csp==X264_CSP_I444 ? 1 : 2, 0, 0 );
            get_plane_ptr( h, src, &pix[2], &stride[2], i_csp==X264_CSP_I444 ? 2 : 1, 0, 0 );
            h->mc.plane_copy( dst->plane[1], dst->i_stride[1], (pixel*)pix[1],
                              stride[1]/sizeof(pixel), h->param.i_width, h->param.i_height );
            h->mc.plane_copy( dst->plane[2], dst->i_stride[2], (pixel*)pix[2],
                              stride[2]/sizeof(pixel), h->param.i_width, h->param.i_height );
        }
    }
    return 0;
}

static void ALWAYS_INLINE pixel_memset( pixel *dst, pixel *src, int len, int size )
{
    uint8_t *dstp = (uint8_t*)dst;
    uint32_t v1 = *src;
    uint32_t v2 = size == 1 ? v1 + (v1 <<  8) : M16( src );
    uint32_t v4 = size <= 2 ? v2 + (v2 << 16) : M32( src );
    int i = 0;
    len *= size;

    /* Align the input pointer if it isn't already */
    if( (intptr_t)dstp & (WORD_SIZE - 1) )
    {
        if( size <= 2 && ((intptr_t)dstp & 3) )
        {
            if( size == 1 && ((intptr_t)dstp & 1) )
                dstp[i++] = v1;
            if( (intptr_t)dstp & 2 )
            {
                M16( dstp+i ) = v2;
                i += 2;
            }
        }
        if( WORD_SIZE == 8 && (intptr_t)dstp & 4 )
        {
            M32( dstp+i ) = v4;
            i += 4;
        }
    }

    /* Main copy loop */
    if( WORD_SIZE == 8 )
    {
        uint64_t v8 = v4 + ((uint64_t)v4<<32);
        for( ; i < len - 7; i+=8 )
            M64( dstp+i ) = v8;
    }
    for( ; i < len - 3; i+=4 )
        M32( dstp+i ) = v4;

    /* Finish up the last few bytes */
    if( size <= 2 )
    {
        if( i < len - 1 )
        {
            M16( dstp+i ) = v2;
            i += 2;
        }
        if( size == 1 && i != len )
            dstp[i] = v1;
    }
}

static void ALWAYS_INLINE plane_expand_border( pixel *pix, int i_stride, int i_width, int i_height, int i_padh, int i_padv, int b_pad_top, int b_pad_bottom, int b_chroma )
{
#define PPIXEL(x, y) ( pix + (x) + (y)*i_stride )
    for( int y = 0; y < i_height; y++ )
    {
        /* left band */
        pixel_memset( PPIXEL(-i_padh, y), PPIXEL(0, y), i_padh>>b_chroma, sizeof(pixel)<<b_chroma );
        /* right band */
        pixel_memset( PPIXEL(i_width, y), PPIXEL(i_width-1-b_chroma, y), i_padh>>b_chroma, sizeof(pixel)<<b_chroma );
    }
    /* upper band */
    if( b_pad_top )
        for( int y = 0; y < i_padv; y++ )
            memcpy( PPIXEL(-i_padh, -y-1), PPIXEL(-i_padh, 0), (i_width+2*i_padh) * sizeof(pixel) );
    /* lower band */
    if( b_pad_bottom )
        for( int y = 0; y < i_padv; y++ )
            memcpy( PPIXEL(-i_padh, i_height+y), PPIXEL(-i_padh, i_height-1), (i_width+2*i_padh) * sizeof(pixel) );
#undef PPIXEL
}

void x264_frame_expand_border( x264_t *h, x264_frame_t *frame, int mb_y )
{
    int pad_top = mb_y == 0;
    int pad_bot = mb_y == h->mb.i_mb_height - (1 << SLICE_MBAFF);
    int b_start = mb_y == h->i_threadslice_start;
    int b_end   = mb_y == h->i_threadslice_end - (1 << SLICE_MBAFF);
    if( mb_y & SLICE_MBAFF )
        return;
    for( int i = 0; i < frame->i_plane; i++ )
    {
        int h_shift = i && CHROMA_H_SHIFT;
        int v_shift = i && CHROMA_V_SHIFT;
        int stride = frame->i_stride[i];
        int width = 16*h->mb.i_mb_width;
        int height = (pad_bot ? 16*(h->mb.i_mb_height - mb_y) >> SLICE_MBAFF : 16) >> v_shift;
        int padh = PADH;
        int padv = PADV >> v_shift;
        // buffer: 2 chroma, 3 luma (rounded to 4) because deblocking goes beyond the top of the mb
        if( b_end && !b_start )
            height += 4 >> (v_shift + SLICE_MBAFF);
        pixel *pix;
        int starty = 16*mb_y - 4*!b_start;
        if( SLICE_MBAFF )
        {
            // border samples for each field are extended separately
            pix = frame->plane_fld[i] + (starty*stride >> v_shift);
            plane_expand_border( pix, stride*2, width, height, padh, padv, pad_top, pad_bot, h_shift );
            plane_expand_border( pix+stride, stride*2, width, height, padh, padv, pad_top, pad_bot, h_shift );

            height = (pad_bot ? 16*(h->mb.i_mb_height - mb_y) : 32) >> v_shift;
            if( b_end && !b_start )
                height += 4 >> v_shift;
            pix = frame->plane[i] + (starty*stride >> v_shift);
            plane_expand_border( pix, stride, width, height, padh, padv, pad_top, pad_bot, h_shift );
        }
        else
        {
            pix = frame->plane[i] + (starty*stride >> v_shift);
            plane_expand_border( pix, stride, width, height, padh, padv, pad_top, pad_bot, h_shift );
        }
    }
}

void x264_frame_expand_border_filtered( x264_t *h, x264_frame_t *frame, int mb_y, int b_end )
{
    /* during filtering, 8 extra pixels were filtered on each edge,
     * but up to 3 of the horizontal ones may be wrong.
       we want to expand border from the last filtered pixel */
    int b_start = !mb_y;
    int width = 16*h->mb.i_mb_width + 8;
    int height = b_end ? (16*(h->mb.i_mb_height - mb_y) >> SLICE_MBAFF) + 16 : 16;
    int padh = PADH - 4;
    int padv = PADV - 8;
    for( int p = 0; p < (CHROMA444 ? 3 : 1); p++ )
        for( int i = 1; i < 4; i++ )
        {
            int stride = frame->i_stride[p];
            // buffer: 8 luma, to match the hpel filter
            pixel *pix;
            if( SLICE_MBAFF )
            {
                pix = frame->filtered_fld[p][i] + (16*mb_y - 16) * stride - 4;
                plane_expand_border( pix, stride*2, width, height, padh, padv, b_start, b_end, 0 );
                plane_expand_border( pix+stride, stride*2, width, height, padh, padv, b_start, b_end, 0 );
            }

            pix = frame->filtered[p][i] + (16*mb_y - 8) * stride - 4;
            plane_expand_border( pix, stride, width, height << SLICE_MBAFF, padh, padv, b_start, b_end, 0 );
        }
}

void x264_frame_expand_border_lowres( x264_frame_t *frame )
{
    for( int i = 0; i < 4; i++ )
        plane_expand_border( frame->lowres[i], frame->i_stride_lowres, frame->i_width_lowres, frame->i_lines_lowres, PADH, PADV, 1, 1, 0 );
}

void x264_frame_expand_border_chroma( x264_t *h, x264_frame_t *frame, int plane )
{
    int v_shift = CHROMA_V_SHIFT;
    plane_expand_border( frame->plane[plane], frame->i_stride[plane], 16*h->mb.i_mb_width, 16*h->mb.i_mb_height>>v_shift,
                         PADH, PADV>>v_shift, 1, 1, CHROMA_H_SHIFT );
}

void x264_frame_expand_border_mod16( x264_t *h, x264_frame_t *frame )
{
    for( int i = 0; i < frame->i_plane; i++ )
    {
        int i_width = h->param.i_width;
        int h_shift = i && CHROMA_H_SHIFT;
        int v_shift = i && CHROMA_V_SHIFT;
        int i_height = h->param.i_height >> v_shift;
        int i_padx = (h->mb.i_mb_width * 16 - h->param.i_width);
        int i_pady = (h->mb.i_mb_height * 16 - h->param.i_height) >> v_shift;

        if( i_padx )
        {
            for( int y = 0; y < i_height; y++ )
                pixel_memset( &frame->plane[i][y*frame->i_stride[i] + i_width],
                              &frame->plane[i][y*frame->i_stride[i] + i_width - 1-h_shift],
                              i_padx>>h_shift, sizeof(pixel)<<h_shift );
        }
        if( i_pady )
        {
            for( int y = i_height; y < i_height + i_pady; y++ )
                memcpy( &frame->plane[i][y*frame->i_stride[i]],
                        &frame->plane[i][(i_height-(~y&PARAM_INTERLACED)-1)*frame->i_stride[i]],
                        (i_width + i_padx) * sizeof(pixel) );
        }
    }
}

void x264_expand_border_mbpair( x264_t *h, int mb_x, int mb_y )
{
    for( int i = 0; i < h->fenc->i_plane; i++ )
    {
        int v_shift = i && CHROMA_V_SHIFT;
        int stride = h->fenc->i_stride[i];
        int height = h->param.i_height >> v_shift;
        int pady = (h->mb.i_mb_height * 16 - h->param.i_height) >> v_shift;
        pixel *fenc = h->fenc->plane[i] + 16*mb_x;
        for( int y = height; y < height + pady; y++ )
            memcpy( fenc + y*stride, fenc + (height-1)*stride, 16*sizeof(pixel) );
    }
}

/* threading */
void x264_frame_cond_broadcast( x264_frame_t *frame, int i_lines_completed )
{
    x264_pthread_mutex_lock( &frame->mutex );
    frame->i_lines_completed = i_lines_completed;
    x264_pthread_cond_broadcast( &frame->cv );
    x264_pthread_mutex_unlock( &frame->mutex );
}

void x264_frame_cond_wait( x264_frame_t *frame, int i_lines_completed )
{
    x264_pthread_mutex_lock( &frame->mutex );
    while( frame->i_lines_completed < i_lines_completed )
        x264_pthread_cond_wait( &frame->cv, &frame->mutex );
    x264_pthread_mutex_unlock( &frame->mutex );
}

void x264_threadslice_cond_broadcast( x264_t *h, int pass )
{
    x264_pthread_mutex_lock( &h->mutex );
    h->i_threadslice_pass = pass;
    if( pass > 0 )
        x264_pthread_cond_broadcast( &h->cv );
    x264_pthread_mutex_unlock( &h->mutex );
}

void x264_threadslice_cond_wait( x264_t *h, int pass )
{
    x264_pthread_mutex_lock( &h->mutex );
    while( h->i_threadslice_pass < pass )
        x264_pthread_cond_wait( &h->cv, &h->mutex );
    x264_pthread_mutex_unlock( &h->mutex );
}

int x264_frame_new_slice( x264_t *h, x264_frame_t *frame )
{
    if( h->param.i_slice_count_max )
    {
        int slice_count;
        if( h->param.b_sliced_threads )
            slice_count = x264_pthread_fetch_and_add( &frame->i_slice_count, 1, &frame->mutex );
        else
            slice_count = frame->i_slice_count++;
        if( slice_count >= h->param.i_slice_count_max )
            return -1;
    }
    return 0;
}

/* list operators */

void x264_frame_push( x264_frame_t **list, x264_frame_t *frame )
{
    int i = 0;
    while( list[i] ) i++;
    list[i] = frame;
}

x264_frame_t *x264_frame_pop( x264_frame_t **list )
{
    x264_frame_t *frame;
    int i = 0;
    assert( list[0] );
    while( list[i+1] ) i++;
    frame = list[i];
    list[i] = NULL;
    return frame;
}

void x264_frame_unshift( x264_frame_t **list, x264_frame_t *frame )
{
    int i = 0;
    while( list[i] ) i++;
    while( i-- )
        list[i+1] = list[i];
    list[0] = frame;
}

x264_frame_t *x264_frame_shift( x264_frame_t **list )
{
    x264_frame_t *frame = list[0];
    int i;
    for( i = 0; list[i]; i++ )
        list[i] = list[i+1];
    assert(frame);
    return frame;
}

void x264_frame_push_unused( x264_t *h, x264_frame_t *frame )
{
    assert( frame->i_reference_count > 0 );
    frame->i_reference_count--;
    if( frame->i_reference_count == 0 )
        x264_frame_push( h->frames.unused[frame->b_fdec], frame );
}

x264_frame_t *x264_frame_pop_unused( x264_t *h, int b_fdec )
{
    x264_frame_t *frame;
    if( h->frames.unused[b_fdec][0] )
        frame = x264_frame_pop( h->frames.unused[b_fdec] );
    else
        frame = x264_frame_new( h, b_fdec );
    if( !frame )
        return NULL;
    frame->b_last_minigop_bframe = 0;
    frame->i_reference_count = 1;
    frame->b_intra_calculated = 0;
    frame->b_scenecut = 1;
    frame->b_keyframe = 0;
    frame->b_corrupt = 0;
    frame->i_slice_count = h->param.b_sliced_threads ? h->param.i_threads : 1;

    memset( frame->weight, 0, sizeof(frame->weight) );
    memset( frame->f_weighted_cost_delta, 0, sizeof(frame->f_weighted_cost_delta) );

    return frame;
}

void x264_frame_push_blank_unused( x264_t *h, x264_frame_t *frame )
{
    assert( frame->i_reference_count > 0 );
    frame->i_reference_count--;
    if( frame->i_reference_count == 0 )
        x264_frame_push( h->frames.blank_unused, frame );
}

x264_frame_t *x264_frame_pop_blank_unused( x264_t *h )
{
    x264_frame_t *frame;
    if( h->frames.blank_unused[0] )
        frame = x264_frame_pop( h->frames.blank_unused );
    else
        frame = x264_malloc( sizeof(x264_frame_t) );
    if( !frame )
        return NULL;
    frame->b_duplicate = 1;
    frame->i_reference_count = 1;
    return frame;
}

void x264_weight_scale_plane( x264_t *h, pixel *dst, intptr_t i_dst_stride, pixel *src, intptr_t i_src_stride,
                              int i_width, int i_height, x264_weight_t *w )
{
    /* Weight horizontal strips of height 16. This was found to be the optimal height
     * in terms of the cache loads. */
    while( i_height > 0 )
    {
        int x;
        for( x = 0; x < i_width-8; x += 16 )
            w->weightfn[16>>2]( dst+x, i_dst_stride, src+x, i_src_stride, w, X264_MIN( i_height, 16 ) );
        if( x < i_width )
            w->weightfn[ 8>>2]( dst+x, i_dst_stride, src+x, i_src_stride, w, X264_MIN( i_height, 16 ) );
        i_height -= 16;
        dst += 16 * i_dst_stride;
        src += 16 * i_src_stride;
    }
}

void x264_frame_delete_list( x264_frame_t **list )
{
    int i = 0;
    if( !list )
        return;
    while( list[i] )
        x264_frame_delete( list[i++] );
    x264_free( list );
}

int x264_sync_frame_list_init( x264_sync_frame_list_t *slist, int max_size )
{
    if( max_size < 0 )
        return -1;
    slist->i_max_size = max_size;
    slist->i_size = 0;
    CHECKED_MALLOCZERO( slist->list, (max_size+1) * sizeof(x264_frame_t*) );
    if( x264_pthread_mutex_init( &slist->mutex, NULL ) ||
        x264_pthread_cond_init( &slist->cv_fill, NULL ) ||
        x264_pthread_cond_init( &slist->cv_empty, NULL ) )
        return -1;
    return 0;
fail:
    return -1;
}

void x264_sync_frame_list_delete( x264_sync_frame_list_t *slist )
{
    x264_pthread_mutex_destroy( &slist->mutex );
    x264_pthread_cond_destroy( &slist->cv_fill );
    x264_pthread_cond_destroy( &slist->cv_empty );
    x264_frame_delete_list( slist->list );
}

void x264_sync_frame_list_push( x264_sync_frame_list_t *slist, x264_frame_t *frame )
{
    x264_pthread_mutex_lock( &slist->mutex );
    while( slist->i_size == slist->i_max_size )
        x264_pthread_cond_wait( &slist->cv_empty, &slist->mutex );
    slist->list[ slist->i_size++ ] = frame;
    x264_pthread_mutex_unlock( &slist->mutex );
    x264_pthread_cond_broadcast( &slist->cv_fill );
}

x264_frame_t *x264_sync_frame_list_pop( x264_sync_frame_list_t *slist )
{
    x264_frame_t *frame;
    x264_pthread_mutex_lock( &slist->mutex );
    while( !slist->i_size )
        x264_pthread_cond_wait( &slist->cv_fill, &slist->mutex );
    frame = slist->list[ --slist->i_size ];
    slist->list[ slist->i_size ] = NULL;
    x264_pthread_cond_broadcast( &slist->cv_empty );
    x264_pthread_mutex_unlock( &slist->mutex );
    return frame;
}
int xClip( int iValue, int imin, int imax )
{
  ROTRS( iValue < imin, imin );
  ROTRS( iValue > imax, imax );
  return iValue;
}
int CeilLog2( int i )
{
  int s = 0; i--;
  while( i > 0 )
  {
    s++;
    i >>= 1;
  }
  return s;
}

void writeCsp(pixel* p, FILE* file, int width, int height,int stride)
{
  int i = 0; 
  for( ;i < height; i++ )
  {
    pixel* buffer = p + i * stride;
    int            wsize  = (int)fwrite( buffer, sizeof(pixel), width, file );
    if( wsize != width )
    {
      fprintf(stderr, "\nERROR: while writing to output file!\n\n");
    }
  }
}




int readColorComponent(pixel *p,pixel *src,int width,int height,int stride,int lines,int src_stride){
  int iMaxPadWidth  = gMin( stride, ( ( width  + 15 ) >> 4 ) << 4 );
  int iMaxPadHeight = gMin( lines,  ( ( height + 31 ) >> 5 ) << 5 );
  int i = 0;
  for( ; i < height; i++ )
  {
    unsigned char* buffer = p + i * stride;
    //int            rsize  = (int)fread( buffer, sizeof(unsigned char), width, file );
    memcpy(buffer,src,width*sizeof(unsigned char));
    src += src_stride;
    int xp = width; 
    for( ;xp < iMaxPadWidth; xp++ )
    {
      buffer[xp] = buffer[xp-1];
    }
  }
  int yp = height;
  for( ; yp < iMaxPadHeight; yp++ )
  {
    unsigned char* buffer  = p + yp * stride;
    unsigned char* bufferX = buffer - stride;
    int xp = 0;
    for( ; xp < stride; xp++ )
    {
      buffer[xp] = bufferX[xp];
    }
  }
  return 0;   
}

void xCopyToImageBuffer( unsigned char* pucSrc, int iWidth, int iHeight, int iStride,DownConvert* cDownConvert )
{
  int* piDes = cDownConvert->m_paiImageBuffer;
  int j = 0;
  for( ; j < iHeight; j++ )
  {
    int i = 0;
    for( ; i < iWidth;  i++ )
    {
      piDes[i] = (int)pucSrc[i];
    }
    piDes   += cDownConvert->m_iImageStride;
    pucSrc  += iStride;
  }
}

void xCopyFromImageBuffer( unsigned char* pucDes, int iWidth, int iHeight, int iStride,DownConvert* cDownConvert  )
{
  int* piSrc = cDownConvert->m_paiImageBuffer;
  int j = 0;
  for( ; j < iHeight; j++ )
  {
    int i = 0;
    for(; i < iWidth;  i++ )
    {
      pucDes[i] = (unsigned char)piSrc[i];
    }
    pucDes  += iStride;
    piSrc   += cDownConvert->m_iImageStride;
  }
}

void xBasicIntraUpsampling( int  iBaseW,   int  iBaseH,   int  iCurrW,   int  iCurrH,
                                    int  iLOffset, int  iTOffset, int  iROffset, int  iBOffset,
                                    int  iShiftX,  int  iShiftY,  int  iScaleX,  int  iScaleY,
                                    int  iOffsetX, int  iOffsetY, int  iAddX,    int  iAddY,
                                    int  iDeltaX,  int  iDeltaY,  int  iYBorder, bool bChromaFilter, int iMargin,DownConvert* cDownConvert)
{

  int filter16_luma[16][4] =
  {
    {  0, 32,  0,  0 },
    { -1, 32,  2, -1 },
    { -2, 31,  4, -1 },
    { -3, 30,  6, -1 },
    { -3, 28,  8, -1 },
    { -4, 26, 11, -1 },
    { -4, 24, 14, -2 },
    { -3, 22, 16, -3 },
    { -3, 19, 19, -3 },
    { -3, 16, 22, -3 },
    { -2, 14, 24, -4 },
    { -1, 11, 26, -4 },
    { -1,  8, 28, -3 },
    { -1,  6, 30, -3 },
    { -1,  4, 31, -2 },
    { -1,  2, 32, -1 }
  };
  int filter16_chroma[16][2] =
  {
    { 32,  0 },
    { 30,  2 },
    { 28,  4 },
    { 26,  6 },
    { 24,  8 },
    { 22, 10 },
    { 20, 12 },
    { 18, 14 },
    { 16, 16 },
    { 14, 18 },
    { 12, 20 },
    { 10, 22 },
    {  8, 24 },
    {  6, 26 },
    {  4, 28 },
    {  2, 30 }
  };

  int iShiftXM4 = iShiftX - 4;
  int iShiftYM4 = iShiftY - 4;

  //========== horizontal upsampling ===========
  {
    int j = 0;
    for( ; j < iBaseH + 2 * iMargin; j++ )
    {
      int* piSrc = &(cDownConvert->m_paiImageBuffer[j*cDownConvert->m_iImageStride]);
      int i = 0;
      for( ; i < iCurrW; i++ )
      {
        if( i < iLOffset || i >= iCurrW - iROffset )
        {
          cDownConvert->m_paiTmp1dBuffer[i] = 128; // set to mid gray
          continue;
        }

        cDownConvert->m_paiTmp1dBuffer[i] = 0;

        int iRefPos16 = (int)( (unsigned int)( ( i - iOffsetX ) * iScaleX + iAddX ) >> iShiftXM4 ) - iDeltaX;
        int iPhase    = iRefPos16 & 15;
        int iRefPos   = iRefPos16 >> 4;

        if( bChromaFilter )
        {
          int k = 0;
          for( ; k < 2; k++ )
          {
            int m = xClip( iRefPos + k, -iMargin, iBaseW - 1 + iMargin ) + iMargin;
            cDownConvert->m_paiTmp1dBuffer[i] += filter16_chroma[iPhase][k] * piSrc[m];
          }
        }
        else
        {
          int k = 0;
          for( ; k < 4; k++ )
          {
            int m = xClip( iRefPos + k - 1, -iMargin, iBaseW - 1 + iMargin ) + iMargin;
            cDownConvert->m_paiTmp1dBuffer[i] += filter16_luma[iPhase][k] * piSrc[m];
          }
        }
      }
      //----- copy row back to image buffer -----
      memcpy( piSrc, cDownConvert->m_paiTmp1dBuffer, iCurrW*sizeof(int) );
    }
  }

  //========== vertical upsampling ===========
  {
    int i = 0;
    for( ; i < iCurrW; i++ )
    {
      int* piSrc = &cDownConvert->m_paiImageBuffer[i];
      int j = -iYBorder;
      for( ; j < iCurrH+iYBorder; j++ )
      {
        if( i < iLOffset            || i >= iCurrW - iROffset           ||
            j < iTOffset - iYBorder || j >= iCurrH - iBOffset + iYBorder  )
        {
          cDownConvert->m_paiTmp1dBuffer[j+iYBorder] = 128; // set to mid gray
          continue;
        }

        cDownConvert->m_paiTmp1dBuffer[j+iYBorder] = 0;

        int iPreShift   = ( j - iOffsetY ) * iScaleY + iAddY;
        int iPostShift  = ( j >= iOffsetY ? (int)( (unsigned int)iPreShift >> iShiftYM4 ) : ( iPreShift >> iShiftYM4 ) );
        int iRefPos16   = iPostShift - iDeltaY;
        int iPhase      = iRefPos16 & 15;
        int iRefPos     = iRefPos16 >> 4;

        if( bChromaFilter )
        {
          int k = 0;
          for( ; k < 2; k++ )
          {
            int m = xClip( iRefPos + k, -iMargin, iBaseH - 1 + iMargin ) + iMargin;
            cDownConvert->m_paiTmp1dBuffer[j+iYBorder] += filter16_chroma[iPhase][k] * piSrc[m*cDownConvert->m_iImageStride];
          }
        }
        else
        {
          int k = 0;
          for( ; k < 4; k++ )
          {
            int m = xClip( iRefPos + k - 1, -iMargin, iBaseH - 1 + iMargin ) + iMargin;
            cDownConvert->m_paiTmp1dBuffer[j+iYBorder] += filter16_luma[iPhase][k] * piSrc[m*cDownConvert->m_iImageStride];
          }
        }
        cDownConvert->m_paiTmp1dBuffer[j+iYBorder] = ( cDownConvert->m_paiTmp1dBuffer[j+iYBorder] + 512 ) >> 10;
      }
      //----- clip and copy back to image buffer -----
      int n = 0;
      for( ; n < iCurrH+2*iYBorder; n++ )
      {
        piSrc[n*cDownConvert->m_iImageStride] = xClip( cDownConvert->m_paiTmp1dBuffer[n], 0, 255 );
      }
    }
  }
}

void xCompIntraUpsampling( ResizeParameters* pcParameters, bool bChroma, bool bBotFlag, bool bVerticalInterpolation, bool bFrameMb, int iMargin ,DownConvert* cDownConvert)
{
  //===== set general parameters =====
  int   iBotField   = ( bBotFlag ? 1 : 0 );
  int   iFactor     = ( !bChroma ? 1 : 2 );
  int   iRefPhaseX  = ( !bChroma ? 0 : pcParameters->m_iRefLayerChromaPhaseX );
  int   iRefPhaseY  = ( !bChroma ? 0 : pcParameters->m_iRefLayerChromaPhaseY );
  int   iPhaseX     = ( !bChroma ? 0 : pcParameters->m_iChromaPhaseX );
  int   iPhaseY     = ( !bChroma ? 0 : pcParameters->m_iChromaPhaseY );
  int   iRefW       = pcParameters->m_iRefLayerFrmWidth   / iFactor;  // reference layer frame width
  int   iRefH       = pcParameters->m_iRefLayerFrmHeight  / iFactor;  // reference layer frame height
  int   iOutW       = pcParameters->m_iScaledRefFrmWidth  / iFactor;  // scaled reference layer frame width
  int   iOutH       = pcParameters->m_iScaledRefFrmHeight / iFactor;  // scaled reference layer frame height
  int   iGlobalW    = pcParameters->m_iFrameWidth         / iFactor;  // current frame width
  int   iGlobalH    = pcParameters->m_iFrameHeight        / iFactor;  // current frame height
  int   iLeftOffset = pcParameters->m_iLeftFrmOffset      / iFactor;  // current left frame offset
  int   iTopOffset  = pcParameters->m_iTopFrmOffset       / iFactor;  // current top  frame offset

  //===== set input/output size =====
  int   iBaseField  = ( pcParameters->m_bRefLayerFrameMbsOnlyFlag ? 0 : 1 );
  int   iCurrField  = ( pcParameters->m_bRefLayerFrameMbsOnlyFlag && pcParameters->m_bFrameMbsOnlyFlag ? 0 : 1 );
  int   iBaseW      = iRefW;
  int   iBaseH      = iRefH      >> iBaseField;
  int   iCurrW      = iGlobalW;
  int   iCurrH      = iGlobalH   >> iCurrField;
  int   iLOffset    = iLeftOffset;
  int   iTOffset    = iTopOffset >> iCurrField;
  int   iROffset    = iCurrW - iLOffset -   iOutW;
  int   iBOffset    = iCurrH - iTOffset - ( iOutH >> iCurrField );
  int   iYBorder    = ( bVerticalInterpolation ? ( bChroma ? 1 : 2 ) : 0 );

  //===== set position calculation parameters =====
  int   iScaledW    = iOutW;
  int   iScaledH    = ( ! pcParameters->m_bRefLayerFrameMbsOnlyFlag || pcParameters->m_bFrameMbsOnlyFlag ? iOutH : iOutH / 2 );
  int   iShiftX     = ( pcParameters->m_iLevelIdc <= 30 ? 16 : 31 - CeilLog2( iRefW ) );
  int   iShiftY     = ( pcParameters->m_iLevelIdc <= 30 ? 16 : 31 - CeilLog2( iRefH ) );
  int   iScaleX     = ( ( (unsigned int)iRefW << iShiftX ) + ( iScaledW >> 1 ) ) / iScaledW;
  int   iScaleY     = ( ( (unsigned int)iRefH << iShiftY ) + ( iScaledH >> 1 ) ) / iScaledH;
  if( ! pcParameters->m_bFrameMbsOnlyFlag || ! pcParameters->m_bRefLayerFrameMbsOnlyFlag )
  {
    if( pcParameters->m_bRefLayerFrameMbsOnlyFlag )
    {
      iPhaseY       = iPhaseY + 4 * iBotField + ( 3 - iFactor );
      iRefPhaseY    = 2 * iRefPhaseY + 2;
    }
    else
    {
      iPhaseY       = iPhaseY    + 4 * iBotField;
      iRefPhaseY    = iRefPhaseY + 4 * iBotField;
    }
  }
  int   iOffsetX    = iLeftOffset;
  int   iOffsetY    = iTopOffset;
  int   iAddX       = ( ( ( iRefW * ( 2 + iPhaseX ) ) << ( iShiftX - 2 ) ) + ( iScaledW >> 1 ) ) / iScaledW + ( 1 << ( iShiftX - 5 ) );
  int   iAddY       = ( ( ( iRefH * ( 2 + iPhaseY ) ) << ( iShiftY - 2 ) ) + ( iScaledH >> 1 ) ) / iScaledH + ( 1 << ( iShiftY - 5 ) );
  int   iDeltaX     = 4 * ( 2 + iRefPhaseX );
  int   iDeltaY     = 4 * ( 2 + iRefPhaseY );
  if( ! pcParameters->m_bFrameMbsOnlyFlag || ! pcParameters->m_bRefLayerFrameMbsOnlyFlag )
  {
    iOffsetY        = iTopOffset / 2;
    iAddY           = ( ( ( iRefH * ( 2 + iPhaseY ) ) << ( iShiftY - 3 ) ) + ( iScaledH >> 1 ) ) / iScaledH + ( 1 << ( iShiftY - 5 ) );
    iDeltaY         = 2 * ( 2 + iRefPhaseY );
  }

  //===== basic interpolation of a frame or a field =====
  xBasicIntraUpsampling ( iBaseW,   iBaseH,   iCurrW,   iCurrH,
                          iLOffset, iTOffset, iROffset, iBOffset,
                          iShiftX,  iShiftY,  iScaleX,  iScaleY,
                          iOffsetX, iOffsetY, iAddX,    iAddY,
                          iDeltaX,  iDeltaY,  iYBorder, bChroma, iMargin,cDownConvert);

  //===== vertical interpolation for second field =====
  //if( bVerticalInterpolation )
//  {
//    xVertIntraUpsampling( iCurrW,   iCurrH,
//                          iLOffset, iTOffset, iROffset, iBOffset,
//                          iYBorder, bBotFlag, bFrameMb, bChroma );
//  }
}


void upsamplingSVC( pixel* pucBufferY,int iStrideY, ResizeParameters* pcParameters,int bBotCoincided,DownConvert* cDownConvert )
{
  int   iBaseW                  =   pcParameters->m_iRefLayerFrmWidth;
  int   iBaseH                  =   pcParameters->m_iRefLayerFrmHeight;
  int   iCurrW                  =   pcParameters->m_iFrameWidth;
  int   iCurrH                  =   pcParameters->m_iFrameHeight;
  bool  bTopAndBottomResampling = ( pcParameters->m_bRefLayerFrameMbsOnlyFlag == false  &&
                                    pcParameters->m_bRefLayerFieldPicFlag     == false  &&
                                    pcParameters->m_bFrameMbsOnlyFlag         == false  &&
                                    pcParameters->m_bFieldPicFlag             == false    );
  bool  bFrameBasedResampling   = ( pcParameters->m_bFrameMbsOnlyFlag         == true   &&
                                    pcParameters->m_bRefLayerFrameMbsOnlyFlag == true     );
  bool  bBotFieldFrameMbsOnly   = ( pcParameters->m_bRefLayerFrameMbsOnlyFlag == true   &&
                                    pcParameters->m_bFieldPicFlag             == true   &&
                                    pcParameters->m_bBotFieldFlag             == true     );
  bool  bVerticalInterpolation  = ( bBotFieldFrameMbsOnly                     == true   ||
                                   (bFrameBasedResampling                     == false  &&
                                    pcParameters->m_bFieldPicFlag             == false   ));
  bool  bFrameMb                = ( bBotFieldFrameMbsOnly                     == false    );
  bool  bCurrBotField           = ( pcParameters->m_bFieldPicFlag             == true   &&
                                    pcParameters->m_bBotFieldFlag             == true     );
  bool  bBotFieldFlag           = ( pcParameters->m_bRefLayerFrameMbsOnlyFlag ?  false
                                  : pcParameters->m_bFieldPicFlag             ?  pcParameters->m_bBotFieldFlag
                                  : pcParameters->m_bRefLayerFieldPicFlag     ?  pcParameters->m_bRefLayerBotFieldFlag
                                  : false );
  int   iBaseField              = ( pcParameters->m_bRefLayerFrameMbsOnlyFlag ?  0 : 1 );
  int   iCurrField              = ( pcParameters->m_bFieldPicFlag             ?  1 : 0 );
  int   iBaseBot                = ( bBotFieldFlag ? 1 : 0 );
  int   iCurrBot                = ( bCurrBotField ? 1 : 0 );

  //==== check bot field coincided parameter for progressive to interlaced resampling =====
  if( pcParameters->m_bRefLayerFrameMbsOnlyFlag && ! pcParameters->m_bFrameMbsOnlyFlag )
  {
    bBotFieldFlag = bBotCoincided;
  }

  //=======================
  //=====   L U M A   =====
  //=======================
    unsigned char* pSrc = pucBufferY + iStrideY * iBaseBot;
    unsigned char* pDes = pucBufferY + iStrideY * iCurrBot;
    xCopyToImageBuffer  ( pSrc,         iBaseW, iBaseH >> iBaseField, iStrideY << iBaseField,cDownConvert);
    xCompIntraUpsampling( pcParameters, false,  bBotFieldFlag,bVerticalInterpolation, bFrameMb,0,cDownConvert);
    xCopyFromImageBuffer( pDes,         iCurrW, iCurrH >> iCurrField, iStrideY << iCurrField ,cDownConvert);
}



//传进参数，DownConvert*   
void xBasicDownsampling(int iBaseW,   int iBaseH,   int iCurrW,   int iCurrH,
                                 int iLOffset, int iTOffset, int iROffset, int iBOffset,
                                 int iShiftX,  int iShiftY,  int iScaleX,  int iScaleY,
                                 int iAddX,    int iAddY,    int iDeltaX,  int iDeltaY,DownConvert* cDownConvert )
{
  const int filter16[8][16][12] =
  {
    { // D = 1
      {   0,   0,   0,   0,   0, 128,   0,   0,   0,   0,   0,   0 },
      {   0,   0,   0,   2,  -6, 127,   7,  -2,   0,   0,   0,   0 },
      {   0,   0,   0,   3, -12, 125,  16,  -5,   1,   0,   0,   0 },
      {   0,   0,   0,   4, -16, 120,  26,  -7,   1,   0,   0,   0 },
      {   0,   0,   0,   5, -18, 114,  36, -10,   1,   0,   0,   0 },
      {   0,   0,   0,   5, -20, 107,  46, -12,   2,   0,   0,   0 },
      {   0,   0,   0,   5, -21,  99,  57, -15,   3,   0,   0,   0 },
      {   0,   0,   0,   5, -20,  89,  68, -18,   4,   0,   0,   0 },
      {   0,   0,   0,   4, -19,  79,  79, -19,   4,   0,   0,   0 },
      {   0,   0,   0,   4, -18,  68,  89, -20,   5,   0,   0,   0 },
      {   0,   0,   0,   3, -15,  57,  99, -21,   5,   0,   0,   0 },
      {   0,   0,   0,   2, -12,  46, 107, -20,   5,   0,   0,   0 },
      {   0,   0,   0,   1, -10,  36, 114, -18,   5,   0,   0,   0 },
      {   0,   0,   0,   1,  -7,  26, 120, -16,   4,   0,   0,   0 },
      {   0,   0,   0,   1,  -5,  16, 125, -12,   3,   0,   0,   0 },
      {   0,   0,   0,   0,  -2,   7, 127,  -6,   2,   0,   0,   0 }
    },
    { // D = 1.5
      {   0,   2,   0, -14,  33,  86,  33, -14,   0,   2,   0,   0 },
      {   0,   1,   1, -14,  29,  85,  38, -13,  -1,   2,   0,   0 },
      {   0,   1,   2, -14,  24,  84,  43, -12,  -2,   2,   0,   0 },
      {   0,   1,   2, -13,  19,  83,  48, -11,  -3,   2,   0,   0 },
      {   0,   0,   3, -13,  15,  81,  53, -10,  -4,   3,   0,   0 },
      {   0,   0,   3, -12,  11,  79,  57,  -8,  -5,   3,   0,   0 },
      {   0,   0,   3, -11,   7,  76,  62,  -5,  -7,   3,   0,   0 },
      {   0,   0,   3, -10,   3,  73,  65,  -2,  -7,   3,   0,   0 },
      {   0,   0,   3,  -9,   0,  70,  70,   0,  -9,   3,   0,   0 },
      {   0,   0,   3,  -7,  -2,  65,  73,   3, -10,   3,   0,   0 },
      {   0,   0,   3,  -7,  -5,  62,  76,   7, -11,   3,   0,   0 },
      {   0,   0,   3,  -5,  -8,  57,  79,  11, -12,   3,   0,   0 },
      {   0,   0,   3,  -4, -10,  53,  81,  15, -13,   3,   0,   0 },
      {   0,   0,   2,  -3, -11,  48,  83,  19, -13,   2,   1,   0 },
      {   0,   0,   2,  -2, -12,  43,  84,  24, -14,   2,   1,   0 },
      {   0,   0,   2,  -1, -13,  38,  85,  29, -14,   1,   1,   0 }
    },
    { // D = 2
      {   2,   0, -10,   0,  40,  64,  40,   0, -10,   0,   2,   0 },
      {   2,   1,  -9,  -2,  37,  64,  42,   2, -10,  -1,   2,   0 },
      {   2,   1,  -9,  -3,  34,  64,  44,   4, -10,  -1,   2,   0 },
      {   2,   1,  -8,  -5,  31,  63,  47,   6, -10,  -2,   3,   0 },
      {   1,   2,  -8,  -6,  29,  62,  49,   8, -10,  -2,   3,   0 },
      {   1,   2,  -7,  -7,  26,  61,  52,  10, -10,  -3,   3,   0 },
      {   1,   2,  -6,  -8,  23,  60,  54,  13, -10,  -4,   3,   0 },
      {   1,   2,  -6,  -9,  20,  59,  56,  15, -10,  -4,   3,   1 },
      {   1,   2,  -5,  -9,  18,  57,  57,  18,  -9,  -5,   2,   1 },
      {   1,   3,  -4, -10,  15,  56,  59,  20,  -9,  -6,   2,   1 },
      {   0,   3,  -4, -10,  13,  54,  60,  23,  -8,  -6,   2,   1 },
      {   0,   3,  -3, -10,  10,  52,  61,  26,  -7,  -7,   2,   1 },
      {   0,   3,  -2, -10,   8,  49,  62,  29,  -6,  -8,   2,   1 },
      {   0,   3,  -2, -10,   6,  47,  63,  31,  -5,  -8,   1,   2 },
      {   0,   2,  -1, -10,   4,  44,  64,  34,  -3,  -9,   1,   2 },
      {   0,   2,  -1, -10,   2,  42,  64,  37,  -2,  -9,   1,   2 }
    },
    { // D = 2.5
      {   0,  -4,  -7,  11,  38,  52,  38,  11,  -7,  -4,   0,   0 },
      {   0,  -4,  -7,   9,  37,  51,  40,  13,  -6,  -7,   2,   0 },
      {   0,  -3,  -7,   8,  35,  51,  41,  14,  -5,  -7,   1,   0 },
      {   0,  -2,  -8,   6,  33,  51,  42,  16,  -5,  -7,   2,   0 },
      {   0,  -2,  -8,   5,  32,  50,  43,  18,  -4,  -8,   2,   0 },
      {   0,  -2,  -8,   4,  30,  50,  45,  19,  -3,  -8,   1,   0 },
      {   0,  -1,  -8,   2,  28,  49,  46,  21,  -2,  -8,   1,   0 },
      {   0,  -1,  -8,   1,  26,  49,  47,  23,  -1,  -8,   0,   0 },
      {   0,   0,  -8,   0,  24,  48,  48,  24,   0,  -8,   0,   0 },
      {   0,   0,  -8,  -1,  23,  47,  49,  26,   1,  -8,  -1,   0 },
      {   0,   1,  -8,  -2,  21,  46,  49,  28,   2,  -8,  -1,   0 },
      {   0,   1,  -8,  -3,  19,  45,  50,  30,   4,  -8,  -2,   0 },
      {   0,   2,  -8,  -4,  18,  43,  50,  32,   5,  -8,  -2,   0 },
      {   0,   2,  -7,  -5,  16,  42,  51,  33,   6,  -8,  -2,   0 },
      {   0,   1,  -7,  -5,  14,  41,  51,  35,   8,  -7,  -3,   0 },
      {   0,   2,  -7,  -6,  13,  40,  51,  37,   9,  -7,  -4,   0 }
    },
    { // D = 3
      {  -2,  -7,   0,  17,  35,  43,  35,  17,   0,  -7,  -5,   2 },
      {  -2,  -7,  -1,  16,  34,  43,  36,  18,   1,  -7,  -5,   2 },
      {  -1,  -7,  -1,  14,  33,  43,  36,  19,   1,  -6,  -5,   2 },
      {  -1,  -7,  -2,  13,  32,  42,  37,  20,   3,  -6,  -5,   2 },
      {   0,  -7,  -3,  12,  31,  42,  38,  21,   3,  -6,  -5,   2 },
      {   0,  -7,  -3,  11,  30,  42,  39,  23,   4,  -6,  -6,   1 },
      {   0,  -7,  -4,  10,  29,  42,  40,  24,   5,  -6,  -6,   1 },
      {   1,  -7,  -4,   9,  27,  41,  40,  25,   6,  -5,  -6,   1 },
      {   1,  -6,  -5,   7,  26,  41,  41,  26,   7,  -5,  -6,   1 },
      {   1,  -6,  -5,   6,  25,  40,  41,  27,   9,  -4,  -7,   1 },
      {   1,  -6,  -6,   5,  24,  40,  42,  29,  10,  -4,  -7,   0 },
      {   1,  -6,  -6,   4,  23,  39,  42,  30,  11,  -3,  -7,   0 },
      {   2,  -5,  -6,   3,  21,  38,  42,  31,  12,  -3,  -7,   0 },
      {   2,  -5,  -6,   3,  20,  37,  42,  32,  13,  -2,  -7,  -1 },
      {   2,  -5,  -6,   1,  19,  36,  43,  33,  14,  -1,  -7,  -1 },
      {   2,  -5,  -7,   1,  18,  36,  43,  34,  16,  -1,  -7,  -2 }
    },
    { // D = 3.5
      {  -6,  -3,   5,  19,  31,  36,  31,  19,   5,  -3,  -6,   0 },
      {  -6,  -4,   4,  18,  31,  37,  32,  20,   6,  -3,  -6,  -1 },
      {  -6,  -4,   4,  17,  30,  36,  33,  21,   7,  -3,  -6,  -1 },
      {  -5,  -5,   3,  16,  30,  36,  33,  22,   8,  -2,  -6,  -2 },
      {  -5,  -5,   2,  15,  29,  36,  34,  23,   9,  -2,  -6,  -2 },
      {  -5,  -5,   2,  15,  28,  36,  34,  24,  10,  -2,  -6,  -3 },
      {  -4,  -5,   1,  14,  27,  36,  35,  24,  10,  -1,  -6,  -3 },
      {  -4,  -5,   0,  13,  26,  35,  35,  25,  11,   0,  -5,  -3 },
      {  -4,  -6,   0,  12,  26,  36,  36,  26,  12,   0,  -6,  -4 },
      {  -3,  -5,   0,  11,  25,  35,  35,  26,  13,   0,  -5,  -4 },
      {  -3,  -6,  -1,  10,  24,  35,  36,  27,  14,   1,  -5,  -4 },
      {  -3,  -6,  -2,  10,  24,  34,  36,  28,  15,   2,  -5,  -5 },
      {  -2,  -6,  -2,   9,  23,  34,  36,  29,  15,   2,  -5,  -5 },
      {  -2,  -6,  -2,   8,  22,  33,  36,  30,  16,   3,  -5,  -5 },
      {  -1,  -6,  -3,   7,  21,  33,  36,  30,  17,   4,  -4,  -6 },
      {  -1,  -6,  -3,   6,  20,  32,  37,  31,  18,   4,  -4,  -6 }
    },
    { // D = 4
      {  -9,   0,   9,  20,  28,  32,  28,  20,   9,   0,  -9,   0 },
      {  -9,   0,   8,  19,  28,  32,  29,  20,  10,   0,  -4,  -5 },
      {  -9,  -1,   8,  18,  28,  32,  29,  21,  10,   1,  -4,  -5 },
      {  -9,  -1,   7,  18,  27,  32,  30,  22,  11,   1,  -4,  -6 },
      {  -8,  -2,   6,  17,  27,  32,  30,  22,  12,   2,  -4,  -6 },
      {  -8,  -2,   6,  16,  26,  32,  31,  23,  12,   2,  -4,  -6 },
      {  -8,  -2,   5,  16,  26,  31,  31,  23,  13,   3,  -3,  -7 },
      {  -8,  -3,   5,  15,  25,  31,  31,  24,  14,   4,  -3,  -7 },
      {  -7,  -3,   4,  14,  25,  31,  31,  25,  14,   4,  -3,  -7 },
      {  -7,  -3,   4,  14,  24,  31,  31,  25,  15,   5,  -3,  -8 },
      {  -7,  -3,   3,  13,  23,  31,  31,  26,  16,   5,  -2,  -8 },
      {  -6,  -4,   2,  12,  23,  31,  32,  26,  16,   6,  -2,  -8 },
      {  -6,  -4,   2,  12,  22,  30,  32,  27,  17,   6,  -2,  -8 },
      {  -6,  -4,   1,  11,  22,  30,  32,  27,  18,   7,  -1,  -9 },
      {  -5,  -4,   1,  10,  21,  29,  32,  28,  18,   8,  -1,  -9 },
      {  -5,  -4,   0,  10,  20,  29,  32,  28,  19,   8,   0,  -9 }
    },
    { // D = 5.5
      {  -8,   7,  13,  18,  22,  24,  22,  18,  13,   7,   2, -10 },
      {  -8,   7,  13,  18,  22,  23,  22,  19,  13,   7,   2, -10 },
      {  -8,   6,  12,  18,  22,  23,  22,  19,  14,   8,   2, -10 },
      {  -9,   6,  12,  17,  22,  23,  23,  19,  14,   8,   3, -10 },
      {  -9,   6,  12,  17,  21,  23,  23,  19,  14,   9,   3, -10 },
      {  -9,   5,  11,  17,  21,  23,  23,  20,  15,   9,   3, -10 },
      {  -9,   5,  11,  16,  21,  23,  23,  20,  15,   9,   4, -10 },
      {  -9,   5,  10,  16,  21,  23,  23,  20,  15,  10,   4, -10 },
      { -10,   5,  10,  16,  20,  23,  23,  20,  16,  10,   5, -10 },
      { -10,   4,  10,  15,  20,  23,  23,  21,  16,  10,   5,  -9 },
      { -10,   4,   9,  15,  20,  23,  23,  21,  16,  11,   5,  -9 },
      { -10,   3,   9,  15,  20,  23,  23,  21,  17,  11,   5,  -9 },
      { -10,   3,   9,  14,  19,  23,  23,  21,  17,  12,   6,  -9 },
      { -10,   3,   8,  14,  19,  23,  23,  22,  17,  12,   6,  -9 },
      { -10,   2,   8,  14,  19,  22,  23,  22,  18,  12,   6,  -8 },
      { -10,   2,   7,  13,  19,  22,  23,  22,  18,  13,   7,  -8 }
    }
  };

  //===== determine filter sets =====
  int iCropW      = iCurrW - iLOffset - iROffset;
  int iCropH      = iCurrH - iTOffset - iBOffset;
  int iVerFilter  = 0;
  int iHorFilter  = 0;
  if      (  4 * iCropH > 15 * iBaseH )   iVerFilter  = 7;
  else if (  7 * iCropH > 20 * iBaseH )   iVerFilter  = 6;
  else if (  2 * iCropH >  5 * iBaseH )   iVerFilter  = 5;
  else if (  1 * iCropH >  2 * iBaseH )   iVerFilter  = 4;
  else if (  3 * iCropH >  5 * iBaseH )   iVerFilter  = 3;//
  else if (  4 * iCropH >  5 * iBaseH )   iVerFilter  = 2;
  else if ( 19 * iCropH > 20 * iBaseH )   iVerFilter  = 1;
  if      (  4 * iCropW > 15 * iBaseW )   iHorFilter  = 7;
  else if (  7 * iCropW > 20 * iBaseW )   iHorFilter  = 6;
  else if (  2 * iCropW >  5 * iBaseW )   iHorFilter  = 5;
  else if (  1 * iCropW >  2 * iBaseW )   iHorFilter  = 4;
  else if (  3 * iCropW >  5 * iBaseW )   iHorFilter  = 3;//
  else if (  4 * iCropW >  5 * iBaseW )   iHorFilter  = 2;
  else if ( 19 * iCropW > 20 * iBaseW )   iHorFilter  = 1;

  int iShiftXM4 = iShiftX - 4;
  int iShiftYM4 = iShiftY - 4;

  //===== horizontal downsampling =====
  {
    int j = 0;
    for(; j < iCurrH; j++ )
    {
      int* piSrc = &(cDownConvert->m_paiImageBuffer[j*cDownConvert->m_iImageStride]);
      int i = 0;
      for( ; i < iBaseW; i++ )
      {
        int iRefPos16 = (int)( (unsigned int)( i * iScaleX + iAddX ) >> iShiftXM4 ) - iDeltaX;
        int iPhase    = iRefPos16  & 15;
        int iRefPos   = iRefPos16 >>  4;

        cDownConvert->m_paiTmp1dBuffer[i] = 0;
        int k = 0; 
        for( ;k < 12; k++ )
        {
          int m = xClip( iRefPos + k - 5, 0, iCurrW - 1 );
          cDownConvert->m_paiTmp1dBuffer[i] += filter16[iHorFilter][iPhase][k] * piSrc[m];
        }
      }
      //--- copy row back to image buffer ---
      memcpy( piSrc, cDownConvert->m_paiTmp1dBuffer, iBaseW*sizeof(int) );
    }
  }

  //===== vertical downsampling =====
  {
    int i = 0; 
    for( ;i < iBaseW; i++ )
    {
      int* piSrc = &cDownConvert->m_paiImageBuffer[i];
      int j = 0;
      for( ; j < iBaseH; j++ )
      {
        int iRefPos16 = (int)( (unsigned int)( j * iScaleY + iAddY ) >> iShiftYM4 ) - iDeltaY;
        int iPhase    = iRefPos16  & 15;
        int iRefPos   = iRefPos16 >>  4;

        cDownConvert->m_paiTmp1dBuffer[j] = 0;
        int k = 0;
        for( ;k < 12; k++ )
        {
          int m = xClip( iRefPos + k - 5, 0, iCurrH - 1 );
          cDownConvert->m_paiTmp1dBuffer[j] += filter16[iVerFilter][iPhase][k] * piSrc[m*cDownConvert->m_iImageStride];
        }
        cDownConvert->m_paiTmp1dBuffer[j] = ( cDownConvert->m_paiTmp1dBuffer[j] + 8192 ) >> 14;
      }
      //--- clip and copy back to image buffer ---
      int n = 0;
      for( ; n < iBaseH; n++ )
      {
        piSrc[n*cDownConvert->m_iImageStride] = xClip( cDownConvert->m_paiTmp1dBuffer[n], 0, 255 );
      }
    }
  }
}

void xCompDownsampling(pixel* p,ResizeParameters* pcParameters, int bChroma, int bBotFlag, int bVerticalDownsampling,DownConvert* cDownConvert)
{
  //===== set general parameters =====
  int   iBotField   = ( bBotFlag ? 1 : 0 );
  int   iFactor     = ( !bChroma ? 1 : 2 );
  int   iRefPhaseX  = ( !bChroma ? 0 : pcParameters->m_iChromaPhaseX );
  int   iRefPhaseY  = ( !bChroma ? 0 : pcParameters->m_iChromaPhaseY );
  int   iPhaseX     = ( !bChroma ? 0 : pcParameters->m_iRefLayerChromaPhaseX );
  int   iPhaseY     = ( !bChroma ? 0 : pcParameters->m_iRefLayerChromaPhaseY );
  int   iRefW       = pcParameters->m_iFrameWidth         / iFactor;  // reference layer frame width
  int   iRefH       = pcParameters->m_iFrameHeight        / iFactor;  // reference layer frame height
  int   iOutW       = pcParameters->m_iScaledRefFrmWidth  / iFactor;  // scaled reference layer frame width
  int   iOutH       = pcParameters->m_iScaledRefFrmHeight / iFactor;  // scaled reference layer frame height
  int   iGlobalW    = pcParameters->m_iRefLayerFrmWidth   / iFactor;  // current frame width
  int   iGlobalH    = pcParameters->m_iRefLayerFrmHeight  / iFactor;  // current frame height
  int   iLeftOffset = pcParameters->m_iLeftFrmOffset      / iFactor;  // current left frame offset
  int   iTopOffset  = pcParameters->m_iTopFrmOffset       / iFactor;  // current  top frame offset

  //===== set input/output size =====
  int   iBaseField  = ( pcParameters->m_bRefLayerFrameMbsOnlyFlag ? 0 : 1 );
  int   iCurrField  = ( pcParameters->m_bRefLayerFrameMbsOnlyFlag && pcParameters->m_bFrameMbsOnlyFlag ? 0 : 1 );
  int   iBaseW      = iRefW;
  int   iBaseH      = iRefH       >> iBaseField;
  int   iCurrW      = iGlobalW;
  int   iCurrH      = iGlobalH    >> iCurrField;
  int   iLOffset    = iLeftOffset;
  int   iTOffset    = iTopOffset  >> iCurrField;
  int   iROffset    = iCurrW - iLOffset -   iOutW;
  int   iBOffset    = iCurrH - iTOffset - ( iOutH >> iCurrField );

  //===== set position calculation parameters =====
  int   iScaledW    = iOutW;
  int   iScaledH    = ( ! pcParameters->m_bRefLayerFrameMbsOnlyFlag || pcParameters->m_bFrameMbsOnlyFlag ? iOutH : iOutH / 2 );
  int   iShiftX     = ( pcParameters->m_iLevelIdc <= 30 ? 16 : 31 - CeilLog2( iScaledW ) );
  int   iShiftY     = ( pcParameters->m_iLevelIdc <= 30 ? 16 : 31 - CeilLog2( iScaledH ) );
  int   iScaleX     = ( ( (unsigned int)iScaledW << iShiftX ) + ( iRefW >> 1 ) ) / iRefW;
  int   iScaleY     = ( ( (unsigned int)iScaledH << iShiftY ) + ( iRefH >> 1 ) ) / iRefH;
  if( ! pcParameters->m_bFrameMbsOnlyFlag || ! pcParameters->m_bRefLayerFrameMbsOnlyFlag )//jump
  {
    if( pcParameters->m_bRefLayerFrameMbsOnlyFlag )
    {
      iPhaseY       = iPhaseY + 4 * iBotField + ( 3 - iFactor );
      iRefPhaseY    = 2 * iRefPhaseY + 2;
    }
    else
    {
      iPhaseY       = iPhaseY    + 4 * iBotField;
      iRefPhaseY    = iRefPhaseY + 4 * iBotField;
    }
  }
  int   iAddX       = ( ( ( iScaledW * ( 2 + iRefPhaseX ) ) << ( iShiftX - 2 ) ) + ( iRefW >> 1 ) ) / iRefW + ( 1 << ( iShiftX - 5 ) );
  int   iAddY       = ( ( ( iScaledH * ( 2 + iRefPhaseY ) ) << ( iShiftY - 2 ) ) + ( iRefH >> 1 ) ) / iRefH + ( 1 << ( iShiftY - 5 ) );
  int   iDeltaX     = 4 * ( 2 + iPhaseX ) - ( iLeftOffset << 4 );
  int   iDeltaY     = 4 * ( 2 + iPhaseY ) - ( iTopOffset  << 4 );
  if( ! pcParameters->m_bFrameMbsOnlyFlag || ! pcParameters->m_bRefLayerFrameMbsOnlyFlag )//jump
  {
    iAddY           = ( ( ( iScaledH * ( 2 + iRefPhaseY ) ) << ( iShiftY - 3 ) ) + ( iRefH >> 1 ) ) / iRefH + ( 1 << ( iShiftY - 5 ) );
    iDeltaY         = 2 * ( 2 + iPhaseY ) - ( iTopOffset  << 3 );
  }

  //===== vertical downsampling to generate a field signal from a progressive frame =====
  if( bVerticalDownsampling )//jump
  {//函数没有调用，所以注释掉 
      ;
    //xVertDownsampling( iCurrW, iCurrH, bBotFlag );
  }
   
  //===== basic downsampling of a frame or field =====  [Linker error] undefined reference to `CeilLog2' 
  xBasicDownsampling( iBaseW,   iBaseH,   iCurrW,   iCurrH,
                      iLOffset, iTOffset, iROffset, iBOffset,
                      iShiftX,  iShiftY,  iScaleX,  iScaleY,
                      iAddX,    iAddY,    iDeltaX,  iDeltaY ,cDownConvert);
}

void downsamplingSVC( pixel* pucBufferY,int iStrideY, ResizeParameters* pcParameters,int bBotCoincided,DownConvert* cDownConvert ){
  int   iBaseW                  =   pcParameters->m_iFrameWidth;//176  
  int   iBaseH                  =   pcParameters->m_iFrameHeight;//128
  int   iCurrW                  =   pcParameters->m_iRefLayerFrmWidth;//352就是输入in:352x240
  int   iCurrH                  =   pcParameters->m_iRefLayerFrmHeight;//240
  int  bTopAndBottomResampling = ( pcParameters->m_bRefLayerFrameMbsOnlyFlag == false  &&
                                    pcParameters->m_bRefLayerFieldPicFlag     == false  &&
                                    pcParameters->m_bFrameMbsOnlyFlag         == false  &&
                                    pcParameters->m_bFieldPicFlag             == false    );//F
  int  bVerticalDownsampling   = ( pcParameters->m_bFrameMbsOnlyFlag         == true   &&
                                    pcParameters->m_bRefLayerFieldPicFlag     == true     );//F
  int  bCurrBotField           = ( pcParameters->m_bFieldPicFlag             == true   &&
                                    pcParameters->m_bBotFieldFlag             == true     );//F
  int  bBotFieldFlag           = ( pcParameters->m_bRefLayerFrameMbsOnlyFlag ?  false
                                  : pcParameters->m_bFieldPicFlag             ?  pcParameters->m_bBotFieldFlag
                                  : pcParameters->m_bRefLayerFieldPicFlag     ?  pcParameters->m_bRefLayerBotFieldFlag
                                  : false );//T
  int   iBaseField              = ( pcParameters->m_bRefLayerFrameMbsOnlyFlag ?  0 : 1 );//0
  int   iCurrField              = ( pcParameters->m_bFieldPicFlag             ?  1 : 0 );//0
  int   iBaseBot                = ( bBotFieldFlag ? 1 : 0 );
  int   iCurrBot                = ( bCurrBotField ? 1 : 0 );

  //==== check bot field coincided parameter for interlaced to progressive resampling =====
  if( pcParameters->m_bRefLayerFrameMbsOnlyFlag && ! pcParameters->m_bFrameMbsOnlyFlag )
  {
    bBotFieldFlag = bBotCoincided;
  }       
  //=======================
  //=====   L U M A   =====
  //=======================
  unsigned char* pSrc = pucBufferY + iStrideY * iCurrBot;
  unsigned char* pDes = pucBufferY + iStrideY * iBaseBot;
  xCopyToImageBuffer  ( pSrc,         iCurrW, iCurrH >> iCurrField, iStrideY << iCurrField ,cDownConvert);
  xCompDownsampling   ( pucBufferY,pcParameters, false,  bBotFieldFlag,        bVerticalDownsampling  ,cDownConvert);
  xCopyFromImageBuffer( pDes,         iBaseW, iBaseH >> iBaseField, iStrideY << iBaseField ,cDownConvert);                                               
}
void cDownConvertInit(int iMaxWidth, int iMaxHeight,int iMaxMargin,DownConvert* cDownConvert){
  iMaxWidth                += 2 * iMaxMargin;
  iMaxHeight               += 2 * iMaxMargin;
  int iPicSize              =   iMaxWidth * iMaxHeight;
  int iMaxDim               = ( iMaxWidth > iMaxHeight ? iMaxWidth : iMaxHeight );
  cDownConvert->m_iImageStride            = iMaxWidth;
  cDownConvert->m_paiImageBuffer          = malloc(sizeof(int)*iPicSize);
  cDownConvert->m_paiTmp1dBuffer          = malloc(sizeof(int)*iMaxDim);

  cDownConvert->m_padFilter               = malloc(sizeof(long)*TMM_TABLE_SIZE);
  cDownConvert->m_aiTmp1dBufferInHalfpel  = malloc(sizeof(int)*iMaxDim);
  cDownConvert->m_aiTmp1dBufferInQ1pel    = malloc(sizeof(int)*iMaxDim);
  cDownConvert->m_aiTmp1dBufferInQ3pel    = malloc(sizeof(int)*iMaxDim);
  cDownConvert->m_paiTmp1dBufferOut       = malloc(sizeof(int)*iMaxDim);                   
}
void resampleFrame( pixel*         p,
               DownConvert*       cDownConvert,
               ResizeParameters  rcRP,
               int                resamplingMethod,
               int                resamplingMode,
               int               resampling,
               int               upsampling,
               int               bSecondInputFrame ,int stride){
  if( upsampling ){
    ResizeParameters cRP = rcRP;
    {
      int iRefVerMbShift        = ( cRP.m_bRefLayerFrameMbsOnlyFlag ? 4 : 5 );
      int iScaledVerShift       = ( cRP.m_bFrameMbsOnlyFlag         ? 1 : 2 );
      int iHorDiv               = ( cRP.m_iRefLayerFrmWidth    <<               1 );
      int iVerDiv               = ( cRP.m_iRefLayerFrmHeight   << iScaledVerShift );
      int iRefFrmW              = ( ( cRP.m_iRefLayerFrmWidth   + ( 1 <<               4 ) - 1 ) >>               4 ) <<               4;  // round to next multiple of 16
      int iRefFrmH              = ( ( cRP.m_iRefLayerFrmHeight  + ( 1 <<  iRefVerMbShift ) - 1 ) >>  iRefVerMbShift ) <<  iRefVerMbShift;  // round to next multiple of 16 or 32 (for interlaced)
      int iScaledRefFrmW        = ( ( cRP.m_iScaledRefFrmWidth  * iRefFrmW + ( iHorDiv >> 1 ) ) / iHorDiv ) <<               1;  // scale and round to next multiple of  2
      int iScaledRefFrmH        = ( ( cRP.m_iScaledRefFrmHeight * iRefFrmH + ( iVerDiv >> 1 ) ) / iVerDiv ) << iScaledVerShift;  // scale and round to next multiple of  2 or  4 (for interlaced)
      cRP.m_iRefLayerFrmWidth   = iRefFrmW;
      cRP.m_iRefLayerFrmHeight  = iRefFrmH;
      cRP.m_iScaledRefFrmWidth  = iScaledRefFrmW;
      cRP.m_iScaledRefFrmHeight = iScaledRefFrmH;
    }
    
    //int stride = ( ( 352  + 15 ) >> 4 ) << 4;
    upsamplingSVC(p,stride,&cRP, resamplingMode == 3,cDownConvert);
    return;
}


//DownSample
  ResizeParameters cRP = rcRP;
  {
    int iRefVerMbShift        = ( cRP.m_bRefLayerFrameMbsOnlyFlag ? 4 : 5 );
    int iScaledVerShift       = ( cRP.m_bFrameMbsOnlyFlag         ? 1 : 2 );
    int iHorDiv               = ( cRP.m_iFrameWidth    <<               1 );
    int iVerDiv               = ( cRP.m_iFrameHeight   << iScaledVerShift );
    int iRefFrmW              = ( ( cRP.m_iFrameWidth   + ( 1 <<               4 ) - 1 ) >>               4 ) <<               4;        // round to next multiple of 16
    int iRefFrmH              = ( ( cRP.m_iFrameHeight  + ( 1 <<  iRefVerMbShift ) - 1 ) >>  iRefVerMbShift ) <<  iRefVerMbShift;        // round to next multiple of 16 or 32 (for interlaced)
    int iScaledRefFrmW        = ( ( cRP.m_iScaledRefFrmWidth  * iRefFrmW + ( iHorDiv >> 1 ) ) / iHorDiv ) <<               1;  // scale and round to next multiple of  2
    int iScaledRefFrmH        = ( ( cRP.m_iScaledRefFrmHeight * iRefFrmH + ( iVerDiv >> 1 ) ) / iVerDiv ) << iScaledVerShift;  // scale and round to next multiple of  2 or  4 (for interlaced)
    cRP.m_iFrameWidth         = iRefFrmW;
    cRP.m_iFrameHeight        = iRefFrmH;
    cRP.m_iScaledRefFrmWidth  = iScaledRefFrmW;
    cRP.m_iScaledRefFrmHeight = iScaledRefFrmH;
  }
  //int stride = ( ( 352/2  + 15 ) >> 4 ) << 4;//需要计算stride； 
  //downsamplingSVC( rcFrame.y.data,  rcFrame.y.stride, rcFrame.u.data,  rcFrame.u.stride, rcFrame.v.data,  rcFrame.v.stride, &cRP, resamplingMode == 3 );
  //downsamplingSVC(p,  rcFrame.y.stride,&cRP, resamplingMode == 3,cDownConvert);
  downsamplingSVC(p,  stride,&cRP, resamplingMode == 3,cDownConvert);
  return;                                                             
}

void x264_frame_expand_layers(x264_t *h,FILE *outputfile,int dst_stride,pixel *src,int src_stride,int win,int hin,int wout,int hout){
	printf ("-------------------x264_frame_expand_layers---------\n");
	pixel* p;
    int maxWidth = gMax(win,wout);
    int maxHeight = gMax(hin,hout);
    int minWidth = gMin(win,wout);
    int minHeight = gMin(hin,hout);
    int width = win;
    int height = hin;  
    int minWRnd16 = ( ( minWidth  + 15 ) >> 4 ) << 4;
    int minHRnd32 = ( ( minHeight + 31 ) >> 5 ) << 5;
    maxWidth      = ( ( maxWidth  * minWRnd16 + ( minWidth  << 4 ) - 1 ) / ( minWidth  << 4 ) ) << 4;
    maxHeight     = ( ( maxHeight * minHRnd32 + ( minHeight << 4 ) - 1 ) / ( minHeight << 4 ) ) << 4;
    int maxwidth  = ( ( maxWidth  + 15 ) >> 4 ) << 4;
    int maxheight = ( ( maxHeight + 15 ) >> 4 ) << 4;
    int size = maxwidth*maxheight;
    p = malloc(sizeof(pixel)*size);
    //CHECKED_MALLOC( p, sizeof(pixel)*size );
    int stride  = maxwidth;
    int lines   = maxheight;
    //readColorComponent(p,src,width,height,stride,lines);
    
    readColorComponent(p,src, width, height, stride, lines, src_stride);
  
     
  	int   resamplingMethod            = 0;
  	int   resamplingMode              = 0;   
  	int   resampling                  = 0;   
  	int   upsampling                  = 0;  
  	//ResizeParameters cRP;
 	cRP.m_bRefLayerFrameMbsOnlyFlag   = true;
  	cRP.m_bFrameMbsOnlyFlag           = true;
    cRP.m_bRefLayerFieldPicFlag       = false;
    cRP.m_bFieldPicFlag               = false;
    cRP.m_bRefLayerBotFieldFlag       = false;
    cRP.m_bBotFieldFlag               = false;
    cRP.m_bRefLayerIsMbAffFrame       = false;
    cRP.m_bIsMbAffFrame               = false;
    cRP.m_iRefLayerChromaPhaseX       = -1;
    cRP.m_iRefLayerChromaPhaseY       = 0;
    cRP.m_iChromaPhaseX               = -1;
    cRP.m_iChromaPhaseY               = 0;
    cRP.m_iRefLayerFrmWidth           = win;
    cRP.m_iRefLayerFrmHeight          = hin;
    cRP.m_iScaledRefFrmWidth          = 0;
    cRP.m_iScaledRefFrmHeight         = 0;
    cRP.m_iFrameWidth                 = wout;
    cRP.m_iFrameHeight                = hout;
    cRP.m_iLeftFrmOffset              = 0;
    cRP.m_iTopFrmOffset               = 0;
    cRP.m_iExtendedSpatialScalability = 0;
    cRP.m_iLevelIdc                   = 0; 
    cRP.m_iScaledRefFrmWidth  = gMax( cRP.m_iRefLayerFrmWidth,  cRP.m_iFrameWidth  );
    cRP.m_iScaledRefFrmHeight = gMax( cRP.m_iRefLayerFrmHeight, cRP.m_iFrameHeight );   
    upsampling  = ( cRP.m_iRefLayerFrmWidth < cRP.m_iFrameWidth ) || ( cRP.m_iRefLayerFrmHeight < cRP.m_iFrameHeight );

    DownConvert *cDownConvert = (DownConvert*)malloc(sizeof(DownConvert));//
    {
      int imaxWidth  = gMax( cRP.m_iRefLayerFrmWidth,  cRP.m_iFrameWidth  );
      int imaxHeight = gMax( cRP.m_iRefLayerFrmHeight, cRP.m_iFrameHeight );
      int iminWidth  = gMin( cRP.m_iRefLayerFrmWidth,  cRP.m_iFrameWidth  );
      int iminHeight = gMin( cRP.m_iRefLayerFrmHeight, cRP.m_iFrameHeight );
      int iminWRnd16 = ( ( iminWidth  + 15 ) >> 4 ) << 4;
      int iminHRnd32 = ( ( iminHeight + 31 ) >> 5 ) << 5;
      imaxWidth      = ( ( imaxWidth  * iminWRnd16 + ( iminWidth  << 4 ) - 1 ) / ( iminWidth  << 4 ) ) << 4;
      imaxHeight     = ( ( imaxHeight * iminHRnd32 + ( iminHeight << 4 ) - 1 ) / ( iminHeight << 4 ) ) << 4;
  	  cDownConvertInit(maxWidth, maxHeight,0,cDownConvert); 
    } 
    resampleFrame(p, cDownConvert, cRP, resamplingMethod, resamplingMode, resampling, upsampling, 1,stride );
    //writeCsp(p,outputFile,wout,hout,stride);
   //writeCsp(pixel* src, pixel* dst, int width, int height,int stride)
 	//h->mc.plane_copy(dst,dst_stride,(pixel*)p,stride, wout, hout);
 	writeCsp(p,outputfile,wout,hout,stride);

}

int xIsInCropWindow(MotionUpsampling * mo_up)
{
  int i_mbaff = (mo_up->m_rc_resize_params->m_bIsMbAffFrame ? 1:0);
  int i_field = (i_mbaff || mo_up->m_rc_resize_params->m_bFieldPicFlag?1:0);
  int i_mb_y0 = (mo_up->i_mby_curr >> i_mbaff) << i_field;
  int i_mb_y1 = i_mb_y0 + i_field;
  ROFRS(mo_up->i_mbx_curr >= mo_up->i_mb_x0_crop_frm && mo_up->i_mbx_curr < mo_up->i_mb_x1_crop_frm,0);
  ROFRS(i_mb_y0 >= mo_up->i_mb_y0_crop_frm && i_mb_y1 < mo_up->i_mb_y1_crop_frm,0);
  
}

int xInitMb(MotionUpsampling* mo_up,int i_mbx_curr,int i_mby_curr,x264_t* h)
{
  mo_up->i_mbx_curr = i_mbx_curr;
  mo_up->i_mby_curr = i_mby_curr;
  mo_up->b_in_crop_window = xIsInCropWindow(mo_up);
  mo_up->b_intraBL = 0;
  //
  mo_up->i_fwd_bwd = 0;
  mo_up->b_res_pred_safe = mo_up->b_in_crop_window;
  //

  mo_up->b_aa_base_intra[0][0] = 0;
  mo_up->b_aa_base_intra[0][1] = 0;
  mo_up->b_aa_base_intra[1][0] = 0;
  mo_up->b_aa_base_intra[1][1] = 0;

  mo_up->i_aaai_ref_idx[0][0][0]  = BLOCK_NOT_PREDICTED;
  mo_up->i_aaai_ref_idx[0][0][1]  = BLOCK_NOT_PREDICTED;
  mo_up->i_aaai_ref_idx[0][1][0]  = BLOCK_NOT_PREDICTED;
  mo_up->i_aaai_ref_idx[0][1][1]  = BLOCK_NOT_PREDICTED;
  mo_up->i_aaai_ref_idx[1][0][0]  = BLOCK_NOT_PREDICTED;
  mo_up->i_aaai_ref_idx[1][0][1]  = BLOCK_NOT_PREDICTED;
  mo_up->i_aaai_ref_idx[1][1][0]  = BLOCK_NOT_PREDICTED;
  mo_up->i_aaai_ref_idx[1][1][1]  = BLOCK_NOT_PREDICTED;

  /*set field mode for SNR scalability*/
  if(mo_up->b_in_crop_window && (mo_up->b_scoeff_pred || mo_up->b_tcoeff_pred))
  {
    int i_field_pic = (mo_up->m_rc_resize_params->m_bFieldPicFlag? 1 :0);
	int i_bot_field = (mo_up->m_rc_resize_params->m_bBotFieldFlag?1:0);
	int i_mbx_base = mo_up->i_mbx_curr - (mo_up->m_rc_resize_params->m_iLeftFrmOffset >> 4);
	int i_mby_base = mo_up->i_mby_curr - ((mo_up->m_rc_resize_params->m_iLeftFrmOffset >> 4)>>i_field_pic);
	int i_mb_stride_base = (mo_up->m_rc_resize_params->m_iRefLayerFrmWidth >> 4) << i_field_pic;
	int i_mb_offset_base = (mo_up->m_rc_resize_params->m_iRefLayerFrmWidth >> 4) * i_bot_field;
	int i_mb_idx_base = i_mb_offset_base + i_mby_base * i_mb_stride_base + i_mbx_base;
 
	int b_field_mb_flag = h->mbBL.field[i_mb_idx_base];
	mo_up->b_curr_field_mb = b_field_mb_flag;
  }
}



int xGetRefLayerPartIdc(MotionUpsampling* mo_up,int iXInsideCurrMB,int iYInsideCurrMb,int* riPartIdc,x264_t* h)
{
   int iBase8x8MbPartIdx = 0;// index of top-left 8x8 block that is covered by the macroblock partition
   int iBase4x4SubMbPartIdx = 0;// index of top-left 4x4 block that is covered by the sub-macroblock partition (inside 8x8 block)
   int iBaseMbIdx = MSYS_INT_MAX;
   int iXInsideBaseMb = MSYS_INT_MAX;
   int iYInsideBaseMb = MSYS_INT_MAX;

   RNOK(xGetRefLayerMb(mo_up,iXInsideCurrMB,iYInsideCurrMb,&iBaseMbIdx ,&iXInsideBaseMb,&iYInsideBaseMb,h));

   
   if(IS_INTRA(h->mbBL.type[iBaseMbIdx]))
   {
   	 *riPartIdc = -1;
	 return m_nOK;
   }

   int iB8x8IdxBase = (( iYInsideBaseMb >> 3) << 1) + ( iXInsideBaseMb >> 3);
   int iB4x4IdxBase = ((( iYInsideBaseMb & 7 ) >> 2)  << 1) + (( iXInsideBaseMb & 7 ) >> 2);
   int eMbModeBase = h->mbBL.mb_mode[iBaseMbIdx];
   int eBlkModeBase = h->mbBL.blk_mode[iBaseMbIdx];
   int b_skip_or_direct = (eMbModeBase == MODE_SKIP && h->sh.i_type == SLICE_TYPE_B);


   if(b_skip_or_direct)
   {  //===== determine macroblock and sub-macroblock partition index for B_Skip and B_Direct_16x16 =====
   	  if(mo_up->i_ref_layer_dqid == 0)
   	  	{
   	  	  iBase8x8MbPartIdx = iB8x8IdxBase;
		  iBase4x4SubMbPartIdx = iB4x4IdxBase;
   	  	}
   }
   else
   {
      const uint16_t aauiNxNPartIdx[6][4] = 
      {
      	   	{0,0,0,0},// MODE_SKIP (P Slice)
      	   	{0,0,0,0},// MODE_16x16     or    BLK_8x8
      	   	{0,0,2,2},// MODE_16x8      or    BLK_8x4
      	   	{0,1,0,1},// MODE_8x16      or    BLK_4x8
      	   	{0,1,2,3},// MODE_8x8       or    BLK_4x4
      	   	{0,1,2,3}// MODE_8x8ref0
      };


	  iBase8x8MbPartIdx = aauiNxNPartIdx[eMbModeBase][iB8x8IdxBase];

	  if(eMbModeBase == MODE_8x8 || eMbModeBase == MODE_8x8ref0)
	  {
	    if(eBlkModeBase = BLK_SKIP)
	    {
	      if(mo_up->i_ref_layer_dqid == 0)
	      {
	        iBase4x4SubMbPartIdx = iB4x4IdxBase;
	      }
	      	
	    }
		else
		{
		  iBase4x4SubMbPartIdx = aauiNxNPartIdx[eBlkModeBase - BLK_8x8 + 1][iB4x4IdxBase];
		}
	  }
   }
   
   //===== this function implements subclause G.8.6.1.1 =====
   int iBase4x4BlkX = ((iBase8x8MbPartIdx & 1) << 1) + (iBase4x4SubMbPartIdx & 1);
   int iBase4x4BlkY = ((iBase8x8MbPartIdx >> 1) << 1) + (iBase4x4SubMbPartIdx >> 1);
   int iBase4x4BlkIdx = (iBase4x4BlkY << 2) + iBase4x4BlkX;

   riPartIdc = (iBaseMbIdx << 4) + iBase4x4BlkIdx;

   return m_nOK;

   
}

int xGetInitialBaseRefIdxAndMv(MotionUpsampling* mo_up,int i4x4BlkX,int i4x4BlkY,int eListIdx,int *riRefIdx,int *rcMv,x264_t* h)
{

}

int xGetRefLayerMb(MotionUpsampling *mo_up,int iXInsideCurrMb,int iYInsideCurrMb,int * riBaseMbIdx,int * riXInsideBaseMb,int *riYInsideBaseMb,x264_t *h)
{
}



int xSetPartIdcArray(MotionUpsampling *mo_up,x264_t * h)
{
  ROFRS(mo_up->b_in_crop_window,m_nOK);

  /*determine all 16 initial partition indices */
  	{
  	  int b_intraBL = 1;
	  for(int iY = 0;iY < 4;iY++)
	  	for(int iX = 0;iX < 4;iX++)
	  	{
	  	  RNOK(xGetRefLayerPartIdc(mo_up,(iX << 2) + 1,(iY << 2) + 1,&(mo_up->i_aai_part_idc[iX][iY]),h));
		  b_intraBL = (b_intraBL && (mo_up->i_aai_part_idc[iX][iY] == -1));

		  if(b_intraBL)
		  {

		   mo_up->mb_mode = INTRA_BL;

		  }	

	  	}

	  ROTRS(b_intraBL,m_nOK);
	  //ROTRS( m_cPosCalc.m_bRSChangeFlag,  Err::m_nOK );

	  //===== replace values of "-1" on a 4x4 block basis =====

	  {
	    for(int iYP = 0; iYP < 2;iYP++)
	    for(int iXP = 0;iXP < 2;iXP++)
	    {

           int b_aaProcI4x4Blk[2][2] = {{0,0},{0,0}};

		   for(int iYS = 0; iYS < 2; iYS++)
		   for(int iXS = 0; iXS < 2; iXS++)
		   {
		     int iYC = (iYP << 1) + iYS;
			 int iXC = (iXP << 1) + iXS;

			 if(mo_up->i_aai_part_idc[iXC][iYC] == -1)
			 {
			    int iYSInv = 1 - iYS;
				int iXSInv = 1 - iXS;
				int iYCInv = (iYP << 1) + iYSInv;
				int iXCInv = (iXP << 1) + iXSInv;
				b_aaProcI4x4Blk[iXS][iYS] = 1;

				if(! b_aaProcI4x4Blk[iXSInv][iYS] && mo_up->i_aai_part_idc[iXCInv][iYC] != -1)
				{
				   mo_up->i_aai_part_idc[iXC][iYC] = mo_up->i_aai_part_idc[iXCInv][iYC];
				}

				else if(! b_aaProcI4x4Blk[iXS][iYSInv] && mo_up->i_aai_part_idc[iXC][iYCInv] != -1)
				{
				   mo_up->i_aai_part_idc[iXC][iYC] = mo_up->i_aai_part_idc[iXC][iYCInv];
				}
				else if(! b_aaProcI4x4Blk[iXSInv][iYSInv] && mo_up->i_aai_part_idc[iXCInv][iYCInv] != -1)
				{
				   mo_up->i_aai_part_idc[iXC][iYC] = mo_up->i_aai_part_idc[iXCInv][iYCInv];
				}
				
			 }
		   }
		  
	    }
	  }
	  
  	}



     //===== replace values of "-1" on an 8x8 block basis =====
    {
      int b_aaProcI8x8Blk[2][2] = {{0,0},{0,0}};

	  for(int iYP = 0; iYP < 2; iYP++)
	  for(int iXP = 0; iXP < 2; iXP++)
	  {
	    int iYPInv = 1 - iYP;
		int iXPInv = 1 - iXP; 
		int iYO = (iYP << 1);
		int iXO = (iXP <<1);
		int iYOInv = (2 - iYP);
		int iXOInv = (2 - iXP);

		if( mo_up->i_aai_part_idc[iXO][iYO] == -1)
		{
		  b_aaProcI8x8Blk[iXP][iYP] = 1;

		  if(! b_aaProcI8x8Blk[iXPInv][iYP] && mo_up->i_aai_part_idc[iXOInv][iYO] != -1)
		  {
		    b_aaProcI8x8Blk[iXO    ][iYO    ] = mo_up->i_aai_part_idc[iXOInv][iYO];
			b_aaProcI8x8Blk[iXO + 1][iYO    ] = mo_up->i_aai_part_idc[iXOInv][iYO];
			b_aaProcI8x8Blk[iXO    ][iYO + 1] = mo_up->i_aai_part_idc[iXOInv][iYO + 1];
			b_aaProcI8x8Blk[iXO + 1][iYO + 1] = mo_up->i_aai_part_idc[iXOInv][iYO + 1];
		  }

		  else if(! b_aaProcI8x8Blk[iXP][iYPInv] && mo_up->i_aai_part_idc[iXO][iYOInv] != -1)
		  {
		    b_aaProcI8x8Blk[iXO    ][iYO    ] = mo_up->i_aai_part_idc[iXO][iYOInv];
			b_aaProcI8x8Blk[iXO + 1][iYO    ] = mo_up->i_aai_part_idc[iXO + 1][iYOInv];
			b_aaProcI8x8Blk[iXO    ][iYO + 1] = mo_up->i_aai_part_idc[iXO][iYOInv];
			b_aaProcI8x8Blk[iXO + 1][iYO + 1] = mo_up->i_aai_part_idc[iXO + 1][iYOInv];
		  }

		  else if(! b_aaProcI8x8Blk[iXPInv][iYPInv] && mo_up->i_aai_part_idc[iXOInv][iYOInv] != -1)
		  {
		    b_aaProcI8x8Blk[iXO    ][iYO    ] = mo_up->i_aai_part_idc[iXOInv][iYOInv];
			b_aaProcI8x8Blk[iXO + 1][iYO    ] = mo_up->i_aai_part_idc[iXOInv][iYOInv];
			b_aaProcI8x8Blk[iXO    ][iYO + 1] = mo_up->i_aai_part_idc[iXOInv][iYOInv];
			b_aaProcI8x8Blk[iXO + 1][iYO + 1] = mo_up->i_aai_part_idc[iXOInv][iYOInv];
		  }
		}
	  }
    }


  return m_nOK;

  
}

int xGetMinRefIdx(int iRefIdxA,int iRefIdxB)
{
  ROTRS(iRefIdxA < 1,iRefIdxB);
  ROTRS(iRefIdxB < 1	,iRefIdxA);

  return X264_MIN(iRefIdxA,iRefIdxB);
  
}

int xMvDiff(int* mvA,int* mvB)
{
  return abs(mvA[0] - mvB[0]) + abs(mvA[1] - mvB[1]);
}

int* xAddMv(int*  mvA,int * mvB)
{
  int cNewMv[2];
  cNewMv[0] = mvA[0] + mvB[0];
  cNewMv[1] = mvA[1] + mvB[1];

  return cNewMv;
}

int xMvCopy(int *mvDes,int *mvSrc)
{
  mvDes[0] = mvSrc[0];
  mvDes[1] = mvSrc[1];
  return m_nOK;
}


int xMvLeftShift(int* mv,short s)
{
  mv[0] = mv[0] << s;
  mv[1] = mv[1] << s;

  return m_nOK;
}

int xMvRightShift(int* mv,short s)
{
  mv[0] = mv[0] >> s;
  mv[1] = mv[1] >> s;

  return m_nOK;
}


int xGetRefIdxAndInitialMvPred(MotionUpsampling *mo_up,int eListIdx,x264_t * h)
{
  //===== get initial predictors for reference indices and motion vectors =====
  	{

	   for(int i4x4BlkY = 0;i4x4BlkY < 4;i4x4BlkY++)
	   for(int i4x4BlkX = 0;i4x4BlkX < 4;i4x4BlkX++)
	   {
	     RNOK(xGetInitialBaseRefIdxAndMv(mo_up,i4x4BlkX,i4x4BlkY,mo_up->i_aai_part_idc[i4x4BlkX][i4x4BlkY],
		 	                                         &(mo_up->i_aai_ref_idx_temp[i4x4BlkX][i4x4BlkY]),&(mo_up->i_aaac_mv[eListIdx][i4x4BlkX][i4x4BlkY]),h));
	   }
	   
  	}

  //===== set reference indices =====
   mo_up->i_aaai_ref_idx[eListIdx][0][0] = mo_up->i_aai_ref_idx_temp[0][0];
   mo_up->i_aaai_ref_idx[eListIdx][0][1] = mo_up->i_aai_ref_idx_temp[0][2];
   mo_up->i_aaai_ref_idx[eListIdx][1][0] = mo_up->i_aai_ref_idx_temp[2][0];
   mo_up->i_aaai_ref_idx[eListIdx][1][1] = mo_up->i_aai_ref_idx_temp[2][2];
   
   // ROTRS( m_cMvScale.m_bRSChangeFlag, Err::m_nOK );

  

   //===== merge reference indices and modify motion vectors accordingly =====
   for( int i8x8BlkY = 0; i8x8BlkY < 2; i8x8BlkY++ )
   for( int i8x8BlkX = 0; i8x8BlkX < 2; i8x8BlkX++ )
   {
	 //----- determine reference indices -----
	 for( int i4x4BlkY = 0; i4x4BlkY < 2; i4x4BlkY++ )
	 for( int i4x4BlkX = 0; i4x4BlkX < 2; i4x4BlkX++ )
	 {
	   int iY  = ( i8x8BlkY << 1 ) + i4x4BlkY;
	   int iX  = ( i8x8BlkX << 1 ) + i4x4BlkX;
	   mo_up->i_aaai_ref_idx[eListIdx][i8x8BlkX][i8x8BlkY] = xGetMinRefIdx(  mo_up->i_aai_ref_idx_temp[eListIdx][i8x8BlkX][i8x8BlkY],  mo_up->i_aai_ref_idx_temp[eListIdx][iX][iY] );
	 }
   
	 //----- update motion vectors -----
	 for( int iYS = 0; iYS < 2; iYS++ )
	 for( int iXS = 0; iXS < 2; iXS++ )
	 {
	   int iY = ( i8x8BlkY << 1 ) + iYS;
	   int iX = ( i8x8BlkX << 1 ) + iXS;
   
	   if( mo_up->i_aaai_ref_idx[eListIdx][i8x8BlkX][i8x8BlkY] != mo_up->i_aai_ref_idx_temp[iX][iY] )
	   {
		 int iYInv = ( i8x8BlkY << 1 ) + 1 - iYS;
		 int iXInv = ( i8x8BlkX << 1 ) + 1 - iXS;
   
		 if( mo_up->i_aaai_ref_idx[eListIdx][i8x8BlkX][i8x8BlkY] == mo_up->i_aai_ref_idx_temp[iXInv][iY] )
		 {
		   xMvCopy( mo_up->i_aaac_mv[eListIdx][iX][iY],mo_up->i_aaac_mv[eListIdx][iXInv][iY]);
		  // mo_up->i_aaac_mv[eListIdx][iX][iY] = mo_up->i_aaac_mv[eListIdx][iXInv][iY];
		 }
		 else if(  mo_up->i_aaai_ref_idx[eListIdx][i8x8BlkX][i8x8BlkY] == mo_up->i_aai_ref_idx_temp[iX][iYInv] )
		 {
		   xMvCopy( mo_up->i_aaac_mv[eListIdx][iX][iY],mo_up->i_aaac_mv[eListIdx][iX][iYInv]);
		   //mo_up->i_aaac_mv[eListIdx][iX][iY] = mo_up->i_aaac_mv[eListIdx][iX][iYInv];
		 }
		 else
		 {
		   ROF( (mo_up->i_aaai_ref_idx[eListIdx][i8x8BlkX][i8x8BlkY] == mo_up->i_aai_ref_idx_temp[eListIdx][iXInv][iYInv]));
		   xMvCopy( mo_up->i_aaac_mv[eListIdx][iX][iY],mo_up->i_aaac_mv[eListIdx][iXInv][iYInv]);
		   //mo_up->i_aaac_mv[eListIdx][iX][iY] = mo_up->i_aaac_mv[eListIdx][iXInv][iYInv];
		 }
	   }
	 }
   }
   
   return m_nOK;

  
}


int xDeriveBlockModeAndUpdateMv(MotionUpsampling *mo_up,int i8x8BlkIdx)
{
  //Int   iAbsMvDiffThreshold = ( m_cMvScale.m_bRSChangeFlag ? 0 : 1 );
  int iAbsMvDiffThreshold = 1;
  int iXO = (i8x8BlkIdx & 1) << 1;
  int iYO = (i8x8BlkIdx >> 1) << 1;
  int b_hor_match = 1;
  int b_ver_match = 1;
  int b_8x8_match = 1;

   //===== unify 8x8 blocks when direct_8x8_inference_flag is equal to 1 =====
   if(mo_up->b_direct8x8_inference &&  /*!m_cMvScale.m_bRSChangeFlag*/ mo_up->e_slice_type == SLICE_TYPE_B)
   {
   	  int iXC = (iXO >> 1) * 3;
	  int iYC = (iYO >> 1) *3;
	  for(int iListIdx = 0; iListIdx < mo_up->i_max_list_idx;iListIdx++)
	  {
	    int* cTmpMv = mo_up->i_aaac_mv[iListIdx][iXC][iYC];
		xMvCopy(mo_up->i_aaac_mv[iListIdx][iXO][iYO],cTmpMv);
		xMvCopy(mo_up->i_aaac_mv[iListIdx][iXO + 1][iYO],cTmpMv);
		xMvCopy(mo_up->i_aaac_mv[iListIdx][iXO][iYO + 1],cTmpMv);
		xMvCopy(mo_up->i_aaac_mv[iListIdx][iXO + 1][iYO + 1],cTmpMv);
		//mo_up->i_aaac_mv[iListIdx][iXO][iYO] = cTmpMv;
		//mo_up->i_aaac_mv[iListIdx][iXO + 1][iYO] = cTmpMv;
		//mo_up->i_aaac_mv[iListIdx][iXO][iYO + 1] = cTmpMv;
		//mo_up->i_aaac_mv[iListIdx][iXO + 1][iYO+ 1] = cTmpMv;
	  }  
   }

   //===== derive partition size =====
   	{
   	  for(int iListIdx = 0; iListIdx < mo_up->i_max_list_idx;iListIdx++)
   	  {
   	    int b_hor1_match = xMvDiff(mo_up->i_aaac_mv[iListIdx][iXO][iYO],mo_up->i_aaac_mv[iListIdx][iXO + 1][iYO]) <= iAbsMvDiffThreshold;
		int b_hor2_match = xMvDiff(mo_up->i_aaac_mv[iListIdx][iXO][iYO + 1],mo_up->i_aaac_mv[iListIdx][iXO + 1][iYO + 1]) <= iAbsMvDiffThreshold;
		int b_ver1_match = xMvDiff(mo_up->i_aaac_mv[iListIdx][iXO][iYO],mo_up->i_aaac_mv[iListIdx][iXO][iYO + 1]) <= iAbsMvDiffThreshold;
		int b_ver2_match = xMvDiff(mo_up->i_aaac_mv[iListIdx][iXO + 1][iYO],mo_up->i_aaac_mv[iListIdx][iXO + 1][iYO + 1]) <= iAbsMvDiffThreshold;
        int b_diag_match = xMvDiff(mo_up->i_aaac_mv[iListIdx][iXO][iYO],mo_up->i_aaac_mv[iListIdx][iXO + 1][iYO + 1]) <= iAbsMvDiffThreshold;
        b_8x8_match = b_8x8_match && b_hor1_match && b_ver1_match && b_diag_match;
		b_hor_match = b_hor_match && b_hor1_match && b_hor2_match;
		b_ver_match = b_ver_match && b_ver1_match && b_ver2_match;
	  }

	  const int aiBlkMode[4] = {BLK_4x4,BLK_8x4,BLK_4x8,BLK_8x8};
	  mo_up->blk_mode[i8x8BlkIdx] = aiBlkMode[ b_8x8_match ? 3: b_hor_match? 1:b_ver_match?2:0];
   	}


   //ROTRS( m_cMvScale.m_bRSChangeFlag,          Err::m_nOK );
   ROTRS(mo_up->blk_mode[i8x8BlkIdx] == BLK_4x4,m_nOK);
   //===== combine motion vectors ===== //===== combine motion vectors =====
   {
     for(int iListIdx = 0;iListIdx < mo_up->i_max_list_idx;iListIdx++)
     {
        int* cNewMv = {0,0}, *cNewMvA = {0,0},*cNewMvB = {0,0},*mv0 = {0,0},*mv1 = {1,1};
       switch(mo_up->blk_mode[i8x8BlkIdx])
       	{
       	  case BLK_8x8:
		  	 //int cNewMv[2] = {0,0};
			 //int mv0[2] = {2,2};
			 cNewMv = xAddMv(mo_up->i_aaac_mv[iListIdx][iXO][iYO],mo_up->i_aaac_mv[iListIdx][iXO + 1][iYO]);
			 cNewMv = xAddMv(cNewMv,mo_up->i_aaac_mv[iListIdx][iXO][iYO + 1]);
			 cNewMv = xAddMv(cNewMv,mo_up->i_aaac_mv[iListIdx][iXO + 1][iYO + 1]);
			 cNewMv = xAddMv(cNewMv,mv0);
             xMvRightShift(cNewMv,2);
			 xMvCopy(mo_up->i_aaac_mv[iListIdx][iXO][iYO],cNewMv);
			 xMvCopy(mo_up->i_aaac_mv[iListIdx][iXO + 1][iYO],cNewMv);
			 xMvCopy(mo_up->i_aaac_mv[iListIdx][iXO][iYO + 1],cNewMv);
			 xMvCopy(mo_up->i_aaac_mv[iListIdx][iXO + 1][iYO + 1],cNewMv);
			 
		  	break;
		  case BLK_8x4:

			//int cNewMvA[2] = {0,0};
			//int cNewMvB[2] = {0,0};
			//int mv0[2] = {1,1};
			cNewMvA = xAddMv(mo_up->i_aaac_mv[iListIdx][iXO][iYO],mo_up->i_aaac_mv[iListIdx][iXO + 1][iYO]);
			cNewMvA = xAddMv(cNewMvA,mv1);
            xMvRightShift(cNewMvA,1);
			
			cNewMvB = xAddMv(mo_up->i_aaac_mv[iListIdx][iXO][iYO + 1],mo_up->i_aaac_mv[iListIdx][iXO + 1][iYO + 1]);
			cNewMvB = xAddMv(cNewMvB,mv1);
			xMvRightShift(cNewMvB,1);
			
			xMvCopy(mo_up->i_aaac_mv[iListIdx][iXO][iYO],cNewMvA);
			xMvCopy(mo_up->i_aaac_mv[iListIdx][iXO + 1][iYO],cNewMvA);
			xMvCopy(mo_up->i_aaac_mv[iListIdx][iXO][iYO + 1],cNewMvB);
			xMvCopy(mo_up->i_aaac_mv[iListIdx][iXO + 1][iYO + 1],cNewMvB);
			
		  	break;
		  case BLK_4x8:
		  //	int cNewMvA[2] = {0,0},cNewMvB[2] = {0,0},mv0[2] = {1,1};
			cNewMvA = xAddMv(mo_up->i_aaac_mv[iListIdx][iXO][iYO],mo_up->i_aaac_mv[iListIdx][iXO][iYO + 1]);
			cNewMvA = xAddMv(cNewMvA,mv1);
            xMvRightShift(cNewMvA,1);
			
			cNewMvB = xAddMv(mo_up->i_aaac_mv[iListIdx][iXO + 1][iYO],mo_up->i_aaac_mv[iListIdx][iXO + 1][iYO + 1]);
			cNewMvB = xAddMv(cNewMvB,mv1);
			xMvRightShift(cNewMvB,1);
			
			xMvCopy(mo_up->i_aaac_mv[iListIdx][iXO][iYO],cNewMvA);
			xMvCopy(mo_up->i_aaac_mv[iListIdx][iXO + 1][iYO],cNewMvB);
			xMvCopy(mo_up->i_aaac_mv[iListIdx][iXO][iYO + 1],cNewMvA);
			xMvCopy(mo_up->i_aaac_mv[iListIdx][iXO + 1][iYO + 1],cNewMvB);
			  break;
		  default:
		  	ROT(1);
			break;

       	}
     }

	 return m_nOK;
   }

   
}

int xMbdataClear(x264_t * h)
{
  h->mbEL1.mb_mode = MODE_SKIP;
  h->mbEL1.b_BL_skip_flag = 0;
  h->mbEL1.i_mb_mode = MODE_SKIP;
  h->mbEL1.i_chroma_pred_mode = 0;
  h->mbEL1.i_cbp_luma = 0;
  h->mbEL1.i_cbp_chroma = 0;
  h->mbEL1.b_residual_pred_flag = 0;
  h->mbEL1.b_transform_8x8 = 0;
  h->mbEL1.b_in_crop_window_flag = 0;
  h->mbEL1.i_blk_mode[0] = h->mbEL1.i_blk_mode[1] = h->mbEL1.i_blk_mode[2]
                             = h->mbEL1.i_blk_mode[3] = BLK_8x8;

/* foreach intra4x4_pred_mode, its value is I_PRED_4x4_DC - BY -MING*/
  for(int i = 0;i < 16; i++)
  {
    
  }

  return 0;
}

int xSetPredMbData(MotionUpsampling *mo_up,x264_t * h)
{
  //=== get MbData reference===
  int iFieldPic = (mo_up->m_rc_resize_params->m_bFieldPicFlag? 1:0);
  int iBotField = (mo_up->m_rc_resize_params->m_bBotFieldFlag?1:0);
  int iMbStride = (mo_up->m_rc_resize_params->m_iFrameWidth >> 4) << iFieldPic;
  int iMbOffset = (mo_up->m_rc_resize_params->m_iFrameWidth >> 4) * iBotField;
  int iMbIdx = iMbOffset = mo_up->i_mby_curr* iMbStride + mo_up->i_mbx_curr;
  int s8x8 = h->mbEL1.i_b8_stride;
  int s4x4 = h->mbEL1.i_b4_stride;
  int iMb8x8Idx = 2*(mo_up->i_mby_curr * s8x8 + mo_up->i_mbx_curr);
  int iMb4x4Idx = 4*(mo_up->i_mby_curr * s4x4 + mo_up->i_mbx_curr);
/*=== reset MbDataStruct data - BY MING*/
  xMbdataClear(h);

  /*=== set motion data(ref idx & motion vectors) == - BY MING*/
  if(! mo_up->b_in_crop_window || mo_up->b_intraBL)
  {
    for(int iListIdx = 0;iListIdx < 2; iListIdx++)
    {
      // apcMotion[iListIdx]->clear( BLOCK_NOT_PREDICTED );
      
    }
  }

  else
  {
    int iListIdx = 0;
	for(; iListIdx < mo_up->i_max_list_idx;iListIdx++)
	{


	 for(int i8x8Idx = 0;i8x8Idx < 4;i8x8Idx++)
	 {
	   int i8x8XIdx = i8x8Idx & 0x11;
	   int i8x8YIdx = i8x8Idx >> 0x01;
	   h->mbEL1.ref[iListIdx][iMb8x8Idx + i8x8YIdx*s8x8 + i8x8XIdx] = mo_up->i_aaai_ref_idx[iListIdx][i8x8XIdx][i8x8YIdx];
	 }
	  


	  for(int i4x4Idx = 0;i4x4Idx < 16; i4x4Idx++)
	  {
	    int i4x4Xidx = i4x4Idx & 0x011;
		int i4x4YIdx = i4x4Idx >> 0x010;
		xMvCopy(h->mbEL1.mv[iListIdx][iMb4x4Idx + i4x4YIdx*s4x4 + i4x4Xidx] , mo_up->i_aaac_mv[iListIdx][i4x4Xidx][i4x4YIdx]);
	  }
	}

	// apcMotion[iListIdx]->clear( BLOCK_NOT_PREDICTED );
  }

   // apcMotion[0]->setFieldMode( m_bCurrFieldMb );
  //  apcMotion[1]->setFieldMode( m_bCurrFieldMb );

  //set general Mb data - BY MING
  h->mbEL1.b_in_crop_window_flag = mo_up->b_in_crop_window;
  h->mbEL1.field[iMbIdx] = mo_up->b_curr_field_mb;
  //rcMbData.setSafeResPred     ( m_bResPredSafe );
  if(mo_up->b_in_crop_window)
  {
    h->mbEL1.mb_mode[iMbIdx] = mo_up->mb_mode;
	//rcMbData.setFwdBwd( (UShort)m_uiFwdBwd );
	for(int  blk = 0;blk < 4; blk++)
	{
	  h->mbEL1.i_blk_mode[blk] =  mo_up->blk_mode[blk];
	}
  }



  if((mo_up->b_scoeff_pred || mo_up->b_tcoeff_pred) && mo_up->b_in_crop_window)
  {
    int iMbXBase = mo_up->i_mbx_curr - (mo_up->m_rc_resize_params->m_iLeftFrmOffset >> 4);
	int iMbYBase = mo_up->i_mby_curr - ((mo_up->m_rc_resize_params->m_iLeftFrmOffset >> 4) >> iFieldPic);
	int     iMbStrideBase = ( mo_up->m_rc_resize_params->m_iRefLayerFrmWidth >> 4 ) << iFieldPic;
    int     iMbOffsetBase = ( mo_up->m_rc_resize_params->m_iRefLayerFrmWidth >> 4 )  * iBotField;
    int     iMbIdxBase    = iMbOffsetBase + iMbYBase * iMbStrideBase + iMbXBase;
   // rcMbData.copyTCoeffs    ( rcMbDataBase );

   h->mbEL1.BL_skip[iMbIdx] = h->mbBL.BL_skip[iMbIdxBase];
   if(mo_up->b_intraBL)
   {
     h->mbEL1.mb_mode[iMbIdx] = h->mbBL.mb_mode[iMbIdxBase];
   }
   if(mo_up->b_tcoeff_pred)
   {
    h->mbEL1.intra16x16_pred_mode[iMbIdx] = h->mbBL.intra16x16_pred_mode[iMbIdxBase];
	h->mbEL1.chroma_pred_mode[iMbIdx] = h->mbBL.chroma_pred_mode[iMbIdxBase];
		// COPY TRANSFORM SIZE
	h->mbEL1.transform8x8_size[iMbIdx] = h->mbBL.transform8x8_size[iMbIdxBase];
	h->mbEL1.cbp[iMbIdx] = h->mbBL.cbp[iMbIdxBase];
	h->mbEL1.qp[iMbIdx] = h->mbBL.qp[iMbIdxBase];
	 //m_ucQp4LF = rcMbData.m_ucQp4LF;
   }

  }
  return m_nOK;
}


/**/
int xUpsampleMotion(MotionUpsampling *mo_up,ResizeParameters* pcResizeParams,int b_field_resampling,
                                                                        int b_residual_pred_check,int i_mv_threshold,x264_t* h)
{
   //ROF(h->sh);
   ROF(pcResizeParams);
   mo_up->b_check_residual_pred = b_residual_pred_check;
   mo_up->b_direct8x8_inference = h->sh.pps->b_transform_8x8_mode;
   mo_up->i_mv_threshold = i_mv_threshold;
   mo_up->b_curr_field_mb = b_field_resampling;
   mo_up->e_slice_type = h->sh.i_type;
   mo_up->i_ref_layer_dqid = h->sh.i_ref_layer_dq_id;
   mo_up->i_max_list_idx = mo_up->e_slice_type == SLICE_TYPE_B ? 2:(mo_up->e_slice_type == SLICE_TYPE_I?0:1);
}







