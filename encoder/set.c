/*****************************************************************************
 * set: header writing
 *****************************************************************************
 * Copyright (C) 2003-2013 x264 project
 *
 * Authors: Laurent Aimar <fenrir@via.ecp.fr>
 *          Loren Merritt <lorenm@u.washington.edu>
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

#include "common/common.h"
#include "set.h"

#define bs_write_ue bs_write_ue_big

// Indexed by pic_struct values
static const uint8_t num_clock_ts[10] = { 0, 1, 1, 1, 2, 2, 3, 3, 2, 3 };
const static uint8_t avcintra_uuid[] = {0xF7, 0x49, 0x3E, 0xB3, 0xD4, 0x00, 0x47, 0x96, 0x86, 0x86, 0xC9, 0x70, 0x7B, 0x64, 0x37, 0x2A};

static void transpose( uint8_t *buf, int w )
{
    for( int i = 0; i < w; i++ )
        for( int j = 0; j < i; j++ )
            XCHG( uint8_t, buf[w*i+j], buf[w*j+i] );
}

static void scaling_list_write( bs_t *s, x264_pps_t *pps, int idx )
{
    const int len = idx<4 ? 16 : 64;
    const uint8_t *zigzag = idx<4 ? x264_zigzag_scan4[0] : x264_zigzag_scan8[0];
    const uint8_t *list = pps->scaling_list[idx];
    const uint8_t *def_list = (idx==CQM_4IC) ? pps->scaling_list[CQM_4IY]
                            : (idx==CQM_4PC) ? pps->scaling_list[CQM_4PY]
                            : (idx==CQM_8IC+4) ? pps->scaling_list[CQM_8IY+4]
                            : (idx==CQM_8PC+4) ? pps->scaling_list[CQM_8PY+4]
                            : x264_cqm_jvt[idx];
    if( !memcmp( list, def_list, len ) )
        bs_write1( s, 0 );   // scaling_list_present_flag
    else if( !memcmp( list, x264_cqm_jvt[idx], len ) )
    {
        bs_write1( s, 1 );   // scaling_list_present_flag
        bs_write_se( s, -8 ); // use jvt list
    }
    else
    {
        int run;
        bs_write1( s, 1 );   // scaling_list_present_flag

        // try run-length compression of trailing values
        for( run = len; run > 1; run-- )
            if( list[zigzag[run-1]] != list[zigzag[run-2]] )
                break;
        if( run < len && len - run < bs_size_se( (int8_t)-list[zigzag[run]] ) )
            run = len;

        for( int j = 0; j < run; j++ )
            bs_write_se( s, (int8_t)(list[zigzag[j]] - (j>0 ? list[zigzag[j-1]] : 8)) ); // delta

        if( run < len )
            bs_write_se( s, (int8_t)-list[zigzag[run]] );
    }
}

void x264_sei_write( bs_t *s, uint8_t *payload, int payload_size, int payload_type )
{
    int i;

    bs_realign( s );

    for( i = 0; i <= payload_type-255; i += 255 )
        bs_write( s, 8, 255 );
    bs_write( s, 8, payload_type-i );

    for( i = 0; i <= payload_size-255; i += 255 )
        bs_write( s, 8, 255 );
    bs_write( s, 8, payload_size-i );

    for( i = 0; i < payload_size; i++ )
        bs_write( s, 8, payload[i] );

    bs_rbsp_trailing( s );
    bs_flush( s );
}

void x264_sps_init( x264_sps_t *sps, int i_id, x264_param_t *param )
{
    int csp = param->i_csp & X264_CSP_MASK;

    sps->i_id = i_id;
	/*sky 2014.08.29 添加 ssps 的宽度高度计算*/
	if(sps->i_nal_type == NAL_SPS)
		{
    			sps->i_mb_width = ( param->i_width + 15 ) / 16;
    			sps->i_mb_height= ( param->i_height + 15 ) / 16;
		}
	else if(sps->i_nal_type == NAL_UNIT_SUBSET_SPS)
		{
			sps->i_mb_width = ( param->i_width + 15 ) / 16 * 2 * (sps->i_id + 1);
    			sps->i_mb_height= ( param->i_height + 15 ) / 16 * 2 * (sps->i_id + 1);
		}
    sps->i_chroma_format_idc = csp >= X264_CSP_I444 ? CHROMA_444 :
                               csp >= X264_CSP_I422 ? CHROMA_422 : CHROMA_420;

    sps->b_qpprime_y_zero_transform_bypass = param->rc.i_rc_method == X264_RC_CQP && param->rc.i_qp_constant == 0;
    if( sps->b_qpprime_y_zero_transform_bypass || sps->i_chroma_format_idc == CHROMA_444 )
        sps->i_profile_idc  = PROFILE_HIGH444_PREDICTIVE;
    else if( sps->i_chroma_format_idc == CHROMA_422 )
        sps->i_profile_idc  = PROFILE_HIGH422;
    else if( BIT_DEPTH > 8 )
        sps->i_profile_idc  = PROFILE_HIGH10;
    else if( param->analyse.b_transform_8x8 || param->i_cqm_preset != X264_CQM_FLAT )
        sps->i_profile_idc  = PROFILE_HIGH;
    else if( param->b_cabac || param->i_bframe > 0 || param->b_interlaced || param->b_fake_interlaced || param->analyse.i_weighted_pred > 0 )
        sps->i_profile_idc  = PROFILE_MAIN;
    else
        sps->i_profile_idc  = PROFILE_BASELINE;

/*	sky 2014.08.28 添加ssps 情况下 i_profile_idc*/
   if( sps ->i_nal_type == NAL_UNIT_SUBSET_SPS)
  	{
		sps ->i_profile_idc = PROFILE_SCALABLE_HIGH_PROFILE;
  	}

    sps->b_constraint_set0  = sps->i_profile_idc == PROFILE_BASELINE;
    /* x264 doesn't support the features that are in Baseline and not in Main,
     * namely arbitrary_slice_order and slice_groups. */

   /*sky 2014.08.28 sps init 添加 ||sps->i_profile_idc == PROFILE_SCALABLE_HIGH_PROFILE;*/
    sps->b_constraint_set1  = (sps->i_profile_idc <= PROFILE_MAIN ||sps->i_profile_idc == PROFILE_SCALABLE_HIGH_PROFILE);

	/* Never set constraint_set2, it is not necessary and not used in real world. */
    sps->b_constraint_set2  = 0;
  
  /*sky 2014.08.21为了和输出码流对应从上述语句改为下述*/
  // sps->b_constraint_set3  = 0;
    sps->b_constraint_set3  = (sps->i_profile_idc <= PROFILE_MAIN ||sps->i_profile_idc == PROFILE_SCALABLE_HIGH_PROFILE);
 
    sps->i_level_idc = param->i_level_idc;
    if( param->i_level_idc == 9 && ( sps->i_profile_idc == PROFILE_BASELINE || sps->i_profile_idc == PROFILE_MAIN ) )
    {
        sps->b_constraint_set3 = 1; /* level 1b with Baseline or Main profile is signalled via constraint_set3 */
        sps->i_level_idc      = 11;
    }
    /* Intra profiles */
    if( param->i_keyint_max == 1 && sps->i_profile_idc > PROFILE_HIGH )
        sps->b_constraint_set3 = 1;

    sps->vui.i_num_reorder_frames = param->i_bframe_pyramid ? 2 : param->i_bframe ? 1 : 0;
    /* extra slot with pyramid so that we don't have to override the
     * order of forgetting old pictures */
    sps->vui.i_max_dec_frame_buffering =
    sps->i_num_ref_frames = X264_MIN(X264_REF_MAX, X264_MAX4(param->i_frame_reference, 1 + sps->vui.i_num_reorder_frames,
                            param->i_bframe_pyramid ? 4 : 1, param->i_dpb_size));
    sps->i_num_ref_frames -= param->i_bframe_pyramid == X264_B_PYRAMID_STRICT;
    if( param->i_keyint_max == 1 )
    {
        sps->i_num_ref_frames = 0;
        sps->vui.i_max_dec_frame_buffering = 0;
    }

  /*sky 2014.08.21，没有看懂邋i_level_idc的取值关系和意义，为了跟码流对上添加如下句子*/
	if(sps->i_nal_type == NAL_UNIT_SUBSET_SPS && sps->i_profile_idc == PROFILE_SCALABLE_HIGH_PROFILE)
		{
			sps ->i_level_idc =  21 ; // 这个 值是根据jsvm编的码流凑出来的
		}
	
    /* number of refs + current frame */
    int max_frame_num = sps->vui.i_max_dec_frame_buffering * (!!param->i_bframe_pyramid+1) + 1;
    /* Intra refresh cannot write a recovery time greater than max frame num-1 */
    if( param->b_intra_refresh )
    {
        int time_to_recovery = X264_MIN( sps->i_mb_width - 1, param->i_keyint_max ) + param->i_bframe - 1;
        max_frame_num = X264_MAX( max_frame_num, time_to_recovery+1 );
    }


	
    sps->i_log2_max_frame_num = 4;
    while( (1 << sps->i_log2_max_frame_num) <= max_frame_num )
        sps->i_log2_max_frame_num++;

    sps->i_poc_type = param->i_bframe || param->b_interlaced ? 0 : 2;
    if( sps->i_poc_type == 0 )
    {
        int max_delta_poc = (param->i_bframe + 2) * (!!param->i_bframe_pyramid + 1) * 2;
        sps->i_log2_max_poc_lsb = 4;
        while( (1 << sps->i_log2_max_poc_lsb) <= max_delta_poc * 2 )
            sps->i_log2_max_poc_lsb++;
    }

    sps->b_vui = 1;

    sps->b_gaps_in_frame_num_value_allowed = 0;
    sps->b_frame_mbs_only = !(param->b_interlaced || param->b_fake_interlaced);
    if( !sps->b_frame_mbs_only )
        sps->i_mb_height = ( sps->i_mb_height + 1 ) & ~1;
    sps->b_mb_adaptive_frame_field = param->b_interlaced;
    sps->b_direct8x8_inference = 1;

    sps->crop.i_left   = param->crop_rect.i_left;
    sps->crop.i_top    = param->crop_rect.i_top;
    sps->crop.i_right  = param->crop_rect.i_right + sps->i_mb_width*16 - param->i_width;
    sps->crop.i_bottom = (param->crop_rect.i_bottom + sps->i_mb_height*16 - param->i_height) >> !sps->b_frame_mbs_only;
    sps->b_crop = sps->crop.i_left  || sps->crop.i_top ||
                  sps->crop.i_right || sps->crop.i_bottom;

    sps->vui.b_aspect_ratio_info_present = 0;
    if( param->vui.i_sar_width > 0 && param->vui.i_sar_height > 0 )
    {
        sps->vui.b_aspect_ratio_info_present = 1;
        sps->vui.i_sar_width = param->vui.i_sar_width;
        sps->vui.i_sar_height= param->vui.i_sar_height;
    }

    sps->vui.b_overscan_info_present = param->vui.i_overscan > 0 && param->vui.i_overscan <= 2;
    if( sps->vui.b_overscan_info_present )
        sps->vui.b_overscan_info = ( param->vui.i_overscan == 2 ? 1 : 0 );

    sps->vui.b_signal_type_present = 0;
    sps->vui.i_vidformat = ( param->vui.i_vidformat >= 0 && param->vui.i_vidformat <= 5 ? param->vui.i_vidformat : 5 );
    sps->vui.b_fullrange = ( param->vui.b_fullrange >= 0 && param->vui.b_fullrange <= 1 ? param->vui.b_fullrange :
                           ( csp >= X264_CSP_BGR ? 1 : 0 ) );
    sps->vui.b_color_description_present = 0;

    sps->vui.i_colorprim = ( param->vui.i_colorprim >= 0 && param->vui.i_colorprim <=  9 ? param->vui.i_colorprim : 2 );
    sps->vui.i_transfer  = ( param->vui.i_transfer  >= 0 && param->vui.i_transfer  <= 15 ? param->vui.i_transfer  : 2 );
    sps->vui.i_colmatrix = ( param->vui.i_colmatrix >= 0 && param->vui.i_colmatrix <= 10 ? param->vui.i_colmatrix :
                           ( csp >= X264_CSP_BGR ? 0 : 2 ) );
    if( sps->vui.i_colorprim != 2 ||
        sps->vui.i_transfer  != 2 ||
        sps->vui.i_colmatrix != 2 )
    {
        sps->vui.b_color_description_present = 1;
    }

    if( sps->vui.i_vidformat != 5 ||
        sps->vui.b_fullrange ||
        sps->vui.b_color_description_present )
    {
        sps->vui.b_signal_type_present = 1;
    }

    /* FIXME: not sufficient for interlaced video */
    sps->vui.b_chroma_loc_info_present = param->vui.i_chroma_loc > 0 && param->vui.i_chroma_loc <= 5;
    if( sps->vui.b_chroma_loc_info_present )
    {
        sps->vui.i_chroma_loc_top = param->vui.i_chroma_loc;
        sps->vui.i_chroma_loc_bottom = param->vui.i_chroma_loc;
    }

    sps->vui.b_timing_info_present = param->i_timebase_num > 0 && param->i_timebase_den > 0;

    if( sps->vui.b_timing_info_present )
    {
        sps->vui.i_num_units_in_tick = param->i_timebase_num;
        sps->vui.i_time_scale = param->i_timebase_den * 2;
        sps->vui.b_fixed_frame_rate = !param->b_vfr_input;
    }

    sps->vui.b_vcl_hrd_parameters_present = 0; // we don't support VCL HRD
    sps->vui.b_nal_hrd_parameters_present = !!param->i_nal_hrd;
    sps->vui.b_pic_struct_present = param->b_pic_struct;

    // NOTE: HRD related parts of the SPS are initialised in x264_ratecontrol_init_reconfigurable

    sps->vui.b_bitstream_restriction = param->i_keyint_max > 1;
    if( sps->vui.b_bitstream_restriction )
    {
        sps->vui.b_motion_vectors_over_pic_boundaries = 1;
        sps->vui.i_max_bytes_per_pic_denom = 0;
        sps->vui.i_max_bits_per_mb_denom = 0;
        sps->vui.i_log2_max_mv_length_horizontal =
        sps->vui.i_log2_max_mv_length_vertical = (int)log2f( X264_MAX( 1, param->analyse.i_mv_range*4-1 ) ) + 1;
    }
 /*sky 2014.08.20 sps init extension*/
 if(sps ->i_nal_type == NAL_SPS)
 	{
		sps ->b_adaptive_tcoeff_level_prediction_flag = 0;
		sps ->b_seq_tcoeff_level_prediction_flag = 0;
 	}
 else 
 	{
		/*这两个值等于输入信息参数 � param

		*	sps ->b_adaptive_tcoeff_level_prediction_flag = 0;
		*	sps ->b_seq_tcoeff_level_prediction_flag = 0;
		*/
		
	
	}
 /*码流检测过*/
		sps ->b_inter_layer_deblocking_present = 0 ;//,码流初值0 
		sps ->b_slice_header_restriction_flag = 0; // jsvm
		sps ->b_svc_vui_parameters_present_flag = 0; // 写死0
		sps ->b_additional_extension2_flag = 0;//jsvm
		sps ->b_chroma_phase_x_plus1_flag = 0;//jsvm
		sps ->i_chroma_phase_y_plus1 = 1;//jsvm
		sps ->b_delta_pic_ord_always_zero_flag = 0;//jsvm
// 扩展结束
}

/*sky 2014.08.29 prifix write*/
void x264_prifix_write(bs_t *s, x264_nal_t *nal)
{
    bs_realign( s );
	
    bs_write1( s , nal ->b_svc_extension);
    bs_write1( s , nal ->b_idr_flag);
    bs_write( s , 6 , nal ->i_priority_id);
    bs_write1(s, nal ->b_no_inter_layer_pred_flag);
    bs_write(s, 3, nal ->i_dependency_id);
    bs_write(s, 4, nal ->i_quality_id);
    bs_write(s, 3, nal ->i_temporal_id);
    bs_write1(s, nal ->b_use_ref_base_pic_flag);
    bs_write1(s, nal ->b_discardable_flag);
    bs_write1(s, nal ->b_output_flag);
    bs_write(s, 2, nal ->i_reserved_three_2bits);
	
  /*这个地方我写死了*/
    bs_write1(s, 0);   // "store_ref_base_pic_flag" 
  /*  if( ( getUseRefBasePicFlag() || m_bStoreRefBasePicFlag ) && ! getIdrFlag() )
    {
      RNOK( m_cDecRefBasePicMarking.write( rcWriteIf ) );
    }*/
    
    bs_write1(s, 0);  //"additional_prefix_nal_unit_extension_flag"
    
    bs_rbsp_trailing( s );
    bs_flush( s );
}

void x264_sps_write( bs_t *s, x264_sps_t *sps )
{
    bs_realign( s );
    bs_write( s, 8, sps->i_profile_idc );
    bs_write1( s, sps->b_constraint_set0 );
    bs_write1( s, sps->b_constraint_set1 );
    bs_write1( s, sps->b_constraint_set2 );
    bs_write1( s, sps->b_constraint_set3 );

    bs_write( s, 4, 0 );    /* reserved */

    bs_write( s, 8, sps->i_level_idc );

    bs_write_ue( s, sps->i_id );


/*sky 2014.08.28 添加了 ||sps->i_profile_idc == PROFILE_SCALABLE_HIGH_PROFILE
		||sps->i_profile_idc == PROFILE_SCALABLE_BASELINE_PROFILE*/
    if( sps->i_profile_idc >= PROFILE_HIGH ||sps->i_profile_idc == PROFILE_SCALABLE_HIGH_PROFILE
		||sps->i_profile_idc == PROFILE_SCALABLE_BASELINE_PROFILE)
    {
        bs_write_ue( s, sps->i_chroma_format_idc );
        if( sps->i_chroma_format_idc == CHROMA_444 )
            bs_write1( s, 0 ); // separate_colour_plane_flag
        bs_write_ue( s, BIT_DEPTH-8 ); // bit_depth_luma_minus8
        bs_write_ue( s, BIT_DEPTH-8 ); // bit_depth_chroma_minus8
        bs_write1( s, sps->b_qpprime_y_zero_transform_bypass );
        bs_write1( s, 0 ); // seq_scaling_matrix_present_flag这个地方给他默认赋 0 了
    }


    bs_write_ue( s, sps->i_log2_max_frame_num - 4 );
    bs_write_ue( s, sps->i_poc_type );
    if( sps->i_poc_type == 0 )
        bs_write_ue( s, sps->i_log2_max_poc_lsb - 4 );
    bs_write_ue( s, sps->i_num_ref_frames );
    bs_write1( s, sps->b_gaps_in_frame_num_value_allowed );
    bs_write_ue( s, sps->i_mb_width - 1 );
    bs_write_ue( s, (sps->i_mb_height >> !sps->b_frame_mbs_only) - 1);
    bs_write1( s, sps->b_frame_mbs_only );
    if( !sps->b_frame_mbs_only )
        bs_write1( s, sps->b_mb_adaptive_frame_field );
    bs_write1( s, sps->b_direct8x8_inference );

    bs_write1( s, sps->b_crop );
    if( sps->b_crop )
    {
        int h_shift = sps->i_chroma_format_idc == CHROMA_420 || sps->i_chroma_format_idc == CHROMA_422;
        int v_shift = sps->i_chroma_format_idc == CHROMA_420;
        bs_write_ue( s, sps->crop.i_left   >> h_shift );
        bs_write_ue( s, sps->crop.i_right  >> h_shift );
        bs_write_ue( s, sps->crop.i_top    >> v_shift );
        bs_write_ue( s, sps->crop.i_bottom >> v_shift );
    }

    bs_write1( s, sps->b_vui );
    if( sps->b_vui )
    {
        bs_write1( s, sps->vui.b_aspect_ratio_info_present );
        if( sps->vui.b_aspect_ratio_info_present )
        {
            int i;
            static const struct { uint8_t w, h, sar; } sar[] =
            {
                // aspect_ratio_idc = 0 -> unspecified
                {  1,  1, 1 }, { 12, 11, 2 }, { 10, 11, 3 }, { 16, 11, 4 },
                { 40, 33, 5 }, { 24, 11, 6 }, { 20, 11, 7 }, { 32, 11, 8 },
                { 80, 33, 9 }, { 18, 11, 10}, { 15, 11, 11}, { 64, 33, 12},
                {160, 99, 13}, {  4,  3, 14}, {  3,  2, 15}, {  2,  1, 16},
                // aspect_ratio_idc = [17..254] -> reserved
                { 0, 0, 255 }
            };
            for( i = 0; sar[i].sar != 255; i++ )
            {
                if( sar[i].w == sps->vui.i_sar_width &&
                    sar[i].h == sps->vui.i_sar_height )
                    break;
            }
            bs_write( s, 8, sar[i].sar );
            if( sar[i].sar == 255 ) /* aspect_ratio_idc (extended) */
            {
                bs_write( s, 16, sps->vui.i_sar_width );
                bs_write( s, 16, sps->vui.i_sar_height );
            }
        }

        bs_write1( s, sps->vui.b_overscan_info_present );
        if( sps->vui.b_overscan_info_present )
            bs_write1( s, sps->vui.b_overscan_info );

        bs_write1( s, sps->vui.b_signal_type_present );
        if( sps->vui.b_signal_type_present )
        {
            bs_write( s, 3, sps->vui.i_vidformat );
            bs_write1( s, sps->vui.b_fullrange );
            bs_write1( s, sps->vui.b_color_description_present );
            if( sps->vui.b_color_description_present )
            {
                bs_write( s, 8, sps->vui.i_colorprim );
                bs_write( s, 8, sps->vui.i_transfer );
                bs_write( s, 8, sps->vui.i_colmatrix );
            }
        }

        bs_write1( s, sps->vui.b_chroma_loc_info_present );
        if( sps->vui.b_chroma_loc_info_present )
        {
            bs_write_ue( s, sps->vui.i_chroma_loc_top );
            bs_write_ue( s, sps->vui.i_chroma_loc_bottom );
        }

        bs_write1( s, sps->vui.b_timing_info_present );
        if( sps->vui.b_timing_info_present )
        {
            bs_write32( s, sps->vui.i_num_units_in_tick );
            bs_write32( s, sps->vui.i_time_scale );
            bs_write1( s, sps->vui.b_fixed_frame_rate );
        }

        bs_write1( s, sps->vui.b_nal_hrd_parameters_present );
        if( sps->vui.b_nal_hrd_parameters_present )
        {
            bs_write_ue( s, sps->vui.hrd.i_cpb_cnt - 1 );
            bs_write( s, 4, sps->vui.hrd.i_bit_rate_scale );
            bs_write( s, 4, sps->vui.hrd.i_cpb_size_scale );

            bs_write_ue( s, sps->vui.hrd.i_bit_rate_value - 1 );
            bs_write_ue( s, sps->vui.hrd.i_cpb_size_value - 1 );

            bs_write1( s, sps->vui.hrd.b_cbr_hrd );

            bs_write( s, 5, sps->vui.hrd.i_initial_cpb_removal_delay_length - 1 );
            bs_write( s, 5, sps->vui.hrd.i_cpb_removal_delay_length - 1 );
            bs_write( s, 5, sps->vui.hrd.i_dpb_output_delay_length - 1 );
            bs_write( s, 5, sps->vui.hrd.i_time_offset_length );
        }

        bs_write1( s, sps->vui.b_vcl_hrd_parameters_present );

        if( sps->vui.b_nal_hrd_parameters_present || sps->vui.b_vcl_hrd_parameters_present )
            bs_write1( s, 0 );   /* low_delay_hrd_flag */

        bs_write1( s, sps->vui.b_pic_struct_present );
        bs_write1( s, sps->vui.b_bitstream_restriction );
        if( sps->vui.b_bitstream_restriction )
        {
            bs_write1( s, sps->vui.b_motion_vectors_over_pic_boundaries );
            bs_write_ue( s, sps->vui.i_max_bytes_per_pic_denom );
            bs_write_ue( s, sps->vui.i_max_bits_per_mb_denom );
            bs_write_ue( s, sps->vui.i_log2_max_mv_length_horizontal );
            bs_write_ue( s, sps->vui.i_log2_max_mv_length_vertical );
            bs_write_ue( s, sps->vui.i_num_reorder_frames );
            bs_write_ue( s, sps->vui.i_max_dec_frame_buffering );
        }
    }
	/*sky 2014.08.28 添加ssps write*/
	 
	  if(sps ->i_nal_type == NAL_UNIT_SUBSET_SPS)
  	{
			if(sps ->i_profile_idc == PROFILE_SCALABLE_BASELINE_PROFILE || sps ->i_profile_idc == PROFILE_SCALABLE_HIGH_PROFILE)
			{

				bs_write1(s, sps ->b_inter_layer_deblocking_present) ;
				bs_write(s, 2, sps  ->i_extended_spatial_scalability) ;
				if( sps->i_chroma_format_idc == CHROMA_420 ||sps->i_chroma_format_idc == CHROMA_422 )
					bs_write ( s , 1 , sps ->b_chroma_phase_x_plus1_flag); 
				if( sps ->i_chroma_format_idc == CHROMA_420 )
					bs_write ( s, 2, sps ->i_chroma_phase_y_plus1);
				if( sps->i_extended_spatial_scalability == 1)   //ess_seq
					{
						if(sps->i_chroma_format_idc > 0)
							{
								bs_write(s, 1, sps ->b_seq_ref_layer_chroma_phase_x_plus1_flag);
								bs_write(s, 2, sps ->i_seq_ref_layer_chroma_phase_y_plus1);
							}
						
						bs_write_se(s, sps ->i_seq_scaled_ref_layer_left_offset);
						bs_write_se(s, sps ->i_seq_scaled_ref_layer_top_offset);
						bs_write_se(s, sps ->i_seq_scaled_ref_layer_right_offset);						
						bs_write_se(s, sps ->i_seq_scaled_ref_layer_bottom_offset );
					}
				bs_write1(s, sps ->b_seq_tcoeff_level_prediction_flag);
				if(sps ->b_seq_tcoeff_level_prediction_flag)
					bs_write1( s, sps ->b_adaptive_tcoeff_level_prediction_flag);
				bs_write1( s, sps ->b_slice_header_restriction_flag ) ;

					/*------------svc yuv扩展 一般标志位为0不进行计算-------------*/
				bs_write1(s, sps ->b_svc_vui_parameters_present_flag) ;
				if(sps ->b_svc_vui_parameters_present_flag)
					{
						//write svc_vui 参数
					}
			}

			bs_write1(s, sps ->b_additional_extension2_flag);
	  }
//扩展结束

    bs_rbsp_trailing( s );
    bs_flush( s );
}

void x264_pps_init( x264_pps_t *pps, int i_id, x264_param_t *param, x264_sps_t *sps )
{
    pps->i_id = i_id;
    pps->i_sps_id = sps->i_id;
    pps->b_cabac = param->b_cabac;

    pps->b_pic_order = !param->b_avcintra_compat && param->b_interlaced;
    pps->i_num_slice_groups = 1;

    pps->i_num_ref_idx_l0_default_active = param->i_frame_reference;
    pps->i_num_ref_idx_l1_default_active = 1;

    pps->b_weighted_pred = param->analyse.i_weighted_pred > 0;
    pps->b_weighted_bipred = param->analyse.b_weighted_bipred ? 2 : 0;

    pps->i_pic_init_qp = param->rc.i_rc_method == X264_RC_ABR || param->b_stitchable ? 26 + QP_BD_OFFSET : SPEC_QP( param->rc.i_qp_constant );
    pps->i_pic_init_qs = 26 + QP_BD_OFFSET;

    pps->i_chroma_qp_index_offset = param->analyse.i_chroma_qp_offset;
    pps->b_deblocking_filter_control = 1;
    pps->b_constrained_intra_pred = param->b_constrained_intra;
    pps->b_redundant_pic_cnt = 0;

    pps->b_transform_8x8_mode = param->analyse.b_transform_8x8 ? 1 : 0;

    pps->i_cqm_preset = param->i_cqm_preset;

    switch( pps->i_cqm_preset )
    {
    case X264_CQM_FLAT:
        for( int i = 0; i < 8; i++ )
            pps->scaling_list[i] = x264_cqm_flat16;
        break;
    case X264_CQM_JVT:
        for( int i = 0; i < 8; i++ )
            pps->scaling_list[i] = x264_cqm_jvt[i];
        break;
    case X264_CQM_CUSTOM:
        /* match the transposed DCT & zigzag */
        transpose( param->cqm_4iy, 4 );
        transpose( param->cqm_4py, 4 );
        transpose( param->cqm_4ic, 4 );
        transpose( param->cqm_4pc, 4 );
        transpose( param->cqm_8iy, 8 );
        transpose( param->cqm_8py, 8 );
        transpose( param->cqm_8ic, 8 );
        transpose( param->cqm_8pc, 8 );
        pps->scaling_list[CQM_4IY] = param->cqm_4iy;
        pps->scaling_list[CQM_4PY] = param->cqm_4py;
        pps->scaling_list[CQM_4IC] = param->cqm_4ic;
        pps->scaling_list[CQM_4PC] = param->cqm_4pc;
        pps->scaling_list[CQM_8IY+4] = param->cqm_8iy;
        pps->scaling_list[CQM_8PY+4] = param->cqm_8py;
        pps->scaling_list[CQM_8IC+4] = param->cqm_8ic;
        pps->scaling_list[CQM_8PC+4] = param->cqm_8pc;
        for( int i = 0; i < 8; i++ )
            for( int j = 0; j < (i < 4 ? 16 : 64); j++ )
                if( pps->scaling_list[i][j] == 0 )
                    pps->scaling_list[i] = x264_cqm_jvt[i];
        break;
    }
	/*sky 2014.08.28 pps init extension*/
	pps->b_base_pred_weight_table_flag = 1; // 这个值暂时赋予 1
}

void x264_pps_write( bs_t *s, x264_sps_t *sps, x264_pps_t *pps )
{
    bs_realign( s );
    bs_write_ue( s, pps->i_id );
    bs_write_ue( s, pps->i_sps_id );

    bs_write1( s, pps->b_cabac );
    bs_write1( s, pps->b_pic_order );
    bs_write_ue( s, pps->i_num_slice_groups - 1 );

    bs_write_ue( s, pps->i_num_ref_idx_l0_default_active - 1 );
    bs_write_ue( s, pps->i_num_ref_idx_l1_default_active - 1 );
    bs_write1( s, pps->b_weighted_pred );
    bs_write( s, 2, pps->b_weighted_bipred );

    bs_write_se( s, pps->i_pic_init_qp - 26 - QP_BD_OFFSET );
    bs_write_se( s, pps->i_pic_init_qs - 26 - QP_BD_OFFSET );
    bs_write_se( s, pps->i_chroma_qp_index_offset );

    bs_write1( s, pps->b_deblocking_filter_control );
    bs_write1( s, pps->b_constrained_intra_pred );
    bs_write1( s, pps->b_redundant_pic_cnt );

    if( pps->b_transform_8x8_mode || pps->i_cqm_preset != X264_CQM_FLAT )
    {
        bs_write1( s, pps->b_transform_8x8_mode );
        bs_write1( s, (pps->i_cqm_preset != X264_CQM_FLAT) );
        if( pps->i_cqm_preset != X264_CQM_FLAT )
        {
            scaling_list_write( s, pps, CQM_4IY );
            scaling_list_write( s, pps, CQM_4IC );
            bs_write1( s, 0 ); // Cr = Cb
            scaling_list_write( s, pps, CQM_4PY );
            scaling_list_write( s, pps, CQM_4PC );
            bs_write1( s, 0 ); // Cr = Cb
            if( pps->b_transform_8x8_mode )
            {
                if( sps->i_chroma_format_idc == CHROMA_444 )
                {
                    scaling_list_write( s, pps, CQM_8IY+4 );
                    scaling_list_write( s, pps, CQM_8IC+4 );
                    bs_write1( s, 0 ); // Cr = Cb
                    scaling_list_write( s, pps, CQM_8PY+4 );
                    scaling_list_write( s, pps, CQM_8PC+4 );
                    bs_write1( s, 0 ); // Cr = Cb
                }
                else
                {
                    scaling_list_write( s, pps, CQM_8IY+4 );
                    scaling_list_write( s, pps, CQM_8PY+4 );
                }
            }
        }
        bs_write_se( s, pps->i_chroma_qp_index_offset );
    }

    bs_rbsp_trailing( s );
    bs_flush( s );
}

void x264_sei_recovery_point_write( x264_t *h, bs_t *s, int recovery_frame_cnt )
{
    bs_t q;
    uint8_t tmp_buf[100];
    bs_init( &q, tmp_buf, 100 );

    bs_realign( &q );

    bs_write_ue( &q, recovery_frame_cnt ); // recovery_frame_cnt
    bs_write1( &q, 1 );   //exact_match_flag 1
    bs_write1( &q, 0 );   //broken_link_flag 0
    bs_write( &q, 2, 0 ); //changing_slice_group 0

    bs_align_10( &q );
    bs_flush( &q );

    x264_sei_write( s, tmp_buf, bs_pos( &q ) / 8, SEI_RECOVERY_POINT );
}

/*sky 2014.09.02 sei_scalability_write*/
int x264_sei_scalability_write(x264_t *h,bs_t *s)
{
	bs_realign( s);

	int sei_type = 24 ; // sei scalability info
	int sei_payload = 29 ; //经验值,sky ,要修改
	
	bs_write(s, 8, sei_type);
	bs_write(s, 8, sei_payload);

	int b_temporal_id_nesting_flag = 0;
	int b_priority_layer_info_present_flag = 0;
	int b_priority_id_setting_flag = 0;
	int i_num_layers_minus1 = h->param.i_layer_number - 1;
	
	bs_write1(s, b_temporal_id_nesting_flag);
	bs_write1(s, b_priority_layer_info_present_flag);
	bs_write1(s, b_priority_id_setting_flag);
	if(i_num_layers_minus1 + 1 <= 0 )
		return -1 ;
	bs_write_ue(s, i_num_layers_minus1);

	int b_sub_pic_layer_flag = 0;
	int b_sub_region_flag = 0;
	int b_iroi_division_info_present_flag = 0;
	int b_proflie_level_info_present_flag = 0;
	int b_bitrate_info_present_flag = 0;
	int b_frm_rate_info_present_flag = 0; // 这个可以再查
	int b_frm_size_info_present_flag = 1; 
	int b_layer_dependency_info_present_flag = 0;  // 这个flag要true
	int b_parameter_sets_info_present_flag = 1;
	int b_bitstream_restriction_info_present_flag = 0 ;
	int b_exact_interlayer_pred_flag = 1;
	int b_layer_conversion_flag = 0 ;
	int b_layer_output_flag = 1;
	
	for(int i = 0 ; i <= i_num_layers_minus1; i++)
		{
			bs_write_ue(s, i );// layer_id
			bs_write(s, 6, 0);//priority_id
			bs_write1(s, 0); // discardable_flag
			bs_write(s, 3, i);//dependency_id
			bs_write(s, 4,0);//quality_id
			bs_write(s, 3,0);//temporal_id
			
			bs_write1(s, b_sub_pic_layer_flag);
			bs_write1(s, b_sub_region_flag);
			bs_write1(s, b_iroi_division_info_present_flag);
			bs_write1(s, b_proflie_level_info_present_flag);
			bs_write1(s, b_bitrate_info_present_flag);
			bs_write1(s, b_frm_rate_info_present_flag);
			bs_write1(s, b_frm_size_info_present_flag);
			bs_write1(s, b_layer_dependency_info_present_flag);
			bs_write1(s, b_parameter_sets_info_present_flag);
			bs_write1(s, b_bitstream_restriction_info_present_flag);
			bs_write1(s, b_exact_interlayer_pred_flag);

		
		        bs_write1(s, b_layer_conversion_flag);
			bs_write1(s, b_layer_output_flag);
		
			if(b_frm_size_info_present_flag || b_iroi_division_info_present_flag)
				{
				   
					bs_write_ue(s, h->sps[i].i_mb_width - 1);
					bs_write_ue(s, h->sps[i].i_mb_height - 1);
					
				}
			if(b_layer_dependency_info_present_flag)
				{
				/*sky 这个地方还要再看*/
					bs_write_ue(s, i-1);
					bs_write_ue(s, 0);
				}
			if(b_parameter_sets_info_present_flag)
				{
					bs_write_ue(s, 0);//"ScalableSEI: num_seq_parameter_set_minus1"		
					bs_write_ue(s, 0);// "ScalableSEI: seq_parameter_set_id_delta"	
					bs_write_ue(s, 0);//"ScalableSEI: num_subset_seq_parameter_set_minus1"
					bs_write_ue(s, 0);//  "ScalableSEI: subset_seq_parameter_set_id_delta"	
					bs_write_ue(s, 0);//"ScalableSEI: num_pic_parameter_set_minus1"
					bs_write_ue(s, 0);//  "ScalableSEI: pic_parameter_set_id_delta"			
				}
		}
	
	return 0;
	
}
int x264_sei_version_write( x264_t *h, bs_t *s )
{
    // random ID number generated according to ISO-11578
    static const uint8_t uuid[16] =
    {
        0xdc, 0x45, 0xe9, 0xbd, 0xe6, 0xd9, 0x48, 0xb7,
        0x96, 0x2c, 0xd8, 0x20, 0xd9, 0x23, 0xee, 0xef
    };
    char *opts = x264_param2string( &h->param, 0 );
    char *payload;
    int length;

    if( !opts )
        return -1;
    CHECKED_MALLOC( payload, 200 + strlen( opts ) );

    memcpy( payload, uuid, 16 );
    sprintf( payload+16, "x264 - core %d%s - H.264/MPEG-4 AVC codec - "
             "Copy%s 2003-2013 - http://www.videolan.org/x264.html - options: %s",
             X264_BUILD, X264_VERSION, HAVE_GPL?"left":"right", opts );
    length = strlen(payload)+1;

    x264_sei_write( s, (uint8_t *)payload, length, SEI_USER_DATA_UNREGISTERED );

    x264_free( opts );
    x264_free( payload );
    return 0;
fail:
    x264_free( opts );
    return -1;
}

void x264_sei_buffering_period_write( x264_t *h, bs_t *s )
{
    x264_sps_t *sps = h->sps;
    bs_t q;
    uint8_t tmp_buf[100];
    bs_init( &q, tmp_buf, 100 );

    bs_realign( &q );
    bs_write_ue( &q, sps->i_id );

    if( sps->vui.b_nal_hrd_parameters_present )
    {
        bs_write( &q, sps->vui.hrd.i_initial_cpb_removal_delay_length, h->initial_cpb_removal_delay );
        bs_write( &q, sps->vui.hrd.i_initial_cpb_removal_delay_length, h->initial_cpb_removal_delay_offset );
    }

    bs_align_10( &q );
    bs_flush( &q );

    x264_sei_write( s, tmp_buf, bs_pos( &q ) / 8, SEI_BUFFERING_PERIOD );
}

void x264_sei_pic_timing_write( x264_t *h, bs_t *s )
{
    x264_sps_t *sps = h->sps;
    bs_t q;
    uint8_t tmp_buf[100];
    bs_init( &q, tmp_buf, 100 );

    bs_realign( &q );

    if( sps->vui.b_nal_hrd_parameters_present || sps->vui.b_vcl_hrd_parameters_present )
    {
        bs_write( &q, sps->vui.hrd.i_cpb_removal_delay_length, h->fenc->i_cpb_delay - h->i_cpb_delay_pir_offset );
        bs_write( &q, sps->vui.hrd.i_dpb_output_delay_length, h->fenc->i_dpb_output_delay );
    }

    if( sps->vui.b_pic_struct_present )
    {
        bs_write( &q, 4, h->fenc->i_pic_struct-1 ); // We use index 0 for "Auto"

        // These clock timestamps are not standardised so we don't set them
        // They could be time of origin, capture or alternative ideal display
        for( int i = 0; i < num_clock_ts[h->fenc->i_pic_struct]; i++ )
            bs_write1( &q, 0 ); // clock_timestamp_flag
    }

    bs_align_10( &q );
    bs_flush( &q );

    x264_sei_write( s, tmp_buf, bs_pos( &q ) / 8, SEI_PIC_TIMING );
}

void x264_sei_frame_packing_write( x264_t *h, bs_t *s )
{
    int quincunx_sampling_flag = h->param.i_frame_packing == 0;
    bs_t q;
    uint8_t tmp_buf[100];
    bs_init( &q, tmp_buf, 100 );

    bs_realign( &q );

    bs_write_ue( &q, 0 );                         // frame_packing_arrangement_id
    bs_write1( &q, 0 );                           // frame_packing_arrangement_cancel_flag
    bs_write ( &q, 7, h->param.i_frame_packing ); // frame_packing_arrangement_type
    bs_write1( &q, quincunx_sampling_flag );      // quincunx_sampling_flag

    // 0: views are unrelated, 1: left view is on the left, 2: left view is on the right
    bs_write ( &q, 6, 1 );                        // content_interpretation_type

    bs_write1( &q, 0 );                           // spatial_flipping_flag
    bs_write1( &q, 0 );                           // frame0_flipped_flag
    bs_write1( &q, 0 );                           // field_views_flag
    bs_write1( &q, h->param.i_frame_packing == 5 && !(h->fenc->i_frame&1) ); // current_frame_is_frame0_flag
    bs_write1( &q, 0 );                           // frame0_self_contained_flag
    bs_write1( &q, 0 );                           // frame1_self_contained_flag
    if ( quincunx_sampling_flag == 0 && h->param.i_frame_packing != 5 )
    {
        bs_write( &q, 4, 0 );                     // frame0_grid_position_x
        bs_write( &q, 4, 0 );                     // frame0_grid_position_y
        bs_write( &q, 4, 0 );                     // frame1_grid_position_x
        bs_write( &q, 4, 0 );                     // frame1_grid_position_y
    }
    bs_write( &q, 8, 0 );                         // frame_packing_arrangement_reserved_byte
    bs_write_ue( &q, 1 );                         // frame_packing_arrangement_repetition_period
    bs_write1( &q, 0 );                           // frame_packing_arrangement_extension_flag

    bs_align_10( &q );
    bs_flush( &q );

    x264_sei_write( s, tmp_buf, bs_pos( &q ) / 8, SEI_FRAME_PACKING );
}

void x264_filler_write( x264_t *h, bs_t *s, int filler )
{
    bs_realign( s );

    for( int i = 0; i < filler; i++ )
        bs_write( s, 8, 0xff );

    bs_rbsp_trailing( s );
    bs_flush( s );
}

void x264_sei_dec_ref_pic_marking_write( x264_t *h, bs_t *s )
{
    x264_slice_header_t *sh = &h->sh_backup;
    bs_t q;
    uint8_t tmp_buf[100];
    bs_init( &q, tmp_buf, 100 );

    bs_realign( &q );

    /* We currently only use this for repeating B-refs, as required by Blu-ray. */
    bs_write1( &q, 0 );                 //original_idr_flag
    bs_write_ue( &q, sh->i_frame_num ); //original_frame_num
    if( !h->sps->b_frame_mbs_only )
        bs_write1( &q, 0 );             //original_field_pic_flag

    bs_write1( &q, sh->i_mmco_command_count > 0 );
    if( sh->i_mmco_command_count > 0 )
    {
        for( int i = 0; i < sh->i_mmco_command_count; i++ )
        {
            bs_write_ue( &q, 1 );
            bs_write_ue( &q, sh->mmco[i].i_difference_of_pic_nums - 1 );
        }
        bs_write_ue( &q, 0 );
    }

    bs_align_10( &q );
    bs_flush( &q );

    x264_sei_write( s, tmp_buf, bs_pos( &q ) / 8, SEI_DEC_REF_PIC_MARKING );
}

int x264_sei_avcintra_umid_write( x264_t *h, bs_t *s )
{
    uint8_t data[512];
    const char *msg = "UMID";
    const int len = 497;

    memset( data, 0xff, len );
    memcpy( data, avcintra_uuid, sizeof(avcintra_uuid) );
    memcpy( data+16, msg, strlen(msg) );

    data[20] = 0x13;
    /* These bytes appear to be some sort of frame/seconds counter in certain applications,
     * but others jump around, so leave them as zero for now */
    data[21] = data[22] = 0;

    data[28] = 0x14;
    data[36] = 0x60;
    data[41] = 0x22; /* Believed to be some sort of end of basic UMID identifier */

    x264_sei_write( &h->out.bs, data, len, SEI_USER_DATA_UNREGISTERED );

    return 0;
}

int x264_sei_avcintra_vanc_write( x264_t *h, bs_t *s, int len )
{
    uint8_t data[6000];
    const char *msg = "VANC";
    if( len > sizeof(data) )
    {
        x264_log( h, X264_LOG_ERROR, "AVC-Intra SEI is too large (%d)\n", len );
        return -1;
    }

    memset( data, 0xff, len );
    memcpy( data, avcintra_uuid, sizeof(avcintra_uuid) );
    memcpy( data+16, msg, strlen(msg) );

    x264_sei_write( &h->out.bs, data, len, SEI_USER_DATA_UNREGISTERED );

    return 0;
}

const x264_level_t x264_levels[] =
{
    { 10,    1485,    99,    396,     64,    175,  64, 64,  0, 2, 0, 0, 1 },
    {  9,    1485,    99,    396,    128,    350,  64, 64,  0, 2, 0, 0, 1 }, /* "1b" */
    { 11,    3000,   396,    900,    192,    500, 128, 64,  0, 2, 0, 0, 1 },
    { 12,    6000,   396,   2376,    384,   1000, 128, 64,  0, 2, 0, 0, 1 },
    { 13,   11880,   396,   2376,    768,   2000, 128, 64,  0, 2, 0, 0, 1 },
    { 20,   11880,   396,   2376,   2000,   2000, 128, 64,  0, 2, 0, 0, 1 },
    { 21,   19800,   792,   4752,   4000,   4000, 256, 64,  0, 2, 0, 0, 0 },
    { 22,   20250,  1620,   8100,   4000,   4000, 256, 64,  0, 2, 0, 0, 0 },
    { 30,   40500,  1620,   8100,  10000,  10000, 256, 32, 22, 2, 0, 1, 0 },
    { 31,  108000,  3600,  18000,  14000,  14000, 512, 16, 60, 4, 1, 1, 0 },
    { 32,  216000,  5120,  20480,  20000,  20000, 512, 16, 60, 4, 1, 1, 0 },
    { 40,  245760,  8192,  32768,  20000,  25000, 512, 16, 60, 4, 1, 1, 0 },
    { 41,  245760,  8192,  32768,  50000,  62500, 512, 16, 24, 2, 1, 1, 0 },
    { 42,  522240,  8704,  34816,  50000,  62500, 512, 16, 24, 2, 1, 1, 1 },
    { 50,  589824, 22080, 110400, 135000, 135000, 512, 16, 24, 2, 1, 1, 1 },
    { 51,  983040, 36864, 184320, 240000, 240000, 512, 16, 24, 2, 1, 1, 1 },
    { 52, 2073600, 36864, 184320, 240000, 240000, 512, 16, 24, 2, 1, 1, 1 },
    { 0 }
};

#define ERROR(...)\
{\
    if( verbose )\
        x264_log( h, X264_LOG_WARNING, __VA_ARGS__ );\
    ret = 1;\
}

int x264_validate_levels( x264_t *h, int verbose )
{
    int ret = 0;
    int mbs = h->sps->i_mb_width * h->sps->i_mb_height;
    int dpb = mbs * h->sps->vui.i_max_dec_frame_buffering;
    int cbp_factor = h->sps->i_profile_idc>=PROFILE_HIGH422 ? 16 :
                     h->sps->i_profile_idc==PROFILE_HIGH10 ? 12 :
                     h->sps->i_profile_idc==PROFILE_HIGH ? 5 : 4;

    const x264_level_t *l = x264_levels;
    while( l->level_idc != 0 && l->level_idc != h->param.i_level_idc )
        l++;

    if( l->frame_size < mbs
        || l->frame_size*8 < h->sps->i_mb_width * h->sps->i_mb_width
        || l->frame_size*8 < h->sps->i_mb_height * h->sps->i_mb_height )
        ERROR( "frame MB size (%dx%d) > level limit (%d)\n",
               h->sps->i_mb_width, h->sps->i_mb_height, l->frame_size );
    if( dpb > l->dpb )
        ERROR( "DPB size (%d frames, %d mbs) > level limit (%d frames, %d mbs)\n",
                h->sps->vui.i_max_dec_frame_buffering, dpb, l->dpb / mbs, l->dpb );

#define CHECK( name, limit, val ) \
    if( (val) > (limit) ) \
        ERROR( name " (%"PRId64") > level limit (%d)\n", (int64_t)(val), (limit) );

    CHECK( "VBV bitrate", (l->bitrate * cbp_factor) / 4, h->param.rc.i_vbv_max_bitrate );
    CHECK( "VBV buffer", (l->cpb * cbp_factor) / 4, h->param.rc.i_vbv_buffer_size );
    CHECK( "MV range", l->mv_range, h->param.analyse.i_mv_range );
    CHECK( "interlaced", !l->frame_only, h->param.b_interlaced );
    CHECK( "fake interlaced", !l->frame_only, h->param.b_fake_interlaced );

    if( h->param.i_fps_den > 0 )
        CHECK( "MB rate", l->mbps, (int64_t)mbs * h->param.i_fps_num / h->param.i_fps_den );

    /* TODO check the rest of the limits */
    return ret;
}
