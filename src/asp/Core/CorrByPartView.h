// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __ASP_CORR_BY_PART_VIEW__
#define __ASP_CORR_BY_PART_VIEW__

#include <vw/Image/ImageViewRef.h>
#include <vw/Stereo/Correlate.h>
#include <vw/Stereo/DisparityMap.h>

template <class ImagePixelT, class MaskPixelT, class PreProcFuncT>
class CorrByPartView : public vw::ImageViewBase<CorrByPartView<ImagePixelT, MaskPixelT, PreProcFuncT> > {

  // Inputs
  vw::ImageViewRef<ImagePixelT> m_left_image, m_right_image;
  vw::ImageViewRef<MaskPixelT> m_left_mask, m_right_mask;
  vw::ImageView<vw::Vector4> m_search_image;
  PreProcFuncT m_preproc_func;
  size_t m_partition_size;

  // Settings
  vw::Vector2i m_kernel_size, m_kernel_pad;
  float m_cross_corr_threshold, m_corr_score_threshold;
  size_t m_cost_blur;
  vw::stereo::CorrelatorType m_correlator_type;

public:
  typedef vw::PixelMask<Vector2f> pixel_type;
  typedef pixel_type result_type;
  typedef vw::ProceduralPixelAccessor<CorrByPartView> pixel_accessor;

  template <class ImageT1, class ImageT2, class MaskT1, class MaskT2, class SearchT>
  CorrByPartView( vw::ImageViewBase<ImageT1> const& left_image,
                  vw::ImageViewBase<ImageT2> const& right_image,
                  vw::ImageViewBase<MaskT1> const& left_mask,
                  vw::ImageViewBase<MaskT2> const& right_mask,
                  vw::ImageViewBase<SearchT> const& search_image,
                  size_t partition_size,
                  PreProcFuncT const& preproc_func ) :
    m_left_image(left_image.impl()), m_right_image(right_image.impl()),
    m_left_mask(left_mask.impl()), m_right_mask(right_mask.impl()),
    m_search_image(search_image.impl()), m_partition_size(partition_size),
    m_preproc_func(preproc_func) {

    // Assertions for operation
    VW_ASSERT((m_left_image.cols() == m_right_image.cols()) &&
              (m_left_image.rows() == m_right_image.rows()),
              ArgumentErr() << "CorrByPartView::CorrByPartView(): input image dimensions do not agree.\n");

    VW_ASSERT((m_left_image.cols() == m_left_mask.cols()) &&
              (m_left_image.rows() == m_left_mask.rows()),
              ArgumentErr() << "CorrByPartView::CorrByPartView(): input image and mask image dimensions do not agree.\n");

    VW_ASSERT((m_left_image.cols() == m_right_mask.cols()) &&
              (m_left_image.rows() == m_right_mask.rows()),
              ArgumentErr() << "CorrByPartView::CorrByPartView(): input image and mask image dimensions do not agree.\n");

    VW_ASSERT((left_image.channels() == 1) && (left_image.impl().planes() == 1) &&
              (right_image.channels() == 1) && (right_image.impl().planes() == 1),
              ArgumentErr() << "CorrByPartView::CorrByPartView(): multi-channel, multi-plane images not supported.\n");

    VW_ASSERT(search_image.channels() == 4 &&
              search_image.planes() == 1,
              ArgumentErr() << "CorrByPartView::CorrByPartView(): Search image is not single plane, Vector4 image.\n");

    VW_ASSERT(m_search_image.cols()*m_partition_size >= m_left_image.cols() &&
              m_search_image.rows()*m_partition_size >= m_left_image.rows(),
              ArgumentErr() << "CorrByPartView::CorrByPartView(): Search image doesn't cover the entire space of input.\n");

    // Set some defaults
    m_kernel_size = m_kernel_pad = vw::Vector2i(1,1);
    m_cross_corr_threshold = 2.0;
    m_corr_score_threshold = 1.3;
    m_cost_blur = 1;
    m_correlator_type = ABS_DIFF_CORRELATOR;
  }

  // Access to settings
  void set_kernel_size(Vector2i size) {
    m_kernel_size = size;
    m_kernel_pad = size/2;
  }
  vw::Vector2i kernel_size() const { return m_kernel_size; }

  void set_correlator_options(size_t cost_blur, stereo::CorrelatorType correlator_type) {
    m_cost_blur = cost_blur;
    m_correlator_type = correlator_type;
  }
  size_t cost_blur() const { return m_cost_blur; }
  stereo::CorrelatorType correlator_type() const { return m_correlator_type; }

  void set_cross_corr_threshold(float threshold) { m_cross_corr_threshold = threshold; }
  float cross_corr_threshold() const { return m_cross_corr_threshold; }

  void set_corr_score_threshold(float threshold) { m_corr_score_threshold = threshold; }
  float corr_score_threshold() const { return m_corr_score_threshold; }

  // Standard ImageView interface methods
  inline int32 cols() const { return m_left_image.cols(); }
  inline int32 rows() const { return m_left_image.rows(); }
  inline int32 planes() const { return 1; }

  inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

  inline pixel_type operator()(vw::int32 /*i*/, vw::int32 /*j*/, vw::int32 /*p*/ = 0) const {
    vw_throw(NoImplErr() << "CorrByPartView::operator()(double i, double j, int32 p) has not been implemented.");
    return pixel_type();
  }

  typedef CropView<ImageView<pixel_type> > prerasterize_type;
  inline prerasterize_type prerasterize(BBox2i bbox) const {
    vw_out(DebugMessage, "asp") << "CorrelatorView: rasterizing image block " << bbox << ".\n";

    // Determine pixels that we touch in search image
  }

  template <class DestT>
  inline void rasterize(DestT const& dest, BBox2i bbox) const {
    vw::rasterize(prerasterize(bbox), dest, bbox);
  }
};

#endif//__ASP_CORR_BY_PART_VIEW__
