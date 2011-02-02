// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#ifndef __ASP_CORR_BY_PART_VIEW__
#define __ASP_CORR_BY_PART_VIEW__

#include <vw/Image/Filter.h>
#include <vw/Image/ImageViewRef.h>
#include <vw/Stereo/OptimizedCorrelator.h>
#include <vw/Stereo/DisparityMap.h>

template <class ImagePixelT, class MaskPixelT, class PreProcFuncT>
class CorrByPartView : public vw::ImageViewBase<CorrByPartView<ImagePixelT, MaskPixelT, PreProcFuncT> > {

  // Inputs
  vw::ImageViewRef<ImagePixelT> m_left_image, m_right_image;
  vw::ImageViewRef<MaskPixelT> m_left_mask, m_right_mask;
  vw::ImageView<vw::Vector4f> m_search_image;
  size_t m_partition_size;
  PreProcFuncT m_preproc_func;

  // Settings
  vw::Vector2i m_kernel_size, m_kernel_pad;
  float m_cross_corr_threshold, m_corr_score_threshold;
  size_t m_cost_blur;
  vw::stereo::CorrelatorType m_correlator_type;

public:
  typedef vw::PixelMask<vw::Vector2f> pixel_type;
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
              vw::ArgumentErr() << "CorrByPartView::CorrByPartView(): input image dimensions do not agree.\n");

    VW_ASSERT((m_left_image.cols() == m_left_mask.cols()) &&
              (m_left_image.rows() == m_left_mask.rows()),
              vw::ArgumentErr() << "CorrByPartView::CorrByPartView(): input image and mask image dimensions do not agree.\n");

    VW_ASSERT((m_left_image.cols() == m_right_mask.cols()) &&
              (m_left_image.rows() == m_right_mask.rows()),
              vw::ArgumentErr() << "CorrByPartView::CorrByPartView(): input image and mask image dimensions do not agree.\n");

    VW_ASSERT((left_image.channels() == 1) && (left_image.impl().planes() == 1) &&
              (right_image.channels() == 1) && (right_image.impl().planes() == 1),
              vw::ArgumentErr() << "CorrByPartView::CorrByPartView(): multi-channel, multi-plane images not supported.\n");

    VW_ASSERT(search_image.channels() == 4 &&
              search_image.impl().planes() == 1,
              vw::ArgumentErr() << "CorrByPartView::CorrByPartView(): Search image is not single plane, Vector4 image.\n");

    VW_ASSERT(m_search_image.cols()*vw::int32(m_partition_size) >= m_left_image.cols() &&
              m_search_image.rows()*vw::int32(m_partition_size) >= m_left_image.rows(),
              vw::ArgumentErr() << "CorrByPartView::CorrByPartView(): Search image doesn't cover the entire space of input.\n");

    // Set some defaults
    m_kernel_size = m_kernel_pad = vw::Vector2i(1,1);
    m_cross_corr_threshold = 2.0;
    m_corr_score_threshold = 1.3;
    m_cost_blur = 1;
    m_correlator_type = vw::stereo::ABS_DIFF_CORRELATOR;
  }

  // Access to settings
  void set_kernel_size(ssize_t size) {
    m_kernel_size = vw::Vector2i(size,size);
    m_kernel_pad = m_kernel_size/2;
  }
  ssize_t kernel_size() const { return m_kernel_size[0]; }

  void set_correlator_options(size_t cost_blur, vw::stereo::CorrelatorType correlator_type) {
    m_cost_blur = cost_blur;
    m_correlator_type = correlator_type;
  }
  size_t cost_blur() const { return m_cost_blur; }
  vw::stereo::CorrelatorType correlator_type() const { return m_correlator_type; }

  void set_cross_corr_threshold(float threshold) { m_cross_corr_threshold = threshold; }
  float cross_corr_threshold() const { return m_cross_corr_threshold; }

  void set_corr_score_threshold(float threshold) { m_corr_score_threshold = threshold; }
  float corr_score_threshold() const { return m_corr_score_threshold; }

  // Standard ImageView interface methods
  inline vw::int32 cols() const { return m_left_image.cols(); }
  inline vw::int32 rows() const { return m_left_image.rows(); }
  inline vw::int32 planes() const { return 1; }

  inline pixel_accessor origin() const { return pixel_accessor( *this, 0, 0 ); }

  inline pixel_type operator()(vw::int32 /*i*/, vw::int32 /*j*/, vw::int32 /*p*/ = 0) const {
    vw_throw(vw::NoImplErr() << "CorrByPartView::operator()(double i, double j, int32 p) has not been implemented.");
    return pixel_type();
  }

  typedef vw::ImageView<pixel_type> prerasterize_type;
  inline prerasterize_type prerasterize(vw::BBox2i bbox) const {
    using vw::ImageView;
    using vw::BBox2i;
    using vw::Vector4f;
    using vw::Vector2i;

    vw_out(vw::DebugMessage, "asp") << "CorrByPartView: rasterizing image block " << bbox << ".\n";

    ImageView<pixel_type> result( bbox.width(), bbox.height() );

    // Determine pixels that we touch in search image
    BBox2i search_image_bbox = bbox / m_partition_size;
    vw_out(vw::DebugMessage, "asp") << "CorrByPartView: search block " << search_image_bbox << ".\n";

    // Rasterizing into result
    for ( ssize_t sj = search_image_bbox.min()[1];
          sj < search_image_bbox.max()[1]; sj++ ) {
      for ( ssize_t si = search_image_bbox.min()[0];
            si < search_image_bbox.max()[0]; si++ ) {
        BBox2i left_crop_bbox( si * m_partition_size,
                                sj * m_partition_size,
                                m_partition_size, m_partition_size );
        BBox2i destination = left_crop_bbox;
        Vector4f srange( m_search_image(si,sj) );
        BBox2i right_crop_bbox( left_crop_bbox.min() + subvector(srange,0,2),
                                left_crop_bbox.max() + subvector(srange,2,2) );

        // The correlator requires the images to be the same size. The
        // search bbox will always be larger than the given left image
        // bbox, so we just make the left bbox the same size as the
        // right bbox.
        left_crop_bbox.max() = left_crop_bbox.min() +
          Vector2i(right_crop_bbox.width(), right_crop_bbox.height());

        // Finally, we must adjust both bounding boxes to account for
        // the size of the kernel itself.
        right_crop_bbox.min() -= m_kernel_pad;
        right_crop_bbox.max() += m_kernel_pad;
        left_crop_bbox.min() -= m_kernel_pad;
        left_crop_bbox.max() += m_kernel_pad;

        vw::stereo::OptimizedCorrelator
          correlator(BBox2i(0,0,srange[2]-srange[0],
                            srange[3]-srange[1]),
                     m_kernel_size[0],
                     m_cross_corr_threshold,
                     m_corr_score_threshold,
                     m_cost_blur, m_correlator_type );
        ImageView<pixel_type> disparity_map =
          disparity_mask( correlator( crop( edge_extend(m_left_image), left_crop_bbox ),
                                      crop( edge_extend(m_right_image), right_crop_bbox ),
                                      m_preproc_func ),
                          crop( edge_extend(m_left_image), left_crop_bbox ),
                          crop( edge_extend(m_right_image), right_crop_bbox ) );

        // Adjust the disparities to be relative to the uncropped
        // pixel location
        // TODO: Make pixel functor, would this be faster?
        for (ssize_t v = 0; v < disparity_map.rows(); ++v)
          for (ssize_t u = 0; u < disparity_map.cols(); ++u)
            if (is_valid(disparity_map(u,v)) )
              remove_mask(disparity_map(u,v)) += subvector(srange,0,2);

        // Assigning back into result
        destination.crop(bbox);
        destination -= bbox.min();
        crop( result, destination ) =
          crop( disparity_map, BBox2i( m_kernel_pad[0], m_kernel_pad[1],
                                       m_partition_size, m_partition_size ) );
      }
    }

    return result;
  }

  template <class DestT>
  inline void rasterize(DestT const& dest, vw::BBox2i bbox) const {
    vw::rasterize(prerasterize(bbox), dest, bbox);
  }
};

template <class ViewT>
vw::ImageView<vw::Vector4f>
compute_search_ranges(vw::ImageViewBase<ViewT> const& map,
                      size_t partition_size ) {
  vw::ImageView<vw::Vector4f> search_range( map.cols()/partition_size,
                                           map.rows()/partition_size);
  for ( ssize_t j = 0; j < search_range.row(); j++ ) {
    for ( ssize_t i = 0; i < search_range.col(); i++ ) {
      ImageViewRef<PixelMask<Vector2f> > crop_disparity =
        crop(map.impl(),BBox2i(i*partition_size,j*partition_size,
                               partition_size,partition_size));
      // TODO: Count valid and get range should happen in the same run
      if ( count_valid_pixels(crop_disparity) > 20*20 ) {
        BBox2f s = get_disparity_range(crop_disparity);
        subvector(search_range(i,j),0,2) = s.min();
        subvector(search_range(i,j),2,2) = s.max();
      } else {
        search_range(i,j) = Vector4f(0,0,0,0);
      }
    }
  }

  // Could grow
}

#endif//__ASP_CORR_BY_PART_VIEW__
