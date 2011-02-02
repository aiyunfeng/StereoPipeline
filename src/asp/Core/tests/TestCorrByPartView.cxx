// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__


#include <gtest/gtest.h>

#include <vw/Image/UtilityViews.h>
#include <vw/Image/Transform.h>
#include <asp/Core/CorrByPartView.h>

#include <boost/random/linear_congruential.hpp>

using namespace vw;
using namespace vw::stereo;

class BasicCorrelationTest : public ::testing::Test {
protected:
  BasicCorrelationTest() {}

  virtual void SetUp() {
    boost::rand48 gen(10);
    image1 = 255*uniform_noise_view( gen, 50, 50 );
    image2 = transform(image1, TranslateTransform(3,3),
                       ZeroEdgeExtension(), NearestPixelInterpolation());
    mask.set_size(50,50);
    fill(mask,PixelMask<uint8>(255));
    search.set_size(2,2);
    fill(search,Vector4(0,0,6,6));
  }

  template <class ViewT, class MViewT, class SViewT, class PreProcT>
  CorrByPartView<typename ViewT::pixel_type,
                 typename MViewT::pixel_type,
                 PreProcT>
  correlate( ImageViewBase<ViewT> const& input1,
             ImageViewBase<ViewT> const& input2,
             ImageViewBase<MViewT> const& mask,
             ImageViewBase<SViewT> const& search,
             PreProcT const& proc) {
    typedef typename ViewT::pixel_type pixel_type;
    typedef typename MViewT::pixel_type mask_type;
    CorrByPartView<pixel_type,mask_type,PreProcT> corr( input1, input2, mask, mask,
                                                        search, 25, proc );
    corr.set_kernel_size( 7 );
    return corr;
  }

  template <class ViewT>
  void check_error( ImageViewBase<ViewT> const& input,
                    float success = 0.9 ) {
    ViewT const& disparity_map = input.impl();
    int count_correct = 0;
    int count_valid = 0;
    for (int j = 0; j < disparity_map.rows(); ++j)
      for (int i = 0; i < disparity_map.cols(); ++i)
        if ( is_valid( disparity_map(i,j) ) ) {
          count_valid++;
          if ( disparity_map(i,j).child() == Vector2f(3,3) )
            count_correct++;
        }
    EXPECT_LT( success, float(count_correct)/float(count_valid) );
  }

  ImageView<uint8> image1, image2;
  ImageView<PixelMask<uint8> > mask;
  ImageView<Vector4> search;
};

TEST_F( BasicCorrelationTest, NullPreprocess ) {
  ImageView<PixelMask<Vector2f> > disparity_map =
    correlate( image1, image2, mask, search,
               NullStereoPreprocessingFilter() );
  check_error( disparity_map, 0.95 );
}

TEST_F( BasicCorrelationTest, SlogPreprocess ) {
  ImageView<PixelMask<Vector2f> > disparity_map =
    correlate( image1, image2, mask, search,
               SlogStereoPreprocessingFilter() );
  check_error( disparity_map, 0.80 );
}

TEST_F( BasicCorrelationTest, LogPreprocess ) {
  ImageView<PixelMask<Vector2f> > disparity_map =
    correlate( image1, image2, mask, search,
               LogStereoPreprocessingFilter() );
  check_error( disparity_map, 0.80 );
}

TEST_F( BasicCorrelationTest, BlurPreprocess ) {
  ImageView<PixelMask<Vector2f> > disparity_map =
    correlate( image1, image2, mask, search,
               BlurStereoPreprocessingFilter() );
  check_error( disparity_map, 0.75 );
}
