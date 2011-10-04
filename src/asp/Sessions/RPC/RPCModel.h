// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#ifndef __STEREO_SESSION_RPC_MODEL_H__
#define __STEREO_SESSION_RPC_MODEL_H__

#include <vw/FileIO/DiskImageResourceGDAL.h>
#include <vw/Camera/CameraModel.h>

namespace asp {

  // TODO: Need to work out a different way to triangulate. Our
  //   standard midpoint method doesn't seem to work.

  class RPCModel : public vw::camera::CameraModel {
    using namespace vw;
    Vector<double,20> m_line_num_coeff, m_line_den_coeff,
      m_sample_num_coeff, m_sample_den_coeff;
    BBox2 m_lonlat_bbox;

    Vector<double,20> calculate_terms( Vector3 const& v );
  public:
    RPCModel( std::string const& filename ) : RPCModel( vw::DiskImageResourceGDAL(filename) ) {}
    RPCModel( vw::DiskImageResourceGDAL const& resource );

    virtual std::string type() const { return "RPC"; }
    virtual ~RPCModel() {}

    // Standard Access Methods
    virtual Vector2 point_to_pixel( Vector3 const& point ) const;
    virtual Vector3 pixel_to_vector( Vector2 const& pix ) const;
    virtual Vector3 camera_center( Vector2 const& pix ) const;

    // Non standard access, but more efficient
    void inverse_transform( Vector2 const& pix, Vector3& point, Vector3& direction );

    // Scaling parameters
    Vector2 line_os; // [OFFSET, SCALE]
    Vector2 sample_os;
    Vector3 lonlatheight_offset;
    Vector3 lonlatheight_scale;
  };

}

#endif//__STEREO_SESSION_RPC_CAMERA_MODEL_H__
