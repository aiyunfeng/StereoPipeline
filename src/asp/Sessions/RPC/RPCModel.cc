// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <asp/Sessions/RPC/RPCModel.h>
#include <gdal.h>

using namespace asp;
using namespace vw;

RPCModel::RPCModel( DiskImageResourceGDAL const& resource  ) {

  boost::shared_ptr<GDALDataset> dataset = resource.get_dataset_ptr();
  if ( !dataset )
    vw_throw( LogicErr() << "RPCModel: Could not read data. No file has been opened." );

  GDALRPCInfo gdal_rpc;
  if ( !GDALExtractRPCInfo( dataset->GetMetadata("RPC"),
                            &gdal_rpc) )
    vw_throw( LogicErr() << "RPCModel: GDAL resource appears not to have RPC metadata." );

  // Copy information over to our data structures.
  lonlatheight_offset = Vector3(gdal_rpc.dfLONG_OFF,
                                gdal_rpc.dfLAT_OFF,
                                gdal_rpc.dfHEIGHT_OFF);
  lonlatheight_scale = Vector3(gdal_rpc.dfLONG_SCALE,
                               gdal_rpc.dfLAT_SCALE,
                               gdal_rpc.dfHEIGHT_SCALE);
  line_os = Vector2(gdal_rpc.dfLINE_OFF, gdal_rpc.dfLINE_SCALE);
  sample_os = Vector2(gdal_rpc.dfSAMP_OFF, gdal_rpc.dfSAMP_SCALE);

  m_line_num_coeff = Vector<double,20>(gdal_rpc.adfLINE_NUM_COEFF);
  m_line_den_coeff = Vector<double,20>(gdal_rpc.adfLINE_DEN_COEFF);
  m_sample_num_coeff = Vector<double,20>(gdal_rpc.adfSAMP_NUM_COEFF);
  m_sample_den_coeff = Vector<double,20>(gdal_rpc.adfSAMP_DEN_COEFF);
  m_lonlat_bbox = BBox2(Vector2(gdal_rpc.dfMIN_LONG, gdal_rpc.dfMIN_LAT),
                        Vector2(gdal_rpc.dfMAX_LONG, gdal_rpc.dfMAX_LAT));
}

Vector<double,20> RPCModel::calculate_terms( Vector3 const& v ) {
  Vector<double,20> result;
  result[0] = 1.0;
  result[1] = v.x();
  result[2] = v.y();
  result[3] = v.z();
  result[4] = v.x() * v.y();
  result[5] = v.x() * v.z();
  result[6] = v.y() * v.z();
  result[7] = v.x() * v.x();
  result[8] = v.y() * v.y();
  result[9] = v.z() * v.z();

  result[10] = v.x() * v.y() * v.z();
  result[11] = v.x() * v.x() * v.x();
  result[12] = v.x() * v.y() * v.y();
  result[13] = v.x() * v.z() * v.z();
  result[14] = v.x() * v.x() * v.y();
  result[15] = v.y() * v.y() * v.y();
  result[16] = v.y() * v.z() * v.z();
  result[17] = v.x() * v.x() * v.z();
  result[18] = v.y() * v.y() * v.z();
  result[19] = v.z() * v.z() * v.z();
  return result;
}

// All of these implementations are largely inspired by the GDAL
// code. We don't use the GDAL code unfortunately because they don't
// make that part of the API available. However I believe this is a
// safe reinterpretation that is safe to distribute.
Vector2 RPCModel::point_to_pixel( Vector3 const& point ) const {

  Vector3 lonlatrad = xyz_to_lon_lat_radius( point );

  Vector<double,20> term =
    calculate_terms( elem_quot(lonlatrad - lonlatheight_offset,
                               lonlatheight_scale) );

  return Vector2( dot_prod(term,m_sample_num_coeff) /
                  dot_prod(term,m_sample_den_coeff) * sample_os[1] + sample_os[0],
                  dot_prod(term,m_line_num_coeff) /
                  dot_prod(term,m_line_den_coeff) * sample_os[1] + sample_os[0] );
}

Vector3 RPCModel::pixel_to_vector( Vector2 const& pix ) const {
  Vector3 point, direction;
  inverse_transform( pix, point, direction );
  return direction;
}

Vector3 RPCModel::camera_center( Vector2 const& pix ) const {
  Vector3 point, direction;
  inverse_transform( pix, point, direction );
  return point;
}

void RPCModel::inverse_transform( Vector2 const& pix, Vector3& point, Vector3& direction ) {
  // Step 1: Calculate an approximate starting point and assume the
  // direction vector is point towards the center of the earth.

  // Step 2: Reduce the reprojection error
}
