// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <asp/Sessions/RPC/RPCModel.h>

#include <vw/FileIO/GdalIO.h>
#include <vw/FileIO/DiskImageResourceGDAL.h>
#include <vw/Cartography.h>

using namespace asp;
using namespace vw;

void asp::RPCModel::initialize( DiskImageResourceGDAL* resource ) {
  // Extract Datum (by means of GeoReference)
  GeoReference georef;
  read_georeference( georef, *resource );
  m_datum = georef.datum();

  // Extract RPC Info
  boost::shared_ptr<GDALDataset> dataset = resource->get_dataset_ptr();
  if ( !dataset )
    vw_throw( LogicErr() << "RPCModel: Could not read data. No file has been opened." );

  GDALRPCInfo gdal_rpc;
  if ( !GDALExtractRPCInfo( dataset->GetMetadata("RPC"),
                            &gdal_rpc) )
    vw_throw( LogicErr() << "RPCModel: GDAL resource appears not to have RPC metadata." );

  // Copy information over to our data structures.
  m_lonlatheight_offset = Vector3(gdal_rpc.dfLONG_OFF,
                                  gdal_rpc.dfLAT_OFF,
                                  gdal_rpc.dfHEIGHT_OFF);
  m_lonlatheight_scale = Vector3(gdal_rpc.dfLONG_SCALE,
                                 gdal_rpc.dfLAT_SCALE,
                                 gdal_rpc.dfHEIGHT_SCALE);
  m_xy_offset = Vector2(gdal_rpc.dfSAMP_OFF,   gdal_rpc.dfLINE_OFF);
  m_xy_scale  = Vector2(gdal_rpc.dfSAMP_SCALE, gdal_rpc.dfLINE_SCALE);

  m_line_num_coeff =   CoeffVec(gdal_rpc.adfLINE_NUM_COEFF);
  m_line_den_coeff =   CoeffVec(gdal_rpc.adfLINE_DEN_COEFF);
  m_sample_num_coeff = CoeffVec(gdal_rpc.adfSAMP_NUM_COEFF);
  m_sample_den_coeff = CoeffVec(gdal_rpc.adfSAMP_DEN_COEFF);
  m_lonlat_bbox = BBox2(Vector2(gdal_rpc.dfMIN_LONG, gdal_rpc.dfMIN_LAT),
                        Vector2(gdal_rpc.dfMAX_LONG, gdal_rpc.dfMAX_LAT));
}

asp::RPCModel::RPCModel( std::string const& filename ) {
  boost::scoped_ptr<DiskImageResourceGDAL> s_ptr( new DiskImageResourceGDAL(filename) );
  initialize( s_ptr.get() );
}

asp::RPCModel::RPCModel( DiskImageResourceGDAL* resource  ) {
  initialize( resource );
}

// All of these implementations are largely inspired by the GDAL
// code. We don't use the GDAL code unfortunately because they don't
// make that part of the API available. However I believe this is a
// safe reinterpretation that is safe to distribute.
Vector2 asp::RPCModel::point_to_pixel( Vector3 const& point ) const {

  Vector3 lonlatrad = cartography::xyz_to_lon_lat_radius( point );

  CoeffVec term =
    calculate_terms( elem_quot(lonlatrad - lonlatheight_offset,
                               lonlatheight_scale) );

  return Vector2( dot_prod(term,m_sample_num_coeff) /
                  dot_prod(term,m_sample_den_coeff) * sample_os[1] + sample_os[0],
                  dot_prod(term,m_line_num_coeff) /
                  dot_prod(term,m_line_den_coeff) * sample_os[1] + sample_os[0] );
}

Vector3 asp::RPCModel::pixel_to_vector( Vector2 const& pix ) const {
  Vector3 point, direction;
  inverse_transform( pix, point, direction );
  return direction;
}

Vector3 asp::RPCModel::camera_center( Vector2 const& pix ) const {
  Vector3 point, direction;
  inverse_transform( pix, point, direction );
  return point;
}

void asp::RPCModel::inverse_transform( Vector2 const& pix, Vector3& point, Vector3& direction ) const {
  // Step 1: Calculate an approximate starting point and assume the
  // direction vector is point towards the center of the earth.

  // Step 2: Reduce the reprojection error
}
