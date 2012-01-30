// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

/// \file stereo_tri_rpc.cc
///
/// Special version of Stereo Triangulation meant only for RPC models.

#include <asp/Tools/stereo.h>

int main( int argc, char* argv[] ) {

  stereo_register_sessions();
  Options opt;
  try {
    handle_arguments( argc, argv, opt );

    // Gather the RPC stereo session
    boost::shared_ptr<StereoSessionRPC> rpc_session =
      boost::dynamic_pointer_cast<StereoSessionRPC>(opt.session);

    if ( !rpc_session )
      vw_throw( ArgumentErr() << "Session doesn't appear to be RPC. Please check input.\n" );

    // Gather the Camera Models and cast them back to RPC
    boost::shared_ptr<camera::CameraModel> base_cam1, base_cam2;
    rpc_session->camera_models(base_cam1, base_cam2);
    boost::shared_ptr<asp::RPCModel> rpc1 =
      boost::dynamic_pointer_cast<RPCModel>(base_cam1);
    boost::shared_ptr<asp::RPCModel> rpc2 =
      boost::dynamic_pointer_cast<RPCModel>(base_cam2);

    if ( !rpc1 || !rpc2 )
      vw_throw( IOErr() << "Failed to recieve RPCModel from file.\n" );

  } ASP_STANDARD_CATCHES;

  return 0;
}
