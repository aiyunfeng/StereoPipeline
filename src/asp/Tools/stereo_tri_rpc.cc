// __BEGIN_LICENSE__
// Copyright (C) 2006-2011 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

/// \file stereo_tri_rpc.cc
///
/// Special version of Stereo Triangulation meant only for RPC models.

#include <asp/Tools/stereo.h>
#include <asp/Sessions/RPC/RPCModel.h>

using namespace vw;

int main( int argc, char* argv[] ) {

  stereo_register_sessions();
  Options opt;
  try {
    handle_arguments( argc, argv, opt );

    // Gather the RPC stereo session
    boost::shared_ptr<asp::StereoSessionRPC> rpc_session =
      boost::dynamic_pointer_cast<asp::StereoSessionRPC>(opt.session);

    if ( !rpc_session )
      vw_throw( ArgumentErr() << "Session doesn't appear to be RPC. Please check input.\n" );

    // Gather the Camera Models and cast them back to RPC
    boost::shared_ptr<camera::CameraModel> base_cam1, base_cam2;
    rpc_session->camera_models(base_cam1, base_cam2);
    boost::shared_ptr<asp::RPCModel> rpc1 =
      boost::dynamic_pointer_cast<asp::RPCModel>(base_cam1);
    boost::shared_ptr<asp::RPCModel> rpc2 =
      boost::dynamic_pointer_cast<asp::RPCModel>(base_cam2);

    if ( !rpc1 || !rpc2 )
      vw_throw( IOErr() << "Failed to recieve RPCModel from file.\n" );

    // I'm still not ready to perform triangulation yet. Here's just
    // some printing to prove that I actually did recieve an RPC
    // model.

    // Triangulation with RPC is a little tricky. Here's the plan.
    // 1.) Solve for the mean XYZ for both cameras by looking at the
    // LLH boxes.
    // 2.) Solve for one ray's intersection at the mean height. This
    // will seed the rigoruos algorithm. I believe this is equation 24
    // in Grodecki 2004. We only need the first iteration's result.
    // 3.) Implement the rigourous triangulation algorithm in LMA with
    // analytical jacobians.

  } ASP_STANDARD_CATCHES;

  return 0;
}
