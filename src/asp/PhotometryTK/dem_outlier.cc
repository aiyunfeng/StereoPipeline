// __BEGIN_LICENSE__
// Copyright (C) 2006-2010 United States Government as represented by
// the Administrator of the National Aeronautics and Space Administration.
// All Rights Reserved.
// __END_LICENSE__

#include <vw/Image.h>
#include <vw/Plate/PlateFile.h>
#include <asp/Core/Macros.h>
using namespace vw;
using namespace vw::platefile;

#include <boost/foreach.hpp>
#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>
namespace po = boost::program_options;

using namespace std;

struct Options {
  // Input
  Url plate_url;
  uint64 mosaic_id;
  uint64 level;
};

void handle_arguments( int argc, char *argv[], Options& opt ) {
  po::options_description general_options("");
  general_options.add_options()
    ("mosaic_id", po::value(&opt.mosaic_id), "transaction id for where the mosaic is to compare against")
    ("level,l", po::value(&opt.level)->default_value(4), "level to calculate the error")
    ("help,h", "Display this help message");

  po::options_description positional("");
  positional.add_options()
    ("plate_url",  po::value(&opt.plate_url),  "Input PTK Url");

  po::positional_options_description positional_desc;
  positional_desc.add("plate_url", 1);

  po::options_description all_options;
  all_options.add(general_options).add(positional);

  po::variables_map vm;
  try {
    po::store( po::command_line_parser( argc, argv ).options(all_options).positional(positional_desc).run(), vm );
    po::notify( vm );
  } catch (po::error &e) {
    vw_throw( ArgumentErr() << "Error parsing input:\n\t"
              << e.what() << general_options );
  }

  std::ostringstream usage;
  usage << "Usage: " << argv[0] << " <plate-url>\n";

  if ( vm.count("help") )
    vw_throw( ArgumentErr() << usage.str() << general_options );
  if ( opt.plate_url.string().empty() )
    vw_throw( ArgumentErr() << "Missing plate file url!\n"
              << usage.str() << general_options );
}

int main( int argc, char *argv[] ) {

  Options opt;
  try {
    handle_arguments( argc, argv, opt );

    // First enforce that we have a local plate file
    if ( opt.plate_url.scheme() != "file" )
      vw_throw( ArgumentErr() << "Dem_outlier only supports local platefiles!\n");

    // Common structures
    std::map<uint64,std::string> id_lookup;
    std::map<uint64,double> id_error;
    std::map<uint64,uint64> id_count;

    // Load up the log file
    std::cout << "Reading log file: " << opt.plate_url.path()+"/plate.log" << "\n";
    std::ifstream log_file( (opt.plate_url.path()+"/plate.log").c_str() );
    std::string line;
    std::getline( log_file, line );
    while ( log_file.good() ) {
      size_t t_start = line.rfind("Transaction");
      size_t s_start = line.rfind("started");
      if ( t_start != std::string::npos &&
           s_start != std::string::npos ) {
        boost::char_separator<char> sep(" ");
        typedef boost::tokenizer<boost::char_separator<char> >
          tokenizer;
        tokenizer tokens( line, sep );
        tokenizer::const_iterator it = tokens.begin();
        for ( uint i = 0; i < 8; i++ )
          it++;
        uint64 id = boost::lexical_cast<uint64>(*it);
        for ( uint i = 0; i < 2; i++ )
          it++;
        if ( id != opt.mosaic_id )
          id_lookup[id] = *it;
      }

      std::getline( log_file, line );
    }
    log_file.close();

    // Print out matches
    std::cout << "\nResults:\n";
    for ( std::map<uint64,std::string>::const_iterator it = id_lookup.begin();
          it != id_lookup.end(); it++ ) {
      std::cout << it->first << "\t" << it->second << "\n";
    }
    std::cout << "\n";

    // Finally load up the plate file and work
    PlateFile platefile( opt.plate_url );
    {
      TerminalProgressCallback tpc(InfoMessage,"Error: ");
      tpc.report_progress(0);
      double inc_amt = 1.0 / double(id_lookup.size());
      BBox2i region( 0, 0, (0x1<<opt.level)-1, (0x1<<opt.level)-1 );
      for ( std::map<uint64,std::string>::const_iterator it = id_lookup.begin();
            it != id_lookup.end(); it++ ) {
        std::list<TileHeader> tiles =
          platefile.search_by_region( opt.level, region, it->first,
                                      it->first+1, 0, true );
        id_error[it->first] = 0;
        id_count[it->first] = 0;

        BOOST_FOREACH( TileHeader tile, tiles ) {
          ImageView<PixelGrayA<float> > subject, reference;
          platefile.read( subject, tile.col(), tile.row(), opt.level,
                          tile.transaction_id(), true );
          platefile.read( reference, tile.col(), tile.row(), opt.level,
                          opt.mosaic_id, true );
          id_error[it->first] += sum_of_channel_values(abs(apply_mask(alpha_to_mask(subject) - alpha_to_mask(reference),0)));
          id_count[it->first] += sum_of_channel_values(select_channel(subject,1));
        }
        tpc.report_incremental_progress( inc_amt );
      }
      tpc.report_finished();
    }

    // Normalizing error
    double mean, stddev;
    {
      TerminalProgressCallback tpc(InfoMessage,"Normalize: ");
      tpc.report_progress(0);
      double inc_amt = 1.0 / double(id_lookup.size());
      MeanAccumulator<double> mean_acc;
      StdDevAccumulator<double> stddev_acc;
      for ( std::map<uint64,std::string>::const_iterator it = id_lookup.begin();
            it != id_lookup.end(); it++ ) {
        id_error[it->first] /= double( id_count[it->first] );
        mean_acc( id_error[it->first] );
        stddev_acc( id_error[it->first] );
        tpc.report_incremental_progress( inc_amt );
      }
      mean = mean_acc.value();
      stddev = stddev_acc.value();
      tpc.report_finished();
    }
    std::cout << "\nMean error:   " << mean << " meters.\n";
    std::cout << "StdDev error: " << stddev << " meters.\n\n";
    std::cout << "Bad ------------------\n";

    

    // Print out who fail me rigorous standards
    for ( std::map<uint64,std::string>::const_iterator it = id_lookup.begin();
          it != id_lookup.end(); it++ ) {
      if ( fabs(id_error[it->first] - mean) >= stddev ) {
        std::cout << it->second << " : [" << fabs(id_error[it->first]-mean) << "]\n";
      }

      // At this point I should search again to find what tiles this image is at.

      // Then see what other images share that same tile

      // Write all combination with this id and ids that are higher. (This avoids repeats)
    }

  } ASP_STANDARD_CATCHES;

  return 0;
}
