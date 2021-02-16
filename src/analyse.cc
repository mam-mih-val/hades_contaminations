#include <iostream>
#include <chrono>
#include <boost/program_options.hpp>

#include <AnalysisTree/TaskManager.hpp>
#include "analysis_task.h"

int main(int n_args, char** args){
  namespace po=boost::program_options;
  if(n_args<2){
    throw std::runtime_error( "Please type \"./acceptance --help\" to get information" );
  }
  std::string file_list;
  std::string output_file{"output.root"};
  std::string efficiency_file{"output.root"};
  int physical_trgger{0};
  int n_events=-1;
  po::options_description options("Options");
  options.add_options()
      ("help,h", "Help screen")
      ("input,i", po::value<std::string>(&file_list),
       "Path to input file list")
      ("output,o", po::value<std::string>(&output_file),
       "output file name")
      ("efficiency,e", po::value<std::string>(&efficiency_file),
       "Path to file with protons efficiency")
      ("physical_trigger,p", po::value<int>(&physical_trgger),
       "Physical trigger number (2 or 3)")
      ("n-events,N", po::value<int>(&n_events),
       "Number of events to process (-1=all)")
      ("start-collisions,s","Selects collisions in START detector");
  po::variables_map vm;
  po::parsed_options parsed = po::command_line_parser(n_args, args).options(options).run();
  po::store(parsed, vm);
  po::notify(vm);
  if (vm.count("help")){
    std::cout << options << std::endl;
    return 0;
  }
  if( physical_trgger != 2 && physical_trgger != 3 && physical_trgger != 0 )
    throw std::runtime_error( R"(Error in physical trigger set value. Only "2" or "3" values are expected)" );
  auto is_in_start = vm.count("start-collisions");
  AnalysisTree::Cuts* evet_cuts;
  std::vector<AnalysisTree::SimpleCut> vector_of_cuts;
  vector_of_cuts.emplace_back(AnalysisTree::SimpleCut({"event_header", "selected_mdc_tracks"}, 2.0, 999.0));
  if( physical_trgger != 0 )
    vector_of_cuts.emplace_back(AnalysisTree::SimpleCut{ {"event_header", "physical_trigger_"+std::to_string(physical_trgger)}, 1 });
  if( !is_in_start )
    vector_of_cuts.emplace_back(AnalysisTree::SimpleCut{ {"event_header", "vtx_z"}, -70.0, -5.0 });
  else
    vector_of_cuts.emplace_back(AnalysisTree::SimpleCut{ {"event_header", "vtx_z"}, -90.0, -70.0 });

  evet_cuts = new AnalysisTree::Cuts( "selected_events", vector_of_cuts );

  AnalysisTree::TaskManager manager({file_list}, {"hades_analysis_tree"});
  manager.SetEventCuts(evet_cuts);
  auto *analysis_task = new AnalysisTree::AnalysisTask;
  analysis_task->InitEffieciencies(efficiency_file);
  manager.AddTask(analysis_task);
  manager.SetOutFileName(output_file);
  manager.Init();
  manager.Run(n_events);
  manager.Finish();
  return 0;
}