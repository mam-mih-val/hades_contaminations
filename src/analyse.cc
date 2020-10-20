#include <iostream>
#include <chrono>

#include <AnalysisTree/TaskManager.hpp>

#include "analysis_task.h"

int main(int n_args, char** args){
  if(n_args<2){
    throw std::runtime_error( "Please use: ./acceptance file.list" );
  }
  std::string list{args[1]};
  AnalysisTree::TaskManager manager({list}, {"hades_analysis_tree"});

  auto *analysis_task = new AnalysisTree::AnalysisTask;

  manager.AddTask(analysis_task);
  manager.SetOutFileName("output.root");
  manager.Init();
  manager.Run(-1);
  manager.Finish();
  return 0;
}