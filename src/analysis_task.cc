//
// Created by mikhail on 6/16/20.
//

#include "analysis_task.h"

namespace AnalysisTree {
void AnalysisTask::Init(std::map<std::string, void *> &branch_map) {
  // linking pointers with branch fields
  event_header_ = static_cast<EventHeader *>(branch_map.at("event_header"));
  mdc_vtx_tracks_ = static_cast<Particles *>(branch_map.at("mdc_vtx_tracks"));
  meta_hits_ = static_cast<HitDetector *>(branch_map.at("meta_hits"));
  wall_hits_ = static_cast<HitDetector *>(branch_map.at("forward_wall_hits"));
  mdc_meta_matching_ = static_cast<Matching *>(branch_map.at("mdc_vtx_tracks2meta_hits"));

  // getting branch configurations, which store information about fields in branches
  auto event_header_config = config_->GetBranchConfig("event_header");
  auto mdc_vtx_tracks_config = config_->GetBranchConfig("mdc_vtx_tracks");
  auto meta_hits_config = config_->GetBranchConfig("meta_hits");
  auto wall_hits_config = config_->GetBranchConfig("forward_wall_hits");

  // linking necessary for analysis fields with enumerator for fast access to them
  fields_id_.insert(std::make_pair(FIELDS::HITS_TOF, event_header_config.GetFieldId("selected_tof_hits")));
  fields_id_.insert(std::make_pair(FIELDS::HITS_RPC, event_header_config.GetFieldId("selected_rpc_hits")));
  fields_id_.insert(std::make_pair(FIELDS::TRACKS_MDC, event_header_config.GetFieldId("selected_mdc_tracks")));
  fields_id_.insert(std::make_pair(FIELDS::FW_SIGNAL, event_header_config.GetFieldId("fw_adc")));
  fields_id_.insert(std::make_pair(FIELDS::PT3, event_header_config.GetFieldId("physical_trigger_3")));
  fields_id_.insert(std::make_pair(FIELDS::PT2, event_header_config.GetFieldId("physical_trigger_2")));
  fields_id_.insert(std::make_pair(FIELDS::CHI_2, mdc_vtx_tracks_config.GetFieldId("chi2")));
  fields_id_.insert(std::make_pair(FIELDS::DCA_XY, mdc_vtx_tracks_config.GetFieldId("dca_xy")));
  fields_id_.insert(std::make_pair(FIELDS::DCA_Z, mdc_vtx_tracks_config.GetFieldId("dca_z")));
  fields_id_.insert(std::make_pair(FIELDS::WALL_RING, wall_hits_config.GetFieldId("ring")));

  // initializing histograms
  pt_rapidity_chi2_ = new TProfile2D( "pt_rapidity_chi2", ";y;p_{T};#chi^{2}", 100, -1.0, 1.0, 100, 0.0, 2.0 );
  pt_rapidity_dca_xy_ = new TProfile2D( "pt_rapidity_dca_xy", ";y;p_{T};DCA_{xy}", 100, -1.0, 1.0, 100, 0.0, 2.0 );
  pt_rapidity_dca_z_ = new TProfile2D( "pt_rapidity_dca_z", ";y;p_{T};DCA_{z}", 100, -1.0, 1.0, 100, 0.0, 2.0 );
  vtx_z_distribution_ = new TH1F( "vtx_z", ";VTX_{z} [mm];counts", 240, -100.0, 20.0 );
  vtx_z_vtx_r_distribution_ = new TH2F("vtx_z_vtx_r", ";VTX_{z} [mm];#sqrt{VTX_{x}^{2}+VTX_{y}^{2}}", 240, -100.0, 20.0, 250, 0.0, 10.0);
  vtx_z_multiplicity_distribution_ = new TH2F("vtx_z_n_tracks", ";VTX_{z} [mm];Hits TOF+RPC", 240, -100.0, 20.0, 250, 0.0, 250.0);
  vtx_x_vtx_y_distribution_ = new TH2F("vtx_x_vtx_y", ";VTX_{x} [mm];VTX_{y} [mm]", 250, -10.0, 10.0, 250, -10.0, 10.0);

  multiplicities_matrix_.insert( std::make_pair(
      MULTIPLICITIES::HITS_TOF, std::map<MULTIPLICITIES, TH2F*>{
                                                            std::pair(MULTIPLICITIES::HITS_RPC, new TH2F( "hits_tof_hits_rpc", ";Hits TOF;Hits RPC", 100, 0, 100.0, 200, 0.0, 200.0 ) ),
                                                            std::pair(MULTIPLICITIES::TRACKS_MDC, new TH2F( "hits_tof_tracks_mdc", ";Hits TOF;Tracks MDC", 100, 0, 100.0, 100, 0.0, 100.0 ) ),
                                                            std::pair(MULTIPLICITIES::FW_ALL_SIGNAL, new TH2F( "hits_tof_fw_all_signal", ";Hits TOF;FW-all signal", 100, 0, 100.0, 200, 0.0, 10000.0 ) ),
                                                            std::pair(MULTIPLICITIES::FW_1_SIGNAL, new TH2F( "hits_tof_fw_1_signal", ";Hits TOF;FW-1 signal", 100, 0, 100.0, 200, 0.0, 10000.0 ) ),
                                                            std::pair(MULTIPLICITIES::FW_2_SIGNAL, new TH2F( "hits_tof_fw_2_signal", ";Hits TOF;FW-2 signal", 100, 0, 100.0, 200, 0.0, 10000.0 ) ),
                                                            std::pair(MULTIPLICITIES::FW_3_SIGNAL, new TH2F( "hits_tof_fw_3_signal", ";Hits TOF;FW-3 signal", 100, 0, 100.0, 200, 0.0, 10000.0 ) ),
                                                                         }));
  multiplicities_matrix_.insert( std::make_pair(
      MULTIPLICITIES::HITS_RPC, std::map<MULTIPLICITIES, TH2F*>{
                                                            std::pair(MULTIPLICITIES::TRACKS_MDC, new TH2F( "hits_rpc_tracks_mdc", ";Hits RPC;Tracks MDC", 200, 0, 200.0, 100, 0.0, 100.0 ) ),
                                                            std::pair(MULTIPLICITIES::FW_ALL_SIGNAL, new TH2F( "hits_rpc_fw_all_signal", ";Hits RPC;FW-all signal", 200, 0, 200.0, 100, 0.0, 10000.0 ) ),
                                                            std::pair(MULTIPLICITIES::FW_1_SIGNAL, new TH2F( "hits_rpc_fw_1_signal", ";Hits TOF;FW-1 signal", 200, 0, 200.0, 200, 0.0, 10000.0 ) ),
                                                            std::pair(MULTIPLICITIES::FW_2_SIGNAL, new TH2F( "hits_rpc_fw_2_signal", ";Hits TOF;FW-2 signal", 200, 0, 200.0, 200, 0.0, 10000.0 ) ),
                                                            std::pair(MULTIPLICITIES::FW_3_SIGNAL, new TH2F( "hits_rpc_fw_3_signal", ";Hits TOF;FW-3 signal", 200, 0, 200.0, 200, 0.0, 10000.0 ) ),
                                                                         }));
  multiplicities_matrix_.insert( std::make_pair(
      MULTIPLICITIES::TRACKS_MDC, std::map<MULTIPLICITIES, TH2F*>{
                                                            std::pair(MULTIPLICITIES::FW_ALL_SIGNAL, new TH2F( "tracks_mdc_fw_signal", ";Tracks MDC;FW signal", 100, 0, 100.0, 100, 0.0, 10000.0 ) ),
                                                            std::pair(MULTIPLICITIES::FW_1_SIGNAL, new TH2F( "tracks_mdc_fw_1_signal", ";Tracks MDC;FW-1 signal", 100, 0, 100.0, 200, 0.0, 10000.0 ) ),
                                                            std::pair(MULTIPLICITIES::FW_2_SIGNAL, new TH2F( "tracks_mdc_fw_2_signal", ";Tracks MDC;FW-2 signal", 100, 0, 100.0, 200, 0.0, 10000.0 ) ),
                                                            std::pair(MULTIPLICITIES::FW_3_SIGNAL, new TH2F( "tracks_mdc_fw_3_signal", ";Tracks MDC;FW-3 signal", 100, 0, 100.0, 200, 0.0, 10000.0 ) ),
                                                                           }));
  multiplicities_matrix_.insert( std::make_pair(
      MULTIPLICITIES::FW_ALL_SIGNAL, std::map<MULTIPLICITIES, TH2F*>{
                                                            std::pair(MULTIPLICITIES::FW_1_SIGNAL, new TH2F( "fw_all_fw_1_signal", ";FW-all signal;FW-1 signal", 200, 0, 10000.0, 200, 0.0, 10000.0 ) ),
                                                            std::pair(MULTIPLICITIES::FW_2_SIGNAL, new TH2F( "fw_all_fw_2_signal", ";FW-all signal;FW-2 signal", 200, 0, 10000.0, 200, 0.0, 10000.0 ) ),
                                                            std::pair(MULTIPLICITIES::FW_3_SIGNAL, new TH2F( "fw_all_fw_3_signal", ";FW-all signal;FW-3 signal", 200, 0, 10000.0, 200, 0.0, 10000.0 ) ),
                                                                           }));
  multiplicities_matrix_.insert( std::make_pair(
      MULTIPLICITIES::FW_1_SIGNAL, std::map<MULTIPLICITIES, TH2F*>{
                                                            std::pair(MULTIPLICITIES::FW_2_SIGNAL, new TH2F( "fw_1_fw_2_signal", ";FW-1 signal;FW-2 signal", 200, 0, 10000.0, 200, 0.0, 10000.0 ) ),
                                                            std::pair(MULTIPLICITIES::FW_3_SIGNAL, new TH2F( "fw_1_fw_3_signal", ";FW-1 signal;FW-3 signal", 200, 0, 10000.0, 200, 0.0, 10000.0 ) ),
                                                                           }));
  multiplicities_matrix_.insert( std::make_pair(
      MULTIPLICITIES::FW_2_SIGNAL, std::map<MULTIPLICITIES, TH2F*>{
                                                            std::pair(MULTIPLICITIES::FW_3_SIGNAL, new TH2F( "fw_2_fw_3_signal", ";FW-2 signal;FW-3 signal", 200, 0, 10000.0, 200, 0.0, 10000.0 ) ),
                                                                           }));
  multiplicities_matrix_track_values_.insert( std::make_pair(
      MULTIPLICITIES::HITS_TOF, std::map<TRACK_VALUES, TH2F*>{
          std::pair(TRACK_VALUES::ERAT, new TH2F( "hits_tof_erat", ";Hits TOF;ERAT", 100, 0, 100.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::PRAT, new TH2F( "hits_tof_prat", ";Hits TOF;PRAT", 100, 0, 100.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_PT, new TH2F( "hits_tof_mean_pt", ";Hits TOF;<p_{T}> [GeV/c]", 100, 0, 100.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_PZ, new TH2F( "hits_tof_mean_pz", ";Hits TOF;<p_{z}> [GeV/c]", 100, 0, 100.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_Y, new TH2F( "hits_tof_mean_y", ";Hits TOF;<y>", 100, 0, 100.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_THETA, new TH2F( "hits_tof_mean_theta", ";Hits TOF;<#Theta> [rad]", 100, 0, 100.0, 200, 0.0, 2.0 ) ),
      }));
  multiplicities_matrix_track_values_.insert( std::make_pair(
      MULTIPLICITIES::TRACKS_MDC, std::map<TRACK_VALUES, TH2F*>{
          std::pair(TRACK_VALUES::ERAT, new TH2F( "tracks_mdc_erat", ";Tracks MDC;ERAR", 100, 0, 100.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::PRAT, new TH2F( "tracks_mdc_prat", ";Tracks MDC;PRAT", 100, 0, 100.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_PT, new TH2F( "tracks_mdc_mean_pt", ";Tracks MDC;<p_{T}> [GeV/c]", 100, 0, 100.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_PZ, new TH2F( "tracks_mdc_mean_pz", ";Tracks MDC;<p_{z}> [GeV/c]", 100, 0, 100.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_Y, new TH2F( "tracks_mdc_mean_y", ";Tracks MDC;<y>", 100, 0, 100.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_THETA, new TH2F( "tracks_mdc_mean_theta", ";Tracks MDC;<#Theta> [rad]", 100, 0, 100.0, 200, 0.0, 2.0 ) ),
      }));
  multiplicities_matrix_track_values_.insert( std::make_pair(
      MULTIPLICITIES::HITS_RPC, std::map<TRACK_VALUES, TH2F*>{
          std::pair(TRACK_VALUES::ERAT, new TH2F( "hits_rpc_erat", ";Hits RPC;ERAR", 200, 0, 200.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::PRAT, new TH2F( "hits_rpc_prat", ";Hits RPC;PRAT", 200, 0, 200.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_PT, new TH2F( "hits_rpc_mean_pt", ";Hits RPC;<p_{T}> [GeV/c]", 200, 0, 200.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_PZ, new TH2F( "hits_rpc_mean_pz", ";Hits RPC;<p_{z}> [GeV/c]", 200, 0, 200.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_Y, new TH2F( "hits_rpc_mean_y", ";Hits RPC;<y>", 200, 0, 200.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_THETA, new TH2F( "hits_rpc_mean_theta", ";Hits RPC;<#Theta> [rad]", 200, 0, 200.0, 200, 0.0, 2.0 ) ),
      }));
  multiplicities_matrix_track_values_.insert( std::make_pair(
      MULTIPLICITIES::FW_ALL_SIGNAL, std::map<TRACK_VALUES, TH2F*>{
          std::pair(TRACK_VALUES::ERAT, new TH2F( "fw_signal_erat", ";FW-all signal;ERAR", 200, 0, 10000.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::PRAT, new TH2F( "fw_signal_prat", ";FW-all signal;PRAT", 200, 0, 10000.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_PT, new TH2F( "fw_signal_mean_pt", ";FW-all Signal;<p_{T}> [GeV/c]", 200, 0, 10000.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_PZ, new TH2F( "fw_signal_mean_pz", ";FW-all Signal;<p_{z}> [GeV/c]", 200, 0, 10000.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_Y, new TH2F( "fw_signal_mean_y", ";FW-all Signal;<y>", 200, 0, 10000.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_THETA, new TH2F( "fw_signal_mean_theta", ";FW-all Signal;<#Theta> [rad]", 200, 0, 10000.0, 200, 0.0, 2.0 ) ),
      }));
  multiplicities_matrix_track_values_.insert( std::make_pair(
      MULTIPLICITIES::FW_1_SIGNAL, std::map<TRACK_VALUES, TH2F*>{
          std::pair(TRACK_VALUES::ERAT, new TH2F( "fw_1_signal_erat", ";FW1 signal;ERAR", 200, 0, 10000.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::PRAT, new TH2F( "fw_1_signal_prat", ";FW1 signal;PRAT", 200, 0, 10000.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_PT, new TH2F( "fw_1_signal_mean_pt", ";FW1 Signal;<p_{T}> [GeV/c]", 200, 0, 10000.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_PZ, new TH2F( "fw_1_signal_mean_pz", ";FW1 Signal;<p_{z}> [GeV/c]", 200, 0, 10000.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_Y, new TH2F( "fw_1_signal_mean_y", ";FW1 Signal;<y>", 200, 0, 10000.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_THETA, new TH2F( "fw_1_signal_mean_theta", ";FW1 Signal;<#Theta> [rad]", 200, 0, 10000.0, 200, 0.0, 2.0 ) ),
      }));
  multiplicities_matrix_track_values_.insert( std::make_pair(
      MULTIPLICITIES::FW_2_SIGNAL, std::map<TRACK_VALUES, TH2F*>{
          std::pair(TRACK_VALUES::ERAT, new TH2F( "fw_2_signal_erat", ";FW2 signal;ERAR", 200, 0, 10000.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::PRAT, new TH2F( "fw_2_signal_prat", ";FW2 signal;PRAT", 200, 0, 10000.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_PT, new TH2F( "fw_2_signal_mean_pt", ";FW2 Signal;<p_{T}> [GeV/c]", 200, 0, 10000.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_PZ, new TH2F( "fw_2_signal_mean_pz", ";FW2 Signal;<p_{z}> [GeV/c]", 200, 0, 10000.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_Y, new TH2F( "fw_2_signal_mean_y", ";FW2 Signal;<y>", 200, 0, 10000.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_THETA, new TH2F( "fw_2_signal_mean_theta", ";FW2 Signal;<#Theta> [rad]", 200, 0, 10000.0, 200, 0.0, 2.0 ) ),
      }));
  multiplicities_matrix_track_values_.insert( std::make_pair(
      MULTIPLICITIES::FW_3_SIGNAL, std::map<TRACK_VALUES, TH2F*>{
          std::pair(TRACK_VALUES::ERAT, new TH2F( "fw_3_signal_erat", ";FW3 signal;ERAR", 200, 0, 10000.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::PRAT, new TH2F( "fw_3_signal_prat", ";FW3 signal;PRAT", 200, 0, 10000.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_PT, new TH2F( "fw_3_signal_mean_pt", ";FW3 Signal;<p_{T}> [GeV/c]", 200, 0, 10000.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_PZ, new TH2F( "fw_3_signal_mean_pz", ";FW3 Signal;<p_{z}> [GeV/c]", 200, 0, 10000.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_Y, new TH2F( "fw_3_signal_mean_y", ";FW3 Signal;<y>", 200, 0, 10000.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_THETA, new TH2F( "fw_3_signal_mean_theta", ";FW3 Signal;<#Theta> [rad]", 200, 0, 10000.0, 200, 0.0, 2.0 ) ),
      }));
  track_values_matrix_.insert(std::make_pair(
      TRACK_VALUES::ERAT, std::map<TRACK_VALUES, TH2F*>{
          std::pair(TRACK_VALUES::PRAT, new TH2F( "erat_prat", ";ERAT;PRAT", 200, 0, 2.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_PT, new TH2F( "erat_mean_pt", ";ERAT;<p_{T}> [GeV/c]", 200, 0, 2.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_PZ, new TH2F( "erat_mean_pz", ";ERAT;<p_{z}> [GeV/c]", 200, 0, 2.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_Y, new TH2F( "erat_mean_y", ";ERAT;<y>", 200, 0, 2.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_THETA, new TH2F( "erat_mean_theta", ";ERAT;<#Theta> [rad]", 200, 0, 2.0, 200, 0.0, 2.0 ) ),
      }));
  track_values_matrix_.insert(std::make_pair(
      TRACK_VALUES::PRAT, std::map<TRACK_VALUES, TH2F*>{
          std::pair(TRACK_VALUES::MEAN_PT, new TH2F( "prat_mean_pt", ";PRAT;<p_{T}> [GeV/c]", 200, 0, 2.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_PZ, new TH2F( "prat_mean_pz", ";PRAT;<p_{z}> [GeV/c]", 200, 0, 2.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_Y, new TH2F( "prat_mean_y", ";PRAT;<y>", 200, 0, 2.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_THETA, new TH2F( "prat_mean_theta", ";PRAT;<#Theta> [rad]", 200, 0, 2.0, 200, 0.0, 2.0 ) ),
      }));
  track_values_matrix_.insert(std::make_pair(
      TRACK_VALUES::MEAN_PT, std::map<TRACK_VALUES, TH2F*>{
          std::pair(TRACK_VALUES::MEAN_PZ, new TH2F( "mean_pt_mean_pz", ";<p_{T}> [GeV/c];<p_{z}> [GeV/c]", 200, 0, 2.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_Y, new TH2F( "mean_pt_mean_y", ";<p_{T}> [GeV/c];<y>", 200, 0, 2.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_THETA, new TH2F( "mean_pt_mean_theta", ";<p_{T}> [GeV/c];<#Theta> [rad]", 200, 0, 2.0, 200, 0.0, 2.0 ) ),
      }));
  track_values_matrix_.insert(std::make_pair(
      TRACK_VALUES::MEAN_PZ, std::map<TRACK_VALUES, TH2F*>{
          std::pair(TRACK_VALUES::MEAN_Y, new TH2F( "mean_pz_mean_y", ";<p_{z}> [GeV/c];<y>", 200, 0, 2.0, 200, 0.0, 2.0 ) ),
          std::pair(TRACK_VALUES::MEAN_THETA, new TH2F( "mean_pz_mean_theta", ";<p_{z}> [GeV/c];<#Theta> [rad]", 200, 0, 2.0, 200, 0.0, 2.0 ) ),
      }));
  track_values_matrix_.insert(std::make_pair(
      TRACK_VALUES::MEAN_Y, std::map<TRACK_VALUES, TH2F*>{
          std::pair(TRACK_VALUES::MEAN_THETA, new TH2F( "mean_y_mean_theta", ";<p_{z}> [GeV/c];<#Theta> [rad]", 200, 0, 2.0, 200, 0.0, 2.0 ) ),
      }));
}

void AnalysisTask::Exec() {
  std::vector<MULTIPLICITIES> mult_keys{
      MULTIPLICITIES::HITS_TOF,
      MULTIPLICITIES::HITS_RPC,
      MULTIPLICITIES::TRACKS_MDC,
      MULTIPLICITIES::FW_ALL_SIGNAL,
      MULTIPLICITIES::FW_1_SIGNAL,
      MULTIPLICITIES::FW_2_SIGNAL,
      MULTIPLICITIES::FW_3_SIGNAL,
  };
  std::vector<TRACK_VALUES> track_val_keys{
      TRACK_VALUES::ERAT,
      TRACK_VALUES::PRAT,
      TRACK_VALUES::MEAN_PT,
      TRACK_VALUES::MEAN_PZ,
      TRACK_VALUES::MEAN_Y,
      TRACK_VALUES::MEAN_THETA,
  };
  std::map<MULTIPLICITIES, int> multiplicities;
  std::map<TRACK_VALUES, float> track_values;
  if( !event_header_->GetField<bool>(fields_id_.at(FIELDS::PT2)) )
    return;
//  if( !(-90 < event_header_->GetField<float>(-3) && event_header_->GetField<float>(-3) < -70) )
//    return;
  multiplicities.insert(std::make_pair(MULTIPLICITIES::HITS_TOF,event_header_->GetField<int>(fields_id_.at(FIELDS::HITS_TOF)))); // getting multiplicity from event header
  multiplicities.insert(std::make_pair(MULTIPLICITIES::HITS_RPC,event_header_->GetField<int>(fields_id_.at(FIELDS::HITS_RPC)))); // getting multiplicity from event header
  multiplicities.insert(std::make_pair(MULTIPLICITIES::TRACKS_MDC,event_header_->GetField<int>(fields_id_.at(FIELDS::TRACKS_MDC)))); // getting multiplicity from event header
  multiplicities.insert(std::make_pair(MULTIPLICITIES::FW_ALL_SIGNAL,event_header_->GetField<int>(fields_id_.at(FIELDS::FW_SIGNAL)))); // getting multiplicity from event header
  auto vtx_x = event_header_->GetVertexX();
  auto vtx_y = event_header_->GetVertexY();
  auto vtx_z = event_header_->GetVertexZ();
  auto vtx_r = sqrt(vtx_x*vtx_x+vtx_y*vtx_y);
  int n_tracks = mdc_vtx_tracks_->GetNumberOfChannels(); // number of tracks in current event

  vtx_z_distribution_->Fill(vtx_z);
  vtx_x_vtx_y_distribution_->Fill(vtx_x,vtx_y);
  vtx_z_vtx_r_distribution_->Fill(vtx_z, vtx_r);
  vtx_z_multiplicity_distribution_->Fill( vtx_z,
                                          event_header_->GetField<int>(fields_id_.at(FIELDS::HITS_TOF))+
                                              event_header_->GetField<int>(fields_id_.at(FIELDS::HITS_RPC)));

  TLorentzVector sum4P{0.0, 0.0, 0.0, 0.0};
  int n_protons=0;
  double erat_x=0;
  double erat_y=0;
  double prat_x=0;
  double prat_y=0;
  double mean_pt=0;
  double mean_pl=0;
  double mean_theta=0;
  double mean_y=0;
  for (size_t i = 0; i < n_tracks; ++i) { // loop over all tracks if current event
    auto track = mdc_vtx_tracks_->GetChannel(i); // getting track from track detector
    int match_meta_hit = mdc_meta_matching_->GetMatchDirect(i); // getting index of matched with track TOF-system hit
    auto hit = meta_hits_->GetChannel(i); // getting matched with track hit in TOF-system
    auto mom4 = track.Get4MomentumByMass( track.GetMass() );
    sum4P+=mom4;
    auto pT = mom4.Pt(); // getting transverse momentum
    auto eta = mom4.Eta(); // getting pseudorapidity
    auto p = mom4.P(); // getting absolute value of momentum
    auto y = mom4.Rapidity()-0.74;
    auto E = mom4.Energy();
    auto theta = mom4.Theta();
    mean_pt+=pT;
    mean_pl+=mom4.Pz();
    mean_theta+=theta;
    mean_y+=mom4.Rapidity();

    erat_x+=E*cos(theta);
    erat_y+=E*sin(theta);
    prat_x+=p*cos(theta);
    prat_y+=p*sin(theta);
    if( track.GetPid()!=2212 ) // protons
      continue;
    auto chi2 = track.GetField<float>(fields_id_.at(FIELDS::CHI_2));
    auto dca_xy = track.GetField<float>(fields_id_.at(FIELDS::DCA_XY));
    auto dca_z = track.GetField<float>(fields_id_.at(FIELDS::DCA_Z));
    // filling distributions
    pt_rapidity_chi2_->Fill(y, pT, chi2);
    pt_rapidity_dca_xy_->Fill(y, pT, fabs(dca_xy));
    pt_rapidity_dca_z_->Fill(y, pT, fabs(dca_z));
    n_protons++;
  }
  auto n_modules = wall_hits_->GetNumberOfChannels();
  float signal_w1=0.0;
  float signal_w2=0.0;
  float signal_w3=0.0;
  for( size_t i=0; i<n_modules; ++i ){
    auto hit = wall_hits_->GetChannel(i);
    auto signal = hit.GetSignal();
    auto ring = hit.GetField<int>(fields_id_.at(FIELDS::WALL_RING));
    if( ring <= 5 )
      signal_w1+=signal;
    if (ring==6 || ring==7)
      signal_w2+=signal;
    if(ring >= 8 && ring <= 10)
      signal_w3+=signal;
  }
  multiplicities.insert(std::make_pair(MULTIPLICITIES::FW_1_SIGNAL,signal_w1)); // getting multiplicity from event header
  multiplicities.insert(std::make_pair(MULTIPLICITIES::FW_2_SIGNAL,signal_w2)); // getting multiplicity from event header
  multiplicities.insert(std::make_pair(MULTIPLICITIES::FW_3_SIGNAL,signal_w3)); // getting multiplicity from event header

  double erat =erat_y/erat_x;
  double prat =prat_y/prat_x;
  sum4P=sum4P*(1.0/ (double) n_tracks);
  mean_pt/= (double) n_tracks;
  mean_pl/= (double) n_tracks;
  mean_y/= (double) n_tracks;
  mean_theta/= (double) n_tracks;
  track_values.insert(std::make_pair( TRACK_VALUES::ERAT, erat ));
  track_values.insert(std::make_pair( TRACK_VALUES::PRAT, prat ));
  track_values.insert(std::make_pair( TRACK_VALUES::MEAN_PT, mean_pt ));
  track_values.insert(std::make_pair( TRACK_VALUES::MEAN_PZ, mean_pl ));
  track_values.insert(std::make_pair( TRACK_VALUES::MEAN_Y, mean_y ));
  track_values.insert(std::make_pair( TRACK_VALUES::MEAN_THETA, mean_theta ));
  for( size_t i=0; i<std::size( mult_keys ); ++i )
    for( size_t j=i+1; j<std::size( mult_keys ); ++j )
      multiplicities_matrix_.at(mult_keys.at(i)).at(mult_keys.at(j))->
          Fill(multiplicities.at( mult_keys.at(i) ),multiplicities.at( mult_keys.at(j) ));
  for( size_t i=0; i<std::size( track_val_keys ); ++i )
    for( size_t j=i+1; j<std::size( track_val_keys ); ++j )
      track_values_matrix_.at(track_val_keys.at(i)).at(track_val_keys.at(j))->
          Fill(track_values.at( track_val_keys.at(i) ),track_values.at( track_val_keys.at(j) ));
  for( size_t i=0; i<std::size( mult_keys ); ++i )
    for( size_t j=0; j<std::size( track_val_keys ); ++j )
      multiplicities_matrix_track_values_.at(mult_keys.at(i)).at(track_val_keys.at(j))->
          Fill(multiplicities.at( mult_keys.at(i) ),track_values.at( track_val_keys.at(j) ));
}

void AnalysisTask::Finish() {
  // Writing histograms to file

  vtx_z_distribution_->Write();
  vtx_x_vtx_y_distribution_->Write();
  vtx_z_vtx_r_distribution_->Write();
  vtx_z_multiplicity_distribution_->Write();
  pt_rapidity_chi2_->Write();
  pt_rapidity_dca_z_->Write();
  pt_rapidity_dca_xy_->Write();


  for( const auto& matrices : multiplicities_matrix_ )
    for( const auto& matrix : matrices.second )
      matrix.second->Write();
  for( const auto& matrices : multiplicities_matrix_track_values_ )
    for( const auto& matrix : matrices.second )
      matrix.second->Write();
  for( const auto& matrices : track_values_matrix_ )
    for( const auto& matrix : matrices.second )
      matrix.second->Write();
}
} // namespace AnalysisTree