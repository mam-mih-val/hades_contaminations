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
  fields_id_.insert(std::make_pair(FIELDS::GEANT_ID, mdc_vtx_tracks_config.GetFieldId("geant_pid")));
  fields_id_.insert(std::make_pair(FIELDS::DCA_Z, mdc_vtx_tracks_config.GetFieldId("dca_z")));
  fields_id_.insert(std::make_pair(FIELDS::WALL_RING, wall_hits_config.GetFieldId("ring")));

  // initializing histograms
  pt_rapidity_chi2_ = new TProfile2D( "pt_rapidity_chi2", ";y;p_{T};#chi^{2}", 100, -1.0, 1.0, 100, 0.0, 2.0 );
  pt_rapidity_dca_xy_ = new TProfile2D( "pt_rapidity_dca_xy", ";y;p_{T};DCA_{xy}", 100, -1.0, 1.0, 100, 0.0, 2.0 );
  pt_rapidity_dca_z_ = new TProfile2D( "pt_rapidity_dca_z", ";y;p_{T};DCA_{z}", 100, -1.0, 1.0, 100, 0.0, 2.0 );
  vtx_z_distribution_ = new TH1F( "vtx_z", ";VTX_{z} [mm];counts", 240, -100.0, 20.0 );
  n_pions_to_all_tracks_ = new TH1F( "n_pions_to_all_tracks", ";N Pions / All tracks (%);counts", 500, 0.0, 100.0 );
  vtx_z_vtx_r_distribution_ = new TH2F("vtx_z_vtx_r", ";VTX_{z} [mm];#sqrt{VTX_{x}^{2}+VTX_{y}^{2}}", 240, -100.0, 20.0, 250, 0.0, 10.0);
  vtx_z_multiplicity_distribution_ = new TH2F("vtx_z_n_tracks", ";VTX_{z} [mm];Hits TOF+RPC", 240, -100.0, 20.0, 250, 0.0, 250.0);
  vtx_x_vtx_y_distribution_ = new TH2F("vtx_x_vtx_y", ";VTX_{x} [mm];VTX_{y} [mm]", 250, -10.0, 10.0, 250, -10.0, 10.0);

  for( const auto& x : multiplicities_axes_ ){
    std::map<MULTIPLICITIES, TH2F*> x_row;
    for(const auto& y : multiplicities_axes_){
      if (x.first == y.first) // to avoid repetitions
        continue;
      x_row.insert( {y.first, Make2DHisto(x.second, y.second)} );
    }
    multiplicities_matrix_.insert({x.first, x_row});
  }
  for( const auto& x : multiplicities_axes_ ){
    std::map<TRACK_VALUES, TH2F*> x_row;
    for(const auto& y : track_values_axes_){
      x_row.insert( {y.first, Make2DHisto(x.second, y.second)} );
    }
    multiplicities_track_values_matrix_.insert({x.first, x_row});
  }
  for( const auto& x : track_values_axes_ ){
    std::map<TRACK_VALUES, TH2F*> x_row;
    for(const auto& y : track_values_axes_){
      if (x.first == y.first) // to avoid repetitions
        continue;
      x_row.insert( {y.first, Make2DHisto(x.second, y.second)} );
    }
    track_values_matrix_.insert({x.first, x_row});
  }

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
      MULTIPLICITIES::N_PIONS,
      MULTIPLICITIES::N_HELIUM
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
  int n_pions = 0;
  int n_helium = 0;
  for (size_t i = 0; i < n_tracks; ++i) { // loop over all tracks if current event
    auto track = mdc_vtx_tracks_->GetChannel(i); // getting track from track detector
    int match_meta_hit = mdc_meta_matching_->GetMatchDirect(i); // getting index of matched with track TOF-system hit
    auto hit = meta_hits_->GetChannel(i); // getting matched with track hit in TOF-system
    auto mom4 = track.Get4MomentumByMass( track.GetMass() );
    auto pid = track.GetPid();
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
    mean_y+=mom4.Rapidity()*mom4.Rapidity();

    erat_x+=E*cos(theta);
    erat_y+=E*sin(theta);
    prat_x+=p*p*cos(theta);
    prat_y+=p*p*sin(theta);
    if( abs(pid) == 211 )
      n_pions++;
    auto geant_pid = track.GetField<int>(fields_id_.at(FIELDS::GEANT_ID));
    if( geant_pid == 47  || geant_pid == 49  )
      n_helium++;
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
  n_pions_to_all_tracks_->Fill( (double) n_pions / (double) n_tracks * 100.0 );
  multiplicities.insert(std::make_pair(MULTIPLICITIES::FW_1_SIGNAL,signal_w1)); // getting multiplicity from event header
  multiplicities.insert(std::make_pair(MULTIPLICITIES::FW_2_SIGNAL,signal_w2)); // getting multiplicity from event header
  multiplicities.insert(std::make_pair(MULTIPLICITIES::FW_3_SIGNAL,signal_w3)); // getting multiplicity from event header
  multiplicities.insert(std::make_pair(MULTIPLICITIES::N_PIONS,n_pions)); // getting multiplicity from event header
  multiplicities.insert(std::make_pair(MULTIPLICITIES::N_HELIUM,n_helium)); // getting multiplicity from event header

  double erat =erat_y/erat_x;
  double prat =prat_y/prat_x;
  sum4P=sum4P*(1.0/ (double) n_tracks);
  mean_pt/= (double) n_tracks;
  mean_pl/= (double) n_tracks;
  mean_y/= (double) n_tracks;
  mean_theta/= (double) n_tracks;
  auto rel_pions = (double) n_pions / (double ) n_tracks * 100.0;
  auto rel_helium = (double) n_helium / (double ) n_tracks * 100.0;
  track_values.insert(std::make_pair( TRACK_VALUES::ERAT, erat ));
  track_values.insert(std::make_pair( TRACK_VALUES::PRAT, prat ));
  track_values.insert(std::make_pair( TRACK_VALUES::MEAN_PT, mean_pt ));
  track_values.insert(std::make_pair( TRACK_VALUES::MEAN_PZ, mean_pl ));
  track_values.insert(std::make_pair( TRACK_VALUES::MEAN_Y, mean_y ));
  track_values.insert(std::make_pair( TRACK_VALUES::MEAN_THETA, mean_theta ));
  track_values.insert(std::make_pair( TRACK_VALUES::REL_AMOUNT_OF_PIONS, rel_pions));
  track_values.insert(std::make_pair( TRACK_VALUES::REL_AMOUNT_OF_HELIUM, rel_helium));
  for( auto x : multiplicities )
    for(auto y : multiplicities){
      if( x.first == y.first )
        continue;
      multiplicities_matrix_.at(x.first).at(y.first)->Fill( x.second, y.second );
    }
  for( auto x : multiplicities )
    for(auto y : track_values){
      multiplicities_track_values_matrix_.at(x.first).at(y.first)->Fill( x.second, y.second );
    }
  for( auto x : track_values )
    for(auto y : track_values){
      if( x.first == y.first )
        continue;
      track_values_matrix_.at(x.first).at(y.first)->Fill( x.second, y.second );
    }
}

void AnalysisTask::Finish() {
  // Writing histograms to file

  vtx_z_distribution_->Write();
  n_pions_to_all_tracks_->Write();
  vtx_x_vtx_y_distribution_->Write();
  vtx_z_vtx_r_distribution_->Write();
  vtx_z_multiplicity_distribution_->Write();
  pt_rapidity_chi2_->Write();
  pt_rapidity_dca_z_->Write();
  pt_rapidity_dca_xy_->Write();


  for( const auto& matrices : multiplicities_matrix_ )
    for( const auto& matrix : matrices.second )
      matrix.second->Write();
  for( const auto& matrices : multiplicities_track_values_matrix_)
    for( const auto& matrix : matrices.second )
      matrix.second->Write();
  for( const auto& matrices : track_values_matrix_ )
    for( const auto& matrix : matrices.second )
      matrix.second->Write();
}
} // namespace AnalysisTree