//
// Created by mikhail on 6/16/20.
//

#ifndef QUALITY_ASSURANCE_SRC_TREE_READER_H_
#define QUALITY_ASSURANCE_SRC_TREE_READER_H_

#include <TChain.h>
#include <TFile.h>
#include <TH3F.h>
#include <TProfile2D.h>

#include <AnalysisTree/EventHeader.hpp>
#include <AnalysisTree/FillTask.hpp>
#include <AnalysisTree/Cuts.hpp>
#include <AnalysisTree/Detector.hpp>
#include <AnalysisTree/Matching.hpp>

namespace AnalysisTree {

struct Axis{
  std::string name;
  std::string title;
  int n_bins;
  double min;
  double max;
};

class AnalysisTask : public FillTask{
public:
 AnalysisTask() = default;
  ~AnalysisTask() override = default;
  void Init( std::map<std::string, void*>& branch_map ) override;
  void Exec() override;
  void Finish() override;
  void InitEffieciencies(const std::string& file_name);
private:
  TH2F* Make2DHisto( Axis first, Axis second ){
    std::string name = first.name + "_" + second.name;
    std::string title = ";" + first.title + ";"+second.title;
    auto histo = new TH2F( name.c_str(), title.c_str(),
                          first.n_bins, first.min, first.max,
                          second.n_bins, second.min, second.max
                          );
    return histo;
  }
  TH3F* Make3DHisto( Axis first, Axis second, Axis third ){
    std::string name = first.name + "_" + second.name+"_"+third.name;
    std::string title = ";" + first.title + ";"+second.title+ ";"+third.title;
    auto histo = new TH3F( name.c_str(), title.c_str(),
                          first.n_bins, first.min, first.max,
                          second.n_bins, second.min, second.max,
                          third.n_bins, third.min, third.max
                          );
    return histo;
  }
 enum class FIELDS { // enumerator to fast access to detectors' fields
    HITS_TOF,         // Hits in TOF-system
    HITS_RPC,         // Hits in RPC-system
    TRACKS_MDC,
    FW_SIGNAL,
    CHI_2,
    DCA_XY,
    DCA_Z,
    WALL_RING,
    PT3,
    PT2,
    GEANT_ID
  };
  enum class MULTIPLICITIES {
   HITS_TOF,
   HITS_RPC,
   TRACKS_MDC,
   FW_ALL_SIGNAL,
   FW_1_SIGNAL,
   FW_2_SIGNAL,
   FW_3_SIGNAL,
 };
  enum class TRACK_VALUES {
    PRAT,
    ERAT,
    MEAN_Y,
    MEAN_PT,
    MEAN_PZ,
    MEAN_THETA,
    MEAN_YCM,
    MEAN_YCM_NO_EFF,
    FW_VS_BW,
    FW_VS_BW_NO_EFF,
  };
  std::map<FIELDS, int> fields_id_; // map to match detectors' fields with enumerator
  std::map<MULTIPLICITIES, Axis> multiplicities_axes_{
      std::pair( MULTIPLICITIES::HITS_TOF, Axis{ "hits_tof", "N hits TOF", 100, 0.0, 100.0 } ),
      std::pair( MULTIPLICITIES::HITS_RPC, Axis{ "hits_rpc", "N hits RPC", 200, 0.0, 200.0 } ),
      std::pair( MULTIPLICITIES::TRACKS_MDC, Axis{ "tracks_mdc", "N tracks MDC", 200, 0.0, 200.0 } ),
      std::pair( MULTIPLICITIES::FW_ALL_SIGNAL, Axis{ "fw_all_signal", "FW-all signal", 200, 0.0, 10000.0 } ),
      std::pair( MULTIPLICITIES::FW_1_SIGNAL, Axis{ "fw1_signal", "FW 1 signal", 200, 0.0, 10000.0 } ),
      std::pair( MULTIPLICITIES::FW_2_SIGNAL, Axis{ "fw2_signal", "FW 2 signal", 200, 0.0, 10000.0 } ),
      std::pair( MULTIPLICITIES::FW_3_SIGNAL, Axis{ "fw3_signal", "FW 3 signal", 200, 0.0, 10000.0 } ),
  };
  std::map<TRACK_VALUES, Axis> track_values_axes_{
      std::pair( TRACK_VALUES::ERAT, Axis{ "erat", "ERAT", 200, 0.0, 2.0 } ),
      std::pair( TRACK_VALUES::PRAT, Axis{ "prat", "PRAT", 200, 0.0, 2.0 } ),
      std::pair( TRACK_VALUES::MEAN_PT, Axis{ "mean_pT", "<p_{T}> [GeV/c]", 200, 0.0, 2.0 } ),
      std::pair( TRACK_VALUES::MEAN_PZ, Axis{ "mean_pz", "<p_{z}> [GeV/c]", 200, 0.0, 2.0 } ),
      std::pair( TRACK_VALUES::MEAN_Y, Axis{ "mean_y", "<y>", 200, 0.0, 2.0 } ),
      std::pair( TRACK_VALUES::MEAN_YCM, Axis{ "mean_protons_ycm", "<y_{cm}>", 200, -1.0, 1.0 } ),
      std::pair( TRACK_VALUES::MEAN_YCM_NO_EFF, Axis{ "no_eff_mean_protons_ycm", "<y_{cm}>", 200, -1.0, 1.0 } ),
      std::pair( TRACK_VALUES::FW_VS_BW, Axis{ "bw_vs_fw_protons", "BW-FW", 60, -30.0, 30.0 } ),
      std::pair( TRACK_VALUES::FW_VS_BW_NO_EFF, Axis{ "no_eff_bw_vs_fw_protons", "BW-FW", 60, -30.0, 30.0 } ),
      std::pair( TRACK_VALUES::MEAN_THETA, Axis{ "mean_theta", "<#theta>", 200, 0.0, 2.0 } ),
  };
  /* pointers to link tree's branches with */
  EventHeader* event_header_{nullptr}; 		// event info
  Particles* mdc_vtx_tracks_{nullptr}; 		// tracks
  HitDetector* meta_hits_{nullptr}; 		// TOF-system
  HitDetector* wall_hits_{nullptr}; 		// FW-system
  Matching* mdc_meta_matching_{nullptr}; 	// matching between tracking system and TOF-system
  TH1F* vtx_z_distribution_;
  TH1F* n_pions_to_all_tracks_;
  TH2F* vtx_z_multiplicity_distribution_;
  TH2F* vtx_z_vtx_r_distribution_;
  TH2F* vtx_x_vtx_y_distribution_;
  std::map<MULTIPLICITIES,std::map<MULTIPLICITIES,TH2F*>>
      multiplicities_matrix_;
  std::map<MULTIPLICITIES,std::map<TRACK_VALUES,TH2F*>>
      multiplicities_track_values_matrix_;
  std::map<TRACK_VALUES,std::map<TRACK_VALUES,TH2F*>>
      track_values_matrix_;
  TProfile2D* pt_rapidity_chi2_;
  TProfile2D* pt_rapidity_dca_xy_;
  TProfile2D* pt_rapidity_dca_z_;
  TH3F* n_tracks_erat_protons_y_;
  TFile* file_efficiency_protons_;
  std::vector<TH2F*> efficiencies_;
};
} // namespace AnalysisTree
#endif // QUALITY_ASSURANCE_SRC_TREE_READER_H_
