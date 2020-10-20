//
// Created by mikhail on 6/16/20.
//

#ifndef QUALITY_ASSURANCE_SRC_TREE_READER_H_
#define QUALITY_ASSURANCE_SRC_TREE_READER_H_

#include <TChain.h>
#include <TH3F.h>
#include <TProfile2D.h>

#include <AnalysisTree/EventHeader.hpp>
#include <AnalysisTree/FillTask.hpp>
#include <AnalysisTree/Cuts.hpp>
#include <AnalysisTree/Detector.hpp>
#include <AnalysisTree/Matching.hpp>

namespace AnalysisTree {
class AnalysisTask : public FillTask{
public:
 AnalysisTask() = default;
  ~AnalysisTask() override = default;
  void Init( std::map<std::string, void*>& branch_map ) override;
  void Exec() override;
  void Finish() override;
private:
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
    PT2
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
    MEAN_THETA
  };
  std::map<FIELDS, int> fields_id_; // map to match detectors' fields with enumerator

  /* pointers to link tree's branches with */
  EventHeader* event_header_{nullptr}; 		// event info
  Particles* mdc_vtx_tracks_{nullptr}; 		// tracks
  HitDetector* meta_hits_{nullptr}; 		// TOF-system
  HitDetector* wall_hits_{nullptr}; 		// FW-system
  Matching* mdc_meta_matching_{nullptr}; 	// matching between tracking system and TOF-system
  TH1F* vtx_z_distribution_;
  TH2F* vtx_z_multiplicity_distribution_;
  TH2F* vtx_z_vtx_r_distribution_;
  TH2F* vtx_x_vtx_y_distribution_;
  std::map<MULTIPLICITIES,std::map<MULTIPLICITIES,TH2F*>>
      multiplicities_matrix_;
  std::map<MULTIPLICITIES,std::map<TRACK_VALUES,TH2F*>>
      multiplicities_matrix_track_values_;
  std::map<TRACK_VALUES,std::map<TRACK_VALUES,TH2F*>>
      track_values_matrix_;
  TProfile2D* pt_rapidity_chi2_;
  TProfile2D* pt_rapidity_dca_xy_;
  TProfile2D* pt_rapidity_dca_z_;
};
} // namespace AnalysisTree
#endif // QUALITY_ASSURANCE_SRC_TREE_READER_H_
