#ifndef GRID_VIEWER_H_INCLUDED
#define GRID_VIEWER_H_INCLUDED
#include <grid.h>

#include <set>

#include <glutils.h>
#include <cpputils.h>

#include <boost/multi_array.hpp>

namespace grid
{
  class datapiece_t ;

  class mscomplex_t;

  class disc_rendata_t
  {
  public:

    cellid_t               cellid;
    uint                   index;

    glutils::renderable_t *ren[GRADDIR_COUNT];
    bool                   show[GRADDIR_COUNT];
    glutils::color_t       color[GRADDIR_COUNT];

    disc_rendata_t(cellid_t c,uint i);
    ~disc_rendata_t();

    void render();
    bool update(mscomplex_t *);

    static void init();
    static void cleanup();
  };

  typedef boost::shared_ptr<glutils::renderable_t> renderable_sp_t;

  typedef std::set<boost::shared_ptr<disc_rendata_t> > disc_rendata_sp_set_t;

  typedef boost::shared_ptr<disc_rendata_t> disc_rendata_ptr_t;

  class octtree_piece_rendata:public configurable_t
  {
  public:

    datapiece_t * dp;

    // set externally to control what is rendered
    bool m_bShowCps[gc_grid_dim+1];
    bool m_bShowAllCps;
    bool m_bShowCpLabels;
    bool m_bShowMsGraph;
    bool m_bShowGrad;
    bool m_bShowCancCps;
    bool m_bShowCancMsGraph;

    // set externally .. cleared by render
    bool m_bNeedUpdateDiscRens;

    renderable_sp_t ren_grad[gc_grid_dim];
    renderable_sp_t ren_cp_labels[gc_grid_dim+1];
    renderable_sp_t ren_cp[gc_grid_dim+1];
    renderable_sp_t ren_cp_conns[gc_grid_dim];
    renderable_sp_t ren_canc_cp_labels[gc_grid_dim+1];
    renderable_sp_t ren_canc_cp[gc_grid_dim+1];
    renderable_sp_t ren_canc_cp_conns[gc_grid_dim];

    glutils::bufobj_ptr_t   cp_loc_bo;

    std::vector<disc_rendata_ptr_t> disc_rds;

    disc_rendata_sp_set_t    active_disc_rens;


    void create_disc_rds();
    void update_active_disc_rens();

    void create_cp_loc_bo();
    void create_cp_rens(const rect_t &roi);
    void create_canc_cp_rens(const rect_t &roi);
    void create_grad_rens(const rect_t &roi);

    void render_msgraph_data(double cp_raise = 0.0,double cp_point_size = 4.0) ;

    void render_dataset_data(double grad_raise = 0.0) ;

    octtree_piece_rendata(datapiece_t *);

    static void init();
    static void cleanup();

    // configurable_t interface
  public:
    virtual data_index_t dim();
    virtual bool exchange_field(const data_index_t &,boost::any &);
    virtual eFieldType exchange_header(const int &,boost::any &);
  };

  class data_manager_t;

  class viewer_t:
      public glutils::renderable_t,
      public configurable_t
  {
  public:
    std::vector<octtree_piece_rendata * >  m_grid_piece_rens;
    cellid_t                               m_size;
    rect_t                                 m_roi;
    double                                 m_scale_factor;
    cellid_t                               m_roi_base_pt;

  public:
    bool                                   m_bShowRoiBB;
    bool                                   m_bRebuildRens;
    bool                                   m_bCenterToRoi;

    double                                 m_cp_raise;
    double                                 m_cp_size;

    data_manager_t *                       m_gdm;
    uint                                   m_rawdata_texture;
    std::string                            m_elevation_filename;

  public:

    viewer_t(data_manager_t * ,std::string ef);

    ~viewer_t();

    // ensure normalization of l and u, l < u , dim in {0,1}
    void set_roi_dim_range_nrm(double l,double u,int dim);

    bool init_rawdata_texture();

  private:

    void build_rens();

    // renderable_t interface
  public:

    void init();

    int  render();

    // configurable_t interface
  public:
    virtual data_index_t dim();
    virtual bool exchange_field(const data_index_t &,boost::any &);
    virtual eFieldType exchange_header(const int &,boost::any &);
  };
}
#endif //VIEWER_H_INCLUDED
