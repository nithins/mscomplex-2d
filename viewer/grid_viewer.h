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

  typedef boost::shared_ptr<glutils::renderable_t> renderable_sp_t;

  struct grid_ren_data_t
  {
    glutils::bufobj_ptr_t m_cell_bo;
    glutils::bufobj_ptr_t m_cell_nrm_bo;
    rect_t                m_roi;
    cellid_t              m_size;
    double                m_scale_factor;
    cellid_t              m_roi_base_pt;
    double                m_cp_raise;
    double                m_cp_size;

    glutils::idx_t cellid_to_index(const cellid_t &c) const
    {
      return (m_size[0]*2-1)*c[1] + c[0];
    }
  };

  class disc_rendata_t
  {
  public:

    cellid_t               cellid;
    uint                   index;

    renderable_sp_t        ren[GRADDIR_COUNT];
    bool                   show[GRADDIR_COUNT];
    glutils::color_t       color[GRADDIR_COUNT];

    disc_rendata_t(cellid_t c,uint i);
    ~disc_rendata_t();

    void render();
    bool update(mscomplex_t *,const grid_ren_data_t& grd);
  };

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
    bool m_bShowGrad[2];
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

    std::vector<disc_rendata_ptr_t> disc_rds;

    disc_rendata_sp_set_t    active_disc_rens;

    void create_disc_rds();
    void update_active_disc_rens(const grid_ren_data_t& grd);

    void create_cp_rens(const grid_ren_data_t& grd);
    void create_canc_cp_rens(const grid_ren_data_t& grd);
    void create_grad_rens(const grid_ren_data_t& grd);

    void render_msgraph_data(const grid_ren_data_t& grd);
    void render_dataset_data(const grid_ren_data_t& grd);

    octtree_piece_rendata(datapiece_t *);

    static void init();

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
    typedef boost::multi_array<cell_fn_t,gc_grid_dim>   varray_t;

  public:
    std::vector<octtree_piece_rendata * >  m_grid_piece_rens;

  public:
    bool                                   m_bShowRoiBB;
    bool                                   m_bRebuildRens;
    bool                                   m_bCenterToRoi;
    bool                                   m_bShowSurface;

    data_manager_t *                       m_gdm;
    std::string                            m_elevation_filename;

    renderable_sp_t                        m_surf_ren;
    grid_ren_data_t                        m_ren_data;

  public:

    void init_surf_ren();

    viewer_t(data_manager_t * ,std::string ef);

    ~viewer_t();

    void set_roi_dim_range_nrm(double l,double u,int dim);

    void set_cp_raise_nrm(double r);

    void set_cp_size_nrm(double n);

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
