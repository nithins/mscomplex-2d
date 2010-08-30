#include <sstream>
#include <cstring>
#include <iostream>

#include <boost/algorithm/string_regex.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>
#include <boost/static_assert.hpp>

#include <GL/glew.h>

#include <glutils.h>
#include <GLSLProgram.h>
#include <logutil.h>

#include <grid_viewer.h>
#include <grid_datamanager.h>
#include <grid_mscomplex.h>
#include <grid_mscomplex_ensure.h>
#include <grid_dataset.h>

#include <shadersources.h>

#define static_assert BOOST_STATIC_ASSERT

glutils::color_t g_grid_cp_colors[grid::gc_grid_dim+1] =
{
  glutils::color_t(0.0,0.0,1.0),
  glutils::color_t(0.0,1.0,0.0),
  glutils::color_t(1.0,0.0,0.0),
};

glutils::color_t g_grid_grad_colors[grid::gc_grid_dim] =
{
  glutils::color_t(0.20,0.20,0.20 ),
  glutils::color_t(0.30,0.30,0.30 ),
};

glutils::color_t g_disc_colors[grid::GRADDIR_COUNT][grid::gc_grid_dim+1] =
{
  {
    glutils::color_t(0.65,0.65,0.65 ),
    glutils::color_t(0.85,0.65,0.75 ),
    glutils::color_t(0.0,0.0,0.0 ),
  },

{
    glutils::color_t(0.0,0.0,0.0 ),
    glutils::color_t(0.65,0.95,0.45 ),
    glutils::color_t(0.15,0.25,0.75 ),
  },
};

glutils::color_t g_grid_cp_conn_colors[grid::gc_grid_dim] =
{
  glutils::color_t(0.65,0.45,0.55 ),
  glutils::color_t(0.45,0.65,0.55 ),
};

glutils::color_t g_roiaabb_color = glutils::color_t(0.85,0.75,0.65);

#ifndef VIEWER_RENDER_AWESOME
double g_max_cp_size  = 8.0;
#else
double g_max_cp_size  = 0.05;
#endif
double g_max_cp_raise = 0.1;

#ifdef VIEWER_RENDER_AWESOME
typedef boost::shared_ptr<GLSLProgram> glsl_program_sp_t;

glsl_program_sp_t s_grad_shader;
glsl_program_sp_t s_sphere_shader;
glsl_program_sp_t s_cylinder_shader;

#endif


namespace grid
{
  void octtree_piece_rendata::init()
  {
#ifdef VIEWER_RENDER_AWESOME
    std::string log;

    s_grad_shader.reset
        (GLSLProgram::createFromSourceStrings
         (grad_vert_glsl,grad_geom_glsl,std::string(),GL_LINES,GL_TRIANGLES));

    s_grad_shader->GetProgramLog(log);

    if(log.size() != 0 )
      throw std::runtime_error("******grad_shader compile error*******\n"+log);

    s_sphere_shader.reset
        (GLSLProgram::createFromSourceStrings
         (sphere_vert_glsl,sphere_geom_glsl,sphere_frag_glsl,GL_POINTS,GL_QUADS));

    s_sphere_shader->GetProgramLog(log);

    if(log.size() != 0 )
      throw std::runtime_error("******sphere_shader compile error*******\n"+log);

    s_cylinder_shader.reset
        (GLSLProgram::createFromSourceStrings
        (cylinder_vert_glsl,cylinder_geom_glsl,cylinder_frag_glsl,GL_LINES,GL_TRIANGLES));

    s_cylinder_shader->GetProgramLog(log);

    if(log.size() != 0 )
      throw std::runtime_error("******cylinder_shader compile error*******\n"+log);
#endif
  }

  viewer_t::viewer_t
      (data_manager_t * gdm,std::string ef):
      m_bRebuildRens(true),m_bShowRoiBB(false),m_bCenterToRoi(false),
      m_bShowSurface(false),
      m_gdm(gdm),m_elevation_filename(ef)
  {

    m_ren_data.m_size         = m_gdm->m_size;

    m_ren_data.m_scale_factor = 0;

    set_cp_raise_nrm(0);

    set_cp_size_nrm(0.5);

    m_ren_data.m_roi          =
        rect_t(cellid_t::zero,(m_ren_data.m_size-cellid_t::one)*2);

    m_ren_data.m_roi_base_pt  =
        ((m_ren_data.m_roi.upper_corner() +  m_ren_data.m_roi.lower_corner())/2);

    for(uint i = 1 ;i < m_gdm->m_pieces.size();++i)
      m_grid_piece_rens.push_back(new octtree_piece_rendata(m_gdm->m_pieces.at(i)));
  }

  viewer_t::~viewer_t()
  {
    for ( uint i = 1 ; i < m_grid_piece_rens.size();i++ )
      delete m_grid_piece_rens[i];

    m_grid_piece_rens.clear();

    glutils::clear();

    delete m_gdm;
  }

  inline glutils::vertex_t cell_to_vertex(cellid_t c,cell_fn_t fn =0)
  {
    return glutils::vertex_t(c[0],fn,c[1]);
  }

  void viewer_t::set_roi_dim_range_nrm(double l,double u,int dim)
  {
    if(!(l<u && 0.0 <= l && u <=1.0 && 0<=dim && dim < gc_grid_dim))
      return;

    rect_t roi = rect_t(cellid_t::zero,(m_ren_data.m_size-cellid_t::one)*2);

    double span = roi[dim].span();

    m_ren_data.m_roi[dim][0]  = (uint)(l*span);
    m_ren_data.m_roi[dim][1]  = (uint)(u*span);

    m_ren_data.m_roi_base_pt  = ((m_ren_data.m_roi.upper_corner() +
                                  m_ren_data.m_roi.lower_corner())/2);
  }

  int viewer_t::render()
  {
    if(m_bRebuildRens)
    {
      build_rens();

      m_bRebuildRens = false;
    }

    glPushAttrib(GL_ENABLE_BIT);

    glEnable(GL_NORMALIZE);

//    glScalef
//        (m_ren_data.m_scale_factor,0.125,
//         m_ren_data.m_scale_factor);

//    if(m_bCenterToRoi)
//      glTranslatef
//          (-m_ren_data.m_roi_base_pt[0],0,
//           -m_ren_data.m_roi_base_pt[1]);
//    else
//      glTranslatef
//          (std::min(-m_ren_data.m_size[0]+1,-1),0,
//           std::min(-m_ren_data.m_size[1]+1,-1));

    glScalef(1.0/(2*m_ren_data.m_size[1]),
             1.0,
             1.0/(2*m_ren_data.m_size[1]));

    glTranslatef(-(2*m_ren_data.m_size[0]-1)/2,0,-(2*m_ren_data.m_size[1]-1)/2);

    if(m_bShowSurface)
    {
      glPushAttrib(GL_ENABLE_BIT|GL_POLYGON_BIT);

      glColor3f(0.65,0.65,0.65);

      m_surf_ren->render();

#ifdef VIEWER_RENDER_AWESOME

      glDisable(GL_LIGHTING);

      glColor3f(0.15,0.15,0.15);

      glPolygonMode ( GL_FRONT_AND_BACK, GL_LINE );

      glPushMatrix();

      glLineWidth(4.0);

      glTranslatef(0,m_ren_data.m_cp_raise/4,0);

      m_surf_ren->render();

      glLineWidth(1.0);

      glPopMatrix();
#endif

      glPopAttrib();
    }

    if(m_bShowRoiBB)
    {
      glPushAttrib(GL_ENABLE_BIT);

      glDisable(GL_LIGHTING);

      glColor3dv(g_roiaabb_color.data());

      glutils::draw_aabb_line
          (cell_to_vertex(m_ren_data.m_roi.lower_corner()),
           cell_to_vertex(m_ren_data.m_roi.upper_corner()));

      glPopAttrib();
    }

    for ( uint i = 0 ; i < m_grid_piece_rens.size();i++ )
    {
      m_grid_piece_rens[i]->render_msgraph_data(m_ren_data);
      m_grid_piece_rens[i]->render_dataset_data(m_ren_data);
    }

    glPopAttrib();
  }

  void viewer_t::set_cp_raise_nrm(double r)
  {
    m_ren_data.m_cp_raise = r*g_max_cp_raise;
  }

  void viewer_t::set_cp_size_nrm(double s)
  {
    m_ren_data.m_cp_size = s*g_max_cp_size;
#ifdef VIEWER_RENDER_AWESOME
    m_ren_data.m_cp_size *= m_ren_data.m_size[1];
#endif
  }

  void viewer_t::build_rens()
  {
    for ( uint i = 0 ; i < m_grid_piece_rens.size();i++ )
    {
      m_grid_piece_rens[i]->create_cp_rens(m_ren_data);
      m_grid_piece_rens[i]->create_canc_cp_rens(m_ren_data);
      m_grid_piece_rens[i]->create_grad_rens(m_ren_data);
    }
  }

  void viewer_t::init()
  {
    glutils::init();

    octtree_piece_rendata::init();

    init_surf_ren();

    m_ren_data.m_scale_factor =
        0.5/ std::max((double) *std::max_element
                      (m_ren_data.m_size.begin(),m_ren_data.m_size.end())-1.0,1.0);

    for ( uint i = 0 ; i < m_grid_piece_rens.size();i++ )
    {
      m_grid_piece_rens[i]->create_disc_rds();
    }
  }

  void viewer_t::init_surf_ren()
  {
    using namespace glutils;

    rect_size_t s = m_gdm->m_size;

    varray_t vert_fns(boost::extents[s[0]][s[1]],boost::fortran_storage_order());

    std::ifstream ifs;

    ifs.open(m_elevation_filename.c_str(),std::ifstream::binary );

    if(ifs.is_open())
    {
      ifs.read ( reinterpret_cast<char *> ( vert_fns.data()),
                 sizeof ( cell_fn_t)*s[0]*s[1]);

      ifs.close();
    }
    else
    {
      memset(vert_fns.data(),0,sizeof ( cell_fn_t)*s[0]*s[1]);
    }

    s = s*2 - cellid_t::one;

    vertex_list_t vlist;

    for (cell_coord_t y = 0; y < s[1];y += 1)
      for (cell_coord_t x = 0; x < s[0];x += 1)
      {
        cell_fn_t fn = 0;
        fn += vert_fns(cellid_t(x + x%2,y + y%2)/2);
        fn += vert_fns(cellid_t(x - x%2,y + y%2)/2);
        fn += vert_fns(cellid_t(x - x%2,y - y%2)/2);
        fn += vert_fns(cellid_t(x + x%2,y - y%2)/2);

        vlist.push_back(vertex_t(x,fn/4.0,y));
      }

    quad_idx_list_t qlist;

    for (cell_coord_t y = 1; y < s[1];y += 2)
      for (cell_coord_t x = 1; x < s[0];x += 2)
      {
        qlist.push_back
            (quad_idx_t
             (m_ren_data.cellid_to_index(cellid_t(x - 1,y - 1)),
              m_ren_data.cellid_to_index(cellid_t(x - 1,y + 1)),
              m_ren_data.cellid_to_index(cellid_t(x + 1,y + 1)),
              m_ren_data.cellid_to_index(cellid_t(x + 1,y - 1)) ));
      }

    normal_list_t nlist;

    compute_vertex_normals(vlist,qlist,nlist);

    for (cell_coord_t y = 0; y < s[1];y += 1)
      for (cell_coord_t x = 0; x < s[0];x += 1)
      {
        normal_t n = normal_t::zero;

        n += nlist[m_ren_data.cellid_to_index(cellid_t(x + x%2,y + y%2))];
        n += nlist[m_ren_data.cellid_to_index(cellid_t(x + x%2,y - y%2))];
        n += nlist[m_ren_data.cellid_to_index(cellid_t(x - x%2,y + y%2))];
        n += nlist[m_ren_data.cellid_to_index(cellid_t(x - x%2,y - y%2))];

        nlist[m_ren_data.cellid_to_index(cellid_t(x,y))] = n/4;
      }

    m_ren_data.m_cell_bo     = make_buf_obj(vlist);
    m_ren_data.m_cell_nrm_bo = make_buf_obj(nlist);
    m_ren_data.m_size        = m_gdm->m_size;

    m_surf_ren.reset
        (create_buffered_quads_ren
         (m_ren_data.m_cell_bo,make_buf_obj(qlist),
          m_ren_data.m_cell_nrm_bo)
         );
  }

  configurable_t::data_index_t viewer_t::dim()
  {
    return data_index_t(11,m_grid_piece_rens.size());
  }
  bool viewer_t::exchange_field(const data_index_t &i,boost::any &v)
  {
    octtree_piece_rendata * dprd = m_grid_piece_rens[i[1]];

    switch(i[0])
    {
    case 0: return s_exchange_data_ro(dprd->dp->label(),v);
    case 1: return s_exchange_data_rw(dprd->m_bShowAllCps,v);
    case 2: return s_exchange_data_rw(dprd->m_bShowCps[0],v);
    case 3: return s_exchange_data_rw(dprd->m_bShowCps[1],v);
    case 4: return s_exchange_data_rw(dprd->m_bShowCps[2],v);
    case 5: return s_exchange_data_rw(dprd->m_bShowCpLabels,v);
    case 6: return s_exchange_data_rw(dprd->m_bShowMsGraph,v);
    case 7: return s_exchange_data_rw(dprd->m_bShowGrad[0],v);
    case 8: return s_exchange_data_rw(dprd->m_bShowGrad[1],v);
    case 9: return s_exchange_data_rw(dprd->m_bShowCancCps,v);
    case 10: return s_exchange_data_rw(dprd->m_bShowCancMsGraph,v);
    }

    throw std::logic_error("unknown index");
  }
  configurable_t::eFieldType viewer_t::exchange_header
      (const int &i,boost::any &v)
  {
    switch(i)
    {

    case 0: v =  std::string("oct tree piece"); return EFT_DATA_RO;
    case 1: v =  std::string("all cps");return EFT_DATA_RW;
    case 2: v =  std::string("minima");return EFT_DATA_RW;
    case 3: v =  std::string("1 saddle");return EFT_DATA_RW;
    case 4: v =  std::string("maxima");return EFT_DATA_RW;
    case 5: v =  std::string("cp labels");return EFT_DATA_RW;
    case 6: v =  std::string("msgraph");return EFT_DATA_RW;
    case 7: v =  std::string("gradient 0 ");return EFT_DATA_RW;
    case 8: v =  std::string("gradient 1 ");return EFT_DATA_RW;
    case 9: v =  std::string("cancelled cps");return EFT_DATA_RW;
    case 10: v =  std::string("cancelled cp msgraph");return EFT_DATA_RW;
    }

    throw std::logic_error("unknown index");

  }

  octtree_piece_rendata::octtree_piece_rendata (datapiece_t * _dp):
      m_bShowAllCps(false),
      m_bShowCpLabels ( false ),
      m_bShowMsGraph ( false ),
      m_bShowCancCps(false),
      m_bShowCancMsGraph(false),
      m_bNeedUpdateDiscRens(false),
      dp(_dp)
  {
    using namespace boost::lambda;

    std::for_each(m_bShowCps,m_bShowCps+gc_grid_dim+1,_1 = false);
    std::for_each(m_bShowGrad,m_bShowGrad+gc_grid_dim,_1 = false);

  }

  void  octtree_piece_rendata::create_cp_rens(const grid_ren_data_t& grd)
  {
    const rect_t roi = grd.m_roi;

    if(dp->msgraph == NULL)
      return;

    std::vector<glutils::point_idx_t>   crit_pt_idxs[gc_grid_dim+1];
    std::vector<glutils::line_idx_t>    crit_conn_idxs[gc_grid_dim];

    for(uint i = 0; i < dp->msgraph->m_cps.size(); ++i)
    {
      if(dp->msgraph->m_cps[i]->is_paired)
        continue;

      cellid_t c = (dp->msgraph->m_cps[i]->cellid);

      if(!roi.contains(c))
        continue;

      uint index = dp->msgraph->m_cps[i]->index;

      if(!dp->msgraph->m_cps[i]->is_paired)
      {
        crit_pt_idxs[index].push_back(grd.cellid_to_index(c));
      }
    }

    for(uint i = 0 ; i < gc_grid_dim+1; ++i)
    {
      ren_cp[i].reset(glutils::create_buffered_points_ren
                      (grd.m_cell_bo,
                       glutils::make_buf_obj(crit_pt_idxs[i])));
    }

    for(uint i = 0 ; i < dp->msgraph->m_cps.size(); ++i)
    {
      if(dp->msgraph->m_cps[i]->isCancelled)
        continue;

      if(dp->msgraph->m_cps[i]->is_paired)
        continue;

      cellid_t c = (dp->msgraph->m_cps[i]->cellid);

      if(!roi.contains(c))
        continue;

      uint index = dp->msgraph->m_cps[i]->index;

      for(conn_iter_t it  = dp->msgraph->m_cps[i]->conn[0].begin();
      it != dp->msgraph->m_cps[i]->conn[0].end(); ++it)
      {
        cellid_t cc =  dp->msgraph->m_cps[*it]->cellid;

        if(!roi.contains(cc)) continue;

        crit_conn_idxs[index-1].push_back
            (glutils::line_idx_t(grd.cellid_to_index(c),grd.cellid_to_index(cc)));
      }
    }

    for(uint i = 0 ; i < gc_grid_dim; ++i)
    {
      ren_cp_conns[i].reset(glutils::create_buffered_lines_ren
                            (grd.m_cell_bo,
                             glutils::make_buf_obj(crit_conn_idxs[i])));
    }

  }

  void  octtree_piece_rendata::create_canc_cp_rens(const grid_ren_data_t& grd)
  {
    if(dp->msgraph == NULL)
      return;

    const rect_t roi = grd.m_roi;

    std::vector<glutils::point_idx_t>   canc_cp_idxs[gc_grid_dim+1];
    std::vector<glutils::line_idx_t>    canc_cp_conn_idxs[gc_grid_dim];

    for(uint i = 0; i < dp->msgraph->m_cps.size(); ++i)
    {
      if(!dp->msgraph->m_cps[i]->is_paired)
        continue;

      cellid_t c = (dp->msgraph->m_cps[i]->cellid);

      if(!roi.contains(c))
        continue;

      uint index = dp->msgraph->m_cps[i]->index;

      canc_cp_idxs[index].push_back(grd.cellid_to_index(c));
    }

    for(uint i = 0 ; i < gc_grid_dim+1; ++i)
    {
      ren_canc_cp[i].reset(
          glutils::create_buffered_points_ren
          (grd.m_cell_bo,glutils::make_buf_obj(canc_cp_idxs[i])));
    }

    for(uint i = 0 ; i < dp->msgraph->m_cps.size(); ++i)
    {
      if(!dp->msgraph->m_cps[i]->is_paired)
        continue;

      cellid_t c = (dp->msgraph->m_cps[i]->cellid);

      if(!roi.contains(c))
        continue;

      uint index = dp->msgraph->m_cps[i]->index;

      for(uint dir = 0 ; dir <2 ;++dir)
      {
        for(conn_iter_t it  = dp->msgraph->m_cps[i]->conn[dir].begin();
        it != dp->msgraph->m_cps[i]->conn[dir].end(); ++it)
        {
          cellid_t cc =  dp->msgraph->m_cps[*it]->cellid;

          if(!roi.contains(cc)) continue;

          canc_cp_conn_idxs[index-(dir^1)].push_back
              (glutils::line_idx_t(grd.cellid_to_index(c),grd.cellid_to_index(cc)));
        }
      }
    }

    for(uint i = 0 ; i < gc_grid_dim; ++i)
    {
      ren_canc_cp_conns[i].reset(glutils::create_buffered_lines_ren
                                 (grd.m_cell_bo,
                                  glutils::make_buf_obj(canc_cp_conn_idxs[i])));
    }

  }

  void octtree_piece_rendata::create_grad_rens(const grid_ren_data_t& grd)
  {
    using namespace glutils;

    if(dp->dataset == NULL)
      return;

    rect_t roi = grd.m_roi,r;

    if(!dp->dataset->get_ext_rect().intersection(roi,r))
      return;

    line_idx_list_t  pair_idxs[gc_grid_dim];

    static_assert(gc_grid_dim == 2 && "defined for 2-manifolds only");

    cellid_t c;

    for(c[1] = r[1][0] ; c[1] <= r[1][1]; ++c[1])
    {
      for(c[0] = r[0][0] ; c[0] <= r[0][1]; ++c[0])
      {
        uint dim = dp->dataset->getCellDim(c);

        if(dp->dataset->isCellPaired(c))
        {
          cellid_t p = dp->dataset->getCellPairId(c);

          if(dp->dataset->isPairOrientationCorrect(c,p))
          {
            pair_idxs[dim].push_back
                (line_idx_t(grd.cellid_to_index(c),grd.cellid_to_index(p)));
          }
        }
      }
    }

    for(uint i = 0 ; i < gc_grid_dim; ++i)

    {
      ren_grad[i].reset(create_buffered_lines_ren
                        (grd.m_cell_bo,
                         make_buf_obj(pair_idxs[i])));
    }
  }

  void octtree_piece_rendata::create_disc_rds()
  {
    if(dp->msgraph == NULL)
      return;

    boost::shared_ptr<disc_rendata_t> sptr;

    for(uint i = 0 ; i < dp->msgraph->m_cps.size();++i)
    {
      critpt_t * cp = dp->msgraph->m_cps[i];

      if(cp->is_paired) continue;

      sptr.reset(new disc_rendata_t(cp->cellid,cp->index));

      disc_rds.push_back(sptr);
    }
  }

  void octtree_piece_rendata::update_active_disc_rens(const grid_ren_data_t& grd)
  {
    if(dp->msgraph == NULL)
      return;

    for(uint i = 0 ; i < disc_rds.size();++i)
    {
      if(disc_rds[i]->update(dp->msgraph,grd))
      {
        active_disc_rens.insert(disc_rds[i]);
      }
      else
      {
        active_disc_rens.erase(disc_rds[i]);
      }
    }
  }

  void octtree_piece_rendata::render_msgraph_data
      (const grid_ren_data_t& grd)
  {
    glPushMatrix();

    glPushAttrib ( GL_ENABLE_BIT );

    glDisable ( GL_LIGHTING );

    glTranslatef(0,grd.m_cp_raise,0);
#ifndef VIEWER_RENDER_AWESOME
    glPointSize ( grd.m_cp_size );

    glEnable(GL_POINT_SMOOTH);
#else
    s_sphere_shader->use();

    s_sphere_shader->sendUniform("g_wc_radius",(float)grd.m_cp_size);
#endif
    for(uint i = 0 ; i < gc_grid_dim+1;++i)
    {
      if(ren_cp[i]&& (m_bShowCps[i]||m_bShowAllCps))
      {
        glColor3dv(g_grid_cp_colors[i].data());

        ren_cp[i]->render();

        if(ren_cp_labels[i] && m_bShowCpLabels)
          ren_cp_labels[i]->render();
      }
    }

#ifdef VIEWER_RENDER_AWESOME
    s_sphere_shader->sendUniform("g_wc_radius",(float)grd.m_cp_size*2/3);
#endif

    if ( m_bShowCancCps)
    {
      for(uint i = 0 ; i < gc_grid_dim+1;++i)
      {
        if(ren_canc_cp[i])
        {
          glColor3dv(g_grid_cp_colors[i].data());

          ren_canc_cp[i]->render();

          if(ren_canc_cp_labels[i] &&
             m_bShowCpLabels)
            ren_canc_cp_labels[i]->render();
        }
      }
    }
#ifdef VIEWER_RENDER_AWESOME
    s_sphere_shader->disable();
#endif
    if (m_bShowMsGraph)
    {
      for(uint i = 0 ; i < gc_grid_dim;++i)
      {
        if(ren_cp_conns[i])
        {
          glColor3dv(g_grid_cp_conn_colors[i].data());

          ren_cp_conns[i]->render();
        }
      }
    }

    if (m_bShowCancMsGraph)
    {
      for(uint i = 0 ; i < gc_grid_dim;++i)
      {
        if(ren_canc_cp_conns[i])
        {
          glColor3dv(g_grid_cp_conn_colors[i].data());

          ren_canc_cp_conns[i]->render();
        }
      }
    }

    glPopAttrib();
    glPopMatrix();
  }

  void octtree_piece_rendata::render_dataset_data(const grid_ren_data_t& grd)
  {
    if(m_bNeedUpdateDiscRens)
    {
      update_active_disc_rens(grd);

      m_bNeedUpdateDiscRens = false;
    }
    glPushMatrix();

    glPushAttrib ( GL_ENABLE_BIT );

    for(disc_rendata_sp_set_t::iterator it = active_disc_rens.begin();
        it != active_disc_rens.end() ; ++it)
    {
      (*it)->render(grd);
    }
    glTranslatef(0,grd.m_cp_raise/2,0);
#ifdef VIEWER_RENDER_AWESOME
    s_grad_shader->use();

    grd.m_cell_nrm_bo->bind_to_normal_pointer();
#endif
    for(uint i = 0 ; i < gc_grid_dim; ++i)
    {
      if(ren_grad[i] && m_bShowGrad[i])
      {
        glColor3dv ( g_grid_grad_colors[i].data() );

        ren_grad[i]->render();
      }
    }
#ifdef VIEWER_RENDER_AWESOME
    grd.m_cell_nrm_bo->unbind_from_normal_pointer();

    s_grad_shader->disable();
#endif
    glPopAttrib();

    glPopMatrix();
  }

  struct random_color_assigner
  {
    disc_rendata_ptr_t m_drd;

    int m_no;

    static const uint MAX_RAND = 256;

    random_color_assigner(disc_rendata_ptr_t drd,int no):m_drd(drd),m_no(no){}

    void operator()()
    {
      for(uint c = 0 ; c < 3 ; ++c)
        m_drd->color[m_no][c] = ((double) (rand()%MAX_RAND))/((double)MAX_RAND);

    }
  };

  configurable_t::data_index_t octtree_piece_rendata::dim()
  {
    return data_index_t(8,disc_rds.size());
  }
  bool octtree_piece_rendata::exchange_field(const data_index_t &i,boost::any &v)
  {
    disc_rendata_ptr_t drd = disc_rds[i[1]];

    switch(i[0])
    {
    case 0:
      return s_exchange_data_ro(drd->cellid.to_string(),v);
    case 1:
      return s_exchange_data_ro((int)drd->index,v);
    case 2:
    case 3:
      {
        bool need_update = false;

        bool is_read     = v.empty();

        need_update =  s_exchange_data_rw(drd->show[i[0]%2],v);

        if(need_update && is_read == false )
          m_bNeedUpdateDiscRens = true;

        return need_update;
      }
    case 4:
    case 5:
      return s_exchange_data_rw(drd->color[i[0]%2],v);
    case 6:
    case 7:
      return s_exchange_action(random_color_assigner(drd,i[0]%2),v);

    };

    throw std::logic_error("invalid index");
  }
  configurable_t::eFieldType octtree_piece_rendata::exchange_header
      (const int &i,boost::any &v)
  {
    switch(i)
    {
    case 0: v = std::string("cellid"); return EFT_DATA_RO;
    case 1: v = std::string("index"); return EFT_DATA_RO;
    case 2: v = std::string("des mfold"); return EFT_DATA_RW;
    case 3: v = std::string("asc mfold"); return EFT_DATA_RW;
    case 4: v = std::string("des mfold color"); return EFT_DATA_RW;
    case 5: v = std::string("asc mfold color"); return EFT_DATA_RW;
    case 6: v = std::string("rand des mflod color"); return EFT_ACTION;
    case 7: v = std::string("rand asc mflod color"); return EFT_ACTION;
    }
    throw std::logic_error("invalid index");
  }

  disc_rendata_t::disc_rendata_t(cellid_t c,uint i):cellid(c),index(i)
  {
    color[0] = g_disc_colors[1][index];
    color[1] = g_disc_colors[0][index];

    show[0] =false; show[1] =false;
  }

  disc_rendata_t::~disc_rendata_t()
  {
    show[0] =false;
    show[1] =false;

  }

  void disc_rendata_t::render(const grid_ren_data_t &grd)
  {
    for(uint dir = 0 ; dir<2;++dir)
    {
      if(show[dir] && ren[dir] != NULL)
      {
#ifdef VIEWER_RENDER_AWESOME
        if(index == 1)
        {
          s_cylinder_shader->use();

          s_cylinder_shader->sendUniform("ug_cylinder_radius",(float)grd.m_cp_size/3);
        }
#endif
        glColor3dv(color[dir].data());

        ren[dir]->render();

#ifdef VIEWER_RENDER_AWESOME
        if(index == 1)
        {
          s_cylinder_shader->disable();

          s_sphere_shader->use();

          s_sphere_shader->sendUniform("g_wc_radius",(float)grd.m_cp_size/3);

          endpt_ren[dir]->render();

          s_sphere_shader->disable();

        }
#endif
      }
    }
  }

  bool disc_rendata_t::update(mscomplex_t *msc,const grid_ren_data_t& grd)
  {

    using namespace glutils;

    critpt_t *cp = msc->m_cps[msc->m_id_cp_map[cellid]];

    for(uint dir = 0 ; dir<2;++dir)
    {
      if(show[dir] == false && this->ren[dir] != NULL ) ren[dir].reset();

      if(show[dir] == false || this->ren[dir] != NULL || msc==NULL) continue;

      std::set<cellid_t> vset;

      for(uint j = 0 ; j < cp->contrib[dir].size();++j)
      {
        critpt_t *cp_contrib = msc->m_cps[cp->contrib[dir][j]];

        if(cp_contrib->index != cp->index)
          throw std::logic_error("contrib and cp must have same idx");

        for(uint i = 0; i < cp_contrib->disc[dir].size(); ++i)
        {
          cellid_t c = cp_contrib->disc[dir][i];

          if(vset.count(c) == 0)
            vset.insert(c);
        }
      }

      if(dir == 0 && cp->index == 2)
      {
        quad_idx_list_t qlist;

        for(std::set<cellid_t>::iterator it = vset.begin();it !=vset.end();++it)
        {
          cellid_t v[4];

          bool add_quad = true;

          quad_idx_t q;

          static const cellid_t v_order[] =
          {cellid_t(-1,-1),cellid_t(-1,+1),cellid_t(+1,+1),cellid_t(+1,-1)};

          for(int i = 0 ;i < 4;++i)
            v[i] = *it + v_order[i];

          for(int i = 0 ;i < 4;++i)
            add_quad &= grd.m_roi.contains(v[i]);

          if(add_quad == false)
            continue;

          for(int i = 0 ;i < 4;++i)
            q[i] = grd.cellid_to_index(v[i]);

          qlist.push_back(q);
        }

        ren[dir].reset(create_buffered_quads_ren
                       (grd.m_cell_bo,make_buf_obj(qlist),grd.m_cell_nrm_bo));
      }
      else if(dir == 0 && cp->index == 1)
      {
        line_idx_list_t llist;

        for(std::set<cellid_t>::iterator it = vset.begin();it !=vset.end();++it)
        {
          cellid_t v[2];

          v[0] = *it - (*it)%2;
          v[1] = *it + (*it)%2;

          if(grd.m_roi.contains(v[0]) && grd.m_roi.contains(v[1]) )
            llist.push_back
                (line_idx_t(grd.cellid_to_index(v[0]),
                            grd.cellid_to_index(v[1])));
        }

        bufobj_ptr_t lbo = make_buf_obj(llist);

        ren[dir].reset(create_buffered_lines_ren
                   (grd.m_cell_bo,lbo));

#ifdef VIEWER_RENDER_AWESOME
        endpt_ren[dir].reset(create_buffered_points_ren
                             (grd.m_cell_bo,recast_buf_obj_num_components(lbo,1)));
#endif
      }

      else if(dir == 1 && cp->index == 1)
      {
        line_idx_list_t llist;

        for(std::set<cellid_t>::iterator it = vset.begin();it !=vset.end();++it)
        {
          cellid_t v[2];

          v[0] = *it - (*it + cellid_t::one)%2;
          v[1] = *it;
          v[2] = *it + (*it + cellid_t::one)%2;

          if(grd.m_roi.contains(v[0]) && grd.m_roi.contains(v[1]) )
            llist.push_back
                (line_idx_t(grd.cellid_to_index(v[0]),
                            grd.cellid_to_index(v[1])));

          if(grd.m_roi.contains(v[1]) && grd.m_roi.contains(v[2]) )
            llist.push_back
                (line_idx_t(grd.cellid_to_index(v[1]),
                            grd.cellid_to_index(v[2])));

        }

        ren[dir].reset(create_buffered_lines_ren
                   (grd.m_cell_bo,make_buf_obj(llist)));

        bufobj_ptr_t lbo = make_buf_obj(llist);

        ren[dir].reset(create_buffered_lines_ren
                   (grd.m_cell_bo,lbo));

#ifdef VIEWER_RENDER_AWESOME
        endpt_ren[dir].reset(create_buffered_points_ren
                             (grd.m_cell_bo,recast_buf_obj_num_components(lbo,1)));
#endif

      }
      else if(dir == 1 && cp->index == 0)
      {
        quad_idx_list_t qlist;

        for(std::set<cellid_t>::iterator it = vset.begin();it !=vset.end();++it)
        {
          static const cellid_t v_order[] =
          {cellid_t(0,0),cellid_t(0,+1),cellid_t(+1,+1),cellid_t(+1,0)};

          for(int j = 0 ; j < 4;++j)
          {
            cellid_t c = *it - v_order[j];

            cellid_t v[4];

            bool add_quad = true;

            quad_idx_t q;

            for(int i = 0 ;i < 4;++i)
              v[i] = c + v_order[i];

            for(int i = 0 ;i < 4;++i)
              add_quad &= grd.m_roi.contains(v[i]);

            if(add_quad == false)
              continue;

            for(int i = 0 ;i < 4;++i)
              q[i] = grd.cellid_to_index(v[i]);

            qlist.push_back(q);
          }

        }

        ren[dir].reset(create_buffered_quads_ren
                       (grd.m_cell_bo,make_buf_obj(qlist),grd.m_cell_nrm_bo));
      }
    }

    return (show[0] || show[1]);
  }

}
