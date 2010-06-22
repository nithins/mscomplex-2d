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

const int s_rawdata_texture_no = 0;


struct s_shader_data_t
{
  enum eSourceType { ST_VERT,ST_GEOM,ST_FRAG,ST_COUNT};

  GLSLProgram * prog;

  const char * source[ST_COUNT];
  const char * header[ST_COUNT];

  uint geom_in;
  uint geom_out;
};

const char * s_is_dual_str[grid::GRADDIR_COUNT]=
{
  "const int is_dual = 0;",
  "const int is_dual = 1;",
};

const char * s_shader_header_regex =
    "//HEADER_REPLACE_BEGIN(.*)//HEADER_REPLACE_END";

s_shader_data_t  s_cell_shaders[grid::GRADDIR_COUNT][grid::gc_grid_dim+1] =
{
  {
    {NULL,{vert_2mfold_glsl,geom_2mfold_glsl,NULL},{NULL,s_is_dual_str[0],NULL},GL_POINTS,GL_TRIANGLES},
    {NULL,{vert_1mfold_glsl,geom_1mfold_glsl,NULL},{NULL,s_is_dual_str[0],NULL},GL_POINTS,GL_LINE_STRIP},
    {NULL,{vert_2mfold_glsl,geom_2mfold_glsl,NULL},{NULL,s_is_dual_str[0],NULL},GL_POINTS,GL_TRIANGLES},
  },
{
    {NULL,{vert_2mfold_glsl,geom_2mfold_glsl,NULL},{NULL,s_is_dual_str[1],NULL},GL_POINTS,GL_TRIANGLES},
    {NULL,{vert_1mfold_glsl,geom_1mfold_glsl,NULL},{NULL,s_is_dual_str[1],NULL},GL_POINTS,GL_LINE_STRIP},
    {NULL,{vert_2mfold_glsl,geom_2mfold_glsl,NULL},{NULL,s_is_dual_str[1],NULL},GL_POINTS,GL_TRIANGLES},
  },
};

glutils::color_t g_grid_cp_colors[grid::gc_grid_dim+1] =
{
  glutils::color_t(0.0,0.0,1.0),
  glutils::color_t(0.0,1.0,0.0),
  glutils::color_t(1.0,0.0,0.0),
};

glutils::color_t g_grid_grad_colors[grid::gc_grid_dim] =
{
  glutils::color_t(0.0,0.5,0.5 ),
  glutils::color_t(0.5,0.0,0.5 ),
};

glutils::color_t g_disc_colors[grid::GRADDIR_COUNT][grid::gc_grid_dim+1] =
{
  {
    glutils::color_t(0.15,0.45,0.35 ),
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
  glutils::color_t(0.0,0.5,0.5 ),
  glutils::color_t(0.5,0.0,0.5 ),
};

glutils::color_t g_roiaabb_color = glutils::color_t(0.85,0.75,0.65);

double g_max_cp_size  = 8.0;
double g_max_cp_raise = 0.025;

GLSLProgram * s_grad_shader = NULL;

namespace grid
{

  glutils::vertex_t cell_to_vertex(cellid_t c)
  {
    return glutils::vertex_t(c[0],0,c[1]);
  }

  glutils::vertex_t cp_to_vertex(mscomplex_t *msc,uint i)
  {
    cellid_t &c = msc->m_cps[i]->cellid;

    return glutils::vertex_t(c[0],msc->m_cps[i]->fn,c[1]);
  }

  void disc_rendata_t::init()
  {
    for(uint k = 0 ;k < GRADDIR_COUNT;++k)
    {
      for(uint j = 0 ;j < gc_grid_dim+1;++j)
      {
        if(s_cell_shaders[k][j].prog != NULL )
          continue;

        std::string shader_source[s_shader_data_t::ST_COUNT];

        for(uint i = 0 ;i < s_shader_data_t::ST_COUNT;++i)
        {
          if(s_cell_shaders[k][j].source[i] != NULL)
            shader_source[i] = s_cell_shaders[k][j].source[i];

          if(s_cell_shaders[k][j].header[i] != NULL)
          {
            boost::replace_regex
                ( shader_source[i],boost::regex(s_shader_header_regex),
                  std::string(s_cell_shaders[k][j].header[i]));
          }
        }

        s_cell_shaders[k][j].prog =
            GLSLProgram::createFromSourceStrings
            (shader_source[0],shader_source[1],shader_source[2],
             s_cell_shaders[k][j].geom_in,s_cell_shaders[k][j].geom_out);

        std::string log;

        s_cell_shaders[k][j].prog->GetProgramLog ( log );

        if(log.size() !=0 )
        {
          std::stringstream ss;

          for(uint i = 0 ;i < s_shader_data_t::ST_COUNT;++i)
          {
            ss<<"shader no"<<i<<"\n";

            ss<<shader_source[i]<<"\n";
          }

          ss<<"shader log ::\n"<<log<<"\n";

          throw std::runtime_error(ss.str());
        }
      }
    }

  }

  void disc_rendata_t::cleanup()
  {
    for(uint k = 0 ;k < GRADDIR_COUNT;++k)
    {
      for(uint j = 0 ;j < gc_grid_dim+1;++j)
      {
        if(s_cell_shaders[k][j].prog == NULL )
          continue;

        delete s_cell_shaders[k][j].prog;

        s_cell_shaders[k][j].prog = NULL;
      }
    }
  }

  void octtree_piece_rendata::init()
  {
    if(s_grad_shader != NULL)
      return;

    s_grad_shader = GLSLProgram::createFromSourceStrings
                    (grad_vert_glsl,grad_geom_glsl,std::string(),
                     GL_LINES,GL_TRIANGLES);

    std::string log;

    s_grad_shader->GetProgramLog(log);

    if(log.size() != 0 )
      throw std::runtime_error("******grad_shader compile error*******\n"+log);

  }

  void octtree_piece_rendata::cleanup()
  {
    if(s_grad_shader == NULL)
      return;

    delete s_grad_shader;

    s_grad_shader = NULL;
  }


  viewer_t::viewer_t
      (data_manager_t * gdm,std::string ef):
      m_size(gdm->m_size),m_scale_factor(0),
      m_bRebuildRens(true),m_bShowRoiBB(false),m_bCenterToRoi(false),
      m_cp_raise(0.0),m_cp_size(0.5),
      m_gdm(gdm)
  {
    m_roi = rect_t(cellid_t::zero,(m_size-cellid_t::one)*2);

    for(uint i = 0 ;i < m_gdm->m_pieces.size();++i)
      m_grid_piece_rens.push_back(new octtree_piece_rendata(m_gdm->m_pieces.at(i)));

    m_roi_base_pt  = ((m_roi.upper_corner() +  m_roi.lower_corner())/2);

    m_elevation_filename = ef;

  }

  viewer_t::~viewer_t()
  {
    for ( uint i = 0 ; i < m_grid_piece_rens.size();i++ )
      delete m_grid_piece_rens[i];

    m_grid_piece_rens.clear();

    glutils::clear();

    disc_rendata_t::cleanup();

    octtree_piece_rendata::cleanup();

    glDeleteTextures( 1, &m_rawdata_texture );

    delete m_gdm;
  }

  void viewer_t::set_roi_dim_range_nrm(double l,double u,int dim)
  {
    if(!(l<u && 0.0 <= l && u <=1.0 && 0<=dim && dim < gc_grid_dim))
      return;

    rect_t roi = rect_t(cellid_t::zero,(m_size-cellid_t::one)*2);

    double span = roi[dim].span();

    m_roi[dim][0]  = (uint)(l*span);
    m_roi[dim][1]  = (uint)(u*span);

    m_roi_base_pt  = ((m_roi.upper_corner() +  m_roi.lower_corner())/2);
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

    glScalef(m_scale_factor,
             0.125,
             m_scale_factor);

    if(m_bCenterToRoi)
      glTranslatef(-m_roi_base_pt[0],0,-m_roi_base_pt[1]);
    else
      glTranslatef(std::min(-m_size[0]+1,-1),0,std::min(-m_size[1]+1,-1));

    if(m_bShowRoiBB)
    {
      glPushAttrib(GL_ENABLE_BIT);

      glDisable(GL_LIGHTING);

      glColor3dv(g_roiaabb_color.data());

      glutils::draw_aabb_line(cell_to_vertex(m_roi.lower_corner()),
                              cell_to_vertex(m_roi.upper_corner()));

      glPopAttrib();
    }

    for ( uint i = 0 ; i < m_grid_piece_rens.size();i++ )
    {
      m_grid_piece_rens[i]->render_msgraph_data
          (m_cp_raise*g_max_cp_raise/0.125,m_cp_size*g_max_cp_size);
    }

    if(m_rawdata_texture != 0 )
    {
      glEnable ( GL_TEXTURE_RECTANGLE_ARB );

      glActiveTexture ( GL_TEXTURE0 + s_rawdata_texture_no );
      glBindTexture ( GL_TEXTURE_RECTANGLE_ARB, m_rawdata_texture );
    }

    for ( uint i = 0 ; i < m_grid_piece_rens.size();i++ )
    {
      m_grid_piece_rens[i]->render_dataset_data(m_cp_raise*g_max_cp_raise/0.125);
    }

    if(m_rawdata_texture != 0 )
    {
      glBindTexture ( GL_TEXTURE_RECTANGLE_ARB, 0);

      glDisable ( GL_TEXTURE_RECTANGLE_ARB );
    }

    glPopAttrib();
  }

  void viewer_t::build_rens()
  {
    for ( uint i = 0 ; i < m_grid_piece_rens.size();i++ )
    {
      m_grid_piece_rens[i]->create_cp_rens(m_roi);
      m_grid_piece_rens[i]->create_canc_cp_rens(m_roi);
      m_grid_piece_rens[i]->create_grad_rens(m_roi);
    }
  }

  void viewer_t::init()
  {
    glutils::init();

    disc_rendata_t::init();

    octtree_piece_rendata::init();

    init_rawdata_texture();

    m_scale_factor = 0.5/ std::max((double) *std::max_element
                                   (m_size.begin(),m_size.end())-1.0,1.0);

    /*turn back face culling off */
    glEnable ( GL_CULL_FACE );

    /*cull backface */
    glCullFace ( GL_BACK );

    /*polymode */
    glPolygonMode ( GL_FRONT, GL_FILL );

    glPolygonMode ( GL_BACK, GL_LINE );

    for ( uint i = 0 ; i < m_grid_piece_rens.size();i++ )
    {
      m_grid_piece_rens[i]->create_cp_loc_bo();
      m_grid_piece_rens[i]->create_disc_rds();
    }
  }

  bool viewer_t::init_rawdata_texture()
  {
    std::ifstream ifs;

    ifs.open(m_elevation_filename.c_str(),std::ifstream::binary );

    if(!ifs.is_open())
      return false;

    cell_fn_t * pData = new cell_fn_t[m_gdm->m_size[0]*m_gdm->m_size[1]];

    ifs.read ( reinterpret_cast<char *> ( pData),
               sizeof ( cell_fn_t)*m_gdm->m_size[0]*m_gdm->m_size[1]);

    ifs.close();

    glGenTextures ( 1, &m_rawdata_texture );

    glBindTexture ( GL_TEXTURE_RECTANGLE_ARB, m_rawdata_texture );

    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexParameteri(GL_TEXTURE_RECTANGLE_ARB, GL_TEXTURE_WRAP_T, GL_CLAMP);

    glTexImage2D(GL_TEXTURE_RECTANGLE_ARB, 0, GL_FLOAT_R32_NV,
                 m_gdm->m_size[0], m_gdm->m_size[1], 0,
                 GL_RED, GL_FLOAT, pData);

    glBindTexture ( GL_TEXTURE_RECTANGLE_ARB, 0 );

    delete []pData;

    return true;
  }

  configurable_t::data_index_t viewer_t::dim()
  {
    return data_index_t(10,m_grid_piece_rens.size());
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
    case 7: return s_exchange_data_rw(dprd->m_bShowGrad,v);
    case 8: return s_exchange_data_rw(dprd->m_bShowCancCps,v);
    case 9: return s_exchange_data_rw(dprd->m_bShowCancMsGraph,v);
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
    case 7: v =  std::string("gradient");return EFT_DATA_RW;
    case 8: v =  std::string("cancelled cps");return EFT_DATA_RW;
    case 9: v =  std::string("cancelled cp msgraph");return EFT_DATA_RW;
    }

    throw std::logic_error("unknown index");

  }

  octtree_piece_rendata::octtree_piece_rendata (datapiece_t * _dp):
      m_bShowAllCps(false),
      m_bShowCpLabels ( false ),
      m_bShowMsGraph ( false ),
      m_bShowGrad ( false ),
      m_bShowCancCps(false),
      m_bShowCancMsGraph(false),
      m_bNeedUpdateDiscRens(false),
      dp(_dp)
  {
    using namespace boost::lambda;

    std::for_each(m_bShowCps,m_bShowCps+gc_grid_dim+1,_1 = false);

  }

  void octtree_piece_rendata::create_cp_loc_bo()
  {
    if(dp->msgraph == NULL)
      return;

    std::vector<glutils::vertex_t>  cp_loc;

    for(uint i = 0; i < dp->msgraph->m_cps.size(); ++i)
    {
      cp_loc.push_back(cp_to_vertex(dp->msgraph,i));
    }

    cp_loc_bo = glutils::make_buf_obj(cp_loc);
  }

  void  octtree_piece_rendata::create_cp_rens(const rect_t & roi)
  {
    if(dp->msgraph == NULL)
      return;

    std::vector<std::string>            crit_labels[gc_grid_dim+1];
    std::vector<glutils::vertex_t>      crit_label_locations[gc_grid_dim+1];
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

      std::stringstream ss;

      ((std::ostream&)ss)<<c;

      if(!dp->msgraph->m_cps[i]->is_paired)
      {
        crit_labels[index].push_back(ss.str());
        crit_label_locations[index].push_back(cp_to_vertex(dp->msgraph,i));
        crit_pt_idxs[index].push_back(i);
      }
    }

    for(uint i = 0 ; i < gc_grid_dim+1; ++i)
    {
      ren_cp_labels[i].reset(glutils::create_buffered_text_ren
                             (crit_labels[i],crit_label_locations[i]));

      ren_cp[i].reset(glutils::create_buffered_points_ren
                      (cp_loc_bo,
                       glutils::make_buf_obj(crit_pt_idxs[i]),
                       glutils::make_buf_obj()));
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
        if(!roi.contains(dp->msgraph->m_cps[*it]->cellid))
          continue;

        crit_conn_idxs[index-1].push_back
            (glutils::line_idx_t(i,*it));
      }
    }

    for(uint i = 0 ; i < gc_grid_dim; ++i)
    {
      ren_cp_conns[i].reset(glutils::create_buffered_lines_ren
                            (cp_loc_bo,
                             glutils::make_buf_obj(crit_conn_idxs[i]),
                             glutils::make_buf_obj()));
    }

  }

  void  octtree_piece_rendata::create_canc_cp_rens(const rect_t & roi)
  {
    if(dp->msgraph == NULL)
      return;

    std::vector<std::string>            canc_cp_labels[gc_grid_dim+1];
    std::vector<glutils::vertex_t>      canc_cp_label_locations[gc_grid_dim+1];
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



      canc_cp_labels[index].push_back(c.to_string());
      canc_cp_label_locations[index].push_back(cp_to_vertex(dp->msgraph,i)) ;
      canc_cp_idxs[index].push_back(i);

    }

    for(uint i = 0 ; i < gc_grid_dim+1; ++i)
    {

      ren_canc_cp_labels[i].reset(glutils::create_buffered_text_ren
                                  (canc_cp_labels[i],canc_cp_label_locations[i]));

      ren_canc_cp[i].reset(glutils::create_buffered_points_ren
                           (cp_loc_bo,
                            glutils::make_buf_obj(canc_cp_idxs[i]),
                            glutils::make_buf_obj()));
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
          if(!roi.contains(dp->msgraph->m_cps[*it]->cellid))
            continue;

          canc_cp_conn_idxs[index-(dir^1)].push_back
              (glutils::line_idx_t(i,*it));
        }
      }
    }

    for(uint i = 0 ; i < gc_grid_dim; ++i)
    {
      ren_canc_cp_conns[i].reset(glutils::create_buffered_lines_ren
                                 (cp_loc_bo,
                                  glutils::make_buf_obj(canc_cp_conn_idxs[i]),
                                  glutils::make_buf_obj()));
    }

  }

  void octtree_piece_rendata::create_grad_rens(const rect_t & roi)
  {
    if(dp->dataset == NULL)
      return;

    rect_t r;
    if(!dp->dataset->get_ext_rect().intersection(roi,r))
      return;

    std::vector<glutils::vertex_t>      cell_locations;
    std::vector<glutils::line_idx_t>    pair_idxs[gc_grid_dim];

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
            cell_locations.push_back(cell_to_vertex(c) );

            cell_locations.push_back(cell_to_vertex(p) );

            pair_idxs[dim].push_back
                (glutils::line_idx_t(cell_locations.size()-2,
                                     cell_locations.size()-1));
          }
        }
      }
    }

    glutils::bufobj_ptr_t cell_bo= glutils::make_buf_obj(cell_locations);

    for(uint i = 0 ; i < gc_grid_dim; ++i)

    {
      ren_grad[i].reset(glutils::create_buffered_lines_ren
                        (cell_bo,
                         glutils::make_buf_obj(pair_idxs[i]),
                         glutils::make_buf_obj()));
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

  void octtree_piece_rendata::update_active_disc_rens()
  {
    if(dp->msgraph == NULL)
      return;

    for(uint i = 0 ; i < disc_rds.size();++i)
    {
      if(disc_rds[i]->update(dp->msgraph))
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
      (double cp_raise,double cp_point_size)
  {
    glPushMatrix();
    glPushAttrib ( GL_ENABLE_BIT );

    glDisable ( GL_LIGHTING );

    glPointSize ( cp_point_size );

    glEnable(GL_POINT_SMOOTH);

    glTranslatef(0,cp_raise,0);

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

  void octtree_piece_rendata::render_dataset_data(double grad_raise)
  {
    if(m_bNeedUpdateDiscRens)
    {
      update_active_disc_rens();
      m_bNeedUpdateDiscRens = false;
    }

    glPushMatrix();
    glPushAttrib ( GL_ENABLE_BIT );

    for(disc_rendata_sp_set_t::iterator it = active_disc_rens.begin();
        it != active_disc_rens.end() ; ++it)
    {
      (*it)->render();
    }

    glTranslatef(0,grad_raise,0);

    s_grad_shader->use();

    if(m_bShowGrad)
    {
      for(uint i = 0 ; i < gc_grid_dim; ++i)
      {
        if(ren_grad[i])
        {
          glColor3dv ( g_grid_grad_colors[i].data() );

          ren_grad[i]->render();
        }
      }
    }

    s_grad_shader->disable();

    glPopAttrib();
    glPopMatrix();

  }

  struct random_color_assigner
  {
    disc_rendata_ptr_t m_drd;

    int m_no;

    static const uint MAX_RAND = 256;

    random_color_assigner(disc_rendata_ptr_t drd,int no):m_drd(drd),m_no(no){};

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

    show[0] =false; ren[0] =NULL;
    show[1] =false; ren[1] =NULL;
  }

  disc_rendata_t::~disc_rendata_t()
  {
    show[0] =false;
    show[1] =false;

    update(NULL);

  }

  void disc_rendata_t::render()
  {
    for(uint dir = 0 ; dir<2;++dir)
    {
      if(show[dir])
      {
        s_cell_shaders[dir][index].prog->use();

        s_cell_shaders[dir][index].prog->sendUniform ( "rawdata_texture",s_rawdata_texture_no );

        glColor3dv(g_grid_cp_colors[index].data());

        glBegin(GL_POINTS);
        glVertex3dv(cell_to_vertex(cellid).data());
        glEnd();

        glColor3dv(color[dir].data());

        ren[dir]->render();

        s_cell_shaders[dir][index].prog->disable();
      }
    }
  }

  bool disc_rendata_t::update(mscomplex_t *msc)
  {

    using namespace boost::lambda;

    for(uint dir = 0 ; dir<2;++dir)
    {
      if(show[dir] && this->ren[dir] == NULL && msc)
      {
        ensure_cellid_critical(msc,cellid);

        critpt_t *cp = msc->m_cps[msc->m_id_cp_map[cellid]];

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

        std::vector<glutils::vertex_t> vlist(vset.size());

        std::transform(vset.begin(),vset.end(),vlist.begin(),cell_to_vertex);

        ren[dir] = glutils::create_buffered_points_ren
                   (glutils::make_buf_obj(vlist),
                    glutils::make_buf_obj(),
                    glutils::make_buf_obj());

      }

      if(!show[dir] && this->ren[dir] != NULL )
      {
        delete ren[dir];

        ren[dir] = NULL;

      }
    }

    return (show[0] || show[1]);
  }

}
