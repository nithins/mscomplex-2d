#ifndef GRID_VIEWER_H_INCLUDED
#define GRID_VIEWER_H_INCLUDED

#include <QGLViewer/qglviewer.h>

namespace glutils
{
  class renderable_t;
}

namespace grid
{

  typedef unsigned char uchar;
  typedef unsigned int  uint;

  class datapiece_t ;

  class grid_piece_rendata
  {
  public:

    datapiece_t * dp;

    // set externally to control what is rendered
    bool m_bShowSurface;
    bool m_bShowCps;
    bool m_bShowCpLabels;
    bool m_bShowMsGraph;
    bool m_bShowGrad;
    bool m_bShowCancCps;
    bool m_bShowCancMsGraph;

    glutils::renderable_t  *ren_surf;
    glutils::renderable_t  *ren_grad[2];
    glutils::renderable_t  *ren_cp_labels[3];
    glutils::renderable_t  *ren_cp[3];
    glutils::renderable_t  *ren_cp_conns[2];
    glutils::renderable_t  *ren_canc_cp_labels[3];
    glutils::renderable_t  *ren_canc_cp[3];
    glutils::renderable_t  *ren_canc_cp_conns[2];

    void create_cp_rens();
    void create_grad_rens();
    void create_surf_ren();
    void render() const ;

    grid_piece_rendata(datapiece_t *);
  };

  class grid_glviewer : public QGLViewer
  {
  public:

    std::vector<grid_piece_rendata * >  m_grid_piece_rens;

    uint              m_size_x;
    uint              m_size_y;

  public:
    grid_glviewer(std::vector<datapiece_t *> * p ,
                  uint size_x,uint size_y);
    ~grid_glviewer();

  protected :
      virtual void draw();
  virtual void init();
  virtual QString helpString() const;
};
}
#endif //VIEWER_H_INCLUDED
