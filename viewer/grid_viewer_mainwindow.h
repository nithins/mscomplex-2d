#ifndef GRID_VIEWER_MAINWINDOW_INCLUDED
#define GRID_VIEWER_MAINWINDOW_INCLUDED

#include <grid.h>

#include <QDialog>
#include <QAbstractItemModel>
#include <QModelIndex>
#include <QVariant>
#include <QItemSelectionModel>
#include <QTreeView>
#include <QSortFilterProxyModel>
#include <QTimer>

#include <QGLViewer/qglviewer.h>

#include <boost/any.hpp>

class configurable_t;

namespace grid
{


  class viewer_t;

  class data_manager_t;

  class glviewer_t : public QGLViewer
  {

  public:

    viewer_t *m_ren;

    bool m_is_recording;
    bool m_bf_cull;
    bool m_wireframe;

    glviewer_t(QWidget *par);

    void setup(data_manager_t *dm,std::string ef);

    ~glviewer_t();

  protected:

    virtual void draw();
    virtual void init();
    virtual QString helpString() const;
    virtual void keyPressEvent(QKeyEvent *e);

  };

  class configurable_item_model;

}

#include <ui_grid_viewer_mainwindow.h>

namespace grid
{

  class viewer_mainwindow:
      public QDialog,
      public Ui::grid_viewer_mainwindow_Dialog
  {

  Q_OBJECT

  public:

    uint                     m_active_otp_idx;

    QSortFilterProxyModel   *m_cp_model_proxy;
    configurable_item_model *m_cp_model;
    configurable_item_model *m_otp_model;
    QTimer                  *m_clear_roi_aabb_timer;

  public:

    viewer_mainwindow(data_manager_t *gdm,std::string ef);

    ~viewer_mainwindow();

    void update_roi_box(double l,double u,uint dim);


  public:

    virtual void showEvent ( QShowEvent * );

  private slots:
    void on_datapiece_view_customContextMenuRequested ( const QPoint &p );

    void on_critpt_view_customContextMenuRequested ( const QPoint &p );

    void on_datapiece_view_activated ( const QModelIndex & index  );

    void on_xroi_spanslider_spanChanged(int l , int u );

    void on_yroi_spanslider_spanChanged(int l , int u );

    void on_cp_point_size_horizontalSlider_valueChanged ( int value );

    void on_cp_point_raise_horizontalSlider_valueChanged ( int value );

    void on_update_roi_pushButton_clicked(bool);

    void on_center_to_roi_checkBox_clicked(bool);

    void on_show_surface_checkBox_clicked(bool);

  private slots:
    void clear_roi_aabb();
  };

  void configurable_ctx_menu
      (configurable_t *c,
       const QModelIndexList & l,
       const QPoint &p);

  class configurable_ctx_menu_sig_collector:public QObject
  {
    Q_OBJECT

  public:

    std::vector<boost::any>  m_vals;
    int                      m_col;
    configurable_t *         m_conf;
    const std::vector<int> & m_rows;

    configurable_ctx_menu_sig_collector
        (configurable_t * conf,
         const std::vector<boost::any> & vals,
         const int & col,
         const std::vector<int> & rows,
         QObject *par):
        m_conf(conf),
        m_vals(vals),
        m_col(col),
        m_rows(rows)
    {setParent(par);}

  private slots:
    void triggered(bool state);
  };

  class configurable_item_model : public QAbstractTableModel
  {
    Q_OBJECT

  public:

    enum eColumnFilter
    {CF_EFT_DATA_RO = 1,
     CF_EFT_DATA_RW = 2,
     CF_EFT_ACTION  = 4};

    configurable_item_model ( configurable_t *conf,QObject *parent = 0 ):
        QAbstractTableModel ( parent ),m_conf(conf)
    {
      setColumnFilter(CF_EFT_DATA_RO|CF_EFT_DATA_RW);
    }

    QVariant data ( const QModelIndex &index, int role ) const;

    QVariant headerData ( int section, Qt::Orientation orientation,
                          int role = Qt::DisplayRole ) const;

    int rowCount ( const QModelIndex &parent = QModelIndex() ) const;

    int columnCount ( const QModelIndex &parent = QModelIndex() ) const;

    void reset_configurable(configurable_t *conf);

    void force_reset(){setColumnFilter(m_column_filter);reset();}

    void setColumnFilter(int columnFilter);

  private:

    int m_column_filter;

    std::vector<int> m_column_idxs;

    configurable_t * m_conf;

  };
}



#endif
