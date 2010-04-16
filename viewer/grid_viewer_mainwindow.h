#ifndef GRID_VIEWER_MAINWINDOW_INCLUDED
#define GRID_VIEWER_MAINWINDOW_INCLUDED

#include <QDialog>
#include <QAbstractItemModel>
#include <QModelIndex>
#include <QVariant>
#include <QItemSelectionModel>

#include <ui_grid_viewer_mainwindow.h>

class grid_glviewer;

class grid_piece_rendata;

class GridDataPiece;

class grid_viewer_mainwindow:
    public QDialog,
    public Ui::grid_viewer_mainwindow_Dialog
{

public:

  grid_glviewer    *m_viewer;

  grid_viewer_mainwindow
      (std::vector<GridDataPiece *> * p ,uint size_x,uint size_y);

  Q_OBJECT

public:

  enum eTreeViewActions
  {
    TVA_SURF,
    TVA_CPS,
    TVA_CPLABELS,
    TVA_GRAPH,
    TVA_GRAD,
    TVA_CANC_CPS,
    TVA_CANC_GRAPH,
  };

  void perform_tva_action ( const eTreeViewActions &,const bool & );
  bool get_tva_state ( const eTreeViewActions & );


private slots:
  void on_datapiece_treeView_customContextMenuRequested ( const QPoint &p );

  void show_surf_toggled ( bool state ) {perform_tva_action ( TVA_SURF,state );}
  void show_cps_toggled ( bool state ) {perform_tva_action ( TVA_CPS,state );}
  void show_cplabels_toggled ( bool state ) {perform_tva_action ( TVA_CPLABELS,state );}
  void show_graph_toggled ( bool state ) {perform_tva_action ( TVA_GRAPH,state );}
  void show_grad_toggled ( bool state ) {perform_tva_action ( TVA_GRAD,state );}
  void show_canc_cps_toggled ( bool state ) {perform_tva_action ( TVA_CANC_CPS,state );}
  void show_canc_graph_toggled ( bool state ) {perform_tva_action ( TVA_CANC_GRAPH,state );}

};

class GridTreeModel : public QAbstractItemModel
{
    Q_OBJECT

  public:

  GridTreeModel ( std::vector<grid_piece_rendata *> *, QObject *parent = 0 );
    ~GridTreeModel();

    QVariant data ( const QModelIndex &index, int role ) const;

    Qt::ItemFlags flags ( const QModelIndex &index ) const;

    QVariant headerData ( int section, Qt::Orientation orientation,
                          int role = Qt::DisplayRole ) const;

    QModelIndex index ( int row, int column,
                        const QModelIndex &parent = QModelIndex() ) const;

    QModelIndex parent ( const QModelIndex &index ) const;

    int rowCount ( const QModelIndex &parent = QModelIndex() ) const;

    int columnCount ( const QModelIndex &parent = QModelIndex() ) const;

    struct tree_item
    {
      std::vector<tree_item *>               children;
      grid_piece_rendata                   * node;
      tree_item                            * parent;

      tree_item ( grid_piece_rendata * _node , tree_item * par );
      tree_item();
      int row();
    };


  private:
    void setupModelData
        ( std::vector<grid_piece_rendata *> *);

    tree_item *m_tree;
};

class RecursiveTreeItemSelectionModel:
    public QItemSelectionModel
{
    Q_OBJECT

  public:
    RecursiveTreeItemSelectionModel ( QAbstractItemModel * m,QTreeView * tv) :
        QItemSelectionModel ( m),m_pTreeView ( tv) {}

  public slots:

    virtual void select ( const QModelIndex &index,
                          QItemSelectionModel::SelectionFlags command );

    virtual void select ( const QItemSelection &selection,
                          QItemSelectionModel::SelectionFlags command );
  private:

    void collect_all_children ( const QModelIndex &index,
                                QModelIndexList &retlist );

    QTreeView *m_pTreeView;
};

#endif
