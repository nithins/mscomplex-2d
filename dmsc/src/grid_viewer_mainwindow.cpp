#include <QMenu>
#include <QTreeView>


#include <grid_viewer.h>
#include <grid_viewer_mainwindow.h>
#include <grid_datamanager.h>


bool &get_tva_item_flag_ref
    (const grid_viewer_mainwindow::eTreeViewActions &action ,
     grid_piece_rendata *gp_rd )
{
  switch (action)
  {
  case grid_viewer_mainwindow::TVA_SURF: return gp_rd->m_bShowSurface;
  case grid_viewer_mainwindow::TVA_CPS: return gp_rd->m_bShowCps;
  case grid_viewer_mainwindow::TVA_CPLABELS: return gp_rd->m_bShowCpLabels;
  case grid_viewer_mainwindow::TVA_GRAPH:return gp_rd->m_bShowMsGraph;
  case grid_viewer_mainwindow::TVA_GRAD:return gp_rd->m_bShowGrad;
  case grid_viewer_mainwindow::TVA_CANC_CPS:return gp_rd->m_bShowCancCps;
  case grid_viewer_mainwindow::TVA_CANC_GRAPH:return gp_rd->m_bShowCancMsGraph;
  }

  throw std::invalid_argument("undefined tva action");

  return gp_rd->m_bShowSurface;
}

bool grid_viewer_mainwindow::get_tva_state ( const eTreeViewActions &action )
{

  uint num_checked_items = 0;
  uint num_unchecked_items = 0;

  QModelIndexList indexes = datapiece_treeView->selectionModel()->selectedIndexes();

  for ( QModelIndexList::iterator ind_it = indexes.begin();ind_it != indexes.end(); ++ind_it )
  {
    GridTreeModel::tree_item *item = static_cast<GridTreeModel::tree_item*> ( ( *ind_it ).internalPointer() );

    if ( get_tva_item_flag_ref(action,item->node) ) ++num_checked_items;
    else ++num_unchecked_items;
  }

  if ( num_checked_items > num_unchecked_items )
    return true;
  else
    return false;
}



void grid_viewer_mainwindow::perform_tva_action ( const eTreeViewActions &action, const bool & state )
{

  QModelIndexList indexes = datapiece_treeView->selectionModel()->selectedIndexes();

  bool need_update = false;

  for ( QModelIndexList::iterator ind_it = indexes.begin();ind_it != indexes.end(); ++ind_it )
  {
    GridTreeModel::tree_item *item = static_cast<GridTreeModel::tree_item*> ( ( *ind_it ).internalPointer() );

    if(state != get_tva_item_flag_ref(action,item->node))
    {
      get_tva_item_flag_ref(action,item->node) = state;
      need_update = true;
    }

  }

  if(need_update )
  {
    m_viewer->updateGL();
  }
}

#define ADD_MENU_ACTION(menu_ptr,action_string,start_state,recv_func_name) \
{\
 QAction * _ama_action  = (menu_ptr)->addAction ( tr(action_string) );\
                          _ama_action->setCheckable ( true );\
                          _ama_action->setChecked ( (start_state));\
                          connect ( _ama_action,SIGNAL ( toggled ( bool ) ),this,SLOT ( recv_func_name(bool) ) );\
                        }\



void grid_viewer_mainwindow::on_datapiece_treeView_customContextMenuRequested ( const QPoint &pos )
{

  QMenu menu;

  ADD_MENU_ACTION ( &menu, "show surface", get_tva_state ( TVA_SURF ), show_surf_toggled );
  ADD_MENU_ACTION ( &menu, "show critical points", get_tva_state ( TVA_CPS ), show_cps_toggled );
  ADD_MENU_ACTION ( &menu, "show critical point labels", get_tva_state ( TVA_CPLABELS ), show_cplabels_toggled );
  ADD_MENU_ACTION ( &menu, "show graph", get_tva_state ( TVA_GRAPH ), show_graph_toggled );
  ADD_MENU_ACTION ( &menu, "show grad", get_tva_state ( TVA_GRAD ), show_grad_toggled );
  ADD_MENU_ACTION ( &menu, "show canc cps", get_tva_state ( TVA_CANC_CPS ), show_canc_cps_toggled );
  ADD_MENU_ACTION ( &menu, "show canc graph", get_tva_state ( TVA_CANC_GRAPH ), show_canc_graph_toggled );
  menu.exec ( datapiece_treeView->mapToGlobal ( pos ) );

}



grid_viewer_mainwindow::grid_viewer_mainwindow
    (std::vector<GridDataPiece *> * p ,
     uint size_x,
     uint size_y)
{
  setupUi (this);

  m_viewer = new grid_glviewer(p,size_x,size_y);

  m_viewer->setParent(glviewer);

  m_viewer->resize(glviewer->size());

  GridTreeModel *model = new GridTreeModel ( &m_viewer->m_grid_piece_rens);

  RecursiveTreeItemSelectionModel * sel_model =
      new RecursiveTreeItemSelectionModel ( model, datapiece_treeView );

  datapiece_treeView->setModel ( model );

  datapiece_treeView->setSelectionModel ( sel_model );

}


GridTreeModel::GridTreeModel ( std::vector<grid_piece_rendata *> * dpList, QObject *parent )
  : QAbstractItemModel ( parent )
{
  setupModelData ( dpList );
}

GridTreeModel::~GridTreeModel()
{
  delete m_tree;
  m_tree = NULL;
}


int GridTreeModel::columnCount ( const QModelIndex &/*parent*/ ) const
{
  return 1;
}

QVariant GridTreeModel::data ( const QModelIndex &index, int role ) const
{
  if ( !index.isValid() )
    return QVariant();

  if ( role != Qt::DisplayRole )
    return QVariant();

  tree_item *item = static_cast<tree_item*> ( index.internalPointer() );

  return QString(item->node->dp->label().c_str());
}

Qt::ItemFlags GridTreeModel::flags ( const QModelIndex &index ) const
{
  if ( !index.isValid() )
    return 0;

  return Qt::ItemIsEnabled | Qt::ItemIsSelectable;
}

QVariant GridTreeModel::headerData ( int /*section*/, Qt::Orientation orientation,
                                 int role ) const
{
  if ( orientation == Qt::Horizontal && role == Qt::DisplayRole )
    return "Data Pieces";

  return QVariant();
}

QModelIndex GridTreeModel::index ( int row, int column, const QModelIndex &parent ) const
{
  if ( !hasIndex ( row, column, parent ) )
    return QModelIndex();

  tree_item *parentItem;

  if ( !parent.isValid() )
    parentItem = m_tree;
  else
    parentItem = static_cast<tree_item*> ( parent.internalPointer() );

  if ( row < ( int ) parentItem->children.size() )
    return createIndex ( row, column, parentItem->children[row] );
  else
    return QModelIndex();

}

QModelIndex GridTreeModel::parent ( const QModelIndex &index ) const
{
  if ( !index.isValid() )
    return QModelIndex();

  tree_item *childItem  = static_cast<tree_item*> ( index.internalPointer() );

  tree_item *parentItem = childItem->parent;

  if ( parentItem == m_tree )
    return QModelIndex();

  return createIndex ( parentItem->row(), 0, parentItem );

}

int GridTreeModel::rowCount ( const QModelIndex &parent ) const
{
  tree_item *parentItem;

  if ( parent.column() > 0 )
    return 0;

  if ( !parent.isValid() )
    parentItem = m_tree;
  else
    parentItem = static_cast<tree_item*> ( parent.internalPointer() );

  return parentItem->children.size();

}

void GridTreeModel::setupModelData( std::vector<grid_piece_rendata *> * dpList)
{
  m_tree = new tree_item();


  for ( std::vector<grid_piece_rendata*>::iterator dp_it =  dpList->begin();
          dp_it != dpList->end(); ++dp_it )
  {
    grid_piece_rendata *dp = *dp_it;

    tree_item *parentItem = m_tree;


    tree_item * dpItem = new tree_item ( dp, parentItem );

    parentItem->children.push_back ( dpItem );

  }
}

GridTreeModel::tree_item::tree_item
    ( grid_piece_rendata * _node , tree_item * par ) :
    node ( _node ),parent ( par )
{
}

GridTreeModel::tree_item::tree_item()
{
  node = NULL;
  parent = NULL;
}

int GridTreeModel::tree_item::row()
{
  return std::find ( parent->children.begin(),parent->children.end(),this )
      - parent->children.begin();
}


void RecursiveTreeItemSelectionModel::select
    ( const QModelIndex &index,
      QItemSelectionModel::SelectionFlags command )
{

  QItemSelectionModel::select ( index,command );

  if ( ! ( command & QItemSelectionModel::Select ) )
  {
    return;
  }

  if ( index.isValid() && !m_pTreeView->isExpanded ( index ) )
  {
    QModelIndexList indexes_to_select;

    collect_all_children ( index,indexes_to_select );

    for ( QModelIndexList::iterator ind_it = indexes_to_select.begin();
          ind_it != indexes_to_select.end();++ind_it )
    {
      QItemSelectionModel::select ( *ind_it,QItemSelectionModel::Select );
    }
  }
}

void RecursiveTreeItemSelectionModel::select
    ( const QItemSelection &selection,
      QItemSelectionModel::SelectionFlags command )
{
  QItemSelectionModel::select ( selection,command );
}

void RecursiveTreeItemSelectionModel::collect_all_children
    ( const QModelIndex &index,QModelIndexList &retlist )
{
  QModelIndexList ind_stack;

  ind_stack.push_back ( index );

  while ( ind_stack.size() !=0 )
  {
    QModelIndex top_ind = ind_stack.front();

    ind_stack.pop_front();

    retlist.push_back ( top_ind );

    uint child_ind = 0;

    QModelIndex child = top_ind.child ( child_ind++, top_ind.column() );

    while ( child.isValid() )
    {
      ind_stack.push_back ( child );

      child = top_ind.child ( child_ind++, top_ind.column() );
    }
  }
}
