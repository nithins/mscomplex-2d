/***************************************************************************
 *   Copyright (C) 2009 by nithin,,,   *
 *   nithin@gauss   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/

#ifndef GRID_DATAMANAGER_H_INCLUDED_
#define GRID_DATAMANAGER_H_INCLUDED_

#include <fstream>
#include <vector>

#include <grid_dataset.h>

struct GridDataPiece
{
  typedef GridDataset::cell_fn_t cell_fn_t;
  typedef GridDataset::mscomplex_t mscomplex_t;
  typedef GridDataset::cell_coord_t cell_coord_t;
  typedef GridDataset::rect_t rect_t;
  typedef GridDataset::rect_size_t rect_size_t;
  typedef GridDataset::cellid_t cellid_t;
  typedef GridDataset::critpt_conn_t conn_t;

  GridDataset *dataset;
  mscomplex_t *msgraph;

  uint level;

  uint m_pieceno;

  GridDataPiece (uint l);

  std::string label();
};

namespace boost
{
  class thread;
}

class GridDataManager
{

  typedef GridDataset::rect_t rect_t;
  typedef GridDataset::cell_coord_t cell_coord_t;
  typedef GridDataset::cellid_t cellid_t;
  typedef GridDataset::rect_point_t rect_point_t;
  typedef GridDataset::rect_size_t rect_size_t;
  typedef std::vector<GridDataPiece *> pieces_list_t;

public:

  pieces_list_t                m_pieces;

  std::string                  m_filename;
  u_int                        m_size_x;
  u_int                        m_size_y;
  u_int                        m_num_levels;
  double                       m_simp_tresh;
  bool                         m_single_threaded_mode;
  bool                         m_use_ocl;
  bool                         m_compute_out_of_core;

  boost::thread **             m_threads;

public:

  uint num_parallel;


  GridDataManager
      ( std::string filename,
        u_int        size_x,
        u_int        size_y,
        u_int        num_levels,
        bool         threaded_mode,
        bool         use_ocl,
        double       simp_tresh,
        bool         compute_out_of_core,
        uint         np);

  virtual ~GridDataManager ();

  void createPieces_quadtree(rect_t r,rect_t e,u_int level );

  void createDataPieces();

  void readDataAndInit(std::ifstream &data_stream,GridDataset::cell_fn_t *,uint start_offset);

  uint getMaxDataBufItems();

  void waitForThreadsInRange(uint,uint);

  void computeMsGraph ( GridDataPiece  * );

  void computeMsGraphInRange(uint ,uint );

  void finalMergeDownPiecesInRange(uint start,uint end);

  void collectManifold( GridDataPiece  * );

  void collectManifoldsInRange(uint start,uint end);

  void writeManifoldsInRange(uint start,uint end);

  void computeSubdomainMsgraphs ();

  void mergePiecesUp( );

  void mergePiecesDown( );

  void collectSubdomainManifolds( );

  void logAllConnections(const std::string &prefix);

  void logAllCancelPairs(const std::string &prefix);

};

#endif
