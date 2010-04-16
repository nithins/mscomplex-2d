/***************************************************************************
 *   Copyright (C) 2009 by Nithin Shivashankar,   *
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

#ifndef DISCRTEMORSEDS_INCLUDED__
#define DISCRTEMORSEDS_INCLUDED__

#include <vector>
#include <map>
#include <set>
#include <iostream>

#include <glutils.h>
#include <quad_edge.h>

enum eGradDirection {GRADIENT_DIR_DOWNWARD,GRADIENT_DIR_UPWARD,GRADIENT_DIR_COUNT};

template <typename id_type>
    class IDiscreteDataset
{
public:
  virtual id_type getCellPairId ( id_type cellid ) const =0;

  virtual bool    ptLt ( id_type cell1,id_type cell2 ) const =0;

  virtual bool    compareCells( id_type ,id_type ) const = 0;

  virtual uint    getCellPoints ( id_type cellId,id_type  *points ) const =0;

  virtual uint    getCellFacets ( id_type cellId,id_type *facets ) const =0;

  virtual uint    getCellCofacets ( id_type cellId,id_type *cofacets ) const =0;

  virtual bool    isPairOrientationCorrect ( id_type cellId, id_type pairId ) const =0;

  virtual bool    isCellMarked ( id_type cellId ) const =0;

  virtual bool    isCellCritical ( id_type cellId ) const =0;

  virtual void    pairCells ( id_type cellId1,id_type cellId2 ) =0;

  virtual void    markCellCritical ( id_type cellId ) =0;

  virtual uint    getCellDim ( id_type cellId ) const = 0;

  virtual uint    getMaxCellDim() const = 0;

  virtual bool    isTrueBoundryCell ( id_type cellId ) const =0;

  virtual bool    isFakeBoundryCell ( id_type cellId ) const =0;

  virtual bool    isCellExterior ( id_type cellId ) const =0;

  // some functions to support debugging

  virtual std::string  getCellFunctionDescription ( id_type pt ) const = 0;

  virtual std::string getCellDescription ( id_type cellid ) const =0;

};

template <typename id_t>
    class MSComplex
{

public:

  struct critical_point;

  typedef std::map<id_t,uint>           id_cp_map_t;
  typedef std::vector<critical_point *> cp_ptr_list_t;


  struct critical_point
  {
    typedef std::multiset<uint> connection_t;
    typedef std::vector<id_t>   disc_t;

    id_t cellid;

    u_int pair_idx;

    bool isCancelled;
    bool isOnStrangulationPath;
    bool isBoundryCancelable;

    critical_point()
    {
      isCancelled           = false;
      isOnStrangulationPath = false;
      isBoundryCancelable   = false;
      pair_idx              = (u_int) -1;
    }

    ~critical_point()
    {
      asc.clear();
      des.clear();

      asc_disc.clear();
      des_disc.clear();
    }


    disc_t asc_disc;
    disc_t des_disc;

    connection_t asc;
    connection_t des;
  };

  cp_ptr_list_t m_cps;
  id_cp_map_t   m_id_cp_map;

  MSComplex()
  {
  }

  ~MSComplex()
  {

    std::for_each(m_cps.begin(),m_cps.end(),delete_ftor<critical_point>);

    m_cps.clear();
    m_id_cp_map.clear();
  }
};

#endif
