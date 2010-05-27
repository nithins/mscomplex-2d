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
#ifndef __GRID_H__
#define __GRID_H__

#include <vector>
#include <algorithm>
#include <iostream>
#include <set>
#include <map>

#include <boost/bind.hpp>

#include <cpputils.h>
#include <logutil.h>

namespace grid
{
  template <typename coord_type>
      class rectangle_complex
  {

  public:

    typedef unsigned int uint;

    struct size_def:public two_tuple_t<coord_type>
    {
      inline coord_type x() const {return (*this)[0];}
      inline coord_type y() const {return (*this)[1];}

      size_def(const coord_type& x,const coord_type& y)
        :two_tuple_t<coord_type>(x,y){}

      static size_def sub ( const size_def & s1,const size_def & s2 )
      {
        return size_def ( s1.x-s2.x,s1.y-s2.y() );
      }
    };

    struct point_def:public two_tuple_t<coord_type>
    {
      inline const coord_type& x() const {return (*this)[0];}
      inline const coord_type& y() const {return (*this)[1];}

      inline coord_type& x(){return (*this)[0];}
      inline coord_type& y(){return (*this)[1];}


      point_def ( const coord_type &x,const coord_type &y ) :
          two_tuple_t<coord_type>(x,y){}

      point_def () :two_tuple_t<coord_type>(-1,-1){}

      static point_def add ( const point_def & p,const size_def & s )
      {
        return point_def ( p.x()+s.x(),p.y()+s.y() );
      }

      static point_def sub ( const point_def & p,const size_def & s )
      {
        return point_def ( p.x()-s.x(),p.y()-s.y() );
      }

      static size_def sub ( const point_def & p1,const point_def & p2 )
      {
        return size_def ( p1.x()-p2.x(),p1.y()-p2.y() );
      }

      point_def add ( const size_def & s ) const
      {
        return add ( *this,s );
      }

      point_def sub ( const size_def & s )
      {
        return sub ( *this,s );
      }

      size_def sub ( const point_def & p ) const
      {
        return sub ( *this,p );
      }

      inline point_def operator/(coord_type s) const
      {
	return point_def((*this)[0]/s,(*this)[1]/s);
      }

      std::string to_string()
      {
        std::stringstream ss;

        ((std::ostream&)ss)<<(*this);

        return ss.str();
      }

    };

    struct rectangle_def
    {
      typedef rectangle_def _SELF_NAME;

      point_def bl;
      point_def tr;


      size_def size() const
      {
        return tr.sub(bl);
      }


      rectangle_def
          (
              const coord_type & start_x,
              const coord_type & start_y,
              const coord_type & end_x,
              const coord_type & end_y,
              const coord_type & bias = 1
                                        )

      {

        bl.x() = std::min ( start_x*bias,end_x*bias );
        bl.y() = std::min ( start_y*bias,end_y*bias );

        tr.x() = std::max ( start_x*bias,end_x*bias );
        tr.y() = std::max ( start_y*bias,end_y*bias );

      }

      rectangle_def
          (
              const point_def &p1,
              const point_def &p2
              )
      {
        bl.x() = std::min ( p1.x(),p2.x() );
        bl.y() = std::min ( p1.y(),p2.y() );

        tr.x() = std::max ( p1.x(),p2.x() );
        tr.y() = std::max ( p1.y(),p2.y() );

      }

      rectangle_def(){}


      bool isInInterior ( const point_def & p ) const
      {
        return
            (
                ( p.x() > bl.x() ) &&
                ( p.y() > bl.y() ) &&
                ( p.x() < tr.x() ) &&
                ( p.y() < tr.y() )
                );

      }

      bool contains ( const point_def & p ) const
      {
        return
            (
                ( p.x() >= bl.x() ) &&
                ( p.y() >= bl.y() ) &&
                ( p.x() <= tr.x() ) &&
                ( p.y() <= tr.y() )
                );
      }

      bool isOnBoundry ( const point_def & p ) const
      {
        return ( contains ( p ) && !isInInterior ( p ) );
      }

      bool contains ( const _SELF_NAME &rec ) const
      {
        return
            (
                ( bl.x() < rec.bl.x() ) &&
                ( bl.y() < rec.bl.y() ) &&
                ( tr.x() > rec.tr.x() ) &&
                ( tr.y() > rec.tr.y() )
                );
      }

      bool intersects ( const _SELF_NAME &rec ) const
      {
        return
            ! (
                ( ( tr.x() ) < ( rec.bl.x() ) ) ||
                ( ( rec.tr.x() ) < ( bl.x() ) ) ||
                ( ( tr.y() ) < ( rec.bl.y() ) ) ||
                ( ( rec.tr.y() ) < ( bl.y() ) )
                );
      }

      bool intersection(const rectangle_def & rect,rectangle_def &i) const
      {
        if(!intersects(rect))
          return false;

        coord_type l = std::max(rect.left(),left());
        coord_type r = std::min(rect.right(),right());
        coord_type b = std::max(rect.bottom(),bottom());
        coord_type t = std::min(rect.top(),top());

        i = rectangle_def(l,b,r,t);
        return true;
      }

      void mul ( const coord_type & b )
      {
        bl.x() *= b;
        bl.y() *= b;
        tr.x() *= b;
        tr.y() *= b;
      }

      void div ( const coord_type & b )
      {
        bl.x() /= b;
        bl.y() /= b;
        tr.x() /= b;
        tr.y() /= b;
      }

      void shrink ( const size_def &s )
      {
        bl = point_def::add ( bl,s );
        tr = point_def::sub ( tr,s );
      }

      void grow ( const size_def &s )
      {
        bl = point_def::sub ( bl,s );
        tr = point_def::add ( tr,s );
      }

      point_def mid() const
      {
        return point_def ( ( bl.x()+tr.x() ) /2, ( bl.y()+tr.y() ) /2 );
      }

      coord_type left() const
      {
        return bl.x();
      }

      coord_type right() const
      {
        return tr.x();
      }

      coord_type top() const
      {
        return tr.y();
      }

      coord_type bottom() const
      {
        return bl.y();
      }

      point_def bottom_left() const
      {
        return point_def ( left(), bottom() );
      }

      point_def top_left() const
      {
        return point_def ( left(),top() );
      }

      point_def top_right() const
      {
        return point_def ( right(),top() );
      }
      point_def bottom_right() const
      {
        return point_def ( right(),bottom() );
      }

      uint getNumVerts()
      {
        return ( tr.x-bl.x+1 ) * ( tr.y-bl.y+1 );
      }

      uint getNumQuads()
      {
        return ( tr.x-bl.x() ) * ( tr.y-bl.y() );
      }

      friend std::ostream& operator<< ( std::ostream& o, const rectangle_def& r )
      {
        return o<<"[ bl="<<r.bl<<" tr= "<<r.tr<<"]";
      }
    };


    typedef rectangle_complex _SELF_NAME;

    rectangle_def m_rec_def;

    std::vector<rectangle_def> m_outside_rects;

    std::vector<rectangle_def> m_inside_rects;

    std::vector<rectangle_def> m_inside_lines;

  public:
    rectangle_complex
        (
            coord_type  start_x,
            coord_type  start_y,
            coord_type  end_x,
            coord_type  end_y
            ) :
        m_rec_def ( start_x,start_y,end_x,end_y,2 )
    {
      init();
    }

    rectangle_complex ( rectangle_def rec )
      :m_rec_def ( rec )
    {
      m_rec_def.mul ( 2 );
      init();
    }

    void init()
    {
      rectangle_def r = m_rec_def;
      r.shrink ( size_def ( 1,1 ) );
      m_inside_rects.push_back ( r );

    }
  };

  const uint gc_grid_dim = 2;

  typedef int16_t                          cell_coord_t;
  typedef float                            cell_fn_t;
  typedef rectangle_complex<cell_coord_t>  rect_cmplx_t;
  typedef rect_cmplx_t::rectangle_def      rect_t;
  typedef rect_cmplx_t::point_def          cellid_t;
  typedef rect_cmplx_t::point_def          rect_point_t;
  typedef rect_cmplx_t::size_def           rect_size_t;
  typedef std::vector<cellid_t>            cellid_list_t;
  typedef unsigned int                     critpt_idx_t;
  typedef std::vector<critpt_idx_t>        critpt_idx_list_t;


  enum eGradientDirection
  {
    GRADDIR_DESCENDING,
    GRADDIR_ASCENDING,
    GRADDIR_COUNT,
  };

}


#endif
