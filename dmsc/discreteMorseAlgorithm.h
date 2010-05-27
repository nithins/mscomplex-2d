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

#ifndef __DISCRETEMORSEALGORITHM_H_INCLUDED
#define __DISCRETEMORSEALGORITHM_H_INCLUDED

#include <queue>
#include <vector>
#include <list>
#include <algorithm>
#include <stack>
#include <cstring>

#include <boost/bind.hpp>

#include <logutil.h>
#include <cpputils.h>

#include <discreteMorseDS.h>

#define LOG_LEVEL 1

#define DEBUG_BEGIN if(LOG_LEVEL > 0) {
#define DEBUG_END }
#define DEBUG_STMT(x) DEBUG_BEGIN x; DEBUG_END;
#define DEBUG_LOG(x) DEBUG_STMT(_LOG(x))

#define DEBUG_BEGIN_V if(LOG_LEVEL > 1) {
#define DEBUG_END_V }
#define DEBUG_STMT_V(x) DEBUG_BEGIN_V x; DEBUG_END_V;
#define DEBUG_LOG_V(x) DEBUG_STMT_V(_LOG(x))

#define DEBUG_BEGIN_VV if(LOG_LEVEL > 2) {
#define DEBUG_END_VV }
#define DEBUG_STMT_VV(x) DEBUG_BEGIN_VV x; DEBUG_END_VV;
#define DEBUG_LOG_VV(x) DEBUG_STMT_VV(_LOG(x))

template <typename id_t,typename cpiter_t>
    void addCriticalPointsToMSComplex
    (
        MSComplex<id_t> *msc,
        cpiter_t begin,
        cpiter_t end
        )
{
  typedef typename MSComplex<id_t>::critical_point cp_t;

  msc->m_cps.resize(end-begin);

  uint pos = 0;

  for ( cpiter_t it = begin;it != end ; ++it )
  {
    cp_t * cp                 = new cp_t;
    msc->m_cps[pos]           = cp;
    cp->cellid                = *it;
    msc->m_id_cp_map.insert ( std::make_pair ( *it,pos ) );

    ++pos;
  }
}

template <typename id_t>
    void cancelPairs ( MSComplex<id_t> *msc,uint cp0_ind,uint cp1_ind )
{

  typedef typename std::multiset<uint> cp_adj_t;

  DEBUG_LOG_V ( "=========================" );
  DEBUG_LOG_V ( "cancelling pairs "<<msc->m_cps[cp0_ind]->cellid<<" "<<msc->m_cps[cp1_ind]->cellid );
  DEBUG_LOG_V ( "-------------------------" );

  if ( msc->m_cps[cp0_ind]->asc.find ( cp1_ind ) == msc->m_cps[cp0_ind]->asc.end() )
  {
    std::swap ( cp0_ind,cp1_ind );
  }


  if ( msc->m_cps[cp0_ind]->asc.find ( cp1_ind ) == msc->m_cps[cp0_ind]->asc.end() )
  {
    _ERROR ( "pair =  "<<msc->m_cps[cp0_ind]->cellid<<" "<<msc->m_cps[cp1_ind]->cellid );
    _ERROR ( "critical points are not connected" );
    return;
  }



  //   std::stringstream ss;
  //
  //   ss<<"asc("<<msc->m_cps[cp0_ind]->cellid<<")="<<"{";
  //
  //   for ( typename cp_adj_t::iterator asc0_it = msc->m_cps[cp0_ind]->asc.begin();asc0_it != msc->m_cps[cp0_ind]->asc.end(); ++asc0_it )
  //   {
  //     ss<<msc->m_cps[*asc0_it]->cellid<<",";
  //   }
  //   ss<<"}";
  //
  //   _LOG ( ss.str() );


  for ( typename cp_adj_t::iterator asc0_it = msc->m_cps[cp0_ind]->asc.begin();asc0_it != msc->m_cps[cp0_ind]->asc.end(); ++asc0_it )
  {
    for ( typename cp_adj_t::iterator des1_it = msc->m_cps[cp1_ind]->des.begin();des1_it != msc->m_cps[cp1_ind]->des.end(); ++des1_it )
    {
      if ( ( *asc0_it == cp1_ind ) || ( *des1_it == cp0_ind ) )
        continue;

      msc->m_cps[ *asc0_it ]->des.insert ( *des1_it );
      msc->m_cps[ *des1_it ]->asc.insert ( *asc0_it );

      DEBUG_LOG_V ( "Connected "<<msc->m_cps[*asc0_it]->cellid<<" to "<<msc->m_cps[*des1_it]->cellid );

      if ( msc->m_cps[ *asc0_it ]->des.count ( *des1_it ) >= 2 )
      {
        msc->m_cps[ *asc0_it ]->isOnStrangulationPath = true;
        msc->m_cps[ *des1_it ]->isOnStrangulationPath = true;

        DEBUG_LOG_V ( "Strangulation by "<<msc->m_cps[*asc0_it]->cellid<<" "<<msc->m_cps[*des1_it]->cellid );

      }
    }
  }

  for ( typename cp_adj_t::iterator asc0_it = msc->m_cps[cp0_ind]->asc.begin();asc0_it != msc->m_cps[cp0_ind]->asc.end(); ++asc0_it )
  {
    if ( *asc0_it == cp1_ind )  continue;

//    DEBUG_LOG_V ( "merging des_manifold of "<< cp1_ind<<" with " << *asc0_it );

//    msc->m_cps[ *asc0_it ]->des_disc.insert ( msc->m_cps[ cp1_ind ]->des_disc.begin(),msc->m_cps[ cp1_ind ]->des_disc.end() );

    DEBUG_LOG_V ( "Removing   "<<cp0_ind<<" from "<< ( *asc0_it ) );

    msc->m_cps[  *asc0_it ]->des.erase ( cp0_ind );
  }

  for ( typename cp_adj_t::iterator des0_it = msc->m_cps[cp0_ind]->des.begin();des0_it != msc->m_cps[cp0_ind]->des.end(); ++des0_it )
  {
    msc->m_cps[ ( *des0_it ) ]->asc.erase ( cp0_ind );

    DEBUG_LOG_V ( "Removing   "<<cp0_ind<<" from "<< ( *des0_it ) );
  }

  for ( typename cp_adj_t::iterator asc1_it = msc->m_cps[cp1_ind]->asc.begin();asc1_it != msc->m_cps[cp1_ind]->asc.end(); ++asc1_it )
  {
    msc->m_cps[ ( *asc1_it ) ]->des.erase ( cp1_ind );

    DEBUG_LOG_V ( "Removing   "<<cp1_ind<<" from "<< ( *asc1_it ) );
  }

  for ( typename cp_adj_t::iterator des1_it = msc->m_cps[cp1_ind]->des.begin();des1_it != msc->m_cps[cp1_ind]->des.end(); ++des1_it )
  {
    if ( ( *des1_it ) == cp0_ind )  continue;

//    DEBUG_LOG_V ( "merging asc_manifold of "<< cp0_ind<<" with " << *des1_it );

//    msc->m_cps[ *des1_it ]->asc_disc.insert ( msc->m_cps[ cp0_ind ]->asc_disc.begin(),msc->m_cps[ cp0_ind ]->asc_disc.end() );

    DEBUG_LOG_V ( "Removing   "<<cp1_ind<<" from "<< ( *des1_it ) );

    msc->m_cps[ ( *des1_it ) ]->asc.erase ( cp1_ind );


  }

  msc->m_cps[cp0_ind]->isCancelled = true;
  //  msc->m_cps[cp0_ind]->asc.clear();
  msc->m_cps[cp0_ind]->des.clear();
  //  msc->m_cps[cp0_ind]->asc_disc.clear();
  //  msc->m_cps[cp0_ind]->des_disc.clear();

  msc->m_cps[cp1_ind]->isCancelled = true;
  msc->m_cps[cp1_ind]->asc.clear();
  //  msc->m_cps[cp1_ind]->des.clear();
  //  msc->m_cps[cp1_ind]->asc_disc.clear();
  //  msc->m_cps[cp1_ind]->des_disc.clear();

  DEBUG_LOG_V ( "=========================" );
}

template <typename id_t>
    void cancelPairs ( MSComplex<id_t> *msc,uint cp0_ind,uint cp1_ind,std::vector<uint> & new_edges )
{

  typedef typename std::multiset<uint> cp_adj_t;

  DEBUG_LOG_V ( "=========================" );
  DEBUG_LOG_V ( "cancelling pairs "<<msc->m_cps[cp0_ind]->cellid<<" "<<msc->m_cps[cp1_ind]->cellid );
  DEBUG_LOG_V ( "-------------------------" );

  if ( msc->m_cps[cp0_ind]->asc.find ( cp1_ind ) == msc->m_cps[cp0_ind]->asc.end() )
  {
    std::swap ( cp0_ind,cp1_ind );
  }


  if ( msc->m_cps[cp0_ind]->asc.find ( cp1_ind ) == msc->m_cps[cp0_ind]->asc.end() )
  {
    _ERROR ( "pair =  "<<msc->m_cps[cp0_ind]->cellid<<" "<<msc->m_cps[cp1_ind]->cellid );
    _ERROR ( "critical points are not connected" );
    return;
  }



  //   std::stringstream ss;
  //
  //   ss<<"asc("<<msc->m_cps[cp0_ind]->cellid<<")="<<"{";
  //
  //   for ( typename cp_adj_t::iterator asc0_it = msc->m_cps[cp0_ind]->asc.begin();asc0_it != msc->m_cps[cp0_ind]->asc.end(); ++asc0_it )
  //   {
  //     ss<<msc->m_cps[*asc0_it]->cellid<<",";
  //   }
  //   ss<<"}";
  //
  //   _LOG ( ss.str() );


  for ( typename cp_adj_t::iterator asc0_it = msc->m_cps[cp0_ind]->asc.begin();asc0_it != msc->m_cps[cp0_ind]->asc.end(); ++asc0_it )
  {
    for ( typename cp_adj_t::iterator des1_it = msc->m_cps[cp1_ind]->des.begin();des1_it != msc->m_cps[cp1_ind]->des.end(); ++des1_it )
    {
      if ( ( *asc0_it == cp1_ind ) || ( *des1_it == cp0_ind ) )
        continue;

      msc->m_cps[ *asc0_it ]->des.insert ( *des1_it );
      msc->m_cps[ *des1_it ]->asc.insert ( *asc0_it );

      DEBUG_LOG_V ( "Connected "<<msc->m_cps[*asc0_it]->cellid<<" to "<<msc->m_cps[*des1_it]->cellid );

      new_edges.push_back(*des1_it );
      new_edges.push_back(*asc0_it );

      if ( msc->m_cps[ *asc0_it ]->des.count ( *des1_it ) >= 2 )
      {
        msc->m_cps[ *asc0_it ]->isOnStrangulationPath = true;
        msc->m_cps[ *des1_it ]->isOnStrangulationPath = true;

        DEBUG_LOG_V ( "Strangulation by "<<msc->m_cps[*asc0_it]->cellid<<" "<<msc->m_cps[*des1_it]->cellid );

      }
    }
  }

  for ( typename cp_adj_t::iterator asc0_it = msc->m_cps[cp0_ind]->asc.begin();asc0_it != msc->m_cps[cp0_ind]->asc.end(); ++asc0_it )
  {
    if ( *asc0_it == cp1_ind )  continue;

//    DEBUG_LOG_V ( "merging des_manifold of "<< cp1_ind<<" with " << *asc0_it );

//    msc->m_cps[ *asc0_it ]->des_disc.insert ( msc->m_cps[ cp1_ind ]->des_disc.begin(),msc->m_cps[ cp1_ind ]->des_disc.end() );

    DEBUG_LOG_V ( "Removing   "<<cp0_ind<<" from "<< ( *asc0_it ) );

    msc->m_cps[  *asc0_it ]->des.erase ( cp0_ind );
  }

  for ( typename cp_adj_t::iterator des0_it = msc->m_cps[cp0_ind]->des.begin();des0_it != msc->m_cps[cp0_ind]->des.end(); ++des0_it )
  {
    msc->m_cps[ ( *des0_it ) ]->asc.erase ( cp0_ind );

    DEBUG_LOG_V ( "Removing   "<<cp0_ind<<" from "<< ( *des0_it ) );
  }

  for ( typename cp_adj_t::iterator asc1_it = msc->m_cps[cp1_ind]->asc.begin();asc1_it != msc->m_cps[cp1_ind]->asc.end(); ++asc1_it )
  {
    msc->m_cps[ ( *asc1_it ) ]->des.erase ( cp1_ind );

    DEBUG_LOG_V ( "Removing   "<<cp1_ind<<" from "<< ( *asc1_it ) );
  }

  for ( typename cp_adj_t::iterator des1_it = msc->m_cps[cp1_ind]->des.begin();des1_it != msc->m_cps[cp1_ind]->des.end(); ++des1_it )
  {
    if ( ( *des1_it ) == cp0_ind )  continue;

//    DEBUG_LOG_V ( "merging asc_manifold of "<< cp0_ind<<" with " << *des1_it );

//    msc->m_cps[ *des1_it ]->asc_disc.insert ( msc->m_cps[ cp0_ind ]->asc_disc.begin(),msc->m_cps[ cp0_ind ]->asc_disc.end() );

    DEBUG_LOG_V ( "Removing   "<<cp1_ind<<" from "<< ( *des1_it ) );

    msc->m_cps[ ( *des1_it ) ]->asc.erase ( cp1_ind );


  }

  msc->m_cps[cp0_ind]->isCancelled = true;
  //  msc->m_cps[cp0_ind]->asc.clear();
  msc->m_cps[cp0_ind]->des.clear();
  //  msc->m_cps[cp0_ind]->asc_disc.clear();
  //  msc->m_cps[cp0_ind]->des_disc.clear();

  msc->m_cps[cp1_ind]->isCancelled = true;
  msc->m_cps[cp1_ind]->asc.clear();
  //  msc->m_cps[cp1_ind]->des.clear();
  //  msc->m_cps[cp1_ind]->asc_disc.clear();
  //  msc->m_cps[cp1_ind]->des_disc.clear();

  DEBUG_LOG_V ( "=========================" );
}

template <typename id_t>
    void uncancel_pairs
    ( MSComplex<id_t>  *msc,
      uint cp1_idx,
      uint cp2_idx)
{

  typedef typename MSComplex<id_t>::critical_point critical_point_t;
  typedef typename MSComplex<id_t>::critical_point::connection_t connection_t;
  typedef typename connection_t::iterator connection_iter_t;
  typedef typename MSComplex<id_t>::id_cp_map_t id_cp_map_t;
  typedef typename MSComplex<id_t>::cp_ptr_list_t cp_ptr_list_t;
  typedef typename std::pair<id_t,id_t> cancel_pair_t;
  typedef typename std::vector<cancel_pair_t> cancel_pair_list_t;

  if( msc->m_cps[cp1_idx]->asc.find ( cp2_idx ) ==
      msc->m_cps[cp1_idx]->asc.end() )
  {
    std::swap ( cp1_idx,cp2_idx );
  }

  connection_t * acdc_conns[] =
  {&msc->m_cps[cp1_idx]->asc,
   &msc->m_cps[cp2_idx]->des};

  uint * cp_idxs[]= {&cp1_idx,&cp2_idx};


  for(uint i = 0 ;i < 2 ; ++i)
  {

    if( acdc_conns[i]->find ( *cp_idxs[(i+1)%2] ) ==
        acdc_conns[i]->end())
    {
      throw std::logic_error("cancellable pair is not connected");
    }

    connection_t new_acdc;

    for(connection_iter_t acdc_conn_it = acdc_conns[i]->begin();
        acdc_conn_it != acdc_conns[i]->end();++acdc_conn_it)
    {
      if(*acdc_conn_it == *cp_idxs[(i+1)%2])
        continue;

      if(msc->m_cps[*acdc_conn_it]->isBoundryCancelable == false)
      {
        new_acdc.insert(*acdc_conn_it);
        continue;
      }

      if(msc->m_cps[*acdc_conn_it]->isCancelled == true)
        throw std::logic_error("*acdc_conn_it should not have been cancelled yet");


      uint conn_cp_pair_idx = msc->m_cps[*acdc_conn_it]->pair_idx;

      if(msc->m_cps[conn_cp_pair_idx]->pair_idx != *acdc_conn_it)
        throw std::logic_error("*acdc_conn_it and its pair dont agree on pairing");

      connection_t * conn_cp_acdc_conns[] =
      {&msc->m_cps[conn_cp_pair_idx]->asc,
       &msc->m_cps[conn_cp_pair_idx]->des};

      for(connection_iter_t conn_cp_acdc_conn_it =
          conn_cp_acdc_conns[i]->begin();
      conn_cp_acdc_conn_it != conn_cp_acdc_conns[i]->end();
      ++conn_cp_acdc_conn_it)
      {
        if(msc->m_cps[*conn_cp_acdc_conn_it]->isBoundryCancelable == true)
          throw std::logic_error("the connected cancelable cp shold not "\
                                 " contain a cancelable cp in its connections");

        new_acdc.insert(*conn_cp_acdc_conn_it);
      }
    }

    acdc_conns[i]->clear();
    acdc_conns[i]->insert(new_acdc.begin(),new_acdc.end());
  }

  msc->m_cps[cp1_idx]->isCancelled = false;
  msc->m_cps[cp2_idx]->isCancelled = false;
}

template <typename id_t>
void print_connections
(std::ostream & os,
 const MSComplex<id_t> &msc,
 const typename MSComplex<id_t>::critical_point::connection_t &conn
 )
{
  typedef typename MSComplex<id_t>::critical_point::connection_t connection_t;
  typedef typename connection_t::iterator conn_iter_t;

  os<<"{ ";
  for(conn_iter_t it = conn.begin(); it != conn.end(); ++it)
  {
    if(msc.m_cps[*it]->isBoundryCancelable)
      os<<"*";
    if(msc.m_cps[*it]->isOnStrangulationPath)
      os<<"-";
    os<<msc.m_cps[*it]->cellid;
    os<<", ";
  }
  os<<"}";
}

template <typename id_t>
void print_connections
(std::ostream & os,const MSComplex<id_t> &msc)
{
  for(uint i = 0 ; i < msc.m_cps.size();++i)
  {
    os<<"des(";
    if(msc.m_cps[i]->isBoundryCancelable)
      os<<"*";
    if(msc.m_cps[i]->isOnStrangulationPath)
      os<<"-";
    os<<msc.m_cps[i]->cellid<<") = ";
    print_connections(os,msc,msc.m_cps[i]->des);
    os<<std::endl;

    os<<"asc(";
    if(msc.m_cps[i]->isBoundryCancelable)
      os<<"*";
    if(msc.m_cps[i]->isOnStrangulationPath)
      os<<"-";
    os<<msc.m_cps[i]->cellid<<") = ";
    print_connections(os,msc,msc.m_cps[i]->asc);
    os<<std::endl;
    os<<std::endl;
  }
}


#endif

