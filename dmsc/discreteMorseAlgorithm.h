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
    void union_complexes_down
    ( MSComplex<id_t>  *msc1,
      MSComplex<id_t>  *msc2 ,
      const std::vector<std::pair<id_t,id_t> > &cancellable_pairs,
      const std::map<id_t,id_t>  &cancel_pair_map)
{
  typedef typename MSComplex<id_t>::critical_point critical_point_t;
  typedef typename MSComplex<id_t>::critical_point::connection_t connection_t;
  typedef typename connection_t::iterator connection_iter_t;
  typedef typename MSComplex<id_t>::id_cp_map_t id_cp_map_t;
  typedef typename id_cp_map_t::iterator id_cp_map_iter_t;
  typedef typename MSComplex<id_t>::cp_ptr_list_t cp_ptr_list_t;
  typedef typename std::pair<id_t,id_t> cancel_pair_t;
  typedef typename std::vector<cancel_pair_t> cancel_pair_list_t;
  typedef typename cancel_pair_list_t::const_reverse_iterator cancel_pair_list_riter_t;
  typedef typename std::map<id_t,id_t> cancel_pair_map_t;
  typedef typename cancel_pair_map_t::iterator cancel_pair_map_iter_t;

  // carry out the uncancellation operation on msc1
  // in reverse order of the cancellation pairs

  for(cancel_pair_list_riter_t canpr_it = cancellable_pairs.rbegin();
  canpr_it != cancellable_pairs.rend() ; ++canpr_it)
  {

    id_t cp1_id = canpr_it->first;
    id_t cp2_id = canpr_it->second;

    if(msc1->m_id_cp_map.find(cp1_id) ==
       msc1->m_id_cp_map.end())
      throw std::logic_error("could not find cancelable cp1 in msc1");

    if(msc1->m_id_cp_map.find(cp2_id) ==
       msc1->m_id_cp_map.end())
      throw std::logic_error("could not find cancelable cp1 in msc2");

    uint cp1_idx = msc1->m_id_cp_map[cp1_id];
    uint cp2_idx = msc1->m_id_cp_map[cp2_id];

    uncancel_pairs(msc1,cp1_idx,cp2_idx,cancel_pair_map);
  }

  // now identify the uncancelled cps with msc2 and copy over the correct
  // connection info.. clear out the old that was there in the first place.
  // create new cp's for any new incidences.
  for(cancel_pair_list_riter_t canpr_it = cancellable_pairs.rbegin();
      canpr_it != cancellable_pairs.rend() ; ++canpr_it)
  {
    id_t cp1_id = canpr_it->first;
    id_t cp2_id = canpr_it->second;

    // the verification of availability of the id's in msc1
    // was already done in the last loop

    uint msc1_cp1_idx = msc1->m_id_cp_map[cp1_id];
    uint msc1_cp2_idx = msc1->m_id_cp_map[cp2_id];

    if(msc2->m_id_cp_map.find(cp1_id) ==
       msc2->m_id_cp_map.end())
      throw std::logic_error("could not find cancelable cp1 in msc1");

    if(msc2->m_id_cp_map.find(cp2_id) ==
       msc2->m_id_cp_map.end())
      throw std::logic_error("could not find cancelable cp2 in msc2");

    uint msc2_cp1_idx = msc2->m_id_cp_map[cp1_id];
    uint msc2_cp2_idx = msc2->m_id_cp_map[cp2_id];

    connection_t * msc1_acdc_conns[] =
    {&msc1->m_cps[msc1_cp1_idx]->asc,&msc1->m_cps[msc1_cp1_idx]->des,
     &msc1->m_cps[msc1_cp2_idx]->asc,&msc1->m_cps[msc1_cp2_idx]->des};

    connection_t * msc2_acdc_conns[] =
    {&msc2->m_cps[msc2_cp1_idx]->asc,&msc2->m_cps[msc2_cp1_idx]->des,
     &msc2->m_cps[msc2_cp2_idx]->asc,&msc2->m_cps[msc2_cp2_idx]->des};

    for(uint i= 0 ; i < 4; ++i)
    {
      msc2_acdc_conns[i]->clear();

      for(connection_iter_t conn_cp_it = msc1_acdc_conns[i]->begin();
          conn_cp_it != msc1_acdc_conns[i]->end(); ++conn_cp_it)
      {
        id_t conn_cp_id = msc1->m_cps[*conn_cp_it]->cellid;

        if(msc1->m_cps[*conn_cp_it]->isBoundryCancelable == true)
          throw std::logic_error("expected only non bc cps");

        critical_point_t* src_cp = msc1->m_cps[*conn_cp_it];

        if(msc2->m_id_cp_map.find(conn_cp_id) ==
           msc2->m_id_cp_map.end())
        {
          // allocate a new cp
          critical_point_t* dest_cp = new critical_point_t;

          // copy over the trivial data
          dest_cp->isCancelled               = src_cp->isCancelled;
          dest_cp->isBoundryCancelable       = src_cp->isBoundryCancelable;
          dest_cp->isOnStrangulationPath     = src_cp->isOnStrangulationPath;
          dest_cp->cellid                    = src_cp->cellid;

          msc2->m_id_cp_map[dest_cp->cellid] = msc2->m_cps.size();
          msc2->m_cps.push_back(dest_cp);
        }

        uint msc2_conn_cp_idx = msc2->m_id_cp_map[conn_cp_id];

        msc2_acdc_conns[i]->insert(msc2_conn_cp_idx);
      }
    }
  }

  // finally for each non bc cp in msc2 adjust connections such that the
  // the bc cp's are cleared out and any connections with the updated list
  // of cps is reflected

  for(uint i = 0 ; i < msc2->m_cps.size();++i)
  {
    critical_point_t * msc2_cp = msc2->m_cps[i];

    if(msc2_cp->isBoundryCancelable == true)
      continue;

    id_cp_map_iter_t msc1_id_cp_it = msc1->m_id_cp_map.find(msc2_cp->cellid);

    if(msc1_id_cp_it == msc1->m_id_cp_map.end())
      throw std::logic_error("missing non bc cp in msc1 present in msc2");

    critical_point_t *msc1_cp = msc1->m_cps[msc1_id_cp_it->second];

    connection_t * msc1_cp_acdc_conns[] = {&msc1_cp->asc,&msc1_cp->des};
    connection_t * msc2_cp_acdc_conns[] = {&msc2_cp->asc,&msc2_cp->des};

    for(uint j = 0 ; j<2; ++j)
    {
      msc2_cp_acdc_conns[j]->clear();

      for(connection_iter_t conn_it = msc1_cp_acdc_conns[j]->begin();
          conn_it != msc1_cp_acdc_conns[j]->end(); ++conn_it)
      {
        id_t conn_cp_id = msc1->m_cps[*conn_it]->cellid;

        id_cp_map_iter_t msc2_conn_cp_it = msc2->m_id_cp_map.find(conn_cp_id);

        if(msc2_conn_cp_it == msc2->m_id_cp_map.end())
          continue;

        msc2_cp_acdc_conns[j]->insert(msc2_conn_cp_it->second);

      }
    }
  }
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

