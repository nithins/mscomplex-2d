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

#include <iostream>
#include <fstream>

#include <boost/thread.hpp>
#include <boost/bind.hpp>
#include <boost/regex.hpp>
#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>

#include <timer.h>

#include <grid_datamanager.h>
#include <discreteMorseAlgorithm.h>


using namespace std;

void GridDataManager::createDataPieces ()
{
  rect_t r(cellid_t(0,0),cellid_t(2*(m_size_x-1),2*(m_size_y-1)));
  rect_t e(cellid_t(0,0),cellid_t(2*(m_size_x-1),2*(m_size_y-1)));

  createPieces_quadtree(r,e,m_num_levels);
  return;
}

const unsigned char split_axes = 1;
const unsigned char scan_axes  = (split_axes+1)%2;

void GridDataManager::createPieces_quadtree(rect_t r,rect_t e,u_int level )
{
  if(level == 0)
  {
    GridDataPiece *dp = new GridDataPiece(m_pieces.size());
    dp->dataset = new GridDataset(r,e);
    dp->msgraph = new GridMSComplex(r,e);
    m_pieces.push_back(dp);

    return;
  }

  u_int dim = split_axes;

  rect_size_t  s = r.size();

  rect_point_t tr1 = r.top_right();
  rect_point_t bl1 = r.bottom_left();

  tr1[dim] -= 2*(s[dim]/4);
  bl1[dim]  = tr1[dim];

  rect_t r1 = rect_t(r.bottom_left(),tr1);
  rect_t r2 = rect_t(bl1,r.top_right());

  rect_point_t tr2 = e.top_right();
  rect_point_t bl2 = e.bottom_left();

  tr2[dim] = tr1[dim] + 2;
  bl2[dim] = bl1[dim] - 2;

  rect_t e1 = rect_t(e.bottom_left(),tr2);
  rect_t e2 = rect_t(bl2,e.top_right());

  createPieces_quadtree(r1,e1,level-1);
  createPieces_quadtree(r2,e2,level-1);
}

void read_msgraph_from_archive(GridDataPiece * dp)
{
  std::string filename("dp_msgraph_");
  filename += dp->label();

  std::ifstream ifs(filename.c_str());

  boost::archive::binary_iarchive ia(ifs);

  dp->msgraph = new GridDataPiece::mscomplex_t;

  ia >> (*dp->msgraph);
}

void write_msgraph_to_archive(GridDataPiece * dp)
{
  std::string filename("dp_msgraph_");
  filename += dp->label();

  std::ofstream ofs(filename.c_str());

  boost::archive::binary_oarchive oa(ofs);

  oa << (*dp->msgraph);

  delete dp->msgraph;
  dp->msgraph = NULL;

}

void GridDataManager::computeMsGraph( GridDataPiece *dp )
{
  if(m_use_ocl != true)
  {
    dp->dataset->work();

    dp->dataset->writeout_connectivity(dp->msgraph);
  }
  else
  {
    dp->dataset->writeout_connectivity_ocl(dp->msgraph);
  }

  dp->dataset->clear_fnref();

  if(m_compute_out_of_core)
  {
    write_msgraph_to_archive(dp);

    dp->dataset->clear();
  }
}

void mergePiecesUp_worker
    ( GridDataPiece  *dp,
      GridDataPiece  *dp1,
      GridDataPiece  *dp2,
      bool is_src_archived,
      bool archive_dest)
{

  if(dp1->level != dp2->level)
      throw std::logic_error("dps must have same level");

  if(is_src_archived == true)
  {
    read_msgraph_from_archive(dp1);
    read_msgraph_from_archive(dp2);
  }

  dp->level         = dp1->level+1;
  dp->msgraph       = GridMSComplex::merge_up(*dp1->msgraph,*dp2->msgraph);

  if(is_src_archived == true)
  {
    delete dp1->msgraph;dp1->msgraph = NULL;
    delete dp2->msgraph;dp2->msgraph = NULL;
  }

  if(archive_dest == true)
  {
    write_msgraph_to_archive(dp);
  }
}


void mergePiecesDown_worker
    ( GridDataPiece  * dp,
      GridDataPiece  * dp1,
      GridDataPiece  * dp2,
      bool is_src_archived,
      bool is_dest_archived,
      bool archive_src,
      bool archive_dest)
{

  if(is_src_archived)
    read_msgraph_from_archive(dp);

  if(is_dest_archived)
  {
    read_msgraph_from_archive(dp1);
    read_msgraph_from_archive(dp2);
  }

  dp->msgraph->merge_down(*dp1->msgraph,*dp2->msgraph);

  if(archive_src)
  {
    write_msgraph_to_archive(dp);
  }

  if(is_src_archived && !archive_src)
  {
    delete dp->msgraph;dp->msgraph = NULL;
  }
  if(archive_dest)
  {
    write_msgraph_to_archive(dp1);
    write_msgraph_to_archive(dp2);
  }
}


void GridDataManager::computeMsGraphInRange(uint start,uint end )
{
  _LOG ( "Begin Gradient/ conn work for pieces from "<<start<<" to "<<end );

  for ( uint i = start ; i < end;i++ )
  {
    GridDataPiece * dp = m_pieces[i];

    dp->dataset->init();

    if(m_use_ocl == true)
      dp->dataset->work_ocl();

    if(m_single_threaded_mode == false)
    {
      _LOG ( "Kicking off thread "<<i-start );


      m_threads[i-start] = new boost::thread
                          ( boost::bind ( &GridDataManager::computeMsGraph,this,dp ) );
    }
    else
    {
      computeMsGraph(dp);
    }

  }
}

void GridDataManager::waitForThreadsInRange(uint start,uint end )
{
  if(m_single_threaded_mode == false)
  {
    _LOG ( "waiting on threads "<<0<<" to "<<end-start );

    for ( uint i = start ; i < end;i++ )
    {
      m_threads[i-start]->join();

      delete m_threads[i-start];

      m_threads[i-start] = NULL;
      _LOG ( "thread "<<i-start<<" joint" );
    }
  }
}

void GridDataManager::mergePiecesUp( )
{
  uint num_leafs = pow(2,m_num_levels);

  for(uint i = 0;i<num_leafs-1;++i)
  {
    m_pieces.push_back(new GridDataPiece(m_pieces.size()));
  }

  uint i_incr = std::min((size_t)num_parallel*2,(m_pieces.size()+1)/2);

  for ( uint i = 0; i<m_pieces.size()-1; i+=i_incr)
  {
    uint threadno = 0;

    if(i_incr+i > m_pieces.size() && i_incr > 1)
      i_incr = i_incr>>1;

    for(int j = i ; j < i+i_incr; j +=2 )
    {

      GridDataPiece *dp1 = m_pieces[j];
      GridDataPiece *dp2 = m_pieces[j+1];
      GridDataPiece *dp  = m_pieces[num_leafs+j/2];

      bool src_archived = m_compute_out_of_core;
      bool archive_dest = m_compute_out_of_core&&(dp1->level != m_num_levels-1) ;


      if(m_single_threaded_mode == false)
      {
        _LOG("Kicking off Merge Up "<<j<<" "<<j+1<<"->"<<num_leafs+j/2);

        m_threads[threadno++] = new boost::thread
                              ( boost::bind ( &mergePiecesUp_worker,dp,dp1,dp2,src_archived,archive_dest ));

      }
      else
      {
        _LOG("Kicking off Merge Up "<<j<<" "<<j+1<<"->"<<num_leafs+j/2);

        mergePiecesUp_worker(dp,dp1,dp2,src_archived,archive_dest);
      }
    }

    waitForThreadsInRange(i,i+i_incr/2);
  }

}

void GridDataManager::mergePiecesDown()
{
  uint num_graphs = (0x01<<(m_num_levels+1))-1;

  uint i_incr = 1;

  for ( int i = 1; i <(0x01<<(m_num_levels-1)); )
  {
    uint threadno = 0 ;

    bool src_archived  = m_compute_out_of_core && (i != 1);
    bool dest_archived = m_compute_out_of_core;

    bool archive_src   = m_compute_out_of_core ;
    bool archive_dest  = m_compute_out_of_core ;

    for ( int j = i; j < i+i_incr; j += 1)
    {

      GridDataPiece *dp  = m_pieces[num_graphs- j];
      GridDataPiece *dp1 = m_pieces[num_graphs- 2*j];
      GridDataPiece *dp2 = m_pieces[num_graphs- 2*j-1];

      if(m_single_threaded_mode == false)
      {
        _LOG("Kicking off Merge Down "<< num_graphs - j<<"->"<<
             num_graphs - 2*j<<" "<<num_graphs - 2*j-1);

        m_threads[threadno++] = new boost::thread
                              ( boost::bind ( &mergePiecesDown_worker,dp,dp1,dp2,src_archived,dest_archived,archive_src,archive_dest ));
      }
      else
      {
        _LOG("Merge Down "<< num_graphs - j<<"->"<<
             num_graphs - 2*j<<" "<<num_graphs - 2*j-1);

        mergePiecesDown_worker(dp,dp1,dp2,src_archived,dest_archived,archive_src,archive_dest);
      }
    }

    waitForThreadsInRange(i,i+i_incr);

    i+=i_incr;

    if(i_incr < num_parallel)
      i_incr <<= 1;
  }

  return;
}

void GridDataManager::finalMergeDownPiecesInRange(uint start,uint end)
{

  uint num_subdomains =  ((0x01)<<m_num_levels);

  bool src_archived  = m_compute_out_of_core && (m_num_levels > 1);
  bool dest_archived = m_compute_out_of_core;

  bool archive_src   = false ;
  bool archive_dest  = false ;

  uint threadno = 0;

  for ( int j = start; j < end; j += 2)
  {

    GridDataPiece *dp  = m_pieces[num_subdomains + j/2];
    GridDataPiece *dp1 = m_pieces[j];
    GridDataPiece *dp2 = m_pieces[j+1];

    if(m_single_threaded_mode == false)
    {
      _LOG("Kicking off Merge Down "<< num_subdomains + j/2<<"->"<<
           j+1<<" "<<j);

      m_threads[threadno++] = new boost::thread
                            ( boost::bind ( &mergePiecesDown_worker,dp,dp1,dp2,src_archived,dest_archived,archive_src,archive_dest ));
    }
    else
    {
      _LOG("Merge Down "<< num_subdomains + j/2<<"->"<<
           j+1<<" "<<j);

      mergePiecesDown_worker(dp,dp1,dp2,src_archived,dest_archived,archive_src,archive_dest);
    }
  }
}

void GridDataManager::collectManifold( GridDataPiece  * dp)
{

  if(m_use_ocl == false)
  {
    dp->dataset->work();
  }

  dp->dataset->postMergeFillDiscs(dp->msgraph);

  dp->dataset->clear_fnref();

  if(m_compute_out_of_core)
  {
    dp->dataset->clear();
  }

  std::stringstream ss;

  ss<<"dp_"/*<<dp->m_pieceno%num_parallel*/<<"_disc_";

  dp->msgraph->write_discs(ss.str());

  if(m_compute_out_of_core)
  {
    delete dp->msgraph;
  }
}

void GridDataManager::collectManifoldsInRange(uint start,uint end)
{

  uint threadno = 0;

  for ( int j = start; j < end; j += 1)
  {
    GridDataPiece *dp  = m_pieces[j];

    if(m_use_ocl)
      dp->dataset->work_ocl(false);

    if(m_single_threaded_mode == false)
    {
      _LOG("Kicking off collect manifolds "<<j );

      m_threads[threadno++] = new boost::thread
                            ( boost::bind ( &GridDataManager::collectManifold,this,dp ));
    }
    else
    {
      _LOG("Collect manifolds "<<j );

      collectManifold(dp);
    }
  }

}


void GridDataManager::collectSubdomainManifolds( )
{

  std::ifstream data_stream;

  GridDataset::cell_fn_t * data_buffer
      = new GridDataset::cell_fn_t[getMaxDataBufItems()];

  uint num_subdomains =  ((0x01)<<m_num_levels);

  int num_pc_per_buf = std::min(num_parallel,num_subdomains);

  for(uint i = 0 ;i < num_subdomains ; i += num_pc_per_buf)
  {
    if(i%2 == 0 && m_num_levels != 0 )
      finalMergeDownPiecesInRange(i,i+num_pc_per_buf);

    // prepare datasets
    readDataAndInit(data_stream,data_buffer,i);

    if(i%2 == 0 && m_num_levels != 0)
      waitForThreadsInRange(i,i+(num_pc_per_buf+1)/2);

    for(uint j = i; j< i+num_pc_per_buf;++j)
    {
      // collect the manifold info into msgraphs
      collectManifoldsInRange(j,j+1);

      // wait for collection to complete
      waitForThreadsInRange(j,j+1);
    }
  }

  data_stream.close();

  delete data_buffer;
}


void GridDataManager::computeSubdomainMsgraphs ( )
{

  std::ifstream data_stream;

  GridDataset::cell_fn_t * data_buffers[2];

  data_buffers[0] = new GridDataset::cell_fn_t[getMaxDataBufItems()];
  data_buffers[1] = new GridDataset::cell_fn_t[getMaxDataBufItems()];

  uint num_subdomains =  ((0x01)<<m_num_levels);

  int num_pc_per_buf = std::min(num_parallel,num_subdomains);

  for(uint i = 0 ;i < num_subdomains ; i += num_pc_per_buf)
  {

    uint active_buffer = (i/num_pc_per_buf)%2;

    readDataAndInit(data_stream,data_buffers[active_buffer],i);

    if(i != 0 )
    {
      waitForThreadsInRange(i-num_pc_per_buf,i);
    }

    computeMsGraphInRange(i,i+num_pc_per_buf);
  }

  waitForThreadsInRange(num_subdomains-num_pc_per_buf,num_subdomains);

  delete data_buffers[0];
  delete data_buffers[1];

  data_stream.close();
}

void GridDataManager::readDataAndInit
    (std::ifstream &data_stream,
     GridDataset::cell_fn_t *buffer,
     uint start_offset)
{

  uint num_subdomains =  ((0x01)<<m_num_levels);

  int num_pc_per_buf = std::min(num_parallel,num_subdomains);

  GridDataset::rect_size_t domain_sz(m_size_x,m_size_y);

  if(start_offset == 0)
    data_stream.open(m_filename.c_str(),fstream::in|fstream::binary );
  else
    data_stream.seekg(-3*domain_sz[scan_axes]*sizeof ( GridDataset::cell_fn_t),std::ios::cur);

  if(data_stream.is_open() == false)
    throw std::runtime_error("unable to read stream");

  GridDataPiece * dp1 = m_pieces[start_offset];

  GridDataPiece * dp2 = m_pieces[start_offset+num_pc_per_buf-1];

  GridDataset::rect_point_t bl = dp1->dataset->get_ext_rect().bottom_left()/2;
  GridDataset::rect_point_t tr = dp2->dataset->get_ext_rect().top_right()/2;

  size_t num_data_items = (domain_sz[scan_axes]*(tr[split_axes] - bl[split_axes]+1));

  if(num_data_items > getMaxDataBufItems())
    throw std::range_error("insufficient buffer space calculted");

  data_stream.read ( reinterpret_cast<char *> ( buffer),
                sizeof ( GridDataset::cell_fn_t)*num_data_items );


  for(uint j = start_offset ; j < start_offset+num_pc_per_buf;++j)
  {
    GridDataPiece * dp = m_pieces[j];

    GridDataset::rect_point_t dp_bl = dp->dataset->get_ext_rect().bottom_left()/2;

    uint data_offset = domain_sz[scan_axes]*(dp_bl[split_axes] - bl[split_axes]);

    dp->dataset->init_fnref(buffer + data_offset);
  }

}

uint GridDataManager::getMaxDataBufItems()
{
  GridDataset::rect_size_t domain_sz(m_size_x,m_size_y);

  uint num_subdomains =  ((0x01)<<m_num_levels);

  int num_pc_per_buf = std::min(num_parallel,num_subdomains);

  uint max_size_split_axis = (domain_sz[split_axes] + num_subdomains-1)/num_subdomains;

  uint max_sz_per_domain = (max_size_split_axis+3)*domain_sz[scan_axes];

  uint buf_ct = max_sz_per_domain*num_pc_per_buf - (num_pc_per_buf-1)*(domain_sz[scan_axes]+1);

  return buf_ct;
}

GridDataManager::GridDataManager
    ( std::string filename,
      u_int        size_x,
      u_int        size_y,
      u_int        num_levels,
      bool         single_threaded_mode,
      bool         use_ocl,
      double       simp_tresh,
      bool         compute_out_of_core,
      uint         np):
    m_filename(filename),
    m_size_x(size_x),
    m_size_y(size_y),
    m_num_levels(num_levels),
    m_single_threaded_mode(single_threaded_mode),
    m_use_ocl(use_ocl),
    m_simp_tresh(simp_tresh),
    m_compute_out_of_core(compute_out_of_core),
    num_parallel(np)
{

  if(num_parallel == 1)
    m_single_threaded_mode = true;

  if (m_single_threaded_mode == true)
    num_parallel = 1;

  m_threads = new boost::thread*[num_parallel];

  createDataPieces();

  if(m_use_ocl)
    GridDataset::init_opencl();

  Timer t;
  t.start();

  _LOG ( "==========================" );
  _LOG ( "Starting Processing Peices" );
  _LOG ( "--------------------------" );

  if(m_num_levels == 0 && m_compute_out_of_core)
    m_compute_out_of_core = false;

  computeSubdomainMsgraphs ();

  _LOG("timer_time = "<<t.getElapsedTimeInMilliSec());

  mergePiecesUp();

  _LOG("timer_time = "<<t.getElapsedTimeInMilliSec());

  if(m_simp_tresh > 0.0)
    m_pieces[m_pieces.size()-1]->msgraph->simplify_un_simplify(m_simp_tresh);

  _LOG("timer_time = "<<t.getElapsedTimeInMilliSec());

  mergePiecesDown();

  _LOG("timer_time = "<<t.getElapsedTimeInMilliSec());

  collectSubdomainManifolds();

  _LOG("timer_time = "<<t.getElapsedTimeInMilliSec());

  _LOG ( "--------------------------" );
  _LOG ( "Finished Processing peices" );
  _LOG ( "==========================" );

  if ( m_compute_out_of_core == true )
  {
    if(m_use_ocl)
      GridDataset::stop_opencl();

    exit(0);
  }

  delete []m_threads;
}

void GridDataManager ::logAllConnections(const std::string &prefix)
{

  for(uint i = 0 ; i <m_pieces.size();++i)
  {
    GridDataPiece *dp = m_pieces[i];

    std::string filename(prefix+dp->label()+string(".txt"));

    ofstream outfile;
    outfile.open(filename.c_str(),  std::ios::out|std::ios::trunc);

    if(outfile.is_open() == false )
    {
      _LOG("failed to open log file");
      break;
    }

    std::stringstream ss;

    print_connections( (ostream&)ss,*dp->msgraph);

    outfile<<ss.str();
  }

}

void GridDataManager::logAllCancelPairs(const std::string &prefix)
{
}

GridDataManager::~GridDataManager()
{
}

GridDataPiece::GridDataPiece (uint pno):
    dataset(NULL),
    msgraph(NULL),
    level(0),
    m_pieceno(pno)
{
}

std::string GridDataPiece::label()
{
  std::stringstream ss;
  ss<<m_pieceno;

  return ss.str();
}
