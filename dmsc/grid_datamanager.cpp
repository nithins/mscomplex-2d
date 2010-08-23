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
#include <logutil.h>

#include <grid_datamanager.h>

#include <grid_mscomplex.h>


namespace grid
{

  using namespace std;

  void data_manager_t::createDataPieces ()
  {
    rect_t r(cellid_t(0,0),(m_size-cellid_t::one)*2);
    rect_t e(cellid_t(0,0),(m_size-cellid_t::one)*2);

    createPieces_quadtree(r,e,m_num_levels);

    for(uint i = 0 ; i < m_num_subdomains; ++i)
      m_pieces.push_back(new datapiece_t);

    std::reverse(m_pieces.begin(),m_pieces.end());

    for(uint i = 0 ; i < m_num_subdomains*2; ++i)
      m_pieces[i]->m_pieceno = i;

    return;
  }

  const unsigned char split_axes = 1;
  const unsigned char scan_axes  = (split_axes+1)%2;

  void data_manager_t::createPieces_quadtree(rect_t r,rect_t e,u_int level )
  {
    if(level == 0)
    {
      datapiece_t *dp = new datapiece_t;
      dp->dataset = new dataset_t(r,e);
      dp->msgraph = new mscomplex_t(r,e);
      m_pieces.push_back(dp);

      return;
    }

    u_int dim = split_axes;

    rect_size_t  s = r.size();

    rect_point_t tr1 = r.upper_corner();
    rect_point_t bl1 = r.lower_corner();

    tr1[dim] -= 2*(s[dim]/4);
    bl1[dim]  = tr1[dim];

    rect_t r1 = rect_t(r.lower_corner(),tr1);
    rect_t r2 = rect_t(bl1,r.upper_corner());

    rect_point_t tr2 = e.upper_corner();
    rect_point_t bl2 = e.lower_corner();

    tr2[dim] = tr1[dim] + 2;
    bl2[dim] = bl1[dim] - 2;

    rect_t e1 = rect_t(e.lower_corner(),tr2);
    rect_t e2 = rect_t(bl2,e.upper_corner());

    createPieces_quadtree(r1,e1,level-1);
    createPieces_quadtree(r2,e2,level-1);
  }

  void read_msgraph_from_archive(datapiece_t * dp)
  {
    std::string filename("dp_msgraph_");
    filename += dp->label();

    std::ifstream ifs(filename.c_str());

    boost::archive::binary_iarchive ia(ifs);

    dp->msgraph = new mscomplex_t;

    ia >> (*dp->msgraph);
  }

  void write_msgraph_to_archive(datapiece_t * dp)
  {
    std::string filename("dp_msgraph_");
    filename += dp->label();

    std::ofstream ofs(filename.c_str());

    boost::archive::binary_oarchive oa(ofs);

    oa << (*dp->msgraph);

    delete dp->msgraph;
    dp->msgraph = NULL;

  }

  void data_manager_t::computeMsGraph( datapiece_t *dp )
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

    if(m_compute_out_of_core)
    {
      write_msgraph_to_archive(dp);
    }
  }

  void mergePiecesUp_worker
      ( datapiece_t  *dp,
        datapiece_t  *dp1,
        datapiece_t  *dp2,
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
    dp->msgraph       = mscomplex_t::merge_up(*dp1->msgraph,*dp2->msgraph);

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
      ( datapiece_t  * dp,
        datapiece_t  * dp1,
        datapiece_t  * dp2,
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


  void data_manager_t::computeMsGraphInRange(uint start,uint end )
  {
    _LOG ( "Begin Gradient/ conn work for pieces from "<<start<<" to "<<end );

    for ( int i = end-1 ; i >= start;i-- )
    {
      datapiece_t * dp = m_pieces[i];

      if(m_use_ocl == true)
        dp->dataset->work_ocl();

      if(m_single_threaded_mode == false)
      {
        _LOG ( "Kicking off thread "<<i );


        m_threads[i] = new boost::thread
                             ( boost::bind ( &data_manager_t::computeMsGraph,this,dp ) );
      }
      else
      {
        computeMsGraph(dp);
      }

    }
  }

  void data_manager_t::waitForPiecesInRange(int start,int end )
  {
    if(m_single_threaded_mode == false)
    {
      if(start < 1) return;

      for ( int i = end-1 ; i >= start;i-- )
      {
        m_threads[i]->join();

        delete m_threads[i];

        m_threads[i] = NULL;
        _LOG ( "thread "<<i<<" joint" );
      }
    }
  }

  void data_manager_t::clearDatasetsInRange(uint start ,uint end)
  {
    if(m_compute_out_of_core == true)
    {
      for ( int i = end-1 ; i >= start;i-- )
      {
        m_pieces[i]->dataset->clear();
        m_pieces[i]->dataset->clear_fnref();
      }
    }
  }

  inline int npot(int val)
  {
    val--;
    val = (val >> 1) | val;
    val = (val >> 2) | val;
    val = (val >> 4) | val;
    val = (val >> 8) | val;
    val = (val >> 16) | val;
    val++;

    return val;
  }

  void data_manager_t::mergePiecesUp( )
  {

    // 1 based indexing.. backwards
    for ( int i = m_num_subdomains ; i>1;)
    {
      int stride = std::min((int)num_parallel,i/2);

      for(int j = 1 ; j <=stride;++j)
      {
        int p   = i - j;
        int c1  = p*2+1;
        int c2  = p*2;

        bool src_archived = m_compute_out_of_core;
        bool archive_dest = m_compute_out_of_core&&(p != 1);

        if(m_single_threaded_mode == false)
        {
          _LOG("Kicking off Merge Up "<<c1<<" "<<c2<<"->"<<p);

          m_threads[p] =
              new boost::thread
              ( boost::bind
                ( &mergePiecesUp_worker,m_pieces[p],m_pieces[c1],m_pieces[c2],
                  src_archived,archive_dest )
                );

        }
        else
        {
          _LOG("Kicking off Merge Up "<<c1<<" "<<c2<<"->"<<p);

          mergePiecesUp_worker
              (m_pieces[p],m_pieces[c1],m_pieces[c2],src_archived,archive_dest);
        }
      }
      i -= stride;

      waitForPiecesInRange(i,i+stride);
    }
  }

  void data_manager_t::mergePiecesDown()
  {
    // 1 based indexing.. forward
    for ( int i = 1; i < m_num_subdomains/2;)
    {
      int stride = std::min((int)num_parallel,i);

      mergeDownPiecesInRange(i,i+stride);

      waitForPiecesInRange(i,i+stride);

      i += stride;
    }
  }

  void data_manager_t::mergeDownPiecesInRange(int start,int end)
  {
    if(start < 1) return;

    for ( int j = start; j < end && j >0 ; j++)
    {

      bool src_archived  = m_compute_out_of_core && (j != 1);
      bool dest_archived = m_compute_out_of_core;

      bool archive_src   = m_compute_out_of_core && (j < m_num_subdomains);
      bool archive_dest  = m_compute_out_of_core && (2*j < m_num_subdomains);


      int p   = j;
      int c1  = 2*p;
      int c2  = 2*p+1;

      if(m_single_threaded_mode == false)
      {
        _LOG("Kicking off Merge Down "<<p<<"->"<<c1<<" "<<c2);

        m_threads[p] =
            new boost::thread
            ( boost::bind
              ( &mergePiecesDown_worker,m_pieces[p],m_pieces[c1],m_pieces[c2],
                src_archived,dest_archived,archive_src,archive_dest )
              );
      }
      else
      {
        _LOG("Kicking off Merge Down "<<p<<"->"<<c1<<" "<<c2);

        mergePiecesDown_worker
            (m_pieces[p],m_pieces[c1],m_pieces[c2],
             src_archived,dest_archived,archive_src,archive_dest);
      }
    }
  }

  void data_manager_t::collectManifold( datapiece_t  * dp)
  {

    if(m_use_ocl == false)
    {
      dp->dataset->work();
    }

    dp->dataset->postMergeFillDiscs(dp->msgraph);

    if(m_save_mfolds_to_file)
    {
      std::stringstream ss;

      ss<<"dp_disc_";

      dp->msgraph->write_discs(ss.str());
    }

    if(m_compute_out_of_core)
    {
      delete dp->msgraph;
    }
  }

  void data_manager_t::collectManifoldsInRange(uint start,uint end)
  {
    for ( int j = start; j < end; j += 1)
    {
      datapiece_t *dp  = m_pieces[j];

      if(m_use_ocl)
        dp->dataset->work_ocl(false);

      if(m_single_threaded_mode == false)
      {
        _LOG("Kicking off collect manifolds "<<j );

        m_threads[j] = new boost::thread
                                ( boost::bind ( &data_manager_t::collectManifold,this,dp ));
      }
      else
      {
        _LOG("Collect manifolds "<<j );

        collectManifold(dp);
      }
    }

  }


  void data_manager_t::collectSubdomainManifolds( )
  {

    std::ifstream data_stream;

    cell_fn_t * data_buffer = new cell_fn_t[getMaxDataBufItems()];

    uint m_num_subdomains =  ((0x01)<<m_num_levels);

    int stride = std::min(num_parallel,m_num_subdomains);

    int i = 2*m_num_subdomains;

    mergeDownPiecesInRange((i-stride)/2,i/2);

    readDataAndInit(data_stream,data_buffer,i);

    waitForPiecesInRange((i-stride)/2,i/2);

    collectManifoldsInRange(i-stride,i);

    waitForPiecesInRange(i-stride,i);

    for(i -= stride;i >m_num_subdomains ; i -= stride)
    {

      mergeDownPiecesInRange((i-stride)/2,i/2);

      clearDatasetsInRange(i,i+stride);

      readDataAndInit(data_stream,data_buffer,i);

      waitForPiecesInRange((i-stride)/2,i/2);

      collectManifoldsInRange(i-stride,i);

      waitForPiecesInRange(i-stride,i);

    }

    clearDatasetsInRange(i,i+stride);

    data_stream.close();

    delete []data_buffer;
  }


  void data_manager_t::computeSubdomainMsgraphs ( )
  {
    std::ifstream data_stream;

    cell_fn_t * data_buffers[2];

    data_buffers[0] = new cell_fn_t[getMaxDataBufItems()];
    data_buffers[1] = new cell_fn_t[getMaxDataBufItems()];

    int stride = std::min(num_parallel,m_num_subdomains);

    int i = 2*m_num_subdomains;

    readDataAndInit(data_stream,data_buffers[(i/stride)%2],i);

    computeMsGraphInRange(i-stride,i);

    for(i -= stride ;i > m_num_subdomains; i -= stride)
    {
      readDataAndInit(data_stream,data_buffers[(i/stride)%2],i);

      waitForPiecesInRange(i,i+stride);

      computeMsGraphInRange(i-stride,i);

      clearDatasetsInRange(i,i+stride);
    }

    waitForPiecesInRange(i,i+stride);

    clearDatasetsInRange(i,i+stride);

    delete []data_buffers[0];
    delete []data_buffers[1];

    data_stream.close();
  }

  void data_manager_t::readDataAndInit
      (std::ifstream &data_stream,cell_fn_t *buffer,uint start_offset)
  {
    int stride = std::min(num_parallel,m_num_subdomains);

    rect_size_t domain_sz(m_size);

    if(start_offset == m_num_subdomains*2)
      data_stream.open(m_filename.c_str(),fstream::in|fstream::binary );
    else
      data_stream.seekg(-3*domain_sz[scan_axes]*sizeof ( cell_fn_t),std::ios::cur);

    if(data_stream.is_open() == false)
      throw std::runtime_error("unable to read stream");

    rect_point_t bl =
        m_pieces[start_offset-1]->dataset->get_ext_rect().lower_corner()/2;

    rect_point_t tr =
        m_pieces[start_offset-stride]->dataset->get_ext_rect().upper_corner()/2;

    size_t num_data_items = (domain_sz[scan_axes]*(tr[split_axes] - bl[split_axes]+1));

    if(num_data_items > getMaxDataBufItems())
      throw std::range_error("insufficient buffer space calculted");

    data_stream.read ( reinterpret_cast<char *> ( buffer),
                       sizeof ( cell_fn_t)*num_data_items );


    for(int j = start_offset-1 ; j >= start_offset-stride;--j)
    {
      datapiece_t * dp = m_pieces[j];

      rect_point_t dp_bl = dp->dataset->get_ext_rect().lower_corner()/2;

      uint data_offset = domain_sz[scan_axes]*(dp_bl[split_axes] - bl[split_axes]);

      dp->dataset->init();

      dp->dataset->init_fnref(buffer + data_offset);
    }
  }

  uint data_manager_t::getMaxDataBufItems()
  {
    rect_size_t domain_sz(m_size);

    uint num_subdomains =  ((0x01)<<m_num_levels);

    int num_pc_per_buf = std::min(num_parallel,num_subdomains);

    uint max_size_split_axis = (domain_sz[split_axes] + num_subdomains-1)/num_subdomains;

    uint max_sz_per_domain = (max_size_split_axis+3)*domain_sz[scan_axes];

    uint buf_ct = max_sz_per_domain*num_pc_per_buf - (num_pc_per_buf-1)*(domain_sz[scan_axes]+1);

    return buf_ct;
  }

  data_manager_t::data_manager_t
      ( std::string  filename,
        cellid_t     size,
        u_int        num_levels,
        bool         threaded_mode,
        bool         use_ocl,
        double       simp_tresh,
        bool         compute_out_of_core,
        uint         np,
        bool         save_mfolds_to_file):
      m_filename(filename),
      m_size(size),
      m_num_levels(num_levels),
      m_single_threaded_mode(threaded_mode),
      m_use_ocl(use_ocl),
      m_simp_tresh(simp_tresh),
      m_compute_out_of_core(compute_out_of_core),
      num_parallel(np),
      m_save_mfolds_to_file(save_mfolds_to_file),
      m_num_subdomains(((0x01)<<num_levels))
  {

    if(num_parallel == 1)
      m_single_threaded_mode = true;

    if (m_single_threaded_mode == true)
      num_parallel = 1;

    m_threads = new boost::thread*[(0x02)<<m_num_levels];

    createDataPieces();

    if(m_use_ocl)
      dataset_t::init_opencl();

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
      m_pieces[1]->msgraph->simplify_un_simplify(m_simp_tresh);

    _LOG("timer_time = "<<t.getElapsedTimeInMilliSec());

    mergePiecesDown();

    _LOG("timer_time = "<<t.getElapsedTimeInMilliSec());

    collectSubdomainManifolds();

    _LOG("timer_time = "<<t.getElapsedTimeInMilliSec());

    _LOG ( "--------------------------" );
    _LOG ( "Finished Processing peices" );
    _LOG ( "==========================" );

    if(m_use_ocl)
      dataset_t::stop_opencl();

    if ( m_compute_out_of_core == true )
      exit(0);

    delete []m_threads;
  }

  void data_manager_t ::logAllConnections(const std::string &prefix)
  {

    for(uint i = 0 ; i <m_pieces.size();++i)
    {
      datapiece_t *dp = m_pieces[i];

      std::string filename(prefix+dp->label()+string(".txt"));

      ofstream outfile;
      outfile.open(filename.c_str(),  std::ios::out|std::ios::trunc);

      if(outfile.is_open() == false )
      {
        _LOG("failed to open log file");
        break;
      }

      std::stringstream ss;

      dp->msgraph->print_connections( (ostream&)ss);

      outfile<<ss.str();
    }

  }

  void data_manager_t::logAllCancelPairs(const std::string &prefix)
  {
  }

  data_manager_t::~data_manager_t()
  {
  }

  datapiece_t::datapiece_t ():
      dataset(NULL),
      msgraph(NULL),
      level(0),
      m_pieceno(-1)
  {
  }

  std::string datapiece_t::label()
  {
    std::stringstream ss;
    ss<<m_pieceno;

    return ss.str();
  }
}
