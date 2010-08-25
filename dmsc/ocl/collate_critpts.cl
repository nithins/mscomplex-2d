#define BLOCK_DIM_X 16
#define BLOCK_DIM_Y 16
#define BLOCK_DIM 128

// #include <common_funcs.cl>

__kernel void collate_cps_initcount(
__read_only   image2d_t  cell_own_img, 
__global unsigned int*   cell_cp_idx_map,
const short2 ext_bl,
const short2 ext_tr
)
{
  short2 c,bb_ext_sz;

  bb_ext_sz.x = ext_tr.x-ext_bl.x;
  bb_ext_sz.y = ext_tr.y-ext_bl.y;

  if(get_global_id(0) > bb_ext_sz.x ||
     get_global_id(1) > bb_ext_sz.y)
   return;

  c.x = get_global_id(0);
  c.y = get_global_id(1);
  
  short2 c_own = read_from_owner_image(c,cell_own_img,ext_bl);
  
  cell_cp_idx_map[c.y*(bb_ext_sz.x+1) + c.x] = ((c.x == c_own.x && c.y == c_own.y)?(1):(0));
} 

__kernel void collate_cps_writeids(
__read_only  image2d_t  cell_own_img,
__global unsigned int*  cell_cp_idx_map,
__global short*         critpt_cellid,
const short2 ext_bl,
const short2 ext_tr
)
{
  short2 c,bb_ext_sz;

  bb_ext_sz.x = ext_tr.x-ext_bl.x;
  bb_ext_sz.y = ext_tr.y-ext_bl.y;

  if(get_global_id(0) > bb_ext_sz.x ||
     get_global_id(1) > bb_ext_sz.y)
   return;

  c.x = get_global_id(0);
  c.y = get_global_id(1);
  
  short2 c_own = read_from_owner_image(c,cell_own_img,ext_bl);

  if(c.x == c_own.x && c.y == c_own.y)
  {
    // write out cellid in the cellid array
    unsigned int idx = cell_cp_idx_map[c.y*(bb_ext_sz.x+1) + c.x];
    critpt_cellid[2*idx + 0] = c.x+ext_bl.x;
    critpt_cellid[2*idx + 1] = c.y+ext_bl.y;
  }
}

__kernel void count_critpt_incidences(
__global short* critpt_cellid,
__read_only image2d_t  cell_own_image,
__global unsigned int* incidence_ct,
unsigned int critpt_ct,
const short2 ext_bl,
const short2 ext_tr
)
{
  int idx = get_global_id(0);
  
  if(idx >= critpt_ct)
    return;

  short2 c,bb_ext_sz;
  c.x = critpt_cellid[2*idx +0]-ext_bl.x;
  c.y = critpt_cellid[2*idx +1]-ext_bl.y;

  bb_ext_sz.x = ext_tr.x-ext_bl.x;
  bb_ext_sz.y = ext_tr.y-ext_bl.y;

  unsigned int incidence_ct_out = 0;

  if(get_cell_dim(c) == 1)
  {
    short2 ncells[4];

    get_cell_facets(c,ncells);
    get_cell_cofacets(c,ncells+2);

    for(int i = 0 ; i < 4;++i)
    {
      if(is_cell_outside_true_boundry(ncells[i],bb_ext_sz) == 1)
        continue;

      short2 n_own = read_from_owner_image(ncells[i],cell_own_image,ext_bl);

      if(is_cell_outside_true_boundry(n_own,bb_ext_sz) == 0)
        incidence_ct_out++;      
    }
  }
  incidence_ct[idx] = incidence_ct_out;
}

__kernel void write_critpt_incidences(
__global short* critpt_cellid,
__read_only image2d_t  cell_own_image,
__global unsigned int* cell_cp_idx_map,
__global unsigned int* incidence_buf_offset,
__global unsigned int* incidence_buf,
unsigned int critpt_ct,
const short2 ext_bl,
const short2 ext_tr
)
{
  int idx = get_global_id(0);

  if(idx >= critpt_ct)
    return;

  short2 c,bb_ext_sz;
  c.x = critpt_cellid[2*idx +0]-ext_bl.x;
  c.y = critpt_cellid[2*idx +1]-ext_bl.y;

  bb_ext_sz.x = ext_tr.x-ext_bl.x;
  bb_ext_sz.y = ext_tr.y-ext_bl.y;
  
  unsigned int out_offset= incidence_buf_offset[idx];

  if(get_cell_dim(c) == 1)
  {
    short2 ncells[4];

    get_cell_facets(c,ncells);
    get_cell_cofacets(c,ncells+2);

    for(int i = 0 ; i < 4;++i)
    {
      if(is_cell_outside_true_boundry(ncells[i],bb_ext_sz) == 1)
        continue;

      short2 n_own = read_from_owner_image(ncells[i],cell_own_image,ext_bl);

      if(is_cell_outside_true_boundry(n_own,bb_ext_sz) == 0)
      {
        unsigned int incident_idx = cell_cp_idx_map[n_own.y*(bb_ext_sz.x+1) + n_own.x];

        incidence_buf[out_offset] = incident_idx;
        
        out_offset++;
      }
    }
  }
}

__kernel void collate_cps_reduce(__global int*  critpt_ct,int n)
{
    // perform first level of reduction,
    // reading from global memory, writing to shared memory
    unsigned int tid = get_local_id(0);
    unsigned int i = get_group_id(0)*(get_local_size(0)*2) + get_local_id(0);
    
    __local int sdata[BLOCK_DIM];

    sdata[tid] = (i < n) ? critpt_ct[i] : 0;
    if (i + get_local_size(0) < n) 
        sdata[tid] += critpt_ct[i+get_local_size(0)];  

    barrier(CLK_LOCAL_MEM_FENCE);

    // do reduction in shared mem
    for(unsigned int s=get_local_size(0)/2; s>0; s>>=1) 
    {
        if (tid < s) 
        {
            sdata[tid] += sdata[tid + s];
        }
        barrier(CLK_LOCAL_MEM_FENCE);
    }

    // write result for this block to global mem 
    if (tid == 0) critpt_ct[get_group_id(0)] = sdata[0];
}