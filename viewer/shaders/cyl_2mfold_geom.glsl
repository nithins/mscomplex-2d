
#version 120
#extension GL_EXT_geometry_shader4 : enable
#extension GL_EXT_gpu_shader4 : enable
#extension GL_ARB_texture_rectangle: enable

//HEADER_REPLACE_BEGIN

const int is_dual = 1;

//HEADER_REPLACE_END

const float small_sz = 0.1;
const float big_sz   = 1.0;

uniform sampler2DRect rawdata_texture;

const float ug_cylinder_radius = 0.05;

varying out vec3  f_wc_p0;
varying out vec3  f_wc_p1;
varying out float f_wc_radius;
varying out vec3  f_wc_prism_crd;

vec3[4] get_quad(vec3 c)
{
  vec3[4] p;  int pos = 0;  ivec3 i = ivec3(0,0,0);ivec3 ex = ivec3(0,0,0); vec3  sz = vec3(0,0,0);
  
  ex.x     = (int(c.x)&1)^is_dual;
  ex.z     = (int(c.z)&1)^is_dual;
  
  sz.x     = (ex.x == 1)?(big_sz):(small_sz);
  sz.z     = (ex.z == 1)?(big_sz):(small_sz);

  for(i[2] = -1 ; i[2] <= 1 ;i[2]+=2)
      for(i[0] = -1 ; i[0] <= 1 ;i[0]+=2)
      {
        p[pos] 	 = c+i*sz;
        p[pos].y = texture2DRect(rawdata_texture, ((c+i*ex).xz)/2).x;
        pos++;
      }

  return p;
}

void draw_cylinder(vec3 cyl_pt1,vec3 cyl_pt2)
{
  float[2] bin01;
  bin01[0] = -1.0;
  bin01[1] =  1.0;

  f_wc_p0      = (gl_ModelViewMatrix*vec4(cyl_pt1,1.0)).xyz;
  f_wc_p1      = (gl_ModelViewMatrix*vec4(cyl_pt2,1.0)).xyz;

  f_wc_radius  = length(gl_ModelViewMatrix*vec4(ug_cylinder_radius,0,0,0));

  vec3  wc_cyl_axis  = f_wc_p1 - f_wc_p0;

  if(wc_cyl_axis.z < 0.0)
    wc_cyl_axis *= -1.0;

  float[3] wc_axes_len;
  wc_axes_len[0] = f_wc_radius;
  wc_axes_len[1] = length(wc_cyl_axis)/2.0;
  wc_axes_len[2] = f_wc_radius;

  mat3 wc_axes;
  wc_axes[1]   = normalize(wc_cyl_axis);
  wc_axes[2]   = vec3(0,0,1);
  if(abs(dot(wc_axes[1],wc_axes[2])) > 0.99)
    wc_axes[2] = vec3(0,1,0);
  wc_axes[0] = normalize(cross(wc_axes[1],wc_axes[2]));
  wc_axes[2] = normalize(cross(wc_axes[0],wc_axes[1]));

  vec3 wc_center = (f_wc_p0+f_wc_p1)/2.0;

  int[3] valid_face_version;

  vec3 wc_eye_dir = normalize(wc_center);

  for(int j = 0 ; j < 3 ; j ++)
  {
    if(dot(wc_axes[(j+2)%3],wc_eye_dir) >0.0) 
      valid_face_version[j] = 0; 
    else 
      valid_face_version[j] = 1;
  }

  f_wc_prism_crd = vec3(0,0,0);

  for(int j = 0 ; j < 3 ; j ++)
  {
    int k = valid_face_version[j];

    vec3 r_dir = wc_axes_len[j]*wc_axes[j];
    vec3 u_dir = bin01[k]*wc_axes_len[(j+1)%3]*wc_axes[(j+1)%3];
    vec3 n_dir = bin01[k]*wc_axes_len[(j+2)%3]*wc_axes[(j+2)%3];	

    for(int l = 0 ; l < 4; l++)
    {
      f_wc_prism_crd    = wc_center + n_dir + bin01[l%2]*r_dir + bin01[(l/2)%2]*u_dir;
      gl_Position       = gl_ProjectionMatrix*vec4(f_wc_prism_crd,1.0); 
      EmitVertex();
    }
    EndPrimitive();
  }

}

void main()
{
  for(int i=0; i< gl_VerticesIn; i++)
  {
    vec3[4] q = get_quad(gl_PositionIn[i].xyz);

    gl_FrontColor     = gl_FrontColorIn[i];
    
    draw_cylinder(q[0],q[1]);
    draw_cylinder(q[1],q[3]);
    draw_cylinder(q[3],q[2]);
    draw_cylinder(q[2],q[0]);
  } 
}

