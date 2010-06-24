
#version 120
#extension GL_EXT_geometry_shader4 : enable
#extension GL_EXT_gpu_shader4 : enable
#extension GL_ARB_texture_rectangle: enable

//HEADER_REPLACE_BEGIN

const int is_dual = 1;

//HEADER_REPLACE_END

const float line_sz   = 1.0;

uniform sampler2DRect rawdata_texture;
uniform vec2 ug_bl;
uniform vec2 ug_tr;

const float ug_cylinder_radius = 0.11;

varying out vec3  f_wc_p0;
varying out vec3  f_wc_p1;
varying out float f_wc_radius;
varying out vec3  f_wc_prism_crd;

vec3[2] get_line(vec3 c)
{
  vec3[2] p; vec3  sz = vec3(0,0,0); 
  
  sz.x   = (((int(c.x)&1)^is_dual) == 1)?(line_sz):(0.0);
  sz.z   = (((int(c.z)&1)^is_dual) == 1)?(line_sz):(0.0);

  p[0]   = c - sz;
  p[1]   = c + sz;

  p[0].xz = max(p[0].xz,ug_bl);
  p[1].xz = max(p[1].xz,ug_bl);

  p[0].xz = min(p[0].xz,ug_tr);
  p[1].xz = min(p[1].xz,ug_tr);

  p[0].y = texture2DRect(rawdata_texture, (p[0].xz)/2).x;
  p[1].y = texture2DRect(rawdata_texture, (p[1].xz)/2).x;

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
    vec3[2] l = get_line(gl_PositionIn[i].xyz);

    gl_FrontColor     = gl_FrontColorIn[i];
    
    draw_cylinder(l[0],l[1]);
  } 
}

