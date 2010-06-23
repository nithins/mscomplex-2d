#version 120
#extension GL_EXT_geometry_shader4 : enable
#extension GL_ARB_texture_rectangle: enable
#extension GL_EXT_gpu_shader4      : enable

varying out vec3  f_wc_quad_coord;
varying out vec3  f_wc_center;
varying out float f_wc_radius;

vec2 bin01 = vec2(-1.0,1.0);

uniform float g_wc_radius = 0.4;

void main ( void )
{
  for(int i=0; i< gl_VerticesIn; i++)
  {
    mat4 mvm       = gl_ModelViewMatrix;

    float mc_sp_r  = g_wc_radius;

    float wc_sp_r  = length(mvm*vec4(mc_sp_r,0,0,0));
    vec3  wc_sp_c  = (mvm*gl_PositionIn[i]).xyz;

    float wc_sp_r2 = wc_sp_r*wc_sp_r;
    float wc_e_dis = length(wc_sp_c); // dist of the eye pt
//     float wc_m     = wc_sp_r2/wc_e_dis;
//     vec3  wc_bb_c  = (1-(wc_m/wc_e_dis))*wc_sp_c;
//     float wc_bb_r  = sqrt(wc_sp_r2 - wc_m*wc_m);

    vec3  wc_bb_c  =  wc_sp_c + normalize(wc_sp_c)*wc_sp_r;
    float wc_bb_r  = (wc_sp_r*wc_e_dis + wc_sp_r2)/sqrt(wc_e_dis*wc_e_dis - wc_sp_r2);

    f_wc_radius     = (wc_sp_r);
    f_wc_center     = (wc_sp_c).xyz;

    mat3 wc_axes;
    wc_axes[2] = normalize(wc_sp_c.xyz);
    wc_axes[1] = vec3(1,0,0);
    wc_axes[0] = normalize(cross(wc_axes[1],wc_axes[2]));
    wc_axes[1] = normalize(cross(wc_axes[2],wc_axes[0]));

    vec4 dc_sp_center  = gl_ProjectionMatrix*vec4(f_wc_center,1.0);

//   a negative radius indicates that the center is out of the dc view frustum
    if(abs(dc_sp_center.z/dc_sp_center.w ) > 1.0)
      f_wc_radius *=  -1.0;

    for (int k = 0 ; k<4 ; k++)
    {

      f_wc_quad_coord = (wc_bb_c + wc_bb_r*(bin01[(k/2)%2]*wc_axes[0] + bin01[k%2]*wc_axes[1]));
      gl_Position     = gl_ProjectionMatrix*(vec4(f_wc_quad_coord,1.0));
      gl_FrontColor   = gl_FrontColorIn[i];

      EmitVertex();
    }

    EndPrimitive();
  }

}
