#version 120
#extension GL_EXT_geometry_shader4 : enable
#extension GL_EXT_gpu_shader4 : enable
#extension GL_ARB_texture_rectangle: enable

const float large_width  = 0.3;
const float small_width  = 0.1;
const float length_ratio = 0.66;
const float u_len_pt = 0;
const float v_len_pt = 1;

vec4 diffuse;
vec4 ambient;
vec4 in_color;

vec3 lightDir;
vec3 halfVector;

varying in vec3 normal[];

void set_front_color_xfm(vec3 n)
{
  vec3  halfV,viewV,ldir;
  float NdotL,NdotHV;
  vec4  color  = ambient;

  

  NdotL = max(dot(n,lightDir),0.0);

  if (NdotL > 0.0)
  {
    halfV = normalize(halfVector);
    NdotHV = max(dot(n,halfV),0.0);
    if(NdotHV > 0.0)
    {
      color += gl_FrontMaterial.specular *gl_LightSource[0].specular *pow(NdotHV,gl_FrontMaterial.shininess);
      color += diffuse * NdotL;
    }
  }
  gl_FrontColor = in_color*color;
}

void set_light_constants()
{
  diffuse    = gl_LightSource[0].diffuse;
  ambient    = gl_LightSource[0].ambient + gl_LightModel.ambient;
  lightDir   = normalize(vec3(gl_LightSource[0].position));
  halfVector = normalize(gl_LightSource[0].halfVector.xyz);
}

vec3[6] get_arrow(vec3 l[2],vec3 n)
{
  vec3 t[6];

  vec3 len_dir = l[1] - l[0];
  vec3 lat_dir = -normalize(cross(len_dir,n));

  l[0] = l[0] + u_len_pt * len_dir;
  l[1] = l[0] + v_len_pt * len_dir;

  vec3 cntr_pt = mix(l[0],l[1],length_ratio);    

  t[0] = l[0];
  t[1] = cntr_pt - small_width*lat_dir;
  t[2] = cntr_pt + small_width*lat_dir;

  t[3] = l[1];
  t[4] = cntr_pt + large_width*lat_dir;
  t[5] = cntr_pt - large_width*lat_dir;  

  return t;
}

void draw_tri(vec4 p1,vec4 p2,vec4 p3,vec3 n)
{
  set_front_color_xfm(n);

  gl_Position = p1; EmitVertex();
  gl_Position = p2; EmitVertex();
  gl_Position = p3; EmitVertex();
  EndPrimitive();
}

void draw_arrow(vec3 c[2],vec3 n)
{
  vec3 mc_t[6] = get_arrow(c,n);
  vec4 wc_t[6];

  for(int i = 0 ; i < 6; ++i)  
    wc_t[i] = (gl_ModelViewProjectionMatrix*vec4(mc_t[i],1.0));

  n = normalize(gl_NormalMatrix*n);

  draw_tri(wc_t[0],wc_t[1],wc_t[2],n);
  draw_tri(wc_t[3],wc_t[4],wc_t[5],n);  
}


void main()
{
  set_light_constants();

  vec3 c[2];

  for(int i=0; i< gl_VerticesIn; i += 2)
  {
    in_color = gl_FrontColorIn[i] + gl_FrontColorIn[i+1]/2;
    
    c[0] = gl_PositionIn[i+0].xyz;
    c[1] = gl_PositionIn[i+1].xyz;

    vec3 n = (normal[i] + normal[i+1])/2;

    draw_arrow(c,n);
  } 
}