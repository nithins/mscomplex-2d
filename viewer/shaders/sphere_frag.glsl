varying in vec3  f_wc_quad_coord;
varying in vec3  f_wc_center;
varying in float f_wc_radius;

#define NUM_LIGHTS 2

void DirectionalLight( in int i,
		       in vec3 normal,
		       inout vec4 ambient,
		       inout vec4 diffuse,
		       inout vec4 specular)
{
  float nDotVP;
  float nDotHV;
  float pf;
  
  nDotVP = max(0.0,dot(normal,normalize(vec3(gl_LightSource[i].position))));
  nDotHV = max(0.0,dot(normal,normalize(vec3(gl_LightSource[i].halfVector))));
  
  if(nDotVP == 0.0)
    pf = 0.0;
  else
    pf = pow(nDotHV,gl_FrontMaterial.shininess);
  
  ambient  += gl_LightSource[i].ambient;
  diffuse  += gl_LightSource[i].diffuse*nDotVP;
  specular += gl_LightSource[i].specular*pf;
}

void PointLight(in int i,
		in vec3 eye,
		in vec3 ecPosition3,
		in vec3 normal,
		inout vec4 ambient,
		inout vec4 diffuse,
		inout vec4 specular)
{
  float nDotVP;
  float nDotHV;
  float pf;
  float attenuation;
  float d;
  vec3  VP;
  vec3  halfVector;
  
  VP = vec3(gl_LightSource[i].position)-ecPosition3;
  
  d  = length(VP);
  
  VP = normalize(VP);
  
  attenuation = 1.0/(gl_LightSource[i].constantAttenuation +
		     gl_LightSource[i].linearAttenuation*d+
		     gl_LightSource[i].quadraticAttenuation*d*d);
		     
  halfVector = normalize(VP+eye);
  
  nDotVP = max(0.0,dot(normal,VP));
  nDotHV = max(0.0,dot(normal,halfVector));
  
  if(nDotVP == 0.0)
    pf = 0.0;
  else
    pf = pow(nDotHV,gl_FrontMaterial.shininess);
  
  ambient  += gl_LightSource[i].ambient*attenuation;
  diffuse  += gl_LightSource[i].diffuse*nDotVP*attenuation;
  specular += gl_LightSource[i].specular*pf*attenuation;  
}

const vec4  uf_border_color              = vec4(0.0,0.0,0.0,1.0);
const float uf_border_color_start_cosine = 0.5;
const float uf_use_border_hinting        = 1.0;

void main()
{

  float a = dot(f_wc_quad_coord,f_wc_quad_coord);
  float b = -2.0*dot(f_wc_quad_coord,f_wc_center);
  float c = dot(f_wc_center,f_wc_center) - f_wc_radius*f_wc_radius;
  float d = b*b - 4.0*a*c;

  if(d <0.0) discard;

  d = sqrt(d);

  if(f_wc_radius < 0.0)
    d = -d;

  vec3  wc_sp_pt   = ((-b -d)/(2.0*a))*f_wc_quad_coord;
  vec3  wc_sp_nrml = normalize(wc_sp_pt - f_wc_center);
  
  
  vec4 dc_sp_pt    = gl_ProjectionMatrix*vec4(wc_sp_pt,1.0);
  float fragDepth  = (dc_sp_pt.z/dc_sp_pt.w +1.0)/2.0 ;

  if(fragDepth <= 0.0 && f_wc_radius > 0.0)
  {
    wc_sp_nrml = vec3(0,0,-1);
    fragDepth  = 0.0001*(1.0-abs(fragDepth));
  }

  if(fragDepth>= 0.0 && f_wc_radius < 0.0)
  {
    wc_sp_nrml = vec3(0,0,-1);
    fragDepth  = 0.0001*(1.0-abs(fragDepth));
  }
  
  vec4 amb  = vec4(0.0);
  vec4 diff = vec4(0.0);
  vec4 spec = vec4(0.0);
  
  for(int i = 0 ; i < NUM_LIGHTS; ++i)
  {
    if(gl_LightSource[i].position.w == 0.0)
      DirectionalLight(i,wc_sp_nrml,amb,diff,spec);
    else
      PointLight(i,vec3(0,0,0),wc_sp_pt,wc_sp_nrml,amb,diff,spec);
  }
  
  vec4 color      = vec4(0);
  vec4 poly_color = gl_Color;
  
  if(uf_use_border_hinting == 1.0)
  {
    vec3 wc_eye_dir = normalize(wc_sp_pt);
    float n_dot_e   = abs(dot(wc_sp_nrml,wc_eye_dir));    
    float alpha     = max(uf_border_color_start_cosine - n_dot_e,0.0)/uf_border_color_start_cosine;     
    poly_color      = mix(gl_Color,uf_border_color,0.75*alpha);    
  }
  
  color += (diff + amb)*poly_color + spec*gl_FrontMaterial.specular;

  gl_FragDepth   = fragDepth;
  gl_FragColor   = color;
}
