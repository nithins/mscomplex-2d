varying in vec3  f_wc_p0;
varying in vec3  f_wc_p1;
varying in float f_wc_radius;
varying in vec3  f_wc_prism_crd;

#define NUM_LIGHTS 1

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

void main()
{

  vec3  v  = f_wc_prism_crd;
  vec3  p  = f_wc_p0;
  vec3  q  = f_wc_p1;
  float r  = f_wc_radius;
  vec3  s  = q-p;

  vec3 vs  = vec3(determinant(mat2(v.yz,s.yz)),
                  determinant(mat2(v.zx,s.zx)),
                  determinant(mat2(v.xy,s.xy)));

  vec3 sp  = vec3(determinant(mat2(s.yz,p.yz)),
                  determinant(mat2(s.zx,p.zx)),
                  determinant(mat2(s.xy,p.xy)));


  float a  = dot(vs,vs);
  float b  = 2.0*dot(vs,sp);
  float c  = dot(sp,sp) - dot(s,s)*r*r ;
  float d  = b*b - 4.0*a*c;

  if(d < 0.0) discard;

  d        = sqrt(d);

  vec3 pt  = ((-b - d)/(2.0*a))*v;
  float t  = -1.0*dot(p-pt,s)/dot(s,s);

  float iv_n = 1.0;

  if( t > 1.0 || t< 0.0) 
  {
    pt     = ((-b + d)/(2.0*a))*v;
    t      = -1.0*dot(p-pt,s)/dot(s,s);
    iv_n   = -1.0;
  }

  if( t > 1.0 || t< 0.0) 
  {
    discard;
  }
  
  vec3 n   = iv_n*normalize(pt - p - t*s);
  
  vec4  dc_pt     = gl_ProjectionMatrix*vec4(pt,1.0); 
  float fragdepth = (dc_pt.z/dc_pt.w +1.0)/2.0;

  vec4 amb  = vec4(0.0);
  vec4 diff = vec4(0.0);
  vec4 spec = vec4(0.0);
  
  for(int i = 0 ; i < NUM_LIGHTS; ++i)
  {
    if(gl_LightSource[i].position.w == 0.0)
      DirectionalLight(i,n,amb,diff,spec);
    else
      PointLight(i,vec3(0,0,0),pt,n,amb,diff,spec);
  }
  
  vec4 color = vec4(0);
  
  color +=  amb*gl_Color;
  color += diff*gl_Color;
  color += spec*gl_FrontMaterial.specular;

  
  gl_FragDepth   =  fragdepth;
  gl_FragColor   =  color;
}
