varying vec3 normal;

void main()
{
  gl_Position   = gl_Vertex;
  gl_FrontColor = gl_Color;
  normal        = gl_Normal;
}