import processing.core.*; 
import processing.data.*; 
import processing.event.*; 
import processing.opengl.*; 

import processing.pdf.*; 
import java.nio.*; 
import processing.core.PMatrix3D; 

import java.util.HashMap; 
import java.util.ArrayList; 
import java.io.File; 
import java.io.BufferedReader; 
import java.io.PrintWriter; 
import java.io.InputStream; 
import java.io.OutputStream; 
import java.io.IOException; 

public class Swirl extends PApplet {

// LecturesInGraphics: vector interpolation
// Template for sketches
// Author: Jarek ROSSIGNAC

//**************************** global variables ****************************
// Skate dancer on moving terrain
float dz=0; // distance to camera. Manipulated with wheel or when 
//float rx=-0.06*TWO_PI, ry=-0.04*TWO_PI;    // view angles manipulated when space pressed but not mouse
float rx=0, ry=0;    // view angles manipulated when space pressed but not mouse
Boolean twistFree=false, animating=true, tracking=false, center=true, gouraud=true, showControlPolygon=false, showNormals=false;
boolean viewpoint=false;
pt Viewer = P();
pt F = P(0,0,0);  // focus point:  the camera is looking at it (moved when 'f or 'F' are pressed
pt Of=P(100,100,0), Ob=P(110,110,0); // red point controlled by the user via mouseDrag : used for inserting vertices ...
pt Vf=P(0,0,0), Vb=P(0,0,0);
float t=1, f=0, s=0;
Boolean animate=true, linear=false, circular=false, beautiful=false, showFrames=true;
float len=60; // length of arrows
//**************************** initialization ****************************
public void setup() {               // executed once at the begining          
  size(900, 900, P3D); // p3D means that we will do 3D graphics
  frameRate(30);             // render 30 frames per second
  smooth();                  // turn on antialiasing
  myFace = loadImage("data/pic.jpg");  // load image from file pic.jpg in folder data *** replace that file with your pic of your own face
  textureMode(NORMAL);
  P.declare();
  P.loadPts("data/pts");
  }

//**************************** display current frame ****************************
public void draw() {      // executed at each frame
  background(white); // clear screen and paints white background
  if(snapPic) beginRecord(PDF,PicturesOutputPath+"/P"+nf(pictureCounter++,3)+".pdf"); 
  

  F0 = F(P.G[0],P.G[1],P.G[2]); 
  F1 = F(P.G[3],P.G[4],P.G[5]);
  E0  = F(P.G[6],P.G[7],P.G[8]);
  // Show user-controlled Frames 
 
  penFill(green,3);  showArrows(F0); // show the start frame F0 as arrows
  penFill(blue,3);  showArrows(F1); // show the second frame F1 (will be hidden by blue iteration)
    
  penFill(magenta,3);   showArrow(E0);  // show magenta decoration arrow 
 
  F1rF0 = F0.invertedOf(F1);
  E0rF0 = F0.invertedOf(E0);
  FF=F(F0); E=F(E0);
  for(int i=0; i<k; i++) {
     FF=FF.of(F1rF0); 
     if(showFrames) {penFill(blue,3); showArrows(FF);}
     E=FF.of(E0rF0); 
     if(showFrames)  {pen(magenta,3,magenta); showArrow(E);}
     Fk = F(FF);
     }

  penFill(red,3);    showArrows(Fk); // show last frame Fk as red arrows
  Ft = F(F0,t,Fk); penFill(black,3);    showArrows(Ft); // show intermediate frame of morph (F0,t,Fk);
  Et=Ft.of(E0rF0); penFill(magenta,3);   showArrow(Et); 
  
  if(animating) {t+=0.01f; if(t>=1) {t=1; animating=false;}} 


  if(snapPic) {endRecord(); snapPic=false;} // end saving a .pdf of the screen

  fill(black); displayHeader();
  if(scribeText && !filming) displayFooter(); // shows title, menu, and my face & name 
  if(filming && (animating || change)) saveFrame("FRAMES/F"+nf(frameCounter++,4)+".tif"); // saves a movie frame 
  change=false; // to avoid capturing movie frames when nothing happens
  }  // end of draw()
  
//**************************** user actions ****************************
public void keyPressed() { // executed each time a key is pressed: sets the "keyPressed" and "key" state variables, 
	if(key=='`') picking=true; 
	  if(key=='?') scribeText=!scribeText;
	  if(key=='!') snapPicture();
	  if(key=='~') filming=!filming;
	  if(key==']') showControlPolygon=!showControlPolygon;
	  if(key=='|') showNormals=!showNormals;
	  if(key=='G') gouraud=!gouraud;
	  if(key=='q') Q.copyFrom(P);
	  if(key=='p') P.copyFrom(Q);
	  if(key=='e') {PtQ.copyFrom(Q);Q.copyFrom(P);P.copyFrom(PtQ);}
	  if(key=='=') {bu=fu; bv=fv;}
	  // if(key=='.') F=P.Picked(); // snaps focus F to the selected vertex of P (easier to rotate and zoom while keeping it in center)
	  if(key=='c') center=!center; // snaps focus F to the selected vertex of P (easier to rotate and zoom while keeping it in center)
	  if(key=='t') tracking=!tracking; // snaps focus F to the selected vertex of P (easier to rotate and zoom while keeping it in center)
	  if(key=='x' || key=='z' || key=='d') P.setPickedTo(pp); // picks the vertex of P that has closest projeciton to mouse
	  if(key=='d') P.deletePicked();
	  if(key=='i') P.insertClosestProjection(Of); // Inserts new vertex in P that is the closeset projection of O
	  if(key=='W') {P.savePts("data/pts"); Q.savePts("data/pts2");}  // save vertices to pts2
	  if(key=='L') {P.loadPts("data/pts"); Q.loadPts("data/pts2");}   // loads saved model
	  if(key=='w') P.savePts("data/pts");   // save vertices to pts
	  if(key=='l') P.loadPts("data/pts"); 
	  if(key=='a') animating=!animating; // toggle animation
	  if(key==',') viewpoint=!viewpoint;
	  if(key=='#') exit();
	  change=true;
  }

public void mouseWheel(MouseEvent event) {dz += event.getAmount(); change=true;}

public void mousePressed() {
   if (!keyPressed) picking=true;
  }
  
public void mouseMoved() {
  if (keyPressed && key==' ') {rx-=PI*(mouseY-pmouseY)/height; ry+=PI*(mouseX-pmouseX)/width;};
  if (keyPressed && key=='s') dz+=(float)(mouseY-pmouseY); // approach view (same as wheel)
  if (keyPressed && key=='v') { //**<01 
      u+=(float)(mouseX-pmouseX)/width;  u=max(min(u,1),0);
      v+=(float)(mouseY-pmouseY)/height; v=max(min(v,1),0); 
      } 
  }
public void mouseDragged() {
  if (!keyPressed) {Of.add(ToIJ(V((float)(mouseX-pmouseX),(float)(mouseY-pmouseY),0))); }
  if (keyPressed && key==CODED && keyCode==SHIFT) {Of.add(ToK(V((float)(mouseX-pmouseX),(float)(mouseY-pmouseY),0)));};
  if (keyPressed && key=='x') P.movePicked(ToIJ(V((float)(mouseX-pmouseX),(float)(mouseY-pmouseY),0))); 
  if (keyPressed && key=='z') P.movePicked(ToK(V((float)(mouseX-pmouseX),(float)(mouseY-pmouseY),0))); 
  if (keyPressed && key=='X') P.moveAll(ToIJ(V((float)(mouseX-pmouseX),(float)(mouseY-pmouseY),0))); 
  if (keyPressed && key=='Z') P.moveAll(ToK(V((float)(mouseX-pmouseX),(float)(mouseY-pmouseY),0))); 
  if (keyPressed && key=='f') { // move focus point on plane
    if(center) F.sub(ToIJ(V((float)(mouseX-pmouseX),(float)(mouseY-pmouseY),0))); 
    else F.add(ToIJ(V((float)(mouseX-pmouseX),(float)(mouseY-pmouseY),0))); 
    }
  if (keyPressed && key=='F') { // move focus point vertically
    if(center) F.sub(ToK(V((float)(mouseX-pmouseX),(float)(mouseY-pmouseY),0))); 
    else F.add(ToK(V((float)(mouseX-pmouseX),(float)(mouseY-pmouseY),0))); 
    }
  }  

//**************************** text for name, title and help  ****************************
String title ="6491 2015 P2: Steady Interpolating of Similarities in 2D", 
       name ="Student: Dingtian Zhang",
       menu="?:(show/hide) help, a: animate, `:snap picture, ~:(start/stop) recording movie frames, Q:quit",
       guide="click and drag to edit, f:showFrames"; // help info

//public void drawObject(pt P, vec V) {
//  beginShape(); 
//    texture(myFace);
//    v(P(P(P,1,V),1,R(V)),0,0);
//    v(P(P(P,1,V),-1,R(V)),1,0);
//    v(P(P(P,-1,V),-1,R(V)),1,1);
//    v(P(P(P,-1,V),1,R(V)),0,1); 
//  endShape(CLOSE);
//  }
// 
//  public void drawObjectS(pt P, vec V) {
//  beginShape(); 
//    v(P(P(P,1,V),0.25f,R(V)));
//    v(P(P(P,1,V),-0.25f,R(V)));
//    v(P(P(P,-1,V),-0.25f,R(V)));
//    v(P(P(P,-1,V),0.25f,R(V))); 
//  endShape(CLOSE);
//  }
  
public float timeWarp(float f) {return sq(sin(f*PI/2));}
/**
 * Frame class
 */
FR F0 = F(), F1 = F(), F1rF0, E0, E0rF0, FF, E, Ft, Et, Fk;
int k=7;
public FR F() {return new FR();}
public FR F(pt A, pt B, pt C) {return new FR(A,B,C);}
public FR F(vec I, vec J, vec K, pt O) {return new FR(I,J, K, O);}
public FR F(vec I, vec J, pt O) {vec K = N(I,J); return new FR(I,J, K, O);}
public FR F(FR F) {return F(F.I,F.J,F.K,F.O);}
public FR F(FR F0, float t, FR F1) {
  float a = angle(F0.I,F1.I); 
  float s = n(F1.I)/n(F0.I);
  vec N = swirlAxis(F0.I, F0.J, F0.K, F1.I, F1.J, F1.K);
  pt G = spiralCenter3D(F0, F1);
  vec unitX = projectedVec(N, F0.O, F0.I);
  unitX = unitX.normalize();
  vec unitY = N(N, unitX);
  vec I = S(pow(s,t),R(F0.I,t*a,unitX, unitY));
  pt O = P(G,pow(s,t),R(V(G,F0.O),t*a, unitX, unitY)); 
  return F(I,J,O); 
  }

public void showArrow(FR F) {F.showArrow();}
public void showArrows(FR F) {F.showArrows();}

class FR { 
  pt O; vec I; vec J; vec K; 
  FR () {O=P(); I=V(1,0,0); J=V(0,1,0); K=V(0,0,1);}
  FR(vec II, vec JJ, vec KK, pt OO) {I=V(II); J=V(JJ); K=V(KK); O=P(OO);}
  FR(pt A, pt B, pt C) {O=P(A); I=V(A,B); J=V(A,C); K=N(I,J); J=N(K,I);}
  public vec of(vec V) {return W(V.x,I,V.y,J,V.z,K);}
  public pt of(pt P) {return P(O,W(P.x,I,P.y,J,P.z,K));}
  public FR of(FR F) {return F(of(F.I),of(F.J),of(F.K),of(F.O));}
  public vec invertedOf(vec V) {return V(det2(V,I)/det2(I,J),det2(V,J)/det2(J,K),det2(V,K)/det2(K,I));}
  public pt invertedOf(pt P) {vec V = V(O,P); return P(det2(V,J)/det2(I,J),det2(V,I)/det2(J,I));}
  public FR invertedOf(FR F) {return F(invertedOf(F.I),invertedOf(F.J), invertedOf(F.K), invertedOf(F.O));}
  public FR showArrow() {show(O,4); arrow(O,I,1); return this;}
  public FR showArrows() {show(O,4); arrow(O,I,1); arrow(O,J,1); return this; }
  }
/******** Editor of an Animated Coons Patch

Implementation steps:
**<01 Manual control of (u,v) parameters. 
**<02 Draw 4 boundary curves CT(u), CB(u), SL(v), CR(v) using proportional Neville
**<03 Compute and show Coons point C(u,v)
**<04 Display quads filed one-by-one for the animated Coons patch
**<05 Compute and show normal at C(u,v) and a ball ON the patch

*/
//**<01: mouseMoved; 'v', draw: uvShow()
float u=0, v=0; 
float bu=0.5f, bv=0.5f;  // back foot
float fu=0.55f, fv=0.55f; // front foot

public void uvShow() { 
  fill(red);
  if(keyPressed && key=='v')  text("u="+u+", v="+v,10,30);
  noStroke(); fill(blue); ellipse(u*width,v*height,5,5); 
  }
/*
0 1 2 3 
11    4
10    5
9 8 7 6
*/
public pt coons(pt[] P, float s, float t) {
  pt Lst = L( L(P[0],s,P[3]), t, L(P[9],s,P[6]) ) ;
  pt Lt = L( N( 0,P[0], 1.f/3,P[1],  2.f/3,P[2],  1,P[3], s) ,t, N(0,P[9], 1.f/3,P[8], 2.f/3,P[7], 1,P[6], s) ) ;
  pt Ls = L( N( 0,P[0], 1.f/3,P[11], 2.f/3,P[10], 1,P[9], t) ,s, N(0,P[3], 1.f/3,P[4], 2.f/3,P[5] ,1,P[6], t) ) ;
  return P(Ls,V(Lst,Lt));
  }
public pt B(pt A, pt B, pt C, float s) {return L(L(A,s,B),s,L(B,s,C)); } 
public pt B(pt A, pt B, pt C, pt D, float s) {return L(B(A,B,C,s),s,B(B,C,D,s)); } 
public pt B(pt A, pt B, pt C, pt D, pt E, float s) {return L(B(A,B,C,D,s),s,B(B,C,D,E,s)); } 
public pt N(float a, pt A, float b, pt B, float t) {return L(A,(t-a)/(b-a),B);}
public pt N(float a, pt A, float b, pt B, float c, pt C, float t) {return N(a,N(a,A,b,B,t),c,N(b,B,c,C,t),t);}
public pt N(float a, pt A, float b, pt B, float c, pt C, float d, pt D, float t) {return N(a,N(a,A,b,B,c,C,t),d,N(b,B,c,C,d,D,t),t);}
public pt N(float a, pt A, float b, pt B, float c, pt C, float d, pt D, float e, pt E, float t) {return N(a,N(a,A,b,B,c,C,d,D,t),e,N(b,B,c,C,d,D,e,E,t),t);}

public void drawBorders(pt[] P){
  float e=0.01f;
  beginShape(); for(float t=0; t<1.001f; t+=e) v(coons(P,0,t)); endShape();
  beginShape(); for(float t=0; t<1.001f; t+=e) v(coons(P,1,t)); endShape();
  beginShape(); for(float t=0; t<1.001f; t+=e) v(coons(P,t,0)); endShape();
  beginShape(); for(float t=0; t<1.001f; t+=e) v(coons(P,t,1)); endShape();
  }
  
public vec CoonsNormal(pt[] P, float u, float v, float e) { 
  vec Tu = V(coons(P,u-e,v),coons(P,u+e,v));
  vec Tv = V(coons(P,u,v-e),coons(P,u,v+e));
  return U(N(Tu,Tv));
  }

public void shadeSurface(pt[] P, float e){ 
  for(float s=0; s<1.001f-e; s+=e) for(float t=0; t<1.001f-e; t+=e) 
  { beginShape(); v(coons(P,s,t)); v(coons(P,s+e,t)); v(coons(P,s+e,t+e)); v(coons(P,s,t+e)); endShape(CLOSE);}
  }

public void shadeSurfaceTextured(pt[] P, float e){ 
  fill(white);
  for(float s=0; s<1.001f-e; s+=e) for(float t=0; t<1.001f-e; t+=e) 
  {beginShape(); texture(myFace); vTextured(coons(P,s,t),s,t); vTextured(coons(P,s+e,t),s+e,t); vTextured(coons(P,s+e,t+e),s+e,t+e); vTextured(coons(P,s,t+e),s,t+e); endShape(CLOSE);}
  }

public void shadeSurfaceGouraud(pt[] P, float e, float ee){ 
  Boolean col=true;
  for(float s=0; s<1.001f-e; s+=e) {
    col=!col;
    for(float t=0; t<1.001f-e; t+=e) {
      if(col) fill(cyan); else fill(magenta); col=!col;
      beginShape(); 
      nv(CoonsNormal(P,s,t,ee)); v(coons(P,s,t)); 
      nv(CoonsNormal(P,s+e,t,ee)); v(coons(P,s+e,t)); 
      nv(CoonsNormal(P,s+e,t+e,ee)); v(coons(P,s+e,t+e)); 
      nv(CoonsNormal(P,s,t+e,ee)); v(coons(P,s,t+e)); 
      endShape(CLOSE);
      }
    }
  }

public void showNormals(pt[] P, float e, float ee){ 
  for(float s=0; s<1.001f-e; s+=e) for(float t=0; t<1.001f-e; t+=e) show(coons(P,s,t),50,CoonsNormal(P,s,t,ee));
  }


public void slide(pt[] P, float e) {
  float nfu=fu, nfv=fv;
  float z=coons(P,fu,fv).z;
  for(float a=0; a<=TWO_PI; a+=PI/20) {
    float nu=fu+e*cos(a), nv=fv+e*sin(a);
    float nz=coons(P,nu,nv).z;
    if(nz<z) {z=nz; nfu=nu; nfv=nv;}
    }
  fu=nfu; fv=nfv;
  }
  
public void attractFront(pt[] P, float e) {
  float nfu=fu, nfv=fv;
  float z=d(coons(P,fu,fv),Of);
  for(float a=0; a<=TWO_PI; a+=PI/20) {
    float nu=fu+e*cos(a), nv=fv+e*sin(a);
    float nz=d(coons(P,nu,nv),Of);
    if(nz<z) {z=nz; nfu=nu; nfv=nv;}
    }
  fu=nfu; fv=nfv;
  }  

public void attractBack(pt[] P, float r, float e) {
  float nbu=bu, nbv=bv;
  pt O = coons(PtQ.G,fu,fv);
  float z = abs(d(coons(PtQ.G,bu,bv),O)-r);
  for(float a=0; a<=TWO_PI; a+=PI/20) {
    float nu=bu+e*cos(a), nv=bv+e*sin(a);
    float nz=abs(d(coons(P,nu,nv),O)-r);
    if(nz<z) {z=nz; nbu=nu; nbv=nv;}
    }
  bu=nbu; bv=nbv;
  }  
  
public void showMan(pt[] P, float h) {
  pt Of = coons(PtQ.G,fu,fv); vec Nf = CoonsNormal(PtQ.G,fu,fv,0.01f); pt Ff = P(Of,5,Nf); pt Kf = P(Of,h,Nf); 
  pt Ob = coons(PtQ.G,bu,bv); vec Nb = CoonsNormal(PtQ.G,bu,bv,0.01f); pt Fb = P(Ob,5,Nb); pt Kb = P(Ob,h,Nb); 
  float d=d(Kf,Kb)/2,b=sqrt(sq(h)-sq(d));
  vec V = V(Nf,Nb); pt B = P(P(Kf,Kb),b,V); 
  pt T = P(B,h,V);
  fill(red); show(Ff,5); collar(Of,V(h,Nf),1,5); show(Kf,5); collar(Kf,V(Kf,B),5,10);
  fill(blue); show(Fb,5); collar(Ob,V(h,Nb),1,5); show(Kb,5); collar(Kb,V(Kb,B),5,10);
  fill(orange); show(B,10); show(T,15); collar(B,V(B,T),10,15);
  }  
  
  


pt PP=P(); // picked point
Boolean  picking=false;

public pt pick(int mX, int mY) { // returns point on visible surface at pixel (mX,My)
  PGL pgl = beginPGL();
  FloatBuffer depthBuffer = ByteBuffer.allocateDirect(1 << 2).order(ByteOrder.nativeOrder()).asFloatBuffer();
  pgl.readPixels(mX, height - mY - 1, 1, 1, PGL.DEPTH_COMPONENT, PGL.FLOAT, depthBuffer);
  float depthValue = depthBuffer.get(0);
  depthBuffer.clear();
  endPGL();

  //get 3d matrices
  PGraphics3D p3d = (PGraphics3D)g;
  PMatrix3D proj = p3d.projection.get();
  PMatrix3D modelView = p3d.modelview.get();
  PMatrix3D modelViewProjInv = proj; modelViewProjInv.apply( modelView ); modelViewProjInv.invert();
  
  float[] viewport = {0, 0, p3d.width, p3d.height};
  float[] normalized = new float[4];
  normalized[0] = ((mX - viewport[0]) / viewport[2]) * 2.0f - 1.0f;
  normalized[1] = ((height - mY - viewport[1]) / viewport[3]) * 2.0f - 1.0f;
  normalized[2] = depthValue * 2.0f - 1.0f;
  normalized[3] = 1.0f;
  
  float[] unprojected = new float[4];
  
  modelViewProjInv.mult( normalized, unprojected );
  return P( unprojected[0]/unprojected[3], unprojected[1]/unprojected[3], unprojected[2]/unprojected[3] );
  }

public pt pick(float mX, float mY, float mZ) { 
  //get 3d matrices
  PGraphics3D p3d = (PGraphics3D)g;
  PMatrix3D proj = p3d.projection.get();
  PMatrix3D modelView = p3d.modelview.get();
  PMatrix3D modelViewProjInv = proj; modelViewProjInv.apply( modelView ); modelViewProjInv.invert();
  float[] viewport = {0, 0, p3d.width, p3d.height};
  float[] normalized = new float[4];
  normalized[0] = ((mX - viewport[0]) / viewport[2]) * 2.0f - 1.0f;
  normalized[1] = ((height - mY - viewport[1]) / viewport[3]) * 2.0f - 1.0f;
  normalized[2] = mZ * 2.0f - 1.0f;
  normalized[3] = 1.0f;
  float[] unprojected = new float[4];
  modelViewProjInv.mult( normalized, unprojected );
  return P( unprojected[0]/unprojected[3], unprojected[1]/unprojected[3], unprojected[2]/unprojected[3] );
  }

public pt viewPoint() {return pick( 0,0, (height/2) / tan(PI/6));}
/*  
in draw, before popMatrix, insert
      
    if(picking) {PP = pick( mouseX, mouseY ); picking=false;} else {fill(yellow); show(PP,3);}

in keyPressed, 

      if(key=='`') picking=true; 

*/
int pp=1; // index of picked vertex
pts P = new pts(); // polyloop in 3D
pts Q = new pts(); // second polyloop in 3D
pts PtQ = new pts(); // inbetweening polyloop L(P,t,Q);
class pts { // class for manipulaitng and sisplaying polyloops
 Boolean loop=true;
 int pv =0, // picked vertex index,
     iv=0, //  insertion vertex index
     nv = 0;  // number of vertices currently used in P
 int maxnv = 16000;                 //  max number of vertices
 pt[] G = new pt [maxnv];           // geometry table (vertices)
  pts() {}
  public pts declare() {for (int i=0; i<maxnv; i++) G[i]=P(); return this;}     // init all point objects
  public pts empty() {nv=0; pv=0; return this;} // resets P so that we can start adding points
  public pts addPt(pt P) { G[nv].setTo(P); pv=nv; nv++;  return this;} // adds a point at the end
  public pts addPt(float x,float y) { G[nv].x=x; G[nv].y=y; pv=nv; nv++; return this;}
  public pts copyFrom(pts Q) {empty(); nv=Q.nv; for (int v=0; v<nv; v++) G[v]=P(Q.G[v]); return this;}
  public pts setToL(pts P, float t, pts Q) { // lerp (linear interpolation betwen P and Q
    empty(); 
    nv=min(P.nv,Q.nv); 
    for (int v=0; v<nv; v++) G[v]=L(P.G[v],t,Q.G[v]); 
    return this;}
  public pts resetOnCircle(int k, float r) { // makes new polyloo[p with k  points on a circle around origin
    empty(); // resert P
    pt C = P(); // center of circle
    for (int i=0; i<k; i++) addPt(R(P(C,V(0,-r,0)),2.f*PI*i/k,C)); // points on z=0 plane
    pv=0; // picked vertex ID is set to 0
    return this;
    } 
  public int idOfVertexWithClosestScreenProjectionTo(pt M) { // for picking a vertex with the mouse
    pp=0; 
    for (int i=1; i<nv; i++) if (d(M,ToScreen(G[i]))<=d(M,ToScreen(G[pp]))) pp=i; 
    return pp;
    }
  public pt closestProjectionOf(pt M) {   // for picking inserting O. Returns projection but also CHANGES iv !!!!
    pt C = P(G[0]); float d=d(M,C);       
    for (int i=1; i<nv; i++) if (d(M,G[i])<=d) {iv=i; C=P(G[i]); d=d(M,C); }  
    for (int i=nv-1, j=0; j<nv; i=j++) { 
       pt A = G[i], B = G[j];
       if(projectsBetween(M,A,B) && disToLine(M,A,B)<d) {d=disToLine(M,A,B); iv=i; C=projectionOnLine(M,A,B);}
       } 
    return C;    
    }
  public pts insertPt(pt P) { // inserts new vertex after vertex with ID iv
    for(int v=nv-1; v>iv; v--) G[v+1].setTo(G[v]); 
     iv++; 
     G[iv].setTo(P);
     nv++; // increments vertex count
     return this;
     }
  public pts insertClosestProjection(pt M) {  
    pt P = closestProjectionOf(M); // also sets iv
    insertPt(P);
    return this;
    }

  public pts deletePicked() {for(int i=pv; i<nv; i++) G[i].setTo(G[i+1]); pv=max(0,pv-1); nv--;  return this;}
  public pts setPt(pt P, int i) { G[i].setTo(P); return this;}
  public pts showPicked() {show(G[pv],13); return this;}
  public pts drawBalls(float r) {for (int v=0; v<nv; v++) show(G[v],r); return this;}
  public pts showPicked(float r) {show(G[pv],r); return this;}
  public pts drawClosedCurve(float r) {for (int v=0; v<nv-1; v++) stub(G[v],V(G[v],G[v+1]),r,r/2);  stub(G[nv-1],V(G[nv-1],G[0]),r,r/2); return this;}
  public pts setPickedTo(int pp) {pv=pp; return this;}
  public pts movePicked(vec V) { G[pv].add(V); return this;}      // moves selected point (index p) by amount mouse moved recently
  public pts moveAll(vec V) {for (int i=0; i<nv; i++) G[i].add(V); return this;};   
  public pt Picked() {return G[pv];} 

public void savePts(String fn) {
  String [] inppts = new String [nv+1];
  int s=0;
  inppts[s++]=str(nv);
  for (int i=0; i<nv; i++) {inppts[s++]=str(G[i].x)+","+str(G[i].y)+","+str(G[i].z);}
  saveStrings(fn,inppts);
  };
  
public void loadPts(String fn) {
  println("loading: "+fn); 
  String [] ss = loadStrings(fn);
  String subpts;
  int s=0;   int comma, comma1, comma2;   float x, y;   int a, b, c;
  nv = PApplet.parseInt(ss[s++]); print("nv="+nv);
  for(int k=0; k<nv; k++) {int i=k+s; float [] xy = PApplet.parseFloat(split(ss[i],",")); G[k].setTo(xy[0],xy[1],xy[2]);}
  pv=0;
  }; 

} // end of pts class
// points, vectors, frames in 3D
class vec { float x=0,y=0,z=0; 
   vec () {}; 
   vec (float px, float py, float pz) {x = px; y = py; z = pz;};
      vec (float px, float py) {x = px; y = py;};
   public vec set (float px, float py, float pz) {x = px; y = py; z = pz; return this;}; 
     public vec setTo(vec V) {x = V.x; y = V.y; z = V.z; return this;}; 
   public vec set (vec V) {x = V.x; y = V.y; z = V.z; return this;}; 
   public vec add(vec V) {x+=V.x; y+=V.y; z+=V.z; return this;};
   public vec add(float s, vec V) {x+=s*V.x; y+=s*V.y; z+=s*V.z; return this;};
   public vec sub(vec V) {x-=V.x; y-=V.y; z-=V.z; return this;};
   public vec mul(float f) {x*=f; y*=f; z*=f; return this;};
   public vec div(float f) {x/=f; y/=f; z/=f; return this;};
   public vec div(int f) {x/=f; y/=f; z/=f; return this;};
   public vec rev() {x=-x; y=-y; z=-z; return this;};
   public float norm() {return(sqrt(sq(x)+sq(y)+sq(z)));}; 
   public vec normalize() {float n=norm(); if (n>0.000001f) {div(n);}; return this;};
   public vec rotate(float a, vec I, vec J) { // Rotate this by angle a parallel in plane (I,J) Assumes I and J are orthogonal
     float x=d(this,I), y=d(this,J); // dot products
     float c=cos(a), s=sin(a); 
     add(x*c-x-y*s,I); add(x*s+y*c-y,J); 
     return this; }; 
   } // end class vec
  
class pt { float x=0,y=0,z=0; 
   pt () {}; 
     pt (float px, float py) {x = px; y = py;};
   pt (float px, float py, float pz) {x = px; y = py; z = pz; };
   public pt set (float px, float py, float pz) {x = px; y = py; z = pz; return this;}; 
   public pt set (pt P) {x = P.x; y = P.y; z = P.z; return this;}; 
   public pt setTo(pt P) {x = P.x; y = P.y; z = P.z; return this;}; 
   public pt setTo(float px, float py, float pz) {x = px; y = py; z = pz; return this;}; 
   public pt add(pt P) {x+=P.x; y+=P.y; z+=P.z; return this;};
   public pt add(vec V) {x+=V.x; y+=V.y; z+=V.z; return this;};
   public pt sub(vec V) {x-=V.x; y-=V.y; z-=V.z; return this;};
   public pt add(float s, vec V) {x+=s*V.x; y+=s*V.y; z+=s*V.z; return this;};
   public pt sub(pt P) {x-=P.x; y-=P.y; z-=P.z; return this;};
   public pt mul(float f) {x*=f; y*=f; z*=f; return this;};
   public pt div(float f) {x/=f; y/=f; z/=f; return this;};
   public pt div(int f) {x/=f; y/=f; z/=f; return this;};
   }
   
// =====  vector functions
public vec V() {return new vec(); };                                                                          // make vector (x,y,z)
public vec V(float x, float y, float z) {return new vec(x,y,z); };                                            // make vector (x,y,z)
public vec V(vec V) {return new vec(V.x,V.y,V.z); };                                                          // make copy of vector V
public vec A(vec A, vec B) {return new vec(A.x+B.x,A.y+B.y,A.z+B.z); };                                       // A+B
public vec A(vec U, float s, vec V) {return V(U.x+s*V.x,U.y+s*V.y,U.z+s*V.z);};                               // U+sV
public vec M(vec U, vec V) {return V(U.x-V.x,U.y-V.y,U.z-V.z);};                                              // U-V
public vec M(vec V) {return V(-V.x,-V.y,-V.z);};                                                              // -V
public vec V(vec A, vec B) {return new vec((A.x+B.x)/2.0f,(A.y+B.y)/2.0f,(A.z+B.z)/2.0f); }                      // (A+B)/2
public vec V(vec A, float s, vec B) {return new vec(A.x+s*(B.x-A.x),A.y+s*(B.y-A.y),A.z+s*(B.z-A.z)); };      // (1-s)A+sB
public vec V(vec A, vec B, vec C) {return new vec((A.x+B.x+C.x)/3.0f,(A.y+B.y+C.y)/3.0f,(A.z+B.z+C.z)/3.0f); };  // (A+B+C)/3
public vec V(vec A, vec B, vec C, vec D) {return V(V(A,B),V(C,D)); };                                         // (A+B+C+D)/4
public vec V(float s, vec A) {return new vec(s*A.x,s*A.y,s*A.z); };                                           // sA
public vec V(float a, vec A, float b, vec B) {return A(V(a,A),V(b,B));}                                       // aA+bB 
public vec V(float a, vec A, float b, vec B, float c, vec C) {return A(V(a,A,b,B),V(c,C));}                   // aA+bB+cC
public vec V(pt P, pt Q) {return new vec(Q.x-P.x,Q.y-P.y,Q.z-P.z);};                                          // PQ
public vec U(vec V) {float n = V.norm(); if (n<0.0000001f) return V(0,0,0); else return V.div(n);};             // V/||V||
public vec U(pt P, pt Q) {return U(V(P,Q));};                                                                 // PQ/||PQ||
public vec U(float x, float y, float z) {return U(V(x,y,z)); };                                               // make vector (x,y,z)
public vec N(vec U, vec V) {return V( U.y*V.z-U.z*V.y, U.z*V.x-U.x*V.z, U.x*V.y-U.y*V.x); };                  // UxV cross product (normal to both)
public vec N(pt A, pt B, pt C) {return N(V(A,B),V(A,C)); };                                                   // normal to triangle (A,B,C), not normalized (proportional to area)
public vec B(vec U, vec V) {return U(N(N(U,V),U)); }        
public vec Normal(vec V) {
  if(abs(V.z)<=min(abs(V.x),abs(V.y))) return V(-V.y,V.x,0); 
  if(abs(V.x)<=min(abs(V.z),abs(V.y))) return V(0,-V.z,V.y);
  return V(V.z,0,-V.x);
  }
//************************************************************************
//**** VECTOR FUNCTIONS
//************************************************************************
//create 
public vec V(pt P) {return new vec(P.x,P.y,P.z); };                                                              // make vector from origin to P

//weighted sum 
public vec W(float s,vec V) {return V(s*V.x,s*V.y,s*V.z);}                                                      // sV
public vec W(vec U, vec V) {return V(U.x+V.x,U.y+V.y,U.z+V.z);}                                                   // U+V 
public vec W(vec U, vec V, vec W) {return V(U.x+V.x+W.x,U.y+V.y+W.y,U.z+V.z+W.z);}                                                   // U+V 
public vec W(vec U,float s,vec V) {return W(U,S(s,V));}                                                   // U+sV
public vec W(float u, vec U, float v, vec V) {return W(S(u,U),S(v,V));}                                   // uU+vV ( Linear combination)
public vec W(float u, vec U, float v, vec V, float w, vec W) {return W(S(u,U),S(v,V),S(w,W));}                                   // uU+vV ( Linear combination)

//transformed 
//public vec R(vec V, vec N, float a) {R(pt P, float a, vec I, vec J, pt G)};          // V rotated by a radians
public vec S(float s,vec V) {return new vec(s*V.x,s*V.y,s*V.z);};                                                  // sV
public vec Reflection(vec V, vec N) { return W(V,-2.f*dot(V,N),N);};                                          // reflection                                                                 // -V

//Interpolation 
public vec L(vec U, vec V, float s) {return new vec(U.x+s*(V.x-U.x),U.y+s*(V.y-U.y),U.z+s*(V.z-U.z));};                      // (1-s)U+sV (Linear interpolation between vectors)
//public vec S(vec U, vec V, float s) {float a = angle(U,V); vec W = R(U,s*a); float u = n(U), v=n(V); return W(pow(v/u,s),W); } // steady interpolation from U to V

//measure 
public float angle(pt A, pt B, pt C) {return  angle(V(B,A),V(B,C)); }                                       // angle <BA,BC>
public float turnAngle(pt A, pt B, pt C) {return  angle(V(A,B),V(B,C)); }                                   // angle <AB,BC> (positive when right turn as seen on screen)
public int toDeg(float a) {return PApplet.parseInt(a*180/PI);}                                                           // convert radians to degrees
public float toRad(float a) {return(a*PI/180);}                                                             // convert degrees to radians 
public float positive(float a) { if(a<0) return a+TWO_PI; else return a;}                                   // adds 2PI to make angle positive

//SLERP
public vec slerp(vec U, float t, vec V) {float a = angle(U,V); float b=sin((1.f-t)*a),c=sin(t*a),d=sin(a); return W(b/d,U,c/d,V); } // UNIT vectors ONLY!

// ===== point functions
public pt P() {return new pt(); };                                                                          // point (x,y,z)
public pt P(float x, float y, float z) {return new pt(x,y,z); };                                            // point (x,y,z)
public pt P(float x, float y) {return new pt(x,y); };                                                       // make point (x,y)
public pt P(pt A) {return new pt(A.x,A.y,A.z); };                                                           // copy of point P
public pt P(pt A, float s, pt B) {return new pt(A.x+s*(B.x-A.x),A.y+s*(B.y-A.y),A.z+s*(B.z-A.z)); };        // A+sAB
public pt L(pt A, float s, pt B) {return new pt(A.x+s*(B.x-A.x),A.y+s*(B.y-A.y),A.z+s*(B.z-A.z)); };        // A+sAB
public pt P(pt A, pt B) {return P((A.x+B.x)/2.0f,(A.y+B.y)/2.0f,(A.z+B.z)/2.0f); }                             // (A+B)/2
public pt P(pt A, pt B, pt C) {return new pt((A.x+B.x+C.x)/3.0f,(A.y+B.y+C.y)/3.0f,(A.z+B.z+C.z)/3.0f); };     // (A+B+C)/3
public pt P(pt A, pt B, pt C, pt D) {return P(P(A,B),P(C,D)); };                                            // (A+B+C+D)/4
public pt P(float s, pt A) {return new pt(s*A.x,s*A.y,s*A.z); };                                            // sA
public pt A(pt A, pt B) {return new pt(A.x+B.x,A.y+B.y,A.z+B.z); };                                         // A+B
public pt P(float a, pt A, float b, pt B) {return A(P(a,A),P(b,B));}                                        // aA+bB 
public pt P(float a, pt A, float b, pt B, float c, pt C) {return A(P(a,A),P(b,B,c,C));}                     // aA+bB+cC 
public pt P(float a, pt A, float b, pt B, float c, pt C, float d, pt D){return A(P(a,A,b,B),P(c,C,d,D));}   // aA+bB+cC+dD
public pt P(pt P, vec V) {return new pt(P.x + V.x, P.y + V.y, P.z + V.z); }                                 // P+V
public pt P(pt P, float s, vec V) {return new pt(P.x+s*V.x,P.y+s*V.y,P.z+s*V.z);}                           // P+sV
public pt P(pt O, float x, vec I, float y, vec J) {return P(O.x+x*I.x+y*J.x,O.y+x*I.y+y*J.y,O.z+x*I.z+y*J.z);}  // O+xI+yJ
public pt P(pt O, float x, vec I, float y, vec J, float z, vec K) {return P(O.x+x*I.x+y*J.x+z*K.x,O.y+x*I.y+y*J.y+z*K.y,O.z+x*I.z+y*J.z+z*K.z);}  // O+xI+yJ+kZ
public void makePts(pt[] C) {for(int i=0; i<C.length; i++) C[i]=P();}
public pt ToScreen(pt P) {return P(screenX(P.x,P.y,P.z),screenY(P.x,P.y,P.z),0);}  // O+xI+yJ+kZ
public pt ToModel(pt P) {return P(modelX(P.x,P.y,P.z),modelY(P.x,P.y,P.z),modelZ(P.x,P.y,P.z));}  // O+xI+yJ+kZ

// ===== mouse
public pt Mouse() {return P(mouseX,mouseY,0);};                                          // current mouse location
public pt Pmouse() {return P(pmouseX,pmouseY,0);};
public vec MouseDrag() {return V(mouseX-pmouseX,mouseY-pmouseY,0);};                     // vector representing recent mouse displacement
public pt ScreenCenter() {return P(width/2,height/2);}                                                        //  point in center of  canvas

// ===== measures
public float d(vec U, vec V) {return U.x*V.x+U.y*V.y+U.z*V.z; };                                            //U*V dot product
public float dot(vec U, vec V) {return U.x*V.x+U.y*V.y+U.z*V.z; };                                            //U*V dot product
public float det2(vec U, vec V) {return -U.y*V.x+U.x*V.y; };                                       // U|V det product
public float det3(vec U, vec V) {return sqrt(d(U,U)*d(V,V) - sq(d(U,V))); };                                       // U|V det product
public float m(vec U, vec V, vec W) {return d(U,N(V,W)); };                                                 // (UxV)*W  mixed product, determinant
public float m(pt E, pt A, pt B, pt C) {return m(V(E,A),V(E,B),V(E,C));}                                    // det (EA EB EC) is >0 when E sees (A,B,C) clockwise
public float n2(vec V) {return sq(V.x)+sq(V.y)+sq(V.z);};                                                   // V*V    norm squared
public float n(vec V) {return sqrt(n2(V));};                                                                // ||V||  norm
public float d(pt P, pt Q) {return sqrt(sq(Q.x-P.x)+sq(Q.y-P.y)+sq(Q.z-P.z)); };                            // ||AB|| distance
public float area(pt A, pt B, pt C) {return n(N(A,B,C))/2; };                                               // area of triangle 
public float volume(pt A, pt B, pt C, pt D) {return m(V(A,B),V(A,C),V(A,D))/6; };                           // volume of tet 
public boolean parallel (vec U, vec V) {return n(N(U,V))<n(U)*n(V)*0.00001f; }                              // true if U and V are almost parallel
public float angle(vec U, vec V) {return acos(d(U,V)/n(V)/n(U)); };                                       // angle(U,V)
public boolean cw(vec U, vec V, vec W) {return m(U,V,W)>0; };                                               // (UxV)*W>0  U,V,W are clockwise
public boolean cw(pt A, pt B, pt C, pt D) {return volume(A,B,C,D)>0; };                                     // tet is oriented so that A sees B, C, D clockwise 
public boolean projectsBetween(pt P, pt A, pt B) {return dot(V(A,P),V(A,B))>0 && dot(V(B,P),V(B,A))>0 ; };
public float disToLine(pt P, pt A, pt B) {return det3(U(A,B),V(A,P)); };
public pt projectionOnLine(pt P, pt A, pt B) {return P(A,dot(V(A,B),V(A,P))/dot(V(A,B),V(A,B)),V(A,B));}

// ===== rotate 
public vec R(vec V) {return V(-V.y,V.x,V.z);} // rotated 90 degrees in XY plane
public pt R(pt P, float a, vec I, vec J, pt G) {float x=d(V(G,P),I), y=d(V(G,P),J); float c=cos(a), s=sin(a); return P(P,x*c-x-y*s,I,x*s+y*c-y,J); }; // Rotated P by a around G in plane (I,J)
public vec R(vec V, float a, vec I, vec J) {float x=d(V,I), y=d(V,J); float c=cos(a), s=sin(a); return A(V,V(x*c-x-y*s,I,x*s+y*c-y,J)); }; // Rotated V by a parallel to plane (I,J)
public pt R(pt Q, pt C, pt P, pt R) { // returns rotated version of Q by angle(CP,CR) parallel to plane (C,P,R)
   vec I0=U(C,P), I1=U(C,R), V=V(C,Q); 
   float c=d(I0,I1), s=sqrt(1.f-sq(c)); 
     if(abs(s)<0.00001f) return Q;
   vec J0=V(1.f/s,I1,-c/s,I0);  
   vec J1=V(-s,I0,c,J0);  
   float x=d(V,I0), y=d(V,J0);
                                //  stroke(red); show(C,400,I0); stroke(blue); show(C,400,I1); stroke(orange); show(C,400,J0); stroke(magenta); show(C,400,J1); noStroke();
   return P(Q,x,M(I1,I0),y,M(J1,J0)); 
  } 
public pt R(pt Q, float a) {float dx=Q.x, dy=Q.y, c=cos(a), s=sin(a); return P(c*dx+s*dy,-s*dx+c*dy,Q.z); };  // Q rotated by angle a around the origin
public pt R(pt Q, float a, pt C) {float dx=Q.x-C.x, dy=Q.y-C.y, c=cos(a), s=sin(a); return P(C.x+c*dx-s*dy, C.y+s*dx+c*dy, Q.z); };  // Q rotated by angle a around point P

// ===== render
public void normal(vec V) {normal(V.x,V.y,V.z);};                                          // changes normal for smooth shading
public void vertex(pt P) {vertex(P.x,P.y,P.z);};                                           // vertex for shading or drawing
public void v(pt P) {vertex(P.x,P.y,P.z);};                                           // vertex for shading or drawing
public void nv(vec N) {normal(N.x,N.y,N.z);};                                           // vertex for shading or drawing
public void vTextured(pt P, float u, float v) {vertex(P.x,P.y,P.z,u,v);};                          // vertex with texture coordinates
public void show(pt P, pt Q) {line(Q.x,Q.y,Q.z,P.x,P.y,P.z); };                       // draws edge (P,Q)
public void show(pt P, vec V) {line(P.x,P.y,P.z,P.x+V.x,P.y+V.y,P.z+V.z); };          // shows edge from P to P+V
public void show(pt P, float d , vec V) {line(P.x,P.y,P.z,P.x+d*V.x,P.y+d*V.y,P.z+d*V.z); }; // shows edge from P to P+dV
public void show(pt A, pt B, pt C) {beginShape(); vertex(A);vertex(B); vertex(C); endShape(CLOSE);};                      // volume of tet 
public void show(pt A, pt B, pt C, pt D) {beginShape(); vertex(A); vertex(B); vertex(C); vertex(D); endShape(CLOSE);};                      // volume of tet 
public void show(pt P, float r) {pushMatrix(); translate(P.x,P.y,P.z); sphere(r); popMatrix();}; // render sphere of radius r and center P
public void show(pt P, float s, vec I, vec J, vec K) {noStroke(); fill(yellow); show(P,5); stroke(red); show(P,s,I); stroke(green); show(P,s,J); stroke(blue); show(P,s,K); }; // render sphere of radius r and center P
public void show(pt P, String s) {text(s, P.x, P.y, P.z); }; // prints string s in 3D at P
public void show(pt P, String s, vec D) {text(s, P.x+D.x, P.y+D.y, P.z+D.z);  }; // prints string s in 3D at P+D
public void showShadow(pt P, float r) {pushMatrix(); translate(P.x,P.y,0); scale(1,1,0.01f); sphere(r); popMatrix();}

public String toText(vec V){ return "("+nf(V.x,1,5)+","+nf(V.y,1,5)+","+nf(V.z,1,5)+")";}
// ==== curve
public void bezier(pt A, pt B, pt C, pt D) {bezier(A.x,A.y,A.z,B.x,B.y,B.z,C.x,C.y,C.z,D.x,D.y,D.z);} // draws a cubic Bezier curve with control points A, B, C, D
public void bezier(pt [] C) {bezier(C[0],C[1],C[2],C[3]);} // draws a cubic Bezier curve with control points A, B, C, D
public pt bezierPoint(pt[] C, float t) {return P(bezierPoint(C[0].x,C[1].x,C[2].x,C[3].x,t),bezierPoint(C[0].y,C[1].y,C[2].y,C[3].y,t),bezierPoint(C[0].z,C[1].z,C[2].z,C[3].z,t)); }
public vec bezierTangent(pt[] C, float t) {return V(bezierTangent(C[0].x,C[1].x,C[2].x,C[3].x,t),bezierTangent(C[0].y,C[1].y,C[2].y,C[3].y,t),bezierTangent(C[0].z,C[1].z,C[2].z,C[3].z,t)); }
public void PT(pt P0, vec T0, pt P1, vec T1) {float d=d(P0,P1)/3;  bezier(P0, P(P0,-d,U(T0)), P(P1,-d,U(T1)), P1);} // draws cubic Bezier interpolating  (P0,T0) and  (P1,T1) 
public void PTtoBezier(pt P0, vec T0, pt P1, vec T1, pt [] C) {float d=d(P0,P1)/3;  C[0].set(P0); C[1].set(P(P0,-d,U(T0))); C[2].set(P(P1,-d,U(T1))); C[3].set(P1);} // draws cubic Bezier interpolating  (P0,T0) and  (P1,T1) 
public vec vecToCubic (pt A, pt B, pt C, pt D, pt E) {return V( (-A.x+4*B.x-6*C.x+4*D.x-E.x)/6, (-A.y+4*B.y-6*C.y+4*D.y-E.y)/6, (-A.z+4*B.z-6*C.z+4*D.z-E.z)/6);}
public vec vecToProp (pt B, pt C, pt D) {float cb=d(C,B);  float cd=d(C,D); return V(C,P(B,cb/(cb+cd),D)); };  

// ==== perspective
public pt Pers(pt P, float d) { return P(d*P.x/(d+P.z) , d*P.y/(d+P.z) , d*P.z/(d+P.z) ); };
public pt InverserPers(pt P, float d) { return P(d*P.x/(d-P.z) , d*P.y/(d-P.z) , d*P.z/(d-P.z) ); };

// ==== intersection
public boolean intersect(pt P, pt Q, pt A, pt B, pt C, pt X)  {return intersect(P,V(P,Q),A,B,C,X); } // if (P,Q) intersects (A,B,C), return true and set X to the intersection point
public boolean intersect(pt E, vec T, pt A, pt B, pt C, pt X) { // if ray from E along T intersects triangle (A,B,C), return true and set X to the intersection point
  vec EA=V(E,A), EB=V(E,B), EC=V(E,C), AB=V(A,B), AC=V(A,C); 
  boolean s=cw(EA,EB,EC), sA=cw(T,EB,EC), sB=cw(EA,T,EC), sC=cw(EA,EB,T); 
  if ( (s==sA) && (s==sB) && (s==sC) ) return false;
  float t = m(EA,AC,AB) / m(T,AC,AB);
  X.set(P(E,t,T));
  return true;
  }
public boolean rayIntersectsTriangle(pt E, vec T, pt A, pt B, pt C) { // true if ray from E with direction T hits triangle (A,B,C)
  vec EA=V(E,A), EB=V(E,B), EC=V(E,C); 
  boolean s=cw(EA,EB,EC), sA=cw(T,EB,EC), sB=cw(EA,T,EC), sC=cw(EA,EB,T); 
  return  (s==sA) && (s==sB) && (s==sC) ;};
public boolean edgeIntersectsTriangle(pt P, pt Q, pt A, pt B, pt C)  {
  vec PA=V(P,A), PQ=V(P,Q), PB=V(P,B), PC=V(P,C), QA=V(Q,A), QB=V(Q,B), QC=V(Q,C); 
  boolean p=cw(PA,PB,PC), q=cw(QA,QB,QC), a=cw(PQ,PB,PC), b=cw(PA,PQ,PC), c=cw(PQ,PB,PQ); 
  return (p!=q) && (p==a) && (p==b) && (p==c);
  }
public float rayParameterToIntersection(pt E, vec T, pt A, pt B, pt C) {vec AE=V(A,E), AB=V(A,B), AC=V(A,C); return - m(AE,AC,AB) / m(T,AC,AB);}
   
public float angleDraggedAround(pt G) {  // returns angle in 2D dragged by the mouse around the screen projection of G
   pt S=P(screenX(G.x,G.y,G.z),screenY(G.x,G.y,G.z),0);
   vec T=V(S,Pmouse()); vec U=V(S,Mouse());
   return atan2(d(R(U),T),d(U,T));
   }
 
public float scaleDraggedFrom(pt G) {pt S=P(screenX(G.x,G.y,G.z),screenY(G.x,G.y,G.z),0); return d(S,Mouse())/d(S,Pmouse()); }

// FANS, CONES, AND ARROWS
public void disk(pt P, vec V, float r) {  
  vec I = U(Normal(V));
  vec J = U(N(I,V));
  disk(P,I,J,r);
  }

public void disk(pt P, vec I, vec J, float r) {
  float da = TWO_PI/36;
  beginShape(TRIANGLE_FAN);
    v(P);
    for(float a=0; a<=TWO_PI+da; a+=da) v(P(P,r*cos(a),I,r*sin(a),J));
  endShape();
  }
  

public void fan(pt P, vec V, float r) {  
  vec I = U(Normal(V));
  vec J = U(N(I,V));
  fan(P,V,I,J,r);
  }

public void fan(pt P, vec V, vec I, vec J, float r) {
  float da = TWO_PI/36;
  beginShape(TRIANGLE_FAN);
    v(P(P,V));
    for(float a=0; a<=TWO_PI+da; a+=da) v(P(P,r*cos(a),I,r*sin(a),J));
  endShape();
  }
  
public void collar(pt P, vec V, float r, float rd) {
  vec I = U(Normal(V));
  vec J = U(N(I,V));
  collar(P,V,I,J,r,rd);
  }
 
public void collar(pt P, vec V, vec I, vec J, float r, float rd) {
  float da = TWO_PI/36;
  beginShape(QUAD_STRIP);
    for(float a=0; a<=TWO_PI+da; a+=da) {v(P(P,r*cos(a),I,r*sin(a),J,0,V)); v(P(P,rd*cos(a),I,rd*sin(a),J,1,V));}
  endShape();
  }

public void cone(pt P, vec V, float r) {fan(P,V,r); disk(P,V,r);}

public void stub(pt P, vec V, float r, float rd) {
  collar(P,V,r,rd); disk(P,V,r); disk(P(P,V),V,rd); 
  }
  
public void arrow(pt P, vec V, float r) {
  stub(P,V(.8f,V),r*2/3,r/3); 
  cone(P(P,V(.8f,V)),V(.2f,V),r); 
  }  

// **************************** PRIMITIVE
public void showFrame(float d) { 
  noStroke(); 
  fill(metal); sphere(d/10);
  fill(blue);  showArrow(d,d/10);
  fill(red); pushMatrix(); rotateY(PI/2); showArrow(d,d/10); popMatrix();
  fill(green); pushMatrix(); rotateX(-PI/2); showArrow(d,d/10); popMatrix();
  }

public void showFan(float d, float r) {
  float da = TWO_PI/36;
  beginShape(TRIANGLE_FAN);
    vertex(0,0,d);
    for(float a=0; a<=TWO_PI+da; a+=da) vertex(r*cos(a),r*sin(a),0);
  endShape();
  }

public void showCollar(float d, float r, float rd) {
  float da = TWO_PI/36;
  beginShape(QUAD_STRIP);
    for(float a=0; a<=TWO_PI+da; a+=da) {vertex(r*cos(a),r*sin(a),0); vertex(rd*cos(a),rd*sin(a),d);}
  endShape();
  }

public void showCone(float d, float r) {showFan(d,r);  showFan(0,r);}

public void showStub(float d, float r, float rd) {
  showCollar(d,r,rd); showFan(0,r);  pushMatrix(); translate(0,0,d); showFan(0,rd); popMatrix();
  }

public void showArrow() {showArrow(1,0.08f);}
 
public void showArrow(float d, float r) {
  float dd=d/5;
  showStub(d-dd,r*2/3,r/3); pushMatrix(); translate(0,0,d-dd); showCone(dd,r); popMatrix();
  }  
  
public void showBlock(float w, float d, float h, float x, float y, float z, float a) {
  pushMatrix(); translate(x,y,h/2); rotateZ(TWO_PI*a); box(w, d, h); popMatrix(); 
  }

//*********** PICK
vec I=V(1,0,0), J=V(0,1,0), K=V(0,0,1); // screen projetions of global model frame

public void computeProjectedVectors() { 
  pt O = ToScreen(P(0,0,0));
  pt A = ToScreen(P(1,0,0));
  pt B = ToScreen(P(0,1,0));
  pt C = ToScreen(P(0,0,1));
  I=V(O,A);
  J=V(O,B);
  K=V(O,C);
  }

public vec ToIJ(vec V) {
 float x = det2(V,J) / det2(I,J);
 float y = det2(V,I) / det2(J,I);
 return V(x,y,0);
 }
 
public vec ToK(vec V) {
 float z = dot(V,K) / dot(K,K);
 return V(0,0,z);
 }
// ************************************ IMAGES & VIDEO 
int pictureCounter=0, frameCounter=0;
Boolean filming=false, change=false;
PImage myFace; // picture of author's face, should be: data/pic.jpg in sketch folder
public void snapPicture() {saveFrame("PICTURES/P"+nf(pictureCounter++,3)+".jpg"); }

// ******************************************COLORS 
int black=0xff000000, white=0xffFFFFFF, // set more colors using Menu >  Tools > Color Selector
   red=0xffFF0000, green=0xff00FF01, blue=0xff0300FF, yellow=0xffFEFF00, cyan=0xff00FDFF, magenta=0xffFF00FB,
   grey=0xff818181, orange=0xffFFA600, brown=0xffB46005, metal=0xffB5CCDE, dgreen=0xff157901;

//************************************************************************ GRAPHICS 
public void pen(int c, float w) {stroke(c); strokeWeight(w); noFill();}
public void pen(int c, float w, int f) {stroke(c); strokeWeight(w); fill(f);}
public void penFill(int c, float w) {stroke(c); strokeWeight(w); fill(c);}
public void showDisk(float x, float y, float r) {ellipse(x,y,r*2,r*2);}

// ******************************** TEXT , TITLE, and USER's GUIDE
Boolean scribeText=true; // toggle for displaying of help text
public void scribe(String S, float x, float y) {fill(0); text(S,x,y); noFill();} // writes on screen at (x,y) with current fill color
public void scribeHeader(String S, int i) {fill(0); text(S,10,20+i*20); noFill();} // writes black at line i
public void scribeHeaderRight(String S) {fill(0); text(S,width-7.5f*S.length(),20); noFill();} // writes black on screen top, right-aligned
public void scribeFooter(String S, int i) {fill(0); text(S,10,height-10-i*20); noFill();} // writes black on screen at line i from bottom
public void scribeAtMouse(String S) {fill(0); text(S,mouseX,mouseY); noFill();} // writes on screen near mouse
public void scribeMouseCoordinates() {fill(black); text("("+mouseX+","+mouseY+")",mouseX+7,mouseY+25); noFill();}

// **************************** FILE SELECTION FOR SAVING AND LOADING MODELS 
//String fileName="data/points";

//String path="data/pts"; 
//void saveToFile(File selection) {
//  if (selection == null) println("Window was closed or the user hit cancel.");
//  else path=selection.getAbsolutePath();
//  println("    save path = "+path);
//  }
//
//void readFromFile(File selection) {
//  if (selection == null) println("Window was closed or the user hit cancel or file not found.");
//  else path=selection.getAbsolutePath();
//  println("    read path = "+path);
//  }
//
//
//void fileSelected(File selection) {
//  if (selection == null) println("Window was closed or the user hit cancel.");
//  else {
//    fileName = selection.getAbsolutePath();
//    println("User selected " + fileName);
//    }
//  }
//
                                 // writes string S at P+V
//************************************************************************
//**** Swirl
//************************************************************************
public pt PtOnSpiral(pt A, pt B, pt C, float t) {
  float a =spiralAngle(A,B,B,C); 
  float s =spiralScale(A,B,B,C);
  pt G = spiralCenter(a, s, A, B); 
  return L(G,pow(s,t),R(B,t*a,G));
  }

public pt spiralPt(pt A, pt G, float s, float a) {return L(G,s,R(A,a,G));}  
public pt spiralPt(pt A, pt G, float s, float a, float t) {return L(G,pow(s,t),R(A,t*a,G));} 
public vec swirlAxis(vec I_0, vec J_0, vec K_0, vec I_1, vec J_1, vec K_1) {
	vec delta_I = I_1.sub(I_0);
	vec delta_J = J_1.sub(J_0);
	vec delta_K = K_1.sub(K_0);
	vec ret = new vec(0,0,0);
	ret.add(N(delta_I, delta_J));
	ret.add(N(delta_J, delta_K));
	ret.add(N(delta_K, delta_I));
	return ret.div(3).normalize();
}

public vec projectedVec(vec N, pt O, vec vector) {
	vec v = W(vector, M(W(dot(N, vector), N)));
	return v;
}


public pt projectedCoord(vec N, pt O, vec unitX, vec OP) {
	vec v = projectedVec(N, O, OP);
	float x = dot(unitX, v);
	float y = det2(unitX, v);
	return new pt(x, y);
}

public pt spiralCenter3D(FR F0, FR F1) { 
	vec N = swirlAxis(F0.I, F0.J, F0.K, F1.I, F1.J, F1.K);
	vec unitX = projectedVec(N, F0.O, F0.I);
	unitX = unitX.normalize();
	pt O_0 = new pt(0,0);
	pt P_0 = projectedCoord(N, F0.O, unitX, F0.I);
	pt O_1 = projectedCoord(N, F0.O, unitX, V(F0.O, F1.O));
	pt P_1 = projectedCoord(N, F0.O, unitX, W(V(F0.O, F1.O), F0.I));
	return spiralCenter(O_0, P_0, O_1, P_1);
  }

public pt spiralCenter(pt A, pt B, pt C, pt D) { // computes center of spiral that takes A to C and B to D
	float a = spiralAngle(A,B,C,D); 
	  float z = spiralScale(A,B,C,D);
	  return spiralCenter(a,z,A,C);
	
  }
public float spiralAngle(pt A, pt B, pt C, pt D) {return angle(V(A,B),V(C,D));}
public float spiralScale(pt A, pt B, pt C, pt D) {return d(C,D)/d(A,B);}
public pt spiralCenter(float a, float z, pt A, pt C) {
  float c=cos(a), s=sin(a);
  float D = sq(c*z-1)+sq(s*z);
  float ex = c*z*A.x - C.x - s*z*A.y;
  float ey = c*z*A.y - C.y + s*z*A.x;
  float x=(ex*(c*z-1) + ey*s*z) / D;
  float y=(ey*(c*z-1) - ex*s*z) / D;
  return P(x,y);
  }
  
//************************************************************************ SAVING INDIVIDUAL IMAGES OF CANVAS 
boolean snapPic=false;
String PicturesOutputPath="data/PDFimages";

public void displayHeader() { // Displays title and authors face on screen
    scribeHeader(title,0); scribeHeaderRight(name); 
    image(myFace, width-myFace.width/2,25,myFace.width/2,myFace.height/2); 
    }
public void displayFooter() { // Displays help text at the bottom
    scribeFooter(guide,1); 
    scribeFooter(menu,0); 
    }
  static public void main(String[] passedArgs) {
    String[] appletArgs = new String[] { "Swirl" };
    if (passedArgs != null) {
      PApplet.main(concat(appletArgs, passedArgs));
    } else {
      PApplet.main(appletArgs);
    }
  }
}
