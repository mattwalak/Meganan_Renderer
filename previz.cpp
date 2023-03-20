#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <math.h>
#include <random>
#include <iostream>
#include <fstream>
#include <float.h>
#include "SETTINGS.h"
#include "skeleton.h"
#include "displaySkeleton.h"
#include "motion.h"
#include "parse_stl.h"
#include "MERSENNE_TWISTER.h"
#include "util.cpp"
#include <sys/time.h>

using namespace std;
using namespace stl;

float PI = 3.1415926535;
int PHONG = 10;
MERSENNE_TWISTER m;

// DEPTH OF FIELD
float aperature = 0.05; // Length of one side of aperature square <--- ANIMATED change in keyframe data plz
int dof_N = 50;
bool dof_on = true;

// MOTION BLUR
bool motionBlur = true;
float motion_N = 10;

// DUST
bool dust = true;
int dust_N = 100;
float dustConst;

bool drawPOI = false;
float downsize = 1; // Percent of resultion to render at

VEC3 ambient;

vector<float> dustKeys = {0,68,83,198,213};
vector<float> dustVals = {1,1,1, 1, 1};
vector<int> dustInterps = {1,1,1};

// Time keyframes for stick man
vector<float> stickKeyframes = {0, 83, 177, 198, 299};
vector<float> stickDataIndexValues = {0, 220, 362, 432, 598};
vector<int> stickInterpolation = {0,0,0,0,0};

// Keyframes for ambient light
vector<float> ambientKeys = {0,68, 83};
vector<VEC3> ambientValues = {VEC3(.01, .01, .01), VEC3(.01, .01, .01), VEC3(0, 0, 0)};
vector<int> ambientInterps = {1,1,1};

// Keyframes for distributed effects
vector<float> dofNKeys = {0, 83, 84};
vector<float> dofNVals = {10, 10, 10};
vector<int> dofNInterps = {1,1,1};

vector<float> motionNKeys = {0, 33,35,50,52, 66,68};
vector<float> motionNVals = {5, 5,5,5,5, 5,5};
vector<int> motionNInterps = {1, 1,1,1,1, 1,1};

// Cam keyframes
vector<float> camKeyFrames = {0, 20, 40, 68, 83, 181, 198, 213};
vector<VEC3> camEyeValues = {VEC3(0, .35, 1.7), VEC3(0, .35, 1.7), VEC3(0, .35, 1.7), VEC3(0, .35, 1.7), VEC3(-4, 5, 10), VEC3(-4, 1.5, 1), VEC3(-4, 1.5, 1), VEC3(-1.3, 2, .75)};
vector<VEC3> camLookValues = {VEC3(1.35, .275, -.55), VEC3(1.35, .275, -.55), VEC3(.3,.333,1.2), VEC3(.3,.333,1.2), VEC3(.4,1,.2), VEC3(.4, .6, 1.2), VEC3(.4, .6, 1.2), VEC3(.6, .7, .9)};
vector<VEC3> camUpValues = {VEC3(0,1,0), VEC3(0,1,0), VEC3(0,1,0), VEC3(0,1,0), VEC3(0,1,0), VEC3(0,1,0), VEC3(0,1,0), VEC3(0,1,0)};
vector<float> aperatureValues = {0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03};
vector<int> camInterpolation = {1,1,1, 1, 1, 1, 1, 1};

// scene geometry
vector<VEC3> sphereCenters;
vector<float> sphereRadii;
vector<VEC3> sphereColors;

vector<VEC3> triangleVert1;
vector<VEC3> triangleVert2;
vector<VEC3> triangleVert3;
vector<VEC3> triangleColors;

vector<VEC3> cylinderV0;
vector<VEC3> cylinderV1;
vector<float> cylinderRadii;
vector<VEC3> cylinderColors;

#include "primitives.cpp"

vector<Primitive *> primitives;
vector<BoundingBox *> bounders;
vector<Light *> lights;

#include "model.cpp"

vector<Model> models;

// Stick-man classes
DisplaySkeleton displayer;    
Skeleton* skeleton;
Motion* motion;

int windowWidth = downsize*1920;
int windowHeight = downsize*1080;

//VEC3 eye(-6, 0.5, 1); ORIGINAL
//VEC3 lookingAt(5, 0.5, 1); ORIGINAL
VEC3 eye(-3, 9, 6);
VEC3 lookingAt(.4, 1, 1.2);
VEC3 up(0,1,0);

void clearGeometry(){
	sphereCenters.clear();
	sphereRadii.clear();
	sphereColors.clear();
	triangleVert1.clear();
	triangleVert2.clear();
	triangleVert3.clear();
	triangleColors.clear();
	cylinderV0.clear();
	cylinderV1.clear();
	cylinderRadii.clear();
	cylinderColors.clear();
}


void buildFloor(){
	// Build floor
	VEC3 color1 = VEC3(.3,.3,.3);
	VEC3 color2 = VEC3(.3,.3,.3);

  float y_plane = .2;
  float z_min = -5;
  float z_max = 5;
  float x_min = -5;
  float x_max = 5;

  vector<Primitive *> prims;

  //int count = 10; Makes every square = 1
  int count = 1;
  for(int i = 0; i < count; i++){
  	for(int j = 0; j < count; j++){
  		float i_1 = (i/(float)count)*(z_max - z_min) + z_min;
  		float i_2 = ((i+1)/(float)count)*(z_max - z_min) + z_min;
  		float j_1 = (j/(float)count)*(x_max - x_min) + x_min;
  		float j_2 = ((j+1)/(float)count)*(x_max - x_min) + x_min;
  		VEC3 corner_1 = VEC3(j_1, y_plane, i_1);
  		VEC3 corner_2 = VEC3(j_1, y_plane, i_2);
  		VEC3 corner_3 = VEC3(j_2, y_plane, i_1);
  		VEC3 corner_4 = VEC3(j_2, y_plane, i_2);

  		VEC3 thisColor;
  		if(((i+j)%2) == 0){
  			thisColor = color1;
  		}else{
  			thisColor = color2;
  		}

  		Tri * t1 = new Tri({corner_1, corner_3, corner_4}, thisColor);
  		Tri * t2 = new Tri({corner_1, corner_4, corner_2}, thisColor);
  		prims.push_back(t1);
  		prims.push_back(t2);
  	}
  }
  bounders.push_back(new BoundingBox(prims));
}

void buildRunway(float frame_num){
	float dotRad = .05;	
	float colRad = .01;
	float y_bottom = .2;
	float y_top = .23;
	float xMid = .45;
	float width = 1.5;
	float zMin = -1.0;
	float zMax = 1.5;
	int numDots = 10; //10;
	VEC3 dotColor = {0,0,1};
	VEC3 colColor = {.5,.5,.5};
	VEC3 lightColor = {0,0,.3};
	vector <Primitive *> prims;

	float intensity = 1.0;
	
	if(frame_num < 90){
		intensity = 0;
	}else if((frame_num >= 90) && (frame_num < 180)){
		// Pulsing blue
		float t = fmod(frame_num, 30.0f)/30.0f;
		if(t <= .5){
			intensity = smerpFloat(0.0f, 1.0f, t*2.0f);
		}else{
			intensity = smerpFloat(1.0f, 0.0f, (t-.5)*2);
		}
	}else if((frame_num >= 180) && (frame_num < 200)){
		intensity = 0;
	}else if((frame_num >= 200) && (frame_num < 255)){
		// Pulsing red
		dotColor = {1,0,0};
		lightColor = {.5,0,0};
		float t = fmod(frame_num, 10.0f)/10.0f;
		if(t <= .5){
			intensity = smerpFloat(0.0f, 1.0f, t*2.0f);
		}else{
			intensity = smerpFloat(1.0f, 0.0f, (t-.5)*2);
		}
	}else{
		dotColor = {1,0,0};
		lightColor = {.5,0,0};
		intensity = 1;
	}

	for(int i = 0; i < numDots; i++){
		float f = ((float)i/numDots);
		Sphere * s1 = new Sphere(VEC3(xMid+width/2.0f, y_top+dotRad, (zMax-zMin)*f + zMin), dotRad, dotColor);
		Sphere * s2 = new Sphere(VEC3(xMid-width/2.0f, y_top+dotRad, (zMax-zMin)*f + zMin), dotRad, dotColor);
		s1->setCastsShadows(false);
		s2->setCastsShadows(false);
		PointLight * l1 = new PointLight(VEC3(xMid+width/2.0f, y_top+dotRad, (zMax-zMin)*f + zMin), lightColor*intensity);
		PointLight * l2 = new PointLight(VEC3(xMid-width/2.0f, y_top+dotRad, (zMax-zMin)*f + zMin), lightColor*intensity);
		lights.push_back(l1);
		lights.push_back(l2);
  	prims.push_back(s1);
  	prims.push_back(s2);
		Cylinder * c1 = new Cylinder({VEC3(xMid+width/2.0f, y_bottom, (zMax-zMin)*f + zMin), VEC3(xMid+width/2.0f, y_top, (zMax-zMin)*f + zMin)}, colRad, colColor);
		Cylinder * c2 = new Cylinder({VEC3(xMid-width/2.0f, y_bottom, (zMax-zMin)*f + zMin), VEC3(xMid-width/2.0f, y_top, (zMax-zMin)*f + zMin)}, colRad, colColor);
		prims.push_back(c1);
		prims.push_back(c2);
	}
	bounders.push_back(new BoundingBox(prims));


}

void buildStickMan(){
	displayer.ComputeBonePositions(DisplaySkeleton::BONES_AND_LOCAL_FRAMES);

  // retrieve all the bones of the skeleton
  vector<MATRIX4>& rotations = displayer.rotations();
  vector<MATRIX4>& scalings  = displayer.scalings();
  vector<VEC4>& translations = displayer.translations();
  vector<float>& lengths     = displayer.lengths();

  // build a sphere list, but skip the first bone, 
  // it's just the origin
  vector<Primitive *> cyls;
  int totalBones = rotations.size();
  for (int x = 1; x < totalBones; x++)
  {
    MATRIX4& rotation = rotations[x];
    MATRIX4& scaling = scalings[x];
    VEC4& translation = translations[x];

    // get the endpoints of the cylinder
    VEC4 leftVertex(0,0,0,1);
    VEC4 rightVertex(0,0,lengths[x],1);

    leftVertex = rotation * scaling * leftVertex + translation;
    rightVertex = rotation * scaling * rightVertex + translation;

    Cylinder * c = new Cylinder({truncate(leftVertex), truncate(rightVertex)}, .05, VEC3(.5,.5,.5));
    cyls.push_back(c);

  }

  BoundingBox * b = new BoundingBox(cyls);
  bounders.push_back(b);

  /*
  Sphere * s1 = new Sphere(b->getMax(), .05, {1,0,0});
  s1->setAcceptsLights(false);
  Sphere * s2 = new Sphere(b->getMin(), .05, {1,0,0});
  s2->setAcceptsLights(false);
  primitives.push_back(s1);
  primitives.push_back(s2);*/
}


void writePPMToFolder(const string& folder, const string& filename, int& xRes, int& yRes, const float* values)
{
  int totalCells = xRes * yRes;
  unsigned char* pixels = new unsigned char[3 * totalCells];
  for (int i = 0; i < 3 * totalCells; i++)
    pixels[i] = values[i];

  system(("mkdir "+folder).c_str());
  FILE *fp;
  fp = fopen((folder + "/" +filename).c_str(), "wb");
  if (fp == NULL)
  {
    cout << " Could not open file \"" << filename.c_str() << "\" for writing." << endl;
    cout << " Make sure you're not trying to write from a weird location or with a " << endl;
    cout << " strange filename. Bailing ... " << endl;
    exit(0);
  }

  fprintf(fp, "P6\n%d %d\n255\n", xRes, yRes);
  fwrite(pixels, 1, totalCells * 3, fp);
  fclose(fp);
  delete[] pixels;
}



//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
void writePPM(const string& filename, int& xRes, int& yRes, const float* values)
{
  int totalCells = xRes * yRes;
  unsigned char* pixels = new unsigned char[3 * totalCells];
  for (int i = 0; i < 3 * totalCells; i++)
    pixels[i] = values[i];

  FILE *fp;
  fp = fopen(filename.c_str(), "wb");
  if (fp == NULL)
  {
    cout << " Could not open file \"" << filename.c_str() << "\" for writing." << endl;
    cout << " Make sure you're not trying to write from a weird location or with a " << endl;
    cout << " strange filename. Bailing ... " << endl;
    exit(0);
  }

  fprintf(fp, "P6\n%d %d\n255\n", xRes, yRes);
  fwrite(pixels, 1, totalCells * 3, fp);
  fclose(fp);
  delete[] pixels;
}


VEC3 starSphere(Ray r){
	VEC3 dir = r.d;
	//cout << "dir = " << dir << endl;
	VEC3 inplane = VEC3(dir[0], dir[1], 0).normalized();
	//cout << "inplane = " << inplane << endl;
	//cout << "dir = " << dir << "; inplane = " << inplane << endl;

	float theta;
	if(isnan(inplane[0])){
		theta = 0;
	}else{
		theta = acos(inplane.dot(VEC3(1,0,0)));
		float zval = inplane.cross(VEC3(1,0,0))[2];
		float sign_theta = -zval/abs(zval);
		if(!isnan(sign_theta) && (sign_theta < 0)) 
			theta = (2*PI) - theta;
	}
	

	float phi = acos(dir.dot(inplane));
	if(isnan(inplane[0])) phi = PI/2.0f;
	float sign_phi = dir[2]/abs(dir[2]);
	//cout << "sign_phi = " << sign_phi << endl;
	if(!isnan(sign_phi) && (sign_phi < 0)) 
		phi = -phi;

	float v = theta/(2.0f*PI);
	float u = (phi + PI/2.0f)/PI;
	cout << "v = "<< v << endl;
	VEC3 color = {0,0,0};
	//color += u*VEC3(1,0,0);
	color += v*VEC3(0,1,0);

	return color;
}


VEC3 getDustColor(Ray r, float t){
	if(!dust)
		return {0,0,0};
	VEC3 colorSum = {0,0,0};
	int samples = floor(t*dust_N);

	for(int i = 0; i < samples; i++){
		float dt = m.rand();
		dt *= t;
		VEC3 testPoint = r.o + dt*r.d;

		for(Light * l : lights){
			// Test to see light
			Ray shadow;
			shadow.o = testPoint;
			shadow.d = (l->getPos() - testPoint).normalized();
			VEC3 lightColor;
			if(l->isVisible(shadow, primitives, bounders, lightColor)){
				colorSum += l->getIntensity(testPoint)*lightColor;
			}
		}
	}
	
	return (colorSum/(20*dust_N));//*dustConst;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
void rayColor(const VEC3& rayPos, const VEC3& rayDir, VEC3& pixelColor) 
{
	// Find closest primitive
	pixelColor = VEC3(0,0,0);
	Ray r;
	r.d = rayDir;
	r.o = rayPos;
  int hitID = -1;
  float tMinFound = FLT_MAX;
  VEC3 hitNormal;
  Primitive *hit = NULL;

  // Intersect solo primitives
	for(int i = 0; i < primitives.size(); i++){
		Primitive * prim = primitives[i];
		float tMin;
		VEC3 thisNorm;
		if(prim->intersectRay(r, tMin, thisNorm)){
			if(tMin < tMinFound){
				hitID = i;
				hit = primitives[i]; // this is new
				tMinFound = tMin;
				hitNormal = thisNorm;
			}
		}
	}

	// Intersect bounded primitives
	for(int i = 0; i < bounders.size(); i++){
		if(bounders[i]->intersectRay(r)){
			// We hit the box!

			vector<Primitive *> primitives = bounders[i]->getPrimitives();
			for(int i = 0; i < primitives.size(); i++){
				Primitive * prim = primitives[i];
				float tMin;
				VEC3 thisNorm;
				if(prim->intersectRay(r, tMin, thisNorm)){

					if(tMin < tMinFound){
						hitID = i;
						hit = primitives[i]; // this is new
						tMinFound = tMin;
						hitNormal = thisNorm;
					}
				}
			}

		}
	}


	VEC3 surfaceColor;
	VEC3 dustColor = {0,0,0};
	if(!hit)
		tMinFound = 20;

	dustColor = getDustColor(r, tMinFound);

	if(hit){
		//cout << "hit something" << endl;
		surfaceColor = hit->getColor(r, tMinFound);
	}else{
		pixelColor += dustColor*1;
		//pixelColor = starSphere(r); // Look to star sphere here
		return;
	}


	// Lighting
	if(!hit->acceptsLights()){
		//cout << "t = " << tMinFound << endl;
		pixelColor = surfaceColor;
		pixelColor += dustColor;
		return;
	}
	VEC3 hitPoint = r.o + (tMinFound-.0001)*r.d;
	VEC3 lightCumulative = {0,0,0};
	for(Light * l : lights){
		// Do the diffuse/specular thing
		Ray shadow;
		shadow.o = hitPoint;
		shadow.d = (l->getPos() - hitPoint).normalized();
		VEC3 r = -shadow.d + 2*(shadow.d.dot(hitNormal))*hitNormal;
		r.normalize();
		VEC3 e = eye - hitPoint;
		e.normalize();
		//cout << "shadow.o = " << shadow.o << "\nshadow.d = " << shadow.d << "\n";
		VEC3 lightColor;
		if(l->isVisible(shadow, primitives, bounders, lightColor)){
			lightCumulative += lightColor*max(0,hitNormal.dot(shadow.d)); // Diffuse
			lightCumulative += lightColor*pow(max(0, e.dot(r)), PHONG); // Specular
			lightCumulative += ambient; // Ambient

			//cout << "light cumulative = \n" << lightColor << "\nshadow.d = \n" << shadow.d << "\nhitNormal = " << hitNormal << endl;
		}
	}
	pixelColor = hadamard(surfaceColor, lightCumulative);
	pixelColor += dustColor;

}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////
// Load up a new motion captured frame
//////////////////////////////////////////////////////////////////////////////////
void setSkeletonsToSpecifiedFrame(int frameIndex)
{
  if (frameIndex < 0)
  {
    printf("Error in SetSkeletonsToSpecifiedFrame: frameIndex %d is illegal.\n", frameIndex);
    exit(0);
  }
  if (displayer.GetSkeletonMotion(0) != NULL)
  {
    int postureID;
    if (frameIndex >= displayer.GetSkeletonMotion(0)->GetNumFrames())
    {
      cout << " We hit the last frame! You might want to pick a different sequence. " << endl;
      postureID = displayer.GetSkeletonMotion(0)->GetNumFrames() - 1;
    }
    else 
      postureID = frameIndex;
    displayer.GetSkeleton(0)->setPosture(* (displayer.GetSkeletonMotion(0)->GetPosture(postureID)));
  }
}


//////////////////////////////////////////////////////////////////////////////////
// Build a list of spheres in the scene
//////////////////////////////////////////////////////////////////////////////////
void buildScene(float frame_num)
{
	primitives.clear();
	bounders.clear();
	lights.clear();

	motion_N = floor(keyframeFloat(motionNKeys, motionNVals, motionNInterps, frame_num));
  dof_N = floor(keyframeFloat(dofNKeys, dofNVals, dofNInterps, frame_num));

	//PointLight * l = new PointLight({1,0,0}, {1,0,0});
  //lights.push_back(l);

	//return;
	// Keyframe motion data and camera
	dustConst = keyframeFloat(dustKeys, dustVals, dustInterps, frame_num);
  eye = keyframeVec(camKeyFrames, camEyeValues, camInterpolation, frame_num);
  lookingAt = keyframeVec(camKeyFrames, camLookValues, camInterpolation, frame_num);
  up = keyframeVec(camKeyFrames, camUpValues, camInterpolation, frame_num);
  aperature = keyframeFloat(camKeyFrames, aperatureValues, camInterpolation, frame_num);
  ambient = keyframeVec(ambientKeys, ambientValues, ambientInterps, frame_num);
	float dataIndex = keyframeFloat(stickKeyframes, stickDataIndexValues, stickInterpolation, frame_num);
  setSkeletonsToSpecifiedFrame(floor(dataIndex));

  

	for(Model m : models){
		m.addToScene(frame_num);
	}


	vector<float> pSpotKeys = {0,83,181};
	vector<VEC3> pSpotPoi = {{.5, .3, -.7}, {.5, .3, -.7}, {.4, .6, .85}};
	vector<int> pSpotInterps = {1,1,1};
	VEC3 pSpotLook = keyframeVec(pSpotKeys, pSpotPoi, pSpotInterps, frame_num);

	SpotLight * megananSpot1 = new SpotLight(VEC3(4.4,8,1.2), VEC3(.4, .3, 1.2), 5.0*PI/180.0, VEC3(1,1,1)); // Meganan spot (y = 8, x = -4)
	SpotLight * megananSpot2 = new SpotLight(VEC3(-3.6,8,1.2), VEC3(.4, .3, 1.2), 5.0*PI/180.0, VEC3(1,1,1)); // Meganan spot (y = 8, x = -4)
	SpotLight * personSpot1 = new SpotLight(VEC3(.5,8,2.2), pSpotLook, 5.0*PI/180.0, VEC3(1,1,1)); // Dude spot



	PointLight * point = new PointLight(VEC3(.45,1,.25), VEC3(.5,0,1));
	//PointLight * point = new PointLight(VEC3(0,4,0), VEC3(.8,0,0)/2);
	lights.push_back(megananSpot1);
	lights.push_back(megananSpot2);
	lights.push_back(personSpot1);
	lights.push_back(point);

	Sphere * s1 = new Sphere(VEC3(.5,2,-.7), .1, VEC3(1,0,1));
  s1->setAcceptsLights(false);
  s1->setCastsShadows(false);
  //primitives.push_back(s1);

  /*Sphere * s2 = new Sphere(VEC3(.3,.333,1.2), .1, VEC3(1,0,0));
  s2->setAcceptsLights(false);
  primitives.push_back(s2);*/


	buildStickMan();
	buildFloor();
	buildRunway(frame_num);

  // MEGANAN placeholder
  /*Sphere * s = new Sphere(VEC3(.4,.3,1.2), .1, VEC3(1,0,1));
  primitives.push_back(s);*/


	if(drawPOI){
		Sphere * poi = new Sphere(lookingAt, .05, VEC3(1,0,0));
		poi->setAcceptsLights(false);
		bounders.push_back(new BoundingBox({poi}));
	}
  
}


//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
void renderImage(int& xRes, int& yRes, const string& filename) 
{
  // allocate the final image
  const int totalCells = xRes * yRes;
  float* ppmOut = new float[3 * totalCells];



  // compute image plane
  const float halfY = (lookingAt - eye).norm() * tan(45.0f / 360.0f * M_PI);
  const float halfX = halfY * 4.0f / 3.0f;

  const VEC3 cameraZ = (lookingAt - eye).normalized();
  const VEC3 cameraX = up.cross(cameraZ).normalized();
  const VEC3 cameraY = cameraZ.cross(cameraX).normalized();



  for (int y = 0; y < yRes; y++) 
    for (int x = 0; x < xRes; x++) 
    {
      // generate the ray, making x-axis go left to right
      const float ratioX = 1.0f - ((xRes - 1) - x) / float(xRes) * 2.0f;
      const float ratioY = 1.0f - y / float(yRes) * 2.0f;
      const VEC3 rayHitImage = lookingAt + 
                               ratioX * halfX * cameraX +
                               ratioY * halfY * cameraY;
      
      VEC3 color;
      // Crazy depth of focus
      if(dof_on){
      	VEC3 dofColorSum = {0,0,0};
	      for(int i = 0; i < dof_N; i++){
	      	VEC3 thisColor;
	      	VEC3 thisEye = eye + (m.rand()*aperature-(aperature/2.0))*cameraX + (m.rand()*aperature-(aperature/2.0))*cameraY;

	      	VEC3 thisDir = (rayHitImage - thisEye).normalized();
	      	rayColor(thisEye, thisDir, thisColor);
	      	dofColorSum += thisColor;
	      	//cout << "eye = " << eye << endl;
	      	//cout << "thisEye = " << thisEye << endl;
	      }
      	color = dofColorSum/dof_N;
      }else{
      	const VEC3 rayDir = (rayHitImage - eye).normalized();
      	rayColor(eye, rayDir, color);
      }

      // set, in final image
      ppmOut[3 * (y * xRes + x)] = clamp(color[0]) * 255.0f;
      ppmOut[3 * (y * xRes + x) + 1] = clamp(color[1]) * 255.0f;
      ppmOut[3 * (y * xRes + x) + 2] = clamp(color[2]) * 255.0f;
    }
  writePPM(filename, xRes, yRes, ppmOut);

  delete[] ppmOut;
}

float * getSceneImage(int& xRes, int& yRes, int& startRow, int& endRow){
	// allocate the final image
  const int totalCells = xRes * yRes;
  float* ppmOut = new float[3 * totalCells]();

  // compute image plane
  const float halfY = (lookingAt - eye).norm() * tan(45.0f / 360.0f * M_PI);
  const float halfX = halfY * 16.0f / 9.0f;

  const VEC3 cameraZ = (lookingAt - eye).normalized();
  const VEC3 cameraX = up.cross(cameraZ).normalized();
  const VEC3 cameraY = cameraZ.cross(cameraX).normalized();

  for (int y = startRow; y < endRow; y++){
    for (int x = 0; x < xRes; x++) 
    {
      // generate the ray, making x-axis go left to right
      const float ratioX = 1.0f - ((xRes - 1) - x) / float(xRes) * 2.0f;
      const float ratioY = 1.0f - y / float(yRes) * 2.0f;
      const VEC3 rayHitImage = lookingAt + 
                               ratioX * halfX * cameraX +
                               ratioY * halfY * cameraY;
      
      VEC3 color;
      // Crazy depth of focus
      if(dof_on){
      	VEC3 dofColorSum = {0,0,0};
	      for(int i = 0; i < dof_N; i++){
	      	VEC3 thisColor;
	      	VEC3 thisEye = eye + (m.rand()*aperature-(aperature/2.0))*cameraX + (m.rand()*aperature-(aperature/2.0))*cameraY;
	      	VEC3 thisDir = (rayHitImage - thisEye).normalized();
	      	rayColor(thisEye, thisDir, thisColor);
	      	dofColorSum += thisColor;
	      }
      	color = dofColorSum/dof_N;
      }else{
      	const VEC3 rayDir = (rayHitImage - eye).normalized();
      	rayColor(eye, rayDir, color);
      }

      // set, in final image
      ppmOut[3 * (y * xRes + x)] = clamp(color[0]) * 255.0f;
      ppmOut[3 * (y * xRes + x) + 1] = clamp(color[1]) * 255.0f;
      ppmOut[3 * (y * xRes + x) + 2] = clamp(color[2]) * 255.0f;
		}
	}
	return ppmOut;
}

void printppm(float * in, int length){
	for(int i = 0; i < length; i++){
		cout << in[i] << ", ";
	}
	cout << endl;
}

float * addPPMs(float * a, float * b, int length){
	float* ppmOut = new float[3 * length];
	for(int i = 0; i < length; i++){
		ppmOut[i] = a[i] + b[i];
	}
	return ppmOut;
}

float * dividePPM(float * a, float div, int length){
	float* ppmOut = new float[3 * length];
	for(int i = 0; i < length; i++){
		ppmOut[i] = a[i]/div;
	}
	return ppmOut;
}


long int getSeconds(){
	struct timeval tp;
	gettimeofday(&tp, NULL);
	long int ms = tp.tv_sec * 1000 + tp.tv_usec / 1000;
	return ms;
}


long int start;

void renderScene(int& xRes, int& yRes, const string& filename, float frame_num, int& startRow, int& endRow){
	const int totalCells = xRes * yRes;
  float* ppmOut = new float[3 * totalCells];

  motion_N = floor(keyframeFloat(motionNKeys, motionNVals, motionNInterps, frame_num));
  dof_N = floor(keyframeFloat(dofNKeys, dofNVals, dofNInterps, frame_num));

	if((motion_N > 0) && motionBlur){
		float* ppmSum = new float[3 * totalCells]; for(int i = 0; i < totalCells; i++) ppmSum[i] = 0;
		for(int i = 0; i < motion_N; i++){
			float dt = m.rand();
			cout << "Builing frame: " << frame_num+dt << endl;
			buildScene(frame_num + dt);
			ppmSum = addPPMs(ppmSum, getSceneImage(xRes, yRes, startRow, endRow), totalCells*3);
			cout << "Time elapsed: " << getSeconds() - start << endl;
			//ppmSum = getSceneImage(xRes, yRes);
		}
		ppmOut = dividePPM(ppmSum, motion_N, totalCells*3);
		//CLAM! (again, just to be real safe)
		for(int i = 0; i < totalCells; i++){
			//if(ppmOut[i] > 255.0f){
			//	cout << "OH SHIT" << endl;
			//}
			ppmOut[i] = clamp(ppmOut[i]/255.0f)*255.0f;
		}
	}else{
		buildScene(frame_num);
		ppmOut = getSceneImage(xRes, yRes, startRow, endRow);
	}
	
	writePPM(filename, xRes, yRes, ppmOut);
	//return ppmOut;
}



//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{

	int startFrame, endFrame;
	int startRow = 0;
	int endRow = windowHeight;
	if(argc == 1){
		startFrame = 185;
		endFrame = 186;
		motionBlur = false;
		dof_on = false;
		dust_N = 25;
	}else if(argc == 2){
		startFrame = atoi(argv[1]);
		endFrame = 300;
	}else if(argc == 4){
		cout << "Please provide a start and end for pixel rows" << endl;
		exit(1);
	}else if(argc == 5){
		startFrame = atoi(argv[1]);
		endFrame = atoi(argv[2]);
		startRow = atoi(argv[3]);
		endRow = atoi(argv[4]);
	}else if(argc == 3){
		startFrame = atoi(argv[1]);
		endFrame = atoi(argv[2]);
	}

	m.seed();

  string skeletonFilename("resources/banana.asf");
  string motionFilename("resources/banana.amc");
  
  // load up skeleton stuff
  skeleton = new Skeleton(skeletonFilename.c_str(), MOCAP_SCALE);
  skeleton->setBasePosture();
  displayer.LoadSkeleton(skeleton);

  // load up the motion
  motion = new Motion(motionFilename.c_str(), MOCAP_SCALE, skeleton);
  displayer.LoadMotion(motion);
  skeleton->setPosture(*(displayer.GetSkeletonMotion(0)->GetPosture(0)));

  // load up models
  Model meganan;
  meganan.Load("resources/meganan.stl");
  meganan.setColor(VEC3(1,1,0));
  meganan.setPosKeys({0, 178, 183}, {VEC3(.4, .2, 1.2), VEC3(.4, .2, 1.2), VEC3(.4, 2.5, 3.9)}, {0, 0});
  meganan.setScaleKeys({0},{VEC3(.01, .01, .01)}, {0});
  meganan.setRotKeys({0, 178, 183}, {VEC3(-PI/2.0,0,0), VEC3(-PI/2.0,0,0), VEC3(2*PI,0,0)}, {0, 0});
  meganan.setActiveKeys({0,183}, {1, 0});
  models.push_back(meganan);

  Model glasses;
  glasses.Load("resources/glasses.stl");
  glasses.setPosKeys({0, 35, 50, 178, 183}, {VEC3(.4, .55, 1.2), VEC3(.4, .55, 1.2), VEC3(.4, .3, 1.2), VEC3(.4, .3, 1.2), VEC3(.4, 3, 3.4)}, {1, 1, 1, 0, 0});
  glasses.setScaleKeys({0}, {VEC3(.02,.02,.02)}, {0});
  glasses.setRotKeys({0, 35, 50, 178, 183}, {VEC3(-PI/2.0,0,-17.0*PI/180.0), VEC3(-PI/2.0,0,-17.0*PI/180.0), VEC3(-PI/2.0,0,-17.0*PI/180.0), VEC3(-PI/2.0,0,-17.0*PI/180.0), VEC3(2.0*PI,0,-17.0*PI/180.0)}, {1, 1, 1, 0, 0});
  glasses.setColor(VEC3(0,0,0));
  glasses.setActiveKeys({0,34, 35, 183}, {0,0,1, 0});
  models.push_back(glasses);

  // Right now we are set up to render individual high quality images
  // Uncomment the other code to render other frames
  // I only did this because I wanted a high quality desktop image of meganan and
  // it was taking way to long to render a single high quality frame and I wanted
  // to parallelize the process
    
    
  /*
  eye = VEC3(-3,.5,0);
  lookingAt = VEC3(1,0,0);
  up = VEC3(0,1,0);
	

	SpotLight * l = new SpotLight(VEC3(1, 1, 0), VEC3(1,0,0), 45.0*PI/180.0, VEC3(1,1,1));
	//PointLight * l = new PointLight(VEC3(1,1,0), VEC3(1,1,1));
	Sphere * s2 = new Sphere({1,.6, 0}, .05, {1,1,0});
	//s2->setAcceptsLights(false);
	primitives.push_back(s2);

	lights.push_back(l);

	Sphere * s = new Sphere({1,1,0}, .1, VEC3(1,0,0));
	s->setAcceptsLights(false);
	s->setCastsShadows(false);
	primitives.push_back(s);

  char buffer[256];
  sprintf(buffer, "./frames/frame.%04i.ppm", 0);
  renderScene(windowWidth, windowHeight, buffer, 0);

  Ray r;
  r.o = {-3,1,0};
  r.d = (l->getPos()-r.o).normalized();
  VEC3 col;
  if(l->isVisible(r, primitives, bounders, col))
  	cout<<"visible"<<endl;
  else
  	cout<<"not visible"<<endl;

  cout << "Rendered " + to_string(0) + " frames" << endl;*/

  // Note we're going 8 frames at a time, otherwise the animation
  // Hahahahahaha good one professor kim

	start = getSeconds();

  cout << "Rendering frames [" << startFrame << ", " << endFrame << ")" << endl;
  cout << "Rendering rows [" << startRow << ", " << endRow << ")" << endl;
  
  for(int frame_num = startFrame; frame_num < endFrame; frame_num++){
    char buffer[256];
    sprintf(buffer, "./frames/frame.%04i.%04i.ppm", frame_num, startRow);
    renderScene(windowWidth, windowHeight, buffer, frame_num, startRow, endRow);
    cout << "Rendered " + to_string(frame_num) + " frames" << endl;




    /*char frame[256];
	  sprintf(frame, "%04i", frame_num);
	  char startRowChar[256];
	  sprintf(startRowChar, "%04i", startRow);
	  char endRowChar[256];
	  sprintf(endRowChar, "%04i", endRow);

	  string numStr = frame;
	  string startRowStr = startRowChar;
	  string endRowStr = endRowChar;

	  string nameStr = startRowStr+"-"+endRowStr+".ppm";
	  string folderStr = "./baked/frame."+numStr;
		float * pixels = renderScene(windowWidth, windowHeight, nameStr, frame_num, startRow, endRow);
		writePPMToFolder(folderStr, nameStr, windowWidth, windowHeight, pixels);
		cout << "Rendered frame to" << folderStr + "/" + nameStr << endl;*/
  }

  cout << "Time elapsed: " << getSeconds() - start << endl;

  return 0;
}
