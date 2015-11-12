#default {
  finish { ambient 0.4 }
}

#include "colors.inc"
#include "textures.inc"

camera {location <-0,-10,-10> sky <0,0,-1> look_at 0}
light_source {<-10,-3,-18>, color White}
background { color rgb <1.0, 1.0, 1.0> }
//background { color rgb Black }

#declare arr = union {
	cylinder {0,<0.5,0,0>, 0.02}
	cone{<0.5,0,0>, 0.04, <0.7,0,0>, 0 }  		
	pigment {color Blue}
	scale <1,2,1>
//	finish{phong 1}
} 

#declare axis = union {
	cylinder {0,<5,0,0>, 0.02}
	cone{<0.5,0,0>, 0.04, <0.7,0,0>, 0 }  		
	pigment {color Black}
	scale <1,2,1>
//	finish{phong 1}
} 

#declare chrome = texture { pigment { color rgb <0.6,0.6,0.6> } finish { specular 0.5 reflection 0.7 roughness 0.005 metallic } }
#declare clearred = texture { pigment { color rgbf <0.75,0,0,0.7> } finish { diffuse 0.6 specular 0.4 } }
#declare greenpaint = texture { pigment { color rgb <0 1 0> } finish { specular 0 reflection 0 diffuse 0.7 } }
#declare offwhitepaint = texture { pigment { color rgb <0.9 0.9 0.8> } finish { specular 0 reflection 0 diffuse 0.7 } }
#declare blackrubber = texture { pigment { color rgb <0.2 0.2 0.2> } finish { specular 0 reflection 0 diffuse 0.7 } }
#declare yellowplastic = texture { pigment { color rgb <0.95 0.95 0> } finish { specular 0.2 diffuse 0.8 } }


#declare push_pin =
union {
  cylinder { <0,0,-0.3> <0,0,0.25> 0.018 texture { chrome } }
  sphere { <0,0,0> 0.05 scale <1,1,0.1 > translate <0,0,-0.3> texture { chrome } }
  cone { <0,0,0.25> 0.018 <0,0,0.375> 0 texture { chrome } }
  difference {
    sphere { <0,0,0.01> 0.17 scale <1,1,0.8> }
    box { <-1,-1,0> <1,1,1> }
  }
  cylinder { <0,0,-0.01> <0,0,-0.4> .09 }
  cone { <0,0,-0.375> .11 <0,0,-0.44> .125 }
  difference {
    sphere { <0,0,0.275> 0.3 }
    box { <-1,-1,0> <1,1,1> }
    translate <0,0,-0.44>
  }
  interior { ior 1.5 }
}

union{	
//box {<0,0,0>, <-10.5,-8,-0.0> pigment { color 1.2*White}}
//box {<0,0,0>, <-10.5,-0.0,-5>pigment { color 1.5*White}}
//box {<0,0,0>, <-0.0,-8,-5> pigment { color 1.5*White}}
// grid coordinates
//cylinder {0,<-8,0,0>, 0.05}
//cylinder {0,<0,-6,0>, 0.05}
//cylinder {0,<0,0,-5>, 0.05}
//cylinder {<-8,0,0>,<-8,-6,0>, 0.03}
//cylinder {<0,-6,0>,<-8,-6,0>, 0.03}
//cylinder {<-8,0,0>,<-8,0,-5>, 0.03}
//cylinder {<0,0,-5>,<-8,0,-5>, 0.03}
//cylinder {<0,-6,0>,<0,-6,-5>, 0.03}
//cylinder {<0,0,-5>,<0,-6,-5>, 0.03}
//union{
//object{axis  rotate<180,0,0>}
//object{axis  rotate<0,180,00>}
//object{axis  rotate<0,0,90>}
//translate<-10,-4,-4>
//}



 union{  
 	sphere{<4.2382688749478120E-003,-4.9855054240550783E-002,0>,0.3 
 	texture{pigment{color 0.5*Magenta}} finish{phong 1}}
 	
 	object { push_pin 
	translate<0,0,-0.4> 
	scale 1.3
 	texture { greenpaint } }

//	object{arr translate < 2.0, 2.0, 0> scale 2.5}
	union{	
		sphere{<0.839759, -0.669981, 0.283580>	0.2}
		sphere{<1.776648, -0.856919, 0.579033>,	0.2 }
		sphere{<2.729676, -1.013058, 0.319497>,	0.2}
		sphere{<3.588657, -1.466057, 0.080869>,	0.2}
		sphere{<4.566778, -1.274397, -0.000029>,	0.2}
		sphere{<5.548086, -1.466831, 0.002032>,	0.2}
		sphere{<6.246263, -0.779688, -0.198925>,	0.2}
		sphere{<6.545694, 0.042727, 0.284785>,	0.2}
		sphere{ <5.579506, 0.044207, 0.026952>,	0.2}
		sphere{<4.859897, -0.206367, -0.620640>,	0.2}
		sphere{<3.917354, -0.428664, -0.371248>,	0.2}
		sphere{<2.964692, -0.196729, -0.567821>,	0.2}
		sphere{<1.970323, -0.296468, -0.531994>,	0.2}
		sphere{ <0.980422, -0.169988, -0.467978>,	0.2}
		texture{pigment{color Cyan}}
		finish{phong 1}
//		normal { bumps 0.4 scale 0.2 }
	}

//		union{	
//		sphere{<0.991352, 0.179832, 0.455388>,	0.2}
//		sphere{<1.978902, 0.337039, 0.449802>,	0.2 }
//		sphere{<2.922226, 0.355975, 0.781135>,	0.2}
//		sphere{<3.907897, 0.424925, 0.627196>,	0.2}
//		sphere{<4.907107, 0.457357, 0.650151>,	0.2}
//		sphere{<5.738312, 0.943304, 0.920254>,	0.2}
//		sphere{<6.635515, 1.041961, 0.489798>,	0.2}
//		sphere{<6.251149, 1.795838, -0.043055>,	0.2}
//		sphere{ <5.595824, 1.372292, -0.668482>,	0.2}
//		sphere{<4.630903, 1.151971, -0.811255>,	0.2}
//		sphere{<3.683918, 1.061149, -0.503080>,	0.2}
//		sphere{<2.727817, 0.898941, -0.259031>,	0.2}
//		sphere{<1.790523, 0.554723, -0.313745>,	0.2}
//		sphere{ <0.793960, 0.635282, -0.294443>,	0.2}
//		texture{pigment{color Black}}
//		finish{phong 1}
////		normal { bumps 0.4 scale 0.2 }
//	}
	union{	
		cylinder {<0.117377, -0.031596, 0.017820>,<0.839759, -0.669981, 0.283580>, 0.1}
		cylinder {<0.839759, -0.669981, 0.283580>,<1.776648, -0.856919, 0.579033>, 0.1}
		cylinder {<1.776648, -0.856919, 0.579033>,<2.729676, -1.013058, 0.319497>, 0.1}
		cylinder {<2.729676, -1.013058, 0.319497>,<3.588657, -1.466057, 0.080869>, 0.1}
		cylinder {<3.588657, -1.466057, 0.080869>,<4.566778, -1.274397, -0.000029>, 0.1}
		cylinder {<4.566778, -1.274397, -0.000029>,<5.548086, -1.466831, 0.002032>, 0.1}
		cylinder {<5.548086, -1.466831, 0.002032>,<6.246263, -0.779688, -0.198925>, 0.1}
		cylinder {<6.246263, -0.779688, -0.198925>,<6.545694, 0.042727, 0.284785>, 0.1}
		cylinder {<6.545694, 0.042727, 0.284785>,<5.579506, 0.044207, 0.026952>, 0.1}
		cylinder {<5.579506, 0.044207, 0.026952>,<4.859897, -0.206367, -0.620640>, 0.1}
		cylinder {<4.859897, -0.206367, -0.620640>,<3.917354, -0.428664, -0.371248>, 0.1}
		cylinder {<3.917354, -0.428664, -0.371248>,<2.964692, -0.196729, -0.567821>, 0.1}
		cylinder {<2.964692, -0.196729, -0.567821>,<1.970323, -0.296468, -0.531994>, 0.1}
		cylinder {<1.970323, -0.296468, -0.531994>,<0.980422, -0.169988, -0.467978>, 0.1}
		cylinder { <0.980422, -0.169988, -0.467978>,<0.117377, -0.031596, 0.017820>, 0.1}

		texture{pigment{color rgb<0.1,0.1,1>}}
		finish{
			ambient .2
			diffuse .6
			specular .75
			roughness .001
			reflection 0.03
			}
		finish{phong 0.3}
	}

//	union{	
//		cylinder {<0.117377, -0.031596, 0.017820>,<0.991352, 0.179832, 0.455388>, 0.1}
//		cylinder {<0.991352, 0.179832, 0.455388>,<1.978902, 0.337039, 0.449802>, 0.1}
//		cylinder {<1.978902, 0.337039, 0.449802>,<2.922226, 0.355975, 0.781135>, 0.1}
//		cylinder {<2.922226, 0.355975, 0.781135>,<3.907897, 0.424925, 0.627196>, 0.1}
//		cylinder {<3.907897, 0.424925, 0.627196>,<4.907107, 0.457357, 0.650151>, 0.1}
//		cylinder {<4.907107, 0.457357, 0.650151>,<5.738312, 0.943304, 0.920254>, 0.1}
//		cylinder {<5.738312, 0.943304, 0.920254>,<6.635515, 1.041961, 0.489798>, 0.1}
//		cylinder {<6.635515, 1.041961, 0.489798>,<6.251149, 1.795838, -0.043055>, 0.1}
//		cylinder {<6.251149, 1.795838, -0.043055>,<5.595824, 1.372292, -0.668482>, 0.1}
//		cylinder {<5.595824, 1.372292, -0.668482>,<4.630903, 1.151971, -0.811255>, 0.1}
//		cylinder {<4.630903, 1.151971, -0.811255>,<3.683918, 1.061149, -0.503080>, 0.1}
//		cylinder {<3.683918, 1.061149, -0.503080>,<2.727817, 0.898941, -0.259031>, 0.1}
//		cylinder {<2.727817, 0.898941, -0.259031>,<1.790523, 0.554723, -0.313745>, 0.1}
//		cylinder {<1.790523, 0.554723, -0.313745>,<0.793960, 0.635282, -0.294443>, 0.1}		cylinder {<0.793960, 0.635282, -0.294443>,<0.117377, -0.031596, 0.017820>, 0.1}
//		texture{pigment{color Orange}}
//		finish{
//			ambient .2
//			diffuse .6
//			specular .75
//			roughness .001
//			reflection 0.1
//			}
////		finish{phong 0.3}
//	}

	translate<-5,-4,-3>
 }
translate<2,0,-2>
}


////box { <-15,0,-15> <15,-1,15> texture { offwhitepaint } }
//cylinder { <-5,0,0> <5,0,0> 0.005 texture { greenpaint } }
//cylinder { <0,0,-5> <0,0,-0.5> 0.005 texture { greenpaint } }
//cylinder { <0,0,5> <0,0,0.4> 0.005 texture { greenpaint } }
//text {
//  ttf "cyrvetic.ttf" "Z" 0.1 0
//  scale 0.25
//  rotate 90*y
//  translate <0.05,0,0.75>
//  texture { greenpaint }
//}
//text {
//  ttf "cyrvetic.ttf" "X" 0.1 0
//  scale 0.25
//  rotate 90*y
//  translate <1.7,0,-0.05>
//  texture { greenpaint }
//}


