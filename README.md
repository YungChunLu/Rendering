# Rendering Course @ NTU, CSIE

## Assignment 1 - Height field
  * [Assignment Link](http://www.csie.ntu.edu.tw/~cyy/courses/rendering/13fall/assignments/proj1/)
  * In this assignment, you are asked to improve pbrt's heightfield class. A heightfield is a 2D array which samples a height function z(x, y). It indicates how high the sampled point (x, y) is and can be used for modeling things such as terrain.

## Assignment 2 - Realistic Camera Model
  * [Assignment Link](http://www.csie.ntu.edu.tw/~cyy/courses/rendering/13fall/assignments/proj2/)
  * Most computer graphics algorithms assume pin-hole camera model. However, it can't capture some important characteristics of lens system used in most modern cameras, usch as depth of field, distortion, vignetting and sptaillay varying exposure. In this project, you will implement the realistic camera model proposed by Kolb et. al. in SIGGRAPH 1995.
  
## Assignment 3 - Environment Lights
  * [Assignment Link](http://www.csie.ntu.edu.tw/~cyy/courses/rendering/13fall/assignments/proj3/)
  * For this project, it's recommend the median cut algorithm proposed by Debevec for its simplicity. This method recursively subdivides the light probe image into equal-energy regions. Finally, a representive point light is created for each region at its energy centroid. In the example of the following image, the Grace Cathedral light probe is subdivided into 64 regions with roughly equal energy values and 64 point lights are created to represent the environment light.
  * ![Alt](http://www.csie.ntu.edu.tw/~cyy/courses/rendering/13fall/assignments/proj3/mediancut.jpg)

## Assignment 4 - Final
  * Since Google released the first explorer wearable device in Feb 2013, more and more people has jumped into this area to develope their ideas about wearable computing devices. Moreover this device could provide a virtual screen which would be floating on the air and combining with the real scene to users. So it inspires me to use the techniques I learned from this course to simulate this amazing feature.
  * ![Alt](https://github.com/ArthurLu/Rendering/blob/master/Final/Result/algorithm.jpg)
