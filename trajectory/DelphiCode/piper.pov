  #include "colors.inc"
  #include "textures.inc"
#include "shapes.inc"
light_source { < 10, 10, -10 >
color White
}             
light_source { < 10, 10, 10 >
color White                     
}                                  
background { color White }
camera {
 right      < -1.33, 0, 0 > 
 up         < 0, 1, 0 >
 direction  < 0, 0, 1 >
 location   < 15, 15, 45 >
 look_at    < 0,0, 20 >
}
cylinder { <0,0,0>,<0,0,100>,           0.1 pigment { color Blue } }
cylinder { <           0.000,0.000,0.000>,
           <0.006,-0.000,1.000>,
           0.1 pigment { color Red } }
cylinder { <           0.006,-0.000,1.000>,
           <0.028,-0.016,2.000>,
           0.1 pigment { color Red } }
cylinder { <           0.028,-0.016,2.000>,
           <0.038,-0.036,2.999>,
           0.1 pigment { color Red } }
cylinder { <           0.038,-0.036,2.999>,
           <0.045,-0.045,3.999>,
           0.1 pigment { color Red } }
cylinder { <           0.045,-0.045,3.999>,
           <0.055,-0.045,4.999>,
           0.1 pigment { color Red } }
cylinder { <           0.055,-0.045,4.999>,
           <0.063,-0.037,5.999>,
           0.1 pigment { color Red } }
cylinder { <           0.063,-0.037,5.999>,
           <0.083,-0.008,6.999>,
           0.1 pigment { color Red } }
cylinder { <           0.083,-0.008,6.999>,
           <0.105,0.013,7.998>,
           0.1 pigment { color Red } }
cylinder { <           0.105,0.013,7.998>,
           <0.077,0.032,8.998>,
           0.1 pigment { color Red } }
cylinder { <           0.077,0.032,8.998>,
           <0.032,0.042,9.996>,
           0.1 pigment { color Red } }
cylinder { <           0.032,0.042,9.996>,
           <-0.010,0.080,10.995>,
           0.1 pigment { color Red } }
cylinder { <           -0.010,0.080,10.995>,
           <0.009,0.119,11.994>,
           0.1 pigment { color Red } }
cylinder { <           0.009,0.119,11.994>,
           <0.044,0.162,12.992>,
           0.1 pigment { color Red } }
cylinder { <           0.044,0.162,12.992>,
           <0.089,0.198,13.991>,
           0.1 pigment { color Red } }
cylinder { <           0.089,0.198,13.991>,
           <0.063,0.463,14.955>,
           0.1 pigment { color Red } }
cylinder { <           0.063,0.463,14.955>,
           <0.047,0.729,15.919>,
           0.1 pigment { color Red } }
cylinder { <           0.047,0.729,15.919>,
           <0.049,0.997,16.882>,
           0.1 pigment { color Red } }
cylinder { <           0.049,0.997,16.882>,
           <0.080,1.276,17.842>,
           0.1 pigment { color Red } }
cylinder { <           0.080,1.276,17.842>,
           <0.114,1.564,18.799>,
           0.1 pigment { color Red } }
cylinder { <           0.114,1.564,18.799>,
           <0.155,1.839,19.759>,
           0.1 pigment { color Red } }
cylinder { <           0.155,1.839,19.759>,
           <0.159,2.147,20.711>,
           0.1 pigment { color Red } }
cylinder { <           0.159,2.147,20.711>,
           <0.117,2.449,21.663>,
           0.1 pigment { color Red } }
cylinder { <           0.117,2.449,21.663>,
           <0.045,2.749,22.614>,
           0.1 pigment { color Red } }
cylinder { <           0.045,2.749,22.614>,
           <-0.053,3.052,23.562>,
           0.1 pigment { color Red } }
cylinder { <           -0.053,3.052,23.562>,
           <-0.133,3.358,24.511>,
           0.1 pigment { color Red } }
cylinder { <           -0.133,3.358,24.511>,
           <-0.295,3.648,25.454>,
           0.1 pigment { color Red } }
cylinder { <           -0.295,3.648,25.454>,
           <-0.471,3.927,26.398>,
           0.1 pigment { color Red } }
cylinder { <           -0.471,3.927,26.398>,
           <-0.692,4.199,27.335>,
           0.1 pigment { color Red } }
cylinder { <           -0.692,4.199,27.335>,
           <-0.906,4.444,28.280>,
           0.1 pigment { color Red } }
cylinder { <           -0.906,4.444,28.280>,
           <-1.312,4.746,29.143>,
           0.1 pigment { color Red } }
cylinder { <           -1.312,4.746,29.143>,
           <-1.718,5.066,29.999>,
           0.1 pigment { color Red } }
cylinder { <           -1.718,5.066,29.999>,
           <-2.110,5.412,30.851>,
           0.1 pigment { color Red } }
cylinder { <           -2.110,5.412,30.851>,
           <-2.521,5.741,31.702>,
           0.1 pigment { color Red } }
cylinder { <           -2.521,5.741,31.702>,
           <-2.631,6.174,32.596>,
           0.1 pigment { color Red } }
cylinder { <           -2.631,6.174,32.596>,
           <-2.719,6.598,33.498>,
           0.1 pigment { color Red } }
cylinder { <           -2.719,6.598,33.498>,
           <-2.766,7.024,34.401>,
           0.1 pigment { color Red } }
cylinder { <           -2.766,7.024,34.401>,
           <-2.809,7.435,35.311>,
           0.1 pigment { color Red } }
cylinder { <           -2.809,7.435,35.311>,
           <-2.776,7.851,36.220>,
           0.1 pigment { color Red } }
cylinder { <           -2.776,7.851,36.220>,
           <-2.762,8.236,37.143>,
           0.1 pigment { color Red } }
cylinder { <           -2.762,8.236,37.143>,
           <-2.807,8.611,38.069>,
           0.1 pigment { color Red } }
cylinder { <           -2.807,8.611,38.069>,
           <-2.840,8.958,39.006>,
           0.1 pigment { color Red } }
cylinder { <           -2.840,8.958,39.006>,
           <-2.881,9.331,39.933>,
           0.1 pigment { color Red } }
cylinder { <           -2.881,9.331,39.933>,
           <-2.814,9.699,40.861>,
           0.1 pigment { color Red } }
cylinder { <           -2.814,9.699,40.861>,
           <-2.896,10.019,41.805>,
           0.1 pigment { color Red } }
cylinder { <           -2.896,10.019,41.805>,
           <-2.560,10.476,42.628>,
           0.1 pigment { color Red } }
cylinder { <           -2.560,10.476,42.628>,
           <-2.343,10.964,43.473>,
           0.1 pigment { color Red } }
cylinder { <           -2.343,10.964,43.473>,
           <-1.907,11.415,44.252>,
           0.1 pigment { color Red } }
cylinder { <           -1.907,11.415,44.252>,
           <-1.690,11.934,45.079>,
           0.1 pigment { color Red } }
cylinder { <           -1.690,11.934,45.079>,
           <-2.019,12.455,45.866>,
           0.1 pigment { color Red } }
cylinder { <           -2.019,12.455,45.866>,
           <-2.641,12.747,46.593>,
           0.1 pigment { color Red } }
cylinder { <           -2.641,12.747,46.593>,
           <-2.977,13.195,47.421>,
           0.1 pigment { color Red } }
cylinder { <           -2.977,13.195,47.421>,
           <-3.646,13.554,48.073>,
           0.1 pigment { color Red } }
cylinder { <           -3.646,13.554,48.073>,
           <-3.844,14.152,48.849>,
           0.1 pigment { color Red } }
cylinder { <           -3.844,14.152,48.849>,
           <-3.349,14.720,49.507>,
           0.1 pigment { color Red } }
