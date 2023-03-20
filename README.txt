https://youtu.be/mjtvSp5ULjw
meganan.mp4

I have implemented depth of field, motion blur, and dust (Custom effect). Depth of field simulates a physical camera with a non-infinetly small aperature, motion blur simulates exposure over time between frames in a physical camera, and dust simulates randomly dispersed dust particles in a room allowing light beams to be seen.

BUILD INSTRUCTIONS
unzip Walak_FinalProject.zip
make using the command 'make all' (Quotes not included)
To render frames [a,b), use './meganan a b' (This will take hours... please don't do this)
To render a representative frame, use './meganan' (No arguments, takes ~1 min on the zoo)
note: frames are output to /frames/frame.xxxx.ppm

Credits:
Banana 3D model by sydd: https://sketchfab.com/3d-models/mario-kart-banana-peel-df7e6bc337d7472ab12fe3ec5004409c
Sunglasses 3D model by printable_models: https://free3d.com/3d-model/sunglasses-v1--803862.html
stl import code by dillonhuff: https://github.com/dillonhuff/stl_parser (Used with license)