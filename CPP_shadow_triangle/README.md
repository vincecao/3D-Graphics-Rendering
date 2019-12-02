Lineng Cao

6298279121

linengca@usc.edu

visual studio 2017 ver 15.9.16

---
Need to switch to "release" instead of "debug" for building

Need Manual change  `valueListShader[5] = (GzPointer)0`, `(GzPointer)(tex_fun)` OR `(GzPointer)(ptex_fun)`in `Application5.cpp` to get texture.