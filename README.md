Tight Isoperimetric Profiles
============================

Prerequisites
-------------
```bash
conda create --name tightiso
conda activate tightiso
conda install -c conda-forge cgal==5.4
conda install pip
pip install matplotlib shapely
```

Compilation
-----------

The standard cmake compilation commands are sufficient:
```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
make
```
This may take a few minutes: CGAL and Boost rely heavily on template
programming, which has slow compilation times.