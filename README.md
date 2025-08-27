
# Forest-Fire-Simulator

A simple forest fire simulator built in **C** using **OpenMP** for parallelism and **SDL3** for graphics.

---

## Features
- Parallelized simulation of forest fire spread using OpenMP
- Visualization with SDL3
- Configurable parameters (forest size, fire spread probability, etc.)
- Lightweight and cross-platform

---

## 🛠Requirements

### Windows
- [MinGW](https://www.mingw-w64.org/) or MSVC (Visual Studio) for compiling
- SDL3 and SDL3_image (already included in `SDL3/` and `SDL3_image/` folders)

### Linux / macOS
Install SDL3 and SDL3_image development libraries:
bash
sudo apt install libsdl3-dev libsdl3-image-dev

Build Instructions
gcc forest_fire.c -o forest_fire.exe -I SDL3/include -L SDL3/lib -lSDL3 -I SDL3_image/include -L SDL3_image/lib -lSDL3_image -fopenmp     #windows

gcc forest_fire.c -o forest_fire -lSDL3 -lSDL3_image -fopenmp    #linux


Run
./forest_fire   # Linux / macOS
forest_fire.exe # Windows






