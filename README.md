# 3D Model Rendering Engine

Project for the course "Computer Graphics", a 3D Model Rendering Engine that makes complex model images out of a 
.init file. 

## Features

- **Simple Initialization:** Define your 3D models effortlessly in the initialization file, specifying parameters such as geometry, lighting, and camera settings.

- **Extensibility:** Easily extend the engine to support additional features and rendering capabilities, adapting it to suit diverse use cases.

- **2D fractals:**

- **3D Bodies:** The following 3D bodies are implemented.
    * Cube
    * Tetrahedron 
    * Octahedron 
    * Icosahedron 
    * Dodecahedron 
    * Cylinder
    * Cone 
    * Sphere 
    * Torus

- **3D fractals:**

- **Linedrawing with cylinders and spheres:** Lines and points of a line drawing are replaced by cylinders and points respectively.

- **Lighting:** Different lighting schemes are implemented.
    * Ambient light
    * Diffuse light (lightsource infinity)
    * Diffuse light (point source)
    * Specular light

- **Shadowing:** Shadowing is only supported for point source.


## Usage

1. **Clone the Repository:**
   ```bash
   git clone https://github.com/sampieters/engine.git

2. **Navigate to the project directory**
   ```bash
   cd engine
   
3. **Build the project**
    * Ensure you have CMake installed
    * Create a build directory
   ```bash
   mkdir build
   cd build
   ```
    * Run CMake to generate build files:
   ```bash
   cmake ..
   ```
    * Build the project
   ```bash
   make
   ```

4. **Run the engine**
   * Execute the engine with your initialization file (replace your_init_file.init with the actual filename):
   ```bash
   ./engine your_init_file.init
   ```

   
