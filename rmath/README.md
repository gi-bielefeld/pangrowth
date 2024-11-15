# Rmath - The R Standalone Math Library

The repository contains the [R standalone Math Library](https://cran.r-project.org/doc/manuals/r-release/R-admin.html#The-standalone-Rmath-library), built to be used with the [Statslabs.Stats](https://github.com/statslabs/stats) C++ library.

## Installation on Ubuntu / macOS
1. Clone the repository.
   ```sh
   git clone git@github.com:statslabs/rmath.git
   ```
2. Configure the project.
   ```sh
   cd rmath
   mkdir build && cd build
   cmake ..
   ```
3. Compile and install the library.
   ```sh
   make
   make install
   ```

## Example program
`demo.c`:
```c
#include <stdio.h>
#include "Rmath.h"

int main() {
  double shape1, shape2, prob;

  printf("Enter first shape parameter: ");
  scanf("%lf",&shape1);

  printf("Enter second shape parameter: ");
  scanf("%lf",&shape2);

  printf("Enter probability level: ");
  scanf("%lf",&prob);

  printf("Critical value is %lf\n",qbeta(prob,shape1,shape2,1,0));

  return 0;
}
```

## Integration of Rmath in your own project
To make the project simple enough, we will create a CMake project for `demo.c`.

1. Make a project folder.
   ```sh
   mkdir example && cd example
   ```

2. Create `demo.cc` and `CMakeLists.txt` in the project folder where file `CMakeLists.txt` should look like:
   ```cmake
   cmake_minimum_required(VERSION 3.0)
   project(example)
   add_executable(example demo.c)

   find_package(Rmath 1.0.0 REQUIRED)
   target_link_libraries(example Rmath::Rmath)
   ```

3. Perform a out-of-source build.
   ```sh
   mkdir build && cd build
   cmake ..
   make
   ```

4. Run the program.
   ```sh
   ./example
   ```

   A quick check on whether the program works:
   ```
   Enter first shape parameter: 1.5
   Enter second shape parameter: 0.3
   Enter probability level: .10
   Critical value is 0.478398
   ```

   This corresponds to the following call in R:
   ```r
   > qbeta(.10,1.5,.3)
   [1] 0.4783981
   ```