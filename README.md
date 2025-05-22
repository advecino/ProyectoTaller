# Initial Orbit Determination (Extended Kalman Filter)

## Descripción
Este proyecto es una traducción de un código original en MATLAB a C++. Implementa un filtro de Kalman extendido (EKF) para la determinación inicial de órbitas a partir de mediciones observacionales.

## Estructura del Proyecto
El código consta de múltiples archivos que han sido adaptados de MATLAB a C++. La estructura del proyecto es la siguiente:

```
Proyecto/
├── src/                    # Código fuente en C++
├── include/                # Archivos de cabecera
├── data/                   # Datos de entrada y archivos de prueba
├── EKF_Test.cpp            # Archivo de pruebas Unitarias
├── EKF_GEOS3.cpp           # Archivo principal
└── README.md               # Este archivo
```


## Autor
Este proyecto ha sido desarrollado por Adrian Vecino Durán

## Requisitos
- Compilador C++ compatible con C++11 o superior (`g++`, `clang++`, etc.)
- Sistema compatible con POSIX (Linux o WSL en Windows recomendado para entorno GCC online)

## Instrucciones iniciales
Descarga el archivo del proyecto.

Descomprimir el contenido en una carpeta local.



## Compilación (GCC en línea o terminal)

Para compilar la **aplicación principal**:

```bash
g++ -std=c++11 -Iinclude src/*.cpp EKF_GEOS3.cpp -o EKF_GEOS3
```
Para compilar los **Test Unitarios**:
```bash
g++ -std=c++11 -Iinclude src/*.cpp EKF_Test.cpp -o EKF_Test
```