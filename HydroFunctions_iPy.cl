//
// HydroFunctions_iPy.cl
//
// Created by Johan Van de Koppel on 31/08/2017.
// Copyright 2017 Johan Van de Koppel. All rights reserved.
//
// Updated by Roeland C. van de Vijsel on 19/09/2023 - Fully checked (OK.)
// Royal Netherlands Institute for Sea Research (NIOZ), Department of Estuarine and Delta Systems, Yerseke, The Netherlands.
//
// This script defines functions that are used in the main numerical model (clPy.ComplexChannelPatterns_...).
// The model simulates saltmarsh development and is used for the following manuscript:
// van de Vijsel, R.C., van Belzen, J., Bouma, T.J., van der Wal, D., Borsje, B.W., Temmerman, S., Cornacchia, L., 
// Gourgue, O., van de Koppel, J. (2023, Nature Communications). Vegetation controls on channel network complexity in coastal 
// wetlands.
// Whenever you use any part of this model code, please make sure to correctly refer to this manuscript and its authors.

////// IMPORTANT INFORMATION //////
// Please note that in this model, X-coordinate increases from left to right and Y-coordinate increases from top to bottom.
// However, in the manuscript, the Y-coordinate increases from left to right and X-coordinate increases from top to bottom.
// Both in this model and in the manuscript, flow component u is in the X-direction and flow component v is in the Y-direction.
// Hence, coordinates X and Y in this model script should be interpreted as Y and X, respectively, in the manuscript.
// Idem, flow components u and v in this model should be interpreted as v and u, respectively, in the manuscript.
// However, this reversal is done consistently throughout the model and the equations are implemented correctly.
// Therefore, the model output is not affected.
///////////////////////////////////

#ifndef HYDROFUNCTIONS_CL
#define HYDROFUNCTIONS_CL

////////////////////////////////////////////////////////////////////////////////
// Apriori prototyping
////////////////////////////////////////////////////////////////////////////////

float d2_dxy2( __global float* );
float d_dx(__global float* );
float d_dy(__global float* );
float dO_dx(__global float*, __global float*, __global float* );
float dO_dy( __global float*, __global float*, __global float* );
float d_uh_dx(__global float*, __global float* );
float d_vh_dy( __global float*, __global float* );
void PeriodicBoundaries( __global float* );
void NeumanBoundaries(__global float* );
void ReflectingBoundaries(__global float*, __global float* );
void PersistentFluxBoundaries(__global float* );
void DirichletBoundaries(__global float*, float );

typedef unsigned int Thread_Id;
#define ON              1
#define OFF             0

////////////////////////////////////////////////////////////////////////////////
// Laplacation operator definition, to calculate diffusive terms
////////////////////////////////////////////////////////////////////////////////

float d2_dxy2(__global float* z)
{
    const float dx = dX;  // Forcing dX to become a float
    const float dy = dY;  // Forcing dY to become a float
    
    // Determine the position at which the current thread is computing
    const size_t current = get_global_id(0);            
    const size_t row     = (size_t)current/(size_t)Grid_Width; // Rounds down to integer, so row = 0 if current < Grid_Width, etc. 
    const size_t column  = (size_t)current%(size_t)Grid_Width;

    // Determine the positions neigboring the current position
    const size_t left    =  row      * Grid_Width + column - 1;
    const size_t right   =  row      * Grid_Width + column + 1;
    const size_t top     = (row - 1) * Grid_Width + column    ;
    const size_t bottom  = (row + 1) * Grid_Width + column    ;
    
    return ( (z[left] + z[right ] - 2.0*z[current])/dx/dx +
             (z[top ] + z[bottom] - 2.0*z[current])/dy/dy );
}

// (OK.)

////////////////////////////////////////////////////////////////////////////////
// Gradient operator definitions, to calculate advective fluxes
////////////////////////////////////////////////////////////////////////////////

float d_dx(__global float* z)
{
    const float dx = dX;  // Forcing dX to become a float
    
    // Determine the position at which the current thread is computing
    const size_t current = get_global_id(0);            
    const size_t row     = (size_t)current/(size_t)Grid_Width; // Rounds down to integer, so row = 0 if current < Grid_Width, etc. 
    const size_t column  = (size_t)current%(size_t)Grid_Width;
    
    // Determine the positions neigboring the current position
    const size_t left    =  row      * Grid_Width + column - 1;
    const size_t right   =  row      * Grid_Width + column + 1;
    
    return ( ( z[right] - z[left] )/2.0/dx );
}

// (OK.)

float d_dy(__global float* z)
{
    const float dy = dY;  // Forcing dY to become a float
    
    // Determine the position at which the current thread is computing
    const size_t current = get_global_id(0);            
    const size_t row     = (size_t)current/(size_t)Grid_Width; // Rounds down to integer, so row = 0 if current < Grid_Width, etc. 
    const size_t column  = (size_t)current%(size_t)Grid_Width;
    
    // Determine the positions neigboring the current position
    const size_t top     = (row - 1) * Grid_Width + column    ;
    const size_t bottom  = (row + 1) * Grid_Width + column    ;
    
    return ( ( z[bottom] - z[top] )/2.0/dy );    
}

// (OK.)

////////////////////////////////////////////////////////////////////////////////
// Gradient operator definitions, to calculate the pressure gradient
////////////////////////////////////////////////////////////////////////////////

float dO_dx(__global float* h, __global float* s,__global float* b)
{
    const float dx = dX;  // Forcing dX to become a float
    
    // Determine the position at which the current thread is computing
    const size_t current = get_global_id(0);            
    const size_t row     = (size_t)current/(size_t)Grid_Width; // Rounds down to integer, so row = 0 if current < Grid_Width, etc. 
    const size_t column  = (size_t)current%(size_t)Grid_Width;
    
    // Determine the positions neigboring the current position
    const size_t left    =  row      * Grid_Width + column - 1;
    const size_t right   =  row      * Grid_Width + column + 1;
    
    return ( ( (h[right]+s[right]+b[right])
             - (h[left ]+s[left ]+b[left ]) ) /2.0/dx );
}

// (OK.)

float dO_dy(__global float* h, __global float* s,__global float* b)
{
    const float dy = dY;  // Forcing dY to become a float
    
    // Determine the position at which the current thread is computing
    const size_t current = get_global_id(0);            
    const size_t row     = (size_t)current/(size_t)Grid_Width; // Rounds down to integer, so row = 0 if current < Grid_Width, etc. 
    const size_t column  = (size_t)current%(size_t)Grid_Width;
    
    // Determine the positions neigboring the current position
    const size_t top     = (row - 1) * Grid_Width + column    ;
    const size_t bottom  = (row + 1) * Grid_Width + column    ;
    
    return ( ( (h[bottom]+s[bottom]+b[bottom])
             - (h[top   ]+s[top   ]+b[top   ]) ) /2.0/dy );
}

// (OK.)

////////////////////////////////////////////////////////////////////////////////
// Definition of the functions that compute the derivatives of uh and vh
////////////////////////////////////////////////////////////////////////////////

float d_uh_dx(__global float* u, __global float* h)
{
    const float dx = dX;  // Forcing dX to become a float
    
    // Determine the position at which the current thread is computing
    const size_t current = get_global_id(0);            
    const size_t row     = (size_t)current/(size_t)Grid_Width; // Rounds down to integer, so row = 0 if current < Grid_Width, etc. 
    const size_t column  = (size_t)current%(size_t)Grid_Width;
    
    // Determine the positions neigboring the current position
    const size_t left    =  row      * Grid_Width + column - 1;
    const size_t right   =  row      * Grid_Width + column + 1;
    
    return ( ( u[right]*h[right] - u[left]*h[left] )/2.0/dx );
}

// (OK.)

float d_vh_dy(__global float* v, __global float* h)
{
    const float dy = dY;  // Forcing dY to become a float
    
    // Determine the position at which the current thread is computing
    const size_t current = get_global_id(0);            
    const size_t row     = (size_t)current/(size_t)Grid_Width; // Rounds down to integer, so row = 0 if current < Grid_Width, etc. 
    const size_t column  = (size_t)current%(size_t)Grid_Width;
    
    // Determine the positions neigboring the current position
    const size_t top     = (row - 1) * Grid_Width + column    ;
    const size_t bottom  = (row + 1) * Grid_Width + column    ;
    
    return ( ( v[bottom]*h[bottom] - v[top]*h[top] )/2.0/dy );
}

// (OK.)

////////////////////////////////////////////////////////////////////////////////
// Periodic Boundary conditions function
////////////////////////////////////////////////////////////////////////////////

void PeriodicBoundaries(__global float* z)
{
    // Determine the position at which the current thread is computing
    const size_t current = get_global_id(0);            
    const size_t row     = (size_t)current/(size_t)Grid_Width; // Rounds down to integer, so row = 0 if current < Grid_Width, etc. 
    const size_t column  = (size_t)current%(size_t)Grid_Width;
    
    if(row==0) // Top boundary
    {
        z[current] = z[ (Grid_Height-2) * Grid_Width + column ];
    }
    else if(row==Grid_Height-1) // Bottom boundary
    {
        z[current] = z[ 1 * Grid_Width + column ];
    }
    else if(column==0) // Left boundary
    {
        z[current] = z[ (row+1) * Grid_Width - 2 ];
    }
    else if(column==Grid_Width-1) // Right boundary
    {
        z[current] = z[ row * Grid_Width + 1 ];
    }
}

// (OK.)

////////////////////////////////////////////////////////////////////////////////
// Neumann Boundary conditions function, having zero flux on the edge
////////////////////////////////////////////////////////////////////////////////

void NeumanBoundaries(__global float* z)
{
    // Determine the position at which the current thread is computing
    const size_t current = get_global_id(0);            
    const size_t row     = (size_t)current/(size_t)Grid_Width; // Rounds down to integer, so row = 0 if current < Grid_Width, etc. 
    const size_t column  = (size_t)current%(size_t)Grid_Width;
    
    if(row==0) // Top boundary
    {
        z[ row * Grid_Width + column ] = z[ 1 * Grid_Width + column ];
    }
    else if(row==Grid_Height-1) // Bottom boundary
    {
        z[ row * Grid_Width + column ] = z[ (Grid_Height-2) * Grid_Width + column ];
    }
    else if(column==0) // Left boundary
    {
        z[ row * Grid_Width + column ] = z[ row * Grid_Width + 1 ];
    }
    else if(column==Grid_Width-1) // Right boundary
    {
        z[ row * Grid_Width + column ] = z[ row * Grid_Width + Grid_Width - 2 ];
    }
}

// (OK.)

////////////////////////////////////////////////////////////////////////////////
// Reflecting Boundary conditions function, for shallow water equations
////////////////////////////////////////////////////////////////////////////////

void ReflectingBoundaries(__global float* u,__global float* v)
{
    // Determine the position at which the current thread is computing
    const size_t current = get_global_id(0);            
    const size_t row     = (size_t)current/(size_t)Grid_Width; // Rounds down to integer, so row = 0 if current < Grid_Width, etc. 
    const size_t column  = (size_t)current%(size_t)Grid_Width;
    
    if(row==0) // Top boundary
    {
        u[row * Grid_Width + column] = u[1*Grid_Width+column];
        v[row * Grid_Width + column] =-v[1*Grid_Width+column];
    }
    else if(row==Grid_Height-1) // Bottom boundary
    {
        u[row * Grid_Width + column] = u[(Grid_Height-2) * Grid_Width+column];
        v[row * Grid_Width + column] =-v[(Grid_Height-2) * Grid_Width+column];
    }
    else if(column==0) // Left boundary
    {
        u[row * Grid_Width + column] =-u[row * Grid_Width + 1];
        v[row * Grid_Width + column] = v[row * Grid_Width + 1];
    }
    else if(column==Grid_Width-1) // Right boundaries
    {
        u[row * Grid_Width + column] =-u[row * Grid_Width + Grid_Width-2];
        v[row * Grid_Width + column] = v[row * Grid_Width + Grid_Width-2];
    }
}

// (OK.)

////////////////////////////////////////////////////////////////////////////////
// Persistent Flux Boundary condition function, extrapolating over the boundaries
////////////////////////////////////////////////////////////////////////////////

void PersistentFluxBoundaries(__global float* z)
{
    // Determine the position at which the current thread is computing
    const size_t current = get_global_id(0);            
    const size_t row     = (size_t)current/(size_t)Grid_Width; // Rounds down to integer, so row = 0 if current < Grid_Width, etc. 
    const size_t column  = (size_t)current%(size_t)Grid_Width;
    
    if(row==0) // Top boundary
    {
        z[current] = 2.0*z[(row+1) * Grid_Width + column] - z[ (row+2) * Grid_Width + column];
    }
    else if(row==Grid_Height-1) // Bottom boundary
    {
        z[current] = 2.0*z[(row-1) * Grid_Width + column] - z[(row-2) * Grid_Width + column];
    }
    else if(column==0) // Left boundary
    {
        z[current] = 2.0*z[row * Grid_Width + column+1] - z[row * Grid_Width + column+2];
    }
    else if(column==Grid_Width-1) // Right boundary
    {
        z[current] = 2.0*z[row * Grid_Width + column-1] - z[row * Grid_Width + column-2];
    }
}

// (OK.)

////////////////////////////////////////////////////////////////////////////////
// Dirichlet Boundary condition function, having fixed values on the edge
////////////////////////////////////////////////////////////////////////////////

void DirichletBoundaries(__global float* z, float Value)
{
    // Determine the position at which the current thread is computing
    const size_t current = get_global_id(0);            
    const size_t row     = (size_t)current/(size_t)Grid_Width; // Rounds down to integer, so row = 0 if current < Grid_Width, etc. 
    const size_t column  = (size_t)current%(size_t)Grid_Width;
    
    if(row==0) // Top boundary
    {
        z[current]=Value;
    }
    else if(row==Grid_Height-1) // Bottom boundary
    {
        z[current]=Value;
    }
    else if(column==0) // Left boundary
    {
        z[current]=Value;
    }
    else if(column==Grid_Width-1) // Right boundary
    {
        z[current]=Value;
    }
}

// (OK.)

#endif