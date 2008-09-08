#include <time.h>
#include <math.h>
#include <iostream>
#include <algorithm>

template <class T>
void linear(T* x_vec, T* y_vec, int len,
            T* new_x_vec, T* new_y_vec, int new_len)
{    
    for (int i=0;i<new_len;i++)
    {
       T new_x = new_x_vec[i];
       int index;
       if (new_x <= x_vec[0])
           index = 0;
       else if (new_x >=x_vec[len-1])
           index = len-2;
       else
       { 
           T* which = std::lower_bound(x_vec, x_vec+len, new_x);
           index = which - x_vec-1;            
       }
             
       if(new_x == x_vec[index])
       {
           // exact value
           new_y_vec[i] = y_vec[index];
       }
       else
       {     
           //interpolate
           double x_lo = x_vec[index];
           double x_hi = x_vec[index+1];
           double y_lo = y_vec[index];
           double y_hi = y_vec[index+1];
           double slope = (y_hi-y_lo)/(x_hi-x_lo);
           new_y_vec[i] = slope * (new_x-x_lo) + y_lo;
       }
    }
}

template <class T>
void loginterp(T* x_vec, T* y_vec, int len,
               T* new_x_vec, T* new_y_vec, int new_len)
{    
    for (int i=0;i<new_len;i++)
    {
        T new_x = new_x_vec[i];
        int index;
        if (new_x <= x_vec[0])
            index = 0;
        else if (new_x >=x_vec[len-1])
            index = len-2;
        else
        { 
            T* which = std::lower_bound(x_vec, x_vec+len, new_x);
            index = which - x_vec-1;            
        }
             
        if(new_x == x_vec[index])
        {
            // exact value
            new_y_vec[i] = y_vec[index];
        }
        else
        {     
            //interpolate
            double x_lo = x_vec[index];
            double x_hi = x_vec[index+1];
            double y_lo = log10(y_vec[index]);
            double y_hi = log10(y_vec[index+1]);
            double slope = (y_hi-y_lo)/(x_hi-x_lo);
            new_y_vec[i] = pow(10.0, (slope * (new_x-x_lo) + y_lo));
        }
    }
}

template <class T>
int block_average_above(T* x_vec, T* y_vec, int len,
                         T* new_x_vec, T* new_y_vec, int new_len)
{    
    int bad_index = -1;
    int start_index = 0;
    T last_y = 0.0;
    T thickness = 0.0;

    for(int i=0;i<new_len;i++)
    {
        T new_x = new_x_vec[i];
        if ((new_x < x_vec[0]) || (new_x > x_vec[len-1]))
        {
            bad_index = i;
            break;
        }
        else if (new_x == x_vec[0])
        {
            // for the first sample, just return the cooresponding y value
            new_y_vec[i] = y_vec[0];             
        }        
        else
        {
            T* which = std::lower_bound(x_vec, x_vec+len, new_x);
            int index = which - x_vec-1;
            
            // calculate weighted average
            
            // Start off with "residue" from last interval in case last x
            // was between to samples.
            T weighted_y_sum = last_y * thickness;
            T thickness_sum = thickness;  
            for(int j=start_index; j<=index; j++)
            {
                    if (x_vec[j+1] < new_x)
                        thickness = x_vec[j+1] - x_vec[j];
                    else
                        thickness = new_x -x_vec[j];
                    weighted_y_sum += y_vec[j] * thickness;
                    thickness_sum += thickness;
            }       
            new_y_vec[i] = weighted_y_sum/thickness_sum; 
            
            // Store the thickness between the x value and the next sample
            // to add to the next weighted average.
            last_y = y_vec[index];
            thickness = x_vec[index+1] - new_x;
            
            // start next weighted average at next sample
            start_index =index+1;
        }
    }
    return bad_index;
}

template <class T>
int window_average(T* x_vec, T* y_vec, int len,
                   T* new_x_vec, T* new_y_vec, int new_len,
                   T width)
{    
    for(int i=0;i<new_len;i++)
    {
        T new_x = new_x_vec[i];
        T bottom = new_x - width/2;
        T top = new_x + width/2;
            
        T* which = std::lower_bound(x_vec, x_vec+len, bottom);
        int bottom_index = which - x_vec;
        if (bottom_index < 0)
        {
            //bottom = x_vec[0];
            bottom_index = 0;
        }
        
        which = std::lower_bound(x_vec, x_vec+len, top);
        int top_index = which - x_vec;
        if (top_index >= len)
        {
            //top = x_vec[len-1];
            top_index = len-1;
        }
        //std::cout << std::endl;
        //std::cout << bottom_index << " " << top_index << std::endl;
        //std::cout << bottom << " " << top << std::endl;
        // calculate weighted average
        T thickness =0.0;
        T thickness_sum =0.0;
        T weighted_y_sum =0.0;
        for(int j=bottom_index; j < top_index; j++)
        {
            thickness = x_vec[j+1] - bottom;
            weighted_y_sum += y_vec[j] * thickness;
            thickness_sum += thickness;
            bottom = x_vec[j+1];
            /*
            std::cout <<  "iter: " << j - bottom_index << " " <<
                     "index: " << j << " " <<
                     "bottom: " << bottom << " " <<
                     "x+1: " << x_vec[j+1] << " " <<
                     "x: " << x_vec[j] << " " << 
                     "y: " << y_vec[j] << " " <<
                     "weighted_sum: " << weighted_y_sum << 
                     "thickness: " << thickness << " " <<
                     "thickness_sum: " << thickness_sum << std::endl;
            */
            //std::cout << x_vec[j] << " ";         
            //std::cout << thickness << " ";         
        }
       
        // last element
        thickness = top - bottom;
        weighted_y_sum += y_vec[top_index] * thickness;
        thickness_sum += thickness;
            /*
            std::cout <<  "iter: last" << " " <<
                     "index: " << top_index << " " <<
                     "x: " << x_vec[top_index] << " " << 
                     "y: " << y_vec[top_index] << " " <<
                     "weighted_sum: " << weighted_y_sum << 
                     "thickness: " << thickness << " " <<
                     "thickness_sum: " << thickness_sum << std::endl;
            */
            //std::cout << x_vec[top_index] << " " <<  thickness_sum << std::endl;   
        new_y_vec[i] = weighted_y_sum/thickness_sum; 
    }
    return -1;
}
