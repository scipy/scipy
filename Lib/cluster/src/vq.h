#ifndef vq_h
#define vq_h
/*
//#define real  float 
//#define scan_format "%f"

#define real double
#define scan_format "%lf"



void vq_obs(real* obs,real* code_book, int Ncodes, int Nfeatures,
			   int& code, real& lowest_dist);

void vq(real* obs,real* code_book, int Nobs, int Ncodes, int Nfeatures,
	    int* codes, real* lowest_dist);
*/
#define BIG 10000.

template<class T>
void tvq_obs(T* obs,T* code_book, int Ncodes, int Nfeatures,
			   int& code, T& lowest_dist)
{
	int i,j,k=0;
	T dist, diff;

	lowest_dist = (T) BIG;
	for(i=0; i < Ncodes; i++)
	{
		dist=0;
		for(j=0; j < Nfeatures; j++)
		{
			diff = code_book[k] - obs[j];
			dist += diff*diff;
			k++;
		}
		dist = (T)sqrt(dist);
		if (dist < lowest_dist)
		{
			code = i;
			lowest_dist = dist;
		}
	}
}

template<class T>
void tvq(T* obs,T* code_book, int Nobs, int Ncodes, int Nfeatures,
	    int* codes, T* lowest_dist)
{
    int i;
	for( i = 0; i < Nobs; i++)
	{		
		tvq_obs<T>(&(obs[i*Nfeatures]),code_book,Ncodes,Nfeatures,
				  codes[i],lowest_dist[i]);
	}
}
#endif