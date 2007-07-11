#ifndef SMOOTHED_AGGREGATION_H
#define SMOOTHED_AGGREGATION_H

#include <iostream>
#include <vector>
#include <iterator>
#include <assert.h>


#define DEBUG


template<class T>
void sa_strong_connections(const int n_row,
			   const T epsilon,
			   const int Ap[], const int Aj[], const T Ax[],
			   std::vector<int> * Sp, std::vector<int> * Sj, std::vector<T> * Sx){
  //Sp,Sj form a CSR representation where the i-th row contains
  //the indices of all the strong connections from node i
  Sp->push_back(0);

  //compute diagonal values
  std::vector<T> diags(n_row);
  for(int i = 0; i < n_row; i++){
    int row_start = Ap[i];
    int row_end   = Ap[i+1];
    for(int jj = row_start; jj < row_end; jj++){
      if(Aj[jj] == i){
	diags[i] = Ax[jj];
	break;
      }
    }    
  }

#ifdef DEBUG
  for(int i = 0; i < n_row; i++){ assert(diags[i] > 0); }
#endif
    


  for(int i = 0; i < n_row; i++){
    int row_start = Ap[i];
    int row_end   = Ap[i+1];

    T eps_Aii = epsilon*epsilon*diags[i];

    for(int jj = row_start; jj < row_end; jj++){
      const int&   j = Aj[jj];
      const T&   Aij = Ax[jj];

      if(i == j){continue;}

      if(Aij*Aij >= eps_Aii * diags[j]){
	Sj->push_back(j);
	Sx->push_back(Aij);
      }
    }
    Sp->push_back(Sj->size());
  }
}


void sa_get_aggregates(const int n_row,
		       const int Ap[], const int Aj[],
		       std::vector<int> * Bj){

  std::vector<int> aggregates(n_row,-1);

  int num_aggregates = 0;

  //Pass #1
  for(int i = 0; i < n_row; i++){
    if(aggregates[i] >= 0){ continue; } //already marked

    const int& row_start = Ap[i];
    const int& row_end   = Ap[i+1];
    

    //Determine whether all neighbors of this node are free (not already aggregates)
    bool free_neighborhood = true;
    for(int jj = row_start; jj < row_end; jj++){
      if(aggregates[Aj[jj]] >= 0){
	free_neighborhood = false;
	break;
      }
    }    
    if(!free_neighborhood){ continue; } //bail out


    //Make an aggregate out of this node and its strong neigbors
    aggregates[i] = num_aggregates;
    for(int jj = row_start; jj < row_end; jj++){
      aggregates[Aj[jj]] = num_aggregates;
    }
    num_aggregates++;
  }



  //Pass #2
  std::vector<int> aggregates_copy(aggregates);
  for(int i = 0; i < n_row; i++){
    if(aggregates[i] >= 0){ continue; } //already marked

    const int& row_start = Ap[i];
    const int& row_end   = Ap[i+1];
    
    for(int jj = row_start; jj < row_end; jj++){
      const int& j = Aj[jj];

      if(aggregates_copy[j] >= 0){
	aggregates[i] = aggregates_copy[j];
	break;
      }
    }    
  }



  //Pass #3
  for(int i = 0; i < n_row; i++){
    if(aggregates[i] >= 0){ continue; } //already marked

    const int& row_start = Ap[i];
    const int& row_end   = Ap[i+1];
    
    aggregates[i] = num_aggregates;

    for(int jj = row_start; jj < row_end; jj++){
      const int& j = Aj[jj];

      if(aggregates[j] < 0){ //unmarked neighbors
	aggregates[j] = num_aggregates;
      }
    }  
    num_aggregates++;
  }


#ifdef DEBUG
  for(int i = 0; i < n_row; i++){ assert(aggregates[i] >= 0 && aggregates[i] < num_aggregates); }
#endif
  
  *Bj = aggregates;  
}






template<class T>
void sa_smoother(const int n_row,
		 const T   omega,
		 const int Ap[], const int Aj[], const T Ax[],
		 const int Sp[], const int Sj[], const T Sx[],
		 std::vector<int> * Bp, std::vector<int> * Bj, std::vector<T> * Bx){


  //compute filtered diagonal
  std::vector<T> diags(n_row,0);
  
  for(int i = 0; i < n_row; i++){
    int row_start = Ap[i];
    int row_end   = Ap[i+1];
    for(int jj = row_start; jj < row_end; jj++){
      diags[i] += Ax[jj];
    }    
  }
  for(int i = 0; i < n_row; i++){
    int row_start = Sp[i];
    int row_end   = Sp[i+1];
    for(int jj = row_start; jj < row_end; jj++){
      diags[i] -= Sx[jj];
    }    
  }
  
#ifdef DEBUG
  for(int i = 0; i < n_row; i++){ assert(diags[i] > 0); }
#endif


  //compute omega Jacobi smoother
  Bp->push_back(0);
  for(int i = 0; i < n_row; i++){
    int row_start = Sp[i];
    int row_end   = Sp[i+1];
    const T row_scale = -omega/diags[i];

    Bx->push_back(1.0);
    Bj->push_back( i );
    
    for(int jj = row_start; jj < row_end; jj++){
      Bx->push_back(row_scale*Sx[jj]);
      Bj->push_back(Sj[jj]);
    }    
    Bp->push_back(Bj->size());
  }
}



#endif
