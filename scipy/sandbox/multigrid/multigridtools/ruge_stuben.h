#ifndef RUGE_STUBEN_H
#define RUGE_STUBEN_H

#include <iostream>
#include <vector>
#include <iterator>
#include <assert.h>


//this will increase the complexity greatly!
//#define DEBUG

enum NodeType {U_NODE,  C_NODE, F_NODE};


template<class T>
void rs_strong_connections(const int n_row,
			   const T theta,
			   const int Ap[], const int Aj[], const T Ax[],
			   std::vector<int> * Sp, std::vector<int> * Sj, std::vector<T> * Sx){
    //Sp,Sj form a CSR representation where the i-th row contains
    //the indices of all the strong connections from node i
    Sp->push_back(0);

    //Compute lambdas for each node
    for(int i = 0; i < n_row; i++){
        T min_offdiagonal = 0.0;

        int row_start = Ap[i];
        int row_end   = Ap[i+1];
        for(int jj = row_start; jj < row_end; jj++){
            min_offdiagonal = std::min(min_offdiagonal,Ax[jj]); //assumes diagonal is positive!
        }

        T threshold = theta*min_offdiagonal;
        for(int jj = row_start; jj < row_end; jj++){
            if(Ax[jj] < threshold){
	            Sj->push_back(Aj[jj]);
	            Sx->push_back(Ax[jj]);
            }
        }

        Sp->push_back(Sj->size());
    }
}




template<class T>
void rs_interpolation(const int n_nodes,
		      const int Ap[], const int Aj[], const T Ax[],
		      const int Sp[], const int Sj[], const T Sx[],
		      const int Tp[], const int Tj[], const T Tx[],
		      std::vector<int> * Bp, std::vector<int> * Bj, std::vector<T> * Bx){
  
  std::vector<int> lambda(n_nodes,0);

  //compute lambdas
  for(int i = 0; i < n_nodes; i++){
    lambda[i] = Tp[i+1] - Tp[i];
  }


  //for each value of lambda, create an interval of nodes with that value
  // ptr - is the first index of the interval
  // count - is the number of indices in that interval
  // index to node - the node located at a given index
  // node to index - the index of a given node
  std::vector<int> interval_ptr(n_nodes,0);
  std::vector<int> interval_count(n_nodes,0);
  std::vector<int> index_to_node(n_nodes);
  std::vector<int> node_to_index(n_nodes);

  for(int i = 0; i < n_nodes; i++){
    interval_count[lambda[i]]++;
  }
  for(int i = 0, cumsum = 0; i < n_nodes; i++){
    interval_ptr[i] = cumsum;
    cumsum += interval_count[i];
    interval_count[i] = 0;
  }
  for(int i = 0; i < n_nodes; i++){
    int lambda_i = lambda[i];
    int index    = interval_ptr[lambda_i]+interval_count[lambda_i];
    index_to_node[index] = i;
    node_to_index[i]     = index;
    interval_count[lambda_i]++;
  }



  

  std::vector<NodeType> NodeSets(n_nodes,U_NODE);

  //Now add elements to C and F, in decending order of lambda
  for(int top_index = n_nodes - 1; top_index > -1; top_index--){
    int i        = index_to_node[top_index];
    int lambda_i = lambda[i];
#ifdef DEBUG
    {
#ifdef DEBUG_PRINT
      std::cout << "top_index " << top_index << std::endl;
      std::cout << "i         " << i << std::endl;
      std::cout << "lambda_i  " << lambda_i << std::endl;

      for(int i = 0; i < n_nodes; i++){
	std::cout << i << "=";
	if(NodeSets[i] == U_NODE)
	  std::cout << "U";
	else if(NodeSets[i] == F_NODE)
	  std::cout << "F";
	else
	  std::cout << "C";
	std::cout << " ";
      }
      std::cout << std::endl;

      std::cout << "node_to_index" << std::endl;
      for(int i = 0; i < n_nodes; i++){
	std::cout << i << "->" << node_to_index[i] << "  ";
      }
      std::cout << std::endl;
      std::cout << "index_to_node" << std::endl;
      for(int i = 0; i < n_nodes; i++){
	std::cout << i << "->" << index_to_node[i] << "  ";
      }
      std::cout << std::endl;

      std::cout << "interval_count ";
      for(int i = 0; i < n_nodes; i++){
	std::cout << interval_count[i] << " ";
      }
      std::cout << std::endl;
#endif

      //make sure arrays are correct
      for(int n = 0; n < n_nodes; n++){
	assert(index_to_node[node_to_index[n]] == n);
      }

      //make sure intervals are reasonable
      int sum_intervals = 0;
      for(int n = 0; n < n_nodes; n++){
	assert(interval_count[n] >= 0);
	if(interval_count[n] > 0){
	  assert(interval_ptr[n] == sum_intervals);
	}
	sum_intervals += interval_count[n];
      }
      assert(sum_intervals == top_index+1);

      
      if(interval_count[lambda_i] <= 0){
	std::cout << "top_index " << top_index << std::endl;
	std::cout << "lambda_i " << lambda_i << std::endl;
	std::cout << "interval_count[lambda_i] " << interval_count[lambda_i] << std::endl;
	std::cout << "top_index " << top_index << std::endl;
	std::cout << "i         " << i << std::endl;
	std::cout << "lambda_i  " << lambda_i << std::endl;
      }
      
      
      for(int n = 0; n <= top_index; n++){
	assert(NodeSets[index_to_node[n]] != C_NODE);
      }
    }
    assert(node_to_index[i] == top_index);
    assert(interval_ptr[lambda_i] + interval_count[lambda_i] - 1 == top_index);
    //max interval should have at least one element
    assert(interval_count[lambda_i] > 0);    
#endif


    //remove i from its interval
    interval_count[lambda_i]--;
    

    if(NodeSets[i] == F_NODE){
      continue;
    } else {
      assert(NodeSets[i] == U_NODE);

      NodeSets[i] = C_NODE;

      //For each j in S^T_i /\ U
      for(int jj = Tp[i]; jj < Tp[i+1]; jj++){
	int j = Tj[jj];

	if(NodeSets[j] == U_NODE){
	  NodeSets[j] = F_NODE;
	  
	  //For each k in S_j /\ U
	  for(int kk = Sp[j]; kk < Sp[j+1]; kk++){
	    int k = Sj[kk];

	    if(NodeSets[k] == U_NODE){	      
	      //move k to the end of its current interval
	      assert(lambda[j] < n_nodes - 1);//this would cause problems!

	      int lambda_k = lambda[k];
	      int old_pos  = node_to_index[k];
	      int new_pos  = interval_ptr[lambda_k] + interval_count[lambda_k] - 1;

	      node_to_index[index_to_node[old_pos]] = new_pos;
	      node_to_index[index_to_node[new_pos]] = old_pos;
	      std::swap(index_to_node[old_pos],index_to_node[new_pos]);
	      
	      //update intervals
	      interval_count[lambda_k]   -= 1;
	      interval_count[lambda_k+1] += 1;
	      interval_ptr[lambda_k+1]    = new_pos;

	      //increment lambda_k
	      lambda[k]++;

#ifdef DEBUG
	      assert(interval_count[lambda_k]   >= 0);
	      assert(interval_count[lambda_k+1] >  0);
	      assert(interval_ptr[lambda[k]] <= node_to_index[k]);
	      assert(node_to_index[k] < interval_ptr[lambda[k]] + interval_count[lambda[k]]);
#endif
	    }
	  }
	}
      }

      //For each j in S_i /\ U
      for(int jj = Sp[i]; jj < Sp[i+1]; jj++){
	int j = Sj[jj];
	if(NodeSets[j] == U_NODE){            //decrement lambda for node j
	  assert(lambda[j] > 0);//this would cause problems!

	  //move j to the beginning of its current interval
	  int lambda_j = lambda[j];
	  int old_pos  = node_to_index[j];
	  int new_pos  = interval_ptr[lambda_j]; 
	      
	  node_to_index[index_to_node[old_pos]] = new_pos;
	  node_to_index[index_to_node[new_pos]] = old_pos;
	  std::swap(index_to_node[old_pos],index_to_node[new_pos]);
	      
	  //update intervals
	  interval_count[lambda_j]   -= 1;
	  interval_count[lambda_j-1] += 1;
	  interval_ptr[lambda_j]     += 1;
	  interval_ptr[lambda_j-1]    = interval_ptr[lambda_j] - interval_count[lambda_j-1];

	  //decrement lambda_j
	  lambda[j]--;

#ifdef DEBUG
	  assert(interval_count[lambda_j]   >= 0);
	  assert(interval_count[lambda_j-1] >  0);
	  assert(interval_ptr[lambda[j]] <= node_to_index[j]);
	  assert(node_to_index[j] < interval_ptr[lambda[j]] + interval_count[lambda[j]]);
#endif
	}
      }
    }
  }




#ifdef DEBUG
  //make sure each f-node has at least one strong c-node neighbor
  for(int i = 0; i < n_nodes; i++){
    if(NodeSets[i] == F_NODE){
      int row_start = Sp[i];
      int row_end   = Sp[i+1];
      bool has_c_neighbor = false;
      for(int jj = row_start; jj < row_end; jj++){
	if(NodeSets[Sj[jj]] == C_NODE){
	  has_c_neighbor = true;
	  break;
	}
      }
      assert(has_c_neighbor);
    }   
  }
#endif

  //Now construct interpolation operator
  std::vector<T> d_k(n_nodes,0);
  std::vector<bool> C_i(n_nodes,0);
  Bp->push_back(0);
  for(int i = 0; i < n_nodes; i++){
    if(NodeSets[i] == C_NODE){
      //interpolate directly
      Bj->push_back(i);
      Bx->push_back(1);      
      Bp->push_back(Bj->size());
    } else {
      //F_NODE
      
      //Step 4
      T d_i = 0; //denominator for this row
      for(int jj = Ap[i]; jj < Ap[i+1]; jj++){ d_i += Ax[jj]; }
      for(int jj = Sp[i]; jj < Sp[i+1]; jj++){ d_i -= Sx[jj]; }
      
      //Create C_i, initialize d_k
      for(int jj = Sp[i]; jj < Sp[i+1]; jj++){ 
	int j = Sj[jj];
	if(NodeSets[j] == C_NODE){
	  C_i[j] = true;
	  d_k[j] = Sx[jj];
	}
      }

      bool Sj_intersects_Ci = true; //in the case that i has no F-neighbors
      for(int jj = Sp[i]; jj < Sp[i+1]; jj++){ //for j in D^s_i
	int    j = Sj[jj];
	T   a_ij = Sx[jj];
	T   a_jl = 0;

	if(NodeSets[j] != F_NODE){continue;}

	//Step 5
	Sj_intersects_Ci = false;

	//compute sum a_jl
	for(int ll = Sp[j]; ll < Sp[j+1]; ll++){
	  if(C_i[Sj[ll]]){
	    Sj_intersects_Ci = true;
	    a_jl += Sx[ll];
	  }	    
	}

	if(!Sj_intersects_Ci){ break; }

	for(int kk = Sp[j]; kk < Sp[j+1]; kk++){
	  int   k = Sj[kk];
	  T  a_jk = Sx[kk];
	  if(C_i[k]){
	    d_k[k] += a_ij*a_jk / a_jl;
	  }	    
	}
      }

      //Step 6
      if(Sj_intersects_Ci){
	for(int jj = Sp[i]; jj < Sp[i+1]; jj++){ 
	  int j = Sj[jj];
	  if(NodeSets[j] == C_NODE){
	    Bj->push_back(j);
	    Bx->push_back(-d_k[j]/d_i);      
	  }
	}	
	Bp->push_back(Bj->size());
      } else { //make i a C_NODE
	NodeSets[i] = C_NODE;
	Bj->push_back(i);
	Bx->push_back(1);      
	Bp->push_back(Bj->size());
      }
      

      //Clear C_i,d_k
      for(int jj = Sp[i]; jj < Sp[i+1]; jj++){ 
	int j = Sj[jj];
	C_i[j] = false;
	d_k[j] = 0;
      }

    }
      
  }

  //for each c-node, determine its index in the coarser lvl
  std::vector<int> cnode_index(n_nodes,-1);
  int n_cnodes = 0;
  for(int i = 0; i < n_nodes; i++){
    if(NodeSets[i] == C_NODE)
      cnode_index[i] = n_cnodes++;
  }
  //map old C indices to coarse indices
  for(std::vector<int>::iterator i = Bj->begin(); i != Bj->end(); i++){
    *i = cnode_index[*i];
  }

  
  
}

#endif
