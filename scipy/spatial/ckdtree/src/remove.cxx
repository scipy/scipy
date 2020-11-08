#include <vector>
#include "ckdtree_decl.h"
#include <math.h>

int
find(ckdtree *self, ckdtreenode *actual_node, const double *point, std::vector<ckdtreenode *> *node_path){
   while(actual_node->split_dim != -1){
        if(*(point + actual_node->split_dim) < actual_node->split)
            actual_node = actual_node->less;
        else
            actual_node = actual_node->greater;

        node_path->push_back(actual_node);
    }
    //now actual_node is the node that can contains point

    for(int idx_indices=actual_node->start_idx; idx_indices<actual_node->end_idx; idx_indices++){
        bool equal = true;
        int idx_data = *(self->raw_indices + idx_indices);

        for(int point_dimension=0; point_dimension<self->m; point_dimension++){
            double *value_of_p = self->raw_data + idx_data*self->m + point_dimension;
            if(*value_of_p != *(point+point_dimension)){
                equal = false;
                break;
            }
        }

        if(equal)
            return idx_indices;
    }

    return -1;
}


int
remove(ckdtree *self, const double *point){
    std::vector<ckdtreenode *> node_path;
    node_path.push_back(self->ctree);
    int found_idx = find(self, self->ctree, point, &node_path);

    if(found_idx != -1){ //there is a point to remove
        self->n--;

        for(int i=found_idx; i<self->n; i++)
            *(self->raw_indices+i) = *(self->raw_indices+i+1);

        std::vector<ckdtreenode *>::iterator path_it = node_path.begin();
        for(;path_it != node_path.end(); ++path_it){
            (*path_it)->children--;

            if((*path_it)->children <= self->leafsize){ //subtree removal
                (*path_it)->split_dim = -1;
                (*path_it)->split = 0;
                (*path_it)->less = 0;
                (*path_it)->greater = 0;
                self->size -= pow(2, node_path.end()-path_it) - 2;
                break;
            }
        }

        std::vector<ckdtreenode>::iterator it = self->tree_buffer->begin();
        for(;it!=self->tree_buffer->end(); ++it){
            it->start_idx -= (it->start_idx > found_idx);
            it->end_idx -= (it->end_idx > found_idx);
        }
    }

    return found_idx;
}