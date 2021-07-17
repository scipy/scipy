#pythran export sort_vertices_of_regions(int[:,:], int list list)
#pythran export sort_vertices_of_regions(intc[:,:], intc list list)
#pythran export sort_vertices_of_regions(int32[:,:], int64 list list)
#pythran export sort_vertices_of_regions(int64[:,:], int64 list list)
def sort_vertices_of_regions(simplices, regions):
    for n in range(0, len(regions)):
        remaining = regions[n][:]
        sorted_vertices = []
        current_simplex = remaining[0]
        for k in simplices[current_simplex]:
            if k != n:
                current_vertex = k
                break
        remaining.remove(current_simplex)
        sorted_vertices.append(current_simplex)
        while remaining:
            for s in remaining:
                if current_vertex in simplices[s]:
                    current_simplex = s
                    break
            for s in simplices[current_simplex]:
                if s != n and s != current_vertex:
                    current_vertex = s
                    break
            remaining.remove(current_simplex)
            sorted_vertices.append(current_simplex)
        regions[n] = list(sorted_vertices)
    
    return regions