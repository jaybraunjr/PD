import numpy as np

def _dfs(graph, start):
    visited, stack = set(), [start]
    while stack:
        vertex = stack.pop()
        if vertex not in visited:
            visited.add(vertex)
            stack.extend(graph[vertex] - visited)
    return visited



def _make_graph(matrix):

    graph = {}
    xis, yis = matrix.shape  # Get the dimensions of the matrix

    for (xi, yi), value in np.ndenumerate(matrix):
        if value == 0:
            continue  
            
        node_index = xi * yis + yi
        neighbor_list = []  
        for dx in [-1, 0, 1]:
            for dy in [-1, 0, 1]:

                x, y = divmod(xi + dx, xis)[1], divmod(yi + dy, yis)[1]

                if matrix[x, y] == 1 and (x, y) != (xi, yi):
                    neighbor_node_index = x * yis + y  
                    neighbor_list.append(neighbor_node_index)  

        graph[node_index] = set(neighbor_list)  

    return graph