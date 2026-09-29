def assemble_truss_stiffness(nodes, elements, E, A, dof_node=2):
    """
    nodes    : array-like (n_nodes, 2)        -> [x, y] per node
    elements : array-like (n_elements, 2) int -> 0-indexed [node_i, node_j]
    E, A     : array-like (n_elements,)       -> per-element properties
    """
    K_global_list = []

    for e, (i, j) in enumerate(elements):
        xi, yi = nodes[i]
        xj, yj = nodes[j]
        dx, dy = xj - xi, yj - yi
        L = np.hypot(dx, dy)
        c, s = dx / L, dy / L

        K_local, _ = local_ESM(E[e], A[e], L)                # 2x2, local axial
        K_global_list.append(K_local)

    K = GSM_assembly(elements, nodes, dof_node, K_global_list)
    return np.asarray(K, dtype=float)
