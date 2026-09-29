import os
import sys
import numpy as np


# ---------------------------------------------------------------- ESM
def local_ESM(E, A, L):
    """Local 2x2 axial stiffness matrix and its 4x4 form (dofs: u1, v1, u2, v2)."""
    if L <= 0:
        raise ValueError(f"Element length must be positive, got L = {L}")

    k = (E * A) / L
    K_local = np.array([[k, -k],
                        [-k, k]])

    k4 = np.zeros((4, 4))
    k4[np.ix_([0, 2], [0, 2])] = K_local   # axial dofs only
    return K_local, k4


# ---------------------------------------------------------------- Rotation
def Rotation(x1, x2, y1, y2):
    """4x4 local->global rotation matrix (block diagonal) and angle theta."""
    theta = np.arctan2(y2 - y1, x2 - x1)
    c, s = np.cos(theta), np.sin(theta)
    R2 = np.array([[c, -s],
                   [s,  c]])
    R = np.zeros((4, 4))
    R[:2, :2] = R2
    R[2:, 2:] = R2
    return R, theta


# ---------------------------------------------------------------- GSM
def assemble_truss_stiffness(nodes, elements, E, A, dof_node=2):
    """
    nodes    : (n_nodes, 2) float  -> [x, y]
    elements : (n_elem, 2)  int    -> 0-based node indices [i, j]
    E, A     : (n_elem,)           -> per-element properties
    """
    n_dof = dof_node * len(nodes)
    K = np.zeros((n_dof, n_dof))

    for e, (i, j) in enumerate(elements):
        x1, y1 = nodes[i]
        x2, y2 = nodes[j]
        L = np.hypot(x2 - x1, y2 - y1)

        _, k4 = local_ESM(E[e], A[e], L)
        R, _ = Rotation(x1, x2, y1, y2)
        k_global = R @ k4 @ R.T            # local -> global

        dofs = np.array([dof_node * i, dof_node * i + 1,
                         dof_node * j, dof_node * j + 1])
        K[np.ix_(dofs, dofs)] += k_global  # scatter into global matrix
    return K


# ---------------------------------------------------------------- Solve
def solve_gsm(K, F, fixed_dofs):
    """Solve K u = F on the free dofs; return displacements and reactions."""
    n_dof = K.shape[0]
    fixed_dofs = np.asarray(fixed_dofs, dtype=int)
    free_dofs = np.setdiff1d(np.arange(n_dof), fixed_dofs)

    u = np.zeros(n_dof)
    u[free_dofs] = np.linalg.solve(K[np.ix_(free_dofs, free_dofs)], F[free_dofs])

    RF = K @ u - F
    return u, RF


# ---------------------------------------------------------------- Post-processing
def strain(U, nodes, elements):
    """Signed axial strain per element (positive = tension)."""
    eps = np.zeros(len(elements))
    for e, (i, j) in enumerate(elements):
        dx, dy = nodes[j] - nodes[i]
        L = np.hypot(dx, dy)
        c, s = dx / L, dy / L
        du = U[2 * j:2 * j + 2] - U[2 * i:2 * i + 2]
        eps[e] = (c * du[0] + s * du[1]) / L   # elongation along the bar axis
    return eps


def postprocess(U, F, K, nodes, elements, E):
    R = K @ U - F
    Strain = strain(U, nodes, elements)
    Stress = E * Strain
    return R, Stress, Strain


# ---------------------------------------------------------------- File input
def read_model(filepath):
    nodes, elements, bcs, loads = {}, {}, {}, {}
    info = ""

    with open(filepath, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line[0] == "[":
                info = line[1:-1].strip().upper()
                continue

            col = [c.strip() for c in line.split(",")]

            if info == "NODES":
                nodes[int(col[0])] = (float(col[1]), float(col[2]))
            elif info == "ELEMENTS":
                elements[int(col[0])] = (int(col[1]), int(col[2]),
                                         float(col[3]),   # E
                                         float(col[4]),   # poisson (unused in truss)
                                         float(col[5]))   # area
            elif info == "BCS":
                bcs[int(col[0])] = (int(col[1]), int(col[2]))
            elif info == "LOADS":
                loads[int(col[0])] = (float(col[1]), float(col[2]))

    return nodes, elements, bcs, loads


# ---------------------------------------------------------------- Main
def main(filepath):
    nodes_d, elements_d, bcs_d, loads_d = read_model(filepath)

    # Map file node ids (any numbering) to 0-based indices
    node_ids = sorted(nodes_d)
    idx = {nid: k for k, nid in enumerate(node_ids)}
    nodes = np.array([nodes_d[n] for n in node_ids], dtype=float)

    elem_ids = sorted(elements_d)
    elements = np.array([[idx[elements_d[e][0]], idx[elements_d[e][1]]] for e in elem_ids])
    E = np.array([elements_d[e][2] for e in elem_ids])
    A = np.array([elements_d[e][4] for e in elem_ids])

    K = assemble_truss_stiffness(nodes, elements, E, A, dof_node=2)

    F = np.zeros(2 * len(nodes))
    for nid, (fx, fy) in loads_d.items():
        F[2 * idx[nid]] += fx
        F[2 * idx[nid] + 1] += fy

    fixed = []
    for nid, (bx, by) in bcs_d.items():      # assumes 1 = fixed, 0 = free
        if bx:
            fixed.append(2 * idx[nid])
        if by:
            fixed.append(2 * idx[nid] + 1)

    U, _ = solve_gsm(K, F, fixed)
    R, Stress, Strain = postprocess(U, F, K, nodes, elements, E)

    print("Displacement vector U:\n", U)
    print("Reaction forces R:\n", R)
    print("Stress in each element:\n", Stress)
    print("Strain in each element:\n", Strain)


if _name_ == "_main_":
    # Resolve the model file relative to this script, not the terminal's cwd
    script_dir = os.path.dirname(os.path.abspath(_file_))
    name = sys.argv[1] if len(sys.argv) > 1 else "model 1.txt"
    path = name if os.path.isabs(name) else os.path.join(script_dir, name)
    if not os.path.exists(path):
        sys.exit(f"Model file not found: {path}")
    main(path)
