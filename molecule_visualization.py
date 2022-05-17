from pymolgen.molecule import Molecule
from pymolgen.molecule_formats import molecule_to_rdkit
import networkx

def plot_molecule(molecule, timeout: float = None, title="Molecule"):
    """
    Show a plot of this molecule.
    """
    from rdkit.Chem import Draw
    import multiprocessing

    def plot_on_thread():
        to_plot = molecule.copy()
        to_plot.hydrogenate()  # Helps with rdkit complaining about unkekulized atoms
        Draw.ShowMol(molecule_to_rdkit(to_plot), size=(1024, 1024), title=title)

    if timeout is None:
        plot_on_thread()
    else:
        p = multiprocessing.Process(target=plot_on_thread)
        p.start()
        time.sleep(timeout)
        p.terminate()


def plot_molecule_graph(molecule):
    """
    Plots the underlying networkx graph of the molecule.
    """
    import matplotlib.pyplot as plt

    pos = networkx.spring_layout(molecule.graph)

    for i, j in molecule.graph.edges:
        xs = pos[i][0], pos[j][0]
        ys = pos[i][1], pos[j][1]
        plt.plot(xs, ys, color="black")
        plt.annotate(str(molecule.graph.edges[i, j]["order"]), (pos[i] + pos[j]) * 0.5)

    for i in molecule.graph.nodes:
        label = str(molecule.graph.nodes[i]["valence"])
        label += f" [{i}" + ("c]" if molecule.is_cyclic(i) else "]")
        plt.annotate(label, pos[i])

    plt.show()
