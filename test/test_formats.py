from pymolgen.molecule_formats import molecule_from_bru
from pymolgen.molecule_visualization import plot_molecule


def test_bru():
	mol = molecule_from_bru("atom H\natom H\nbond 0 1 1")
	assert mol.atom_count == 2
	assert len(mol.graph.edges) == 1
	assert mol.graph.nodes[0]["valence"] == 1.0
	assert mol.graph.nodes[1]["valence"] == 1.0
	assert 0 in mol.graph[1]
	assert 1 in mol.graph[0]