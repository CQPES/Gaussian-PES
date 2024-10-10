import sys
from typing import Union

import numpy as np

PERIODIC_TABLE: list[str] = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
    "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
    "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
    "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk",
    "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh",
    "Fl", "Mc", "Lv", "Ts", "Og",
]

# constants
_ANG2BOHR: float = 1.8897259886


class GauDriver:
    layer: str
    input_file: str
    output_file: str
    msg_file: str
    fchk_file: str
    mat_el_file: str

    natom: int
    derivs: int
    charge: int
    spin: int

    sysmbols: list[str]
    atoms_nuc: list[int]
    coords: np.ndarray

    def __init__(self) -> None:
        (layer, InputFile, OutputFile, MsgFile,
         FChkFile, MatElFile) = sys.argv[1:]

        self.layer = layer
        self.input_file = InputFile
        self.output_file = OutputFile
        self.msg_file = MsgFile
        self.fchk_file = FChkFile
        self.mat_el_file = MatElFile

        coords: list[list[float]] = []
        atoms_nuc: list[int] = []
        sysmbols: list[str] = []
        with open(self.input_file, "r") as f:
            (natom, derivs, charge, spin) = \
                [int(x) for x in f.readline().split()]

            self.natom = natom
            self.derivs = derivs
            self.charge = charge
            self.spin = spin
            for _ in range(self.natom):
                arr = f.readline().split()
                nuc = int(arr[0])
                atoms_nuc.append(nuc)
                sysmbols.append(PERIODIC_TABLE[nuc-1])
                coord = [(float(x) / _ANG2BOHR) for x in arr[1:4]]
                coords.append(coord)

        self.atoms_nuc = atoms_nuc
        self.coords = np.array(coords)

    def write(
        self,
        energy: float,
        gradients: Union[np.array, None],
        force_constants: Union[np.array, None],
    ) -> None:
        contents = []

        # energy and dipole-monent
        line = f"{energy:20.12E}" + f"{0:20.12E}" * 3
        contents.append(line)

        # gradient on atom
        if self.derivs in [1, 2]:
            for gradient_on_atom in gradients:
                line = "".join([f"{g:20.12E}" for g in gradient_on_atom])
                contents.append(line)

        # polarizability, dipole derivatives, & force constants
        if self.derivs == 2:
            empty_line = f"{0:20.12E}" * 3

            # polarizability
            for _ in range(2):
                contents.append(empty_line)

            # dipole derivatives
            for _ in range(3 * self.natom):
                contents.append(empty_line)

            # force constants
            for component in force_constants.reshape((-1, 3)):
                line = "".join([f"{c:20.12E}" for c in component])
                contents.append(line)

        contents = "\n".join(contents)

        with open(self.output_file, "w") as f:
            f.write(contents)

    def xyz(self, comment: str = "") -> str:
        xyz_str = f"{self.natom}\n{comment}\n"
        for sym, coord in zip(self.sysmbols, self.coords):
            xyz_str += f"  {sym}{coord[0]:>23.17e}{coord[1]:>23.17e}{coord[2]:>23.17e}\n"
        return xyz_str

    def pyscf_atom(self) -> str:
        atom_str = ""
        for nuc, coord in zip(self.atoms_nuc, self.coords):
            atom_str += f"{nuc}{coord[0]:>23.17e}{coord[1]:>23.17e}{coord[2]:>23.17e};"
        return atom_str
