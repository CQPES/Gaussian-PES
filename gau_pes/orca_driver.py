from typing import Optional

from sys import argv
import numpy as np
import re

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

class OrcaDeriver:

    inp_file: str
    xyz_file: str
    engrad_file: str
    pointcharges_file: str

    charge: int
    multiplicity: int
    NCores: int
    do_gradient: int

    natom: Optional[int] = None
    xyz_comment: Optional[str] = None
    sysmbols: Optional[list[str]] = None
    atoms_nuc: Optional[list[int]] = None
    coords: Optional[np.ndarray] = None

    def __init__(self, read_xyz: bool = False) -> None:
        self.inp_file = argv[1]

        extinp = open(self.inp_file).readlines()

        self.xyz_file = extinp[0].strip()
        self.engrad_file = self.xyz_file.replace(".xyz", ".engrad")
        (self.charge, self.multiplicity, self.NCores, self.do_gradient) = \
            [int(x) for x in extinp[1:]]
        
        if read_xyz:
            xyzfile = open(self.xyz_file).readlines()
            self.natom = int(xyzfile[0])
            self.xyz_comment = xyzfile[1]

            coords: list[list[float]] = []
            sysmbols: list[str] = []
            atoms_nuc: list[int] = []
            for line_str in xyzfile[2:]:
                (sym, x, y, z) = re.split(" +", line_str.strip())
                sysmbols.append(sym)
                atoms_nuc.append(PERIODIC_TABLE.index(sym))
                coords.append([float(x), float(y), float(z)])

            self.sysmbols = sysmbols
            self.atoms_nuc = atoms_nuc
            self.coords = np.array(coords)
        
    def write(
        self,
        energy: float,
        gradients: Optional[np.array],
    ) -> None:
        header_str = f"{self.natom}\n{energy}"
        if gradients == None:
            open(self.engrad_file, "w").write(header_str)
        else:
            np.savetxt(self.engrad_file, X=gradients.reshape(-1), header=header_str, comments="")

    def xyz(self) -> str:
        return open(self.xyz_file).read()
    
    def pyscf_atom(self) -> str:
        return self.xyz_file
    