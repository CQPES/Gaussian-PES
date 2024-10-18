from __future__ import annotations
from typing import Optional
from sys import argv
from re import split

import numpy as np

from gau_pes.data import PERIODIC_TABLE


class Orca6Driver:
    # ORCA input and output file
    inp_file: str
    xyz_file: str
    engrad_file: str
    pointcharges_file: str

    # ORCA input data
    charge: int
    multiplicity: int
    NCores: int
    do_gradient: bool
    
    # external infomation
    do_hessian: bool = False

    natom: int = None
    xyz_comment: str = None
    sysmbols: list[str] = None
    atoms_nuc: list[int] = None
    coords: np.ndarray = None

    @staticmethod
    def from_stdio() -> Orca6Driver:
        inp_file = argv[1]
        return Orca6Driver(inp_file)

    def __init__(self, inp_file: str) -> None:
        # ORCA input file, see ORCA Manual 6.0 sec 7.26.6
        self.inp_file = inp_file
        extinp = open(self.inp_file).readlines()
        extinp = [x.split("#")[0] for x in extinp]

        # first line is the xyz file name
        self.xyz_file = extinp[0].strip()
        # get the output file name
        self.engrad_file = self.xyz_file.replace(".xyz", ".engrad")
        # remainder line
        (self.charge, self.multiplicity, self.NCores) = \
            [int(x) for x in extinp[1:4]]
        self.do_gradient = bool(extinp[4])

        # parase xyz file
        xyzfile = open(self.xyz_file).readlines()
        self.natom = int(xyzfile[0])
        self.xyz_comment = xyzfile[1]

        coords: list[list[float]] = []
        sysmbols: list[str] = []
        atoms_nuc: list[int] = []
        for line_str in xyzfile[2:]:
            (sym, x, y, z) = split(" +", line_str.strip())
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
        if gradients is None:
            open(self.engrad_file, "w").write(header_str)
        else:
            np.savetxt(self.engrad_file, X=gradients.reshape(-1),
                       header=header_str, comments="")

    def xyz(self) -> str:
        return open(self.xyz_file).read()

    def pyscf_atom(self) -> str:
        return self.xyz_file
