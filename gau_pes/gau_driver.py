import sys
from typing import Union

import numpy as np

# constants
_ANG2BOHR = 1.8897259886


class GauDriver:
    def __init__(self) -> None:
        (layer, InputFile, OutputFile, MsgFile,
         FChkFile, MatElFile) = sys.argv[1:]

        self.layer = layer
        self.input_file = InputFile
        self.output_file = OutputFile
        self.msg_file = MsgFile
        self.fchk_file = FChkFile
        self.mat_el_file = MatElFile

        coords = []
        with open(self.input_file, "r") as f:
            (natom, derivs, charge, spin) = \
                [int(x) for x in f.readline().split()]

            self.natom = natom
            self.derivs = derivs
            self.charge = charge
            self.spin = spin

            for _ in range(self.natom):
                arr = f.readline().split()
                coord = [(float(x) / _ANG2BOHR) for x in arr[1:4]]
                coords.append(coord)

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
